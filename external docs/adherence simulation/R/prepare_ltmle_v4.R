###############################################################################
# prepare_ltmle_v4.R — LTMLE helpers for categorised adherence in V4
#
# Uses 3-level adherence categories (low / mid / high) encoded as two
# binary dummy treatment nodes per block (proven pattern from V3):
#   Acat2_t = 1{PDC_cat == "mid"}
#   Acat3_t = 1{PDC_cat == "high"}
#   Reference (both 0) = "low"
#
# Node ordering per block:  [C_t] → L_t → Acat2_t → Acat3_t → Y_t
#
# Nuisance models use EXPLICIT simplified Qform/gform to avoid
# rank-deficiency from the default full-history formulas.
###############################################################################

if (!requireNamespace("ltmle", quietly = TRUE))
  message("Install the ltmle package: install.packages('ltmle')")


# ── Prepare wide data for ltmle ─────────────────────────────────────────────
make_ltmle_data_v4 <- function(df, K = NULL, include_censoring = FALSE) {

  if (is.null(K)) K <- attr(df, "K")
  if (is.null(K)) stop("Supply K or use data from simulate_v4_data().")

  n <- nrow(df)

  # Start with baseline
  out <- data.frame(W0 = df$W0, A0 = df$A0)

  for (t in seq_len(K)) {
    if (include_censoring) {
      cname <- paste0("C_", t)
      if (cname %in% names(df)) {
        out[[cname]] <- df[[cname]]
      } else {
        out[[cname]] <- rep(1L, n)
      }
    }

    # Time-varying confounder
    out[[paste0("L_", t)]] <- df[[paste0("L_", t)]]

    # Dummy-coded treatment nodes
    pdc_int <- df[[paste0("PDC_int_", t)]]
    out[[paste0("Acat2_", t)]] <- as.integer(!is.na(pdc_int) & pdc_int == 2L)
    out[[paste0("Acat3_", t)]] <- as.integer(!is.na(pdc_int) & pdc_int == 3L)
    out[[paste0("Acat2_", t)]][is.na(pdc_int)] <- NA_integer_
    out[[paste0("Acat3_", t)]][is.na(pdc_int)] <- NA_integer_

    # Outcome
    out[[paste0("Y_", t)]] <- df[[paste0("Y_", t)]]
  }

  # Enforce absorbing state: once Y_t = 1, subsequent L/A nodes → NA, Y → 1
  for (t in seq_len(K)) {
    if (t < K) {
      failed_by_t <- !is.na(out[[paste0("Y_", t)]]) & out[[paste0("Y_", t)]] == 1L
      for (s in (t + 1):K) {
        if (include_censoring) {
          out[[paste0("C_", s)]][failed_by_t] <- NA_integer_
        }
        out[[paste0("L_", s)]][failed_by_t]     <- NA_real_
        out[[paste0("Acat2_", s)]][failed_by_t]  <- NA_integer_
        out[[paste0("Acat3_", s)]][failed_by_t]  <- NA_integer_
        out[[paste0("Y_", s)]][failed_by_t]      <- 1L
      }
    }
  }

  out
}


# ── Node lists for ltmle ────────────────────────────────────────────────────
make_ltmle_nodes_v4 <- function(K, include_censoring = FALSE) {
  Lnodes <- c("W0", "A0")
  Anodes <- character()
  Cnodes <- character()
  Ynodes <- character()

  for (t in seq_len(K)) {
    if (include_censoring) {
      Cnodes <- c(Cnodes, paste0("C_", t))
    }
    Lnodes <- c(Lnodes, paste0("L_", t))
    Anodes <- c(Anodes, paste0("Acat2_", t), paste0("Acat3_", t))
    Ynodes <- c(Ynodes, paste0("Y_", t))
  }

  list(Lnodes = Lnodes, Anodes = Anodes, Cnodes = Cnodes, Ynodes = Ynodes)
}


# ── Explicit Qform and gform ────────────────────────────────────────────────
# Simplified formulas avoid rank deficiency from default full-history models.
# Q-model for Y_t: depends only on W0, current L_t, current Acat dummies
# g-model for Acat2_t / Acat3_t: depends on W0, current L_t, previous category
make_ltmle_formulas_v4 <- function(K) {
  Qform <- NULL
  gform <- NULL

  q_list <- character()
  g_list <- character()

  for (t in seq_len(K)) {
    L_t    <- paste0("L_", t)
    A2_t   <- paste0("Acat2_", t)
    A3_t   <- paste0("Acat3_", t)
    Y_t    <- paste0("Y_", t)

    # Q-model: Y_t ~ W0 + L_t + Acat2_t + Acat3_t
    q_list[Y_t] <- sprintf("Q.kplus1 ~ W0 + %s + %s + %s", L_t, A2_t, A3_t)

    # g-model: condition on W0, current L_t, and previous category (if t > 1)
    if (t == 1) {
      g_rhs <- sprintf("W0 + %s", L_t)
    } else {
      A2_prev <- paste0("Acat2_", t - 1)
      A3_prev <- paste0("Acat3_", t - 1)
      g_rhs <- sprintf("W0 + %s + %s + %s", L_t, A2_prev, A3_prev)
    }
    g_list[A2_t] <- sprintf("%s ~ %s", A2_t, g_rhs)
    g_list[A3_t] <- sprintf("%s ~ %s + %s", A3_t, g_rhs, A2_t)
  }

  list(Qform = q_list, gform = g_list)
}


# ── Deterministic g-function ────────────────────────────────────────────────
# If Acat2_t = 1, then P(Acat3_t = 1) = 0.
make_deterministic_g_v4 <- function(K) {
  a3_to_a2 <- list()
  for (t in seq_len(K)) {
    a3_to_a2[[paste0("Acat3_", t)]] <- paste0("Acat2_", t)
  }

  function(data, current.node, nodes) {
    # current.node is a column INDEX into the data frame.
    # 'nodes' is a LIST of index vectors (A, C, L, Y, ...), NOT a character vector.
    # The correct way to get the column name is names(data)[current.node].
    node_name <- names(data)[current.node]

    if (node_name %in% names(a3_to_a2)) {
      acat2_col <- a3_to_a2[[node_name]]
      if (acat2_col %in% names(data)) {
        is_mid <- !is.na(data[[acat2_col]]) & data[[acat2_col]] == 1L
        if (any(is_mid)) {
          return(list(
            is.deterministic = is_mid,
            prob1 = 0   # scalar: P(Acat3=1) = 0 when Acat2=1
          ))
        }
      }
    }
    NULL
  }
}


# ── Build abar for static category regimes ──────────────────────────────────
make_static_regimes_v4 <- function(K) {
  list(
    always_high = rep(c(0L, 1L), K),   # Acat2=0, Acat3=1 at each block
    always_mid  = rep(c(1L, 0L), K),   # Acat2=1, Acat3=0 at each block
    always_low  = rep(c(0L, 0L), K)    # both 0 at each block (reference)
  )
}


# ── Fit LTMLE for one static regime, within one drug stratum ────────────────
fit_ltmle_v4_static <- function(df,
                                drug,
                                regime       = "always_high",
                                K            = 4,
                                include_censoring = FALSE,
                                sl_lib       = "glm",
                                verbose      = TRUE) {

  # Subset to one drug
  idx <- df$A0 == drug
  df_sub <- df[idx, , drop = FALSE]

  # Prepare ltmle data
  ldata <- make_ltmle_data_v4(df_sub, K = K,
                               include_censoring = include_censoring)
  # Remove A0 from ltmle data (it's constant within stratum)
  ldata$A0 <- NULL

  nodes <- make_ltmle_nodes_v4(K, include_censoring)
  nodes$Lnodes <- setdiff(nodes$Lnodes, "A0")

  regimes <- make_static_regimes_v4(K)
  abar <- regimes[[regime]]
  if (is.null(abar)) stop("Unknown regime: ", regime)

  det_g <- make_deterministic_g_v4(K)

  # Explicit simplified formulas
  forms <- make_ltmle_formulas_v4(K)

  # Fit
  tryCatch({
    fit <- ltmle::ltmle(
      data       = ldata,
      Anodes     = nodes$Anodes,
      Cnodes     = if (include_censoring) nodes$Cnodes else NULL,
      Lnodes     = nodes$Lnodes,
      Ynodes     = nodes$Ynodes,
      abar       = abar,
      Qform      = forms$Qform,
      gform      = forms$gform,
      SL.library = sl_lib,
      estimate.time = FALSE,
      survivalOutcome = (K > 1),
      deterministic.g.function = det_g,
      variance.method = "ic"
    )

    # Extract treatment-specific mean (risk under regime)
    summ     <- summary(fit)
    tmle_est <- as.numeric(fit$estimates["tmle"])
    iptw_est <- as.numeric(fit$estimates["iptw"])
    se       <- summ$treatment$std.dev

    if (verbose) {
      cat(sprintf("  LTMLE drug=%d %s: est=%.4f (SE=%.4f)\n",
                  drug, regime, tmle_est, se))
    }

    list(
      estimate  = tmle_est,
      se        = se,
      ci_lo     = tmle_est - 1.96 * se,
      ci_hi     = tmle_est + 1.96 * se,
      iptw      = iptw_est,
      regime    = regime,
      drug      = drug,
      fit       = fit,
      converged = TRUE
    )
  }, error = function(e) {
    msg <- conditionMessage(e)
    if (verbose) message("LTMLE failed for drug=", drug, " regime=", regime,
                         ": ", msg)
    list(
      estimate = NA_real_, se = NA_real_,
      ci_lo = NA_real_, ci_hi = NA_real_,
      iptw = NA_real_,
      regime = regime, drug = drug, fit = NULL, converged = FALSE,
      error = msg
    )
  })
}


# ── Convenience: fit all static regimes for both drugs ──────────────────────
fit_ltmle_v4_all <- function(df, K = 4, include_censoring = FALSE,
                             sl_lib = "glm", verbose = FALSE) {

  regimes <- c("always_high", "always_mid", "always_low")
  drugs   <- 0:1

  results <- list()
  for (d in drugs) {
    for (reg in regimes) {
      res <- fit_ltmle_v4_static(
        df, drug = d, regime = reg, K = K,
        include_censoring = include_censoring,
        sl_lib = sl_lib, verbose = verbose
      )
      results[[length(results) + 1]] <- data.frame(
        drug       = d,
        drug_label = ifelse(d == 0, "Drug B", "Drug A"),
        regime     = reg,
        estimate   = res$estimate,
        se         = res$se,
        ci_lo      = res$ci_lo,
        ci_hi      = res$ci_hi,
        iptw       = res$iptw,
        converged  = res$converged,
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, results)
}


# ── Support / positivity diagnostics ────────────────────────────────────────
diagnose_ltmle_support_v4 <- function(df, K = NULL) {

  if (is.null(K)) K <- attr(df, "K")
  if (is.null(K)) stop("Supply K or use data from simulate_v4_data().")

  cat("\n=== LTMLE Support Diagnostics (V4) ===\n\n")

  cat("── Category frequencies by block and drug ──\n")
  for (t in seq_len(K)) {
    cat_col <- paste0("PDC_cat_", t)
    if (!cat_col %in% names(df)) next
    tab <- table(Drug = ifelse(df$A0 == 1, "Drug A", "Drug B"),
                 Category = df[[cat_col]], useNA = "ifany")
    cat(sprintf("\nBlock %d:\n", t))
    print(tab)
    cat("  Row proportions:\n")
    print(round(prop.table(tab, margin = 1), 3))
  }

  cat("\n── Subjects following static regimes for ALL blocks ──\n")
  for (d in 0:1) {
    idx <- df$A0 == d
    n_d <- sum(idx)
    dlabel <- ifelse(d == 0, "Drug B", "Drug A")

    for (regime in c("always_high", "always_mid", "always_low")) {
      target <- switch(regime,
                       always_high = "high",
                       always_mid  = "mid",
                       always_low  = "low")
      follows <- rep(TRUE, n_d)
      for (t in seq_len(K)) {
        cat_col <- paste0("PDC_cat_", t)
        vals <- as.character(df[[cat_col]][idx])
        follows <- follows & !is.na(vals) & vals == target
      }
      n_follow <- sum(follows)
      cat(sprintf("  %s, %s: %d / %d (%.1f%%)\n",
                  dlabel, regime, n_follow, n_d, 100 * n_follow / n_d))
    }
  }

  cat("\n── PDC summary by block and drug ──\n")
  for (t in seq_len(K)) {
    pdc_col <- paste0("PDC_", t)
    if (!pdc_col %in% names(df)) next
    for (d in 0:1) {
      vals <- df[[pdc_col]][df$A0 == d]
      vals <- vals[!is.na(vals)]
      dlabel <- ifelse(d == 0, "Drug B", "Drug A")
      cat(sprintf("  Block %d, %s: mean=%.3f, sd=%.3f, min=%.3f, max=%.3f, %%<0.05=%.1f%%, %%>0.95=%.1f%%\n",
                  t, dlabel, mean(vals), sd(vals), min(vals), max(vals),
                  100 * mean(vals < 0.05), 100 * mean(vals > 0.95)))
    }
  }

  invisible(NULL)
}
