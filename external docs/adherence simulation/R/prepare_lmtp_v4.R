###############################################################################
# prepare_lmtp_v4.R — LMTP helpers for continuous PDC analyses in V4
#
# Uses the lmtp package in wide format with:
#   - continuous PDC treatment nodes
#   - optional censoring nodes
#   - survival-type block outcomes (absorbing 0/1)
#   - mtp = TRUE for modified treatment policies
#   - default GLM learners and folds = 1 for speed
#
# IMPORTANT: lmtp with outcome_type = "survival" returns P(survival).
#   Risk = 1 - P(survival).  All functions in this file return RISK.
#
# Two main intervention strategies:
#   A. Additive feasible shift: PDC_t → min(PDC_t + delta, 1)
#   B. Quantile-based adherence matching: equalize PDC distributions across drugs
###############################################################################

if (!requireNamespace("lmtp", quietly = TRUE))
  message("Install the lmtp package: install.packages('lmtp')")


# ── Safe extraction helpers for lmtp fit objects ────────────────────────────
# lmtp >= 1.5.0 uses ife::influence_func_estimate (S7 object).
# Use lmtp::tidy(fit) which returns data.frame(estimate, std.error, conf.low, conf.high).
# The 'estimate' column is P(survival); we convert to risk = 1 - P(survival).

safe_lmtp_extract <- function(fit) {
  td <- tryCatch(lmtp::tidy(fit), error = function(e) {
    # Fallback for older lmtp versions
    tryCatch(generics::tidy(fit), error = function(e2) NULL)
  })

  if (!is.null(td) && is.data.frame(td) && nrow(td) > 0) {
    # lmtp returns P(survival); convert to risk
    return(list(
      risk    = 1 - td$estimate[1],
      se      = td$std.error[1],
      ci_lo   = 1 - td$conf.high[1],  # flip for risk
      ci_hi   = 1 - td$conf.low[1]
    ))
  }

  # Last resort: try legacy API
  theta <- if (!is.null(fit$theta) && length(fit$theta) > 0) fit$theta else NA_real_
  se    <- if (!is.null(fit$standard_error) && length(fit$standard_error) > 0) fit$standard_error else NA_real_
  list(
    risk  = 1 - as.numeric(theta),
    se    = as.numeric(se),
    ci_lo = NA_real_,
    ci_hi = NA_real_
  )
}


# ── Prepare wide data for lmtp ──────────────────────────────────────────────
# Returns a list with: data, trt_nodes, outcome_nodes, baseline, tv_confounders,
#                       cens_nodes (if applicable)
make_lmtp_data_v4 <- function(df,
                              K                 = NULL,
                              drug              = NULL,
                              include_censoring = FALSE) {

  if (is.null(K)) K <- attr(df, "K")
  if (is.null(K)) stop("Supply K or use data from simulate_v4_data().")

  # Optionally subset to one drug
  if (!is.null(drug)) {
    df <- df[df$A0 == drug, , drop = FALSE]
  }

  n <- nrow(df)

  # Build a clean wide data frame for lmtp
  # Column order matters: baseline, then for each t: L_t, [C_t], PDC_t, Y_t
  out <- data.frame(id = seq_len(n), W0 = df$W0, A0 = df$A0)

  trt_nodes  <- character(K)
  out_nodes  <- character(K)
  tv_conf    <- vector("list", K)
  cens_nodes <- character(0)

  for (t in seq_len(K)) {
    # Time-varying confounder
    lname <- paste0("L_", t)
    out[[lname]] <- df[[lname]]
    tv_conf[[t]] <- lname

    # Censoring node (lmtp convention: 1 = observed, 0 = censored)
    if (include_censoring) {
      cname <- paste0("C_", t)
      if (cname %in% names(df)) {
        out[[cname]] <- df[[cname]]
      } else {
        out[[cname]] <- rep(1L, n)
      }
      cens_nodes <- c(cens_nodes, cname)
    }

    # Treatment: continuous PDC
    pdc_name <- paste0("PDC_", t)
    out[[pdc_name]] <- df[[pdc_name]]
    trt_nodes[t] <- pdc_name

    # Outcome: absorbing failure
    yname <- paste0("Y_", t)
    out[[yname]] <- df[[yname]]
    out_nodes[t] <- yname
  }

  # Enforce absorbing state for lmtp survival: once Y_t = 1,
  # all subsequent Y must be 1 and post-event nodes can be NA
  for (t in seq_len(K)) {
    if (t < K) {
      failed <- !is.na(out[[out_nodes[t]]]) & out[[out_nodes[t]]] == 1L
      for (s in (t + 1):K) {
        out[[tv_conf[[s]]]][failed]  <- NA_real_
        out[[trt_nodes[s]]][failed]  <- NA_real_
        if (include_censoring && length(cens_nodes) >= s) {
          out[[cens_nodes[s]]][failed] <- NA_integer_
        }
        out[[out_nodes[s]]][failed]  <- 1L
      }
    }
  }

  # Baseline covariates
  baseline <- c("W0", "A0")
  if (!is.null(drug)) {
    # If stratified, A0 is constant → drop from baseline
    baseline <- "W0"
    out$A0 <- NULL
  }

  list(
    data         = out,
    trt          = trt_nodes,
    outcome      = out_nodes,
    baseline     = baseline,
    time_vary    = tv_conf,
    cens         = if (include_censoring) cens_nodes else NULL,
    K            = K,
    drug         = drug
  )
}


# ── Additive shift function factory ─────────────────────────────────────────
# Returns a shift function: PDC_t → min(PDC_t + delta, 1)
make_additive_shift_v4 <- function(delta = 0.10) {
  force(delta)
  function(data, trt) {
    pmin(data[[trt]] + delta, 1.0)
  }
}


# ── Natural (no shift) function ─────────────────────────────────────────────
# Returns natural treatment (identity)
natural_shift_v4 <- function(data, trt) {
  data[[trt]]
}


# ── Core LMTP fitting wrapper ───────────────────────────────────────────────
# Returns list(estimate=RISK, se, ci_lo, ci_hi, converged, fit, ...)
fit_lmtp_core_v4 <- function(prep, shift_fn, mtp = TRUE,
                             folds = 1,
                             learners_trt = "SL.glm",
                             learners_outcome = "SL.glm",
                             estimator = "lmtp_tmle",
                             verbose = FALSE) {

  est_fn <- switch(estimator,
                   lmtp_tmle = lmtp::lmtp_tmle,
                   lmtp_sdr  = lmtp::lmtp_sdr,
                   stop("Unknown estimator: ", estimator))

  tryCatch({
    fit <- est_fn(
      data             = prep$data,
      trt              = prep$trt,
      outcome          = prep$outcome,
      baseline         = prep$baseline,
      time_vary        = prep$time_vary,
      cens             = prep$cens,
      shift            = shift_fn,
      mtp              = mtp,
      outcome_type     = "survival",
      folds            = folds,
      learners_trt     = learners_trt,
      learners_outcome = learners_outcome
    )

    ext <- safe_lmtp_extract(fit)

    list(
      estimate  = ext$risk,
      se        = ext$se,
      ci_lo     = ext$ci_lo,
      ci_hi     = ext$ci_hi,
      fit       = fit,
      converged = TRUE
    )
  }, error = function(e) {
    if (verbose) message("LMTP failed: ", conditionMessage(e))
    list(
      estimate = NA_real_, se = NA_real_,
      ci_lo = NA_real_, ci_hi = NA_real_,
      fit = NULL, converged = FALSE,
      error = conditionMessage(e)
    )
  })
}


# ── Fit LMTP with additive shift ────────────────────────────────────────────
fit_lmtp_shift_v4 <- function(df, delta = 0.10, drug = NULL, K = 4,
                              include_censoring = FALSE, folds = 1,
                              learners_trt = "SL.glm",
                              learners_outcome = "SL.glm",
                              estimator = "lmtp_tmle", verbose = FALSE) {

  prep <- make_lmtp_data_v4(df, K = K, drug = drug,
                             include_censoring = include_censoring)
  shift_fn <- make_additive_shift_v4(delta)

  res <- fit_lmtp_core_v4(prep, shift_fn, mtp = TRUE, folds = folds,
                           learners_trt = learners_trt,
                           learners_outcome = learners_outcome,
                           estimator = estimator, verbose = verbose)
  res$delta <- delta
  res$drug  <- drug
  res
}


# ── Fit LMTP for natural (observed) treatment ───────────────────────────────
fit_lmtp_natural_v4 <- function(df, drug = NULL, K = 4,
                                include_censoring = FALSE, folds = 1,
                                learners_trt = "SL.glm",
                                learners_outcome = "SL.glm",
                                estimator = "lmtp_tmle", verbose = FALSE) {

  prep <- make_lmtp_data_v4(df, K = K, drug = drug,
                             include_censoring = include_censoring)

  res <- fit_lmtp_core_v4(prep, natural_shift_v4, mtp = TRUE, folds = folds,
                           learners_trt = learners_trt,
                           learners_outcome = learners_outcome,
                           estimator = estimator, verbose = verbose)
  res$drug <- drug
  res
}


# ── Quantile-matching shift function factory ────────────────────────────────
make_quantile_match_shift_v4 <- function(df, from_drug, to_drug, K = NULL) {

  if (is.null(K)) K <- attr(df, "K")
  if (is.null(K)) stop("Supply K")

  ecdf_list    <- list()
  quantile_list <- list()

  for (t in seq_len(K)) {
    pdc_col <- paste0("PDC_", t)
    ecdf_list[[t]]    <- list()
    quantile_list[[t]] <- list()
    for (d in 0:1) {
      vals <- df[[pdc_col]][df$A0 == d]
      vals <- vals[!is.na(vals)]
      ecdf_list[[t]][[d + 1]]    <- ecdf(vals)
      quantile_list[[t]][[d + 1]] <- vals
    }
  }

  trt_nodes <- paste0("PDC_", seq_len(K))

  function(data, trt) {
    t_idx <- match(trt, trt_nodes)
    if (is.na(t_idx)) return(data[[trt]])

    pdc_vals <- data[[trt]]
    n <- length(pdc_vals)
    shifted <- numeric(n)

    F_from      <- ecdf_list[[t_idx]][[from_drug + 1]]
    target_vals <- sort(quantile_list[[t_idx]][[to_drug + 1]])

    for (i in seq_len(n)) {
      if (is.na(pdc_vals[i])) {
        shifted[i] <- NA_real_
        next
      }
      q <- F_from(pdc_vals[i])
      shifted[i] <- quantile(target_vals, probs = q, type = 1)
    }
    shifted
  }
}


# ── Fit LMTP drug comparison under matched adherence ────────────────────────
fit_lmtp_matched_drug_contrast_v4 <- function(
    df, reference_drug = 0, K = 4, include_censoring = FALSE,
    folds = 1, learners_trt = "SL.glm", learners_outcome = "SL.glm",
    estimator = "lmtp_tmle", verbose = FALSE) {

  results <- list()
  other_drug <- 1 - reference_drug

  for (d in c(reference_drug, other_drug)) {
    if (d == reference_drug) {
      shift_fn <- natural_shift_v4
    } else {
      shift_fn <- make_quantile_match_shift_v4(
        df, from_drug = d, to_drug = reference_drug, K = K
      )
    }

    prep <- make_lmtp_data_v4(df, K = K, drug = d,
                               include_censoring = include_censoring)

    res <- fit_lmtp_core_v4(prep, shift_fn, mtp = TRUE, folds = folds,
                             learners_trt = learners_trt,
                             learners_outcome = learners_outcome,
                             estimator = estimator, verbose = verbose)

    results[[length(results) + 1]] <- data.frame(
      drug       = d,
      drug_label = ifelse(d == 0, "Drug B", "Drug A"),
      estimate   = res$estimate,
      se         = res$se,
      ci_lo      = res$ci_lo,
      ci_hi      = res$ci_hi,
      converged  = res$converged,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, results)
  rownames(out) <- NULL

  ref_est <- out$estimate[out$drug == reference_drug]
  oth_est <- out$estimate[out$drug == other_drug]
  rd <- oth_est - ref_est

  ref_se <- out$se[out$drug == reference_drug]
  oth_se <- out$se[out$drug == other_drug]
  se_rd  <- sqrt(ref_se^2 + oth_se^2)

  attr(out, "contrast") <- data.frame(
    label = sprintf("%s - %s (under %s adherence)",
                    ifelse(other_drug == 1, "Drug A", "Drug B"),
                    ifelse(reference_drug == 0, "Drug B", "Drug A"),
                    ifelse(reference_drug == 0, "Drug B", "Drug A")),
    rd    = rd, se = se_rd,
    ci_lo = rd - 1.96 * se_rd,
    ci_hi = rd + 1.96 * se_rd,
    stringsAsFactors = FALSE
  )

  out
}
