###############################################################################
# diagnostics_v4.R — Targeted diagnostic figures and tables for V4
#
# Provides:
#   - fig_support_heatmap_v4()     [FIGURE A: positivity heatmap]
#   - fig_horizon_validation_v4()  [FIGURE B: horizon-by-horizon LTMLE]
#   - fig_pdc_by_drug_block_v4()   [FIGURE C: PDC distributions]
#   - fig_shift_distribution_v4()  [FIGURE D: observed vs shifted PDC]
#   - tbl_support_detail_v4()      [TABLE A: effective support diagnostics]
#   - tbl_estimand_alignment_v4()  [TABLE B: oracle midpoint vs uniform]
#   - fit_ltmle_horizons_v4()      [horizon sweep with status tracking]
###############################################################################

if (!requireNamespace("ggplot2", quietly = TRUE))
  stop("ggplot2 required for diagnostics")


# ── FIGURE A: Static regime support heatmap ─────────────────────────────────
fig_support_heatmap_v4 <- function(df, K = NULL) {
  if (is.null(K)) K <- attr(df, "K")
  library(ggplot2)

  rows <- data.frame(drug = character(), regime = character(),
                     pct = numeric(), stringsAsFactors = FALSE)
  for (d in 0:1) {
    dl <- ifelse(d == 0, "Drug B", "Drug A")
    idx <- df$A0 == d
    nd <- sum(idx)
    for (tgt in c("high", "mid", "low")) {
      follows <- rep(TRUE, nd)
      for (t in seq_len(K)) {
        v <- as.character(df[[paste0("PDC_cat_", t)]][idx])
        follows <- follows & !is.na(v) & v == tgt
      }
      rows <- rbind(rows, data.frame(
        drug = dl, regime = paste0("Always ", tgt),
        pct = 100 * sum(follows) / nd, stringsAsFactors = FALSE))
    }
  }

  rows$regime <- factor(rows$regime,
                        levels = c("Always low", "Always mid", "Always high"))
  rows$drug <- factor(rows$drug, levels = c("Drug A", "Drug B"))
  rows$label <- sprintf("%.1f%%", rows$pct)
  rows$warn <- ifelse(rows$pct < 5, "< 5%", "")

  ggplot(rows, aes(x = drug, y = regime, fill = pct)) +
    geom_tile(colour = "white", linewidth = 1.5) +
    geom_text(aes(label = label), size = 5, fontface = "bold") +
    geom_text(aes(label = warn), vjust = 2.2, size = 3, colour = "firebrick") +
    scale_fill_gradient2(low = "#D73027", mid = "#FEE08B", high = "#1A9850",
                         midpoint = 10, limits = c(0, NA),
                         name = "% following\nall 4 blocks") +
    labs(x = NULL, y = NULL,
         title = "Static Regime Support (Positivity Diagnostic)") +
    theme_minimal(base_size = 14) +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 13))
}


# ── Horizon-by-horizon LTMLE sweep with status tracking ────────────────────
fit_ltmle_horizons_v4 <- function(df, drug, regime = "always_low",
                                  K_max = 4, sl_lib = "glm") {
  results <- data.frame(
    horizon = integer(), estimate = numeric(), se = numeric(),
    oracle_u = numeric(), oracle_m = numeric(),
    bias_u = numeric(), converged = logical(),
    error_msg = character(), stringsAsFactors = FALSE
  )

  for (h in seq_len(K_max)) {
    oracle_u <- oracle_static_category_v4(
      category = sub("always_", "", regime), drug = drug, K = h,
      n_oracle = 100000, seed = 7777 + h)
    oracle_m <- oracle_static_category_midpoint_v4(
      category = sub("always_", "", regime), drug = drug, K = h,
      n_oracle = 100000, seed = 8888 + h)

    res <- fit_ltmle_v4_static(df, drug = drug, regime = regime,
                                K = h, sl_lib = sl_lib, verbose = FALSE)

    est <- if (res$converged) res$estimate else NA_real_
    se  <- if (res$converged) res$se else NA_real_
    bu  <- if (!is.na(est)) est - oracle_u else NA_real_
    emsg <- if (!is.null(res$error)) res$error else ""

    results <- rbind(results, data.frame(
      horizon = h, estimate = est, se = se,
      oracle_u = oracle_u, oracle_m = oracle_m,
      bias_u = bu, converged = res$converged,
      error_msg = emsg, stringsAsFactors = FALSE))
  }
  results
}


# ── FIGURE B: Horizon validation plot ───────────────────────────────────────
fig_horizon_validation_v4 <- function(horizon_results, title = NULL) {
  library(ggplot2)
  hr <- horizon_results

  if (is.null(title))
    title <- "LTMLE Estimate vs Oracle by Horizon"

  p <- ggplot(hr, aes(x = horizon)) +
    geom_line(aes(y = oracle_u), colour = "black", linetype = "dashed",
              linewidth = 0.8) +
    geom_point(aes(y = oracle_u), colour = "black", size = 2, shape = 4, stroke = 1.2)

  # Plot successful estimates
  ok <- hr[hr$converged & !is.na(hr$estimate), ]
  fail <- hr[!hr$converged | is.na(hr$estimate), ]

  if (nrow(ok) > 0) {
    p <- p +
      geom_point(data = ok, aes(y = estimate), colour = "#2166AC", size = 3.5) +
      geom_errorbar(data = ok, aes(ymin = estimate - 1.96 * se,
                                    ymax = estimate + 1.96 * se),
                    colour = "#2166AC", width = 0.12)
  }
  if (nrow(fail) > 0) {
    p <- p +
      geom_point(data = fail, aes(y = oracle_u), colour = "firebrick",
                 size = 5, shape = 4, stroke = 2) +
      geom_text(data = fail, aes(y = oracle_u, label = "FAILED"),
                vjust = -1.5, colour = "firebrick", size = 3, fontface = "bold")
  }

  p + scale_x_continuous(breaks = hr$horizon, labels = paste0("Q", hr$horizon)) +
    scale_y_continuous(labels = scales::percent_format(1)) +
    labs(x = "Horizon (blocks)", y = "Cumulative Failure Risk",
         title = title,
         subtitle = "Blue = LTMLE estimate (95% CI), Black X = oracle, Red X = failed fit") +
    theme_minimal(base_size = 13)
}


# ── FIGURE C: PDC distribution by drug and block ───────────────────────────
fig_pdc_by_drug_block_v4 <- function(df, K = NULL) {
  if (is.null(K)) K <- attr(df, "K")
  library(ggplot2)

  long <- do.call(rbind, lapply(seq_len(K), function(t) {
    data.frame(block = paste0("Q", t),
               drug = ifelse(df$A0 == 1, "Drug A", "Drug B"),
               PDC = df[[paste0("PDC_", t)]], stringsAsFactors = FALSE)
  }))
  long <- long[!is.na(long$PDC), ]
  long$drug <- factor(long$drug, levels = c("Drug A", "Drug B"))

  ggplot(long, aes(x = drug, y = PDC, fill = drug)) +
    geom_violin(alpha = 0.5, trim = FALSE) +
    geom_boxplot(width = 0.15, outlier.size = 0.5, alpha = 0.7) +
    facet_wrap(~ block, nrow = 1) +
    scale_y_continuous(labels = scales::percent_format(1)) +
    scale_fill_manual(values = c("Drug A" = "#2166AC", "Drug B" = "#B2182B")) +
    labs(x = NULL, y = "PDC", fill = "Drug",
         title = "Adherence (PDC) Distribution by Drug and Quarter") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom")
}


# ── FIGURE D: Observed vs shifted PDC ──────────────────────────────────────
fig_shift_distribution_v4 <- function(df, delta = 0.10, K = NULL) {
  if (is.null(K)) K <- attr(df, "K")
  library(ggplot2)

  # Use block 1 PDC for illustration
  pdc_obs <- df$PDC_1[!is.na(df$PDC_1)]
  pdc_shift <- pmin(pdc_obs + delta, 1.0)

  long <- rbind(
    data.frame(PDC = pdc_obs, Policy = "Observed", stringsAsFactors = FALSE),
    data.frame(PDC = pdc_shift,
               Policy = sprintf("+%g%% shift", delta * 100),
               stringsAsFactors = FALSE))
  long$Policy <- factor(long$Policy, levels = c("Observed", sprintf("+%g%% shift", delta * 100)))

  ggplot(long, aes(x = PDC, fill = Policy)) +
    geom_density(alpha = 0.5, linewidth = 0.6) +
    scale_x_continuous(labels = scales::percent_format(1)) +
    scale_fill_manual(values = c("Observed" = "#D73027",
                                  setNames("#1A9850", sprintf("+%g%% shift", delta * 100)))) +
    labs(x = "PDC (Block 1)", y = "Density", fill = NULL,
         title = sprintf("Additive +%g%% PDC Shift: Observed vs Shifted", delta * 100)) +
    theme_minimal(base_size = 13) + theme(legend.position = "bottom")
}


# ── TABLE A: Effective support diagnostics ──────────────────────────────────
tbl_support_detail_v4 <- function(df, K = NULL, threshold = 5) {
  if (is.null(K)) K <- attr(df, "K")

  rows <- list()
  for (d in 0:1) {
    dl <- ifelse(d == 0, "Drug B", "Drug A")
    idx <- df$A0 == d
    nd <- sum(idx)
    for (tgt in c("high", "mid", "low")) {
      # Per-block category frequency
      block_pcts <- numeric(K)
      for (t in seq_len(K)) {
        v <- as.character(df[[paste0("PDC_cat_", t)]][idx])
        block_pcts[t] <- 100 * mean(v == tgt, na.rm = TRUE)
      }

      # Full regime followers
      follows <- rep(TRUE, nd)
      for (t in seq_len(K)) {
        v <- as.character(df[[paste0("PDC_cat_", t)]][idx])
        follows <- follows & !is.na(v) & v == tgt
      }
      n_follow <- sum(follows)
      pct_follow <- 100 * n_follow / nd

      # Events among followers
      events <- sum(df$Y_final[idx][follows], na.rm = TRUE)

      flag <- if (pct_follow < threshold) "SPARSE" else "OK"

      rows[[length(rows) + 1]] <- data.frame(
        Drug = dl, Regime = paste0("always_", tgt),
        Block_pcts = paste(sprintf("%.0f%%", block_pcts), collapse = " / "),
        N_follow = n_follow, Pct_follow = sprintf("%.1f%%", pct_follow),
        Events = events, Flag = flag,
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}


# ── TABLE B: Estimand alignment — midpoint vs uniform oracle ────────────────
tbl_estimand_alignment_v4 <- function(K = 4, n_oracle = 100000, seed = 9999) {
  rows <- list()
  for (d in 0:1) {
    dl <- ifelse(d == 0, "Drug B", "Drug A")
    for (cat_name in c("low", "mid", "high")) {
      ou <- oracle_static_category_v4(cat_name, drug = d, K = K,
                                       n_oracle = n_oracle,
                                       seed = seed + d * 100 + match(cat_name, c("low","mid","high")))
      om <- oracle_static_category_midpoint_v4(cat_name, drug = d, K = K,
                                                n_oracle = n_oracle,
                                                seed = seed + d * 100 + match(cat_name, c("low","mid","high")) + 50)
      rows[[length(rows) + 1]] <- data.frame(
        Drug = dl, Category = cat_name,
        Oracle_uniform = ou, Oracle_midpoint = om,
        Difference = ou - om,
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}
