###############################################################################
# 04_mc_phase2_comparison.R — Slide-ready method comparison MC (Phase 2)
#
# Purpose: generate bias-comparison figure for the talk
# Default: n=2000, B=300
# Methods: naive, naive g-comp, LTMLE (GLM), LMTP (+10% shift)
# Estimands: always_low, always_high, shift
# Output: summary CSV + dot-whisker bias figure
###############################################################################

cat("=== MC Phase 2: Method Comparison ===\n")

source(file.path("R", "sim_v4_dgp.R"))
source(file.path("R", "oracle_v4.R"))
source(file.path("R", "prepare_ltmle_v4.R"))
source(file.path("R", "prepare_lmtp_v4.R"))

for (pkg in c("ltmle", "lmtp"))
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(sprintf("Package '%s' required", pkg))

args <- commandArgs(trailingOnly = TRUE)
B      <- if (length(args) >= 1) as.integer(args[1]) else 300L
n_sim  <- 2000L
K      <- 4L
seed0  <- 20260501L

out_dir <- file.path("results", "mc_phase2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Oracle ──────────────────────────────────────────────────────────────────
cat("Computing oracle...\n")
p <- default_params_v4()
oracle <- list()
for (d in 0:1) {
  dl <- ifelse(d == 0, "Drug_B", "Drug_A")
  for (cn in c("low", "high")) {
    oracle[[paste0("cat_", cn, "_", dl)]] <-
      oracle_static_category_v4(cn, drug = d, K = K, n_oracle = 100000, seed = seed0 + d*10 + match(cn,c("low","high")))
  }
  oracle[[paste0("shift_", dl)]] <-
    oracle_shift_v4(0.10, drug = d, K = K, n_oracle = 100000, seed = seed0 + 100 + d)
}

# ── MC loop ─────────────────────────────────────────────────────────────────
store <- data.frame(rep = integer(), estimand = character(), method = character(),
                    estimate = numeric(), stringsAsFactors = FALSE)

t0 <- Sys.time()
for (b in seq_len(B)) {
  if (b %% 25 == 0 || b == 1) cat(sprintf("  Rep %d/%d\n", b, B))
  dat <- simulate_v4_data(n = n_sim, K = K, seed = seed0 + b * 17)

  pdc_cols <- paste0("PDC_", 1:K)
  dat$mean_pdc <- rowMeans(dat[, pdc_cols])
  dat$strat <- cut(dat$mean_pdc, c(0, p$pdc_cut_low, p$pdc_cut_high, 1),
                   labels = c("low","mid","high"), include.lowest = TRUE)

  for (d in 0:1) {
    dl <- ifelse(d == 0, "Drug_B", "Drug_A")
    idx <- dat$A0 == d

    # Naive: always_low and always_high
    for (cn in c("low", "high")) {
      nv <- mean(dat$Y_final[idx & dat$strat == cn], na.rm = TRUE)
      store <- rbind(store, data.frame(
        rep = b, estimand = paste0("cat_", cn, "_", dl),
        method = "naive", estimate = nv, stringsAsFactors = FALSE))
    }

    # LTMLE
    for (reg in c("always_low", "always_high")) {
      cn <- sub("always_", "", reg)
      res <- fit_ltmle_v4_static(dat, drug = d, regime = reg, K = K,
                                  sl_lib = "glm", verbose = FALSE)
      if (res$converged) {
        store <- rbind(store, data.frame(
          rep = b, estimand = paste0("cat_", cn, "_", dl),
          method = "ltmle", estimate = res$estimate, stringsAsFactors = FALSE))
      }
    }

    # LMTP shift
    res <- fit_lmtp_shift_v4(dat, delta = 0.10, drug = d, K = K, folds = 1)
    if (res$converged) {
      store <- rbind(store, data.frame(
        rep = b, estimand = paste0("shift_", dl),
        method = "lmtp_shift", estimate = res$estimate, stringsAsFactors = FALSE))
    }
  }

  if (b %% 50 == 0) saveRDS(store, file.path(out_dir, "phase2_raw.rds"))
}
elapsed <- difftime(Sys.time(), t0, units = "mins")
cat(sprintf("MC loop done in %.1f min\n", as.numeric(elapsed)))

# ── Summary + figure ────────────────────────────────────────────────────────
saveRDS(store, file.path(out_dir, "phase2_raw.rds"))

# Compute bias per rep
store$truth <- sapply(store$estimand, function(e) {
  v <- oracle[[e]]; if (is.null(v)) NA_real_ else v
})
store$bias <- store$estimate - store$truth

summary_rows <- list()
for (est in unique(store$estimand)) {
  for (m in unique(store$method[store$estimand == est])) {
    sub <- store[store$estimand == est & store$method == m & !is.na(store$bias), ]
    if (nrow(sub) == 0) next
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      estimand = est, method = m,
      truth = sub$truth[1],
      mean_est = mean(sub$estimate), bias = mean(sub$bias),
      emp_se = sd(sub$estimate), n_reps = nrow(sub),
      stringsAsFactors = FALSE)
  }
}
mc_summary <- do.call(rbind, summary_rows)
write.csv(mc_summary, file.path(out_dir, "phase2_summary.csv"), row.names = FALSE)

# Bias figure
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  pal <- c("naive" = "#D73027", "ltmle" = "#4575B4", "lmtp_shift" = "#1A9850")
  fig <- ggplot(store[!is.na(store$bias), ],
                aes(x = estimand, y = bias, fill = method)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = pal) +
    labs(x = "Estimand", y = "Bias (estimate - oracle)", fill = "Method",
         title = sprintf("MC Bias Comparison (B=%d, n=%d)", B, n_sim)) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(out_dir, "phase2_bias_boxplot.pdf"), fig, width = 10, height = 6)
  cat("Figure saved to", file.path(out_dir, "phase2_bias_boxplot.pdf"), "\n")
}

cat("\n=== Phase 2 Summary ===\n")
print(mc_summary)
