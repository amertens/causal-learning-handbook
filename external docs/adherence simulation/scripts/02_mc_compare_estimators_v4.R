###############################################################################
# 02_mc_compare_estimators_v4.R — Monte Carlo comparison of estimators
#
# Repeats simulation B times to assess:
#   - bias, empirical SE, RMSE, coverage
#   - naïve vs LTMLE vs LMTP
#
# Default fast settings: n=1500, B=50
# Usage:
#   Rscript scripts/02_mc_compare_estimators_v4.R          # B=50
#   Rscript scripts/02_mc_compare_estimators_v4.R 200      # B=200
###############################################################################

cat("=== V4 Monte Carlo Estimator Comparison ===\n\n")

# ── Source helpers ──────────────────────────────────────────────────────────
source(file.path("R", "sim_v4_dgp.R"))
source(file.path("R", "oracle_v4.R"))
source(file.path("R", "prepare_ltmle_v4.R"))
source(file.path("R", "prepare_lmtp_v4.R"))

for (pkg in c("ltmle", "lmtp")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(sprintf("Package '%s' required. install.packages('%s')", pkg, pkg))
}

# ── Configuration ───────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
B <- if (length(args) >= 1) as.integer(args[1]) else 50L

n_sim      <- 1500L
K          <- 4L
n_oracle   <- 50000L
seed_base  <- 20260301L
sl_lib     <- "glm"
lmtp_folds <- 1L
delta      <- 0.10

cat(sprintf("Config: n=%d, K=%d, B=%d, oracle_n=%d\n\n", n_sim, K, B, n_oracle))

# ── Output directory ────────────────────────────────────────────────────────
out_dir <- file.path("results", "mc_v4")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# ══════════════════════════════════════════════════════════════════════════════
# Step 1: Compute oracle truth (once)
# ══════════════════════════════════════════════════════════════════════════════
cat("── Computing oracle truth ──\n")

oracle <- list()

# Static category regimes by drug
for (d in 0:1) {
  dlabel <- ifelse(d == 0, "Drug_B", "Drug_A")
  for (cat_name in c("high", "mid", "low")) {
    key <- paste0("cat_", cat_name, "_", dlabel)
    oracle[[key]] <- oracle_static_category_v4(
      category = cat_name, drug = d, K = K, n_oracle = n_oracle,
      seed = seed_base + d * 100 + match(cat_name, c("high", "mid", "low"))
    )
    cat(sprintf("  %s always_%s: %.4f\n", dlabel, cat_name, oracle[[key]]))
  }
}

# Shift by drug
for (d in 0:1) {
  dlabel <- ifelse(d == 0, "Drug_B", "Drug_A")
  key <- paste0("shift_", dlabel)
  oracle[[key]] <- oracle_shift_v4(
    delta = delta, drug = d, K = K, n_oracle = n_oracle,
    seed = seed_base + 500 + d
  )
  cat(sprintf("  %s +%.0f%% shift: %.4f\n", dlabel, delta * 100, oracle[[key]]))
}

# Matched adherence drug comparison
matched <- oracle_matched_adherence_v4(
  reference_drug = 0, K = K, n_oracle = n_oracle, seed = seed_base + 600
)
oracle[["matched_Drug_B"]] <- matched$risk_ref
oracle[["matched_Drug_A"]] <- matched$risk_other
oracle[["matched_RD"]]     <- matched$rd
cat(sprintf("  Matched RD (Drug A - Drug B under Drug B adh): %.4f\n", matched$rd))

cat("\n")


# ══════════════════════════════════════════════════════════════════════════════
# Step 2: Monte Carlo loop
# ══════════════════════════════════════════════════════════════════════════════
cat("── Starting MC loop ──\n")

# Define estimands we track
estimand_names <- c(
  "cat_high_Drug_A", "cat_high_Drug_B",
  "cat_low_Drug_A",  "cat_low_Drug_B",
  "shift_Drug_A",    "shift_Drug_B",
  "matched_RD"
)

# Storage
mc_store <- data.frame(
  rep       = integer(),
  estimand  = character(),
  method    = character(),
  estimate  = numeric(),
  se        = numeric(),
  ci_lo     = numeric(),
  ci_hi     = numeric(),
  converged = logical(),
  stringsAsFactors = FALSE
)

t_start <- Sys.time()

for (b in seq_len(B)) {
  rep_seed <- seed_base + b * 13
  cat(sprintf("  Rep %d/%d (seed=%d)...", b, B, rep_seed))
  rep_t0 <- Sys.time()

  # Simulate data
  dat <- simulate_v4_data(n = n_sim, K = K, seed = rep_seed)
  p   <- attr(dat, "dgp_params")

  # ── Naïve estimates ────────────────────────────────────────────────────
  pdc_cols <- paste0("PDC_", 1:K)
  dat$mean_pdc <- rowMeans(dat[, pdc_cols], na.rm = TRUE)
  dat$pdc_strat <- cut(dat$mean_pdc,
                       breaks = c(0, p$pdc_cut_low, p$pdc_cut_high, 1),
                       labels = c("low", "mid", "high"),
                       include.lowest = TRUE)

  for (d in 0:1) {
    dlabel <- ifelse(d == 0, "Drug_B", "Drug_A")
    idx <- dat$A0 == d
    for (cat_name in c("high", "low")) {
      sidx <- idx & dat$pdc_strat == cat_name
      est <- mean(dat$Y_final[sidx], na.rm = TRUE)
      n_s <- sum(sidx, na.rm = TRUE)
      se <- sqrt(est * (1 - est) / max(n_s, 1))
      mc_store <- rbind(mc_store, data.frame(
        rep = b, estimand = paste0("cat_", cat_name, "_", dlabel),
        method = "naive", estimate = est, se = se,
        ci_lo = est - 1.96 * se, ci_hi = est + 1.96 * se,
        converged = TRUE, stringsAsFactors = FALSE
      ))
    }
  }

  # ── LTMLE estimates ────────────────────────────────────────────────────
  for (d in 0:1) {
    dlabel <- ifelse(d == 0, "Drug_B", "Drug_A")
    for (regime in c("always_high", "always_low")) {
      cat_name <- sub("always_", "", regime)
      res <- fit_ltmle_v4_static(dat, drug = d, regime = regime,
                                  K = K, sl_lib = sl_lib)
      mc_store <- rbind(mc_store, data.frame(
        rep = b, estimand = paste0("cat_", cat_name, "_", dlabel),
        method = "ltmle", estimate = res$estimate, se = res$se,
        ci_lo = res$ci_lo, ci_hi = res$ci_hi,
        converged = res$converged, stringsAsFactors = FALSE
      ))
    }
  }

  # ── LMTP shift estimates ──────────────────────────────────────────────
  for (d in 0:1) {
    dlabel <- ifelse(d == 0, "Drug_B", "Drug_A")
    res <- fit_lmtp_shift_v4(dat, delta = delta, drug = d, K = K,
                              folds = lmtp_folds)
    mc_store <- rbind(mc_store, data.frame(
      rep = b, estimand = paste0("shift_", dlabel),
      method = "lmtp_shift", estimate = res$estimate, se = res$se,
      ci_lo = res$ci_lo, ci_hi = res$ci_hi,
      converged = res$converged, stringsAsFactors = FALSE
    ))
  }

  # ── LMTP matched adherence comparison ─────────────────────────────────
  matched_res <- tryCatch({
    mres <- fit_lmtp_matched_drug_contrast_v4(
      dat, reference_drug = 0, K = K, folds = lmtp_folds
    )
    contrast <- attr(mres, "contrast")
    data.frame(
      rep = b, estimand = "matched_RD", method = "lmtp_matched",
      estimate = contrast$rd, se = contrast$se,
      ci_lo = contrast$ci_lo, ci_hi = contrast$ci_hi,
      converged = TRUE, stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      rep = b, estimand = "matched_RD", method = "lmtp_matched",
      estimate = NA_real_, se = NA_real_,
      ci_lo = NA_real_, ci_hi = NA_real_,
      converged = FALSE, stringsAsFactors = FALSE
    )
  })
  mc_store <- rbind(mc_store, matched_res)

  rep_time <- as.numeric(difftime(Sys.time(), rep_t0, units = "secs"))
  cat(sprintf(" %.1fs\n", rep_time))

  # Save checkpoint every 10 reps
  if (b %% 10 == 0 || b == B) {
    saveRDS(mc_store, file.path(out_dir, "mc_raw_v4.rds"))
    cat(sprintf("    [Checkpoint saved: %d reps]\n", b))
  }
}

total_time <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
cat(sprintf("\nMC loop complete. Total time: %.1f minutes\n\n", total_time))


# ══════════════════════════════════════════════════════════════════════════════
# Step 3: Summarise results
# ══════════════════════════════════════════════════════════════════════════════
cat("── MC Summary ──\n\n")

summary_rows <- list()

for (est_name in unique(mc_store$estimand)) {
  for (method in unique(mc_store$method[mc_store$estimand == est_name])) {
    sub <- mc_store[mc_store$estimand == est_name &
                      mc_store$method == method &
                      mc_store$converged, ]
    if (nrow(sub) == 0) next

    truth <- oracle[[est_name]]
    if (is.null(truth)) truth <- NA_real_

    mean_est <- mean(sub$estimate, na.rm = TRUE)
    bias     <- mean_est - truth
    emp_se   <- sd(sub$estimate, na.rm = TRUE)
    rmse     <- sqrt(bias^2 + emp_se^2)
    mean_se  <- mean(sub$se, na.rm = TRUE)
    se_ratio <- mean_se / emp_se

    # Coverage
    if (!is.na(truth) && all(!is.na(sub$ci_lo)) && all(!is.na(sub$ci_hi))) {
      coverage <- mean(sub$ci_lo <= truth & sub$ci_hi >= truth, na.rm = TRUE)
    } else {
      coverage <- NA_real_
    }

    n_valid <- sum(sub$converged)

    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      estimand  = est_name,
      method    = method,
      truth     = truth,
      mean_est  = mean_est,
      bias      = bias,
      emp_se    = emp_se,
      rmse      = rmse,
      mean_se   = mean_se,
      se_ratio  = se_ratio,
      coverage  = coverage,
      n_valid   = n_valid,
      stringsAsFactors = FALSE
    )
  }
}

mc_summary <- do.call(rbind, summary_rows)
rownames(mc_summary) <- NULL

# Print formatted
cat(sprintf("  %-22s %-14s %7s %7s %8s %7s %7s %7s %5s\n",
            "Estimand", "Method", "Truth", "Mean", "Bias", "EmpSE", "RMSE", "Cov95", "N"))
cat(paste(rep("-", 105), collapse = ""), "\n")
for (i in seq_len(nrow(mc_summary))) {
  r <- mc_summary[i, ]
  cat(sprintf("  %-22s %-14s %7s %7.4f %+8.4f %7.4f %7.4f %7s %5d\n",
              r$estimand, r$method,
              ifelse(is.na(r$truth), "   --  ", sprintf("%.4f", r$truth)),
              r$mean_est, r$bias, r$emp_se, r$rmse,
              ifelse(is.na(r$coverage), "  --  ", sprintf("%.3f", r$coverage)),
              r$n_valid))
}

# Save summary
write.csv(mc_summary, file.path(out_dir, "mc_summary_v4.csv"), row.names = FALSE)
saveRDS(mc_store,   file.path(out_dir, "mc_raw_v4.rds"))
saveRDS(oracle,     file.path(out_dir, "oracle_v4.rds"))

cat(sprintf("\nResults saved to %s/\n", out_dir))


# ══════════════════════════════════════════════════════════════════════════════
# Step 4: Didactic interpretation
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── Interpretation ──\n\n")

cat("Expected story from these results:\n")
cat("  1. NAIVE estimates are biased — they confound adherence with\n")
cat("     underlying health status (time-dependent confounding).\n")
cat("  2. LTMLE on categorised adherence is approximately unbiased when\n")
cat("     the GLM models are well-specified (low-dimensional DGP).\n")
cat("  3. LTMLE for 'always_low' may show higher variance or instability\n")
cat("     due to emerging positivity problems (few subjects naturally\n")
cat("     maintain low adherence across all blocks).\n")
cat("  4. LMTP shift estimates are stable because the +10% PDC shift\n")
cat("     is a feasible modification of the observed treatment process.\n")
cat("  5. LMTP matched-adherence contrast isolates the drug effect by\n")
cat("     equalising adherence distributions via quantile matching.\n")
cat("\n")

# Convergence report
cat("── Convergence report ──\n")
conv_tab <- table(mc_store$method, mc_store$converged)
print(conv_tab)

cat("\nDone.\n")
