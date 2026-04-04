###############################################################################
# 01_single_run_examples_v4.R — Fast end-to-end demo of all estimators
#
# 1. Simulate a clean dataset (no censoring)
# 2. Run naïve adherence-stratified analyses (biased)
# 3. Run GLM-based LTMLE on categorised adherence
# 4. Run LMTP shift analyses on continuous PDC
# 5. Run LMTP matched-adherence drug comparison
# 6. Compare estimates to oracle truth
###############################################################################

cat("=== V4 Single-Run Examples ===\n\n")

# ── Source helpers ──────────────────────────────────────────────────────────
source(file.path("R", "sim_v4_dgp.R"))
source(file.path("R", "oracle_v4.R"))
source(file.path("R", "prepare_ltmle_v4.R"))
source(file.path("R", "prepare_lmtp_v4.R"))

# Check required packages
for (pkg in c("ltmle", "lmtp")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(sprintf("Package '%s' is required. Install with: install.packages('%s')", pkg, pkg))
}

# ── Settings ────────────────────────────────────────────────────────────────
n_sim   <- 2000
K       <- 4
seed    <- 2024
n_oracle <- 30000

cat(sprintf("Settings: n=%d, K=%d, seed=%d\n\n", n_sim, K, seed))


# ══════════════════════════════════════════════════════════════════════════════
# Step 1: Simulate clean data
# ══════════════════════════════════════════════════════════════════════════════
cat("── Step 1: Simulating data ──\n")
dat <- simulate_v4_data(n = n_sim, K = K, seed = seed, include_censoring = FALSE)
cat(sprintf("  n = %d, cumulative failure = %.1f%%\n\n",
            nrow(dat), 100 * mean(dat$Y_final)))


# ══════════════════════════════════════════════════════════════════════════════
# Step 2: Oracle truth
# ══════════════════════════════════════════════════════════════════════════════
cat("── Step 2: Computing oracle truth ──\n")

oracle_results <- list()

# Natural risks
for (d in 0:1) {
  dlabel <- ifelse(d == 0, "Drug B", "Drug A")
  r <- oracle_natural_v4(drug = d, K = K, n_oracle = n_oracle, seed = seed + 100)
  oracle_results[[paste0("natural_", dlabel)]] <- r
  cat(sprintf("  Natural %s: %.4f\n", dlabel, r))
}

# Static category regimes
for (d in 0:1) {
  dlabel <- ifelse(d == 0, "Drug B", "Drug A")
  for (cat_name in c("high", "mid", "low")) {
    r <- oracle_static_category_v4(
      category = cat_name, drug = d, K = K, n_oracle = n_oracle,
      seed = seed + 200 + d * 10 + match(cat_name, c("high", "mid", "low"))
    )
    oracle_results[[paste0(cat_name, "_", dlabel)]] <- r
    cat(sprintf("  %s, always_%s: %.4f\n", dlabel, cat_name, r))
  }
}

# Shift (+10%)
for (d in 0:1) {
  dlabel <- ifelse(d == 0, "Drug B", "Drug A")
  r <- oracle_shift_v4(delta = 0.10, drug = d, K = K,
                       n_oracle = n_oracle, seed = seed + 300 + d)
  oracle_results[[paste0("shift10_", dlabel)]] <- r
  cat(sprintf("  %s, +10%% shift: %.4f\n", dlabel, r))
}

# Matched adherence
matched <- oracle_matched_adherence_v4(
  reference_drug = 0, K = K, n_oracle = n_oracle, seed = seed + 400
)
oracle_results[["matched_Drug B"]] <- matched$risk_ref
oracle_results[["matched_Drug A"]] <- matched$risk_other
oracle_results[["matched_RD"]]     <- matched$rd
cat(sprintf("  Matched (Drug B adh): Drug B=%.4f, Drug A=%.4f, RD=%.4f\n",
            matched$risk_ref, matched$risk_other, matched$rd))

cat("\n")


# ══════════════════════════════════════════════════════════════════════════════
# Step 3: Naïve analyses (biased)
# ══════════════════════════════════════════════════════════════════════════════
cat("── Step 3: Naïve adherence-stratified analysis ──\n")
cat("  (Simply stratifying by observed mean PDC — ignores time-dependent confounding)\n\n")

pdc_cols <- paste0("PDC_", 1:K)
dat$mean_pdc <- rowMeans(dat[, pdc_cols], na.rm = TRUE)
p <- attr(dat, "dgp_params")
dat$pdc_strat <- cut(dat$mean_pdc,
                     breaks = c(0, p$pdc_cut_low, p$pdc_cut_high, 1),
                     labels = c("low", "mid", "high"),
                     include.lowest = TRUE)

naive_results <- list()
cat("  Drug | Stratum |  Risk  |  N\n")
cat("  -----|---------|--------|----\n")
for (d in 0:1) {
  dlabel <- ifelse(d == 0, "Drug B", "Drug A")
  idx <- dat$A0 == d
  for (s in c("high", "mid", "low")) {
    sidx <- idx & dat$pdc_strat == s
    risk <- mean(dat$Y_final[sidx], na.rm = TRUE)
    n_s  <- sum(sidx, na.rm = TRUE)
    naive_results[[paste0(s, "_", dlabel)]] <- risk
    cat(sprintf("  %s | %s  | %.4f | %d\n", dlabel, s, risk, n_s))
  }
}
cat("\n")


# ══════════════════════════════════════════════════════════════════════════════
# Step 4: LTMLE on categorised adherence
# ══════════════════════════════════════════════════════════════════════════════
cat("── Step 4: LTMLE (GLM-based, categorised adherence) ──\n")

ltmle_res <- fit_ltmle_v4_all(dat, K = K, sl_lib = "glm", verbose = TRUE)
print(ltmle_res[, c("drug_label", "regime", "estimate", "se", "ci_lo", "ci_hi", "converged")])
cat("\n")


# ══════════════════════════════════════════════════════════════════════════════
# Step 5: LMTP shift on continuous PDC
# ══════════════════════════════════════════════════════════════════════════════
cat("── Step 5: LMTP additive shift (+10% PDC) ──\n")

lmtp_shift_results <- list()
for (d in 0:1) {
  dlabel <- ifelse(d == 0, "Drug B", "Drug A")
  cat(sprintf("  Fitting %s...\n", dlabel))

  # Natural risk
  nat_res <- fit_lmtp_natural_v4(dat, drug = d, K = K, folds = 1)
  cat(sprintf("    Natural: %.4f (SE=%.4f)\n", nat_res$estimate, nat_res$se))

  # Shift
  shift_res <- fit_lmtp_shift_v4(dat, delta = 0.10, drug = d, K = K, folds = 1)
  cat(sprintf("    +10%% shift: %.4f (SE=%.4f)\n", shift_res$estimate, shift_res$se))

  lmtp_shift_results[[dlabel]] <- list(natural = nat_res, shift = shift_res)
}
cat("\n")


# ══════════════════════════════════════════════════════════════════════════════
# Step 6: LMTP matched-adherence drug comparison
# ══════════════════════════════════════════════════════════════════════════════
cat("── Step 6: LMTP matched-adherence drug comparison ──\n")
cat("  (Comparing drugs under Drug B's adherence distribution)\n")

matched_res <- fit_lmtp_matched_drug_contrast_v4(
  dat, reference_drug = 0, K = K, folds = 1
)
print(matched_res)
contrast <- attr(matched_res, "contrast")
cat(sprintf("\n  Risk difference: %.4f (95%% CI: %.4f to %.4f)\n",
            contrast$rd, contrast$ci_lo, contrast$ci_hi))
cat("\n")


# ══════════════════════════════════════════════════════════════════════════════
# Step 7: Summary comparison table
# ══════════════════════════════════════════════════════════════════════════════
cat("══════════════════════════════════════════════════════\n")
cat("  SUMMARY: Estimates vs Oracle Truth\n")
cat("══════════════════════════════════════════════════════\n\n")

# Build comparison table
comp <- data.frame(
  estimand = character(),
  drug     = character(),
  oracle   = numeric(),
  naive    = numeric(),
  ltmle    = numeric(),
  lmtp     = numeric(),
  stringsAsFactors = FALSE
)

# Static category regimes (LTMLE vs naive vs oracle)
for (d in 0:1) {
  dlabel <- ifelse(d == 0, "Drug B", "Drug A")
  for (cat_name in c("high", "mid", "low")) {
    regime_name <- paste0("always_", cat_name)
    oracle_val <- oracle_results[[paste0(cat_name, "_", dlabel)]]
    naive_val  <- naive_results[[paste0(cat_name, "_", dlabel)]]
    ltmle_row  <- ltmle_res[ltmle_res$drug == d & ltmle_res$regime == regime_name, ]
    ltmle_val  <- if (nrow(ltmle_row) > 0) ltmle_row$estimate else NA

    comp <- rbind(comp, data.frame(
      estimand = paste0("always_", cat_name),
      drug     = dlabel,
      oracle   = round(oracle_val, 4),
      naive    = round(naive_val, 4),
      ltmle    = round(ltmle_val, 4),
      lmtp     = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
}

# Shift estimands (LMTP vs oracle)
for (d in 0:1) {
  dlabel <- ifelse(d == 0, "Drug B", "Drug A")
  oracle_val <- oracle_results[[paste0("shift10_", dlabel)]]
  lmtp_val <- lmtp_shift_results[[dlabel]]$shift$estimate

  comp <- rbind(comp, data.frame(
    estimand = "+10% shift",
    drug     = dlabel,
    oracle   = round(oracle_val, 4),
    naive    = NA_real_,
    ltmle    = NA_real_,
    lmtp     = round(lmtp_val, 4),
    stringsAsFactors = FALSE
  ))
}

# Print
cat("  Estimand        | Drug    | Oracle | Naïve  | LTMLE  | LMTP\n")
cat("  ----------------|---------|--------|--------|--------|------\n")
for (i in seq_len(nrow(comp))) {
  r <- comp[i, ]
  cat(sprintf("  %-16s| %-8s| %s | %s | %s | %s\n",
              r$estimand, r$drug,
              ifelse(is.na(r$oracle), "  --  ", sprintf("%.4f", r$oracle)),
              ifelse(is.na(r$naive),  "  --  ", sprintf("%.4f", r$naive)),
              ifelse(is.na(r$ltmle),  "  --  ", sprintf("%.4f", r$ltmle)),
              ifelse(is.na(r$lmtp),   "  --  ", sprintf("%.4f", r$lmtp))))
}

# Matched adherence comparison
cat("\n  Drug comparison under matched adherence (Drug B distribution):\n")
cat(sprintf("    Oracle RD: %.4f\n", oracle_results[["matched_RD"]]))
cat(sprintf("    LMTP   RD: %.4f (SE=%.4f)\n", contrast$rd, contrast$se))

cat("\n── Bias summary ──\n")
cat("  Naïve bias = naive - oracle (expected: biased due to confounding)\n")
cat("  LTMLE bias = ltmle - oracle (expected: small if models correct)\n")
for (i in seq_len(nrow(comp))) {
  r <- comp[i, ]
  if (!is.na(r$naive) && !is.na(r$oracle)) {
    cat(sprintf("  %s %s: naive bias = %+.4f\n",
                r$drug, r$estimand, r$naive - r$oracle))
  }
  if (!is.na(r$ltmle) && !is.na(r$oracle)) {
    cat(sprintf("  %s %s: LTMLE bias = %+.4f\n",
                r$drug, r$estimand, r$ltmle - r$oracle))
  }
  if (!is.na(r$lmtp) && !is.na(r$oracle)) {
    cat(sprintf("  %s %s: LMTP  bias = %+.4f\n",
                r$drug, r$estimand, r$lmtp - r$oracle))
  }
}

# Support diagnostics
cat("\n── Support diagnostics ──\n")
diagnose_ltmle_support_v4(dat, K = K)

cat("\n\nDone. Single-run examples complete.\n")
