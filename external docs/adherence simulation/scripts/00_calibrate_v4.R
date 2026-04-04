###############################################################################
# 00_calibrate_v4.R — Calibration & diagnostics for the V4 DGP
#
# Run this script repeatedly while tuning parameters in default_params_v4().
# It simulates one moderately large dataset and prints/plots diagnostics
# to verify the qualitative features of the DGP.
###############################################################################

# ── Source DGP ──────────────────────────────────────────────────────────────
source(file.path("R", "sim_v4_dgp.R"))
source(file.path("R", "oracle_v4.R"))

# ── Settings ────────────────────────────────────────────────────────────────
n_cal  <- 5000
K      <- 4
seed   <- 42

# Override any params here for quick tuning:
# my_params <- default_params_v4()
# my_params$alpha0 <- 1.4
# dat <- simulate_v4_data(n = n_cal, K = K, seed = seed, params = my_params)

dat <- simulate_v4_data(n = n_cal, K = K, seed = seed)
p   <- attr(dat, "dgp_params")

cat("\n====================================================\n")
cat("  V4 DGP CALIBRATION REPORT\n")
cat("  n =", n_cal, "  K =", K, "  seed =", seed, "\n")
cat("====================================================\n")


# ── 1. Overall summaries ────────────────────────────────────────────────────
cat("\n── 1. Sample size by drug ──\n")
print(table(Drug = ifelse(dat$A0 == 1, "Drug A", "Drug B")))

cat("\n── 2. Cumulative failure by drug ──\n")
for (d in 0:1) {
  dlabel <- ifelse(d == 0, "Drug B", "Drug A")
  idx <- dat$A0 == d
  risk <- mean(dat$Y_final[idx], na.rm = TRUE)
  cat(sprintf("  %s: %.1f%% failed by block %d\n", dlabel, 100 * risk, K))
}


# ── 2. PDC distributions by drug and block ──────────────────────────────────
cat("\n── 3. PDC summary by drug and block ──\n")
for (t in 1:K) {
  pdc_col <- paste0("PDC_", t)
  for (d in 0:1) {
    vals <- dat[[pdc_col]][dat$A0 == d]
    vals <- vals[!is.na(vals)]
    dlabel <- ifelse(d == 0, "Drug B", "Drug A")
    cat(sprintf("  Block %d, %s: mean=%.3f  sd=%.3f  min=%.3f  max=%.3f\n",
                t, dlabel, mean(vals), sd(vals), min(vals), max(vals)))
  }
}


# ── 3. Category frequencies ────────────────────────────────────────────────
cat("\n── 4. Category proportions by drug and block ──\n")
for (t in 1:K) {
  cat_col <- paste0("PDC_cat_", t)
  tab <- table(Drug = ifelse(dat$A0 == 1, "Drug A", "Drug B"),
               Category = dat[[cat_col]], useNA = "ifany")
  cat(sprintf("\n  Block %d:\n", t))
  print(round(prop.table(tab, margin = 1), 3))
}


# ── 4. Cumulative failure by drug × adherence stratum ──────────────────────
cat("\n── 5. Cumulative failure by drug × mean PDC stratum ──\n")
# Compute mean PDC across blocks
pdc_cols <- paste0("PDC_", 1:K)
dat$mean_pdc <- rowMeans(dat[, pdc_cols], na.rm = TRUE)
dat$pdc_strat <- cut(dat$mean_pdc,
                     breaks = c(0, p$pdc_cut_low, p$pdc_cut_high, 1),
                     labels = c("low", "mid", "high"),
                     include.lowest = TRUE)

for (d in 0:1) {
  dlabel <- ifelse(d == 0, "Drug B", "Drug A")
  idx <- dat$A0 == d
  tab <- tapply(dat$Y_final[idx], dat$pdc_strat[idx], mean, na.rm = TRUE)
  n_tab <- tapply(dat$Y_final[idx], dat$pdc_strat[idx], length)
  cat(sprintf("\n  %s:\n", dlabel))
  for (s in names(tab)) {
    cat(sprintf("    %s: risk=%.3f (n=%d)\n", s, tab[s], n_tab[s]))
  }
}


# ── 5. Check qualitative patterns ──────────────────────────────────────────
cat("\n── 6. Qualitative pattern check ──\n")

# (a) Adherence differs by drug
mean_pdc_A <- mean(dat$mean_pdc[dat$A0 == 1], na.rm = TRUE)
mean_pdc_B <- mean(dat$mean_pdc[dat$A0 == 0], na.rm = TRUE)
cat(sprintf("  (a) Mean PDC: Drug A=%.3f, Drug B=%.3f  → %s\n",
            mean_pdc_A, mean_pdc_B,
            ifelse(abs(mean_pdc_A - mean_pdc_B) > 0.02,
                   "PASS (drugs differ)", "CHECK (drugs too similar)")))

# (b) Poor health predicts lower future adherence
# Correlate L_1 with PDC_2
if (K >= 2) {
  r <- cor(dat$L_1, dat$PDC_2, use = "complete.obs")
  cat(sprintf("  (b) Corr(L_1, PDC_2) = %.3f  → %s\n", r,
              ifelse(r < -0.05, "PASS (sicker → lower PDC)", "CHECK")))
}

# (c) Poor adherence worsens later health
if (K >= 2) {
  r <- cor(dat$PDC_1, dat$L_2, use = "complete.obs")
  cat(sprintf("  (c) Corr(PDC_1, L_2) = %.3f  → %s\n", r,
              ifelse(r < -0.05, "PASS (low PDC → worse health)", "CHECK")))
}

# (d) Lower adherence increases failure
low_risk  <- mean(dat$Y_final[dat$pdc_strat == "low"], na.rm = TRUE)
high_risk <- mean(dat$Y_final[dat$pdc_strat == "high"], na.rm = TRUE)
cat(sprintf("  (d) Failure: low_adh=%.3f, high_adh=%.3f  → %s\n",
            low_risk, high_risk,
            ifelse(low_risk > high_risk + 0.02,
                   "PASS (low PDC → more failure)", "CHECK")))

# (e) Drug difference small at high PDC, larger at low PDC
high_A <- mean(dat$Y_final[dat$A0 == 1 & dat$pdc_strat == "high"], na.rm = TRUE)
high_B <- mean(dat$Y_final[dat$A0 == 0 & dat$pdc_strat == "high"], na.rm = TRUE)
low_A  <- mean(dat$Y_final[dat$A0 == 1 & dat$pdc_strat == "low"], na.rm = TRUE)
low_B  <- mean(dat$Y_final[dat$A0 == 0 & dat$pdc_strat == "low"], na.rm = TRUE)
diff_high <- abs(high_A - high_B)
diff_low  <- abs(low_A - low_B)
cat(sprintf("  (e) Drug diff at high PDC: |%.3f - %.3f| = %.3f\n",
            high_A, high_B, diff_high))
cat(sprintf("      Drug diff at low PDC:  |%.3f - %.3f| = %.3f\n",
            low_A, low_B, diff_low))
cat(sprintf("      → %s\n",
            ifelse(diff_low > diff_high + 0.01,
                   "PASS (larger drug diff at low PDC)",
                   "CHECK (drug diff pattern not as expected)")))


# ── 6. Support diagnostics ─────────────────────────────────────────────────
cat("\n── 7. Support diagnostics ──\n")

# Proportion following static regimes for all blocks
for (d in 0:1) {
  idx <- dat$A0 == d
  n_d <- sum(idx)
  dlabel <- ifelse(d == 0, "Drug B", "Drug A")

  for (target in c("high", "mid", "low")) {
    follows <- rep(TRUE, n_d)
    for (t in 1:K) {
      cat_col <- paste0("PDC_cat_", t)
      vals <- as.character(dat[[cat_col]][idx])
      follows <- follows & !is.na(vals) & vals == target
    }
    pct <- 100 * sum(follows) / n_d
    warning_flag <- ifelse(pct < 1, " *** SPARSE ***", "")
    cat(sprintf("  %s, always_%s: %d / %d (%.1f%%)%s\n",
                dlabel, target, sum(follows), n_d, pct, warning_flag))
  }
}

# Fraction of PDC near boundaries
cat("\n  Fraction of PDC values near boundaries:\n")
for (t in 1:K) {
  pdc_col <- paste0("PDC_", t)
  vals <- dat[[pdc_col]][!is.na(dat[[pdc_col]])]
  cat(sprintf("    Block %d: PDC<0.10: %.2f%%,  PDC>0.95: %.2f%%\n",
              t, 100 * mean(vals < 0.10), 100 * mean(vals > 0.95)))
}


# ── 7. Block-specific failure rates ────────────────────────────────────────
cat("\n── 8. Block-specific new failure rates ──\n")
for (t in 1:K) {
  y_col <- paste0("Y_", t)
  if (t == 1) {
    new_fail <- dat[[y_col]]
  } else {
    prev_col <- paste0("Y_", t - 1)
    new_fail <- dat[[y_col]] == 1 & dat[[prev_col]] == 0
  }
  for (d in 0:1) {
    idx <- dat$A0 == d
    at_risk <- if (t == 1) idx else idx & dat[[prev_col]] == 0
    n_risk <- sum(at_risk, na.rm = TRUE)
    n_new  <- sum(new_fail[at_risk], na.rm = TRUE)
    dlabel <- ifelse(d == 0, "Drug B", "Drug A")
    cat(sprintf("  Block %d, %s: %d new failures / %d at risk (%.1f%%)\n",
                t, dlabel, n_new, n_risk,
                ifelse(n_risk > 0, 100 * n_new / n_risk, 0)))
  }
}


# ── 8. Plots (if interactive) ──────────────────────────────────────────────
if (interactive()) {
  cat("\n── 9. Generating plots ──\n")

  # PDC density by drug and block
  old_par <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  for (t in 1:K) {
    pdc_col <- paste0("PDC_", t)
    pdc_A <- dat[[pdc_col]][dat$A0 == 1 & !is.na(dat[[pdc_col]])]
    pdc_B <- dat[[pdc_col]][dat$A0 == 0 & !is.na(dat[[pdc_col]])]

    plot(density(pdc_B, from = 0, to = 1), col = "blue", lwd = 2,
         main = sprintf("Block %d: PDC by drug", t),
         xlab = "PDC", ylim = c(0, 5))
    lines(density(pdc_A, from = 0, to = 1), col = "red", lwd = 2)
    abline(v = c(p$pdc_cut_low, p$pdc_cut_high), lty = 2, col = "grey40")
    legend("topleft", c("Drug B", "Drug A"), col = c("blue", "red"),
           lwd = 2, cex = 0.7)
  }
  par(old_par)

  cat("  Plots generated in graphics device.\n")
}


# ── 9. Quick oracle check ──────────────────────────────────────────────────
cat("\n── 10. Quick oracle truth estimates ──\n")
cat("  (Using n_oracle = 20000 for speed)\n\n")

for (d in 0:1) {
  dlabel <- ifelse(d == 0, "Drug B", "Drug A")
  r_nat <- oracle_natural_v4(drug = d, K = K, n_oracle = 20000, seed = 7777)
  cat(sprintf("  %s natural risk: %.3f\n", dlabel, r_nat))
}

for (d in 0:1) {
  dlabel <- ifelse(d == 0, "Drug B", "Drug A")
  for (cat in c("high", "mid", "low")) {
    r <- oracle_static_category_v4(
      category = cat, drug = d, K = K, n_oracle = 20000,
      seed = 8000 + d * 10 + match(cat, c("high", "mid", "low"))
    )
    cat(sprintf("  %s, always_%s: risk = %.3f\n", dlabel, cat, r))
  }
}

for (d in 0:1) {
  dlabel <- ifelse(d == 0, "Drug B", "Drug A")
  r <- oracle_shift_v4(delta = 0.10, drug = d, K = K,
                       n_oracle = 20000, seed = 9000 + d)
  cat(sprintf("  %s, +10%% shift: risk = %.3f\n", dlabel, r))
}

cat("\nCalibration complete.\n")
