###############################################################################
# 03_mc_phase1_validation.R — Quick validation MC (Phase 1)
#
# Purpose: validate 1-block and 4-block LTMLE + LMTP after code fixes
# Default: n=1500, B=100
# Estimands: Drug A/B always_low, always_high, +10% shift
# Output: CSV summary + bias figure
###############################################################################

cat("=== MC Phase 1: Quick Validation ===\n")

source(file.path("R", "sim_v4_dgp.R"))
source(file.path("R", "oracle_v4.R"))
source(file.path("R", "prepare_ltmle_v4.R"))
source(file.path("R", "prepare_lmtp_v4.R"))

for (pkg in c("ltmle", "lmtp"))
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(sprintf("Package '%s' required", pkg))

# ── Config ──────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
B      <- if (length(args) >= 1) as.integer(args[1]) else 100L
n_sim  <- 1500L
K      <- 4L
seed0  <- 20260401L

out_dir <- file.path("results", "mc_phase1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Oracle (precompute once) ────────────────────────────────────────────────
cat("Computing oracle...\n")
oracle <- list()
for (d in 0:1) {
  dl <- ifelse(d == 0, "Drug_B", "Drug_A")
  for (cn in c("low", "high")) {
    oracle[[paste0("cat_", cn, "_", dl)]] <-
      oracle_static_category_v4(cn, drug = d, K = K, n_oracle = 100000, seed = seed0 + d * 10 + match(cn, c("low","high")))
  }
  oracle[[paste0("shift_", dl)]] <-
    oracle_shift_v4(delta = 0.10, drug = d, K = K, n_oracle = 100000, seed = seed0 + 100 + d)
}
cat("Oracle done.\n")

# ── MC loop ─────────────────────────────────────────────────────────────────
store <- data.frame(rep = integer(), estimand = character(), method = character(),
                    estimate = numeric(), se = numeric(),
                    converged = logical(), stringsAsFactors = FALSE)

t0 <- Sys.time()
for (b in seq_len(B)) {
  if (b %% 10 == 0 || b == 1) cat(sprintf("  Rep %d/%d\n", b, B))
  dat <- simulate_v4_data(n = n_sim, K = K, seed = seed0 + b * 13)

  # LTMLE
  for (d in 0:1) {
    dl <- ifelse(d == 0, "Drug_B", "Drug_A")
    for (reg in c("always_low", "always_high")) {
      cn <- sub("always_", "", reg)
      res <- fit_ltmle_v4_static(dat, drug = d, regime = reg, K = K,
                                  sl_lib = "glm", verbose = FALSE)
      store <- rbind(store, data.frame(
        rep = b, estimand = paste0("cat_", cn, "_", dl), method = "ltmle",
        estimate = res$estimate, se = res$se,
        converged = res$converged, stringsAsFactors = FALSE))
    }
  }

  # LMTP shift
  for (d in 0:1) {
    dl <- ifelse(d == 0, "Drug_B", "Drug_A")
    res <- fit_lmtp_shift_v4(dat, delta = 0.10, drug = d, K = K, folds = 1)
    store <- rbind(store, data.frame(
      rep = b, estimand = paste0("shift_", dl), method = "lmtp_shift",
      estimate = res$estimate, se = res$se,
      converged = res$converged, stringsAsFactors = FALSE))
  }

  if (b %% 25 == 0) saveRDS(store, file.path(out_dir, "phase1_raw.rds"))
}
elapsed <- difftime(Sys.time(), t0, units = "mins")
cat(sprintf("MC loop done in %.1f min\n", as.numeric(elapsed)))

# ── Summary ─────────────────────────────────────────────────────────────────
summary_rows <- list()
for (est_name in unique(store$estimand)) {
  for (meth in unique(store$method[store$estimand == est_name])) {
    sub <- store[store$estimand == est_name & store$method == meth & store$converged, ]
    if (nrow(sub) == 0) next
    truth <- oracle[[est_name]]
    if (is.null(truth)) truth <- NA_real_
    mn <- mean(sub$estimate, na.rm = TRUE)
    bias <- mn - truth
    emp_se <- sd(sub$estimate, na.rm = TRUE)
    rmse <- sqrt(bias^2 + emp_se^2)
    mn_se <- mean(sub$se, na.rm = TRUE)
    cov <- if (!is.na(truth)) mean(sub$estimate - 1.96*sub$se <= truth &
                                    sub$estimate + 1.96*sub$se >= truth) else NA
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      estimand = est_name, method = meth, truth = truth,
      mean_est = mn, bias = bias, emp_se = emp_se, rmse = rmse,
      mean_se = mn_se, coverage = cov,
      n_converged = nrow(sub), n_total = B, stringsAsFactors = FALSE)
  }
}
mc_summary <- do.call(rbind, summary_rows)
write.csv(mc_summary, file.path(out_dir, "phase1_summary.csv"), row.names = FALSE)
saveRDS(store, file.path(out_dir, "phase1_raw.rds"))

cat("\n=== Phase 1 Summary ===\n")
print(mc_summary[, c("estimand","method","truth","mean_est","bias","emp_se","coverage","n_converged")])
cat(sprintf("\nSaved to %s/\n", out_dir))
