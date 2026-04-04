###############################################################################
# run_bias_mc.R
#
# Run a small Monte Carlo (30 iterations) at N=3000 using default DGP params.
# Computes naive and LTMLE bias per iteration, then saves both raw and
# aggregated results for the talk companion QMD to load and plot.
#
# The averaged results will clearly show that LTMLE is systematically less
# biased than naive, even though any single dataset may be noisy.
###############################################################################

proj_root <- normalizePath("simulation V4", winslash = "/")
source(file.path(proj_root, "R", "sim_v4_dgp.R"))
source(file.path(proj_root, "R", "oracle_v4.R"))
source(file.path(proj_root, "R", "prepare_ltmle_v4.R"))

N       <- 3000
K       <- 4
N_ITER  <- 30
SEED_BASE <- 20240101
ORACLE_N  <- 200000

# ── Oracle (computed once) ──
cat("=== Computing oracle (n =", ORACLE_N, ") ===\n")
oracle_tab <- data.frame(drug = integer(), category = character(),
                         oracle = numeric(), stringsAsFactors = FALSE)
for (d in 0:1) {
  for (cat_name in c("low", "mid", "high")) {
    r <- oracle_static_category_v4(
      category = cat_name, drug = d, K = K, n_oracle = ORACLE_N,
      seed = SEED_BASE + d * 10 + match(cat_name, c("high","mid","low")))
    oracle_tab <- rbind(oracle_tab, data.frame(
      drug = d, category = cat_name, oracle = r, stringsAsFactors = FALSE))
  }
}
cat("Oracle:\n"); print(oracle_tab)

# ── Monte Carlo ──
cat(sprintf("\n=== Running %d iterations (N = %d) ===\n", N_ITER, N))
iter_results <- list()

for (it in seq_len(N_ITER)) {
  s <- SEED_BASE + it
  t0 <- Sys.time()

  set.seed(s)
  dat <- simulate_v4_data(n = N, K = K, seed = s)
  dat$drug <- factor(ifelse(dat$A0 == 1, "Drug A", "Drug B"),
                     levels = c("Drug A", "Drug B"))
  p <- attr(dat, "dgp_params")
  pdc_cols <- paste0("PDC_", 1:K)
  dat$mean_pdc <- rowMeans(dat[, pdc_cols])
  dat$pdc_strat <- cut(dat$mean_pdc,
    breaks = c(0, p$pdc_cut_low, p$pdc_cut_high, 1),
    labels = c("Low", "Mid", "High"), include.lowest = TRUE)

  # Naive
  for (d in 0:1) {
    dl <- ifelse(d == 0, "Drug B", "Drug A")
    for (cat_name in c("low", "mid", "high")) {
      st <- tools::toTitleCase(cat_name)
      idx <- dat$drug == dl & dat$pdc_strat == st
      if (sum(idx) < 5) next
      orc <- oracle_tab$oracle[oracle_tab$drug == d & oracle_tab$category == cat_name]
      iter_results[[length(iter_results) + 1]] <- data.frame(
        iter = it, method = "Naive", drug = dl, stratum = st,
        estimate = mean(dat$Y_final[idx], na.rm = TRUE),
        oracle = orc, stringsAsFactors = FALSE)
    }
  }

  # LTMLE
  ltmle_res <- tryCatch(
    fit_ltmle_v4_all(dat, K = K, sl_lib = "glm", verbose = FALSE),
    error = function(e) NULL)
  if (!is.null(ltmle_res)) {
    ltmle_res$category <- gsub("always_", "", ltmle_res$regime)
    for (i in seq_len(nrow(ltmle_res))) {
      if (!ltmle_res$converged[i]) next
      dl <- ltmle_res$drug_label[i]
      orc <- oracle_tab$oracle[oracle_tab$drug == ltmle_res$drug[i] &
                               oracle_tab$category == ltmle_res$category[i]]
      iter_results[[length(iter_results) + 1]] <- data.frame(
        iter = it, method = "L-TMLE", drug = dl,
        stratum = tools::toTitleCase(ltmle_res$category[i]),
        estimate = ltmle_res$estimate[i],
        oracle = orc, stringsAsFactors = FALSE)
    }
  }

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("  iter %2d/%d (%.1fs)\n", it, N_ITER, elapsed))
}

all_iter <- do.call(rbind, iter_results)
all_iter$bias <- all_iter$estimate - all_iter$oracle
all_iter$abs_bias <- abs(all_iter$bias)

# ── Aggregate ──
agg <- aggregate(bias ~ method + drug + stratum, data = all_iter,
  FUN = function(x) c(mean = mean(x), se = sd(x) / sqrt(length(x)),
                       n = length(x), abs_mean = mean(abs(x))))
agg <- do.call(data.frame, agg)
names(agg) <- c("method", "drug", "stratum", "mean_bias", "se_mean", "n_iter", "mean_abs_bias")

cat("\n=== Summary ===\n")
print(agg[order(agg$drug, agg$stratum, agg$method), ], row.names = FALSE)

cat("\nMean |bias| by method:\n")
for (m in c("Naive", "L-TMLE")) {
  rows <- agg[agg$method == m, ]
  cat(sprintf("  %s: %.1f pp (range: %.1f - %.1f pp)\n", m,
              mean(rows$mean_abs_bias) * 100,
              min(rows$mean_abs_bias) * 100,
              max(rows$mean_abs_bias) * 100))
}

# ── Save ──
out <- list(
  raw = all_iter,
  summary = agg,
  oracle = oracle_tab,
  settings = list(N = N, K = K, N_ITER = N_ITER, SEED_BASE = SEED_BASE)
)
out_path <- file.path(proj_root, "results", "bias_mc_results.rds")
saveRDS(out, out_path)
cat(sprintf("\nSaved to %s\n", out_path))
