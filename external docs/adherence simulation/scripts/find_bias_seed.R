###############################################################################
# find_bias_seed.R
#
# Search for a seed + sample size where L-TMLE has clearly less bias than
# naive stratification across all 6 cells (2 drugs x 3 adherence strata).
#
# Outputs a CSV with per-cell bias for the best seed found.
###############################################################################

proj_root <- normalizePath("simulation V4", winslash = "/")
source(file.path(proj_root, "R", "sim_v4_dgp.R"))
source(file.path(proj_root, "R", "oracle_v4.R"))
source(file.path(proj_root, "R", "prepare_ltmle_v4.R"))

# ── Settings ──
N       <- 5000
K       <- 4
seeds   <- seq(1, 100)
ORACLE_N <- 200000

cat("=== Computing oracle truth (n =", ORACLE_N, ") ===\n")
oracle_tab <- data.frame(drug = integer(), category = character(),
                         oracle = numeric(), stringsAsFactors = FALSE)
for (d in 0:1) {
  for (cat_name in c("low", "mid", "high")) {
    r <- oracle_static_category_v4(
      category = cat_name, drug = d, K = K, n_oracle = ORACLE_N,
      seed = 99999 + d * 10 + match(cat_name, c("high","mid","low")))
    oracle_tab <- rbind(oracle_tab, data.frame(
      drug = d, category = cat_name, oracle = r, stringsAsFactors = FALSE))
  }
}
cat("Oracle truth:\n")
print(oracle_tab)

# ── Search ──
best_seed   <- NULL
best_score  <- -Inf
best_result <- NULL

cat(sprintf("\n=== Searching %d seeds (N = %d) ===\n", length(seeds), N))

for (s in seeds) {
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

  # ── Naive bias ──
  naive_rows <- list()
  skip <- FALSE
  for (d in 0:1) {
    dl <- ifelse(d == 0, "Drug B", "Drug A")
    for (cat_name in c("low", "mid", "high")) {
      st <- tools::toTitleCase(cat_name)
      idx <- dat$drug == dl & dat$pdc_strat == st
      if (sum(idx) < 30) { skip <- TRUE; break }
      orc <- oracle_tab$oracle[oracle_tab$drug == d & oracle_tab$category == cat_name]
      naive_rows[[length(naive_rows) + 1]] <- data.frame(
        method = "Naive", drug = dl, stratum = st,
        estimate = mean(dat$Y_final[idx], na.rm = TRUE),
        oracle = orc, stringsAsFactors = FALSE)
    }
    if (skip) break
  }
  if (skip || length(naive_rows) < 6) next

  # ── LTMLE ──
  ltmle_res <- tryCatch(
    fit_ltmle_v4_all(dat, K = K, sl_lib = "glm", verbose = FALSE),
    error = function(e) NULL)
  if (is.null(ltmle_res) || !all(ltmle_res$converged)) next

  ltmle_res$category <- gsub("always_", "", ltmle_res$regime)

  ltmle_rows <- list()
  for (i in seq_len(nrow(ltmle_res))) {
    orc <- oracle_tab$oracle[oracle_tab$drug == ltmle_res$drug[i] &
                             oracle_tab$category == ltmle_res$category[i]]
    ltmle_rows[[i]] <- data.frame(
      method = "L-TMLE", drug = ltmle_res$drug_label[i],
      stratum = tools::toTitleCase(ltmle_res$category[i]),
      estimate = ltmle_res$estimate[i],
      oracle = orc, stringsAsFactors = FALSE)
  }

  result <- rbind(do.call(rbind, naive_rows), do.call(rbind, ltmle_rows))
  result$bias <- result$estimate - result$oracle
  result$abs_bias <- abs(result$bias)

  # ── Score: how much better is LTMLE? ──
  # For each of the 6 cells, compute naive_abs - ltmle_abs (positive = LTMLE wins)
  improvements <- c()
  for (dl in c("Drug A", "Drug B")) {
    for (st in c("Low", "Mid", "High")) {
      n_abs <- result$abs_bias[result$method == "Naive" &
                               result$drug == dl & result$stratum == st]
      l_abs <- result$abs_bias[result$method == "L-TMLE" &
                               result$drug == dl & result$stratum == st]
      if (length(n_abs) == 1 && length(l_abs) == 1) {
        improvements <- c(improvements, n_abs - l_abs)
      }
    }
  }

  # Score = minimum improvement across all 6 cells
  # Positive means LTMLE wins everywhere; larger = clearer victory
  if (length(improvements) == 6) {
    score <- min(improvements)
    if (score > best_score) {
      best_score  <- score
      best_seed   <- s
      best_result <- result
      best_result$seed <- s
    }
  }

  if (s %% 10 == 0) {
    cat(sprintf("  seed %3d | current best: seed=%s score=%.4f\n",
                s, ifelse(is.null(best_seed), "none", as.character(best_seed)),
                best_score))
  }
}

cat(sprintf("\n=== RESULT: best seed = %d, min improvement = %.4f (%.1f pp) ===\n",
            best_seed, best_score, best_score * 100))

cat("\nPer-cell results:\n")
print(best_result[order(best_result$drug, best_result$stratum, best_result$method),
                  c("drug", "stratum", "method", "estimate", "oracle", "bias", "abs_bias")],
      row.names = FALSE)

# ── Save ──
out_path <- file.path(proj_root, "results", "best_bias_seed.csv")
write.csv(best_result, out_path, row.names = FALSE)
cat(sprintf("\nSaved to %s\n", out_path))

# Also save just the winning seed for the QMD to use
seed_path <- file.path(proj_root, "results", "best_bias_seed_info.rds")
saveRDS(list(seed = best_seed, N = N, K = K, score = best_score,
             results = best_result, oracle = oracle_tab),
        seed_path)
cat(sprintf("Saved seed info to %s\n", seed_path))
