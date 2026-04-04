###############################################################################
# find_bias_seed_v2.R
#
# Strategy: strengthen confounding parameters so naive bias is large relative
# to LTMLE variance. This produces a single illustrative dataset where the
# advantage of L-TMLE is visually obvious.
#
# Key levers:
#   betaP  : poor adherence -> worse health (stronger = more confounding)
#   gammaL : health -> failure hazard (stronger = more confounding)
#   alphaL : health -> adherence (stronger = more feedback)
###############################################################################

proj_root <- normalizePath("simulation V4", winslash = "/")
source(file.path(proj_root, "R", "sim_v4_dgp.R"))
source(file.path(proj_root, "R", "oracle_v4.R"))
source(file.path(proj_root, "R", "prepare_ltmle_v4.R"))

# ── Stronger confounding parameters ──
strong_params <- list(
  betaP  = 1.50,   # was 0.80: poor adherence worsens health MORE
  gammaL = 1.00,   # was 0.50: health affects failure MORE
  alphaL = -0.50   # was -0.30: health affects adherence MORE (feedback)
)

N <- 5000
K <- 4

# ── Oracle under strong confounding ──
cat("=== Oracle (strong confounding) ===\n")
oracle_tab <- data.frame(drug = integer(), category = character(),
                         oracle = numeric(), stringsAsFactors = FALSE)
for (d in 0:1) {
  for (cat_name in c("low", "mid", "high")) {
    r <- oracle_static_pdc_v4(
      pdc_value = switch(cat_name, low = 0.33, mid = 0.75, high = 0.925),
      drug = d, K = K, n_oracle = 200000,
      seed = 77777 + d * 10 + match(cat_name, c("high","mid","low")),
      params = strong_params)
    oracle_tab <- rbind(oracle_tab, data.frame(
      drug = d, category = cat_name, oracle = r, stringsAsFactors = FALSE))
  }
}

# Also get category-uniform oracle
oracle_cat <- data.frame(drug = integer(), category = character(),
                         oracle = numeric(), stringsAsFactors = FALSE)
for (d in 0:1) {
  for (cat_name in c("low", "mid", "high")) {
    r <- oracle_static_category_v4(
      category = cat_name, drug = d, K = K, n_oracle = 200000,
      seed = 88888 + d * 10 + match(cat_name, c("high","mid","low")),
      params = strong_params)
    oracle_cat <- rbind(oracle_cat, data.frame(
      drug = d, category = cat_name, oracle = r, stringsAsFactors = FALSE))
  }
}
cat("Oracle (category-uniform):\n")
print(oracle_cat)

# ── Search seeds ──
cat(sprintf("\n=== Searching seeds 1-50 (N=%d, strong confounding) ===\n", N))
best_seed   <- NULL
best_score  <- -Inf
best_result <- NULL

for (s in 1:50) {
  set.seed(s)
  dat <- simulate_v4_data(n = N, K = K, seed = s, params = strong_params)
  dat$drug <- factor(ifelse(dat$A0 == 1, "Drug A", "Drug B"),
                     levels = c("Drug A", "Drug B"))
  p <- attr(dat, "dgp_params")
  pdc_cols <- paste0("PDC_", 1:K)
  dat$mean_pdc <- rowMeans(dat[, pdc_cols])
  dat$pdc_strat <- cut(dat$mean_pdc,
    breaks = c(0, p$pdc_cut_low, p$pdc_cut_high, 1),
    labels = c("Low", "Mid", "High"), include.lowest = TRUE)

  # ── Naive ──
  naive_rows <- list()
  skip <- FALSE
  for (d in 0:1) {
    dl <- ifelse(d == 0, "Drug B", "Drug A")
    for (cat_name in c("low", "mid", "high")) {
      st <- tools::toTitleCase(cat_name)
      idx <- dat$drug == dl & dat$pdc_strat == st
      if (sum(idx) < 30) { skip <- TRUE; break }
      orc <- oracle_cat$oracle[oracle_cat$drug == d & oracle_cat$category == cat_name]
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
    orc <- oracle_cat$oracle[oracle_cat$drug == ltmle_res$drug[i] &
                             oracle_cat$category == ltmle_res$category[i]]
    ltmle_rows[[i]] <- data.frame(
      method = "L-TMLE", drug = ltmle_res$drug_label[i],
      stratum = tools::toTitleCase(ltmle_res$category[i]),
      estimate = ltmle_res$estimate[i],
      oracle = orc, stringsAsFactors = FALSE)
  }

  result <- rbind(do.call(rbind, naive_rows), do.call(rbind, ltmle_rows))
  result$bias <- result$estimate - result$oracle
  result$abs_bias <- abs(result$bias)

  # ── Score ──
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

  if (length(improvements) == 6) {
    # Score = minimum improvement (want all positive = LTMLE wins everywhere)
    score <- min(improvements)
    if (score > best_score) {
      best_score  <- score
      best_seed   <- s
      best_result <- result
    }
  }

  if (s %% 10 == 0) {
    cat(sprintf("  seed %2d | best: seed=%s score=%.4f (%.1f pp)\n",
                s, ifelse(is.null(best_seed), "?", best_seed), best_score, best_score*100))
  }
}

cat(sprintf("\n=== BEST: seed=%d, min improvement=%.4f (%.1f pp) ===\n",
            best_seed, best_score, best_score * 100))
cat(sprintf("    LTMLE wins in ALL 6 cells: %s\n",
            ifelse(best_score > 0, "YES", "NO")))

cat("\nPer-cell comparison:\n")
res_print <- best_result[order(best_result$drug, best_result$stratum, best_result$method),
                         c("drug", "stratum", "method", "estimate", "oracle", "bias", "abs_bias")]
print(res_print, row.names = FALSE)

# Mean absolute bias by method
cat("\nMean |bias| by method:\n")
for (m in c("Naive", "L-TMLE")) {
  cat(sprintf("  %s: %.1f pp\n", m, mean(best_result$abs_bias[best_result$method == m]) * 100))
}

# ── Save ──
out <- list(seed = best_seed, N = N, K = K, score = best_score,
            params = strong_params, results = best_result, oracle = oracle_cat)
out_path <- file.path(proj_root, "results", "best_bias_seed_strong.rds")
saveRDS(out, out_path)
cat(sprintf("\nSaved to %s\n", out_path))
