###############################################################################
# 07_mc_bias_study.R — Monte Carlo bias comparison for symposium
#
# Three estimator classes:
#   1. Static categorical L-TMLE (3-category: low/mid/high)
#   2. Shift LMTP (+10% PDC)
#   3. Naive adherence-stratified conditioning
#
# Three oracle variants for categorical L-TMLE:
#   A. Uniform within category
#   B. Empirical within category (PRIMARY)
#   C. Grid/midpoint
#
# Usage:
#   Rscript scripts/07_mc_bias_study.R          # pilot: B=50
#   Rscript scripts/07_mc_bias_study.R 200      # full run
###############################################################################

cat("=== V4 MC Bias Study ===\n\n")

source(file.path("R", "sim_v4_dgp.R"))
source(file.path("R", "oracle_v4.R"))
source(file.path("R", "prepare_ltmle_v4.R"))
source(file.path("R", "prepare_lmtp_v4.R"))

for (pkg in c("ltmle", "lmtp"))
  if (!requireNamespace(pkg, quietly = TRUE)) stop(paste("Need", pkg))

# ── Config ──────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
B      <- if (length(args) >= 1) as.integer(args[1]) else 50L
n_sim  <- 2000L
K      <- 4L
seed0  <- 20260701L
delta  <- 0.10

out_dir <- file.path("results", "mc_bias_study")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Config: B=%d, n=%d, K=%d, delta=%.2f\n\n", B, n_sim, K, delta))

# ── Step 1: Precompute oracles (once, large n) ─────────────────────────────
cat("── Computing oracles ──\n")

# Generate a large reference dataset for empirical oracle
ref_data <- simulate_v4_data(n = 20000, K = K, seed = seed0 - 1)
p <- default_params_v4()

oracle <- list()
for (d in 0:1) {
  dl <- ifelse(d == 0, "Drug_B", "Drug_A")
  for (cn in c("low", "mid", "high")) {
    s_base <- seed0 + d * 100 + match(cn, c("low","mid","high"))

    oracle[[paste0(cn, "_", dl, "_uniform")]] <-
      oracle_static_category_v4(cn, drug = d, K = K, n_oracle = 200000,
                                seed = s_base)
    oracle[[paste0(cn, "_", dl, "_empirical")]] <-
      oracle_static_category_empirical_v4(cn, drug = d, K = K, n_oracle = 200000,
                                          seed = s_base + 10, ref_data = ref_data)
    oracle[[paste0(cn, "_", dl, "_midpoint")]] <-
      oracle_static_category_midpoint_v4(cn, drug = d, K = K, n_oracle = 200000,
                                         seed = s_base + 20)

    cat(sprintf("  %s %s: uniform=%.4f empirical=%.4f midpoint=%.4f\n",
                dl, cn,
                oracle[[paste0(cn, "_", dl, "_uniform")]],
                oracle[[paste0(cn, "_", dl, "_empirical")]],
                oracle[[paste0(cn, "_", dl, "_midpoint")]]))
  }

  # Shift oracle
  oracle[[paste0("shift_", dl)]] <-
    oracle_shift_v4(delta, drug = d, K = K, n_oracle = 200000,
                    seed = seed0 + 500 + d)
  cat(sprintf("  %s shift+%.0f%%: %.4f\n", dl, delta*100, oracle[[paste0("shift_", dl)]]))

  # Natural oracle (for naive bias computation)
  oracle[[paste0("natural_", dl)]] <-
    oracle_natural_v4(drug = d, K = K, n_oracle = 200000, seed = seed0 + 600 + d)
}

cat("\n── Oracle alignment table ──\n")
align_rows <- list()
for (d in 0:1) {
  dl <- ifelse(d == 0, "Drug_B", "Drug_A")
  for (cn in c("low", "mid", "high")) {
    align_rows[[length(align_rows) + 1]] <- data.frame(
      drug = dl, category = cn,
      oracle_uniform  = oracle[[paste0(cn, "_", dl, "_uniform")]],
      oracle_empirical = oracle[[paste0(cn, "_", dl, "_empirical")]],
      oracle_midpoint = oracle[[paste0(cn, "_", dl, "_midpoint")]],
      stringsAsFactors = FALSE)
  }
}
oracle_align <- do.call(rbind, align_rows)
oracle_align$diff_emp_unif <- oracle_align$oracle_empirical - oracle_align$oracle_uniform
print(oracle_align)
write.csv(oracle_align, file.path(out_dir, "oracle_alignment.csv"), row.names = FALSE)

cat("\n")

# ── Step 2: MC loop ────────────────────────────────────────────────────────
cat("── Starting MC loop ──\n")

store <- data.frame(
  rep = integer(), drug = character(), category = character(),
  method = character(), estimate = numeric(), se = numeric(),
  converged = logical(), support_pct = numeric(),
  stringsAsFactors = FALSE
)

t0 <- Sys.time()
for (b in seq_len(B)) {
  if (b %% 10 == 0 || b == 1) cat(sprintf("  Rep %d/%d\n", b, B))
  dat <- simulate_v4_data(n = n_sim, K = K, seed = seed0 + b * 17)

  # Precompute mean PDC strata for naive
  pdc_cols <- paste0("PDC_", 1:K)
  dat$mean_pdc <- rowMeans(dat[, pdc_cols])
  dat$strat <- cut(dat$mean_pdc, c(0, p$pdc_cut_low, p$pdc_cut_high, 1),
                   labels = c("low", "mid", "high"), include.lowest = TRUE)

  for (d in 0:1) {
    dl <- ifelse(d == 0, "Drug_B", "Drug_A")
    idx <- dat$A0 == d
    nd <- sum(idx)

    # ── Naive: condition on mean PDC stratum ──
    for (cn in c("low", "mid", "high")) {
      nv <- mean(dat$Y_final[idx & dat$strat == cn], na.rm = TRUE)
      n_in <- sum(idx & dat$strat == cn, na.rm = TRUE)
      store <- rbind(store, data.frame(
        rep = b, drug = dl, category = cn, method = "naive",
        estimate = nv, se = NA_real_, converged = TRUE,
        support_pct = 100 * n_in / nd, stringsAsFactors = FALSE))
    }

    # ── LTMLE: static category regimes ──
    for (cn in c("low", "mid", "high")) {
      regime <- paste0("always_", cn)
      res <- fit_ltmle_v4_static(dat, drug = d, regime = regime,
                                  K = K, sl_lib = "glm", verbose = FALSE)
      # Support: fraction following sustained regime
      follows <- rep(TRUE, nd)
      for (t in seq_len(K)) {
        v <- as.character(dat[[paste0("PDC_cat_", t)]][idx])
        follows <- follows & !is.na(v) & v == cn
      }
      supp <- 100 * sum(follows) / nd

      store <- rbind(store, data.frame(
        rep = b, drug = dl, category = cn, method = "ltmle",
        estimate = res$estimate, se = res$se,
        converged = res$converged, support_pct = supp,
        stringsAsFactors = FALSE))
    }

    # ── LMTP: +10% shift ──
    res_s <- fit_lmtp_shift_v4(dat, delta = delta, drug = d, K = K, folds = 1)
    store <- rbind(store, data.frame(
      rep = b, drug = dl, category = "shift", method = "lmtp_shift",
      estimate = res_s$estimate, se = res_s$se,
      converged = res_s$converged, support_pct = NA_real_,
      stringsAsFactors = FALSE))
  }

  # Checkpoint
  if (b %% 25 == 0 || b == B) {
    saveRDS(store, file.path(out_dir, "mc_raw.rds"))
  }
}
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("MC loop done in %.1f min\n\n", elapsed))

# ── Step 3: Summarize ──────────────────────────────────────────────────────
cat("── Summarizing ──\n")

# Map each (drug, category, method) to its oracle
get_oracle <- function(drug, category, method) {
  if (method == "lmtp_shift") return(oracle[[paste0("shift_", drug)]])
  # For ltmle, use empirical oracle (PRIMARY)
  if (method == "ltmle") return(oracle[[paste0(category, "_", drug, "_empirical")]])
  # For naive, use empirical oracle (same estimand target for comparison)
  if (method == "naive") return(oracle[[paste0(category, "_", drug, "_empirical")]])
  NA_real_
}

summary_rows <- list()
for (dl in c("Drug_A", "Drug_B")) {
  for (cn in unique(store$category)) {
    for (m in unique(store$method[store$category == cn])) {
      sub <- store[store$drug == dl & store$category == cn &
                     store$method == m & store$converged, ]
      if (nrow(sub) == 0) next

      truth <- get_oracle(dl, cn, m)
      mn <- mean(sub$estimate, na.rm = TRUE)
      bias <- mn - truth
      emp_se <- sd(sub$estimate, na.rm = TRUE)
      rmse <- sqrt(bias^2 + emp_se^2)
      mn_se <- mean(sub$se, na.rm = TRUE)
      cov <- if (!is.na(truth) && all(!is.na(sub$se)))
        mean(sub$estimate - 1.96*sub$se <= truth &
               sub$estimate + 1.96*sub$se >= truth) else NA_real_
      mn_supp <- mean(sub$support_pct, na.rm = TRUE)
      n_conv <- nrow(sub)
      n_fail <- sum(store$drug == dl & store$category == cn &
                      store$method == m & !store$converged)

      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        drug = dl, category = cn, method = m,
        oracle = truth, mean_est = mn, bias = bias,
        emp_se = emp_se, rmse = rmse, mean_se = mn_se,
        coverage = cov, mean_support = mn_supp,
        n_converged = n_conv, n_failed = n_fail,
        stringsAsFactors = FALSE)
    }
  }
}
mc_summary <- do.call(rbind, summary_rows)

# Also add oracle columns for all three variants (for ltmle rows)
mc_summary$oracle_uniform <- NA_real_
mc_summary$oracle_empirical <- NA_real_
mc_summary$oracle_midpoint <- NA_real_
for (i in seq_len(nrow(mc_summary))) {
  r <- mc_summary[i, ]
  if (r$method %in% c("ltmle", "naive") && r$category != "shift") {
    mc_summary$oracle_uniform[i] <- oracle[[paste0(r$category, "_", r$drug, "_uniform")]]
    mc_summary$oracle_empirical[i] <- oracle[[paste0(r$category, "_", r$drug, "_empirical")]]
    mc_summary$oracle_midpoint[i] <- oracle[[paste0(r$category, "_", r$drug, "_midpoint")]]
  }
}

write.csv(mc_summary, file.path(out_dir, "mc_summary.csv"), row.names = FALSE)
saveRDS(store, file.path(out_dir, "mc_raw.rds"))
saveRDS(oracle, file.path(out_dir, "oracle_all.rds"))

# ── Step 4: Print ──────────────────────────────────────────────────────────
cat("\n=== MC BIAS SUMMARY (primary oracle = empirical) ===\n\n")
cat(sprintf("  %-8s %-8s %-12s %7s %7s %+8s %7s %7s %5s %5s\n",
            "Drug","Cat","Method","Oracle","Mean","Bias","EmpSE","RMSE","Cov","N"))
cat(paste(rep("-", 95), collapse = ""), "\n")
for (i in seq_len(nrow(mc_summary))) {
  r <- mc_summary[i, ]
  cat(sprintf("  %-8s %-8s %-12s %7.4f %7.4f %+8.4f %7.4f %7.4f %5s %5d\n",
              r$drug, r$category, r$method,
              r$oracle, r$mean_est, r$bias, r$emp_se, r$rmse,
              ifelse(is.na(r$coverage), "  -- ", sprintf("%.2f", r$coverage)),
              r$n_converged))
}

# Support diagnostics
cat("\n=== SUPPORT DIAGNOSTICS (static regimes) ===\n\n")
supp_diag <- store[store$method == "ltmle" & !is.na(store$support_pct), ]
supp_agg <- aggregate(support_pct ~ drug + category, data = supp_diag,
                      FUN = function(x) c(mean = mean(x), min = min(x),
                                           pct_below_5 = mean(x < 5) * 100))
supp_agg <- do.call(data.frame, supp_agg)
names(supp_agg) <- c("Drug", "Category", "Mean_support", "Min_support", "Pct_below_5pct")
print(supp_agg)
write.csv(supp_agg, file.path(out_dir, "support_diagnostics.csv"), row.names = FALSE)

cat(sprintf("\nResults saved to %s/\n", out_dir))

# ── Step 5: Bias figure ───────────────────────────────────────────────────
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  # Compute per-rep bias using empirical oracle
  store$truth <- mapply(get_oracle, store$drug, store$category, store$method)
  store$bias <- store$estimate - store$truth

  plot_data <- store[store$converged & !is.na(store$bias), ]
  plot_data$method <- factor(plot_data$method,
                             levels = c("naive", "ltmle", "lmtp_shift"))
  plot_data$drug_label <- ifelse(plot_data$drug == "Drug_A", "Drug A", "Drug B")

  fig <- ggplot(plot_data, aes(x = category, y = bias, fill = method)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
    facet_wrap(~ drug_label) +
    scale_fill_manual(values = c("naive" = "#D73027", "ltmle" = "#4575B4",
                                  "lmtp_shift" = "#1A9850"),
                      labels = c("Naive", "L-TMLE", "LMTP (+10%)")) +
    labs(x = "Adherence category / estimand", y = "Bias (estimate - oracle)",
         fill = "Method",
         title = sprintf("MC Bias Comparison (B=%d, n=%d)", B, n_sim),
         subtitle = "Primary oracle: empirical within-category") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1),
          legend.position = "bottom")

  ggsave(file.path(out_dir, "mc_bias_boxplot.pdf"), fig, width = 10, height = 6)
  ggsave(file.path(out_dir, "mc_bias_boxplot.png"), fig, width = 10, height = 6, dpi = 150)
  cat("Bias figure saved.\n")
}
