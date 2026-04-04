###############################################################################
# 05_mc_phase3_positivity_stress.R — Positivity stress-test MC (Phase 3)
#
# Purpose: show why LMTP is useful when static regimes become poorly supported
# Approach: vary the adherence intercept to create better vs worse support
# Default: 3 scenarios, n=2000, B=200 per scenario
# Output: support diagnostics, LTMLE/LMTP bias by scenario, figure
###############################################################################

cat("=== MC Phase 3: Positivity Stress Test ===\n")

source(file.path("R", "sim_v4_dgp.R"))
source(file.path("R", "oracle_v4.R"))
source(file.path("R", "prepare_ltmle_v4.R"))
source(file.path("R", "prepare_lmtp_v4.R"))
source(file.path("R", "diagnostics_v4.R"))

for (pkg in c("ltmle", "lmtp"))
  if (!requireNamespace(pkg, quietly = TRUE)) stop(sprintf("Package '%s' required", pkg))

args <- commandArgs(trailingOnly = TRUE)
B      <- if (length(args) >= 1) as.integer(args[1]) else 200L
n_sim  <- 2000L
K      <- 4L
seed0  <- 20260601L

out_dir <- file.path("results", "mc_phase3")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Scenarios: vary adherence intercept ─────────────────────────────────────
# Higher alpha0 → more high-adherence patients → better support for always_high,
#   worse for always_low
# Lower alpha0 → more spread → better support for always_low,
#   worse for always_high
scenarios <- list(
  easy   = list(label = "Easy support (alpha0=1.0)",  alpha0 = 1.0),
  medium = list(label = "Medium (alpha0=1.3, default)", alpha0 = 1.3),
  hard   = list(label = "Hard support (alpha0=1.6)",  alpha0 = 1.6)
)

# ── MC loop by scenario ────────────────────────────────────────────────────
all_results <- list()

for (sc_name in names(scenarios)) {
  sc <- scenarios[[sc_name]]
  cat(sprintf("\n--- Scenario: %s ---\n", sc$label))

  my_params <- default_params_v4()
  my_params$alpha0 <- sc$alpha0

  # Oracle for this scenario
  oracle_low <- oracle_static_category_v4("low", drug = 1, K = K,
    n_oracle = 100000, seed = seed0 + match(sc_name, names(scenarios)), params = my_params)
  oracle_shift <- oracle_shift_v4(0.10, drug = 1, K = K,
    n_oracle = 100000, seed = seed0 + 100 + match(sc_name, names(scenarios)), params = my_params)

  store <- data.frame(rep = integer(), method = character(),
                      estimate = numeric(), converged = logical(),
                      stringsAsFactors = FALSE)

  for (b in seq_len(B)) {
    if (b %% 50 == 0) cat(sprintf("  Rep %d/%d\n", b, B))
    dat <- simulate_v4_data(n = n_sim, K = K, seed = seed0 + b * 19 + match(sc_name, names(scenarios)) * 1000,
                            params = my_params)

    # LTMLE always_low Drug A
    res_l <- fit_ltmle_v4_static(dat, drug = 1, regime = "always_low", K = K,
                                  sl_lib = "glm", verbose = FALSE)
    store <- rbind(store, data.frame(rep = b, method = "ltmle_low",
      estimate = res_l$estimate, converged = res_l$converged, stringsAsFactors = FALSE))

    # LMTP shift Drug A
    res_s <- fit_lmtp_shift_v4(dat, delta = 0.10, drug = 1, K = K, folds = 1)
    store <- rbind(store, data.frame(rep = b, method = "lmtp_shift",
      estimate = res_s$estimate, converged = res_s$converged, stringsAsFactors = FALSE))
  }

  # Summarise
  for (m in c("ltmle_low", "lmtp_shift")) {
    sub <- store[store$method == m & store$converged & !is.na(store$estimate), ]
    truth <- if (m == "ltmle_low") oracle_low else oracle_shift
    if (nrow(sub) > 0) {
      all_results[[length(all_results) + 1]] <- data.frame(
        scenario = sc$label, method = m, truth = truth,
        mean_est = mean(sub$estimate), bias = mean(sub$estimate) - truth,
        emp_se = sd(sub$estimate), n_converged = nrow(sub),
        conv_rate = sprintf("%.0f%%", 100 * nrow(sub) / B),
        stringsAsFactors = FALSE)
    }
  }

  # Support diagnostic for one dataset
  dat_diag <- simulate_v4_data(n = n_sim, K = K, seed = seed0, params = my_params)
  supp <- tbl_support_detail_v4(dat_diag, K = K)
  cat("  Support (Drug A):\n")
  print(supp[supp$Drug == "Drug A", c("Regime", "Pct_follow", "Flag")])
}

# ── Save results ────────────────────────────────────────────────────────────
phase3 <- do.call(rbind, all_results)
write.csv(phase3, file.path(out_dir, "phase3_summary.csv"), row.names = FALSE)

cat("\n=== Phase 3 Summary ===\n")
print(phase3)

# ── Figure ──────────────────────────────────────────────────────────────────
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  phase3$method_label <- ifelse(phase3$method == "ltmle_low",
                                "LTMLE (always low)", "LMTP (+10% shift)")

  fig <- ggplot(phase3, aes(x = scenario, y = abs(bias), fill = method_label)) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    geom_text(aes(label = conv_rate), position = position_dodge(0.8),
              vjust = -0.5, size = 3) +
    scale_fill_manual(values = c("LTMLE (always low)" = "#4575B4",
                                  "LMTP (+10% shift)" = "#1A9850")) +
    labs(x = "Scenario", y = "|Bias|", fill = "Method",
         title = "Positivity Stress Test: LTMLE vs LMTP",
         subtitle = "Labels show convergence rate") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 15, hjust = 1))
  ggsave(file.path(out_dir, "phase3_stress_test.pdf"), fig, width = 9, height = 5)
  cat("Figure saved\n")
}
