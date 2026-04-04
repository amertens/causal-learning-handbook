###############################################################################
# 06_ltmle_fine_grid.R — LTMLE on fine 10%-width PDC categories
#
# Runs LTMLE for each drug on 7 adherence categories:
#   [0.30,0.40), [0.40,0.50), [0.50,0.60), [0.60,0.70),
#   [0.70,0.80), [0.80,0.90), [0.90,1.00]
#
# Uses a larger dataset (n=10000) for adequate support in sparse bins.
# Saves results to results/ltmle_fine_grid.rds for use in reports.
#
# Usage:
#   Rscript scripts/06_ltmle_fine_grid.R          # default n=10000
#   Rscript scripts/06_ltmle_fine_grid.R 15000    # larger n
###############################################################################

cat("=== LTMLE Fine Grid (10%-width PDC bins) ===\n\n")

source(file.path("R", "sim_v4_dgp.R"))
source(file.path("R", "oracle_v4.R"))

if (!requireNamespace("ltmle", quietly = TRUE))
  stop("Package 'ltmle' required")
library(ltmle)

# ── Config ──────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
N_SIM <- if (length(args) >= 1) as.integer(args[1]) else 10000L
K     <- 4L
SEED  <- 2024L

cat(sprintf("n = %d, K = %d, seed = %d\n\n", N_SIM, K, SEED))

# ── Define fine bins ────────────────────────────────────────────────────────
bin_breaks <- seq(0.30, 1.00, by = 0.10)
n_bins     <- length(bin_breaks) - 1  # 7 bins
bin_labels <- sprintf("[%.0f%%,%.0f%%)", bin_breaks[-length(bin_breaks)] * 100,
                      bin_breaks[-1] * 100)
bin_labels[n_bins] <- sub("\\)", "]", bin_labels[n_bins])  # last bin is closed
bin_midpoints <- (bin_breaks[-length(bin_breaks)] + bin_breaks[-1]) / 2

cat("Bins:", paste(bin_labels, collapse = ", "), "\n")
cat("Midpoints:", paste(bin_midpoints, collapse = ", "), "\n\n")

# Number of dummy treatment nodes per block = n_bins - 1 (reference = lowest bin)
n_dummies <- n_bins - 1  # 6 dummies per block

# ── Simulate data ──────────────────────────────────────────────────────────
dat <- simulate_v4_data(n = N_SIM, K = K, seed = SEED)
cat(sprintf("Data: n=%d, failure=%.1f%%\n\n", nrow(dat), 100 * mean(dat$Y_final)))

# ── Helper: categorise PDC into fine bins ───────────────────────────────────
categorise_fine <- function(pdc) {
  # Returns integer 1..n_bins; values below 0.30 go to bin 1, above 1.0 to last
  out <- findInterval(pdc, bin_breaks, rightmost.closed = TRUE, left.open = FALSE)
  out[out < 1] <- 1L
  out[out > n_bins] <- n_bins
  as.integer(out)
}

# ── Build LTMLE data with fine-bin dummy encoding ───────────────────────────
build_ltmle_fine <- function(df_sub, K, n_bins) {
  n <- nrow(df_sub)
  n_dummies <- n_bins - 1

  out <- data.frame(W0 = df_sub$W0)

  for (t in seq_len(K)) {
    # L node
    out[[paste0("L_", t)]] <- df_sub[[paste0("L_", t)]]

    # Categorise PDC
    pdc_vals <- df_sub[[paste0("PDC_", t)]]
    bin_int  <- categorise_fine(pdc_vals)

    # Create n_dummies dummy nodes: Abin2_t, Abin3_t, ..., Abin{n_bins}_t
    for (d in 2:n_bins) {
      dname <- sprintf("Abin%d_%d", d, t)
      out[[dname]] <- as.integer(!is.na(bin_int) & bin_int == d)
      out[[dname]][is.na(bin_int)] <- NA_integer_
    }

    # Outcome
    out[[paste0("Y_", t)]] <- df_sub[[paste0("Y_", t)]]
  }

  # Enforce absorbing state
  for (t in seq_len(K)) {
    if (t < K) {
      failed <- !is.na(out[[paste0("Y_", t)]]) & out[[paste0("Y_", t)]] == 1L
      for (s in (t + 1):K) {
        out[[paste0("L_", s)]][failed] <- NA_real_
        for (d in 2:n_bins) {
          out[[sprintf("Abin%d_%d", d, s)]][failed] <- NA_integer_
        }
        out[[paste0("Y_", s)]][failed] <- 1L
      }
    }
  }

  out
}

# ── Node lists ──────────────────────────────────────────────────────────────
make_nodes_fine <- function(K, n_bins) {
  Lnodes <- "W0"
  Anodes <- character()
  Ynodes <- character()

  for (t in seq_len(K)) {
    Lnodes <- c(Lnodes, paste0("L_", t))
    for (d in 2:n_bins) {
      Anodes <- c(Anodes, sprintf("Abin%d_%d", d, t))
    }
    Ynodes <- c(Ynodes, paste0("Y_", t))
  }

  list(Lnodes = Lnodes, Anodes = Anodes, Ynodes = Ynodes)
}

# ── Formulas ────────────────────────────────────────────────────────────────
make_formulas_fine <- function(K, n_bins) {
  q_list <- character()
  g_list <- character()

  for (t in seq_len(K)) {
    L_t <- paste0("L_", t)
    Y_t <- paste0("Y_", t)
    A_dummies_t <- sprintf("Abin%d_%d", 2:n_bins, t)

    # Q: Y_t ~ W0 + L_t + all current dummies
    q_list[Y_t] <- sprintf("Q.kplus1 ~ W0 + %s + %s",
                           L_t, paste(A_dummies_t, collapse = " + "))

    # g: each dummy ~ W0 + L_t + preceding dummies at same block
    for (di in seq_along(A_dummies_t)) {
      dname <- A_dummies_t[di]
      if (di == 1) {
        # First dummy: condition on W0 + L_t (+ previous block dummies if t>1)
        if (t == 1) {
          g_list[dname] <- sprintf("%s ~ W0 + %s", dname, L_t)
        } else {
          prev_dummies <- sprintf("Abin%d_%d", 2:n_bins, t - 1)
          g_list[dname] <- sprintf("%s ~ W0 + %s + %s", dname, L_t,
                                   paste(prev_dummies, collapse = " + "))
        }
      } else {
        # Later dummies: add preceding same-block dummies
        preceding <- A_dummies_t[1:(di - 1)]
        if (t == 1) {
          g_list[dname] <- sprintf("%s ~ W0 + %s + %s", dname, L_t,
                                   paste(preceding, collapse = " + "))
        } else {
          prev_dummies <- sprintf("Abin%d_%d", 2:n_bins, t - 1)
          g_list[dname] <- sprintf("%s ~ W0 + %s + %s + %s", dname, L_t,
                                   paste(prev_dummies, collapse = " + "),
                                   paste(preceding, collapse = " + "))
        }
      }
    }
  }

  list(Qform = q_list, gform = g_list)
}

# ── Deterministic g-function ────────────────────────────────────────────────
# At each block t, at most one dummy can be 1. If any earlier dummy is 1,
# all later dummies must be 0.
make_det_g_fine <- function(K, n_bins) {
  # Build lookup: for each dummy node, which earlier dummies block it
  blockers <- list()
  for (t in seq_len(K)) {
    dummies_t <- sprintf("Abin%d_%d", 2:n_bins, t)
    for (di in 2:length(dummies_t)) {
      blockers[[ dummies_t[di] ]] <- dummies_t[1:(di - 1)]
    }
  }

  function(data, current.node, nodes) {
    node_name <- names(data)[current.node]
    if (node_name %in% names(blockers)) {
      bl_cols <- blockers[[node_name]]
      is_blocked <- rep(FALSE, nrow(data))
      for (bc in bl_cols) {
        if (bc %in% names(data)) {
          is_blocked <- is_blocked | (!is.na(data[[bc]]) & data[[bc]] == 1L)
        }
      }
      if (any(is_blocked)) {
        return(list(is.deterministic = is_blocked, prob1 = 0))
      }
    }
    NULL
  }
}

# ── Build abar for a target bin ─────────────────────────────────────────────
# target_bin: integer 1..n_bins
make_abar_fine <- function(target_bin, K, n_bins) {
  # At each block: dummy d is 1 iff d == target_bin, else 0
  block_abar <- integer(n_bins - 1)
  if (target_bin >= 2) {
    block_abar[target_bin - 1] <- 1L
  }
  # All 0 if target_bin == 1 (reference category)
  rep(block_abar, K)
}

# ── Main loop: fit LTMLE for each drug × bin ────────────────────────────────
nodes  <- make_nodes_fine(K, n_bins)
forms  <- make_formulas_fine(K, n_bins)
det_g  <- make_det_g_fine(K, n_bins)

out_dir <- file.path("results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

results <- data.frame(
  drug = integer(), drug_label = character(),
  bin = integer(), bin_label = character(), midpoint = numeric(),
  estimate = numeric(), se = numeric(),
  oracle = numeric(), bias = numeric(),
  converged = logical(), error_msg = character(),
  n_in_bin = integer(),
  stringsAsFactors = FALSE
)

for (drug in 0:1) {
  dlabel <- ifelse(drug == 0, "Drug B", "Drug A")
  cat(sprintf("── %s ──\n", dlabel))

  # Subset to drug
  df_sub <- dat[dat$A0 == drug, , drop = FALSE]
  ldata  <- build_ltmle_fine(df_sub, K, n_bins)

  # Support check
  cat("  Support (block 1 bin counts):\n")
  pdc1 <- df_sub$PDC_1
  bins1 <- categorise_fine(pdc1)
  for (b in seq_len(n_bins)) {
    cat(sprintf("    Bin %d %s: n=%d (%.1f%%)\n", b, bin_labels[b],
                sum(bins1 == b, na.rm = TRUE),
                100 * mean(bins1 == b, na.rm = TRUE)))
  }

  for (b in seq_len(n_bins)) {
    cat(sprintf("  Fitting bin %d (%s, midpoint=%.2f)...", b, bin_labels[b], bin_midpoints[b]))

    abar <- make_abar_fine(b, K, n_bins)

    # Oracle at midpoint
    oracle_val <- oracle_static_pdc_v4(bin_midpoints[b], drug = drug, K = K,
                                       n_oracle = 50000,
                                       seed = SEED + drug * 1000 + b)

    res <- tryCatch({
      fit <- ltmle::ltmle(
        data       = ldata,
        Anodes     = nodes$Anodes,
        Lnodes     = nodes$Lnodes,
        Ynodes     = nodes$Ynodes,
        abar       = abar,
        Qform      = forms$Qform,
        gform      = forms$gform,
        SL.library = "glm",
        estimate.time = FALSE,
        survivalOutcome = TRUE,
        deterministic.g.function = det_g,
        variance.method = "ic"
      )
      summ <- summary(fit)
      est  <- as.numeric(fit$estimates["tmle"])
      se   <- summ$treatment$std.dev
      list(estimate = est, se = se, converged = TRUE, error_msg = "")
    }, error = function(e) {
      list(estimate = NA_real_, se = NA_real_, converged = FALSE,
           error_msg = conditionMessage(e))
    })

    est_str <- if (res$converged) sprintf("%.4f", res$estimate) else "FAILED"
    cat(sprintf(" %s (oracle=%.4f)\n", est_str, oracle_val))

    # Count subjects in this bin at block 1
    n_in <- sum(bins1 == b, na.rm = TRUE)

    results <- rbind(results, data.frame(
      drug = drug, drug_label = dlabel,
      bin = b, bin_label = bin_labels[b], midpoint = bin_midpoints[b],
      estimate = res$estimate, se = res$se,
      oracle = oracle_val, bias = res$estimate - oracle_val,
      converged = res$converged, error_msg = res$error_msg,
      n_in_bin = n_in,
      stringsAsFactors = FALSE))
  }
}

# ── Save ────────────────────────────────────────────────────────────────────
saveRDS(results, file.path(out_dir, "ltmle_fine_grid.rds"))
write.csv(results, file.path(out_dir, "ltmle_fine_grid.csv"), row.names = FALSE)

cat("\n=== Summary ===\n")
cat(sprintf("Converged: %d/%d\n", sum(results$converged), nrow(results)))
print(results[, c("drug_label", "bin_label", "midpoint", "estimate", "oracle", "bias", "converged")])

cat(sprintf("\nSaved to %s/ltmle_fine_grid.rds\n", out_dir))
