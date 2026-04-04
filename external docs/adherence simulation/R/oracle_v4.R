###############################################################################
# oracle_v4.R — Oracle truth computation via large Monte Carlo from known DGP
#
# Provides "ground truth" estimands by simulating from the interventional
# distribution with very large n.  Used to assess bias of estimators.
###############################################################################

source_dgp <- function() {
  # Locate DGP relative to this file or working directory
  candidates <- c(
    file.path("R", "sim_v4_dgp.R"),
    file.path("simulation V4", "R", "sim_v4_dgp.R"),
    "sim_v4_dgp.R"
  )
  for (f in candidates) {
    if (file.exists(f)) { source(f, local = FALSE); return(invisible()) }
  }
  stop("Cannot find sim_v4_dgp.R. Source it manually before calling oracle functions.")
}


# ── Oracle 1: Static continuous PDC within a drug ───────────────────────────
# Fix PDC to a constant value at every block for a given drug.
# Returns risk = P(failure by block K) under this deterministic intervention.
oracle_static_pdc_v4 <- function(pdc_value,
                                 drug,
                                 K        = 4,
                                 n_oracle = 50000,
                                 seed     = 9999,
                                 params   = NULL) {
  if (!exists("simulate_v4_interventional", mode = "function")) source_dgp()

  intervene <- function(nat_pdc, t, A0, W0, L_t) {
    rep(pdc_value, length(nat_pdc))
  }

  res <- simulate_v4_interventional(
    n = n_oracle, K = K, seed = seed,
    drug = drug, intervene_pdc = intervene, params = params
  )
  res$risk
}


# ── Oracle 2a: Static category regime — UNIFORM within category ─────────────
# "Always high", "always mid", "always low" — draw PDC uniformly from the
# category range at every block (stochastic within-category intervention).
oracle_static_category_v4 <- function(category,
                                      drug,
                                      K        = 4,
                                      n_oracle = 50000,
                                      seed     = 9999,
                                      params   = NULL) {
  if (!exists("simulate_v4_interventional", mode = "function")) source_dgp()

  p <- default_params_v4()
  if (!is.null(params)) for (nm in names(params)) p[[nm]] <- params[[nm]]

  lo <- switch(category,
               low  = 0.01,
               mid  = p$pdc_cut_low,
               high = p$pdc_cut_high)
  hi <- switch(category,
               low  = p$pdc_cut_low - 0.001,
               mid  = p$pdc_cut_high - 0.001,
               high = 1.00)

  intervene <- function(nat_pdc, t, A0, W0, L_t) {
    runif(length(nat_pdc), lo, hi)
  }

  res <- simulate_v4_interventional(
    n = n_oracle, K = K, seed = seed,
    drug = drug, intervene_pdc = intervene, params = params
  )
  res$risk
}


# ── Oracle 2b: Static category regime — MIDPOINT of category ───────────────
# Fix PDC to the midpoint of each category at every block.
# Midpoints: low = 0.33, mid = 0.75, high = 0.925
oracle_static_category_midpoint_v4 <- function(category,
                                               drug,
                                               K        = 4,
                                               n_oracle = 50000,
                                               seed     = 9999,
                                               params   = NULL) {
  if (!exists("simulate_v4_interventional", mode = "function")) source_dgp()

  p <- default_params_v4()
  if (!is.null(params)) for (nm in names(params)) p[[nm]] <- params[[nm]]

  midpoint <- switch(category,
                     low  = (0.01 + p$pdc_cut_low) / 2,
                     mid  = (p$pdc_cut_low + p$pdc_cut_high) / 2,
                     high = (p$pdc_cut_high + 1.0) / 2)

  oracle_static_pdc_v4(midpoint, drug = drug, K = K,
                       n_oracle = n_oracle, seed = seed, params = params)
}


# ── Oracle 2c: Static category regime — EMPIRICAL within category ───────────
# At each block, draw PDC from the empirical observed distribution within the
# target category. This best matches what the categorical L-TMLE targets:
# the intervention "set category = c" maps to the within-category distribution
# actually observed in the data.
#
# ref_data: a data frame from simulate_v4_data() used to build the empirical
#           distribution. If NULL, generates one internally.
oracle_static_category_empirical_v4 <- function(category,
                                                drug,
                                                K        = 4,
                                                n_oracle = 50000,
                                                seed     = 9999,
                                                params   = NULL,
                                                ref_data = NULL) {
  if (!exists("simulate_v4_interventional", mode = "function")) source_dgp()

  p <- default_params_v4()
  if (!is.null(params)) for (nm in names(params)) p[[nm]] <- params[[nm]]

  # Build empirical PDC pools by block from reference data
  if (is.null(ref_data)) {
    set.seed(seed + 5555)
    ref_data <- simulate_v4_data(n = 10000, K = K, seed = seed + 5555, params = params)
  }

  # Collect PDC values within the target category, by block
  pdc_pools <- vector("list", K)
  for (t in seq_len(K)) {
    pdc_col <- paste0("PDC_", t)
    cat_col <- paste0("PDC_cat_", t)
    vals <- ref_data[[pdc_col]][as.character(ref_data[[cat_col]]) == category &
                                  !is.na(ref_data[[pdc_col]])]
    if (length(vals) < 10) {
      warning(sprintf("oracle_empirical: only %d obs in category '%s' at block %d",
                      length(vals), category, t))
    }
    pdc_pools[[t]] <- vals
  }

  intervene <- function(nat_pdc, t, A0, W0, L_t) {
    pool <- pdc_pools[[t]]
    sample(pool, length(nat_pdc), replace = TRUE)
  }

  res <- simulate_v4_interventional(
    n = n_oracle, K = K, seed = seed,
    drug = drug, intervene_pdc = intervene, params = params
  )
  res$risk
}


# ── Oracle 3: Feasible additive shift ───────────────────────────────────────
# Increase each person's natural PDC by +delta, capped at 1.
# This is the estimand targeted by LMTP additive shift.
oracle_shift_v4 <- function(delta    = 0.10,
                            drug     = NULL,
                            K        = 4,
                            n_oracle = 50000,
                            seed     = 9999,
                            params   = NULL) {
  if (!exists("simulate_v4_interventional", mode = "function")) source_dgp()

  intervene <- function(nat_pdc, t, A0, W0, L_t) {
    pmin(nat_pdc + delta, 1.0)
  }

  res <- simulate_v4_interventional(
    n = n_oracle, K = K, seed = seed,
    drug = drug, intervene_pdc = intervene, params = params
  )
  res$risk
}


# ── Oracle 4: Natural (observational) risk by drug ─────────────────────────
oracle_natural_v4 <- function(drug     = NULL,
                              K        = 4,
                              n_oracle = 50000,
                              seed     = 9999,
                              params   = NULL) {
  if (!exists("simulate_v4_interventional", mode = "function")) source_dgp()

  res <- simulate_v4_interventional(
    n = n_oracle, K = K, seed = seed,
    drug = drug, intervene_pdc = NULL, params = params
  )
  res$risk
}


# ── Oracle 5: Matched-adherence drug comparison ─────────────────────────────
# Compare Drug A vs Drug B under the SAME adherence distribution.
# Approach: simulate both drugs with natural DGP, then swap adherence.
# "Under Drug B adherence": both drugs receive Drug B's natural PDC process.
#
# Implemented by:
#   1. Simulate a large dataset for each drug (natural).
#   2. Record the natural PDC sequences for Drug B.
#   3. Re-simulate Drug A using Drug B's PDC sequences (replacing its natural PDC).
#   4. Compare risks.
oracle_matched_adherence_v4 <- function(reference_drug = 0,
                                        K             = 4,
                                        n_oracle      = 50000,
                                        seed          = 9999,
                                        params        = NULL) {
  if (!exists("simulate_v4_interventional", mode = "function")) source_dgp()

  p <- default_params_v4()
  if (!is.null(params)) for (nm in names(params)) p[[nm]] <- params[[nm]]

  set.seed(seed)

  # --- Step 1: simulate reference drug naturally to get PDC sequences ---
  ref_seed <- seed
  W0_ref <- rnorm(n_oracle, 0, 1)
  A0_ref <- rep(as.integer(reference_drug), n_oracle)
  u_ref  <- rnorm(n_oracle, 0, p$sigma_u)
  mu_ref <- p$alpha0 + p$alphaA * A0_ref + p$alphaW * W0_ref + u_ref

  # Record PDC sequences from reference drug
  PDC_ref_seq <- matrix(NA_real_, n_oracle, K)
  L_prev <- rep(NA_real_, n_oracle)
  PDC_prev <- rep(NA_real_, n_oracle)
  eta_prev <- rep(NA_real_, n_oracle)

  for (t in seq_len(K)) {
    noise_L <- rnorm(n_oracle, 0, p$sigma_L)
    if (t == 1) {
      L_t <- p$beta0 + p$betaW * W0_ref + p$betaA * A0_ref + noise_L
    } else {
      L_t <- p$beta0 + p$betaW * W0_ref + p$betaA * A0_ref +
        p$betaLag * L_prev + p$betaP * (1 - PDC_prev) + noise_L
    }
    noise_e <- rnorm(n_oracle, 0, p$sigma_e)
    if (t == 1) {
      eta_t <- mu_ref + noise_e
    } else {
      eta_t <- (1 - p$rho) * mu_ref + p$rho * eta_prev +
        p$alphaL * L_prev + noise_e
    }
    PDC_t <- plogis(eta_t)
    PDC_ref_seq[, t] <- PDC_t

    L_prev   <- L_t
    PDC_prev <- PDC_t
    eta_prev <- eta_t
  }

  # --- Step 2: risk under reference drug with its own adherence ---
  risk_ref <- oracle_static_drug_with_pdc_seq(
    drug = reference_drug, pdc_seq = PDC_ref_seq,
    W0 = W0_ref, K = K, params = p, seed = seed + 1
  )

  # --- Step 3: risk under other drug with REFERENCE drug's adherence ---
  other_drug <- 1 - reference_drug
  risk_other <- oracle_static_drug_with_pdc_seq(
    drug = other_drug, pdc_seq = PDC_ref_seq,
    W0 = W0_ref, K = K, params = p, seed = seed + 2
  )

  list(
    risk_ref        = risk_ref,
    risk_other      = risk_other,
    rd              = risk_other - risk_ref,
    reference_drug  = reference_drug,
    reference_label = ifelse(reference_drug == 0, "Drug B", "Drug A"),
    other_label     = ifelse(reference_drug == 0, "Drug A", "Drug B")
  )
}


# ── Helper: compute risk for a given drug using pre-specified PDC sequences ─
oracle_static_drug_with_pdc_seq <- function(drug, pdc_seq, W0, K, params, seed) {
  p <- params
  n <- length(W0)
  A0 <- rep(as.integer(drug), n)

  set.seed(seed)
  failed <- rep(FALSE, n)
  Y_mat  <- matrix(0L, n, K)

  L_prev   <- rep(NA_real_, n)
  PDC_prev <- rep(NA_real_, n)


  for (t in seq_len(K)) {
    # L_t depends on the intervened PDC (from pdc_seq)
    noise_L <- rnorm(n, 0, p$sigma_L)
    if (t == 1) {
      L_t <- p$beta0 + p$betaW * W0 + p$betaA * A0 + noise_L
    } else {
      L_t <- p$beta0 + p$betaW * W0 + p$betaA * A0 +
        p$betaLag * L_prev + p$betaP * (1 - PDC_prev) + noise_L
    }

    PDC_t <- pdc_seq[, t]

    # Failure
    at_risk <- !failed
    if (any(at_risk)) {
      lp <- p$gamma0 +
        p$gammaL * L_t[at_risk] +
        p$gammaW * W0[at_risk] +
        p$gammaP * (1 - PDC_t[at_risk]) +
        p$gammaA * A0[at_risk] +
        p$gammaAP * A0[at_risk] * (1 - PDC_t[at_risk])
      prob <- plogis(lp)
      new_fail <- rbinom(sum(at_risk), 1, prob) == 1L
      idx <- which(at_risk)
      Y_mat[idx[new_fail], t] <- 1L
      failed[idx[new_fail]] <- TRUE
    }
    if (t > 1) Y_mat[Y_mat[, t - 1] == 1L, t] <- 1L

    L_prev   <- L_t
    PDC_prev <- PDC_t
  }

  mean(Y_mat[, K] == 1L)
}


# ── Convenience: compute a grid of oracle values ────────────────────────────
oracle_grid_v4 <- function(K = 4, n_oracle = 50000, seed = 9999, params = NULL) {

  if (!exists("default_params_v4", mode = "function")) source_dgp()

  cat("Computing oracle grid (this may take a minute)...\n")

  results <- list()

  # Natural risks
  for (d in 0:1) {
    r <- oracle_natural_v4(drug = d, K = K, n_oracle = n_oracle,
                           seed = seed, params = params)
    results[[length(results) + 1]] <- data.frame(
      estimand = "natural", drug = d,
      drug_label = ifelse(d == 0, "Drug B", "Drug A"),
      detail = "observed", risk = r, stringsAsFactors = FALSE
    )
  }

  # Static category regimes
  for (d in 0:1) {
    for (cat in c("low", "mid", "high")) {
      r <- oracle_static_category_v4(
        category = cat, drug = d, K = K,
        n_oracle = n_oracle, seed = seed + d * 10 + match(cat, c("low","mid","high")),
        params = params
      )
      results[[length(results) + 1]] <- data.frame(
        estimand = "static_category", drug = d,
        drug_label = ifelse(d == 0, "Drug B", "Drug A"),
        detail = cat, risk = r, stringsAsFactors = FALSE
      )
    }
  }

  # Additive shifts
  for (d in 0:1) {
    for (delta in c(0.05, 0.10, 0.15)) {
      r <- oracle_shift_v4(delta = delta, drug = d, K = K,
                           n_oracle = n_oracle,
                           seed = seed + 100 + d * 10 + delta * 100,
                           params = params)
      results[[length(results) + 1]] <- data.frame(
        estimand = "shift", drug = d,
        drug_label = ifelse(d == 0, "Drug B", "Drug A"),
        detail = sprintf("+%.0f%%", delta * 100), risk = r,
        stringsAsFactors = FALSE
      )
    }
  }

  # Matched-adherence comparison
  matched <- oracle_matched_adherence_v4(
    reference_drug = 0, K = K, n_oracle = n_oracle,
    seed = seed + 200, params = params
  )
  results[[length(results) + 1]] <- data.frame(
    estimand = "matched_adherence", drug = 0,
    drug_label = "Drug B (ref)", detail = "under Drug B adherence",
    risk = matched$risk_ref, stringsAsFactors = FALSE
  )
  results[[length(results) + 1]] <- data.frame(
    estimand = "matched_adherence", drug = 1,
    drug_label = "Drug A", detail = "under Drug B adherence",
    risk = matched$risk_other, stringsAsFactors = FALSE
  )

  grid <- do.call(rbind, results)
  rownames(grid) <- NULL

  cat("Oracle grid complete.\n")
  grid
}
