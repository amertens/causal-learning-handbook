###############################################################################
# sim_v4_dgp.R — Simplified longitudinal DGP for HIV adherence simulation V4
#
# Generates 90-day block data with:
#   - W0: baseline severity (continuous)
#   - A0: baseline regimen indicator (0 = Drug B, 1 = Drug A)
#   - L_t: time-varying health confounder (continuous, per block)
#   - PDC_t: continuous proportion of days covered in [0,1] (per block)
#   - PDC_cat_t: categorised adherence (low / mid / high)
#   - C_t: optional censoring indicator (1 = observed, 0 = censored)
#   - Y_t: absorbing virologic-failure indicator
#
# Causal structure per block (t = 1, …, K):
#   W0, A0 → L_t → PDC_t → Y_t
#   with feedback: PDC_{t-1} → L_t → PDC_t  (time-dependent confounding)
#
# Design goals:
#   - simple GLM-based LTMLE ≈ unbiased in clean setting
#   - static always-high / always-low regimes show positivity strain
#   - LMTP shift interventions remain feasible
#   - Drug A is more forgiving at low adherence; drugs similar at high PDC
###############################################################################

# ── Default parameter list ──────────────────────────────────────────────────
default_params_v4 <- function() {
  list(
    # ── Adherence process (logit-normal, mean-reverting AR(1)) ──
    #   mu_i = alpha0 + alphaA*A0 + alphaW*W0 + u_i    (subject mean)
    #   eta_1 = mu_i + e_1
    #   eta_t = (1-rho)*mu_i + rho*eta_{t-1} + alphaL*L_{t-1} + e_t
    #   PDC_t = plogis(eta_t)
    alpha0   = 1.30,     # intercept → plogis(1.3) ≈ 0.79 baseline mean
    alphaA   = -0.35,    # Drug A shifts adherence lower
    alphaW   = -0.20,    # worse baseline health → lower adherence
    alphaL   = -0.30,    # higher L (sicker) → lower adherence
    rho      = 0.45,     # AR(1) persistence on logit scale
    sigma_u  = 0.50,     # between-subject random-intercept SD
    sigma_e  = 0.65,     # within-subject innovation SD

    # ── Time-varying covariate L_t ──
    #   L_1 = beta0 + betaW*W0 + betaA*A0 + noise
    #   L_t = beta0 + betaW*W0 + betaA*A0 + betaLag*L_{t-1}
    #         + betaP*(1 - PDC_{t-1}) + noise
    beta0    = 0.0,      # intercept (centered)
    betaW    = 0.30,     # baseline severity → worse health
    betaA    = 0.0,      # no direct drug effect on L
    betaLag  = 0.35,     # persistence of L
    betaP    = 0.80,     # poor adherence worsens future health
    sigma_L  = 0.50,     # noise SD

    # ── Failure hazard ──
    #   logit P(Y_t = 1 | not yet failed)
    #     = gamma0 + gammaL*L_t + gammaW*W0
    #       + gammaP*(1-PDC_t) + gammaA*A0 + gammaAP*A0*(1-PDC_t)
    gamma0   = -3.80,    # intercept (low baseline hazard)
    gammaL   = 0.50,     # sicker → higher hazard
    gammaW   = 0.30,     # baseline severity → higher hazard
    gammaP   = 2.50,     # poor adherence → higher hazard
    gammaA   = 0.0,      # no main drug effect (same at high PDC)
    gammaAP  = -1.20,    # Drug A × poor adherence (protective = forgiving)

    # ── Censoring (optional) ──
    cens0    = -4.50,    # intercept (low baseline censoring)
    censL    = 0.30,     # sicker → more LTFU
    censW    = 0.10,     # severity → slightly more LTFU

    # ── PDC category cutpoints ──
    pdc_cut_low  = 0.65,
    pdc_cut_high = 0.85
  )
}


# ── Main simulation function ────────────────────────────────────────────────
simulate_v4_data <- function(n             = 2000,
                             K             = 4,
                             seed          = 123,
                             include_censoring = FALSE,
                             return_long   = FALSE,
                             params        = NULL) {
  # Packages
  if (!requireNamespace("stats", quietly = TRUE))
    stop("Base R stats package required")


  # Merge user-supplied params over defaults
  p <- default_params_v4()
  if (!is.null(params)) {
    for (nm in names(params)) p[[nm]] <- params[[nm]]
  }

  set.seed(seed)

  # ── Baseline variables ──────────────────────────────────────────────────
  W0 <- rnorm(n, 0, 1)                     # baseline severity

  A0 <- rbinom(n, 1, 0.5)                  # regimen: 0 = Drug B, 1 = Drug A

  # Subject-specific adherence random intercept
  u_i <- rnorm(n, 0, p$sigma_u)

  # Subject-specific adherence long-run mean on logit scale
  mu_i <- p$alpha0 + p$alphaA * A0 + p$alphaW * W0 + u_i

  # ── Storage ─────────────────────────────────────────────────────────────
  # Matrices: rows = subjects, cols = blocks
  L_mat       <- matrix(NA_real_, n, K)
  PDC_mat     <- matrix(NA_real_, n, K)
  PDC_cat_mat <- matrix(NA_character_, n, K)
  PDC_int_mat <- matrix(NA_integer_, n, K)
  Y_mat       <- matrix(0L, n, K)
  C_mat       <- matrix(1L, n, K)        # 1 = observed by default
  eta_mat     <- matrix(NA_real_, n, K)  # latent adherence propensity

  failed   <- rep(FALSE, n)              # absorbing failure flag
  censored <- rep(FALSE, n)              # absorbing censoring flag

  for (t in seq_len(K)) {

    # ── Time-varying covariate L_t ──────────────────────────────────────
    noise_L <- rnorm(n, 0, p$sigma_L)
    if (t == 1) {
      L_mat[, t] <- p$beta0 + p$betaW * W0 + p$betaA * A0 + noise_L
    } else {
      L_mat[, t] <- p$beta0 + p$betaW * W0 + p$betaA * A0 +
        p$betaLag * L_mat[, t - 1] +
        p$betaP * (1 - PDC_mat[, t - 1]) + noise_L
    }

    # Enforce NA after censoring
    if (include_censoring) {
      L_mat[censored, t] <- NA_real_
    }

    # ── Continuous adherence PDC_t ──────────────────────────────────────
    noise_e <- rnorm(n, 0, p$sigma_e)
    if (t == 1) {
      eta_mat[, t] <- mu_i + noise_e
    } else {
      eta_mat[, t] <- (1 - p$rho) * mu_i +
        p$rho * eta_mat[, t - 1] +
        p$alphaL * L_mat[, t - 1] + noise_e
      # note: L_{t-1} from the *current* block's L generation uses t-1
      # but here alphaL uses L_{t-1} which was generated in the prior block
      # At t>1 we use L from the PREVIOUS block to affect adherence
      # (L_t for the current block was just generated above, but we use
      # L_{t-1} for the adherence model to maintain temporal ordering:
      #   L_{t-1} → PDC_t, and also L_t → PDC_t is not in the model)
      # Actually: re-read the spec more carefully.
      # The spec says alphaL*L_{t-1}. L_{t} is generated first in each block,
      # then PDC_t. But the adherence model uses L_{t-1}, the PRIOR block's L.
      # This is correct: within a block, L_t and PDC_t are generated, but

      # adherence depends on the PRIOR block's health status, not the current.
      # The current L_t then affects the failure hazard for this block.
    }

    PDC_mat[, t] <- plogis(eta_mat[, t])

    # Enforce NA after censoring
    if (include_censoring) {
      PDC_mat[censored, t] <- NA_real_
      eta_mat[censored, t] <- NA_real_
    }

    # ── Categorise PDC ──────────────────────────────────────────────────
    pdc_vals <- PDC_mat[, t]
    cat_vec  <- rep(NA_character_, n)
    int_vec  <- rep(NA_integer_, n)
    obs <- !is.na(pdc_vals)
    cat_vec[obs & pdc_vals <  p$pdc_cut_low]  <- "low"
    cat_vec[obs & pdc_vals >= p$pdc_cut_low &
              pdc_vals <  p$pdc_cut_high] <- "mid"
    cat_vec[obs & pdc_vals >= p$pdc_cut_high] <- "high"
    int_vec[cat_vec == "low"]  <- 1L
    int_vec[cat_vec == "mid"]  <- 2L
    int_vec[cat_vec == "high"] <- 3L
    PDC_cat_mat[, t] <- cat_vec
    PDC_int_mat[, t] <- int_vec

    # ── Virologic failure Y_t (absorbing) ───────────────────────────────
    at_risk <- !failed & !censored
    if (any(at_risk)) {
      lp_fail <- p$gamma0 +
        p$gammaL * L_mat[at_risk, t] +
        p$gammaW * W0[at_risk] +
        p$gammaP * (1 - PDC_mat[at_risk, t]) +
        p$gammaA * A0[at_risk] +
        p$gammaAP * A0[at_risk] * (1 - PDC_mat[at_risk, t])

      prob_fail <- plogis(lp_fail)
      new_fail  <- rbinom(sum(at_risk), 1, prob_fail) == 1L
      idx_risk  <- which(at_risk)
      Y_mat[idx_risk[new_fail], t] <- 1L
      failed[idx_risk[new_fail]]   <- TRUE
    }

    # Carry forward failure (absorbing)
    if (t > 1) {
      still_failed_before <- Y_mat[, t - 1] == 1L & !censored
      Y_mat[still_failed_before, t] <- 1L
    }

    # Enforce NA after censoring
    if (include_censoring) {
      Y_mat[censored, t] <- NA_integer_
    }

    # ── Optional censoring C_t ──────────────────────────────────────────
    # C_t = 1 means subject is still observed at the *next* block.
    # Censoring is generated after Y_t; if censored, next block is NA.
    if (include_censoring && t < K) {
      at_risk_cens <- !failed & !censored
      if (any(at_risk_cens)) {
        lp_cens <- p$cens0 +
          p$censL * L_mat[at_risk_cens, t] +
          p$censW * W0[at_risk_cens]
        prob_cens <- plogis(lp_cens)
        new_cens  <- rbinom(sum(at_risk_cens), 1, prob_cens) == 1L
        idx_cens  <- which(at_risk_cens)
        C_mat[idx_cens[new_cens], t] <- 0L
        censored[idx_cens[new_cens]] <- TRUE
      }
    }
  }

  # ── Assemble wide data frame ──────────────────────────────────────────
  df <- data.frame(id = seq_len(n), W0 = W0, A0 = A0)

  for (t in seq_len(K)) {
    df[[paste0("L_", t)]]       <- L_mat[, t]
    df[[paste0("PDC_", t)]]     <- PDC_mat[, t]
    df[[paste0("PDC_cat_", t)]] <- factor(PDC_cat_mat[, t],
                                          levels = c("low", "mid", "high"),
                                          ordered = TRUE)
    df[[paste0("PDC_int_", t)]] <- PDC_int_mat[, t]
    if (include_censoring) {
      df[[paste0("C_", t)]]     <- C_mat[, t]
    }
    df[[paste0("Y_", t)]]       <- Y_mat[, t]
  }

  # Final cumulative outcome: 1 if failed at any block
  df$Y_final <- as.integer(rowSums(Y_mat, na.rm = TRUE) > 0)

  # Store params as attribute for reproducibility
  attr(df, "dgp_params") <- p
  attr(df, "K")           <- K

  # ── Optional long format ──────────────────────────────────────────────
  if (return_long) {
    long_list <- vector("list", K)
    for (t in seq_len(K)) {
      long_list[[t]] <- data.frame(
        id       = seq_len(n),
        block    = t,
        W0       = W0,
        A0       = A0,
        L        = L_mat[, t],
        PDC      = PDC_mat[, t],
        PDC_cat  = PDC_cat_mat[, t],
        PDC_int  = PDC_int_mat[, t],
        C        = C_mat[, t],
        Y        = Y_mat[, t],
        failed   = as.integer(failed),
        stringsAsFactors = FALSE
      )
    }
    long_df <- do.call(rbind, long_list)
    long_df <- long_df[order(long_df$id, long_df$block), ]
    rownames(long_df) <- NULL
    attr(df, "long") <- long_df
  }

  df
}


# ── Helper: categorise a PDC vector ─────────────────────────────────────────
categorise_pdc_v4 <- function(pdc, params = NULL) {
  p <- if (is.null(params)) default_params_v4() else params
  cut_low  <- p$pdc_cut_low
  cut_high <- p$pdc_cut_high
  out <- rep(NA_character_, length(pdc))
  out[!is.na(pdc) & pdc <  cut_low]                  <- "low"
  out[!is.na(pdc) & pdc >= cut_low & pdc < cut_high] <- "mid"
  out[!is.na(pdc) & pdc >= cut_high]                  <- "high"
  factor(out, levels = c("low", "mid", "high"), ordered = TRUE)
}


# ── Helper: simulate interventional data from the DGP ───────────────────────
# Used by oracle functions. Allows replacing the natural PDC process
# with an arbitrary intervention function.
#
# intervene_pdc: function(natural_pdc, t, A0, W0, L_t) → intervened PDC
#   If NULL, uses natural PDC (observational).
simulate_v4_interventional <- function(n             = 50000,
                                       K             = 4,
                                       seed          = 999,
                                       drug          = NULL,
                                       intervene_pdc = NULL,
                                       params        = NULL) {
  p <- default_params_v4()
  if (!is.null(params)) {
    for (nm in names(params)) p[[nm]] <- params[[nm]]
  }

  set.seed(seed)

  # Baseline
  W0 <- rnorm(n, 0, 1)
  if (is.null(drug)) {
    A0 <- rbinom(n, 1, 0.5)
  } else {
    A0 <- rep(as.integer(drug), n)
  }

  u_i  <- rnorm(n, 0, p$sigma_u)
  mu_i <- p$alpha0 + p$alphaA * A0 + p$alphaW * W0 + u_i

  L_prev   <- rep(NA_real_, n)
  PDC_prev <- rep(NA_real_, n)
  eta_prev <- rep(NA_real_, n)
  failed   <- rep(FALSE, n)
  Y_mat    <- matrix(0L, n, K)

  for (t in seq_len(K)) {
    # L_t
    noise_L <- rnorm(n, 0, p$sigma_L)
    if (t == 1) {
      L_t <- p$beta0 + p$betaW * W0 + p$betaA * A0 + noise_L
    } else {
      L_t <- p$beta0 + p$betaW * W0 + p$betaA * A0 +
        p$betaLag * L_prev + p$betaP * (1 - PDC_prev) + noise_L
    }

    # Natural PDC
    noise_e <- rnorm(n, 0, p$sigma_e)
    if (t == 1) {
      eta_t <- mu_i + noise_e
    } else {
      eta_t <- (1 - p$rho) * mu_i + p$rho * eta_prev +
        p$alphaL * L_prev + noise_e
    }
    nat_pdc <- plogis(eta_t)

    # Apply intervention if provided
    if (!is.null(intervene_pdc)) {
      PDC_t <- intervene_pdc(nat_pdc, t, A0, W0, L_t)
      PDC_t <- pmin(pmax(PDC_t, 0), 1)
    } else {
      PDC_t <- nat_pdc
    }

    # Failure hazard
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

    # Carry forward
    if (t > 1) {
      Y_mat[Y_mat[, t - 1] == 1L, t] <- 1L
    }

    L_prev   <- L_t
    PDC_prev <- PDC_t

    # Under an MTP, the "natural value at time t" arises from the
    # post-intervention history.  The AR(1) term feeds back the *realised*
    # (shifted) adherence on the latent scale so that future natural PDC
    # reflects the intervention that already occurred.
    # When no intervention is applied, PDC_t == nat_pdc so this is a no-op.
    if (!is.null(intervene_pdc)) {
      eta_prev <- qlogis(pmin(pmax(PDC_t, 1e-4), 1 - 1e-4))
    } else {
      eta_prev <- eta_t
    }
  }

  # Return cumulative failure by end of follow-up
  Y_final <- as.integer(Y_mat[, K] == 1L)
  list(Y_final = Y_final, risk = mean(Y_final), A0 = A0, W0 = W0)
}
