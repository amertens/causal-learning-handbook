# Simulation V4 — Simplified Longitudinal Adherence Simulation

## What changed relative to V3

V3 used a complex EHR-inspired DGP with many baseline covariates (age, sex, race, comorbidity, CD4, etc.), multiple time-varying confounders (CD4, symptoms, depression, substance use), informative monitoring, treatment switching, and a two-arm oral-vs-injectable comparison with pharmacokinetic decay modelling.

**V4 is intentionally simpler.** It strips down to the minimum structure needed to demonstrate the core methodological lessons:

| Feature | V3 | V4 |
|---------|----|----|
| Baseline covariates | 6+ | 1 (`W0`) |
| Time-varying confounders | 2-4 per block | 1 (`L_t`) |
| Treatment | Oral PDC + LA timeliness | Continuous PDC only |
| Regimen comparison | Oral vs. injectable | Drug A vs. Drug B (both oral) |
| Monitoring model | Informative | None |
| Switching | Yes | No |
| Pharmacokinetics | Exponential decay for LA | None |
| Default sample size | 2000-5000 | 1500-2000 |
| Default MC reps | 200 | 50 |
| LMTP learners | SuperLearner | `SL.glm` |

## Why V4 is intentionally simpler

The purpose is pedagogical. V4 is designed so that:

1. **GLM-based LTMLE works well** in the clean setting — because the DGP is low-dimensional and the treatment/confounder relationships are smooth and parametric.
2. **Positivity problems emerge naturally** for extreme static regimes (always-low adherence) — demonstrating a real limitation of categorical LTMLE.
3. **LMTP shift interventions remain feasible** — because they modify the observed treatment rather than forcing implausible static values.
4. **The code is readable** — each file is self-contained with clear variable names and comments.

## Why 90-day PDC is generated directly

V3 simulated monthly adherence and aggregated to 90-day blocks. V4 generates PDC directly at the 90-day block level. This avoids unnecessary intermediate steps and makes the causal structure transparent: one treatment node per block, one confounder per block, one outcome per block.

## Why one baseline and one time-varying confounder

The core methodological lessons do not require high-dimensional confounding. One baseline confounder (`W0`, representing clinical severity) and one time-varying confounder (`L_t`, representing current health status) are sufficient to create:

- **Time-dependent confounding**: `PDC_{t-1}` → `L_t` → `PDC_t` → `Y_t`
- **Confounding bias in naïve analyses**
- **A setting where GLM-based LTMLE can succeed**

Adding more covariates would slow execution and obscure the key lessons without changing the qualitative conclusions.

## Main estimators

### 1. Naïve adherence-stratified analysis
Stratify by observed mean PDC and compare failure rates. This is biased because it conditions on a post-treatment variable affected by confounders.

### 2. LTMLE on categorised adherence
Uses the `ltmle` package with 3-level categorical treatment (low/mid/high PDC) encoded as two binary dummy nodes per block. Static regimes like "always high" and "always low" are compared. GLM nuisance models are used by default.

### 3. LMTP on continuous PDC
Uses the `lmtp` package with continuous PDC as the treatment. Two intervention types:
- **Additive shift**: increase each person's PDC by +10%, capped at 1.0
- **Quantile-matched drug comparison**: map one drug's adherence distribution to another's via empirical quantile matching

## Scientific lessons

1. **Naïve adherence-stratified analyses are biased** due to time-dependent confounding (sicker patients have both lower adherence and higher failure risk).

2. **LTMLE recovers approximately unbiased estimates** for static adherence regimes when the DGP is low-dimensional and GLM models are well-specified.

3. **Static regimes create positivity problems**: "always low adherence" is rarely observed in practice, especially over multiple blocks. LTMLE estimates for these regimes become unstable.

4. **LMTP avoids the worst positivity problems** because modified treatment policies depend on the natural treatment value and only shift it by a feasible amount.

5. **Drug comparisons require matching adherence distributions** to isolate the pharmacological effect from differential adherence patterns.

## Why GLM-based LTMLE succeeds in the clean setting

The DGP has:
- **One continuous baseline covariate** and **one continuous time-varying confounder** — easily captured by linear models
- **Smooth logistic relationships** between covariates, treatment, and outcome
- **Adequate support** for common regimes (always-high, always-mid) at moderate sample sizes
- **Treatment categories broad enough** (3 levels) that each has reasonable mass

This means a correctly specified GLM can approximate the true nuisance functions well, and the TMLE targeting step corrects remaining bias.

## Why LMTP is appealing here

Modified treatment policies (MTPs) estimate the effect of *realistic modifications* of the observed treatment process. Key advantages:

- **Feasibility**: a +10% PDC shift is something a realistic intervention might achieve, whereas "always maintain PDC > 0.85 for 12 months" may be implausible for many patients.
- **Support preservation**: the shifted treatment value remains in the support of the observed treatment distribution (by construction), avoiding the positivity violations that plague static regimes.
- **Continuous treatment**: LMTP works directly with continuous PDC, avoiding the information loss from categorisation.
- **Drug comparison**: quantile matching creates a fair comparison by equalising adherence distributions, isolating the pharmacological drug effect.

The `lmtp` package supports all of these interventions in longitudinal wide data with `mtp = TRUE`.

## Scripts — run order

```bash
# From the simulation V4/ directory:

# 1. Calibrate DGP parameters (tune and inspect)
Rscript scripts/00_calibrate_v4.R

# 2. Single-run demo of all estimators
Rscript scripts/01_single_run_examples_v4.R

# 3. Monte Carlo comparison (default B=50, fast)
Rscript scripts/02_mc_compare_estimators_v4.R

# 3b. Monte Carlo with more reps
Rscript scripts/02_mc_compare_estimators_v4.R 200
```

## Dependencies

- R (>= 4.0)
- `ltmle` — for categorical LTMLE
- `lmtp` — for continuous LMTP
- `SuperLearner` — required by both (but we use `SL.glm` only)

Install with:
```r
install.packages(c("ltmle", "lmtp", "SuperLearner"))
```

## Monte Carlo plan (symposium-targeted)

Three phases of MC simulation, each producing a compact CSV and one slide-ready figure:

| Phase | Script | Purpose | Default |
|-------|--------|---------|---------|
| 1 | `03_mc_phase1_validation.R` | Validate 1- and 4-block LTMLE + LMTP | n=1500, B=100 |
| 2 | `04_mc_phase2_comparison.R` | Slide-ready bias comparison (naive vs LTMLE vs LMTP) | n=2000, B=300 |
| 3 | `05_mc_phase3_positivity_stress.R` | Show LMTP stability under positivity stress | n=2000, B=200 x 3 scenarios |

```bash
# Phase 1: quick validation
Rscript scripts/03_mc_phase1_validation.R 100

# Phase 2: method comparison (slower)
Rscript scripts/04_mc_phase2_comparison.R 300

# Phase 3: positivity stress test
Rscript scripts/05_mc_phase3_positivity_stress.R 200
```

## File structure

```
simulation V4/
├── R/
│   ├── sim_v4_dgp.R            # Core DGP + defaults + interventional simulator
│   ├── oracle_v4.R             # Oracle truth via large MC (uniform + midpoint)
│   ├── prepare_ltmle_v4.R      # LTMLE data prep, fitting, explicit Qform/gform
│   ├── prepare_lmtp_v4.R       # LMTP data prep, shift/matching, safe extraction
│   └── diagnostics_v4.R        # Targeted figures + tables for slides
├── scripts/
│   ├── 00_calibrate_v4.R       # DGP calibration & diagnostics
│   ├── 01_single_run_examples_v4.R  # End-to-end demo
│   ├── 02_mc_compare_estimators_v4.R  # General Monte Carlo study
│   ├── 03_mc_phase1_validation.R     # Phase 1: quick validation
│   ├── 04_mc_phase2_comparison.R     # Phase 2: slide-ready comparison
│   └── 05_mc_phase3_positivity_stress.R  # Phase 3: positivity stress test
├── manuscript/
│   ├── v4_descriptive_report.qmd     # Single-dataset descriptive summary
│   ├── ltmle_debug_v4.qmd           # LTMLE bias diagnosis
│   └── talk_results_companion_v4.qmd # Presentation figure/table bank
├── slides/
│   └── forgiveness_talk_v4.qmd       # RevealJS symposium deck
└── README.md
```
