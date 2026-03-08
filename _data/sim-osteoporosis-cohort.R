# Shared Simulated Dataset: Osteoporosis Treatment Cohort
# -------------------------------------------------------
# This script generates the simulated cohort used throughout the handbook.
# Source this file in any chapter that needs the data:
#   source("_data/sim-osteoporosis-cohort.R")
#
# The dataset simulates a claims-database cohort comparing denosumab (A=1)
# vs zoledronic acid (A=0) for 3-year MI risk among osteoporosis patients.

set.seed(2025)
n <- 4000

# Baseline covariates
age <- rnorm(n, 75, 6)
cvd <- rbinom(n, 1, plogis(0.12 * (age - 70)))

# Treatment assignment (confounded by age and CVD history)
A <- rbinom(n, 1, plogis(-1 + 0.09 * (age - 70) + 1.5 * cvd))

# Binary outcome: 3-year MI (causal risk difference ~ 0.08)
Y <- rbinom(n, 1, plogis(-2 + 0.5 * A + 0.10 * (age - 70) + 1.0 * cvd))

# Time-to-event outcome (for survival analyses)
baseline_hazard <- 0.05
log_hr <- 0.4 * A + 0.05 * (age - 70) + 1 * cvd
event_time <- rexp(n, rate = baseline_hazard * exp(log_hr))
censor_time <- 3
time <- pmin(event_time, censor_time)
event <- as.integer(event_time <= censor_time)

# Assemble the analysis dataset
dat <- data.frame(
  age       = age,
  cvd       = cvd,
  A         = A,
  Y         = Y,
  time      = time,
  event     = event
)

# Clean up intermediate objects
rm(baseline_hazard, log_hr, event_time, censor_time)
