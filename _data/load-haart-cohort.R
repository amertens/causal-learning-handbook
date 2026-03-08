# Shared Data Loader: HIV Antiretroviral Therapy Cohort
# -----------------------------------------------------
# Source this file in any chapter that needs the running case study:
#   source("_data/load-haart-cohort.R")
#
# The dataset is the longitudinal HIV antiretroviral therapy cohort
# used in the marginal structural model literature (Van der Wal & Geskus, 2011).
# It contains 1,200 simulated HIV-positive patients observed in 100-day
# intervals from HIV seroconversion, with time-varying CD4 counts,
# HAART initiation, mortality, and dropout.
#
# Variables:
#   patient   - patient ID
#   tstart    - interval start (days since seroconversion)
#   fuptime   - interval end (days since seroconversion)
#   haartind  - HAART initiated at interval end (0/1)
#   event     - death at interval end (0/1)
#   sex       - sex (0 = male, 1 = female)
#   age       - age at follow-up start (years)
#   cd4.sqrt  - square root of CD4 count at interval end
#   endtime   - patient's final observed time
#   dropout   - censored at interval end (0/1)

# Load from bundled .rda file
load("_data/haartdat.rda")

# Also create a point-treatment (baseline) snapshot for chapters
# that need a simpler cross-sectional version of the data
haart_baseline <- haartdat |>
  dplyr::group_by(patient) |>
  dplyr::summarise(
    sex = dplyr::first(sex),
    age = dplyr::first(age),
    cd4_baseline = dplyr::first(cd4.sqrt)^2,  # back-transform to CD4 count
    cd4_sqrt_baseline = dplyr::first(cd4.sqrt),
    ever_haart = max(haartind),
    died = max(event),
    dropout = max(dropout),
    total_fup_days = max(fuptime),
    .groups = "drop"
  )
