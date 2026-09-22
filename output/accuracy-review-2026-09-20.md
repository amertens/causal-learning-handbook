# Accuracy review of the Causal Learning Handbook (2026-09-20)

Scope: every chapter listed in `_quarto.yml`, checked for (1) code that runs and output that matches the prose, (2) bibliography and citation accuracy, (3) statistical claims, and (4) whether the chapter is developed enough to sit on the live site. "Fixed" means the change is in the working tree as of this review; "Open" means a decision or work is needed.

Method: full static render (no errors); every executed output cell scanned for errors/warnings/NA/Inf; every `eval: false` chunk parsed and each called function checked against the installed packages; numeric prose claims compared with the rendered output; cohort facts re-derived from `ipw::haartdat`; bibliography entries checked against the literature from memory (flagged, not changed, where I could not be sure).

## 1. Findings that were wrong and are now fixed

| Chapter | Problem | Fix |
|---|---|---|
| 02-01 G-computation | The "watch it fail" demo did not show what it claimed. In its DGP a main-terms logistic model gives an essentially unbiased marginal risk difference (bias about 0.004 at n = 200,000); the 0.196 vs 0.215 gap shown was sampling noise, and the "correct" model landed at 0.199. | New DGP where treatment is concentrated in older patients and the benefit in younger ones (n = 20,000): main-terms bias is 15 to 30% of the effect, the correctly specified model is unbiased. Prose rewritten and now states honestly that outcome-model misspecification bias is usually moderate and grows with confounding-driven extrapolation. |
| 02-01 G-computation | "Prevented fraction" was computed with the population-attributable-fraction formula, giving a negative number that then needed a sign warning (the reviewer's item 6.c). | Now the actual prevented fraction, (P(Y(0)=1) − P(Y=1)) / P(Y(0)=1), positive for a protective exposure; callout explains when to use each formula. |
| 03-08 Clean room | Prose said the crude risk difference is "large and positive" and that the adjusted estimators "reverse the sign". On the cohort the crude RD is −0.010 and the adjusted estimates are ≈ −0.014; nothing reverses. | Paragraph rewritten with the inline computed value. |
| 03-08 Clean room | The `#>` outputs shown for `cleanTMLE` (287 events, 23.9% mortality, RD −0.092) cannot come from this cohort (31 deaths, 2.6%). And none of the 27 workflow functions the chapter calls (`attach_estimand`, `checkpoint_*`, `tmle_candidate`, `run_plasmode_feasibility`, `gate_all`, `authorize_outcome_analysis`, `run_iptw_workflow`, …) exist in the installed `cleanTMLE` 0.3.0, whose exports are `define_candidates`, `select_candidate`, `stress_test`, `estimate_effect`, `design_report`, `mask_outcome`/`unmask_outcome`, `negative_control_ladder`, … | Added a warning callout stating that the package code targets the v0.1 interface, is not executed, and that the shown output is illustrative. **Open:** the chapter needs rewriting against the current API or archiving (see §4). |
| 03-09 / 03-09b L-TMLE | The single simulated realization (n = 1,000, seed 2026) did not support the text: naive 0.128 vs g-computation 0.120 vs L-TMLE 0.122 vs truth 0.094, so all estimators looked equally off. A 200-replicate Monte Carlo confirms the claim is right (naive bias +0.025, g-computation bias −0.001); the draw was unlucky. | n raised to 5,000 in both chapters (same seed; at n = 5,000 naive = 0.124, g-computation = 0.093, truth = 0.094) and the Monte Carlo result stated in the prose. |
| 02-04 Super Learner | The expanded-library output shows `SL.ranger.mtry4_All` with cross-validated risk `Inf` and no explanation. | Annotation explains it (log-likelihood loss is infinite for a learner that predicts exactly 0/1 against the opposite outcome; the metalearner gives it zero weight). |
| 03-10 Mediation | Installation note omitted `mlr3extralearners`, which `crumble` requires for the `earth`/`ranger` learners (the fit errors without it). | Note updated. |
| Citations | `robins1994estimation` (Robins, Rotnitzky & Zhao on AIPW/missing data) was cited as the source for g-computation (02-01), IPTW/MSMs (02-02, 03-07) and described as "foundational paper on marginal structural models". `diaz2020lmtp` (a Biostatistics commentary on TMLE vs DML) was cited for the LMTP framework (03-03b). | Re-pointed to `robins1986gformula`, `robins2000msm`/`cole2008ipw`, and `diaz2021lmtp`. |
| book.bib | `wei2023simultaneous` had a placeholder arXiv id (`2306.00000`) and a made-up author list. | Replaced with Wei, Petersen, van der Laan, Zheng, Wu & Wang (2023), *Biometrics* 79(3):1934–1946, which is the simultaneous-subgroup TMLE paper the text describes. |
| book.bib | `smith2023systematic` combined the title of the TMLE systematic review with the author list and venue of a different Smith et al. paper (the 2022 *Statistics in Medicine* computational tutorial). | Corrected to Smith, Phillips, Luque-Fernandez & Maringe (2023), *Annals of Epidemiology* 86:34–48.e28. Please confirm the volume/pages. |
| book.bib | `vanderlaan2024atmle` and `vanderLaan2025atmle` were the same paper; the second key was unused. | Duplicate removed. |

## 2. Code: what runs, what does not, what is unverifiable

Everything with `eval: true` executes without errors or warnings in the render. The problems are all in the 84 `eval: false` chunks:

- **03-08** (12 of 18 chunks not executed): API mismatch above. Nothing in the package walkthrough can be run.
- **03-10** (14 of 36): `crumble` 0.1.2 and the `torch` backend are installed, so the fits could in principle run. I tried the chapter's recanting-twins fit on the chapter's own simulated data: one run returned all-zero estimates, a second failed inside torch (`Tensor_slice: Unsupported index`). The results table is labelled as saved output from a representative run; it cannot currently be reproduced. The `mma::weight_behavior` example was not tested.
- **03-06** (9 of 17): the `ES.cvtmle()` example runs as written (2 s with V = 3); it could be switched to `eval: true`. The `atmle` package is not installed and `nonparametric()`/`atmle()` argument names could not be checked.
- **03-04b** (6 of 9): the whole worked example is non-executed; `InterXshift` is not installed.
- **03-05** (33 of 51): template code; parses cleanly, but four templates use the `moments` package, which is not installed. Add `install.packages("moments")` or use base R for skewness/kurtosis.
- **03-03b**, **02-04**, **02-05**, **03-04**: non-executed chunks parse and call only existing functions.

## 3. Content claims to verify (not changed)

- 01-02: "In 2006, several observational studies reported that statin users had substantially lower mortality … cancer, trauma, and even infections." Suissa (2008) is cited but the year and the specific outcomes are not from that paper as far as I can tell. Suggest "In the mid-2000s" and dropping "trauma" unless you have a source.
- 03-04b: "A 2023 systematic review identified 89 observational studies that applied TMLE." I believe the Smith et al. review reports a different count (81). Check the abstract.
- 03-11: "continuous-time TMLE gave 7 to 29% efficiency gains" and "discrete-time TMLE with non-optimized bin widths sometimes did worse than the unadjusted estimators" are attributed to Chen et al. (2023, `concrete`); verify the numbers against the paper.
- 02-05: "Asymptotically, CV-TMLE is never worse than standard TMLE, and it is strictly better when the nuisance models rely on aggressive ML." Both are efficient when standard TMLE is valid; the honest claim is that CV-TMLE stays valid in cases where standard TMLE does not. Wording only.
- book.bib entries I could not verify: `talbot2025tmle` (marked "In press"), `basten2025additive`, `cashin2025target` (author list is "Cashin and others"; the TARGET guideline's first author is Hansford in the versions I know), `muntner2024staging`, `gatto2023spifd2`. `hernan2020whatif` and `hernan2025whatif` both cite *What If*; pick one.
- 03-14 cites `[@hernan2016target]` for the prescribing-preference instrument; that paper is about target trial emulation, not instrumental variables. Brookhart et al. (2006, *Epidemiology* 17(3):268–275, "Evaluating short-term drug effects using a physician-specific prescribing preference as an instrumental variable") is the standard reference.

## 4. Chapters to archive, trim, or hold back

Recommendation, in order of confidence:

1. **03-08 Staged and clean-room designs: hold back or rewrite.** The concept and the base-R stages are sound, but two thirds of the code addresses a package interface that no longer exists, and the sample outputs were not produced from the data on the page. Either rewrite the twelve chunks against `cleanTMLE` 0.3.0 and run them, or move the chapter to an "in development" appendix until then. The base-R diagnostics could be folded into 03-05 in the meantime.
2. **03-04b Interactions in RWE: trim by half or archive.** 2,377 lines, six of nine chunks not executed, an "Abstract" section, two Python sections, and heavy overlap with 03-04 on estimands, positivity in subgroups, multiplicity, and reporting. Its distinctive contributions (the four-cell joint-intervention estimand, the additive/multiplicative scale discussion, the preventive-exposure recoding, the Layer 2 `ltmle` sketch) would fit in 600 to 800 lines with one executed worked example. As it stands it reads as a manuscript pasted into the book.
3. **03-06 RWD meets RCT: keep as "in development" but make it honest.** Switch the two `ES.cvtmle()` chunks to `eval: true` (they run in seconds), drop or clearly caveat the `atmle` block (package not installable from CRAN, functions unverified), and remove the duplicate question sets ("Check Your Understanding" plus "Additional Questions" plus "Answers to Additional Questions").
4. **03-10 Mediation: keep, but decide about the table.** Either pin working `crumble`/`torch` versions and run the fits live, or present the saved-output table with the caveat moved up to the table caption. The conceptual material is good.
5. **03-09b L-TMLE by hand: keep.** Clearly marked optional, code is correct, and with n = 5,000 the demonstration now supports the text.
6. **04-03 Synthesis: keep, it is short.** It repeats the case study's framing but that is its job.
7. **02-04b Highly adaptive lasso** is already out of `_quarto.yml`; move the file to an `archive/` folder so it stops showing up in `_workshop-src/` builds, along with `03-11b`, `03-11c`, `03-12`, `03-13` leftovers.

Everything else (short course, case study, Parts I to III, positivity, diagnostics, applied LMTP, effect modification, time-to-event, identification beyond exchangeability, appendices) is executed, internally consistent, and I found no substantive errors beyond those listed above.
