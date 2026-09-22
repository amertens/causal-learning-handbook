# Paper Polishing Report: Causal Learning Handbook

*Constructive-advisor review of all 24 rendered chapters. Adapted from the single-paper paper-polishing skill to a pedagogical Quarto textbook: conference concerns (novelty/SOTA/ablations/baseline-fairness) were dropped; correctness of exposition, pedagogy, presentation, figures/tables, and scholarship were emphasized. R-code execution was NOT re-audited here (a separate accuracy audit already fixed those bugs); this pass targets prose, math exposition, structure, and citations. Full per-chapter reports were generated for every chapter; this document consolidates them.*

## Overall Assessment

**Recommendation: Minor Revision (book-wide).** The handbook is in strong shape. Across 24 chapters the modal verdict is "minor polish": the exposition is correct, well-motivated, and often better than published tutorials (the estimand-vs-model framing, the "watch it fail" simulations, and the honest treatment of untestable assumptions are recurring strengths). No chapter has a fatal correctness error in its prose. What holds the book back from "ready to ship" is a set of **systemic, easily-batched issues** (heading hierarchy, acronym expansion, figure/table cross-referencing, hardcoded section numbers, and missing foundational citations) plus a **short list of genuine per-chapter exposition errors** worth fixing before release. Two chapters need more than polish: **03-10 (mediation)** renders no estimates (its worked example is entirely `eval:false`), and **03-06 (hybrid designs)** contains an internally contradictory A-TMLE description.

## Cross-Cutting Themes (fix once, apply everywhere)

These recur in most chapters and are the highest-leverage fixes.

1. **Heading hierarchy — multiple H1s per chapter.** Many chapters use level-1 `#` for every section, so the chapter has several H1s and a broken sidebar/TOC nesting. Affected: 01-01, 02-01, 02-02, 02-03, 03-07, 03-09, 03-10, and others. *Two chapters' reviewers noted the multi-`#` pattern appears to be house style (03-03b, 03-09b).* **Decide once, book-wide:** either one H1 (title) + `##` sections everywhere, or confirm the current convention and make it uniform.

2. **Acronyms not expanded on first use.** Nearly universal: HAART, ATE, ATT, TMLE, AIPW, IPTW, IPCW, LMTP, MACE, RMST, PFS, RWE, CKD, eGFR, SAP, ESS, HAL, AUC, GEE, DAG, LATE, 2SLS. Expand each at first occurrence in each chapter (a reader may enter mid-book).

3. **Figures/tables labeled but not `@`-cross-referenced (or unlabeled).** Most chapters define `#| label:` figures but refer to them in prose as "the figure below," and several key tables have no `#| tbl-cap`/label at all. This wastes Quarto numbering and violates the handbook's own `@sec-`/`@fig-`/`@tbl-` convention (per project memory). Systemic.

4. **Hardcoded "Section N" / "Table N" references.** 02-03, 03-04, 03-04b, 03-05, 02-02 refer to sections/tables by manual number against unlabeled or manually numbered headings. These silently desync on reordering. Convert to `{#sec-}`/`{#tbl-}` labels + `@` refs.

5. **Missing foundational citations (the big one).** The same canonical references recur as gaps, and several keys are absent from `book.bib`. Highest-frequency: **Robins (1986)** g-formula; **Bang & Robins (2005)** sequential regression; **Robins, Hernán & Brumback (2000)** MSMs; **Cole & Hernán (2008)** IP weights/truncation; **Greenland, Robins & Pearl (1999)** non-collapsibility; **Suissa (2008)** immortal time; **Hernán & Robins (2016)** target trial; **VanderWeele & Ding (2017)** E-value; **Rothman (1976)** synergy/RERI; **Athey, Tibshirani & Wager (2019)** grf; **Wolpert (1992)/Breiman (1996)** stacking; **van der Laan & Dudoit (2003)** CV oracle; **Lipsitch et al. (2010)** negative controls; **Miao/Tchetgen Tchetgen (2018)** proximal; **Manski (1990)** bounds; **Aalen-Johansen (1978)** and **Fine & Gray (1999)**; **Colnet et al. (2024)**, **Dahabreh et al.**, **Bareinboim & Pearl (2016)** for data fusion/transport. See `missing_citations.json` for the per-location list. A focused "add missing .bib entries + inline cites" pass would close most of these.

6. **Redundancy in the long chapters.** 03-04b, 03-05, 03-06, 02-04, 02-03, 03-03b restate key points 3-6 times each. Targeted trims (stated below per chapter) would tighten without loss.

7. **Reading-time estimates run high** in several chapters (03-08, 03-09, 03-14, 03-04b, 04-03, 05-01). Sanity-check against actual length.

## Book-Wide Prioritized Checklist

**Critical / correctness of exposition (fix before release):**
- [ ] 03-04:86 — restore `\text{age<75}` in the CATE defining equation; literal TAB characters currently break the LaTeX (`ext{...}`).
- [ ] 03-08 Stage 1b — per-arm event counts (209|78 vs 824|376) let the reader recover the crude risk difference, defeating the chapter's outcome-blind guarantee; show pooled events + per-arm N only, or justify explicitly. Also fix the crude-effect **sign** (line 638 says treated have higher mortality; the illustrative counts imply lower).
- [ ] 03-06 @sec-atmle — A-TMLE is described as two incompatible methods (convex combination of estimators vs. pooled-minus-HAL-bias). Delete the `sum alpha_k Psi_k` framing; keep the bias-correction description.
- [ ] 03-10 assumption 6 — relabel "cross-world independence" as "no treatment-induced mediator-outcome confounding"; they are distinct. Also change "doubly robust" to "multiply robust" (two-of-three nuisance condition).
- [ ] 02-03:271 — the variance sentence contradicts its formula: `(1/n)Σ D_i^2` IS the EIF variance, not "EIF variance / n". State Var(psi_hat) = sigma^2/n separately (the code `sqrt(var(EIC)/n)` is correct).
- [ ] 02-05 — ES-CV-TMLE is the **Experiment-Selector** CV-TMLE (a data-fusion method, sibling of A-TMLE), not "an ensemble selector for the targeting step"; and "CV-TMLE removes the Donsker requirement" should not imply arbitrarily flexible ML is always valid (the nuisance rate condition remains).
- [ ] 04-03:37 — "g-methods ... modeling the treatment mechanism" is false for g-computation (it models the outcome and standardizes); reword. Also the "collider bias" claim (line 35) is overadjustment/mediator-blocking, not collider bias as set up.
- [ ] 03-05 — excess-kurtosis threshold (">3" should be ">0" for *excess*; reconcile with the ">5-6" callout and the raw `moments::kurtosis` output); replace the two non-rigorous displayed equations (the "·(terms)" DR identity at ~1103 and the bias formula at ~2206); delete the paragraph duplicated verbatim at ~2590 and ~2620.

**Major / clarity, pedagogy, scholarship:**
- [ ] 00-short-course — IPTW is billed as one of the estimators that should "land close together" but no IPTW estimate is ever computed; add a Hájek/stabilized estimate or soften the "N ways" framing (and reconcile "three ways" vs "four ways").
- [ ] 03-10 — render at least the simulated RT/RI fits (or a saved results table); the chapter currently promises a "headline number" and a three-way comparison but shows no estimate (all `crumble` chunks `eval:false`).
- [ ] 05-01 — distinguish dynamic treatment regime (the CD4<350 rule) from modified treatment policy (depends on the natural treatment value); fix the mis-attributed SDR/cross-fit citation (should be Díaz et al., not the Super Learner paper); either implement or explicitly scope-out the censoring node the prose claims to handle.
- [ ] 03-07 — pick one canonical estimand (always-vs-never used in the demo vs. the vs-observed estimand in the notation section); unify the intermediate confounder's name (L1 in some DAGs, L2 in others); align the intro DAG with the DGP.
- [ ] 03-04b — reconcile the `Anodes = c("A","E")` prose with the working `c("E","A")` code; state the temporal-ordering rule once; trim ~25% redundancy.
- [ ] 03-03b:237 — the drug's *direct* adherence effect (`alphaA = -0.35`) means the adherence gap is not pure confounding; reword. Fix the notation clash in the forgiveness-curve equation (uses A-bar/D vs. the chapter's PDC/A0).
- [ ] 01-02 — the roadmap "mapping" table contradicts its own caption (rows are shifted, so stages misalign across formulations); re-align or reword.
- [ ] 02-04 — soften the finite-sample "SL is at least as good as the best learner" to asymptotic, and note the honest comparison is CV.SuperLearner (the in-sample `Z %*% alpha` comparison is guaranteed to favor SL); fix the Brier-score range ("0 to 0.25" — 0.25 is the trivial-predictor score, not the max); merge the duplicated loss-function sections.
- [ ] 01-01 — "cumulative mortality risk over follow-up" overstates what the logistic g-computation computes (it ignores time-at-risk/censoring); interpolate the real value where `$-X$ per 100` placeholders currently render literally; deliver or drop the "collider bias" learning objective.
- [ ] 03-11 — add primary citations for Aalen-Johansen, Fine & Gray, survtmle, and separable effects (named but uncited; keys missing from `book.bib`).
- [ ] 03-14 — honor the scope note's promise of DiD/RDD/ITS further reading (or reword); add IV and proximal DAGs; add at least a bounds pointer (Manski).
- [ ] 03-06 — separate data-fusion-for-precision (ES-CV-TMLE, A-TMLE target the trial ATE) from transportability-for-generalization (PATE); unify the `@vanderlaan2024atmle`/`@vanderLaan2025atmle` key.

**Minor (batchable):** acronym expansions; `@fig-`/`@tbl-` references and missing captions; heading demotions; reading-time sanity checks; redundancy trims. Details in the per-chapter summaries.

## Per-Chapter Summaries

Each entry: verdict, then the critical/major items. (Minor, figure, and citation detail is captured in the themes above and `missing_citations.json`. The full per-chapter report for `index`/`00-short-course` is saved in `reports/00-front-matter.md` as a depth sample; the rest are summarized here.)

- **index.qmd** — *Ready (minor).* Verify `.Rmd` vs `.qmd` link targets resolve; the `anotated_bibliography` filename reads as a typo; DOI placeholder pending.
- **00-short-course.qmd** — *Minor-to-moderate.* IPTW promised but never estimated ("three/four ways" inconsistent); add method citations; duplicate glossary rows.
- **00-case-study.qmd** — *Minor.* Add a formal estimand equation and fix a concrete time horizon (τ) for the running example; state the longitudinal identification assumptions (esp. positivity for the dynamic rule); reconcile "during interval" vs "at end of interval" ordering; unify the two citation mechanisms.
- **01-01-regression-vs-causal.qmd** — *Minor.* "Cumulative mortality risk" overstates the logistic g-comp; `$-X$` placeholders render literally; cross-reference the figure/table; deliver-or-drop the collider-bias objective; fix "IPTW" expansion (of *treatment*).
- **01-02-causal-roadmap.qmd** — *Minor.* Roadmap mapping table contradicts its caption; cross-reference `fig-roadmap-dag`; remove the thrice-stated causal-gap passage; the risk-difference-vs-baseline-prevalence claim is loose.
- **02-01-gcomputation.qmd** — *Minor.* State the three identification assumptions up front (currently only a parenthetical); reference `fig-gcomp-counterfactual`; add Robins (1986); trim the PAF sign digression.
- **02-02-iptw.qmd** — *Minor.* Hardcoded "Section 8" ref; heading hierarchy; add MSM (Robins-Hernán-Brumback 2000) and Cole & Hernán (2008); balancing-score sentence conflates balance with unconfoundedness.
- **02-03-doubly-robust-tmle.qmd** — *Minor.* **Variance equation off by n (line 271, correctness).** Hardcoded "Section 7"; heading hierarchy; efficiency vs. targeting conflated in the IF callout; add efficiency-bound and DR-origin citations.
- **02-04-superlearner.qmd** — *Minor.* Finite-sample "at least as good" overstated + in-sample comparison favors SL by construction; Brier range wrong; duplicated loss-function sections; convex-constraint under-specified; add stacking (Wolpert/Breiman) and CV-oracle (van der Laan & Dudoit) citations.
- **02-05-advanced-tmle.qmd** — *Minor.* **ES-CV-TMLE misdescribed; CV-TMLE "removes all conditions" overstated (correctness).** Relocate the AIPW cross-check; note cross-fit AIPW = DML.
- **03-02-positivity-overlap.qmd** — *Minor.* IPTW estimator mislabeled "Horvitz-Thompson" (it is Hájek); add Westreich & Cole (2010), Crump et al. (2009), Cole & Hernán (2008); ESS cross-ref points away from this chapter.
- **03-07-longitudinal-td-confounding.qmd** — *Minor.* Two estimands both called canonical; L1/L2 naming inconsistent across DAGs; intro DAG contradicts the DGP; add Robins (1986), Hernán-Brumback-Robins (2000), Daniel et al. (2013).
- **03-09-illustrated-ltmle.qmd** — *Minor.* "Static vs dynamic regimes" set up but never delivered; "cohort" reads as real data but is simulated; add Bang & Robins (2005) for the backwards-regression core.
- **03-09b-ltmle-by-hand.qmd** — *Minor.* The efficient influence curve (the payoff) lives only in code comments — promote it to a displayed equation and define it; cite Bang & Robins and van der Laan & Gruber; explain *why* the naive/g-comp estimate is biased (standardizes over observed L1).
- **03-03b-adherence-shift-methods.qmd** — *Minor.* Cohort interpretation contradicts the DGP's direct drug-on-adherence effect (`alphaA`); notation clash in the central equation; DAG omits `A0 -> PDC_t` edges; expand LMTP; add Díaz & van der Laan (2012), Young et al.
- **03-05-advanced-diagnostics.qmd** — *Minor (but several correctness nits).* Excess-kurtosis threshold; two non-rigorous displayed equations; verbatim-duplicated paragraph; ~30 hardcoded section refs; flagship methods cited only in back-matter; consider splitting the ~2800-line chapter.
- **03-08-staged-clean-room.qmd** — *Minor.* **Per-arm event counts leak the outcome-blind association (central-concept, correctness); crude-effect sign contradicts the illustrative numbers.** Add Hernán & Robins (2016); render a static bias-severity figure; `LTMPS` typo (→ LMTP); expand SAP/ESS/IPCW.
- **05-01-applied-lmtp-workflow.qmd** — *Minor.* Dynamic regime vs. MTP conflation; SDR/cross-fit citation mis-attributed; censoring claimed but not implemented; cross-reference the dose-response figure; unit notation (cells/µL vs cells/mm³).
- **03-04-effect-modification.qmd** — *Minor.* **Broken `\text{}` in the CATE definition (line 86, correctness).** Notation drift (V/X/X_g/W; Y^1 vs Y(1)); te-VIM citation mis-keyed; add grf (Athey-Tibshirani-Wager 2019) and metalearners (Künzel et al.); hardcoded section refs.
- **03-04b-interactions-rwe.qmd** — *Minor.* Length/redundancy (~2370 lines, trim ~25%); Anode-ordering prose vs. code inconsistency; add Rothman (1976); modifier notation V vs M; hardcoded section/table numbers.
- **03-10-tmle-mediation.qmd** — *Major polish.* **Worked example renders no estimate (all `eval:false`); cross-world-independence mislabeled; "doubly" vs "multiply" robust.** Heading hierarchy; add Robins & Greenland (1992), Pearl (2001), VanderWeele-Vansteelandt-Robins (2014).
- **03-11-time-to-event-competing-risks.qmd** — *Minor.* Foundational estimators (Aalen-Johansen, Fine-Gray, survtmle, separable effects) named but uncited (keys missing from bib); RMST/all-cause-vs-CIF notation nits; expand MACE/RMST/PFS/AIPW.
- **03-14-identification-beyond-exchangeability.qmd** — *Minor.* Scope note promises DiD/RDD/ITS reading it never gives; no bounds coverage; zero figures (add IV and proximal DAGs); define LATE; thin inline citations.
- **04-03-case-study-synthesis.qmd** — *Minor.* **"g-computation models the treatment mechanism" is false; collider-bias claim unsupported (correctness).** Inconsistent Sources list; add MSM/g-formula/E-value citations.

## How to Act on This

A high-leverage sequence: (1) a single **book-wide mechanical pass** for the five systemic themes (headings, acronyms, `@`-refs, hardcoded numbers, reading times); (2) a **`book.bib` + inline-citation pass** driven by `missing_citations.json`; (3) the **per-chapter correctness fixes** in the Critical checklist; (4) decide on 03-10's rendering and 03-06's A-TMLE rewrite. Steps 1-2 alone would move most chapters from "minor polish" to "ready."
