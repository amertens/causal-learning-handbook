# Polishing report: index.qmd

## Assessment
Ready to ship (with minor polish) - Short, well-organized preface; defines its three headline acronyms on first use and orients the reader across the six parts. Only small link-hygiene and consistency items remain.

## Strengths
- Clear "Who is this book for?" scoping (assumes R + regression, not potential outcomes) sets an honest prerequisite contract.
- Two-path "How to use this book" (Short Course vs front-to-back) plus per-part chapter tables give excellent navigational structure.

## Critical issues
None.

## Major issues
None. (Preface; substance lives in the chapters.)

## Minor issues
1. Mixed/possibly-broken link extensions - internal links point to `.Rmd` targets (`other_resources.Rmd`, `anotated_bibliography.Rmd`) while the rest use `.qmd`. Verify they resolve.
2. Misspelled filename in a user-facing link - `anotated_bibliography.Rmd` ("anotated", one n). Rename to `annotated_bibliography.qmd` or note intentionality.
3. DOI placeholder "DOI: *pending*" - leave a visible TODO so the Zenodo mint isn't forgotten.
4. Part-numbering vs file-numbering mismatch (03-07, 03-09, 03-09b, 03-03b) invites editing mistakes; consider a note that file prefixes are legacy and order is authoritative in `_quarto.yml`.

## Figures & tables
None.

## Citation gaps
None (preface defers method citations to chapters).

---

# Polishing report: 00-short-course.qmd

## Assessment
Minor-to-moderate polish - Pedagogically excellent (each method has a "why it exists" lead-in, a runnable estimate, a plain-language reading, and a worked exercise). Main issue: IPTW is promised as an estimator but never produces an estimate, colliding with the "watch the numbers land together" framing.

## Strengths
- The immortal-time-bias "teaching simplification" callout models intellectual honesty and motivates Module 3.
- Module 3's dual-role explanation of time-varying CD4 (confounder for the later decision, mediator of the earlier treatment) is one of the clearest statements of why g-methods are needed at this level.

## Critical issues
None. (g-computation formula, three assumptions, clever covariate H = A/g - (1-A)/(1-g), fluctuation-via-offset, and EIC-based SE all correct.)

## Major issues
1. IPTW billed as an estimator but never yields an estimate - the IPTW section fits the propensity model and plots overlap but computes no weighted ATE, yet the framing promises estimates that "land close together." Add a Hajek/stabilized IPTW estimate, or soften every "N ways" framing.
2. "Four ways / three ways" count is inconsistent between the Module 2 intro and the at-a-glance callout. Reconcile after fixing #1.
3. Headline methods (Roadmap, TMLE, Super Learner, MSM/IPTW, longitudinal TMLE, LMTP, E-values) introduced with no inline citations. Add one per method's first mention.

## Minor issues
1. Redundant glossary rows (propensity score twice; ATE and "risk difference" defined as the same quantity on two rows).
2. Inline figures (`fig-sc-dag`, `fig-sc-overlap`, `fig-sc-long-dag`, `fig-sc-positivity`) have labels but are never `@fig-`-referenced.
3. Exercise 1.1 `glimpse` solution names columns not yet shown to the reader; ensure names match.
4. Call the standardization identity "g-formula (Robins) / standardization formula" once.

## Figures & tables
- `fig-sc-overlap` and `fig-sc-positivity` visualize the same Module-2 propensity scores; note in Module 4 that it revisits the same `ps`, and cross-reference `@fig-sc-overlap`.
- `fig-sc-dag` - anchor with `@fig-sc-dag` in the Step 2 sentence.

## Citation gaps
- Roadmap: Petersen & van der Laan (2014); target trial: Hernan & Robins (2016).
- TMLE: van der Laan & Rubin (2006); van der Laan & Rose Targeted Learning (2011).
- Super Learner: van der Laan, Polley & Hubbard (2007).
- g-formula: Robins (1986); MSMs/IPTW: Robins, Hernan & Brumback (2000).
- Longitudinal TMLE: van der Laan & Gruber (2012).
- Modified treatment policies / lmtp: Diaz et al. (2023); Haneuse & Rotnitzky (2013).
- E-values: VanderWeele & Ding (2017); immortal time: Suissa (2008).
