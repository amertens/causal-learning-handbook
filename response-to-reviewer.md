---
title: "Response to reviewer (Carrie Nielson) — Causal Learning Handbook drafts"
date: 2026-05-20
---

This document responds to the 26 comments left by Carrie Nielson on the docx drafts dated 2026-05-07 / 2026-05-11. Comments are organized by source file in the order they appear. Each entry shows the anchored text, the reviewer's comment, and how it was addressed. Items marked **FLAGGED** are not yet implemented and need a decision before we make the change.

---

## File: `000-ReadMe.docx` → `index.qmd`

### 1. "lmtp" — define or add a glossary
- **Anchor:** "…longitudinal causal inference, and lmtp in a single chapter."
- **Comment:** Define here or add a "glossary" of abbreviations and terms.
- **Resolution:** Expanded inline to "the `lmtp` package (**L**ongitudinal **M**odified **T**reatment **P**olicies)" on first mention. I also expanded the related abbreviations encountered later in the same paragraph and in the Part-by-Part tables (L-TMLE, LMTP, AIPW, HTE, CATE, IOSW, DGP — see items below).
- **Resolved (added 2026-05-20, second pass).** Created [`glossary.qmd`](glossary.qmd), grouped by *Causal framework and study design*, *Identification assumptions*, *Estimation*, and *Clinical / pharmacoepidemiology*. Registered as the first appendix in [`_quarto.yml`](_quarto.yml) and linked from the appendix table in [`index.qmd`](index.qmd). The "Who is this book for?" section also points readers at the glossary.

### 2. "regression is not enough"
- **Anchor:** "Part I builds the conceptual foundations: why regression is not enough…"
- **Comment:** Clarify — why it's not enough for what purpose?
- **Resolution:** Edited to "why standard regression is not enough **for causal questions**".

### 3. "L-" (L-TMLE)
- **Anchor:** "…time-varying confounding, L-TMLE, and LMTP."
- **Comment:** Define first.
- **Resolution:** Edited to "longitudinal TMLE (L-TMLE), and Longitudinal Modified Treatment Policies (LMTP)".

### 4. "table-2 fallacy" — plainer language
- **Anchor:** "Why regression coefficients are not causal effects, the counterfactual framework, and the table-2 fallacy."
- **Comment:** Put this in simpler terms here?
- **Resolution:** Rewrote the row description to: "Why a regression coefficient for a treatment is not, on its own, a causal effect; the counterfactual ('what if everyone had been treated?') framework; and the 'Table 2 fallacy' — interpreting every adjusted coefficient in a regression output as a causal effect for its own variable."

### 5. "predict" — predict what?
- **Anchor:** "Outcome modeling and standardization: predict under interventions and average."
- **Comment:** Predict what?
- **Resolution:** Rewrote to: "fit a model for the outcome, predict each person's outcome under each treatment value, and average across the sample."

### 6. "reweighting" — why *re*-weighting?
- **Anchor:** "Propensity-score reweighting, stabilized weights, truncation, and diagnostics."
- **Comment:** Why *re*weighting?
- **Resolution:** Rewrote the row to drop the "re-" and to spell out what IPTW does up front: "Inverse-probability-of-treatment weighting: weight each individual by the inverse of their probability of receiving the treatment they actually got, so the weighted sample mimics one in which treatment was randomized. Covers stabilized weights, truncation, and diagnostics."

### 7. "AIPW" — define first
- **Anchor:** "AIPW and TMLE: combining outcome and treatment models for double robustness."
- **Comment:** Define first or link to glossary.
- **Resolution:** Expanded to "Augmented Inverse Probability Weighting (AIPW) and Targeted Maximum Likelihood Estimation (TMLE)".

### 8. "breaks" — softer phrasing
- **Anchor:** "Why treatment-confounder feedback breaks standard regression and motivates g-methods."
- **Comment:** "is not addressed by"?
- **Resolution:** Adopted reviewer's wording: "Why treatment–confounder feedback **is not addressed by** standard regression and motivates g-methods."

### 9. "distribution-shift" — distribution of what?
- **Anchor:** "Comparing drugs under realistic adherence with distribution-shift interventions."
- **Comment:** Distribution of what?
- **Resolution:** Rewrote to specify: "interventions that shift the *distribution of adherence* (rather than setting adherence to a fixed value)."

### 10. "DGP" — define first
- **Anchor:** "Forgiveness curves, the adherence DGP, and four ways to shift an adherence distribution."
- **Comment:** Define first or link to glossary.
- **Resolution:** Expanded to "the adherence **data-generating process (DGP)**". I also merged this row with the duplicate-targeting adherence row immediately above it — both pointed at `03-03b-adherence-shift-methods.qmd`, which was clearly a bug.

### 11. "HTE" — define first
- **Anchor:** "Effect Modification and HTE" (chapter title row).
- **Comment:** Define HTE first or link to glossary.
- **Resolution:** Expanded chapter-link text to "Effect Modification and Heterogeneous Treatment Effects (HTE)".

### 12. "CATE" — define
- **Anchor:** "CATE estimation, TMLE-based subgroup analyses, and causal forests."
- **Comment:** Define.
- **Resolution:** Expanded to "Estimating the conditional average treatment effect (CATE)…".

### 13. "IOSW" — define
- **Anchor:** "Trial-sample vs target-population effects and IOSW/TMLE transportability methods."
- **Comment:** Define.
- **Resolution:** Rewrote to: "transportability methods based on inverse odds of sampling weights (IOSW) and TMLE."

### 14. "Pre-specified analysis for regulatory pharmacoepidemiology" — reframe
- **Anchor:** Chapter title row in Part VI of the drafted README.
- **Comment:** Reframe as guarding against investigator bias and improving transparency. Regulatory studies (and most rigorous studies) are prespecified. Clean room adds the formal staging and review process to help guard against accusations of data dredging and selective reporting.
- **Resolved (added 2026-05-20, second pass).** Two edits:
    1. **Primary home** — added a new subsection *Pre-specification and the clean-room workflow* to [`01-02-causal-roadmap.qmd`](01-02-causal-roadmap.qmd), positioned between Step 7 (Interpretation) and the Running Example. This is the natural home because the chapter already has a learning objective committed to it ("Explain why pre-specifying an analysis plan before data access reduces researcher degrees of freedom") that was previously not delivered on. The new subsection adopts Carrie's framing essentially verbatim: pre-specification "guards against investigator bias and improves transparency" and "well-designed observational studies (and essentially all confirmatory clinical trials) have always been pre-specified"; the clean room is presented as the *formal staging and external review* that builds on pre-specification, "guarding regulatory and high-stakes RWE analyses against accusations of data dredging and selective reporting." A short, practical paragraph translates the principle for readers doing non-regulatory work (write a dated, version-controlled analysis plan; lock it before pulling outcomes). I also added a one-clause forward reference from the outcome-blind simulation bullet in Step 6 so readers find the new subsection.
    2. **Forward-reference fix** — the closing sentence at [`03-06-rwd-meets-rct-hybrid-designs.qmd:591`](03-06-rwd-meets-rct-hybrid-designs.qmd:591) used to promise a "next chapter" on TMLE clean rooms (the deleted 03-08). Rewrote it to absorb Carrie's framing and point back to the new 1.2 subsection: hybrid designs raise the stakes for analytic discipline; the clean-room workflow guards them against data dredging and selective reporting; regulatory submissions formalize this through staged review. No new chapter promised.
- I did **not** edit [`03-05a-assumption-diagnostics.qmd:497`](03-05a-assumption-diagnostics.qmd) — its existing "credibility established before the results are known" framing is already very close to Carrie's and a third re-statement would be redundant.

### 15. "longitudinal methods are needed" — needed for what?
- **Anchor:** "What naive analysis misses, how the Causal Roadmap clarifies the question, and why longitudinal methods are needed."
- **Comment:** Needed for what?
- **Resolution:** Rewrote to: "…and why longitudinal **g-methods are needed to recover an unbiased estimate of the HAART-on-mortality effect under treatment–confounder feedback**."

### 16. "How to use this book" — add a "who is this book for?" section
- **Anchor:** "0.1 How to use this book"
- **Comment:** Add a section on "who is this book for?" (statistician, pharmacoepi, journal / regulatory reviewer?)
- **Resolved (added 2026-05-20, second pass).** Inserted a "Who is this book for?" section in [`index.qmd`](index.qmd) immediately after the author/contact paragraph and before "How to use this book". The inserted text commits to the audience laid out in the draft (primary: pharmacoepidemiologists, biostatisticians, epidemiologists; secondary: reviewers and methods-oriented graduate students), states that prior comfort with R and basic regression is assumed but potential-outcomes notation and g-methods are not, points readers at the new Glossary appendix, and recommends Hernán & Robins, *Causal Inference: What If* as the foundational text. If the audience claim is wrong (e.g., you'd rather lead with "biostatisticians" or include health economists), say which audience to add/remove and I'll adjust — but the text is now in the book.

---

## File: `00-case-study.docx` → `00-case-study.qmd`

### 17. Opening — frame by question, not cohort
- **Anchor:** "…a longitudinal HIV antiretroviral therapy cohort from the marginal structural model literature…"
- **Comment:** I'd frame the use case as the causal question, not the cohort. So the use case is "an investigation of whether HAART initiation affects mortality in people with HIV" or similar?
- **Resolution:** Rewrote the opening paragraph to lead with the bolded causal question ("does initiating HAART reduce mortality among people living with HIV, and if so, by how much?") and present the cohort as the dataset used to answer it.

### 18. "initiating" — initiation vs use
- **Anchor:** "…does **initiating** HAART reduce mortality among HIV-positive patients…"
- **Comment:** Initiation is used here, but then time-dependent confounding is described. So it seems the original question was HAART use?
- **Resolution:** Added a clarifying sentence right after the question: "Although this sounds like a question about HAART 'use,' it is really a question about *strategies for when to start* HAART, because therapy is initiated (and then sustained) in response to evolving clinical status. The longitudinal, time-dependent nature of the initiation decision is what creates the confounding structure we have to address." This explicitly reconciles "initiation" with the subsequent treatment-confounder-feedback discussion.

### 19. "conflates the drug's benefit" — expand
- **Anchor:** "…conflates the drug's benefit with baseline severity."
- **Comment:** Perhaps expand this. A naïve comparison makes HAART patient outcomes appear worse than ART(?) simply due to the severity of HAART initiators' infection at baseline.
- **Resolution:** Expanded with: "In this dataset that means HAART initiators look *worse* than non-initiators in a crude analysis, simply because they had lower CD4 counts and worse prognosis at the moment they started therapy." (I read the reviewer's "ART(?)" as a transcription glitch — the live contrast in this cohort is HAART-initiators vs HAART-non-initiators, not HAART vs older mono/dual ART; let me know if I misread.)

### 20. "vague answers" — interpretability framing
- **Anchor:** "…vague causal questions produce vague answers."
- **Comment:** I'd say they can produce precise answers but can be uninterpretable relative to the original causal question.
- **Resolution:** Replaced with: "…vague causal questions produce numerically precise answers that cannot be interpreted relative to any well-defined contrast."

### 21. "haartdat" — prospective or simulated?
- **Anchor:** "…patients observed from HIV seroconversion in person-period (long) format, with one row per patient per 100-day interval."
- **Comment:** Was this a prospective cohort study? (The interval seems too standard to be from RWD.)
- **Resolution:** Added an explicit statement that the dataset is simulated: "It is a *simulated* cohort (not a prospective study) built to reproduce the structure of real HIV cohorts used in the marginal structural model literature." Also added a note that the regular 100-day grid is a simulation artefact and that real RWD would have irregular visit times that themselves need modelling.

### 22. "summary" → summary of baseline characteristics
- **Anchor:** Section heading "3.3 Cohort summary"
- **Comment:** of baseline characteristics
- **Resolution:** Renamed heading to "Cohort summary of baseline characteristics".

### 23. "chapters that use point-treatment examples" — confusing phrasing
- **Anchor:** "We also create a baseline snapshot (one row per patient) for chapters that use point-treatment examples:"
- **Comment:** This phrase is a bit confusing. I guess it means in some chapters, we will only consider baseline treatment initiation?
- **Resolution:** Rewrote to: "Several chapters (g-computation, IPTW, point-treatment TMLE) work with a simpler point-treatment version of the question, contrasting *ever* versus *never* initiating HAART using only baseline confounders. For those chapters we collapse the longitudinal data into a baseline snapshot with one row per patient."

### 24. "when clinically indicated" — what triggers initiation?
- **Anchor:** "…a strategy of initiating HAART when clinically indicated, compared with a strategy of deferring HAART"
- **Comment:** So each time the CD4 count drops below a certain level, HAART initiation is indicated? Or are there other indicators?
- **Resolution (partial — see flag).** Added a clarifying paragraph after the question: "'Clinically indicated' here means a CD4-driven initiation rule consistent with the guidance in force during the cohort era — for concreteness, *initiate HAART at the first interval at which CD4 falls below a fixed threshold (e.g., 350 cells/µL)* and continue thereafter."
- **FLAGGED — threshold choice.** I picked "350 cells/µL" as a representative historical guideline value, hedged with "e.g.,". The `haartdat` simulation is described in Van der Wal & Geskus (2011) but the *exact* CD4 threshold the simulation uses to trigger initiation isn't surfaced in the case-study chapter, and I couldn't run R in this environment to check empirical CD4-at-initiation. Please confirm whether 350 is the right number, or substitute the actual threshold the simulation uses (200, 350, 500, or a sliding rule). If you'd rather not commit to a number in the case study, the safe alternative is to drop the parenthetical entirely and keep just "a CD4-driven initiation rule" — let me know.

### 25. "depending" — ITT vs per-protocol
- **Anchor:** "Intention-to-treat or per-protocol, depending on the question"
- **Comment:** Clarify which question aligns with which analysis? Or refer to a later estimand table with the plan for each treatment policy and intercurrent event?
- **Resolution:** Replaced the row's right-hand cell with an explicit per-protocol commitment and a forward reference: "Per-protocol: we estimate the effect of *adhering to* each sustained strategy, since the scientific question is about treatment strategies rather than randomized assignment. An ITT analogue is not natural in an observational target-trial emulation because there is no baseline randomization to be 'intended.' Chapter 1.2 discusses the ITT-vs-per-protocol choice as part of the Causal Roadmap, and Chapter 3.7 walks through how the per-protocol effect is identified from observational data via g-methods."

### 26. "deferring HAART" — does that include never-initiate?
- **Anchor:** "…compared with a strategy of deferring HAART"
- **Comment:** Is never starting HAART also among the compared strategies?
- **Resolution (with a small editorial decision — flag if you disagree).** I rewrote the question to make the comparator explicit: "…compared with a strategy of **never initiating HAART during follow-up**." I also updated the corresponding "Treatment strategies" row in the target trial table to read: "(1) Initiate HAART at the first interval at which CD4 falls below the contemporaneous clinical threshold (e.g., 350 cells/µL) and continue thereafter, vs. (2) **never initiate HAART during follow-up**." This static "always-when-indicated vs never" contrast matches the marginal-structural-model literature on this cohort and aligns with the ever/never point-treatment contrast used in Chapters 2.1–2.3. Chapters 3.7 and 3.9 explore richer dynamic regimes that vary the CD4 threshold.
- **FLAGGED — confirm.** If your intent in the chapter is a *dynamic* contrast (e.g., "initiate at CD4 < 350 vs. initiate at CD4 < 200") rather than a static "always-when-indicated vs never," tell me and I'll switch the framing back to two dynamic regimes.

---

## Other files in the draft set

`00-short-course.docx`, `01-01-regression-vs-causal.docx`, `01-02-causal-roadmap.docx`, `02-01-gcomputation.docx`, `02-02-iptw.docx`, `02-03-doubly-robust-tmle.docx`, `02-04-superlearner.docx`, `02-05-advanced-tmle.docx` — **no reviewer comments found** in the docx XML. Nothing to address in those files in this pass.

---

## Open items requiring your input

Status of the five items flagged in the first pass:

1. **Glossary appendix** (item 1). ✅ Done — [`glossary.qmd`](glossary.qmd) created and registered.
2. **Pre-specified analysis / clean-room reframing** (item 14). ✅ Done — new subsection in [`01-02-causal-roadmap.qmd`](01-02-causal-roadmap.qmd) and tweak to the closing sentence of [`03-06-rwd-meets-rct-hybrid-designs.qmd`](03-06-rwd-meets-rct-hybrid-designs.qmd).
3. **"Who is this book for?" section in `index.qmd`** (item 16). ✅ Done — inserted in [`index.qmd`](index.qmd).
4. **CD4 initiation threshold in the case-study question** (item 24). ✅ Done (recommended Option A applied — see *Third-pass resolution* below).
5. **Dynamic vs static comparator** (item 26). ✅ Done (recommended Option (i) applied — see *Third-pass resolution* below).

---

### Third-pass resolution of items 4 and 5

Both clarifications were added to a single rewritten paragraph in [`00-case-study.qmd`](00-case-study.qmd), placed immediately after the bolded primary causal question and before the "sustained treatment strategies" paragraph. The block now reads (in three paragraphs):

1. **Threshold definition (item 4, paragraph 1).** Keeps "*initiate HAART at the first interval at which CD4 falls below a fixed threshold (e.g., 350 cells/µL)*" as the concrete intervention regime. No change to the numerical anchor.
2. **Intervention regime vs DGP (item 4, paragraph 2 — new).** Explicit one-paragraph distinction: the threshold defines the *intervention regime we want to evaluate*, **not** the policy that generated the observed data. The `haartdat` simulation's actual initiation process is a probabilistic, CD4/age/sex-driven hazard model. The causal analysis estimates what would have happened under the cleaner threshold-based strategy despite the noisier observed initiation process. This forecloses the "is 350 the simulation's rule?" misreading and makes the intervention/DGP distinction available for use later in the book.
3. **Dynamic-vs-static framing bridge (item 5, paragraph 3 — rewritten).** Explicitly notes (a) that the static "never initiate" comparator is the longitudinal analogue of the *ever vs never* contrast used in Chapters [2.1](02-01-gcomputation.qmd)–[2.3](02-03-doubly-robust-tmle.qmd) on the baseline snapshot; (b) that the dynamic-vs-static framing is *deliberate* — it is what makes a single regression insufficient and what motivates the g-methods in Part III; and (c) that the clinically most relevant dynamic-vs-dynamic contrast (e.g., CD4 < 350 vs CD4 < 200) is developed in Chapters [3.7](03-07-longitudinal-td-confounding.qmd) and [3.9](03-09-illustrated-ltmle.qmd). Readers who want that comparison have a clear forward pointer.

The target-trial table on the next page already carries the same threshold language ("e.g., 350 cells/µL") and the same "never initiate HAART during follow-up" comparator, so it remains consistent with the rewritten paragraph and was not edited further. All flagged items from the first two passes are now closed.

---

### Deep-dive on item 4: CD4 initiation threshold

**Current state.** The case-study question in [`00-case-study.qmd:150`](00-case-study.qmd) and the corresponding target-trial table row both read "initiate HAART at the first interval at which CD4 falls below the contemporaneous clinical threshold (e.g., 350 cells/µL) and continue thereafter." I picked 350 as a placeholder and hedged with "e.g.,".

**Why this is not a small decision.**

The threshold defines the *dynamic regime under intervention* — it is the strategy we are estimating the counterfactual mortality risk under. Changing the number changes the causal estimand, the positivity profile, and the numerical answer reported in every subsequent chapter that touches the longitudinal HAART question. Specifically:

- **Estimand identity.** "Risk of death under (initiate at CD4 < 350)" and "Risk of death under (initiate at CD4 < 200)" are *different* causal quantities. They cannot be interchanged; one cannot be inferred from the other without re-running the analysis.
- **Positivity.** A higher threshold (e.g., 500) makes more patients "indicated" at baseline and therefore moves more person-time into the treated arm of the contrast. A lower threshold (e.g., 200) makes fewer patients indicated and stresses positivity in the *deferred* arm. Some thresholds may produce empirical positivity violations in this cohort that others don't.
- **Clinical interpretability.** Each threshold corresponds to a real-world guideline era. Choosing one ties the case study to a particular historical context.

**The four realistic options.**

| Option | What goes in the chapter | Pros | Cons |
|---|---|---|---|
| **A. Keep "e.g., 350 cells/µL"** | Status quo. Hedged single example. | Concrete enough to anchor reader intuition; aligns with mid-2000s WHO/DHHS guidance, which is roughly when the marginal-structural-model literature this cohort comes from was written. | Risk that a reader treats it as the actual threshold the simulation uses (it isn't — see below). |
| **B. Drop the threshold from the question entirely** | "…a CD4-driven initiation rule" without naming a number. The specific threshold is introduced only when Chapter 3.7/3.9 actually estimates effects under it. | Lowest commitment; defers the choice to the chapter that has to make it anyway. | Reader loses a concrete mental anchor; "clinically indicated" stays slightly abstract. |
| **C. Use multiple thresholds throughout** | Frame the case study as a *family* of dynamic-regime contrasts (e.g., 200 vs 350 vs 500). | Scientifically the strongest framing; matches the modern *When to Start* literature and lets the longitudinal chapters show how dynamic-regime contrasts work. | Largest revision. Means rewriting the target trial table, the "primary causal question" block, and the Synthesis chapter to handle multiple estimands. Probably not justified just by Carrie's comment. |
| **D. Match the simulation's actual data-generating process** | Look up (or infer empirically from the data) the rule used by the Van der Wal & Geskus simulation, and write *that* rule into the case study. | The most defensible choice: the regime under intervention matches the regime under the simulation, so the "naive vs g-methods" comparisons later have a known ground truth. | Requires reading the simulation code (or fitting `P(haartind | CD4, sex, age)` empirically) to identify the rule. The simulation actually uses a *probabilistic* model of initiation, not a hard threshold, so the answer is "no single threshold — the model is logistic in CD4." Writing that honestly is a small extra paragraph. |

**Important factual point I should have stated in the first pass.** The `haartdat` simulation in `ipw` does *not* initiate HAART by a hard threshold rule. It uses a probabilistic model in which the hazard of HAART initiation increases as CD4 falls (and depends on sex and age). There is no "the threshold is X" in the data-generating process. So *any* concrete threshold in the case study is describing the **intervention strategy we are evaluating**, not the **policy under which the observed data were generated**. These are conceptually distinct objects and it's worth being explicit about that in the chapter, regardless of which option above you pick.

**My recommendation.** Option **A** (status quo) for the case-study chapter, but with one extra sentence explicitly distinguishing the *intervention regime* (a threshold rule the analyst picks) from the *observed initiation process* (a probabilistic model in this simulation). If you'd rather not name a number at all, switch to **B**. I would not pick **C** unless you want to substantially rework the synthesis chapter. **D** is the most rigorous but only worth doing if the synthesis chapter is going to compare estimated effects against a known truth.

**Action requested.** Tell me which of A/B/C/D you want, and (if A) whether 350 is the right anchor number or whether you'd prefer 200 or 500. I'll apply the change once you decide.

---

### Deep-dive on item 5: dynamic vs static comparator

**Current state.** The case-study question in [`00-case-study.qmd:150`](00-case-study.qmd) compares two strategies:

> initiate HAART at the first interval at which CD4 falls below the threshold and continue thereafter
> vs.
> never initiate HAART during follow-up

The first arm is a *dynamic* regime (its action depends on the time-varying CD4 covariate); the second is a *static* regime (always 0). I picked this hybrid because it's the contrast that's already present in the point-treatment chapters (ever vs never) and it produces a single, intuitive risk difference.

**Why this is not a small decision either.** Static and dynamic comparators don't just differ in difficulty — they answer materially different scientific questions, and they exercise different parts of the handbook's machinery.

**The three realistic framings.**

| Framing | Contrast | Scientific question it answers | Where it fits in the handbook |
|---|---|---|---|
| **Static vs static** | Always treat at time zero vs never treat. | "What if we forced everyone onto HAART from seroconversion vs forced no one onto it?" | Matches the *ever vs never* comparison in Chapters 2.1–2.3 exactly. Cleanest pedagogy for point-treatment methods. **Weakness:** "Treat everyone immediately at seroconversion" is not a clinical strategy anyone proposed; it is a useful technical benchmark but not a clinically interpretable target. |
| **Dynamic vs static** (current) | Treat when CD4 falls below threshold, vs never treat. | "What if everyone followed the CD4-driven guideline vs no one was ever treated?" | Bridges the point-treatment chapters and the longitudinal chapters. The dynamic arm is what makes treatment-confounder feedback necessary to handle (you can't estimate this with a single regression). The static comparator keeps the contrast interpretable and concrete. **Weakness:** Half-and-half is intellectually a little awkward — it exercises dynamic-regime machinery on one arm and not the other. |
| **Dynamic vs dynamic** | Treat at CD4 < 350 vs treat at CD4 < 200 (or similar). | "What if guidelines recommended early initiation vs late initiation?" | The clinically most interesting and policy-relevant contrast. Matches the *When to Start Consortium* literature. **Weakness:** Two dynamic arms is harder to teach and asks the reader to track two threshold rules at once. Less clean for the synthesis-chapter narrative of "naive analysis fails → g-methods recover the right number" because both arms involve the same machinery. |

**Why the choice matters for downstream chapters.**

- **Chapter 3.7 (Time-Dependent Confounding) and Chapter 3.9 (Longitudinal TMLE)** are the chapters that actually do the longitudinal estimation. If the case-study contrast is static-vs-static or dynamic-vs-static, those chapters mostly demonstrate how to handle the dynamic-arm complications; if it's dynamic-vs-dynamic, both arms are dynamic and the contrast is between two MSM-style regimes.
- **The synthesis chapter ([`04-03-case-study-synthesis.qmd`](04-03-case-study-synthesis.qmd))** tells the punchline of the running case study. Whichever framing you pick, the synthesis number changes. A reader's "aha" moment from the synthesis depends on the contrast being clinically interpretable, which slightly favors dynamic-vs-static (current) or dynamic-vs-dynamic over static-vs-static.
- **Point-treatment chapters 2.1–2.3** currently use *ever vs never* as a static-vs-static contrast on the baseline snapshot. None of the three framings above for the *longitudinal* case-study question conflict with that — the point-treatment chapters operate on a different, simpler version of the data — but it's worth being explicit that the point-treatment contrast and the longitudinal contrast are *different estimands*. If you pick dynamic-vs-static for the longitudinal case study, the chapter should explicitly say "the ever-vs-never contrast in Chapters 2.1–2.3 is a simplified, static-vs-static stand-in for this dynamic-vs-static contrast."

**My recommendation.** Keep the current **dynamic vs static** framing as the primary case-study contrast and add a single sentence to the case-study chapter explicitly noting (a) that the point-treatment chapters use a simpler static-vs-static (ever vs never) version of the same question, and (b) that Chapters 3.7 and 3.9 also evaluate a richer dynamic-vs-dynamic contrast (e.g., CD4 < 350 vs CD4 < 200) for readers who want the clinically most relevant question. This preserves the synthesis-chapter narrative arc, keeps the case study concrete, and forecloses the "is this the same as ever-vs-never?" confusion.

**Action requested.** Confirm whether you want me to:
- (i) keep the current dynamic-vs-static framing and add the bridging sentence described above (my recommendation), or
- (ii) switch the primary case-study contrast to static-vs-static (always-at-time-zero vs never) — simpler and matches Chapters 2.1–2.3 exactly, or
- (iii) switch to dynamic-vs-dynamic — most clinically interesting but requires reworking the synthesis chapter.

Once you say which, I'll make the corresponding edits.
