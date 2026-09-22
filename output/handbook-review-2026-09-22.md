# Handbook review: content and organization

**Date:** 2026-09-22
**Scope:** all 33 files listed in `_quarto.yml` (28 substantive chapters plus index and appendices)
**Method:** measured, not eyeballed. Counts below come from parsing the sources: prose words exclude code blocks, maths counts separate display equations from inline expressions, and scaffolding is detected per chapter.

## The book in numbers

| Measure | Value |
|---|---|
| Prose words | 89,516 |
| Lines of code | 6,711 |
| Code chunks | 331, of which **52 (16%) are `eval: false`** |
| Display equations | 81 |
| Inline maths expressions | 968 |
| Figures | 52, plus 5 Mermaid diagrams |
| **Visual density** | **0.64 figures per 1,000 prose words** |
| Interactive (WebR) cells | 53, confined to **4 chapters** |
| Chapters with learning objectives | 11 of 28 |
| Chapters with check-your-understanding | 15 of 28 |
| Chapters with a "what comes next" | 9 of 28 |
| Two largest chapters | 03-05 (12,690 words) and 03-11 (8,699) = **24% of all prose** |

---

# Reader 1: the new student who is afraid of the maths

This reader wants pictures, things to click, and a path that tells them where they are. Here is what they meet.

## What already works

The running HAART cohort threads through every chapter, so the reader never has to learn a new dataset to learn a new method. The 67 collapsible check-your-understanding questions are genuinely good: they ask for judgment rather than recall, and the answers are paragraphs rather than a letter. `code-fold: show` means code is present but escapable. Callout boxes carry the warnings. The short course is a real on-ramp.

## The five things that will lose them

### 1. The interactive promise is broken

The sidebar carries a prominent **"Interactive (WebR)"** button on every page. It points at a site last rebuilt on **12 June 2026**, more than three months ago. That build:

- is missing the four newest chapters (estimand, design and time zero, missing data, optimal regimes),
- still contains the mediation chapter that was removed from the book,
- predates the entire reviewer-notes pass, the prose humanization, and every accuracy fix since.

Worse, even when rebuilt, only **4 of 28 chapters** contain any interactive cell at all: g-computation (8), IPTW (10), doubly robust TMLE (11), and Super Learner (24). A student who clicks "Interactive" and opens the roadmap, the case study, positivity, or any Part III chapter finds a static page with an interactive-looking button in the corner.

This is the single highest-value fix for this reader, and it is two separate fixes: rebuild the site, and spread `interactive: true` beyond Part II.

### 2. Visual density is a third of what a visual learner needs

0.64 figures per 1,000 words, against roughly 2 to 3 for a textbook built around pictures. The distribution is worse than the average:

| Chapter | Prose words | Figures |
|---|---|---|
| `01-02-causal-roadmap` | 5,231 | 1 figure + 1 diagram |
| `02-03-doubly-robust-tmle` | 3,306 | 1 |
| `02-05-advanced-tmle` | 3,508 | 1 |
| `03-04b-interactions-rwe` | 3,607 | 0 figures, 1 diagram |
| `index` | 1,072 | 0 |
| `04-03-case-study-synthesis` | 931 | 0 |

The roadmap chapter is the conceptual spine of the book and carries one figure across 5,231 words. Doubly robust TMLE, the moment the whole method clicks or does not, has one figure and 75 inline maths expressions.

### 3. No learning objectives anywhere in Part II

Every chapter in **Part II: Core Estimators** (g-computation, IPTW, doubly robust TMLE, Super Learner) opens without a statement of what the reader will be able to do. These are exactly the four chapters where a nervous reader most needs to be told the destination before the derivation starts. Part I, Part IV, and most of Part VI have objectives; the core does not.

### 4. The reading path breaks 19 times

Only 9 of 28 chapters end with a "what comes next" (a further chapter had one under a different capitalization, and 7 of the files without one are the index and appendices, which do not need it). A reader who finishes L-TMLE concepts, or the adherence chapter, or the applied LMTP workflow, reaches the bottom of the page with no idea where to go. The sidebar is the only navigation, and the sidebar is a flat list of 32 items.

### 5. The densest maths in the book sits in the main flow

`03-09b-ltmle-by-hand` has the highest maths density of any chapter (10.4 on the combined display-plus-inline measure, against a median near 2.5), no learning objectives, and no check-your-understanding. The index calls it "optional"; the sidebar does not, so a reader working through Part III in order walks straight into it.

## What to do, in order

1. **Rebuild the workshop site** (`Rscript _build-workshop.R && quarto render --profile workshop`). One command, removes a three-month-old inconsistency. Then rebuild it as part of the release routine, not occasionally.
2. **Add `interactive: true` to one or two chunks in every chapter that has none.** Target the chunk where changing a number teaches something: the positivity threshold in 03-02, the grace period in 01-04, the misclassification sensitivity in 03-15, the CD4 initiation threshold in 05-01.
3. **Add learning objectives to the four Part II chapters.** Copy the format already used in Part I. Half a day.
4. **Add "what comes next" to the 19 chapters missing it.** Two or three sentences each, pointing with `@sec-` refs.
5. **Put a visual map on the landing page.** One diagram showing the six parts and the dependency between them would do more for orientation than the current three paragraphs of prose.
6. **Raise the figure count where the concepts are hardest**, not uniformly. The three places a picture would earn the most: the targeting step in 02-03 (the fluctuation moving an initial fit toward the truth, which is a natural animation), the roadmap's seven steps as a single diagram in 01-02, and the weight distribution in 02-02.
7. **Mark `03-09b` as optional in the sidebar**, either by moving it to an appendix or by retitling it so the sidebar says so.

---

# Reader 2: the causal inference expert

This reader is checking whether the book is trustworthy, whether it says anything they do not already know, and whether they would cite it or assign it.

## What earns their respect

The estimand-first framing is now genuinely strong: specifying the estimand before the estimator, with a worked demonstration that the choice of intercurrent-event strategy moves the answer further than any estimator choice would. The clean-room chapter is the most distinctive thing in the book and has no real equivalent in the published literature. Ground-truth simulations appear throughout, so claims are checkable. Package versions are pinned where the results depend on them. And the book reports negative results honestly: the learned treatment rule that loses to treat-all out of sample, the CATE rule indistinguishable from the blanket policy, the g-computation estimate that lands near but not on the truth.

## The five things that will lose them

### 1. The theory is asserted, never derived, even once

| Concept | Chapters mentioning it |
|---|---|
| Efficient influence function | 11 |
| Cross-fitting | 12 |
| Asymptotic linearity / root-n | 4 |
| Donsker conditions | 2 (plus glossary and bibliography) |
| Second-order remainder | 2 |
| Semiparametric efficiency bound | 2 |
| Non-regularity | 1 |

The EIF is invoked in eleven chapters and the von Mises expansion that explains why it is the right object appears nowhere. An expert reading this concludes the book knows the vocabulary rather than the theory. **One appendix** deriving the expansion for the ATE, showing where the second-order remainder comes from, why Donsker or cross-fitting is needed to control the empirical-process term, and why that yields root-n inference, would change how the whole book reads. It costs perhaps 2,500 words and it is the highest-value addition for this reader.

### 2. A sixth of the code does not run, concentrated in two chapters

| Chapter | Chunks | `eval: false` | Dead |
|---|---|---|---|
| `03-05-advanced-diagnostics` | 51 | 33 | **65%** |
| `03-06-rwd-meets-rct-hybrid-designs` | 17 | 9 | **53%** |
| Whole book | 331 | 52 | 16% |

The diagnostics chapter is the longest in the book and two thirds of its code is display-only. An expert reads a non-executing chunk as a chunk that was never run. Some of this is unavoidable (missing packages), but the earlier review established that the `ES.cvtmle()` example in 03-06 runs in about two seconds and could simply be switched on.

### 3. Every demonstration is a single draw

With one partial exception, the book never reports bias and coverage over repeated simulation. A methods handbook is judged on exactly that: does the estimator cover at 95%, and does the naive one fail to? The machinery is already in place in a dozen chapters; what is missing is the loop. **One figure showing coverage across estimators over 500 replications** would carry more weight with this reader than another worked example.

### 4. Chapter 03-14 promises three literatures and delivers 1,560 words

Instrumental variables, negative controls, and proximal inference in 1,560 prose words with **zero display equations**. Proximal inference cannot be stated at that length, let alone motivated. Either expand it substantially or split it and cut proximal inference back to a pointer.

### 5. Chapter 03-05 is a 12,690-word grab-bag

It absorbed the archived sensitivity-analysis chapter and now holds assumption diagnostics, E-values, quantitative bias analysis, negative controls, nuisance-model diagnostics, positivity diagnostics, calibration, power, and the bootstrap. At 14% of the book's prose it has no narrative and an expert cannot navigate it. Splitting it into "diagnostics for a fitted analysis" and "sensitivity analysis for unverifiable assumptions" would restore both.

## Remaining topic gaps

From the earlier scoping review, still open after the four new chapters: **SUTVA and interference** appear once in the entire book, which is untenable in a text covering vaccines or anything with spillover; the **highly adaptive lasso** remains archived, which is a conspicuous hole in a book whose title is targeted learning; and panel methods (difference-in-differences, synthetic control) are absent, which is defensible but should be stated as a scope decision rather than left silent.

---

# Organization

**What is sound.** Six parts with a clear arc from foundations to estimators to longitudinal to diagnostics to application to special topics. Part I now carries the question-framing chapters, which is where they belong. Part IV's sequence of diagnostics, then missing data, then the clean-room workflow reads well.

**Three problems.**

1. **Part VI carries a blanket "In development" warning** over six chapters, four of which (effect modification, interactions, optimal regimes, time-to-event) are as finished as anything in the book. The warning is doing the opposite of its job: it tells the reader to discount mature material. Move the warning onto the two chapters that need it.

2. **The capstone is 931 words.** `04-03-case-study-synthesis` has to carry the payoff of the entire book, and it has no figures, no code, and six short sections. A reader arriving there after 89,000 words deserves the whole argument assembled: one table of every estimate the book produced on the HAART cohort, what each one assumed, and where they disagree.

3. **The sidebar is a flat list of 32 entries** with no indication of difficulty, dependency, or which chapters are optional. Both readers suffer from this, for opposite reasons.

---

# Priority list

| # | Action | Serves | Status |
|---|---|---|---|
| 1 | Rebuild the WebR workshop site and add it to the release routine | Student | **Done.** Build script rewritten as a self-contained project; it also converted 3 cells that the old script silently skipped |
| 2 | Point at the semiparametric theory rather than deriving it | Expert | **Done.** Bickel et al., Tsiatis, van der Vaart, Chernozhukov et al., Hines et al., Kennedy added at the places the EIF and the regularity conditions are invoked |
| 3 | Learning objectives for the four Part II chapters | Student | **Done**, plus six other chapters that lacked them (26 of 28 chapters now have them) |
| 4 | "What comes next" in the chapters missing it | Student | **Done.** 12 substantive chapters gained one; the other 7 without one are the index and appendices |
| 5 | Coverage-over-replications figure for the headline estimators | Expert | **Done.** 500 replications in 02-03: g-computation 0% coverage, IPTW 98%, AIPW 97%, TMLE 96% |
| 6 | Spread `interactive: true` beyond Part II | Student | **Done.** 4 chapters to 9; the new cells are self-contained and use only packages WebR can load |
| 7 | Split 03-05 into diagnostics and sensitivity analysis | Both | Open. Still 12,690 words |
| 8 | Move the "In development" warning off the four mature Part VI chapters | Both | Open |
| 9 | Expand the synthesis chapter into a real capstone with an estimate table | Both | Open. It now has learning objectives and a forward pointer, but is still 931 words |
| 10 | Turn on the executable examples that can run | Expert | **Partly done.** 03-06 from 53% to 35% dead, 03-05 from 65% to 59%, book from 16% to 14%. The rest need `moments` or are template snippets referencing undefined objects |
| 11 | Visual map of the book on the landing page | Student | **Done.** `fig-book-map` |
| 12 | A section on SUTVA and interference | Expert | Open |
| 13 | Expand or split 03-14 (1,560 words for IV, negative controls, and proximal) | Expert | Open |

## What changed on 2026-09-22, after this review

The student-facing items and the citation-only version of item 2 were done the same day. Three things surfaced while doing them:

- **The workshop build was broken in two ways, not one.** Besides being three months stale, `_build-workshop.R` only looked at the first option line of each chunk, so any chunk whose `interactive: true` sat below a `label:` or `fig-cap:` was silently left static. Three cells in Part II were affected. The script is now a self-contained project generator: it writes `_workshop-src/_quarto.yml` from the book's own chapter list, so the interactive site cannot drift from the book again.
- **The coverage simulation changed what it demonstrates.** The first version showed every estimator as nearly unbiased, because a quadratic term omitted from the outcome model does not bias the marginal effect much when treatment depends only linearly on the covariate. Making treatment depend on the squared term too produced the intended result: g-computation overstates the effect by 84% and its intervals never cover.
- **Most of 03-05's non-executing code cannot be revived.** Of its 33 dead chunks, 6 run clean in sequence and only 3 of those produce output; the rest need the uninstalled `moments` package or are template snippets that reference objects the chapter never creates. Splitting the chapter (item 7) is the real fix.
