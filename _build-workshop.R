#!/usr/bin/env Rscript
# _build-workshop.R
#
# Generates the workshop version of the book: a self-contained Quarto project
# in _workshop-src/ whose tagged chunks run in the browser via WebR.
#
# Usage:
#   Rscript _build-workshop.R
#   quarto render _workshop-src
#
# Chunks tagged with  #| interactive: true  anywhere in their option block
# become {webr-r} interactive cells. All other chunks remain {r} and are
# executed by knitr as usual. The original source files are NOT modified.
#
# The generated _workshop-src/_quarto.yml merges the book structure from
# _quarto.yml with the workshop settings from _quarto-workshop.yml and writes
# its output to ../docs/workshop, so the interactive site is built by one
# command and cannot drift from the book's chapter list.

library(fs)

workshop_dir <- "_workshop-src"

if (dir_exists(workshop_dir)) dir_delete(workshop_dir)
dir_create(workshop_dir)

# ---- 1. Copy chapter sources, converting the tagged chunks -----------------

# A {webr-r} cell is not executed by knitr: it only runs when the reader presses
# Run. So converting a chunk in place would strip its objects out of the knitr
# session and break every later chunk that uses them. Each tagged chunk is
# therefore emitted twice: an editable {webr-r} cell for the reader, and a
# hidden {r} copy that still runs so the rest of the chapter computes.
convert_interactive <- function(lines) {
  out <- character()
  converted <- 0L
  n <- length(lines)
  i <- 1L
  while (i <= n) {
    if (!grepl("^```\\{r", lines[i])) {
      out <- c(out, lines[i]); i <- i + 1L; next
    }
    open_line <- lines[i]
    j <- i + 1L
    opts <- character()
    while (j <= n && grepl("^#\\|", lines[j])) { opts <- c(opts, lines[j]); j <- j + 1L }
    body_start <- j
    while (j <= n && !grepl("^```\\s*$", lines[j])) j <- j + 1L
    body <- if (body_start <= j - 1L) lines[body_start:(j - 1L)] else character()

    if (any(grepl("^#\\|\\s*interactive:\\s*true\\s*$", opts))) {
      keep <- opts[!grepl("^#\\|\\s*(interactive|label|code-fold|include|echo):", opts)]
      out <- c(out,
               "```{webr-r}", keep, body, "```",
               "",
               "```{r}", "#| include: false", body, "```")
      converted <- converted + 1L
    } else {
      out <- c(out, open_line, opts, body, "```")
    }
    i <- j + 1L
  }
  list(lines = out, converted = converted)
}

src_files <- c(dir_ls(".", glob = "*.qmd", recurse = FALSE),
               dir_ls(".", glob = "*.Rmd", recurse = FALSE))

total <- 0L
for (f in src_files) {
  res <- convert_interactive(readLines(f, warn = FALSE))
  writeLines(res$lines, path(workshop_dir, path_file(f)))
  if (res$converted > 0L) {
    message("[workshop] ", path_file(f), ": ", res$converted, " interactive cell(s)")
    total <- total + res$converted
  }
}

# ---- 2. Copy the supporting files the chapters need ------------------------

for (f in c("book.bib", "style.css", "packages.bib")) {
  if (file_exists(f)) file_copy(f, path(workshop_dir, f), overwrite = TRUE)
}
for (d in c("_data", "_extensions")) {
  if (dir_exists(d)) dir_copy(d, path(workshop_dir, d), overwrite = TRUE)
}

# Seed the workshop project with the book's freeze cache. Chapters with no
# interactive cells are copied verbatim, so their cached results still match
# and knitr skips them; the rest re-execute only the chunks that stayed {r}.
if (dir_exists("_freeze")) dir_copy("_freeze", path(workshop_dir, "_freeze"), overwrite = TRUE)

# ---- 3. Write a self-contained _quarto.yml ---------------------------------
# Take the book structure from _quarto.yml, then apply the workshop overrides
# and repoint the output directory.

base <- readLines("_quarto.yml", warn = FALSE)
over <- readLines("_quarto-workshop.yml", warn = FALSE)

# Drop the profile block from the base config: the generated project is
# rendered directly, not through a profile.
prof_start <- grep("^profile:", base)
if (length(prof_start)) {
  nxt <- grep("^[a-z]", base)
  nxt <- nxt[nxt > prof_start[1]]
  base <- base[-(prof_start[1]:(if (length(nxt)) nxt[1] - 1L else length(base)))]
}

base <- sub("^  output-dir: docs.*$", "  output-dir: ../docs/workshop", base)
base <- sub('^  title: "Causal Learning Handbook"$',
            '  title: "Causal Learning Handbook (Interactive Workshop)"', base)
base <- sub("^(\\s*)site-url: .*$",
            "\\1site-url: https://amertens.github.io/causal-learning-handbook/workshop/",
            base)

# Inside the workshop, the sidebar button should lead back to the static book
# rather than to the page the reader is already on.
base <- sub("^(\\s*-\\s*)icon: play-circle-fill$", "\\1icon: book", base)
base <- sub('^(\\s*)text: "Interactive \\(WebR\\)"$', '\\1text: "Static version"', base)
base <- sub("^(\\s*)href: https://amertens.github.io/causal-learning-handbook/workshop/$",
            "\\1href: https://amertens.github.io/causal-learning-handbook/", base)

# Append only the workshop-only top-level sections. `format:` already exists in
# the base config, and YAML rejects a duplicate key, so its overrides are
# applied to the base block instead of appended.
keep <- FALSE
extra <- character()
for (ln in over) {
  if (grepl("^(filters|webr):", ln)) keep <- TRUE
  else if (grepl("^[a-z]", ln) && !grepl("^\\s", ln)) keep <- FALSE
  if (keep) extra <- c(extra, ln)
}

# Code is shown unfolded in the workshop so the editable cells are visible.
base <- sub("^(\\s*)code-fold: .*$", "\\1code-fold: false", base)

writeLines(c(base, "", extra), path(workshop_dir, "_quarto.yml"))

message("\n", total, " interactive cell(s) across the book")
message("Workshop project written to: ", workshop_dir, "/")
message("Now run: quarto render ", workshop_dir)
