#!/usr/bin/env Rscript
# _build-workshop.R
#
# Generates the workshop version of the book by copying all .qmd files
# and converting tagged chunks from {r} to {webr-r}.
#
# Usage:
#   Rscript _build-workshop.R
#   quarto render --profile workshop
#
# Chunks tagged with  #| interactive: true  become {webr-r} interactive cells.
# All other chunks remain {r} and are executed by knitr as usual.
# The original files are NOT modified.

library(fs)

workshop_dir <- "_workshop-src"

# Clean and recreate
if (dir_exists(workshop_dir)) dir_delete(workshop_dir)
dir_create(workshop_dir)

# Copy all .qmd files
qmd_files <- dir_ls(".", glob = "*.qmd", recurse = FALSE)

for (f in qmd_files) {
  lines <- readLines(f, warn = FALSE)
  modified <- FALSE
  i <- 1L

  while (i <= length(lines)) {
    if (grepl("^```\\{r", lines[i]) && i < length(lines)) {
      if (grepl("^#\\|\\s*interactive:\\s*true", lines[i + 1])) {
        lines[i] <- sub("\\{r", "{webr-r", lines[i])
        lines[i + 1] <- ""  # blank out the tag (keep line count stable)
        modified <- TRUE
      }
    }
    i <- i + 1L
  }

  out_path <- path(workshop_dir, path_file(f))
  writeLines(lines, out_path)
  if (modified) message("[workshop] Converted: ", f)
}

message("\nWorkshop sources written to: ", workshop_dir, "/")
message("Now run: quarto render --profile workshop")
