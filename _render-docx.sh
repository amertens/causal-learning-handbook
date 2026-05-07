#!/bin/bash
# _render-docx.sh
# Render every handbook chapter to a separate .docx file in docs-docx/.
# Intended for collaborator review with Track Changes / Comments in Word.
#
# Strategy:
#   The repo is configured as a Quarto book (type: book in _quarto.yml),
#   which means `quarto render <file>.qmd --to docx` produces one big
#   book.docx instead of a per-chapter docx. To get per-chapter output,
#   we temporarily stash _quarto.yml so each chapter renders as a
#   standalone document, then restore it on exit.
#
# Usage:
#   ./_render-docx.sh              # render all chapters
#   ./_render-docx.sh chapter.qmd  # render a single chapter (e.g. for testing)

set -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

QUARTO="/c/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.exe"
OUT_DIR="docs-docx"
mkdir -p "$OUT_DIR"

# Restore _quarto.yml on any exit, including Ctrl-C and errors.
restore_quarto_yml() {
  if [[ -f _quarto.yml.docx-bak ]]; then
    mv -f _quarto.yml.docx-bak _quarto.yml
    echo "    (restored _quarto.yml)"
  fi
}
trap restore_quarto_yml EXIT INT TERM

# Stash the book-project config so chapters render as standalone documents.
if [[ -f _quarto.yml ]]; then
  cp _quarto.yml _quarto.yml.docx-bak
  rm -f _quarto.yml
fi

# Chapter list, in the order they appear in _quarto.yml
ALL_CHAPTERS=(
  "index.qmd"
  "00-short-course.qmd"
  "00-case-study.qmd"
  "01-01-regression-vs-causal.qmd"
  "01-02-causal-roadmap.qmd"
  "02-01-gcomputation.qmd"
  "02-02-iptw.qmd"
  "02-03-doubly-robust-tmle.qmd"
  "02-04-superlearner.qmd"
  "03-07-longitudinal-td-confounding.qmd"
  "03-09-illustrated-ltmle.qmd"
  "03-09b-ltmle-by-hand.qmd"
  "03-03b-adherence-shift-methods.qmd"
  "03-02-positivity-overlap.qmd"
  "03-05-advanced-diagnostics.qmd"
  "03-12-sensitivity-analysis.qmd"
  "05-01-applied-lmtp-workflow.qmd"
  "03-04-effect-modification.qmd"
  "03-04b-interactions-rwe.qmd"
  "03-10-tmle-mediation.qmd"
  "03-11-time-to-event-competing-risks.qmd"
  "02-05-advanced-tmle.qmd"
  "03-06-rwd-meets-rct-hybrid-designs.qmd"
  "03-06b-data-fusion-tutorial.qmd"
  "03-13-transportability-generalizability.qmd"
  "03-08-tmle-clean-room.qmd"
  "04-03-case-study-synthesis.qmd"
  "anotated_bibliography.Rmd"
  "other_resources.Rmd"
)

# Allow rendering a single chapter when passed as an argument
if [[ $# -gt 0 ]]; then
  CHAPTERS=("$@")
else
  CHAPTERS=("${ALL_CHAPTERS[@]}")
fi

FAILED=()
SUCCEEDED=()

for f in "${CHAPTERS[@]}"; do
  if [[ ! -f "$f" ]]; then
    echo "SKIP: $f (file not found)"
    FAILED+=("$f (file not found)")
    continue
  fi

  echo ""
  echo "==> Rendering $f"

  base="${f%.*}"

  # Snapshot which sibling _files/ and _cache/ dirs already existed,
  # so we only clean up the ones Quarto creates for THIS render.
  pre_files_existed=0
  pre_cache_existed=0
  [[ -d "${base}_files" ]] && pre_files_existed=1
  [[ -d "${base}_cache" ]] && pre_cache_existed=1

  # Render to docx as a standalone document (no project context).
  # --bibliography is passed explicitly because we hid _quarto.yml.
  if "$QUARTO" render "$f" --to docx \
       --bibliography book.bib \
       --metadata toc:true \
       --metadata number-sections:true \
       --metadata toc-depth:3; then
    if [[ -f "${base}.docx" ]]; then
      mv -f "${base}.docx" "$OUT_DIR/"
      SUCCEEDED+=("$f")
      echo "    OK -> $OUT_DIR/${base}.docx"
    else
      echo "    WARN: render reported success but ${base}.docx was not found"
      FAILED+=("$f (output file missing)")
    fi
  else
    echo "    FAIL: render error"
    FAILED+=("$f (render error)")
  fi

  # Clean up Quarto's per-chapter intermediates if they did not already exist.
  if [[ $pre_files_existed -eq 0 && -d "${base}_files" ]]; then
    rm -rf "${base}_files"
  fi
  if [[ $pre_cache_existed -eq 0 && -d "${base}_cache" ]]; then
    rm -rf "${base}_cache"
  fi
  # If they did exist before, remove only the docx-specific subdirs we added.
  rm -rf "${base}_files/figure-docx" "${base}_cache/docx" 2>/dev/null || true
done

echo ""
echo "=== Summary ==="
echo "Succeeded: ${#SUCCEEDED[@]}"
echo "Failed:    ${#FAILED[@]}"
if [[ ${#FAILED[@]} -gt 0 ]]; then
  printf '  - %s\n' "${FAILED[@]}"
fi
echo "Output:    $OUT_DIR/"
