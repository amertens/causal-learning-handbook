#!/bin/bash
# _render-workshop.sh
# Renders the full interactive workshop book with WebR.
#
# Strategy:
#   1. Back up original .qmd files
#   2. Convert {r} chunks tagged #| interactive: true to {webr-r}
#   3. Convert remaining {r} chunks to {webr-r} (so knitr state is consistent)
#   4. Fix source() calls to use download.file() from GitHub
#   5. Remove inline R expressions that reference webr-computed variables
#   6. Render with the workshop profile
#   7. Restore originals

set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

QUARTO="/c/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.exe"

PART2_FILES="02-01-gcomputation.qmd 02-02-iptw.qmd 02-03-doubly-robust-tmle.qmd 02-04-tmle-teaching-examples.qmd 02-05-superlearner.qmd"

echo "=== Step 1: Back up originals ==="
mkdir -p _originals_backup
for f in *.qmd; do
  cp "$f" "_originals_backup/$f"
done

echo "=== Step 2: Convert Part II chunks to {webr-r} ==="
python3 << 'PYEOF'
import os, glob

files = glob.glob("02-0*.qmd")
for fname in files:
    with open(fname, "r", encoding="utf-8") as f:
        lines = f.readlines()
    new_lines = []
    converted = 0
    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()
        if stripped.startswith("```{r") and not stripped.startswith("```{webr-r"):
            skip = False
            for j in range(1, min(5, len(lines) - i)):
                ns = lines[i+j].strip() if i+j < len(lines) else ""
                if "eval: false" in ns or "eval=FALSE" in ns or "include: false" in ns:
                    skip = True
                    break
                if not ns.startswith("#|"):
                    break
            if not skip:
                new_lines.append(line.replace("{r", "{webr-r", 1))
                converted += 1
                # Remove #| interactive: true line if present
                if i+1 < len(lines) and "#| interactive: true" in lines[i+1]:
                    i += 2
                    continue
            else:
                new_lines.append(line)
        else:
            new_lines.append(line)
        i += 1
    with open(fname, "w", encoding="utf-8") as f:
        f.writelines(new_lines)
    print(f"  {fname}: {converted} chunks converted")
PYEOF

echo "=== Step 3: Fix source() for WebR ==="
for f in $PART2_FILES; do
  sed -i 's|source("_data/load-haart-cohort.R")|download.file("https://raw.githubusercontent.com/amertens/causal-learning-handbook/main/_data/load-haart-cohort.R", "load-haart-cohort.R")\nsource("load-haart-cohort.R")|g' "$f" 2>/dev/null
done

echo "=== Step 4: Fix inline R in 02-04 ==="
sed -i "s|\`r round(abs(tmle_ate)\*100, 1)\`|several|g" 02-04-tmle-teaching-examples.qmd 2>/dev/null
sed -i "s|\`r round(tmle_ci\[1\]\*100, 1)\`, \`r round(tmle_ci\[2\]\*100, 1)\`|see results table above|g" 02-04-tmle-teaching-examples.qmd 2>/dev/null

echo "=== Step 5: Render workshop ==="
"$QUARTO" render --profile workshop

echo "=== Step 6: Restore originals ==="
for f in _originals_backup/*.qmd; do
  base=$(basename "$f")
  cp "$f" "$base"
done
rm -rf _originals_backup

echo "=== Done ==="
echo "Static:   amertens.github.io/causal-learning-handbook/"
echo "Workshop: amertens.github.io/causal-learning-handbook/workshop/"
