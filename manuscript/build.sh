#!/bin/sh
# Build the review manuscript and supplement; no model fitting is performed.
set -eu
cd "$(dirname "$0")"
PYTHON="${AFT_PYTHON:-python3}"
command -v pdflatex >/dev/null 2>&1 || { echo "pdflatex is required (TeX Live or MacTeX)." >&2; exit 1; }
command -v "$PYTHON" >/dev/null 2>&1 || { echo "Python 3 is required for cross-references." >&2; exit 1; }
mkdir -p build_logs
for pass in 1 2 3; do
  for name in main supplement; do
    pdflatex -interaction=nonstopmode -halt-on-error "$name.tex" > "build_logs/$name-pass$pass.txt"
    "$PYTHON" export_crossrefs.py
  done
done
for pass in 1 2 3; do
  pdflatex -interaction=nonstopmode -halt-on-error complete.tex > "build_logs/complete-pass$pass.txt"
done
cp supplement.pdf ESM_1.pdf
printf '\nBuilt main.pdf, ESM_1.pdf, and complete.pdf.\n'
