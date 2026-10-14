#!/usr/bin/env bash
# Build poster.pdf from poster.tex. Figures: run `Rscript R/make_figures.R` first.
set -euo pipefail
cd "$(dirname "$0")"
# standalone components (diagrams) -> components/*.pdf
for f in components/*.tex; do
  (cd components && pdflatex -interaction=nonstopmode -halt-on-error "$(basename "$f")" >/dev/null) \
    || { tail -40 "${f%.tex}.log"; exit 1; }
done
pdflatex -interaction=nonstopmode -halt-on-error poster.tex >/dev/null || { tail -40 poster.log; exit 1; }
pdflatex -interaction=nonstopmode -halt-on-error poster.tex >/dev/null
echo "built $(pwd)/poster.pdf"
