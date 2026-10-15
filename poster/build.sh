#!/usr/bin/env bash
# Build poster.pdf and supplemental_slides.pdf. Figures: run `Rscript R/make_figures.R`,
# `python scripts/poster_additive.py`, `python scripts/poster_eftedal_predictions.py` and
# `Rscript R/poster_figures.R` first.
set -euo pipefail
cd "$(dirname "$0")"
# standalone components (diagrams) -> components/*.pdf
for f in components/*.tex; do
  (cd components && pdflatex -interaction=nonstopmode -halt-on-error "$(basename "$f")" >/dev/null) \
    || { tail -40 "${f%.tex}.log"; exit 1; }
done
pdflatex -interaction=nonstopmode -halt-on-error poster.tex >/dev/null || { tail -40 poster.log; exit 1; }
pdflatex -interaction=nonstopmode -halt-on-error poster.tex >/dev/null
for i in 1 2; do
  pdflatex -interaction=nonstopmode -halt-on-error supplemental_slides.tex >/dev/null \
    || { tail -40 supplemental_slides.log; exit 1; }
done
echo "built $(pwd)/poster.pdf and $(pwd)/supplemental_slides.pdf"
