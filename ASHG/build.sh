#!/usr/bin/env bash
# Build the ASHG 2026 materials: poster/poster.pdf, supplemental_slides/supplemental_slides.pdf and
# QnA/QnA.pdf. Inputs (results/, untracked) come first, from the repo root:
#   python scripts/01_replicate_eftedal.py ... 08_param_corr.py   (see the top-level README)
#   python scripts/poster_additive.py && python scripts/poster_eftedal_predictions.py
#   Rscript R/make_figures.R && Rscript R/poster_figures.R
#   python ASHG/QnA/make_tables.py
set -euo pipefail
cd "$(dirname "$0")"
latex() {   # latex <dir> <file.tex> [passes]
  (cd "$1" && for _ in $(seq "${3:-2}"); do
     pdflatex -interaction=nonstopmode -halt-on-error "$2" >/dev/null || { tail -40 "${2%.tex}.log"; exit 1; }
   done)
}
for f in poster/components/*.tex; do latex poster/components "$(basename "$f")" 1; done
latex poster poster.tex
latex supplemental_slides supplemental_slides.tex
latex QnA QnA.tex 3   # table of contents and cross-references
echo "built $(pwd)/poster/poster.pdf, $(pwd)/supplemental_slides/supplemental_slides.pdf and $(pwd)/QnA/QnA.pdf"
