#!/usr/bin/env bash
# Build the ASHG 2026 materials: poster/poster.pdf, supplemental_slides/supplemental_slides.pdf and
# QnA/QnA.pdf. Inputs (results/, untracked) come first, from the repo root:
#   python scripts/01_replicate_eftedal.py ... 08_param_corr.py   (see the top-level README)
#   python scripts/poster_additive.py && python scripts/poster_eftedal_predictions.py
#   Rscript R/make_figures.R && Rscript R/poster_figures.R
#   python ASHG/QnA/make_tables.py
set -euo pipefail
cd "$(dirname "$0")"

# ---- guard: these materials are archived as presented (tag ASHG2026; results snapshot at tag
# ASHG2026-results, branch archive/ashg2026-results). Rebuilding from newer code or newer results
# would silently change the poster's numbers. Rebuild faithfully from the tag instead:
#   git worktree add ../ashg2026 ASHG2026
#   git show ASHG2026-results  # then copy that branch's results/ into ../ashg2026/results/
#   bash ../ashg2026/ASHG/build.sh
# Pass --force to build anyway (e.g. a deliberate post-ASHG version).
TAG=ASHG2026
if [[ "${1:-}" != "--force" ]]; then
  REPO="$(git rev-parse --show-toplevel)"
  problems=()
  if ! git -C "$REPO" diff --quiet "$TAG" -- src scripts R 2>/dev/null; then
    problems+=("the code (src/, scripts/, R/) differs from tag $TAG")
  fi
  shopt -s nullglob
  stamps=("$REPO"/results/PROVENANCE/*.txt)
  if (( ${#stamps[@]} == 0 )); then
    problems+=("results/ has no provenance stamps, so which code produced it is unknown")
  fi
  for s in "${stamps[@]}"; do
    c=$(sed -n 's/^commit=//p' "$s"); d=$(sed -n 's/^dirty=//p' "$s")
    if [[ "$d" == 1 ]] || ! git -C "$REPO" diff --quiet "$c" "$TAG" -- src scripts R 2>/dev/null; then
      problems+=("$(basename "$s" .txt) results were made by code other than tag $TAG (commit ${c:0:7}${d/1/, uncommitted changes})")
    fi
  done
  if (( ${#problems[@]} )); then
    echo "WARNING: refusing to rebuild the ASHG 2026 materials:" >&2
    printf '  - %s\n' "${problems[@]}" >&2
    echo "Rebuild from tag $TAG with the ASHG2026-results snapshot (see the top of ASHG/build.sh), or pass --force." >&2
    exit 1
  fi
fi
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
