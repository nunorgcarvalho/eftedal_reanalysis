# ASHG 2026

Materials for the poster *Relatives' similarity in test score performance is explained by a wide range of
genetic and environmental parameters* (Nuno R. G. Carvalho and Alexander Gusev), ASHG 2026 Annual Meeting,
Montréal.

| folder | what |
|--------|------|
| [`poster/`](poster/) | The poster ([PDF](poster/poster.pdf), 72 × 36 in): `poster.tex`, the pedigree diagrams in `components/`, and the logo script (logos themselves are not tracked). |
| [`supplemental_slides/`](supplemental_slides/) | The supplemental slides ([PDF](supplemental_slides/supplemental_slides.pdf)) linked from the poster's QR code: every core model, cross-validation, and clarifications to the poster. |
| [`QnA/`](QnA/) | Preparation notes ([PDF](QnA/QnA.pdf)): likely questions with short answers, figures and tables. AI-written; see the disclosure at the top of the PDF. |

## Archived as presented

These materials are frozen as presented at ASHG 2026:
- **Code and built PDFs:** tag [`ASHG2026`](../../../tree/ASHG2026).
- **The `results/` they were built from**, plus the Python/R environment: tag `ASHG2026-results`
  (branch `archive/ashg2026-results`).

The analysis has moved on since then (Phase 2: sex-composition and family-side multipliers;
t-based jackknife intervals). So `build.sh` **refuses to rebuild** when the code or the results
differ from the tag; pass `--force` to override. To rebuild faithfully, check out the tag in its
own worktree, copy in that snapshot's `results/`, and run `build.sh` there.

## Building

Everything is generated from the analysis results (`results/`, not tracked). From the repository root:

```bash
python scripts/01_replicate_eftedal.py   # ... through scripts/08_param_corr.py (see the top-level README)
python scripts/poster_additive.py
python scripts/poster_eftedal_predictions.py
Rscript R/make_figures.R
Rscript R/poster_figures.R               # poster and slide figures -> results/figures/poster/
python ASHG/QnA/make_tables.py           # Q&A tables and numbers -> results/qna/
bash ASHG/build.sh                       # -> the three PDFs
```

The poster and slide figures are drawn by `R/poster_figures.R`; the Q&A's numbers all come from
`ASHG/QnA/make_tables.py`, so none is typed by hand.
