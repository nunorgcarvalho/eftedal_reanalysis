# ASHG 2026

Materials for the poster *Relatives' similarity in test score performance is explained by a wide range of
genetic and environmental parameters* (Nuno R. G. Carvalho and Alexander Gusev), ASHG 2026 Annual Meeting,
Montréal.

| folder | what |
|--------|------|
| [`poster/`](poster/) | The poster ([PDF](poster/poster.pdf), 72 × 36 in): `poster.tex`, the pedigree diagrams in `components/`, and the logo script (logos themselves are not tracked). |
| [`supplemental_slides/`](supplemental_slides/) | The supplemental slides ([PDF](supplemental_slides/supplemental_slides.pdf)) linked from the poster's QR code: every core model, cross-validation, and clarifications to the poster. |
| [`QnA/`](QnA/) | Preparation notes ([PDF](QnA/QnA.pdf)): likely questions with short answers, figures and tables. AI-written; see the disclosure at the top of the PDF. |

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
