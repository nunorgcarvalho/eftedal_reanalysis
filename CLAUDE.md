# eftedal_reanalysis — notes for agents

- **Fitting in Python, figures in R.** `src/eftedal/` does all estimation and writes tidy
  CSVs to `results/`; `R/` only reads those CSVs and draws ggplot2 figures. The figures
  should keep the look of Eftedal et al. Fig. 2 (and the NHGRI 2026 poster).
- **Environment:** `source ~/alkes_nuno/.venv_py13/bin/activate`, with the package installed
  editable (`pip install -e .`). R: `module load gcc/14.2.0 R/4.4.2`.
- **Tests:** `pytest` from the repo root. Reproductions of published and legacy numbers
  are regression tests, so keep them passing.
- **Untracked by design:** `archive/` (the 2025 R code and write-ups; read-only reference,
  never commit it), `scratch/` (personal notebooks) and `results/` (regenerable).
- **Data:** `data/eftedal_osf/` holds the authors' public correlation files, which are
  aggregate statistics only. Note the duplicated 2nd-cousin row described in
  `data/README.md`.
- **Conventions:** every analysis choice (weighting, RSS scale, cohort, m–r handling,
  step-sibling formula, ...) is a field of `FitSpec`, never a hard-coded branch in a
  script. Keep the "legacy" options that reproduce earlier results next to the
  corrected defaults.
