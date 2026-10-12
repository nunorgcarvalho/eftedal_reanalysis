# eftedal_reanalysis

A reanalysis of the relative-pair correlations published by Eftedal et al. (2025,
*PNAS*) for Norwegian national test scores. Each of the 82 categories of relatives has
one observed phenotypic correlation. We fit models that predict those correlations
from genetic relatedness under assortative mating, with optional nonadditive (GxG) and
shared-environment components, and ask how sensitive the conclusions are to modelling
and fitting choices.

## Layout

| path | what |
|------|------|
| `data/` | Input data: the authors' public correlation files and our relationship-type annotation (see `data/README.md`). |
| `src/eftedal/` | Python package: data loading, covariance components, model specification, fitting, uncertainty. |
| `scripts/` | Python entry points that run the analyses and write tables to `results/`. |
| `R/` | ggplot2 figures drawn from the tables in `results/`. |
| `tests/` | pytest suite. |
| `results/` | Generated output (not tracked). |

## Setup

```bash
source ~/alkes_nuno/.venv_py13/bin/activate      # any Python >= 3.11 with numpy/scipy/pandas
pip install -e .
module load gcc/14.2.0 R/4.4.2                    # R with tidyverse, data.table, ggnewscale, ggrepel
```

## Running

```bash
pytest
python scripts/01_replicate_eftedal.py
python scripts/02_replicate_legacy.py
python scripts/03_sensitivity.py
python scripts/04_uncertainty.py
Rscript R/make_figures.R
```
