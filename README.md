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

## Model

The predicted correlation for category $k$ is a sum of components

$$\rho_k = V_A K_k(m, r) + V_{GxG} K_k^{(2)}(m, r) + \textstyle\sum_j V_j I_{jk} + V_{fam}\, r_f^{f_k},$$

where $K$ is the additive-genetic sharing under assortative mating at equilibrium (Fisher;
Nagylaki's half-sibling multiplier; Eftedal et al. Table S9), $K^{(2)}$ its GxG analogue
(taken as $K^2$), $I_{jk}$ indicators for shared household / cousin environments, and $f_k$ a
"family degree" for a decaying family environment. Parameters are fitted by weighted least
squares on the correlation or log-correlation scale.

Every choice is a field of `eftedal.model.FitSpec`:

| field | options |
|-------|---------|
| `components` | any of `A`, `GxG`, `house1`, `house2`, `house3`, `housepart`, `cousins1`, `cousins2`, `family` (named sets in `MODELS`) |
| `cohort` | `ALL`, `ALL_but_adopted`, `ALL_but_adopted_UZ` (Eftedal "bio + in-law"), `bio_only` |
| `scale` | `lin`, `log` |
| `weights` | `invSqErr` ($1/SE^2$ on either scale, as Eftedal et al.), `invvar` (inverse sampling variance on the scale used; $\rho^2/SE^2$ on the log scale by the delta method), `uniform`, `hybrid` / `hybrid_invvar` (equal total weight per relationship type) |
| `mr` | `free` ($m$, $r$ separate), `equilibrium` ($m = V_A r$, primary phenotypic assortment), `equal` ($r = m$) |
| `stepsib` | `eftedal` ($m((1+m)/2)^2$) or `corrected` ($m((1+r)/2)^2$) |
| `uz`, `mz_frac` | relatives through twins of unknown zygosity: MZ/DZ `mixture` (default, fraction 0.5) or the averaged degree in the `exponent` |
| `fixed`, `exclude_types`, `drop_duplicate`, `cocousin_family_degree` | fixed parameter values, leave-out sets, data checks |
| `log_zero`, `adoptive` | handling of unloggable rows and of adoptive pairs; `legacy*` values reproduce the 2025 analysis |

Uncertainty (`eftedal.inference`): sandwich and quasi (residual-variance-scaled) Wald
intervals, parametric bootstrap, and, for replication only, the inverse Hessian of the RSS.

## Running

```bash
pytest
python scripts/01_replicate_eftedal.py
python scripts/02_replicate_legacy.py
python scripts/03_sensitivity.py
python scripts/04_uncertainty.py
Rscript R/make_figures.R
```
