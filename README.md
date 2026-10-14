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

$$\rho_k = \big[V_A K_k(m, r) + V_{GxG} K_k^{(2)}(m, r) + V_{house}(h^{full}_k + \pi h^{part}_k) + V_{fam}\, r_f^{f_k}\big](1 + \delta M_k),$$

where $K$ is the additive-genetic sharing under assortative mating at equilibrium (Fisher;
Nagylaki's half-sibling multiplier; Eftedal et al. Table S9), $K^{(2)}$ its GxG analogue
(taken as $K^2$), a household environment shared fully by full-time and by a fraction $\pi$ by
part-time sibling pairs, a decaying family environment over the "family degree" $f_k$, and a
multiplicative maternal excess $\delta$ for categories related only through mothers. Parameters are fitted by weighted least
squares on the correlation or log-correlation scale.

Every choice is a field of `eftedal.model.FitSpec`:

| field | options |
|-------|---------|
| `components` | any of `A`, `GxG`, `house`, `family`, `maternal` (core set: `CORE_MODELS`, 8 models without the maternal term; `CORE_MODELS_WITH_MATERNAL`, 16), plus the 2025 indicator components `house1`, `house2`, `house3`, `housepart`, `cousins1`, `cousins2` (`MODELS_2025`) |
| `cohort` | `ALL`, `ALL_but_adopted`, `ALL_but_adopted_UZ` (Eftedal "bio + in-law"), `bio_only` |
| `scale` | `lin`, `log` |
| `weights` | `invSqErr` ($1/SE^2$ on either scale, as Eftedal et al.), `invvar` (inverse sampling variance on the scale used; $\rho^2/SE^2$ on the log scale by the delta method), `uniform`, `hybrid` / `hybrid_invvar` (equal total weight per relationship type), `hybrid_block` (equal weight per cv block) |
| `mr` | `free` ($m$, $r$ separate), `equilibrium` ($m = V_A r$, primary phenotypic assortment), `equal` ($r = m$) |
| `stepsib` | `eftedal` ($m((1+m)/2)^2$) or `corrected` ($m((1+r)/2)^2$) |
| `uz`, `mz_frac` | relatives through twins of unknown zygosity: MZ/DZ `mixture` (default, fraction 0.5) or the averaged degree in the `exponent` |
| `fixed`, `exclude_types`, `drop_duplicate`, `cocousin_family_degree` | fixed parameter values, leave-out sets, data checks |
| `sequential_base`, `sequential_carry` | sequential estimation: fit the base model first and hold the carried parameters (default `V_A`) fixed |
| `log_zero`, `adoptive` | handling of unloggable rows and of adoptive pairs; `legacy*` values reproduce the 2025 analysis |

Main settings (`presets.MAIN`, `presets.main_spec`): linear scale, $m = V_A r$, corrected
step-sibling kernel, co-cousin family degree 2, $\sum V \le 1$; weightings `hybrid_block` (primary) and `invvar`.

Uncertainty (`eftedal.inference`, `eftedal.cv`): parametric bootstrap, block jackknife over the
10 cv blocks, full-Hessian sandwich, and (for replication only) the inverse Hessian of the RSS.
Model comparison: block cross-validation (`eftedal.cv`).

## Running

```bash
pytest
python scripts/01_replicate_eftedal.py     # Eftedal et al. (2025), exactly
python scripts/02_fit_models.py            # 8 core models x {hybrid_block, invvar}
python scripts/03_cross_validation.py      # block cross-validation
python scripts/04_uncertainty.py           # bootstrap + block jackknife (~1 h on one core)
python scripts/05_robustness.py            # robustness rows, bio vs in-law contrast, phenotypes
python scripts/06_sequential.py            # joint vs sequential (A first) estimation
python scripts/07_profile.py               # profile over V_A
Rscript R/make_figures.R                   # -> results/figures/main/
```

`docs/codebase_overview.html` explains the model, estimation, optimisers, uncertainty and
cross-validation in detail. `scripts/exploratory/` and `R/exploratory_figures.R` hold the first,
exploratory round of analyses.
