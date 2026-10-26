# eftedal_reanalysis

A reanalysis of the relative-pair correlations published by Eftedal et al. (2025,
*PNAS*) for Norwegian national test scores. Each of the 82 categories of relatives has
one observed phenotypic correlation. We fit models that predict those correlations
from genetic relatedness under assortative mating, with optional nonadditive (GxG) and
shared-environment components, and ask how sensitive the conclusions are to modelling
and fitting choices.

## ASHG 2026 poster and supplemental slides

> **Here from the QR code on the poster?**
>
> ### [Poster (PDF)](https://github.com/nunorgcarvalho/eftedal_reanalysis/blob/main/ASHG/poster/poster.pdf)
> ### [Supplemental slides (PDF)](https://github.com/nunorgcarvalho/eftedal_reanalysis/blob/main/ASHG/supplemental_slides/supplemental_slides.pdf)
>
> *Relatives' similarity in test score performance is explained by a wide range of genetic and
> environmental parameters.* Nuno R. G. Carvalho and Alexander Gusev, ASHG 2026 Annual Meeting,
> Montréal. The supplemental slides add a model-by-model view of every fit, the cross-validation
> results, and clarifications to the poster. All ASHG materials are in [`ASHG/`](ASHG/).

## AI Disclosure
I used Claude (Opus 5.5) extensively to help me with this project, and almost
everything in this repository was ultimately written by Claude.
However, I was responsible for deciding which analyses to run and was constantly
in the loop regarding their implementation and provided direct feedback about
improving it.
With the exception of this disclosure, this README was written by Claude, as well
as the figures, tables, poster(s), slides, notebooks, and more.
However, the text on the poster was written by me.

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

The predicted correlation for relationship category $k$ is a sum of components

$$\rho_k = V_A\, K_k(m, r) + V_{GxG}\, (1/4)^{n_k} + V_C\, c^{f_k},$$

- **A**, additive genetic effects under assortative mating at equilibrium: $K_k$ is the
  additive-genetic sharing of the pair (Fisher; Nagylaki's half-sibling multiplier; Eftedal et al.
  Table S9), with mates' genetic correlation $m = V_A r$ for a phenotypic correlation $r$;
- **GxG**, pairwise gene-gene interactions: $(1/4)^{n_k}$ is the square of the random-mating
  kernel $(1/2)^{n_k}$ for genetic degree $n_k$, not affected by assortative mating (0 for in-laws
  and adoptees);
- **C**, shared environment, decaying by a factor $c$ per "family degree" $f_k$. It splits into a
  within-sibling part $V_C(1-c)$ and a beyond-sibling part $V_C\, c$.

The nine core models (`CORE_MODELS`) are A, A→GxG, A+GxG, A→C, A+C, A+GxG→C, A+C→GxG, A→GxG+C
and A+GxG+C: "+" means fitted jointly, "→" fitted in steps, with the earlier components'
parameters held fixed. In the code, C is the `family` component (`V_family`, decay `r_f` = $c$),
and model names use F for it.

**Phase 2 (current; after ASHG 2026).** Every core and joint model also fits the sex-composition
and family-side multipliers its components allow (`FitSpec.modifiers`):
- **ρ_g** scales opposite-sex pairs' genetic covariance (A × ρ_g, GxG × ρ_g²).
- **ρ_e** scales opposite-sex pairs' shared environment C.
- **κ** redistributes C between the maternal (1 + κ) and paternal (1 − κ) side of a type.

Notes:
- **Two-step models** estimate each multiplier in the step where its component enters.
- **Phase 1** (no multipliers) is the `core_p1` set. The Additive model of the Eftedal-vs-revised
  comparison stays Phase 1.
- **Jackknife CIs** use t₉ (10 blocks).
- **X-linked kernel:** a sensitivity analysis (`x_share`, set `core_x`), derived in
  `theory/x_kernel/`.
- **ASHG 2026 materials** were Phase 1 with 1.96 intervals, and are frozen at tag `ASHG2026` (see
  `ASHG/README.md`).

Parameters are fitted by weighted least squares on the correlation scale, in two ways:

| | Marginal Pair (`invvar`) | Marginal Type (`hybrid_block`) |
|---|---|---|
| unit | each pair of relatives | each of 10 relative types |
| weights | $1/SE^2$ | equal total weight per relative type; $1/SE^2$ within a type |
| 95% CI | parametric bootstrap | leave-one-type-out jackknife |

Every analysis choice is a field of `eftedal.model.FitSpec`:

| field | options |
|-------|---------|
| `components` | `A`, `GxG`, `family` (C) |
| `cohort` | `ALL`, `ALL_but_adopted`, `ALL_but_adopted_UZ` (Eftedal "bio + in-law"), `bio_only` |
| `scale` | `lin`, `log` |
| `weights` | `invSqErr` ($1/SE^2$ on either scale, as Eftedal et al.), `invvar` (inverse sampling variance on the scale used; $\rho^2/SE^2$ on the log scale by the delta method), `hybrid_block` (Marginal Type) |
| `mr` | `free` ($m$, $r$ separate, as Eftedal et al.), `equilibrium` ($m = V_A r$) |
| `stepsib` | `eftedal` ($m((1+m)/2)^2$) or `corrected` ($m((1+r)/2)^2$) |
| `uz`, `mz_frac` | relatives through twins of unknown zygosity: an MZ/DZ `mixture` (fraction 0.5) |
| `sequential_base`, `sequential_carry` | two-step ("→") estimation |

Further fields, used for replication and sensitivity checks, are documented in `FitSpec`.

Main settings (`presets.MAIN`, `presets.main_spec`): linear scale, $m = V_A r$, corrected
step-sibling kernel, GxG kernel $(1/4)^n$, $\sum V \le 1$; fit types Marginal Type (primary) and
Marginal Pair.

Uncertainty (`eftedal.inference`, `eftedal.cv`): parametric bootstrap, block jackknife over the
10 relative types, and (for replication only) the inverse Hessian of the RSS. Model comparison:
leave-one-type-out cross-validation (`eftedal.cv`).

## Running

```bash
pytest
python scripts/01_replicate_eftedal.py     # Eftedal et al. (2025), exactly
python scripts/02_corrections.py           # Eftedal -> ours step by step; bio-only vs in-laws
python scripts/03_models.py                # core models: fits, cross-validation, jackknife
python scripts/04_bootstrap.py             # parametric bootstrap of the core models (~1 h)
python scripts/05_robustness.py            # robustness rows; in-law contrast for every test
python scripts/06_profile.py               # profile over V_A
python scripts/07_decay.py                 # component decay across relative types (no fitting)
python scripts/08_param_corr.py            # parameter sampling correlations (needs 03)
Rscript R/make_figures.R                   # -> results/figures/

# ASHG 2026 poster, supplemental slides and Q&A notes (see ASHG/README.md)
python scripts/poster_additive.py          # additive model by estimator and cohort
python scripts/poster_eftedal_predictions.py
Rscript R/poster_figures.R                 # -> results/figures/poster/
python ASHG/QnA/make_tables.py             # -> results/qna/
bash ASHG/build.sh                         # -> the PDFs in ASHG/
```

`docs/codebase_overview.html` explains the model, estimation, optimisers, uncertainty and
cross-validation in detail. `scripts/exploratory/` and `R/exploratory_figures.R` hold the first,
exploratory round of analyses.
