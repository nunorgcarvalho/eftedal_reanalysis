# Is K^2 the right GxG kernel under assortative mating?

Report: `gxg_kernel.pdf` (source `gxg_kernel.tex`). Short answer: K^2 is exact for
*scale* epistasis (GxG as a quadratic function of the additive value) when the curvature is weak.
It is wrong for *random* pairwise epistasis: there the GxG kernel stays near a^2, the squared
pedigree relatedness, plus a small term from sorting on GxG itself.

| file | role |
|---|---|
| `simulate.py` | allele-level forward simulation: AM on Y = A + GxG + E; relative classes; GxG kernels |
| `ag_path_model.py` | linearized GxG transmission as a pathMgr path model (equilibrium + per-class predictions) |
| `compare.py` | simulation vs predictions -> `results/comparison.csv` |
| `realistic.py` | predictions at illustrative test-score-like parameters -> `results/realistic.csv` |
| `make_tables.py`, `figure.R` | tables and figure included by the report |

Reproduce (one core, ~1 h):
```bash
python simulate.py --arch random  --reps 3 --gens 20 --out results/strong
python simulate.py --arch random  --reps 3 --gens 20 --noise 1.2 --out results/moderate
python simulate.py --arch aligned --reps 3 --gens 20 --VG 0.0001 --VE 0.5499 --out results/strong
python simulate.py --arch aligned --reps 3 --gens 20 --VG 0.0001 --VE 0.5499 --noise 1.2 --out results/moderate
python simulate.py --arch aligned --reps 3 --gens 20 --VG 0.03 --VE 0.52 --out results/strong_curved
python simulate.py --arch aligned --reps 3 --gens 20 --VG 0.03 --VE 0.52 --noise 1.2 --out results/moderate_curved
python compare.py && python realistic.py && python make_tables.py && Rscript figure.R
pdflatex gxg_kernel.tex   # twice
```
`results/` is not tracked; the report PDF is.
