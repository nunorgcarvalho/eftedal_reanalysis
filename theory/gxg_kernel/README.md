# Is K^2 the right GxG kernel under assortative mating?

Report: `gxg_kernel.pdf` (source `gxg_kernel.tex`). Short answer: K^2 is exact for
*scale* epistasis (GxG as a quadratic function of the additive value) when the curvature is weak.
It is wrong for *random* pairwise epistasis: there the GxG kernel stays near a^2, the squared
pedigree relatedness, plus a small term from sorting on GxG itself. Linked (cis) interactions
are the exception: a tightly linked pair is one haplotype locus, with kernel (a + delta)/2 under
random mating, and AM scales it further (section "Linked interacting loci").

| file | role |
|---|---|
| `simulate.py` | allele-level forward simulation: AM on Y = A + GxG + E; relative classes; GxG kernels |
| `ag_path_model.py` | linearized GxG transmission as a pathMgr path model (equilibrium + per-class predictions) |
| `compare.py` | simulation vs predictions -> `results/comparison.csv` |
| `realistic.py` | predictions at illustrative test-score-like parameters -> `results/realistic.csv` |
| `make_tables.py`, `figure.R` | tables and figure included by the report |
| `digenic.py` | linkage, random mating: exact digenic coefficients psi(c) by gene dropping |
| `ag_linked_model.py` | linkage + AM: cis/trans GxG per recombination bin as a pathMgr model |
| `simulate_linked.py` | haplotype-level simulation with a genetic map and five interaction architectures |
| `compare_linked.py`, `realistic_linked.py`, `make_tables_linked.py`, `figure_linked.R` | linkage section outputs (`results/linked/`) |

Reproduce (one core, ~1 h):
```bash
python simulate.py --arch random  --reps 3 --gens 20 --out results/strong
python simulate.py --arch random  --reps 3 --gens 20 --noise 1.2 --out results/moderate
python simulate.py --arch aligned --reps 3 --gens 20 --VG 0.0001 --VE 0.5499 --out results/strong
python simulate.py --arch aligned --reps 3 --gens 20 --VG 0.0001 --VE 0.5499 --noise 1.2 --out results/moderate
python simulate.py --arch aligned --reps 3 --gens 20 --VG 0.03 --VE 0.52 --out results/strong_curved
python simulate.py --arch aligned --reps 3 --gens 20 --VG 0.03 --VE 0.52 --noise 1.2 --out results/moderate_curved
python compare.py && python realistic.py && python make_tables.py && Rscript figure.R
# linkage section (one core, ~1.5 h): run_linked.sh runs the 12 simulate_linked.py runs
bash run_linked.sh
python digenic.py && python realistic_linked.py && python compare_linked.py
python make_tables_linked.py && Rscript figure_linked.R
pdflatex gxg_kernel.tex   # twice
```
`results/` is not tracked; the report PDF is.
