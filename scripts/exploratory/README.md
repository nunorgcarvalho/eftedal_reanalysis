# Exploratory scripts (2026-10-12)

The first round of analyses after the rebuild: grids over fitting settings (scale x weights
x m-r handling x step-sibling kernel) for the 2025 model set, leave-one-type-out refits,
uncertainty comparisons (sandwich / quasi / inverse Hessian / bootstrap) and the additive
model across all phenotypes. They informed the analysis plan; the main pipeline is
`scripts/0*.py`. Outputs go to `results/sensitivity`, `results/uncertainty` and
`results/phenotypes`.
