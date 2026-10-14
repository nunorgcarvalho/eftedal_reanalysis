#!/bin/bash
source ~/alkes_nuno/.venv_py13/bin/activate
cd /n/groups/price/nuno/eftedal_reanalysis/theory/gxg_kernel
for reg in rm:-1 strong:0.6; do
  name=${reg%%:*}; noise=${reg##*:}
  for spec in "uniform 0" "beta 0" "beta_1mc 0" "beta_prox 0.001" "beta_prox 0.05" "beta_prox_dir 0.001"; do
    set -- $spec
    python simulate_linked.py --arch $1 --lam $2 --reps 3 --gens 20 --noise $noise --out results/linked/$name
  done
done
touch results/linked/all.done
