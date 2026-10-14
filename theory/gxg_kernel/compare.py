"""Compare simulated relative kernels with the predictions.

For each regime directory under results/ and each architecture, per replicate: take the
simulation's measured base variances (generation 0) and mean realized mate correlation, solve the
equilibrium and compute the path-model predictions (ag_path_model.py), then summarize across
replicates. Writes results/comparison.csv and results/equilibrium.csv.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from ag_path_model import solve_equilibrium, predict

R = Path(__file__).parent / "results"
V_E = 0.25            # simulate.py default


def main():
    comp, eqrows = [], []
    for regime in sorted(p for p in R.iterdir() if p.is_dir()):
        for arch in ("random", "aligned"):
            fc, fp = regime / f"{arch}_classes.csv", regime / f"{arch}_population.csv"
            if not fc.exists():
                continue
            cls, pop = pd.read_csv(fc), pd.read_csv(fp)
            for rep, pr in pop.groupby("rep"):
                g0 = pr[pr.gen == 0].iloc[0]
                last = pr[pr.gen == pr.gen.max()].iloc[0]
                rho = pr["rho_y"].mean()
                eq = solve_equilibrium(g0.V_A, g0.V_G, V_E, rho)
                eqrows.append(dict(regime=regime.name, arch=arch, rep=rep, rho_y=rho,
                                   V_A_sim=last.V_A, V_A_pred=eq.V_A, V_G_sim=last.V_G,
                                   V_G_pred=eq.V_G, c_sim=last.c_AG, c_pred=eq.c))
                pred = predict(eq).rename(columns={"K_A": "K_A_path", "K_G": "K_G_path",
                                                   "cross": "cross_path"})
                sim = cls[cls.rep == rep]
                comp.append(sim.merge(pred, on="cls").assign(regime=regime.name, arch=arch))
            print(regime.name, arch, flush=True)
    c = pd.concat(comp, ignore_index=True)
    c["K_A_sim_sq"] = c["K_A"] ** 2
    c.to_csv(R / "comparison_reps.csv", index=False)
    agg = (c.groupby(["regime", "arch", "cls"])
             .agg(n_reps=("rep", "nunique"), n_pairs=("n_pairs", "mean"),
                  K_A_sim=("K_A", "mean"), K_A_path=("K_A_path", "mean"),
                  K_G_sim=("K_G", "mean"), K_G_sim_se=("K_G", lambda x: x.std(ddof=1) / np.sqrt(len(x))),
                  K_G_path=("K_G_path", "mean"), K_A_sq=("K_A_sim_sq", "mean"), a_sq=("a_sq", "mean"),
                  cross_sim=("cross", "mean"), cross_path=("cross_path", "mean"))
             .reset_index())
    agg.to_csv(R / "comparison.csv", index=False)
    pd.DataFrame(eqrows).to_csv(R / "equilibrium.csv", index=False)
    with pd.option_context("display.width", 250):
        print(agg.round(4).to_string(index=False))
        print(pd.DataFrame(eqrows).round(4).to_string(index=False))


if __name__ == "__main__":
    main()
