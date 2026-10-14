"""Linked-loci simulations vs predictions -> results/linked/comparison.csv.

Per regime x architecture (bins and mate correlation averaged over replicates; base variances
from generation 0, which is at linkage equilibrium):
    a2         squared pedigree relatedness (unlinked, random mating)
    K_A_sq     simulated additive kernel squared (the current K^2 kernel)
    linked_rm  exact random-mating prediction with linkage, sum_b V_b psi(c_b) / V_G
               (= ag_linked_model at rho = 0; identical to digenic.py)
    linked_am  ag_linked_model at the measured rho (linkage + AM), equilibrium
    linked_am_T  the same after the simulated number of AM generations (founders of the class
               pedigrees at generation T-2); tight pairs approach equilibrium on a 1/c time scale
    unlinked_am  ag_path_model at the measured rho (AM, ignoring linkage)
"""

from pathlib import Path

import numpy as np
import pandas as pd

import ag_linked_model as alm
import ag_path_model as apm

R = Path(__file__).parent / "results" / "linked"
V_E = 0.25
ORDER = ["parent-offspring", "full sibs", "half sibs", "grandparent", "avuncular", "first cousins",
         "step sibs", "co-cousins"]


def bins_for(fb, V_G_gen0):
    b = pd.read_csv(fb)
    b = b[b.V0 > 0].groupby("bin").agg(c=("c", "mean"), V0=("V0", "mean")).reset_index()
    b = b[b.V0 > 1e-4 * b.V0.sum()]
    return b.c.to_numpy(), b.V0.to_numpy() * V_G_gen0 / b.V0.sum()


def main():
    comp, eqs = [], []
    for reg in sorted(p for p in R.iterdir() if p.is_dir()):
        for fc in sorted(reg.glob("*_classes.csv")):
            tag = fc.name.removesuffix("_classes.csv")
            cls = pd.read_csv(fc)
            pop = pd.read_csv(reg / f"{tag}_population.csv")
            g0 = pop[pop.gen == 0].mean()
            last = pop[pop.gen == pop.gen.max()].mean()
            rho = 0.0 if reg.name == "rm" else pop.rho_y.mean()
            c, V0 = bins_for(reg / f"{tag}_bins.csv", g0.V_G)
            eq_am = alm.solve_equilibrium(g0.V_A, c, V0, V_E, rho)
            eq_rm = alm.solve_equilibrium(g0.V_A, c, V0, V_E, 0.0)
            T = int(pop.gen.max()) + 1                          # generations of mating simulated
            eq_T = alm.solve_equilibrium(g0.V_A, c, V0, V_E, rho, iters=T - 2)
            eq_un = apm.solve_equilibrium(g0.V_A, g0.V_G, V_E, rho)
            p_am = alm.predict(eq_am, ORDER).set_index("cls")
            p_rm = alm.predict(eq_rm, ORDER).set_index("cls")
            p_T = p_rm if rho == 0 else alm.predict(eq_T, ORDER).set_index("cls")
            p_un = apm.predict(eq_un, ORDER).set_index("cls")
            agg = (cls.groupby("cls").agg(n_pairs=("n_pairs", "mean"), K_A=("K_A", "mean"),
                                          K_G=("K_G", "mean"),
                                          K_G_se=("K_G", lambda x: x.std(ddof=1) / np.sqrt(len(x))),
                                          cross=("cross", "mean")))
            agg = agg[agg.index.isin(ORDER)]
            agg["K_A_sq"] = (cls.assign(k2=cls.K_A ** 2).groupby("cls").k2.mean())
            agg["a2"] = [apm.IBD[k] ** 2 for k in agg.index]
            agg["linked_rm"] = p_rm.K_G
            agg["linked_am"] = p_am.K_G
            agg["linked_am_T"] = p_T.K_G
            agg["unlinked_am"] = p_un.K_G
            agg["cross_linked_am"] = p_am.cross
            comp.append(agg.reset_index().assign(regime=reg.name, arch=tag))
            share_tight = V0[c < 0.01].sum() / V0.sum()
            eqs.append(dict(regime=reg.name, arch=tag, rho_y=rho, V_A_sim=last.V_A, V_A_pred=eq_am.V_A, V_A_pred_T=eq_T.V_A,
                            V_G_sim=last.V_G, V_G_pred=eq_am.V_G, V_G_pred_T=eq_T.V_G, V_G_gen0=g0.V_G,
                            share_tight=share_tight,
                            share_unlinked=V0[c > 0.4999].sum() / V0.sum(),
                            skew_A=last.skew_A))
            print(reg.name, tag, flush=True)
    out = pd.concat(comp, ignore_index=True)
    out.to_csv(R / "comparison.csv", index=False)
    pd.DataFrame(eqs).to_csv(R / "equilibrium.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_rows", 500):
        print(out.round(4).to_string(index=False))
        print(pd.DataFrame(eqs).round(4).to_string(index=False))


if __name__ == "__main__":
    main()
