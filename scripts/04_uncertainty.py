"""Parameter uncertainty for the chosen models (the top cross-validation cluster plus the
additive reference) under the main settings. Writes results/uncertainty_main/:

intervals.csv    95% intervals per model x weighting x parameter:
                   param_boot x1     parametric bootstrap (rho* ~ N(rho, SE^2), refit; sampling
                                     noise only)
                   param_boot x1.3   same with SEs inflated 1.3x (pair overlap between categories)
                   block jackknife   refits leaving out each of the 10 cv blocks (sensitivity to
                                     which relationship types exist; misspecification-aware)
                   sandwich          full-Hessian sandwich Wald interval (cross-check only)
boot_reps.csv    bootstrap replicates
boot_corr.csv    correlations between parameters across bootstrap replicates (identifiability)
jackknife.csv    the leave-one-block-out estimates themselves
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from eftedal import cv
from eftedal import data as D
from eftedal.fit import fit
from eftedal.inference import parametric_bootstrap, percentile_table, wald_table
from eftedal.model import PARAM_BOUNDS
from eftedal.presets import MAIN_WEIGHTINGS, main_spec

CHOSEN = ["A", "A+GxG", "A+F", "A+GxG+F", "A+F+mat"]
REPORT = ["V_A", "V_GxG", "V_family", "r_f", "delta", "m", "r", "c"]


def main(out: Path, n_boot: int, n_boot_infl: int, n_starts: int):
    out.mkdir(parents=True, exist_ok=True)
    df = D.load("MRE8")
    blocks = {b: D.rels_in_blocks(df, b) for b in D.CV_BLOCKS}
    ints, reps, corrs, jks = [], [], [], []
    for k, (w, m) in enumerate([(w, m) for w in MAIN_WEIGHTINGS for m in CHOSEN]):
        s = main_spec(m, w)
        res = fit(df, s, n_starts=n_starts)
        est = dict(res.params)
        est["c"] = res.summary()["c"]
        est = {p: v for p, v in est.items() if p in REPORT}
        tag = dict(model=m, weights=w)

        blocks_out = [wald_table(res, methods=("sandwich",))]
        for infl, nb in ((1.0, n_boot), (1.3, n_boot_infl)):
            b = parametric_bootstrap(df, res, n_boot=nb, seed=100 + k, se_inflation=infl, n_starts=3)
            reps.append(b.assign(se_inflation=infl, **tag))
            blocks_out.append(percentile_table(b, est, method=f"param_boot x{infl:g}"))
            if infl == 1.0:
                c = b[[p for p in REPORT if p in b and b[p].std() > 0]].corr()
                c = c.stack().rename("corr").reset_index().rename(columns={"level_0": "p1", "level_1": "p2"})
                corrs.append(c[c.p1 < c.p2].assign(**tag))

        _, fp, _, _ = cv.cross_validate(df, s, n_starts=n_starts, folds=blocks)
        fp = fp.assign(**tag)
        if "m" in fp and "r" in fp:
            fp["c"] = (1 + 2 * fp["m"] + fp["r"] * fp["m"]) / (1 + fp["m"]) ** 2
        jks.append(fp)
        blocks_out.append(cv.jackknife(fp, est, [p for p in est if p in fp]))

        t = pd.concat(blocks_out, ignore_index=True)
        t = t[t["param"].isin(REPORT)].assign(**tag)
        lb = t["param"].map(lambda p: PARAM_BOUNDS.get(p, (0.0, 1.0))[0])
        ub = t["param"].map(lambda p: PARAM_BOUNDS.get(p, (0.0, 1.0))[1])
        t["spans_bound"] = ((t["lo"] <= lb + 1e-3) | (t["hi"] >= ub - 1e-3)).where(t["param"] != "c")
        ints.append(t)
        print(w, m, "done", flush=True)
    pd.concat(ints, ignore_index=True).to_csv(out / "intervals.csv", index=False)
    pd.concat(reps, ignore_index=True).to_csv(out / "boot_reps.csv", index=False)
    pd.concat(corrs, ignore_index=True).to_csv(out / "boot_corr.csv", index=False)
    pd.concat(jks, ignore_index=True).to_csv(out / "jackknife.csv", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=D.REPO_ROOT / "results" / "uncertainty_main")
    ap.add_argument("--n-boot", type=int, default=500)
    ap.add_argument("--n-boot-infl", type=int, default=200)
    ap.add_argument("--n-starts", type=int, default=20)
    a = ap.parse_args()
    main(a.out, a.n_boot, a.n_boot_infl, a.n_starts)
