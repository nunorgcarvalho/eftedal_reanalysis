"""Supplement: parametric-bootstrap intervals for the core models (sampling noise only:
rho* ~ N(rho, SE^2), every step refitted), to set against the block-jackknife intervals of
scripts/03_models.py. Writes results/bootstrap/:

intervals.csv    95% percentile intervals per model x weighting x parameter (incl. V_sum)
reps.csv         the replicates
corr.csv         correlations between parameters across replicates (identifiability)
"""

import argparse
from pathlib import Path

import pandas as pd

from eftedal import data as D
from eftedal.fit import fit
from eftedal.inference import parametric_bootstrap, percentile_table
from eftedal.model import CORE_MODELS
from eftedal.presets import MAIN_WEIGHTINGS, main_spec

PARAMS = ["V_A", "V_GxG", "V_house", "V_family", "V_sum", "pi", "r_f", "m", "r"]


def main(out: Path, n_boot: int, n_starts: int):
    out.mkdir(parents=True, exist_ok=True)
    df = D.load("MRE8")
    ints, reps, corrs = [], [], []
    for k, (w, m) in enumerate([(w, m) for w in MAIN_WEIGHTINGS for m in CORE_MODELS]):
        res = fit(df, main_spec(m, w), n_starts=n_starts)
        b = parametric_bootstrap(df, res, n_boot=n_boot, seed=100 + k, n_starts=3)
        b["V_sum"] = b[[c for c in b if c.startswith("V_")]].sum(axis=1)
        est = {p: v for p, v in res.params.items() if p in PARAMS}
        est["V_sum"] = res.summary()["V_sum"]
        tag = dict(model=m, weights=w)
        ints.append(percentile_table(b, est, method="param_boot").assign(**tag))
        reps.append(b.assign(**tag))
        c = b[[p for p in PARAMS if p in b and p != "V_sum" and b[p].std() > 0]].corr()
        c = c.stack().rename("corr").reset_index().rename(columns={"level_0": "p1", "level_1": "p2"})
        corrs.append(c[c.p1 < c.p2].assign(**tag))
        print(w, m, flush=True)
    pd.concat(ints, ignore_index=True).to_csv(out / "intervals.csv", index=False)
    pd.concat(reps, ignore_index=True).to_csv(out / "reps.csv", index=False)
    pd.concat(corrs, ignore_index=True).to_csv(out / "corr.csv", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=D.REPO_ROOT / "results" / "bootstrap")
    ap.add_argument("--n-boot", type=int, default=400)
    ap.add_argument("--n-starts", type=int, default=12)
    a = ap.parse_args()
    main(a.out, a.n_boot, a.n_starts)
