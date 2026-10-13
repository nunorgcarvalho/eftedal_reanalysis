"""Uncertainty analyses for the key specs. Writes results/uncertainty/:

cis.csv              95% intervals per key spec and parameter: Wald (sandwich, quasi,
                     inverse Hessian of the RSS as in Eftedal et al.'s code) and parametric
                     bootstrap (with and without inflating the SEs by 1.3 for pair overlap)
boot_reps.csv        bootstrap replicates (every parameter)
hessian_scaling.csv  the Hessian-of-RSS 'SE' of the same estimate when the weights are
                     multiplied by a constant (which leaves the estimate unchanged)
profile_r.csv        loss profile over a grid of fixed r (other parameters refitted)
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from eftedal import data as D
from eftedal.fit import fit
from eftedal.inference import (covariance, parametric_bootstrap, percentile_table,
                               wald_table)
from eftedal.model import Model
from eftedal.presets import key_specs, spec


def intervals(df, n_boot, n_starts):
    cis, reps = [], []
    for k, s in enumerate(key_specs()):
        res = fit(df, s, n_starts=n_starts)
        w = wald_table(res)
        est = {p: res.params[p] for p in res.params}
        est["c"] = res.summary()["c"]
        blocks = [w]
        for infl in (1.0, 1.3):
            b = parametric_bootstrap(df, res, n_boot=n_boot, seed=k, se_inflation=infl)
            b["se_inflation"] = infl
            reps.append(b.assign(key_id=k, key_label=s.label))
            blocks.append(percentile_table(b, est, method=f"param_boot x{infl:g}"))
        t = pd.concat(blocks, ignore_index=True)
        t.insert(0, "key_label", s.label)
        t.insert(0, "key_id", k)
        t["chi2_per_df"] = res.stats["chi2_per_df"]
        cis.append(t)
        print(s.label, "done", flush=True)
    return pd.concat(cis, ignore_index=True), pd.concat(reps, ignore_index=True)


def hessian_scaling(df):
    """Same estimate, weights multiplied by a constant: the inverse-Hessian 'SE' scales by
    1/sqrt(constant), so its value depends on an arbitrary normalisation; the sandwich
    does not."""
    rows = []
    for s in (key_specs()[0], key_specs()[4]):
        res = fit(df, s)
        w0 = res.model.w.copy()
        for scale_name, factor in [("as used", 1.0), ("sum to 1", 1 / w0.sum()),
                                   ("mean 1", len(w0) / w0.sum()), ("x1e6", 1e6), ("x1e-6", 1e-6)]:
            res.model.w = w0 * factor
            for meth in ("hessian_rss", "sandwich"):
                se = np.sqrt(np.clip(np.diag(covariance(res, meth)), 0, None))
                for p, v, e in zip(res.free, res.theta, se):
                    rows.append({"spec": s.label, "weights": scale_name, "factor": factor,
                                 "method": meth, "param": p, "estimate": v, "se": e})
        res.model.w = w0
    return pd.DataFrame(rows)


def profile_r(df, grid=np.round(np.linspace(0, 1, 21), 3)):
    specs = [spec("M01", "ALL_but_adopted_UZ", st) for st in ("log-invSqErr", "log-invvar", "lin-invvar")]
    specs += [spec("M01", "ALL_but_adopted_UZ", "log-invvar", stepsib="corrected")]
    rows = []
    for s in specs:
        best = fit(df, s)
        for r in grid:
            res = fit(df, s.with_(fixed=(("r", float(r)),)), n_starts=4)
            rows.append({"spec": s.label, "r": r, "loss": res.loss, "loss_best": best.loss,
                         "df_best": best.stats["n_rows_in_loss"] - len(best.free),
                         "V_A": res.params["V_A"], "m": res.params["m"],
                         "implied_r_equilibrium": res.params["m"] / res.params["V_A"]})
    return pd.DataFrame(rows)


def main(out: Path, n_boot: int, n_starts: int):
    out.mkdir(parents=True, exist_ok=True)
    df = D.load("MRE8")
    hessian_scaling(df).to_csv(out / "hessian_scaling.csv", index=False)
    profile_r(df).to_csv(out / "profile_r.csv", index=False)
    cis, reps = intervals(df, n_boot, n_starts)
    cis.to_csv(out / "cis.csv", index=False)
    reps.to_csv(out / "boot_reps.csv", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=D.REPO_ROOT / "results" / "uncertainty")
    ap.add_argument("--n-boot", type=int, default=500)
    ap.add_argument("--n-starts", type=int, default=12)
    a = ap.parse_args()
    main(a.out, a.n_boot, a.n_starts)
