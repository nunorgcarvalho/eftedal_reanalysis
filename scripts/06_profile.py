"""Profile over V_A: hold V_A at each value of a grid, refit every other parameter, and record
the loss and the refitted components; the joint fit and the two-step core model (additive part
frozen at the A-only fit) are marked. Shows how flat the loss is along the direction in which
V_A trades off against V_GxG / V_family. Writes results/profile/profile_VA.csv with, per point:

loss_rel   loss / minimum loss over the grid and the joint fit
qlr        quasi-likelihood ratio (loss - min) / (min / df): the loss increase in units of the
           residual variance; ~chi2_1 under a quasi-likelihood view, so 3.84 is a rough 95% cut-off
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from eftedal import data as D
from eftedal.fit import fit
from eftedal.presets import MAIN_WEIGHTINGS, main_spec

# joint model -> its two-step core counterpart (None: no two-step version)
MODELS = {"A+GxG": "A>GxG", "A+GxG+H": "A>GxG+H", "A+GxG+F": "A>GxG+F", "A+F": None}
GRID = np.round(np.arange(0.20, 0.901, 0.025), 3)


def main(out: Path, n_starts: int):
    out.mkdir(parents=True, exist_ok=True)
    df = D.load("MRE8")
    rows = []
    for w in MAIN_WEIGHTINGS:
        for m, two_step in MODELS.items():
            joint = fit(df, main_spec(m, w), n_starts=n_starts)
            seq = fit(df, main_spec(two_step, w), n_starts=n_starts) if two_step else None
            pts = [("grid", v) for v in GRID] + [("joint", joint.params["V_A"])]
            if seq is not None:
                pts.append(("two-step", seq.params["V_A"]))
            for kind, v in pts:
                res = joint if kind == "joint" else seq if kind == "two-step" else fit(
                    df, main_spec(m, w).with_(fixed=(("V_A", float(v)),)), n_starts=n_starts,
                    starts=[np.delete(joint.theta, joint.free.index("V_A"))])
                d = {k: res.params.get(k) for k in ("V_A", "V_GxG", "V_family", "r_f", "m", "r")}
                d.update(model=m, weights=w, point=kind, loss=res.loss,
                         df_resid=res.stats["n_rows_in_loss"] - len(joint.free),
                         sum_V=sum(x for k, x in res.params.items() if k.startswith("V_")))
                rows.append(d)
            print(w, m, flush=True)
    p = pd.DataFrame(rows)
    best = p.groupby(["model", "weights"])["loss"].transform("min")
    p["loss_rel"] = p["loss"] / best
    p["qlr"] = (p["loss"] - best) / (best / p["df_resid"])
    p.to_csv(out / "profile_VA.csv", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=D.REPO_ROOT / "results" / "profile")
    ap.add_argument("--n-starts", type=int, default=8)
    a = ap.parse_args()
    main(a.out, a.n_starts)
