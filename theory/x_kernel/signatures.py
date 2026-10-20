"""Predicted X-only signatures vs the observed MRE8 correlations.

X adds V_A * s_X * D_k(m) to the correlation of category k (D = K^X - K^auto, x_kernel_grid.csv).
Contrasts, computed identically for predictions and data:
    ss - os            per type and side: mean(ssf, ssm) - os; then the mean over sides
    maternal - paternal per type and sex composition: mat - pat; then the mean over sex compositions
Observed standard errors assume independent category estimates.
Writes signatures.csv and prints a summary.
"""

from pathlib import Path

import warnings

import numpy as np
import pandas as pd

from eftedal.data import load

HERE = Path(__file__).parent
warnings.simplefilter("ignore", pd.errors.PerformanceWarning)
TYPES_SS = ["dz", "sibs", "half sibs", "1st cousins", "half cousins", "2nd cousins", "step sibs", "cocousins"]
TYPES_MP = ["half sibs", "1st cousins", "half cousins", "2nd cousins", "1st cousins t/uz", "cocousins"]


def contrasts(df, val, var):
    """df: one row per category with columns val (value) and var (sampling variance)."""
    out = []
    g = df.set_index(["rel", "side", "sex"])
    for rel in TYPES_SS:
        sides = sorted(df[df.rel == rel].side.unique())
        ests, vs = [], []
        for sd in sides:
            r = g.loc[(rel, sd)]
            ests.append((r.loc["ssf", val] + r.loc["ssm", val]) / 2 - r.loc["os", val])
            vs.append((r.loc["ssf", var] + r.loc["ssm", var]) / 4 + r.loc["os", var])
        out.append(dict(contrast="ss - os", rel=rel, est=np.mean(ests), se=np.sqrt(np.sum(vs)) / len(sides)))
    for rel in TYPES_MP:
        ests, vs = [], []
        for sx in ("ssf", "ssm", "os"):
            ests.append(g.loc[(rel, "maternal", sx), val] - g.loc[(rel, "paternal", sx), val])
            vs.append(g.loc[(rel, "maternal", sx), var] + g.loc[(rel, "paternal", sx), var])
        out.append(dict(contrast="maternal - paternal", rel=rel, est=np.mean(ests), se=np.sqrt(np.sum(vs)) / 3))
    return pd.DataFrame(out)


def main():
    d = load("MRE8")
    d = d[~d.suspect_duplicate.astype(bool)] if False else d       # keep all rows (no 2nd-cousin 'other' contrast uses it alone)
    obs = contrasts(d.assign(v=d.se ** 2), "corr", "v").rename(columns={"est": "observed", "se": "obs_se"})
    grid = pd.read_csv(HERE / "x_kernel_grid.csv")
    rows = [obs]
    for m in (0.0, 0.43):
        gm = grid[np.isclose(grid.m, m)].assign(v=0.0)
        p = contrasts(gm, "D", "v")[["contrast", "rel", "est"]]
        rows.append(p.rename(columns={"est": f"D_contrast_m{m}"}))
    xk = pd.read_csv(HERE / "x_kinship.csv").assign(v=0.0)
    rows.append(contrasts(xk.assign(D=xk.KX_NDC - xk.a), "D", "v")[["contrast", "rel", "est"]]
                .rename(columns={"est": "D_contrast_NDC_m0.0"}))
    out = rows[0]
    for r in rows[1:]:
        out = out.merge(r, on=["contrast", "rel"])
    for VA in (0.4, 0.7):
        for s in (0.02, 0.05):
            out[f"pred_VA{VA}_s{s}"] = VA * s * out["D_contrast_m0.43"]
    out.round(5).to_csv(HERE / "signatures.csv", index=False)
    with pd.option_context("display.width", 220):
        print(out.round(4).to_string(index=False))
    # per sex composition, half sibs and 1st cousins: maternal - paternal
    g = grid[np.isclose(grid.m, 0.43)].set_index(["rel", "side", "sex"])
    o = d.set_index(["rel", "side", "sex"])
    print("\nmaternal - paternal by sex composition (observed; predicted at V_A 0.5, s 0.02, m 0.43):")
    for rel in ("half sibs", "1st cousins", "2nd cousins", "half cousins"):
        for sx in ("ssf", "ssm", "os"):
            ob = o.loc[(rel, "maternal", sx), "corr"] - o.loc[(rel, "paternal", sx), "corr"]
            se = np.hypot(o.loc[(rel, "maternal", sx), "se"], o.loc[(rel, "paternal", sx), "se"])
            pr = 0.5 * 0.02 * (g.loc[(rel, "maternal", sx), "D"] - g.loc[(rel, "paternal", sx), "D"])
            print(f"  {rel:13s} {sx}: observed {ob:+.4f} (SE {se:.4f})   X predicted {pr:+.4f}")


if __name__ == "__main__":
    main()
