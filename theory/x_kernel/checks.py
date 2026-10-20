"""Validation checks quoted in x_kernel.tex. Writes tex/checks.tex (one macro per number).

    PYTHONPATH=<worktree>/src python checks.py
"""

from pathlib import Path

import numpy as np
import pandas as pd

from eftedal.data import load
from pedigrees import configs
from x_kernel import eq_for_m, k_total

HERE = Path(__file__).parent


def kx(rel, side, sex, m, s, VA0=0.45, VE=0.45):
    e0, e1 = eq_for_m(m, VA0=VA0, s=0.0, VE=VE), eq_for_m(m, VA0=VA0, s=s, VE=VE)
    out = 0.0
    for c in configs(rel, side, sex):
        out += c.w * (k_total(c.ped, c.i, c.j, e1) - (1 - s) * k_total(c.ped, c.i, c.j, e0)) / s
    return out


def fmt(v):
    """Errors below the 6 decimals stored in the CSVs."""
    return "$<10^{-6}$" if v < 5e-7 else f"{v:.1e}"


def main():
    xk = pd.read_csv(HERE / "x_kinship.csv")
    grid = pd.read_csv(HERE / "x_kernel_grid.csv")
    g0 = grid[grid.m == 0].merge(xk[["rel", "side", "sex", "KX_FDC", "a"]], on=["rel", "side", "sex"])
    rm_err = float(max((g0.KX - g0.KX_FDC).abs().max(), (g0.Kauto - g0.a).abs().max()))
    sib = grid[(grid.rel == "sibs") & (grid.sex == "os")]
    fisher_err = float((sib.Kauto - (1 + sib.m) / 2).abs().max())
    s_vals = [kx("sibs", "other", "os", 0.43, s) for s in (0.005, 0.02, 0.05)]
    # h2 sensitivity at m = 0.43: base (V_A0, V_E) = (0.45, 0.45) vs (0.70, 0.20), all categories
    cats = load("MRE8")[["rel", "side", "sex"]].drop_duplicates()
    ref = grid[np.isclose(grid.m, 0.43)].set_index(["rel", "side", "sex"])
    h2_err = 0.0
    for _, c in cats.iterrows():
        if not configs(c.rel, c.side, c.sex):
            continue
        d_alt = kx(c.rel, c.side, c.sex, 0.43, 0.02, 0.70, 0.20)
        e0 = eq_for_m(0.43, VA0=0.70, s=0.0, VE=0.20)
        ka = sum(cf.w * k_total(cf.ped, cf.i, cf.j, e0) for cf in configs(c.rel, c.side, c.sex))
        h2_err = max(h2_err, abs((d_alt - ka) - ref.loc[(c.rel, c.side, c.sex), "D"]))
    macros = {
        "chkRM": fmt(rm_err), "chkFisher": fmt(fisher_err),
        "chkSa": f"{s_vals[0]:.4f}", "chkSb": f"{s_vals[1]:.4f}", "chkSc": f"{s_vals[2]:.4f}",
        "chkHtwo": f"{h2_err:.3f}", "chkHtwoRho": f"{0.7 * 0.02 * h2_err:.4f}",
    }
    (HERE / "tex" / "checks.tex").write_text(
        "".join(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in macros.items()))
    print(macros)


if __name__ == "__main__":
    main()
