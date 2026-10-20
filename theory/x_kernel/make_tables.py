"""LaTeX tables and sentences for x_kernel.tex, from the CSVs written by compute.py, grid.py and
signatures.py. Writes tex/*.tex.

    python make_tables.py
"""

from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).parent
OUT = HERE / "tex"
OUT.mkdir(exist_ok=True)

TYPES = ["mz", "uz", "dz", "sibs", "half sibs", "1st cousins t/uz", "1st cousins", "half cousins",
         "2nd cousins t/uz", "2nd cousins"]
SIDES = ["maternal", "paternal", "other"]
f3 = lambda v: f"{v:.3f}".replace("-0.000", "0.000")


def kernel_table(xk):
    rows = []
    for rel in TYPES:
        d = xk[xk.rel == rel]
        for side in [s for s in SIDES if s in set(d.side)]:
            g = d[d.side == side].set_index("sex")
            cell = lambda col, sx: f3(g.loc[sx, col]) if sx in g.index else "--"
            label = rel if side == "other" and rel in ("mz", "uz", "dz", "sibs") else f"{rel}, {side}"
            rows.append(" & ".join([label, f3(g.iloc[0].a)]
                                   + [cell("KX_FDC", s) for s in ("ssf", "ssm", "os")]
                                   + [cell("KX_NDC", s) for s in ("ssf", "ssm", "os")]
                                   + [f"{g.iloc[0].mt_share:.2f}"]) + r" \\")
    (OUT / "tab_kernels.tex").write_text(
        r"\begin{tabular}{lr|rrr|rrr|r}\toprule" "\n"
        r" & & \multicolumn{3}{c|}{FDC (default)} & \multicolumn{3}{c|}{NDC} & \\" "\n"
        r"category & $a$ & ssf & ssm & os & ssf & ssm & os & mt\\\midrule" "\n"
        + "\n".join(rows) + "\n" r"\bottomrule\end{tabular}" "\n")


def am_table(grid):
    sel = [("sibs", "other"), ("half sibs", "maternal"), ("half sibs", "paternal"),
           ("1st cousins", "maternal"), ("1st cousins", "paternal"), ("1st cousins", "other"),
           ("2nd cousins", "maternal"), ("2nd cousins", "paternal"),
           ("step sibs", "other"), ("cocousins", "maternal"), ("cocousins", "paternal")]
    rows = []
    for rel, side in sel:
        for sx in ("ssf", "ssm", "os"):
            g = grid[(grid.rel == rel) & (grid.side == side) & (grid.sex == sx)].set_index("m")
            rows.append(f"{rel}, {side} & {sx} & " + " & ".join(
                f"{f3(g.loc[m, 'KX'])} & {f3(g.loc[m, 'D'])}" for m in (0.0, 0.2, 0.43)) + r" \\")
    (OUT / "tab_am.tex").write_text(
        r"\begin{tabular}{ll|rr|rr|rr}\toprule" "\n"
        r" & & \multicolumn{2}{c|}{$m=0$} & \multicolumn{2}{c|}{$m=0.2$} & \multicolumn{2}{c}{$m=0.43$}\\" "\n"
        r"category & sex & $K^X$ & $D$ & $K^X$ & $D$ & $K^X$ & $D$\\\midrule" "\n"
        + "\n".join(rows) + "\n" r"\bottomrule\end{tabular}" "\n")


def approx_table(xk):
    bio = xk[xk.a > 0]
    rows = []
    for m in (0.2, 0.43):
        for name, col in (("constant $D$: $K^X(0)+K^{auto}(m)-a$", f"KX_add_m{m}"),
                          ("AM factor: $K^X(0)\\,K^{auto}(m)/a$", f"KX_mult_m{m}")):
            e = (bio[col] - bio[f"KX_m{m}"]).abs()
            rows.append(f"{m} & {name} & {e.median():.3f} & {e.max():.3f} & {0.5 * 0.02 * e.max():.4f}" + r" \\")
    (OUT / "tab_approx.tex").write_text(
        r"\begin{tabular}{rl|rr|r}\toprule" "\n"
        r"$m$ & approximation & median $|\Delta K^X|$ & max $|\Delta K^X|$ & max $|\Delta\rho|$\\\midrule" "\n"
        + "\n".join(rows) + "\n" r"\bottomrule\end{tabular}" "\n")


def signature_table(sg):
    rows = []
    lab = {"ss - os": "ss $-$ os", "maternal - paternal": "mat $-$ pat"}
    for _, r in sg.iterrows():
        lo, hi = r["pred_VA0.4_s0.02"], r["pred_VA0.7_s0.02"]
        up = r["pred_VA0.7_s0.05"]
        if abs(r.observed) < 2 * r.obs_se:
            share = "n.s."
        elif np.sign(hi) != np.sign(r.observed):
            share = "opp. sign"
        else:
            pct = 100 * hi / r.observed
            share = f"{pct:.0f}\\%" if pct >= 1 else "$<$1\\%"
        sg4 = lambda v: f"{v:+.4f}".replace("-0.0000", "0.0000").replace("+0.0000", "0.0000")
        rows.append(f"{lab[r.contrast]} & {r.rel} & {r.observed:+.4f} & {r.obs_se:.4f} & "
                    f"{sg4(lo)} -- {sg4(hi)} & {sg4(up)} & {share}" + r" \\")
    (OUT / "tab_signatures.tex").write_text(
        r"\begin{tabular}{ll|rr|rr|r}\toprule" "\n"
        r" & & \multicolumn{2}{c|}{observed} & \multicolumn{2}{c|}{X only, $m=0.43$} & \\" "\n"
        r"contrast & type & estimate & SE & $s_X{=}0.02$, $V_A$ 0.4--0.7 & $s_X{=}0.05$, $V_A{=}0.7$ & explained\\\midrule" "\n"
        + "\n".join(rows) + "\n" r"\bottomrule\end{tabular}" "\n")


def main():
    xk = pd.read_csv(HERE / "x_kinship.csv")
    grid = pd.read_csv(HERE / "x_kernel_grid.csv")
    sg = pd.read_csv(HERE / "signatures.csv")
    kernel_table(xk)
    am_table(grid)
    approx_table(xk)
    signature_table(sg)
    print(sorted(p.name for p in OUT.glob("*.tex")))


if __name__ == "__main__":
    main()
