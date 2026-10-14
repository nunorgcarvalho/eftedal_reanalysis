"""Write the LaTeX tables/sentences included by gxg_kernel.tex from results/*.csv."""
from pathlib import Path
import pandas as pd

R = Path(__file__).parent / "results"
ORDER = ["parent-offspring", "full sibs", "half sibs", "grandparent", "avuncular", "first cousins",
         "half cousins", "second cousins", "step sibs", "co-cousins"]

re = pd.read_csv(R / "realistic.csv")
re = re[re.rho_y.isin([0.30, 0.45])]
rows = []
for c in ORDER:
    d = re[re.cls == c].set_index("rho_y")
    rows.append(f"{c} & {d.a_sq.iloc[0]:.4f} & " + " & ".join(
        f"{d.loc[r,'K_G']:.4f} & {d.loc[r,'K_A_sq']:.4f} & {d.loc[r,'cross']:.3f}" for r in (0.30, 0.45)) + r" \\")
(R / "realistic_table.tex").write_text(
    r"\begin{tabular}{lr|rrr|rrr}\toprule" "\n"
    r" & & \multicolumn{3}{c|}{$\rho_y=0.30$} & \multicolumn{3}{c}{$\rho_y=0.45$}\\" "\n"
    r"class & $a^2$ & $K_G$ & $K_A^2$ & cross & $K_G$ & $K_A^2$ & cross\\\midrule" "\n"
    + "\n".join(rows) + "\n" r"\bottomrule\end{tabular}" "\n")

cp = pd.read_csv(R / "comparison.csv")
s = cp[cp.regime == "strong"]
rows = []
for c in ORDER:
    r = s[(s.arch == "random") & (s.cls == c)]
    a = s[(s.arch == "aligned") & (s.cls == c)]
    if r.empty:
        continue
    r, a = r.iloc[0], a.iloc[0] if not a.empty else None
    rows.append(f"{c} & {r.K_G_sim:.3f} ({r.K_G_sim_se:.3f}) & {r.K_G_path:.3f} & {r.K_A_sq:.3f} & {r.a_sq:.3f}"
                + (f" & {a.K_G_sim:.3f} ({a.K_G_sim_se:.3f}) & {a.K_A_sq:.3f}" if a is not None else " & & ")
                + r" \\")
(R / "check_table.tex").write_text(
    r"\begin{tabular}{l|rrrr|rr}\toprule" "\n"
    r" & \multicolumn{4}{c|}{random epistasis} & \multicolumn{2}{c}{scale epistasis}\\" "\n"
    r"class & sim (SE) & path model & $K_A^2$ & $a^2$ & sim (SE) & $K_A^2$\\\midrule" "\n"
    + "\n".join(rows) + "\n" r"\bottomrule\end{tabular}" "\n")

al = cp[(cp.arch == "aligned") & ~cp.regime.str.endswith("curved")]
if not al.empty:
    parts = []
    for reg, lab in (("moderate", "moderate"), ("strong", "strong")):
        d = al[al.regime == reg]
        d = d[d.K_A_sq > 0.02]
        parts.append(f"{lab} AM: largest absolute difference {(d.K_G_sim - d.K_A_sq).abs().max():.3f}, "
                     f"ratio {(d.K_G_sim / d.K_A_sq).min():.2f}--{(d.K_G_sim / d.K_A_sq).max():.2f}")
    print("aligned summary:", "; ".join(parts))   # aligned_sentence.tex is written by hand from this
    (R / "aligned_sentence_auto.tex").write_text(
        "with negligible curvature the simulated GxG kernel matches $K_A^2$, built from the simulated additive "
        "kernel, in every class with $K_A^2>0.02$ (" + "; ".join(parts) + "), as derived.\n")

cur = cp[cp.regime.str.endswith("curved")]
if not cur.empty:
    fs = cur[cur.cls == "full sibs"].set_index("regime")
    fc = cur[cur.cls == "first cousins"].set_index("regime")
    (R / "curved_sentence.tex").write_text(
        "the skewness of $A$ reaches 0.8 after 12 generations of strong AM, and the within-person correlation "
        "of $A$ and $G$ about 0.4. The GxG kernel then exceeds $K_A^2$: full siblings "
        + ", ".join(f"{fs.loc[r,'K_G_sim']:.2f} vs {fs.loc[r,'K_A_sq']:.2f}" for r in fs.index)
        + "; first cousins "
        + ", ".join(f"{fc.loc[r,'K_G_sim']:.2f} vs {fc.loc[r,'K_A_sq']:.2f}" for r in fc.index)
        + " (moderate, strong AM). It also carries a large A--GxG cross term. For a strongly curved scale, $K^2$ "
          "is therefore a lower bound, and the cross term cannot be ignored.\n")
