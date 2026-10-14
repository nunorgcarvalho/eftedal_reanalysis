"""LaTeX tables/sentences for the linkage section of gxg_kernel.tex, from results/linked/*.csv."""
from pathlib import Path

import numpy as np
import pandas as pd

R = Path(__file__).parent / "results" / "linked"
ORDER = ["parent-offspring", "full sibs", "half sibs", "grandparent", "avuncular", "first cousins",
         "half cousins", "second cousins", "step sibs", "co-cousins"]
DELTA = {"full sibs": 0.25}                          # P(both alleles IBD), non-inbred

# ---- psi(c) at random mating -----------------------------------------------------------------
psi = pd.read_csv(R / "psi.csv")
CS = [0.0, 0.01, 0.1, 0.5]
rows = []
for k in ORDER:
    d = psi[psi.cls == k].set_index("c")
    a = d.a.iloc[0]
    lim = (a + DELTA.get(k, 0.0)) / 2
    ratio = f"{lim / a ** 2:.0f}" if a > 0 else "--"
    rows.append(f"{k} & {a:.4f} & " + " & ".join(f"{d.loc[c, 'psi']:.4f}" for c in CS) + f" & {ratio}" + r" \\")
(R / "psi_table.tex").write_text(
    r"\begin{tabular}{lr|rrrr|r}\toprule" "\n"
    r"class & $a$ & $c=0$ & $0.01$ & $0.1$ & $0.5$ (${}=a^2$) & $\psi(0)/a^2$\\\midrule" "\n"
    + "\n".join(rows) + "\n" r"\bottomrule\end{tabular}" "\n")

# ---- realistic parameters: c x rho ---------------------------------------------------------
re = pd.read_csv(R / "realistic.csv")
rows = []
for k in ORDER:
    d = re[re.cls == k].set_index(["c", "rho_y"])
    rows.append(f"{k} & " + " & ".join(
        f"{d.loc[(c, 0.0), 'K_G']:.4f} & {d.loc[(c, 0.45), 'K_G']:.4f}" for c in (0.0, 0.1, 0.5)) + r" \\")
(R / "realistic_table.tex").write_text(
    r"\begin{tabular}{l|rr|rr|rr}\toprule" "\n"
    r" & \multicolumn{2}{c|}{$c=0$} & \multicolumn{2}{c|}{$c=0.1$} & \multicolumn{2}{c}{$c=\frac12$}\\" "\n"
    r"class & RM & $\rho_y{=}0.45$ & RM & $\rho_y{=}0.45$ & RM & $\rho_y{=}0.45$\\\midrule" "\n"
    + "\n".join(rows) + "\n" r"\bottomrule\end{tabular}" "\n")

# ---- simulation agreement sentence -----------------------------------------------------------
cp = pd.read_csv(R / "comparison.csv")
eq = pd.read_csv(R / "equilibrium.csv")
undirected = cp[~cp.arch.str.contains("dir")]
dev = (undirected.K_G - undirected.linked_am_T).abs()
z = dev / undirected.K_G_se.clip(lower=1e-4)
dir_ = cp[cp.arch.str.contains("dir") & (cp.regime == "strong")].set_index("cls")
st = cp[cp.regime == "strong"]
k2ratio = (st[~st.arch.str.contains("dir") & st.cls.isin(["full sibs", "half sibs", "first cousins"])]
           .assign(r=lambda x: x.K_A_sq / x.K_G).r)
(R / "sim_sentence.tex").write_text(
    f"Across the five undirected architectures, both regimes and all eight classes, the linkage path model "
    f"is within ${dev.max():.3f}$ of the simulated GxG correlation (median absolute error ${dev.median():.4f}$; "
    f"{(z > 3).sum()} of {len(z)} cells beyond 3 replicate standard errors, which are themselves estimated "
    f"from three replicates). $K^2$ overstates the sibling, half-sibling and first-cousin GxG correlation "
    f"under strong AM by ${k2ratio.min():.1f}$--${k2ratio.max():.1f}\\times$. For the directional "
    f"within-gene architecture (3d) under strong AM the path model, which assumes interaction effects "
    f"unrelated in sign to the additive effects, gives "
    + ", ".join(f"{k} ${dir_.loc[k, "linked_am_T"]:.3f}$ vs.\\ simulated ${dir_.loc[k, 'K_G']:.3f}$"
                for k in ("full sibs", "half sibs", "first cousins")) + ".\n")
print((R / "sim_sentence.tex").read_text())
print(eq.round(4).to_string(index=False))
