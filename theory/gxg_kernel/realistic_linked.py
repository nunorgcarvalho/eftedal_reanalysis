"""Linked path-model kernels at the illustrative parameters of realistic.py (V_A0=0.40, V_G0=0.20,
V_E=0.40), with every interacting pair at one recombination fraction c, for several mate
correlations. Writes results/linked/realistic.csv."""
from pathlib import Path
import pandas as pd
import ag_linked_model as alm

out = Path(__file__).parent / "results" / "linked"
out.mkdir(parents=True, exist_ok=True)
rows = []
for c in (0.0, 0.01, 0.1, 0.5):
    for rho in (0.0, 0.30, 0.45):
        eq = alm.solve_equilibrium(0.40, [c], [0.20], 0.40, rho)
        rows.append(alm.predict(eq).assign(c=c, rho_y=rho, V_A=eq.V_A, V_G=eq.V_G, V_Y=eq.V_Y))
d = pd.concat(rows, ignore_index=True)
d.to_csv(out / "realistic.csv", index=False)
with pd.option_context("display.width", 200):
    print(d.pivot_table(index="cls", columns=["c", "rho_y"], values="K_G", sort=False).round(4).to_string())
    print(d.pivot_table(index="cls", columns=["c", "rho_y"], values="cross", sort=False).round(4).to_string())
    print(d.groupby(["c", "rho_y"])[["V_A", "V_G"]].first().round(4))
