"""Path-model kernels at illustrative, test-score-like parameters (base V_A0=0.40, V_G0=0.20,
V_E=0.40) for several mate correlations, incl. random mating (rho_y=0, where K_G must equal a^2).
Writes results/realistic.csv."""
from pathlib import Path
import pandas as pd
from ag_path_model import solve_equilibrium, predict

rows = []
for rho in (0.0, 0.30, 0.45):
    eq = solve_equilibrium(0.40, 0.20, 0.40, rho)
    p = predict(eq).assign(rho_y=rho, V_A=eq.V_A, V_G=eq.V_G, c=eq.c, V_Y=eq.V_Y)
    rows.append(p)
out = pd.concat(rows, ignore_index=True)
out.to_csv(Path(__file__).parent / "results" / "realistic.csv", index=False)
with pd.option_context("display.width", 200):
    print(out.round(4).to_string(index=False))
