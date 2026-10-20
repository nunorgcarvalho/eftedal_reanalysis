"""Exact X kernel on a grid of the genotypic mate correlation m (FDC, s = 0.02), for interpolation
by the Phase-2 model. Writes x_kernel_grid.csv with, per category and m,
    Kauto   autosomal kernel of the same path model (no X)
    KX      X kernel, defined by K^total = (1 - s) Kauto + s KX
    D       KX - Kauto: the X correction, K^total = Kauto + s D
Also checks that the kernels depend on m only, not on h^2 (two base variance settings).

    PYTHONPATH=<worktree>/src python grid.py
"""

from pathlib import Path

import numpy as np
import pandas as pd

from eftedal.data import load
from pedigrees import configs
from x_kernel import eq_for_m, k_total

HERE = Path(__file__).parent
S = 0.02
GRID = np.round(np.sort(np.append(np.arange(0.0, 0.601, 0.05), 0.43)), 2)   # 0.43 = Eftedal et al. estimate


def kernels(cats, m, VA0, VE):
    e0, e1 = eq_for_m(m, VA0=VA0, s=0.0, VE=VE), eq_for_m(m, VA0=VA0, s=S, VE=VE)
    rows = []
    for _, c in cats.iterrows():
        ka = kx = 0.0
        for cfg in configs(c.rel, c.side, c.sex):
            a = k_total(cfg.ped, cfg.i, cfg.j, e0)
            t = k_total(cfg.ped, cfg.i, cfg.j, e1)
            ka += cfg.w * a
            kx += cfg.w * (t - (1 - S) * a) / S
        rows.append(dict(rel=c.rel, side=c.side, sex=c.sex, m=m, Kauto=ka, KX=kx, D=kx - ka))
    return pd.DataFrame(rows)


def main():
    cats = load("MRE8")[["rel", "side", "sex"]].drop_duplicates().reset_index(drop=True)
    out = pd.concat([kernels(cats, m, 0.45, 0.45) for m in GRID], ignore_index=True)
    out.round(6).to_csv(HERE / "x_kernel_grid.csv", index=False)
    alt = kernels(cats, 0.43, 0.70, 0.20)                     # h^2 ~ 0.78 instead of 0.5
    ref = out[out.m == 0.43].reset_index(drop=True)
    print("max |D(h2 high) - D(h2 0.5)| at m=0.43:", float((alt.D - ref.D).abs().max()))


if __name__ == "__main__":
    main()
