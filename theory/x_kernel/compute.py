"""Compute the X-linked kernel for every (rel, side, sex) category.

Writes
    x_kinship.csv   one row per category: X kinship and kernels at random mating (FDC, NDC),
                    mtDNA / Y line-sharing probabilities,
                    exact K^X under AM from the sex-aware path model at m = 0.2 and 0.43 (FDC,
                    s = 0.02), the autosomal kernel from the same model, and two approximations
    x_configs.csv   every pedigree configuration with its weight (what the category averages over)

    PYTHONPATH=<worktree>/src python compute.py
"""

from pathlib import Path

import numpy as np
import pandas as pd

from eftedal.data import load
from pedigrees import configs
from x_kernel import eq_for_m, k_total, kinship, kx_rm, line_share

HERE = Path(__file__).parent
MS = (0.2, 0.43)
S = 0.02


def main():
    cats = load("MRE8")[["rel", "side", "sex"]].drop_duplicates().reset_index(drop=True)
    eqs = {m: (eq_for_m(m, s=0.0), eq_for_m(m, s=S)) for m in MS}
    rows, crow = [], []
    for _, c in cats.iterrows():
        cf = configs(c.rel, c.side, c.sex)
        r = dict(rel=c.rel, side=c.side, sex=c.sex, n_configs=len(cf))
        if not cf:                                            # adoptive: unrelated
            r.update(phiX=0.0, a=0.0, KX_FDC=0.0, KX_NDC=0.0, mt_share=0.0, y_share=0.0)
            for m in MS:
                r.update({f"Kauto_m{m}": 0.0, f"KX_m{m}": 0.0})
            rows.append(r)
            continue
        acc = dict(phiX=0.0, a=0.0, KX_FDC=0.0, KX_NDC=0.0, mt_share=0.0, y_share=0.0)
        accm = {m: [0.0, 0.0] for m in MS}
        for cfg in cf:
            phx = kinship(cfg.ped, cfg.i, cfg.j, x=True)
            a = 2 * kinship(cfg.ped, cfg.i, cfg.j, x=False)
            kf, kn = kx_rm(cfg.ped, cfg.i, cfg.j, "FDC"), kx_rm(cfg.ped, cfg.i, cfg.j, "NDC")
            mt, y = line_share(cfg.ped, cfg.i, cfg.j, "mt"), line_share(cfg.ped, cfg.i, cfg.j, "y")
            for k, v in (("phiX", phx), ("a", a), ("KX_FDC", kf), ("KX_NDC", kn), ("mt_share", mt), ("y_share", y)):
                acc[k] += cfg.w * v
            cr = dict(rel=c.rel, side=c.side, sex=c.sex, config=cfg.label, weight=cfg.w,
                      phiX=phx, a=a, KX_FDC=kf, KX_NDC=kn, mt_share=mt, y_share=y)
            for m in MS:
                e0, e1 = eqs[m]
                ka = k_total(cfg.ped, cfg.i, cfg.j, e0)
                kt = k_total(cfg.ped, cfg.i, cfg.j, e1)
                kx = (kt - (1 - S) * ka) / S
                accm[m][0] += cfg.w * ka
                accm[m][1] += cfg.w * kx
                cr.update({f"Kauto_m{m}": ka, f"KX_m{m}": kx})
            crow.append(cr)
        r.update(acc)
        for m in MS:
            r.update({f"Kauto_m{m}": accm[m][0], f"KX_m{m}": accm[m][1]})
        rows.append(r)
        print(c.rel, c.side, c.sex, len(cf), flush=True)
    out = pd.DataFrame(rows)
    for m in MS:
        a0 = out.a.where(out.a > 0)
        out[f"KX_mult_m{m}"] = (out.KX_FDC * out[f"Kauto_m{m}"] / a0).fillna(0.0)   # RM X kernel x auto AM factor
        out[f"KX_add_m{m}"] = out.KX_FDC + out[f"Kauto_m{m}"] - out.a              # RM X kernel + auto AM excess
    out.round(6).to_csv(HERE / "x_kinship.csv", index=False)
    pd.DataFrame(crow).round(6).to_csv(HERE / "x_configs.csv", index=False)
    for m in MS:
        bio = out[out.a > 0]
        e_add = (bio[f"KX_add_m{m}"] - bio[f"KX_m{m}"]).abs()
        e_mul = (bio[f"KX_mult_m{m}"] - bio[f"KX_m{m}"]).abs()
        print(f"m={m}: |additive approx error| max {e_add.max():.3f} median {e_add.median():.3f}; "
              f"|multiplicative| max {e_mul.max():.3f} median {e_mul.median():.3f}")


if __name__ == "__main__":
    main()
