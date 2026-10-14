"""Exact random-mating GxG kernel for a pair of loci at recombination fraction c.

For pairwise (A x A) epistasis on standardized genotypes in a non-inbred, linkage-equilibrium base
population, the covariance between relatives i and j contributed by the pair (k, l) is
gamma_kl^2 * psi_ij(c_kl), where

    psi_ij(c) = E[a_ij(k) a_ij(l)],   a_ij(k) = realized relatedness at locus k (twice kinship),

the "digenic" identity coefficient (Cockerham 1956). Unlinked loci (c = 1/2) give psi = a^2 for
these relationships; complete linkage (c = 0) turns the pair into one haplotype locus and gives
psi = (a + delta)/2, half additive-like and half dominance-like (delta = P(both alleles IBD)).

This module gene-drops two loci through the class pedigrees of ag_path_model.py and returns psi(c)
on a grid (Monte Carlo, but with closed forms checked in __main__).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ag_path_model import class_pedigrees


def psi_class(ped, i, j, c, R=400_000, rng=None):
    """Monte Carlo E[a_ij(k) a_ij(l)] for loci k, l at recombination fraction c."""
    rng = rng or np.random.default_rng(1)
    lab = {}                                    # key -> (2 alleles at k, 2 alleles at l), shape (R, 2)
    nxt = 0
    for key, ind in ped.individuals.items():    # insertion order is generation order
        if ind.is_founder:
            ak = np.tile([nxt, nxt + 1], (R, 1))
            lab[key] = (ak, ak.copy())          # founder allele labels, same haplotype labels at k and l
            nxt += 2
            continue
        gk, gl = [], []
        for par in (ind.maternal, ind.paternal):
            pk, pl = lab[par]
            s = rng.integers(0, 2, R)                          # haplotype transmitted at k
            t = s ^ (rng.random(R) < c)                         # at l: switch w.p. c
            gk.append(pk[np.arange(R), s])
            gl.append(pl[np.arange(R), t])
        lab[key] = (np.stack(gk, 1), np.stack(gl, 1))

    def a(x, y):                                 # realized relatedness = (# IBD allele pairs)/2
        return sum((x[:, u] == y[:, v]).astype(float) for u in range(2) for v in range(2)) / 2

    a_k = a(lab[i][0], lab[j][0])
    a_l = a(lab[i][1], lab[j][1])
    return float(np.mean(a_k * a_l)), float(np.mean(a_k))


def psi_table(cs=(0.0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5), R=400_000):
    rows = []
    for name, (ped, i, j) in class_pedigrees().items():
        for c in cs:
            p, a = psi_class(ped, i, j, c, R)
            rows.append(dict(cls=name, c=c, psi=p, a=a))
    return pd.DataFrame(rows)


def interp(table: pd.DataFrame, cls: str, c):
    t = table[table.cls == cls].sort_values("c")
    return np.interp(c, t.c, t.psi)


if __name__ == "__main__":
    from pathlib import Path
    t = psi_table(cs=(0.0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5))
    out = Path(__file__).parent / "results" / "linked"
    out.mkdir(parents=True, exist_ok=True)
    t.to_csv(out / "psi.csv", index=False)
    w = t.pivot(index="cls", columns="c", values="psi")
    print(w.round(4).to_string())
    # closed forms: s = (1-c)^2 + c^2
    s = lambda c: (1 - c) ** 2 + c ** 2
    for c in (0.0, 0.1, 0.5):
        print(f"c={c}: FS {s(c)/4+1/8:.4f}  HS {s(c)/8:.4f}  GP {(1-c)/8:.4f}  PO 0.2500")
