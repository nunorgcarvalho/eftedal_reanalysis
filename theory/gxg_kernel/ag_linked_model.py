"""Path model for A + pairwise GxG under primary phenotypic AM with *linked* interacting loci.

Extends ag_path_model.py. With alleles centred and LE in the base, a pair (k, l) contributes
gamma_kl x_k x_l = gamma_kl (cis + trans), where cis = h1k h1l + h2k h2l (products within a
haplotype) and trans = h1k h2l + h2k h1l (products across the two haplotypes); each carries half
of the pair's variance. Over one meiosis a parent passes on one of its cis products w.p. 1 - c and
one of its trans products w.p. c, so, per recombination bin b (pairs with recombination fraction
c_b, base variance V_b0):

    Gc_b,o = (1 - c_b)/2 (Gc_b,m + Gc_b,f) + c_b/2 (Gt_b,m + Gt_b,f) + s_cb
    Gt_b,o = B_b,mf / 4 + s_tb                 (products across the two parental gametes)
    Var[B_b] = 2 V_b0,   Var[s_cb] = V_b0/2 - V_b0 ((1-c_b)^2 + c_b^2)/4,   Var[s_tb] = 3 V_b0/8.

c_b = 1/2 recovers ag_path_model.py (G transmitted with 1/4 plus the couple term); c_b = 0 makes
the cis part a haplotype-additive value (transmitted with 1/2, like A) and the trans part a
dominance-like value (shared only through the couple). At random mating the model reproduces the
exact digenic coefficients of digenic.py. Mates are matched on Y = A + sum_b (Gc_b + Gt_b) + E.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
import sympy as sp

from pathmgr.core.model import Model
from pathmgr.core.ram import RAMEngine
from pathmgr.core.units import Units

from ag_path_model import class_pedigrees, IBD


@dataclass
class LinkedEq:
    Sigma: np.ndarray          # equilibrium covariance of (A, Gc_1..B, Gt_1..B)
    c: np.ndarray              # bin recombination fractions
    V0: np.ndarray             # bin base variances
    V_A0: float
    V_E: float
    rho: float

    @property
    def B(self):
        return len(self.c)

    @property
    def names(self):
        return ["A"] + [f"Gc{b}" for b in range(self.B)] + [f"Gt{b}" for b in range(self.B)]

    @property
    def V_Y(self):
        return float(self.Sigma.sum() + self.V_E)

    @property
    def V_A(self):
        return float(self.Sigma[0, 0])

    @property
    def V_G(self):
        return float(self.Sigma[1:, 1:].sum())


def _T(c):
    B = len(c)
    T = np.zeros((1 + 2 * B, 1 + 2 * B))
    T[0, 0] = 0.5
    for b in range(B):
        T[1 + b, 1 + b] = (1 - c[b]) / 2
        T[1 + b, 1 + B + b] = c[b] / 2
    return T


def _resid(V_A0, c, V0):
    s = (1 - c) ** 2 + c ** 2
    return np.concatenate([[V_A0 / 2], V0 / 2 - V0 * s / 4, 3 * V0 / 8])


def solve_equilibrium(V_A0, c, V0, V_E, rho, iters=5000, tol=1e-14) -> LinkedEq:
    c, V0 = np.asarray(c, float), np.asarray(V0, float)
    B = len(c)
    S = np.diag(np.concatenate([[V_A0], V0 / 2, V0 / 2]))      # base (LE) covariance
    T = _T(c)
    D = np.diag(_resid(V_A0, c, V0))
    D[1 + B:, 1 + B:] += np.diag(2 * V0 / 16)                  # couple term B_b / 4
    for _ in range(iters):
        u = S.sum(1)
        C = rho * np.outer(u, u) / (S.sum() + V_E)             # Cov(z_m, z_f)
        S_new = T @ (2 * S + C + C.T) @ T.T + D
        if np.max(np.abs(S_new - S)) < tol:
            S = S_new
            break
        S = S_new
    return LinkedEq(S, c, V0, V_A0, V_E, rho)


def build(ped, eq: LinkedEq) -> Model:
    m = Model("A + linked GxG under AM", units=Units.unstandardized())
    f = lambda v: sp.Float(v, 30)
    B, names = eq.B, eq.names
    resid = _resid(eq.V_A0, eq.c, eq.V0)
    T = _T(eq.c)
    children = {}
    for k, ind in ped.individuals.items():
        if not ind.is_founder:
            children.setdefault((ind.maternal, ind.paternal), []).append(k)
    for (mo, fa) in children:
        for b in range(B):
            m.add_var(f"B{b}_{mo}_{fa}", latent=True)
            m.add_variance(f"B{b}_{mo}_{fa}", f(2 * eq.V0[b]))
    for k, ind in ped.individuals.items():
        for v in names + ["E"]:
            m.add_var(f"{v}_{k}", latent=True)
        m.add_var(f"Y_{k}")
        for v in names + ["E"]:
            m.add_path(f"{v}_{k}", f"Y_{k}", 1)
        m.add_variance(f"E_{k}", f(eq.V_E))
        if ind.is_founder:
            for a in range(len(names)):
                m.add_variance(f"{names[a]}_{k}", f(eq.Sigma[a, a]))
                for b in range(a + 1, len(names)):
                    if abs(eq.Sigma[a, b]) > 1e-15:
                        m.add_cov(f"{names[a]}_{k}", f"{names[b]}_{k}", f(eq.Sigma[a, b]))
        else:
            mo, fa = ind.maternal, ind.paternal
            for a in range(len(names)):
                for b in range(len(names)):
                    if T[a, b] != 0:
                        for par in (mo, fa):
                            m.add_path(f"{names[b]}_{par}", f"{names[a]}_{k}", f(T[a, b]))
                m.add_variance(f"{names[a]}_{k}", f(resid[a]))
            for b in range(B):
                m.add_path(f"B{b}_{mo}_{fa}", f"Gt{b}_{k}", sp.Rational(1, 4))
    for cpl in ped.couples:
        m.add_copath(f"Y_{cpl.maternal}", f"Y_{cpl.paternal}", correlation=f(eq.rho), process=cpl.key)
    return m


def predict(eq: LinkedEq, classes=None) -> pd.DataFrame:
    rows = []
    B = eq.B
    for name, (ped, i, j) in class_pedigrees().items():
        if classes and name not in classes:
            continue
        e = RAMEngine(build(ped, eq))
        G = lambda k: [f"Gc{b}_{k}" for b in range(B)] + [f"Gt{b}_{k}" for b in range(B)]
        cov = lambda xs, ys: sum(float(sp.N(e.cov(x, y, form="raw"))) for x in xs for y in ys)
        vA_i, vA_j = cov([f"A_{i}"], [f"A_{i}"]), cov([f"A_{j}"], [f"A_{j}"])
        vG_i, vG_j = cov(G(i), G(i)), cov(G(j), G(j))
        sA, sG = np.sqrt(vA_i * vA_j), np.sqrt(vG_i * vG_j)
        K_A = cov([f"A_{i}"], [f"A_{j}"]) / sA
        rows.append(dict(cls=name, K_A=K_A, K_G=cov(G(i), G(j)) / sG,
                         cross=(cov([f"A_{i}"], G(j)) + cov(G(i), [f"A_{j}"])) / 2 / np.sqrt(sA * sG),
                         K_A_sq=K_A ** 2, a_sq=IBD[name] ** 2,
                         stationary_check=max(abs(vA_i / eq.V_A - 1), abs(vG_j / eq.V_G - 1))))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    import time
    from digenic import psi_table, interp
    import ag_path_model as apm
    t0 = time.time()
    c = np.array([0.0005, 0.05, 0.5])
    V0 = np.array([0.1, 0.1, 0.1])
    # random mating: must equal the V0-weighted digenic coefficients
    eq = solve_equilibrium(0.45, c, V0, 0.25, 0.0)
    p = predict(eq)
    tab = psi_table(cs=sorted(set([0.0, 0.0005, 0.05, 0.5])), R=200_000)
    p["digenic"] = [np.sum(V0 * np.array([interp(tab, k, ci) for ci in c])) / V0.sum() for k in p.cls]
    print(p.round(4).to_string(index=False))
    # c = 1/2 everywhere under AM: must equal the unlinked model
    eq1 = solve_equilibrium(0.45, [0.5], [0.3], 0.25, 0.85)
    old = apm.predict(apm.solve_equilibrium(0.45, 0.3, 0.25, 0.85))
    new = predict(eq1)
    print(pd.DataFrame(dict(cls=new.cls, K_G_new=new.K_G, K_G_old=old.K_G, K_A_new=new.K_A, K_A_old=old.K_A)).round(5).to_string(index=False))
    print(f"{time.time() - t0:.1f}s")
