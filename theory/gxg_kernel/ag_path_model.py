"""Path-model predictions for additive (A) and pairwise-epistatic (G) values of relatives under
primary phenotypic assortative mating at equilibrium, built with pathMgr.

The derivation (see gxg_kernel.tex) shows that, for interaction effects independent of the additive
effects, G is transmitted *linearly* to leading order in 1/M:

    G_o = (G_m + G_f)/4 + B_mf/4 + s_G,      Var[B_mf] = 2 V_G0,  Var[s_G] = 3 V_G0 / 4,
    A_o = (A_m + A_f)/2 + s_A,               Var[s_A] = V_A0 / 2,

where B_mf = sum gamma_kl (x_mk x_fl + x_fk x_ml) is a couple-specific cross-parent interaction value
shared by the couple's children, and V_A0, V_G0 are the base (linkage-equilibrium) variances. Mates
are matched on Y = A + G + E (one co-path per couple, declared by the phenotypic correlation rho_y),
which induces Cov[A_m, G_f] and hence a within-person Cov[A, G] = c at equilibrium.

This module (1) solves the equilibrium (V_A, V_G, c), (2) builds the stationary model over a pedigree
for each relationship class, and (3) returns the kernels the simulation measures:
    K_A = Cov[A_i, A_j]/V_A,  K_G = Cov[G_i, G_j]/V_G,
    cross = (Cov[A_i, G_j] + Cov[G_i, A_j]) / (2 sqrt(V_A V_G)).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
import sympy as sp

from pathmgr.core.model import Model
from pathmgr.core.ram import RAMEngine
from pathmgr.core.units import Units
from pathmgr.genetics.pedigree import Pedigree


@dataclass
class Equilibrium:
    V_A: float
    V_G: float
    c: float
    V_E: float
    rho: float
    V_A0: float
    V_G0: float

    @property
    def V_Y(self):
        return self.V_A + self.V_G + 2 * self.c + self.V_E


def solve_equilibrium(V_A0, V_G0, V_E, rho, iters=2000) -> Equilibrium:
    """Fixed point of the generation recursion (derivation, eq. recursions)."""
    V_A, V_G, c = V_A0, V_G0, 0.0
    for _ in range(iters):
        V_Y = V_A + V_G + 2 * c + V_E
        aa = rho * (V_A + c) ** 2 / V_Y
        gg = rho * (V_G + c) ** 2 / V_Y
        ag = rho * (V_A + c) * (V_G + c) / V_Y
        V_A, V_G, c = ((V_A + aa) / 2 + V_A0 / 2,
                       (V_G + gg) / 8 + V_G0 / 8 + 3 * V_G0 / 4,
                       (c + ag) / 4)
    return Equilibrium(V_A, V_G, c, V_E, rho, V_A0, V_G0)


def build(ped: Pedigree, eq: Equilibrium) -> Model:
    m = Model("A + GxG under AM", units=Units.unstandardized())
    f = lambda v: sp.Float(v, 30)
    children = {}
    for k, ind in ped.individuals.items():
        if not ind.is_founder:
            children.setdefault((ind.maternal, ind.paternal), []).append(k)
    for (mo, fa), kids in children.items():
        m.add_var(f"B_{mo}_{fa}", latent=True)
        m.add_variance(f"B_{mo}_{fa}", f(2 * eq.V_G0))
    for k, ind in ped.individuals.items():
        for v in ("A", "G", "E"):
            m.add_var(f"{v}_{k}", latent=True)
        m.add_var(f"Y_{k}")
        for v in ("A", "G", "E"):
            m.add_path(f"{v}_{k}", f"Y_{k}", 1)
        m.add_variance(f"E_{k}", f(eq.V_E))
        if ind.is_founder:
            m.add_variance(f"A_{k}", f(eq.V_A))
            m.add_variance(f"G_{k}", f(eq.V_G))
            m.add_cov(f"A_{k}", f"G_{k}", f(eq.c))
        else:
            mo, fa = ind.maternal, ind.paternal
            m.add_path(f"A_{mo}", f"A_{k}", sp.Rational(1, 2))
            m.add_path(f"A_{fa}", f"A_{k}", sp.Rational(1, 2))
            m.add_variance(f"A_{k}", f(eq.V_A0 / 2))                 # segregation s_A
            m.add_path(f"G_{mo}", f"G_{k}", sp.Rational(1, 4))
            m.add_path(f"G_{fa}", f"G_{k}", sp.Rational(1, 4))
            m.add_path(f"B_{mo}_{fa}", f"G_{k}", sp.Rational(1, 4))
            m.add_variance(f"G_{k}", f(3 * eq.V_G0 / 4))             # segregation s_G
    for cpl in ped.couples:
        m.add_copath(f"Y_{cpl.maternal}", f"Y_{cpl.paternal}", correlation=f(eq.rho),
                     process=cpl.key)
    return m


# ---- relationship-class pedigrees (founders at generation 0) ---------------------------------
def _ped(spec):
    """spec: list of (key, generation, mother, father); couples inferred from children, plus
    explicit ('mate', mother, father) entries for childless unions."""
    p = Pedigree()
    for s in spec:
        if s[0] != "mate":
            p.add(s[0], s[1], s[2], s[3])
    seen = set()
    for s in spec:
        if s[0] == "mate":
            pair = (s[1], s[2])
        elif s[2] is not None:
            pair = (s[2], s[3])
        else:
            continue
        if pair not in seen:
            seen.add(pair)
            p.mate(*pair)
    return p


def class_pedigrees():
    """{class: (pedigree, proband_i, proband_j)}."""
    P = {}
    P["parent-offspring"] = (_ped([("m", 0, None, None), ("f", 0, None, None), ("c", 1, "m", "f")]), "m", "c")
    P["full sibs"] = (_ped([("m", 0, None, None), ("f", 0, None, None),
                            ("a", 1, "m", "f"), ("b", 1, "m", "f")]), "a", "b")
    P["half sibs"] = (_ped([("p", 0, None, None), ("x", 0, None, None), ("y", 0, None, None),
                            ("a", 1, "p", "x"), ("b", 1, "p", "y")]), "a", "b")
    P["step sibs"] = (_ped([("fA", 0, None, None), ("mA", 0, None, None), ("fB", 0, None, None),
                            ("mB", 0, None, None), ("a", 1, "fA", "mA"), ("b", 1, "fB", "mB"),
                            ("mate", "fB", "mA")]), "a", "b")
    P["grandparent"] = (_ped([("g1", 0, None, None), ("g2", 0, None, None), ("p", 1, "g1", "g2"),
                              ("o", 1, None, None), ("c", 2, "p", "o")]), "g1", "c")
    P["avuncular"] = (_ped([("g1", 0, None, None), ("g2", 0, None, None), ("p", 1, "g1", "g2"),
                            ("u", 1, "g1", "g2"), ("o", 1, None, None), ("c", 2, "p", "o")]), "u", "c")
    P["first cousins"] = (_ped([("g1", 0, None, None), ("g2", 0, None, None),
                                ("p1", 1, "g1", "g2"), ("p2", 1, "g1", "g2"),
                                ("o1", 1, None, None), ("o2", 1, None, None),
                                ("c1", 2, "p1", "o1"), ("c2", 2, "p2", "o2")]), "c1", "c2")
    P["co-cousins"] = (_ped([("g1", 0, None, None), ("g2", 0, None, None), ("g3", 0, None, None),
                             ("g4", 0, None, None), ("s1", 1, "g1", "g2"), ("p1", 1, "g1", "g2"),
                             ("t1", 1, "g3", "g4"), ("p2", 1, "g3", "g4"),
                             ("o1", 1, None, None), ("o2", 1, None, None),
                             ("mate", "s1", "t1"),
                             ("c1", 2, "p1", "o1"), ("c2", 2, "p2", "o2")]), "c1", "c2")
    P["half cousins"] = (_ped([("g", 0, None, None), ("x", 0, None, None), ("y", 0, None, None),
                               ("h1", 1, "g", "x"), ("h2", 1, "g", "y"),
                               ("o1", 1, None, None), ("o2", 1, None, None),
                               ("c1", 2, "h1", "o1"), ("c2", 2, "h2", "o2")]), "c1", "c2")
    P["second cousins"] = (_ped([("g1", 0, None, None), ("g2", 0, None, None),
                                 ("p1", 1, "g1", "g2"), ("p2", 1, "g1", "g2"),
                                 ("o1", 1, None, None), ("o2", 1, None, None),
                                 ("q1", 2, "p1", "o1"), ("q2", 2, "p2", "o2"),
                                 ("r1", 2, None, None), ("r2", 2, None, None),
                                 ("c1", 3, "q1", "r1"), ("c2", 3, "q2", "r2")]), "c1", "c2")
    return P


# pedigree kinship coefficient a = 2*phi (IBD-sharing; no AM)
IBD = {"parent-offspring": 1 / 2, "full sibs": 1 / 2, "half sibs": 1 / 4, "step sibs": 0.0,
       "grandparent": 1 / 4, "avuncular": 1 / 4, "first cousins": 1 / 8, "co-cousins": 0.0,
       "half cousins": 1 / 16, "second cousins": 1 / 32}


def predict(eq: Equilibrium, classes=None) -> pd.DataFrame:
    rows = []
    for name, (ped, i, j) in class_pedigrees().items():
        if classes and name not in classes:
            continue
        e = RAMEngine(build(ped, eq))
        num = lambda x, y: float(sp.N(e.cov(x, y, form="raw")))
        vA_i, vA_j = num(f"A_{i}", f"A_{i}"), num(f"A_{j}", f"A_{j}")
        vG_i, vG_j = num(f"G_{i}", f"G_{i}"), num(f"G_{j}", f"G_{j}")
        sA, sG = np.sqrt(vA_i * vA_j), np.sqrt(vG_i * vG_j)
        K_A = num(f"A_{i}", f"A_{j}") / sA
        rows.append(dict(cls=name, K_A=K_A, K_G=num(f"G_{i}", f"G_{j}") / sG,
                         cross=(num(f"A_{i}", f"G_{j}") + num(f"G_{i}", f"A_{j}")) / 2 / np.sqrt(sA * sG),
                         K_A_sq=K_A ** 2, a_sq=IBD[name] ** 2,
                         stationary_check=max(abs(vA_i / eq.V_A - 1), abs(vG_j / eq.V_G - 1))))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    eq = solve_equilibrium(0.45, 0.30, 0.25, 0.85)
    print(eq)
    print(predict(eq).round(4).to_string(index=False))
