"""X-linked additive kernel for the Eftedal et al. (2025) relative categories.

Random mating (exact): X kinship phi^X by the standard recursion (a male's X comes from his mother;
a female's from both parents), on the pedigree configurations of pedigrees.py.

Kernels, normalised by the sex-averaged X additive variance so that they are on the scale of the
autosomal kernel (K^X_ii averages to 1 over the sexes):
    full dosage compensation (FDC; default):  male X variance = 2 x female  ->  K^X = (4/3) phi^X
    no dosage compensation (NDC):             male X variance = 1/2 x female -> K^X = (4/1.5) phi^X (1/2)^{#males}

Assortative mating: a pathMgr model with two additive components per person, autosomal (Aa) and
X-linked (X), with sex-specific X transmission (FDC: son <- mother 1, son <- father 0, daughter <-
each parent 1/2), mates matched on Y = Aa + X + E. With s = the sex-averaged base X share of V_A,
K^X(m) is defined by the decomposition the Phase-2 model uses,
    K^total(m) = (1 - s) K^auto(m) + s K^X(m),
with K^total from the model with X and K^auto from the same model without X, both at the same
genotypic mate correlation m.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
import sympy as sp
from scipy.optimize import brentq

from pathmgr.core.model import Model
from pathmgr.core.ram import RAMEngine
from pathmgr.core.units import Units

from pedigrees import F, M, Ped


# ---- random mating: kinship recursions --------------------------------------------------------
def _depths(ped: Ped):
    d = {}

    def depth(k):
        if k not in d:
            _, mo, fa = ped.ind[k]
            d[k] = 0 if mo is None else 1 + max(depth(mo), depth(fa))
        return d[k]
    for k in ped.ind:
        depth(k)
    return d


def kinship(ped: Ped, a: str, b: str, x: bool) -> float:
    d = _depths(ped)

    @lru_cache(maxsize=None)
    def phi(a, b):
        if a == b:
            if not x:
                return 0.5
            return 1.0 if ped.sex(a) == M else 0.5
        if d[a] < d[b]:
            a, b = b, a
        sa, mo, fa = ped.ind[a]
        if mo is None:
            return 0.0                       # a is a founder and b is not its descendant
        if x and sa == M:
            return phi(mo, b)
        return 0.5 * (phi(mo, b) + phi(fa, b))
    return phi(a, b)


def line_share(ped: Ped, i: str, j: str, line: str) -> float:
    """1 if i and j carry the same mtDNA (line='mt': same maternal-line founder) or the same Y
    (line='y': both male, same paternal-line founder), else 0."""
    def founder(k):
        while True:
            _, mo, fa = ped.ind[k]
            nxt = mo if line == "mt" else fa
            if nxt is None:
                return k
            k = nxt
    if line == "y" and not (ped.sex(i) == M and ped.sex(j) == M):
        return 0.0
    return float(founder(i) == founder(j))


def kx_rm(ped: Ped, i: str, j: str, dc: str = "FDC") -> float:
    phx = kinship(ped, i, j, x=True)
    if dc == "FDC":
        return 4 / 3 * phx
    males = (ped.sex(i) == M) + (ped.sex(j) == M)
    return 4 / 1.5 * phx * 0.5 ** males


# ---- assortative mating: sex-aware equilibrium and path model ---------------------------------
def transmission(dc: str):
    """(T[child_sex][parent_role], residual factors) for the (Aa, X) vector."""
    if dc == "FDC":
        T = {F: {"mo": np.diag([.5, .5]), "fa": np.diag([.5, .5])},
             M: {"mo": np.diag([.5, 1.]), "fa": np.diag([.5, 0.])}}
        resid_x = {F: 0.25, M: 1.0}          # x V_Xf0
        vxm = 2.0                            # male base X variance / female
    else:
        T = {F: {"mo": np.diag([.5, .5]), "fa": np.diag([.5, 1.])},
             M: {"mo": np.diag([.5, .5]), "fa": np.diag([.5, 0.])}}
        resid_x = {F: 0.25, M: 0.25}
        vxm = 0.5
    return T, resid_x, vxm


class XEq:
    """Equilibrium of A = Aa + X under primary phenotypic AM (Gaussian assortment on Y)."""

    def __init__(self, VA0, s, VE, rho, dc="FDC", iters=4000):
        self.T, rx, vxm = transmission(dc)
        Vx_bar = s * VA0
        self.VXf0 = 2 * Vx_bar / (1 + vxm)
        self.VXm0 = vxm * self.VXf0
        self.Va0 = (1 - s) * VA0
        self.VE, self.rho, self.dc, self.s = VE, rho, dc, s
        self.D = {sx: np.diag([self.Va0 / 2, rx[sx] * self.VXf0]) for sx in (F, M)}
        S = {F: np.diag([self.Va0, self.VXf0]), M: np.diag([self.Va0, self.VXm0])}
        u = np.ones(2)
        for _ in range(iters):
            C = self.mate_cov(S)
            new = {}
            for sx in (F, M):
                Tm, Tf = self.T[sx]["mo"], self.T[sx]["fa"]
                new[sx] = Tm @ S[F] @ Tm.T + Tf @ S[M] @ Tf.T + Tm @ C @ Tf.T + Tf @ C.T @ Tm.T + self.D[sx]
            if max(np.abs(new[k] - S[k]).max() for k in S) < 1e-15:
                S = new
                break
            S = new
        self.S = S
        C = self.mate_cov(S)
        self.m = float(u @ C @ u / np.sqrt((u @ S[F] @ u) * (u @ S[M] @ u)))
        self.VA_bar = float((u @ S[F] @ u + u @ S[M] @ u) / 2)

    def mate_cov(self, S):
        u = np.ones(2)
        vyf, vym = u @ S[F] @ u + self.VE, u @ S[M] @ u + self.VE
        return self.rho * np.outer(S[F] @ u, S[M] @ u) / np.sqrt(vyf * vym)


def eq_for_m(m, VA0=0.45, s=0.02, VE=0.45, dc="FDC") -> XEq:
    """Equilibrium whose genotypic mate correlation is m (rho found by root finding)."""
    if m == 0:
        return XEq(VA0, s, VE, 0.0, dc)
    rho = brentq(lambda r: XEq(VA0, s, VE, r, dc).m - m, 1e-6, 0.999, xtol=1e-12)
    return XEq(VA0, s, VE, rho, dc)


def build(ped: Ped, eq: XEq) -> Model:
    mdl = Model("A + X under AM", units=Units.unstandardized())
    f = lambda v: sp.Float(float(v), 30)
    for k, (sx, mo, fa) in ped.ind.items():
        for v in ("Aa", "X", "E"):
            mdl.add_var(f"{v}_{k}", latent=True)
        mdl.add_var(f"Y_{k}")
        for v in ("Aa", "X", "E"):
            mdl.add_path(f"{v}_{k}", f"Y_{k}", 1)
        mdl.add_variance(f"E_{k}", f(eq.VE))
        if mo is None:
            S = eq.S[sx]
            mdl.add_variance(f"Aa_{k}", f(S[0, 0]))
            mdl.add_variance(f"X_{k}", f(S[1, 1]))
            if abs(S[0, 1]) > 1e-16:
                mdl.add_cov(f"Aa_{k}", f"X_{k}", f(S[0, 1]))
        else:
            for role, par in (("mo", mo), ("fa", fa)):
                T = eq.T[sx][role]
                for a, va in enumerate(("Aa", "X")):
                    if T[a, a] != 0:
                        mdl.add_path(f"{va}_{par}", f"{va}_{k}", f(T[a, a]))
            mdl.add_variance(f"Aa_{k}", f(eq.D[sx][0, 0]))
            mdl.add_variance(f"X_{k}", f(eq.D[sx][1, 1]))
    for n, (mo, fa) in enumerate(ped.unions):
        mdl.add_copath(f"Y_{mo}", f"Y_{fa}", correlation=f(eq.rho), process=f"u{n}")
    return mdl


def k_total(ped: Ped, i: str, j: str, eq: XEq) -> float:
    """Cov(A_i, A_j) / sex-averaged V_A, with A = Aa + X."""
    e = RAMEngine(build(ped, eq))
    c = sum(float(sp.N(e.cov(f"{a}_{i}", f"{b}_{j}", form="raw")))
            for a in ("Aa", "X") for b in ("Aa", "X"))
    return c / eq.VA_bar
