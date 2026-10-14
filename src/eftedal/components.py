"""Covariance components: each maps parameters to a predicted correlation for every row.

The predicted correlation for a category of relative is a sum of components, each a
variance parameter times a "sharing" coefficient that depends on the relationship:

    rho_k = V_A * K_k(m, r) + V_GxG * K2_k(m, r) + sum_j V_j * I_jk + V_family * r_f^f_k

where K is the additive-genetic kernel under assortative mating at equilibrium
(Fisher / Nagylaki; Eftedal et al. Table S9) and K2 its GxG analogue.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------------------
# Genetic kernels
# --------------------------------------------------------------------------------------


def _transmission(m, n):
    """((1+m)/2)^n, the AM-adjusted genetic correlation after n degrees of relatedness."""
    return ((1.0 + m) / 2.0) ** n


def half_multiplier(m, r):
    """Nagylaki's multiplier for half sibs and their descendants, (1+2m+rm)/(1+m)^2."""
    return (1.0 + 2.0 * m + r * m) / (1.0 + m) ** 2


class KernelData:
    """Row masks and degrees needed by `genetic_kernels`, extracted once from the data."""

    def __init__(self, data: pd.DataFrame):
        s = data["genetic_structure"].to_numpy()
        n = data["n"].to_numpy(dtype=float)
        self.size = len(data)
        self.masks = {k: s == k for k in ("full", "full_uz", "half", "step", "cocousin", "none")}
        self.n = {k: n[v] for k, v in self.masks.items()}
        self.ado_cocousin = data["rel"].to_numpy()[self.masks["none"]] == "adoptive cocousins"


def genetic_kernels(data, m: float, r: float, *, stepsib: str = "eftedal",
                    uz: str = "mixture", mz_frac: float = 0.5,
                    adoptive: str = "zero", gxg: str = "am2") -> tuple[np.ndarray, np.ndarray]:
    """Return (K, K2): additive-genetic and GxG sharing coefficients per row.

    gxg: 'am2' : K2 = K^2, the square of the assortative-mating kernel (incl. the half-sib
                 multiplier and the in-law kernels)
         'a2'  : K2 = a^2, the square of the random-mating kernel a = (1/2)^n: half relatives
                 as full ones of the same degree, in-laws and adoptive pairs 0; independent of m, r

    `data` is the annotated table or a `KernelData` built from it.

    K is the expected genotypic correlation of the pair; the phenotypic contribution is
    V_A * K. The GxG kernel is taken as the square of the additive one (two loci must both
    be shared), so K2 = K^2, except for classes linked through twins of unknown zygosity,
    where both kernels are averaged over the MZ/DZ alternatives.

    stepsib:
        'eftedal'   : m * ((1+m)/2)^2     (Eftedal et al. Table S9)
        'corrected' : m * ((1+r)/2)^2     (relative_covariance v1: every mating step after
                                          the first along a mating chain costs r, not m)
    uz:
        'mixture'   : mz_frac * K(n - 1/2) + (1 - mz_frac) * K(n + 1/2)  (Eftedal's predictions)
        'exponent'  : K(n) with the averaged n plugged into the exponent (legacy 2025 code)
    adoptive:
        'zero'      : adoptive pairs share nothing genetically
        'legacy_n99': 2025 code, which coded adoptive pairs as n = 99 without zeroing them,
                      i.e. K = ((1+m)/2)^99 (times m for adoptive co-cousins). Negligible
                      unless m is close to 1.
    """
    kd = data if isinstance(data, KernelData) else KernelData(data)
    M, N = kd.masks, kd.n
    K = np.zeros(kd.size)
    K2 = np.zeros(kd.size)

    K[M["full"]] = _transmission(m, N["full"])
    K2[M["full"]] = K[M["full"]] ** 2

    if uz == "mixture":
        lo, hi = _transmission(m, N["full_uz"] - 0.5), _transmission(m, N["full_uz"] + 0.5)
        K[M["full_uz"]] = mz_frac * lo + (1 - mz_frac) * hi
        K2[M["full_uz"]] = mz_frac * lo ** 2 + (1 - mz_frac) * hi ** 2
    elif uz == "exponent":
        K[M["full_uz"]] = _transmission(m, N["full_uz"])
        K2[M["full_uz"]] = K[M["full_uz"]] ** 2
    else:
        raise ValueError(f"unknown uz mode {uz!r}")

    K[M["half"]] = _transmission(m, N["half"]) * half_multiplier(m, r)
    K2[M["half"]] = K[M["half"]] ** 2

    if stepsib == "eftedal":
        K[M["step"]] = m * _transmission(m, N["step"])
    elif stepsib == "corrected":
        K[M["step"]] = m * ((1.0 + r) / 2.0) ** 2
    else:
        raise ValueError(f"unknown stepsib mode {stepsib!r}")
    K2[M["step"]] = K[M["step"]] ** 2

    K[M["cocousin"]] = m * _transmission(m, N["cocousin"])
    K2[M["cocousin"]] = K[M["cocousin"]] ** 2

    if adoptive == "legacy_n99":
        K[M["none"]] = _transmission(m, 99.0) * np.where(kd.ado_cocousin, m, 1.0)
        K2[M["none"]] = K[M["none"]] ** 2
    elif adoptive != "zero":
        raise ValueError(f"unknown adoptive mode {adoptive!r}")

    if gxg == "a2":
        K2 = np.zeros(kd.size)
        K2[M["full"]] = 0.25 ** N["full"]
        K2[M["half"]] = 0.25 ** N["half"]
        K2[M["full_uz"]] = mz_frac * 0.25 ** (N["full_uz"] - 0.5) + (1 - mz_frac) * 0.25 ** (N["full_uz"] + 0.5)
    elif gxg != "am2":
        raise ValueError(f"unknown gxg kernel {gxg!r}")

    return K, K2


# --------------------------------------------------------------------------------------
# Shared-environment indicators
# --------------------------------------------------------------------------------------


def household_share(data: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """(full_time, part_time) household indicators.

    'by_side' (half sibs) counts as full-time if maternal and part-time if paternal:
    children of separated parents mostly live with their mother.
    """
    h = data["household"].to_numpy()
    maternal = data["side"].to_numpy() == "maternal"
    full_time = (h == "full") | ((h == "by_side") & maternal)
    part_time = (h == "partial") | ((h == "by_side") & ~maternal)
    return full_time.astype(float), part_time.astype(float)


def env_indicator(data: pd.DataFrame, name: str) -> np.ndarray:
    """Indicator (or weight) for the simple shared-environment components."""
    full_time, part_time = household_share(data)
    cc = data["cousin_class"].to_numpy()
    if name == "house1":      # any sibling-type pair (incl. step, paternal half, adoptive)
        return full_time + part_time
    if name == "house2":      # only pairs likely to share a household full-time
        return full_time
    if name == "house3":      # blend: part-time pairs count one half
        return full_time + 0.5 * part_time
    if name == "housepart":   # part-time pairs only (use together with house2)
        return part_time
    if name == "cousins1":
        return (cc == 1).astype(float)
    if name == "cousins2":
        return (cc == 2).astype(float)
    raise ValueError(f"unknown indicator component {name!r}")


def maternal_indicator(data: pd.DataFrame) -> np.ndarray:
    """1 for categories related only through mothers (side == 'maternal'), else 0.

    Paternal, opposite-sex-linked ('other') and both-parent pairs (twins, full sibs) are 0:
    for both-parent pairs a maternal excess cannot be separated from household or genetic
    sharing, so the maternal parameter is a contrast identified within half sibs, cousins,
    half cousins, 2nd cousins and co-cousins.
    """
    return (data["side"].to_numpy() == "maternal").astype(float)


INDICATOR_COMPONENTS = ("house1", "house2", "house3", "housepart", "cousins1", "cousins2")
GENETIC_COMPONENTS = ("A", "GxG")
# 'house': V_house * (full_time + pi * part_time), with the part-time share pi estimated
# 'maternal': multiplicative, rho_k * (1 + delta * maternal_k); not a variance component
ALL_COMPONENTS = GENETIC_COMPONENTS + INDICATOR_COMPONENTS + ("house", "family", "maternal")
MODIFIER_COMPONENTS = ("maternal",)
