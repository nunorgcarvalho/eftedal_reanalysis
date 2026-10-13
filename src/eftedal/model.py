"""Model specification (`FitSpec`) and the `Model` that turns parameters into predictions
and a loss.

Every analysis choice lives in `FitSpec`; scripts only build specs and call `fit()`.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Mapping

import numpy as np
import pandas as pd

from . import components as C
from .data import cohort_mask

PARAM_ORDER = (
    "V_A", "V_GxG", "V_house", "V_house1", "V_house2", "V_house3", "V_housepart",
    "V_cousins1", "V_cousins2", "V_family", "m", "r", "r_f", "pi", "delta",
)

# Starting values (the 2025 starting values, plus the new parameters).
DEFAULT_START = {
    "V_A": 0.5, "V_GxG": 0.25, "V_house": 0.1, "V_house1": 0.1, "V_house2": 0.1,
    "V_house3": 0.1, "V_housepart": 0.1, "V_cousins1": 0.1, "V_cousins2": 0.1,
    "V_family": 0.1, "m": 0.5, "r": 0.5, "r_f": 0.5, "pi": 0.5, "delta": 0.0,
}

# Box bounds; [0, 1] unless listed. delta may be negative so that a paternal excess is
# allowed: the maternal contrast is tested two-sided.
PARAM_BOUNDS = {"delta": (-0.5, 1.0)}

WEIGHT_SCHEMES = ("invSqErr", "invvar", "uniform", "hybrid", "hybrid_invvar", "hybrid_block")
MR_MODES = ("free", "equilibrium", "equal")

# Named component sets used in the 2025 analysis (NHGRI 2026 poster).
MODELS_2025 = {
    "M01": ("A",),
    "M02": ("A", "GxG"),
    "M03a": ("A", "house1"),
    "M03b": ("A", "house2"),
    "M04a": ("A", "house2", "cousins1", "cousins2"),
    "M04b": ("A", "family"),
    "M05a": ("A", "GxG", "house2", "cousins1", "cousins2"),
    "M05b": ("A", "GxG", "family"),
}
MODELS = MODELS_2025   # backwards-compatible name


def _core_models() -> dict[str, tuple[str, ...]]:
    """The a priori model set: genetic {A, A+GxG} x environment {none, household(pi),
    decaying family, household + family} x {without, with maternal inflation delta}."""
    out = {}
    for g in (("A",), ("A", "GxG")):
        for e in ((), ("house",), ("family",), ("house", "family")):
            for mat in ((), ("maternal",)):
                comps = g + e + mat
                name = "+".join({"A": "A", "GxG": "GxG", "house": "H", "family": "F",
                                 "maternal": "mat"}[c] for c in comps)
                out[name] = comps
    return out


CORE_MODELS = _core_models()

@dataclass(frozen=True)
class FitSpec:
    """Everything that defines one fit.

    components : which covariance components are summed (see `components.ALL_COMPONENTS`)
    cohort     : which categories of relatives the loss is computed over (`data.COHORTS`)
    scale      : 'lin' -> residuals on the correlation scale; 'log' -> on ln(correlation)
    weights    : 'invSqErr'      1/SE^2 on whatever scale (Eftedal et al.)
                 'invvar'        1/Var on the scale used: 1/SE^2 for 'lin', corr^2/SE^2 for
                                 'log' (delta method: Var(ln r) ~ SE^2 / r^2)
                 'uniform'       1
                 'hybrid'        every relationship type gets equal total weight, split within
                                 the type in proportion to 1/SE^2 (2025 analysis)
                 'hybrid_invvar' as 'hybrid' but split within type by the scale's 1/Var
                 'hybrid_block'  as 'hybrid' but balanced over cv blocks (near-duplicate types,
                                 e.g. cousins and cousins via UZ twins, share one block's weight)
    mr         : 'free'          m and r estimated separately (Eftedal et al.)
                 'equilibrium'   r estimated, m = V_A * r (primary phenotypic assortment at
                                 equilibrium: corr(g_mate1, g_mate2) = h^2 * r)
                 'equal'         r = m
    stepsib    : 'eftedal' | 'corrected' step-sibling genetic kernel (`components`)
    uz         : 'mixture' | 'exponent' handling of relatives via unknown-zygosity twins
    mz_frac    : assumed fraction of same-sex twins that are MZ (for uz='mixture')
    log_zero   : how rows that cannot be logged are handled on the log scale:
                 'structural' drop rows with corr <= 0 and rows every component structurally
                              predicts to be 0 (e.g. adoptive rows in a genetic-only model)
                 'legacy'     drop rows with corr <= 1e-4, and rows with pred <= 1e-4 at
                              the current parameters (2025 code; makes the loss
                              parameter-dependent). Unlike the 2025 code, only structurally
                              zero rows may be dropped for a small prediction, which removes
                              the 2025 loss's degenerate minimum at V_A = 0.
                 'legacy_raw' exactly the 2025 rule, degenerate minimum included (only for
                              reproducing 2025 fits with a local search from the 2025 start)
    adoptive   : 'zero' | 'legacy_n99' genetic sharing of adoptive pairs (`components`)
    cocousin_family_degree : override f for biological co-cousins (None = annotation table)
    drop_duplicate : drop the duplicated 2nd-cousin row (see data/README.md)
    exclude_types  : relationship types (`rel`) removed from the cohort (leave-one-out)
    max_total_variance : if set, constrain the sum of all variance components V_* to be at
                     most this (1 = phenotypic variance; leaves V_E >= 0). None = box bounds only.
    fixed      : parameters held at fixed values, e.g. (('r', 0.47),)
    """

    components: tuple[str, ...] = ("A",)
    cohort: str = "ALL"
    scale: str = "lin"
    weights: str = "invvar"
    mr: str = "free"
    stepsib: str = "eftedal"
    uz: str = "mixture"
    mz_frac: float = 0.5
    log_zero: str = "structural"
    adoptive: str = "zero"
    cocousin_family_degree: float | None = None
    drop_duplicate: bool = False
    exclude_types: tuple[str, ...] = ()
    max_total_variance: float | None = None
    fixed: tuple[tuple[str, float], ...] = ()
    label: str = field(default="", compare=False)

    def __post_init__(self):
        bad = set(self.components) - set(C.ALL_COMPONENTS)
        if bad:
            raise ValueError(f"unknown components {bad}")
        if self.weights not in WEIGHT_SCHEMES:
            raise ValueError(f"unknown weights {self.weights!r}")
        if self.mr not in MR_MODES:
            raise ValueError(f"unknown mr mode {self.mr!r}")
        if self.scale not in ("lin", "log"):
            raise ValueError(f"unknown scale {self.scale!r}")
        if self.mr == "equilibrium" and "A" not in self.components:
            raise ValueError("mr='equilibrium' needs an additive component")

    def with_(self, **kw) -> "FitSpec":
        return replace(self, **kw)

    def as_dict(self) -> dict:
        d = {k: getattr(self, k) for k in (
            "cohort", "scale", "weights", "mr", "stepsib", "uz", "mz_frac", "log_zero", "adoptive",
            "cocousin_family_degree", "drop_duplicate", "max_total_variance")}
        d["exclude_types"] = ";".join(self.exclude_types)
        d["components"] = "+".join(self.components)
        d["fixed"] = ";".join(f"{k}={v:g}" for k, v in self.fixed)
        d["label"] = self.label
        return d


def required_params(spec: FitSpec) -> list[str]:
    """All parameters the model's prediction depends on (free, fixed or derived)."""
    p = [f"V_{c}" for c in spec.components if c not in C.MODIFIER_COMPONENTS]
    if any(c in C.GENETIC_COMPONENTS for c in spec.components):
        p += ["m", "r"]
    if "family" in spec.components:
        p += ["r_f"]
    if "house" in spec.components:
        p += ["pi"]
    if "maternal" in spec.components:
        p += ["delta"]
    return [q for q in PARAM_ORDER if q in p]


class Model:
    """A `FitSpec` bound to a data table: predictions, weights and loss."""

    MIN_CORR = 1e-4     # legacy log-scale cut-off
    LOG_FLOOR = 1e-10   # floor for predictions on the log scale ('structural')

    def __init__(self, data: pd.DataFrame, spec: FitSpec, corr: np.ndarray | None = None):
        self.spec = spec
        self.data = data.copy()
        if spec.cocousin_family_degree is not None:
            sel = self.data["rel"] == "cocousins"
            self.data.loc[sel, "family_degree"] = spec.cocousin_family_degree
        # `corr` lets bootstrap replicates swap in resampled correlations
        self.corr = data["corr"].to_numpy(dtype=float) if corr is None else np.asarray(corr, float)
        self.se = data["se"].to_numpy(dtype=float)

        self.cohort = cohort_mask(self.data, spec.cohort, spec.drop_duplicate)
        if spec.exclude_types:
            self.cohort &= ~self.data["rel"].isin(spec.exclude_types).to_numpy()
        self.indicators = {c: C.env_indicator(self.data, c)
                           for c in spec.components if c in C.INDICATOR_COMPONENTS}
        self.f = self.data["family_degree"].to_numpy(dtype=float)
        self.kd = C.KernelData(self.data)
        self.house_full, self.house_part = C.household_share(self.data)
        self.mat = C.maternal_indicator(self.data)

        fixed = dict(spec.fixed)
        req = required_params(spec)
        derived = set()
        if spec.mr == "equilibrium":
            derived.add("m")
        elif spec.mr == "equal":
            derived.add("r")
        self.fixed = {k: v for k, v in fixed.items() if k in req}
        self.free = [p for p in req if p not in derived and p not in self.fixed]
        self.required = req
        self.lb = np.array([PARAM_BOUNDS.get(q, (0.0, 1.0))[0] for q in self.free])
        self.ub = np.array([PARAM_BOUNDS.get(q, (0.0, 1.0))[1] for q in self.free])
        self.var_idx = np.array([i for i, q in enumerate(self.free) if q.startswith("V_")], dtype=int)
        self.var_fixed = sum(v for k, v in self.fixed.items() if k.startswith("V_"))

        self._set_rows()

    # ---- rows & weights --------------------------------------------------------------

    def _structurally_zero(self, adoptive_zero: bool | None = None) -> np.ndarray:
        """Rows every component predicts to be exactly 0, whatever the parameters."""
        if adoptive_zero is None:
            adoptive_zero = self.spec.adoptive == "zero"
        zero = np.ones(len(self.data), dtype=bool)
        s = self.data["genetic_structure"].to_numpy()
        for c in self.spec.components:
            if c in C.GENETIC_COMPONENTS:
                zero &= (s == "none") & adoptive_zero
            elif c in C.INDICATOR_COMPONENTS:
                zero &= self.indicators[c] == 0
            elif c == "house":
                zero &= (self.house_full + self.house_part) == 0
            elif c == "family":
                zero &= False
        return zero

    def _set_rows(self):
        rows = self.cohort.copy()
        if self.spec.scale == "log" and self.spec.log_zero == "structural":
            rows &= (self.corr > 0) & ~self._structurally_zero()
        self.rows = rows
        # rows the legacy log-scale rule may drop: those with no genetically/environmentally
        # shared component (adoptive pairs in a genetic-only model). Other rows are always
        # kept; otherwise the legacy loss has a degenerate minimum of 0 at V_A -> 0, where
        # every prediction falls below the cut-off and every row is dropped.
        self._legacy_droppable = self._structurally_zero(adoptive_zero=True)
        self.w = self._weights(self.cohort)

    def _weights(self, mask: np.ndarray) -> np.ndarray:
        """Weights for every row (0 outside `mask`)."""
        se, corr, scale = self.se, self.corr, self.spec.scale
        inv_se2 = 1.0 / se ** 2
        if scale == "log":
            with np.errstate(invalid="ignore"):
                inv_var = np.where(corr > 0, corr ** 2 / se ** 2, 0.0)
        else:
            inv_var = inv_se2
        scheme = self.spec.weights
        if scheme == "invSqErr":
            w = inv_se2
        elif scheme == "invvar":
            w = inv_var
        elif scheme == "uniform":
            w = np.ones_like(se)
        elif scheme in ("hybrid", "hybrid_invvar", "hybrid_block"):
            base = inv_var if scheme == "hybrid_invvar" else inv_se2
            w = np.zeros_like(se)
            col = "cv_block" if scheme == "hybrid_block" else "broad_type"
            types = self.data[col].to_numpy()
            for t in np.unique(types[mask]):
                sel = mask & (types == t)
                tot = base[sel].sum()
                if tot > 0:
                    w[sel] = base[sel] / tot
        w = np.where(mask, w, 0.0)
        return w

    # ---- parameters ------------------------------------------------------------------

    @property
    def constrained(self) -> bool:
        return self.spec.max_total_variance is not None and len(self.var_idx) > 0

    def variance_slack(self, theta) -> float:
        """max_total_variance - sum of variance components (>= 0 when feasible)."""
        return self.spec.max_total_variance - self.var_fixed - float(np.sum(np.asarray(theta)[self.var_idx]))

    def make_feasible(self, theta, margin: float = 1e-3) -> np.ndarray:
        theta = np.clip(np.asarray(theta, float), self.lb, self.ub)
        if self.constrained and self.variance_slack(theta) < 0:
            room = max(self.spec.max_total_variance - self.var_fixed - margin, 0.0)
            tot = theta[self.var_idx].sum()
            theta[self.var_idx] *= room / tot
        return theta

    def start(self) -> np.ndarray:
        return np.array([DEFAULT_START[p] for p in self.free])

    def unpack(self, theta) -> dict:
        p = dict(zip(self.free, np.asarray(theta, float)))
        p.update(self.fixed)
        if self.spec.mr == "equilibrium":
            p["m"] = p["V_A"] * p["r"]
        elif self.spec.mr == "equal":
            p["r"] = p["m"]
        return p

    # ---- predictions & loss ----------------------------------------------------------

    def predict(self, p: Mapping[str, float]) -> np.ndarray:
        spec, pred = self.spec, np.zeros(len(self.data))
        if any(c in C.GENETIC_COMPONENTS for c in spec.components):
            K, K2 = C.genetic_kernels(self.kd, p["m"], p["r"], stepsib=spec.stepsib,
                                      uz=spec.uz, mz_frac=spec.mz_frac,
                                      adoptive=spec.adoptive)
            if "A" in spec.components:
                pred += p["V_A"] * K
            if "GxG" in spec.components:
                pred += p["V_GxG"] * K2
        for c, ind in self.indicators.items():
            pred += p[f"V_{c}"] * ind
        if "house" in spec.components:
            pred += p["V_house"] * (self.house_full + p["pi"] * self.house_part)
        if "family" in spec.components:
            pred += p["V_family"] * p["r_f"] ** self.f
        if "maternal" in spec.components:
            pred = pred * (1.0 + p["delta"] * self.mat)
        return pred

    def components(self, p: Mapping[str, float]) -> pd.DataFrame:
        """Per-row contribution of each component (columns sum to `predict(p)`)."""
        spec, out = self.spec, {}
        if any(c in C.GENETIC_COMPONENTS for c in spec.components):
            K, K2 = C.genetic_kernels(self.kd, p["m"], p["r"], stepsib=spec.stepsib,
                                      uz=spec.uz, mz_frac=spec.mz_frac,
                                      adoptive=spec.adoptive)
            if "A" in spec.components:
                out["A"] = p["V_A"] * K
            if "GxG" in spec.components:
                out["GxG"] = p["V_GxG"] * K2
        for c, ind in self.indicators.items():
            out[c] = p[f"V_{c}"] * ind
        if "house" in spec.components:
            out["house"] = p["V_house"] * (self.house_full + p["pi"] * self.house_part)
        if "family" in spec.components:
            out["family"] = p["V_family"] * p["r_f"] ** self.f
        if "maternal" in spec.components:
            # the inflation applied to the sum of the other components
            base = np.sum(list(out.values()), axis=0) if out else np.zeros(len(self.data))
            out["maternal"] = base * p["delta"] * self.mat
        return pd.DataFrame(out)

    def residuals(self, pred: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """(residual on the fitting scale, rows used); residual is 0 outside rows used."""
        if self.spec.scale == "lin":
            rows = self.rows
            res = np.where(rows, pred - self.corr, 0.0)
        else:
            if self.spec.log_zero == "legacy_raw":
                rows = self.cohort & (self.corr > self.MIN_CORR) & (pred > self.MIN_CORR)
            elif self.spec.log_zero == "legacy":
                rows = (self.cohort & (self.corr > self.MIN_CORR)
                        & ((pred > self.MIN_CORR) | ~self._legacy_droppable))
            else:
                rows = self.rows
            with np.errstate(divide="ignore", invalid="ignore"):
                res = np.where(rows, np.log(np.maximum(pred, self.LOG_FLOOR))
                               - np.log(np.where(rows, self.corr, 1.0)), 0.0)
        return res, rows

    def loss(self, theta) -> float:
        res, _ = self.residuals(self.predict(self.unpack(theta)))
        return float(np.sum(self.w * res ** 2))
