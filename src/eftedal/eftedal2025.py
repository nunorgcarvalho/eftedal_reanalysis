"""Replication of the analysis code published with Eftedal et al. (2025)
("analysis_ beyond additive.R"), so their numbers can be reproduced in the same pipeline.

Their Fisherian model is fitted on the log scale:
    ln(rho) = b0 + b1 * n + b2 * d_half + b3 * d_inlaw
with h2 = exp(b0), (1+m)/2 = exp(b1), c = exp(b2) (half-relative multiplier) and m = exp(b3).
Each log-correlation is weighted by 1/SE^2 of the *untransformed* correlation.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize

from .data import cohort_mask
from .inference import numerical_hessian, numerical_jacobian


def twin_ace(df: pd.DataFrame) -> dict:
    """Classical twin design from precision-weighted MZ and DZ correlations."""
    def wmean(rel):
        d = df[df["rel"] == rel]
        w = 1 / d["se"] ** 2
        return float(np.sum(w * d["corr"]) / np.sum(w))
    rmz, rdz = wmean("mz"), wmean("dz")
    return {"r_mz": rmz, "r_dz": rdz, "A": 2 * (rmz - rdz), "C": 2 * rdz - rmz, "E": 1 - rmz}


def predict(df: pd.DataFrame, h2: float, m: float, c: float, mz_frac: float = 0.5) -> np.ndarray:
    """Their predicted correlations (their R script's `dt$pr`; adoptive pairs predicted 0)."""
    n = df["n"].to_numpy(dtype=float)
    s = df["genetic_structure"].to_numpy()
    t = lambda k: ((1 + m) / 2) ** k
    pred = np.zeros(len(df))
    pred[s == "full"] = h2 * t(n[s == "full"])
    u = s == "full_uz"
    pred[u] = mz_frac * h2 * t(n[u] - 0.5) + (1 - mz_frac) * h2 * t(n[u] + 0.5)
    pred[s == "half"] = h2 * c * t(n[s == "half"])
    inlaw = (s == "step") | (s == "cocousin")
    pred[inlaw] = h2 * m * t(n[inlaw])
    return pred


def weighted_r2(obs, pred, w) -> float:
    rss = np.sum(w * (obs - pred) ** 2)
    tss = np.sum(w * (obs - np.sum(w * obs) / np.sum(w)) ** 2)
    return float(1 - rss / tss)


@dataclass
class EftedalFit:
    cohort: str
    h2: float
    m: float
    c: float
    r_implied: float          # r implied by c and m: c = (1+2m+rm)/(1+m)^2
    R2: float                 # their R^2: linear scale, 1/SE^2 weights, over the cohort
    R2_log: float             # R^2 of the log-scale weighted regression itself
    b: np.ndarray             # parameters on their optimisation scale
    ci: pd.DataFrame | None = None


def implied_r(m, c):
    return (c * (1 + m) ** 2 - 1 - 2 * m) / m


def fit_bio_only(df: pd.DataFrame) -> EftedalFit:
    """Bio-only model: weighted linear regression of ln(corr) on n and d_half.

    m is identified from the slope alone (no in-law rows); 95% CIs are the usual lm
    t-intervals (residual variance estimated from the data), transformed as they did.
    """
    d = df[cohort_mask(df, "bio_only")]
    y = np.log(d["corr"].to_numpy())
    X = np.column_stack([np.ones(len(d)), d["n"].to_numpy(float),
                         (d["genetic_structure"] == "half").to_numpy(float)])
    w = 1 / d["se"].to_numpy() ** 2
    sw = np.sqrt(w)
    b, *_ = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)
    res = y - X @ b
    dof = len(y) - X.shape[1]
    sigma2 = np.sum(w * res ** 2) / dof
    cov = sigma2 * np.linalg.inv(X.T @ (w[:, None] * X))
    se = np.sqrt(np.diag(cov))
    q = stats.t.ppf(0.975, dof)
    lo, hi = b - q * se, b + q * se
    h2, m, c = np.exp(b[0]), 2 * np.exp(b[1]) - 1, np.exp(b[2])
    ci = pd.DataFrame({
        "param": ["h2", "m", "c"],
        "estimate": [h2, m, c],
        "lo": [np.exp(lo[0]), 2 * np.exp(lo[1]) - 1, np.exp(lo[2])],
        "hi": [np.exp(hi[0]), 2 * np.exp(hi[1]) - 1, np.exp(hi[2])],
        "method": "lm t-interval",
    })
    pred = predict(d, h2, m, c)
    R2 = weighted_r2(d["corr"].to_numpy(), pred, w)
    R2_log = 1 - np.sum(w * res ** 2) / np.sum(w * (y - np.sum(w * y) / np.sum(w)) ** 2)
    return EftedalFit("bio_only", h2, m, c, implied_r(m, c), R2, R2_log, b, ci)


def _inlaw_design(d: pd.DataFrame):
    s = d["genetic_structure"]
    return (np.log(d["corr"].to_numpy()), d["n"].to_numpy(float),
            (s == "half").to_numpy(float), s.isin(["step", "cocousin"]).to_numpy(float))


def _inlaw_obj(b, y, n, half, inlaw, w):
    """Their objective: b = (b0, ril = ln m, half = ln c); the slope is ln((1+m)/2)."""
    b0, ril, bh = b
    slope = np.log((np.exp(ril) + 1) / 2)
    pred = b0 + slope * n + bh * half + ril * inlaw
    return np.sum(w * (y - pred) ** 2)


def _minimize_inlaw(y, n, half, inlaw, w, x0=(-0.5, -0.5, 0.0), quick=False):
    f = lambda b: _inlaw_obj(b, y, n, half, inlaw, w)
    if quick:   # warm-started bootstrap replicates: the objective is smooth
        return minimize(f, np.array(x0), method="BFGS", options={"gtol": 1e-9}).x
    best = minimize(f, np.array(x0), method="Nelder-Mead",
                    options={"xatol": 1e-12, "fatol": 1e-12, "maxiter": 20000})
    r = minimize(f, best.x, method="BFGS", options={"gtol": 1e-10})
    return r.x if r.fun <= best.fun else best.x


def fit_bio_inlaw(df: pd.DataFrame, weights: np.ndarray | None = None) -> EftedalFit:
    """Bio + in-law model (their headline fit): biological relatives, step sibs and co-cousins,
    excluding adoptive pairs and classes linked through twins of unknown zygosity."""
    d = df[cohort_mask(df, "ALL_but_adopted_UZ")]
    y, n, half, inlaw = _inlaw_design(d)
    w = 1 / d["se"].to_numpy() ** 2 if weights is None else weights
    b = _minimize_inlaw(y, n, half, inlaw, w)
    h2, m, c = np.exp(b)
    pred = predict(d, h2, m, c)
    R2 = weighted_r2(d["corr"].to_numpy(), pred, 1 / d["se"].to_numpy() ** 2)
    yhat = b[0] + np.log((m + 1) / 2) * n + b[2] * half + b[1] * inlaw
    R2_log = 1 - np.sum(w * (y - yhat) ** 2) / np.sum(w * (y - np.sum(w * y) / np.sum(w)) ** 2)
    return EftedalFit("ALL_but_adopted_UZ", h2, m, c, implied_r(m, c), R2, R2_log, b)


def bootstrap_bio_inlaw(df: pd.DataFrame, n_boot: int = 5000, seed: int = 1) -> pd.DataFrame:
    """Their bootstrap: resample the categories with probability proportional to 1/SE^2,
    refit, and take percentile intervals of exp(b) (= h2, m, c)."""
    d = df[cohort_mask(df, "ALL_but_adopted_UZ")].reset_index(drop=True)
    y, n, half, inlaw = _inlaw_design(d)
    w = 1 / d["se"].to_numpy() ** 2
    rng = np.random.default_rng(seed)
    p = w / w.sum()
    x0 = _minimize_inlaw(y, n, half, inlaw, w)
    out = np.empty((n_boot, 3))
    for i in range(n_boot):
        idx = rng.choice(len(d), size=len(d), replace=True, p=p)
        out[i] = _minimize_inlaw(y[idx], n[idx], half[idx], inlaw[idx], w[idx], x0=x0, quick=True)
    est = np.exp(x0)
    q = np.exp(np.quantile(out, [0.025, 0.975], axis=0))
    return pd.DataFrame({"param": ["h2", "m", "c"], "estimate": est, "lo": q[0], "hi": q[1],
                         "method": "case bootstrap (prob ~ 1/SE^2)"})


def hessian_ci_bio_inlaw(df: pd.DataFrame) -> pd.DataFrame:
    """Their alternative CIs: inverse Hessian of the weighted RSS (treated as a covariance),
    delta-method to exp(b), +/- 2 SE.

    The inverse Hessian of a sum of squares is not a covariance matrix: for 1/SE^2 weights
    the Gauss-Newton covariance is 2 * H^-1, and with any other weight normalisation the
    result scales arbitrarily. Kept only to reproduce their numbers.
    """
    d = df[cohort_mask(df, "ALL_but_adopted_UZ")]
    y, n, half, inlaw = _inlaw_design(d)
    w = 1 / d["se"].to_numpy() ** 2
    b = _minimize_inlaw(y, n, half, inlaw, w)
    H = numerical_hessian(lambda x: _inlaw_obj(x, y, n, half, inlaw, w), b)
    cov = np.linalg.inv(H)
    J = numerical_jacobian(np.exp, b)
    se = np.sqrt(np.diag(J @ cov @ J.T))
    est = np.exp(b)
    return pd.DataFrame({"param": ["h2", "m", "c"], "estimate": est, "se": se,
                         "lo": est - 2 * se, "hi": est + 2 * se, "method": "inverse Hessian of RSS"})
