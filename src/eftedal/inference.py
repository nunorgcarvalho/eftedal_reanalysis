"""Uncertainty for fitted models.

The observed correlations are treated as independent with known sampling variances SE^2
(ignoring pair overlap between categories; Eftedal et al. estimate CIs are <30% too
narrow). With weights W and fitting-scale residuals e(theta), the estimator minimises
e' W e. Under that sampling model:

* `sandwich`  - A^-1 B A^-1 with A = half the Hessian of the loss, B = J'W S W J and
                S = diag(Var of each observation on the fitting scale). Valid for any weights,
                and it targets the weighted-projection ("pseudo-true") parameter even when the
                model is misspecified. Only sampling noise in the correlations is counted.
                (`sandwich_gn` uses the Gauss-Newton A = J'WJ instead, which drops the
                residual-curvature term and can be far off under misspecification.)
* `quasi`     - lm/nls convention: sigma^2 (J'WJ)^-1 with sigma^2 = RSS/(n - p). Treats the
                weights as relative precisions and absorbs lack of fit into the noise. With
                inverse-variance weights it equals the model-based covariance times chi2/df.
* `param_boot`- parametric bootstrap: redraw every correlation from N(corr, SE^2), refit,
                percentile intervals. Same target as the sandwich, but respects the [0, 1]
                bounds and nonlinearity.
* `hessian_rss` - inverse Hessian of the weighted RSS, as in Eftedal et al.'s code. NOT a
                covariance (off by 2 for 1/SE^2 weights; arbitrary for other weightings),
                kept for replication only.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from .fit import FitResult, fit


def numerical_jacobian(f, x, h=1e-6):
    x = np.asarray(x, float)
    f0 = np.atleast_1d(f(x))
    J = np.empty((len(f0), len(x)))
    for j in range(len(x)):
        e = np.zeros_like(x)
        e[j] = h * max(1.0, abs(x[j]))
        J[:, j] = (np.atleast_1d(f(x + e)) - np.atleast_1d(f(x - e))) / (2 * e[j])
    return J


def numerical_hessian(f, x, h=1e-4):
    x = np.asarray(x, float)
    k = len(x)
    H = np.empty((k, k))
    steps = h * np.maximum(1.0, np.abs(x))
    for i in range(k):
        for j in range(i, k):
            ei = np.zeros(k); ei[i] = steps[i]
            ej = np.zeros(k); ej[j] = steps[j]
            H[i, j] = H[j, i] = (f(x + ei + ej) - f(x + ei - ej) - f(x - ei + ej)
                                 + f(x - ei - ej)) / (4 * steps[i] * steps[j])
    return H


def _transformed_pred(res: FitResult, theta) -> np.ndarray:
    m = res.model
    pred = m.predict(m.unpack(theta))
    if m.spec.scale == "log":
        return np.log(np.maximum(pred, m.LOG_FLOOR))
    return pred


def _obs_var(res: FitResult) -> np.ndarray:
    m = res.model
    if m.spec.scale == "log":
        with np.errstate(divide="ignore"):
            return np.where(m.corr > 0, m.se ** 2 / np.maximum(m.corr, 1e-300) ** 2, np.inf)
    return m.se ** 2


def _jacobian(res: FitResult) -> np.ndarray:
    """d(prediction on fitting scale)/d(free params), one-sided at active bounds."""
    m, th = res.model, res.theta
    J = np.empty((len(m.data), len(th)))
    base = _transformed_pred(res, th)
    for j in range(len(th)):
        h = 1e-6
        up, dn = th.copy(), th.copy()
        if th[j] + h <= m.ub[j] and th[j] - h >= m.lb[j]:
            up[j] += h; dn[j] -= h
            J[:, j] = (_transformed_pred(res, up) - _transformed_pred(res, dn)) / (2 * h)
        elif th[j] + h <= m.ub[j]:
            up[j] += h
            J[:, j] = (_transformed_pred(res, up) - base) / h
        else:
            dn[j] -= h
            J[:, j] = (base - _transformed_pred(res, dn)) / h
    return J


def covariance(res: FitResult, method: str = "sandwich") -> np.ndarray:
    """Covariance of the free parameters (see module docstring for `method`)."""
    m = res.model
    rows = res.rows
    J = _jacobian(res)[rows]
    w = m.w[rows]
    A = J.T @ (w[:, None] * J)
    Ainv = np.linalg.pinv(A)
    if method in ("sandwich", "sandwich_gn"):
        # The estimate solves J'W e = 0, so d(theta)/d(obs) = H^-1 J'W with H half the full
        # Hessian of the loss. Gauss-Newton (H ~ J'WJ) drops the residual-curvature term,
        # which is not negligible when the model is misspecified (large residuals).
        S = _obs_var(res)[rows]
        B = J.T @ ((w ** 2 * S)[:, None] * J)
        if method == "sandwich":
            Ainv = np.linalg.pinv(0.5 * numerical_hessian(m.loss, res.theta))
        return Ainv @ B @ Ainv
    if method == "quasi":
        e, _ = m.residuals(res.pred)
        dof = rows.sum() - len(res.theta)
        sigma2 = np.sum(w * e[rows] ** 2) / dof
        return sigma2 * Ainv
    if method == "hessian_rss":
        H = numerical_hessian(m.loss, res.theta)
        return np.linalg.pinv(H)
    raise ValueError(method)


def wald_table(res: FitResult, methods=("sandwich", "sandwich_gn", "quasi", "hessian_rss"),
               level: float = 0.95) -> pd.DataFrame:
    """Wald intervals (estimate +/- z SE) for the free parameters, clipped to the bounds."""
    z = stats.norm.ppf(0.5 + level / 2)
    out = []
    for meth in methods:
        cov = covariance(res, meth)
        se = np.sqrt(np.clip(np.diag(cov), 0, None))
        for p, est, s, lb, ub in zip(res.free, res.theta, se, res.model.lb, res.model.ub):
            out.append({"param": p, "estimate": est, "se": s, "lo": max(est - z * s, lb),
                        "hi": min(est + z * s, ub), "lo_raw": est - z * s, "hi_raw": est + z * s,
                        "method": meth})
    return pd.DataFrame(out)


def parametric_bootstrap(data: pd.DataFrame, res: FitResult, n_boot: int = 500, seed: int = 0,
                         se_inflation: float = 1.0, n_starts: int = 2) -> pd.DataFrame:
    """Refit to correlations redrawn from N(corr, (se_inflation*SE)^2). Returns one row per
    replicate with every parameter (free, fixed and derived) plus R2_lin."""
    rng = np.random.default_rng(seed)
    corr0 = data["corr"].to_numpy(float)
    se = data["se"].to_numpy(float) * se_inflation
    reps = []
    for b in range(n_boot):
        corr = corr0 + rng.normal(size=len(corr0)) * se
        r = fit(data, res.spec, n_starts=n_starts, seed=seed + 1000 + b, starts=[res.theta],
                corr=corr, polish=False)
        row = {k: v for k, v in r.params.items()}
        row["c"] = (1 + 2 * r.params["m"] + r.params["r"] * r.params["m"]) / (1 + r.params["m"]) ** 2 \
            if "m" in r.params else np.nan
        row["rep"] = b
        row["loss"] = r.loss
        reps.append(row)
    return pd.DataFrame(reps)


def percentile_table(reps: pd.DataFrame, estimate: dict, level: float = 0.95,
                     method: str = "param_boot") -> pd.DataFrame:
    a = (1 - level) / 2
    out = []
    for p in estimate:
        if p in reps:
            out.append({"param": p, "estimate": estimate[p], "se": reps[p].std(),
                        "lo": reps[p].quantile(a), "hi": reps[p].quantile(1 - a), "method": method})
    return pd.DataFrame(out)
