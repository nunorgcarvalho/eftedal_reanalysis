"""Fitting a `Model`: multi-start bounded optimisation of the weighted least-squares loss."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from .model import FitSpec, Model, PARAM_ORDER


@dataclass
class FitResult:
    spec: FitSpec
    model: Model
    theta: np.ndarray            # free parameters at the optimum
    params: dict                 # all parameters (free, fixed and derived)
    loss: float                  # weighted RSS on the fitting scale (spec weights, unnormalised)
    pred: np.ndarray             # predicted correlation for all 82 rows
    rows: np.ndarray             # rows that entered the loss
    n_starts: int
    n_converged_to_best: int     # how many starts reached (within tol) the best loss
    stats: dict = field(default_factory=dict)

    @property
    def free(self):
        return self.model.free

    def at_bound(self, tol=1e-6) -> list[str]:
        return [p for p, t in zip(self.free, self.theta)
                if t <= self.model.lb[self.free.index(p)] + tol
                or t >= self.model.ub[self.free.index(p)] - tol]

    def summary(self) -> dict:
        d = {"spec_label": self.spec.label, **self.spec.as_dict()}
        for p in PARAM_ORDER:
            d[p] = self.params.get(p, np.nan)
        d["c"] = half_multiplier_or_nan(self.params)
        d["V_sum"] = sum(v for k, v in self.params.items() if k.startswith("V_"))
        d.update(self.stats)
        d["loss"] = self.loss
        d["n_free"] = len(self.free)
        d["free"] = ",".join(self.free)
        d["at_bound"] = ",".join(self.at_bound())
        d["n_starts"] = self.n_starts
        d["n_converged_to_best"] = self.n_converged_to_best
        return d

    def predictions(self) -> pd.DataFrame:
        df = self.model.data[["row", "rel", "side", "sex", "corr", "se", "group"]].copy()
        df["pred"] = self.pred
        df["in_cohort"] = self.model.cohort
        df["in_loss"] = self.rows
        df["weight"] = self.model.w
        comps = self.model.components(self.params)
        for c in comps:
            df[f"comp_{c}"] = comps[c].to_numpy()
        return df


def half_multiplier_or_nan(p: dict) -> float:
    if "m" in p and "r" in p:
        m, r = p["m"], p["r"]
        return (1 + 2 * m + r * m) / (1 + m) ** 2
    return np.nan


def fit_stats(model: Model, pred: np.ndarray) -> dict:
    """Goodness-of-fit summaries comparable across specs.

    R2_fit  : 1 - RSS/TSS on the fitting scale with the fitting weights over the rows used
              (the 2025 analysis' R^2)
    R2_lin  : the same with linear-scale 1/SE^2 weights over the whole cohort (Eftedal et
              al.'s R^2; a common yardstick across weighting schemes and scales)
    chi2    : sum over the cohort of ((obs - pred) / SE)^2; if the model were right and the
              correlations independent this would be ~ chi^2 with df = n_cohort - n_free
    """
    res, rows = model.residuals(pred)
    w = model.w
    if model.spec.scale == "log":
        y = np.log(np.where(rows, model.corr, 1.0))
    else:
        y = model.corr
    if np.sum(w[rows]) > 0:
        ybar = np.sum(w[rows] * y[rows]) / np.sum(w[rows])
        rss = np.sum(w[rows] * res[rows] ** 2)
        tss = np.sum(w[rows] * (y[rows] - ybar) ** 2)
    else:   # degenerate legacy fit with every row dropped
        rss, tss = np.nan, np.nan

    c = model.cohort
    wl = 1.0 / model.se[c] ** 2
    yl = model.corr[c]
    rss_l = np.sum(wl * (pred[c] - yl) ** 2)
    tss_l = np.sum(wl * (yl - np.sum(wl * yl) / np.sum(wl)) ** 2)
    n = int(c.sum())
    df = n - len(model.free)
    return {
        "R2_fit": 1 - rss / tss,
        "R2_lin": 1 - rss_l / tss_l,
        "chi2": rss_l,
        "df": df,
        "chi2_per_df": rss_l / df if df > 0 else np.nan,
        "n_cohort": n,
        "n_rows_in_loss": int(rows.sum()),
    }


def _random_starts(model: Model, n: int, rng: np.random.Generator) -> list[np.ndarray]:
    k = len(model.free)
    return [model.make_feasible(model.lb + (model.ub - model.lb) * rng.uniform(0.02, 0.98, size=k))
            for _ in range(n)]


def _local_search(model: Model, x0: np.ndarray, scale: float) -> tuple[np.ndarray, float]:
    """Local search from x0: L-BFGS-B, or SLSQP when the total-variance constraint is on."""
    f = lambda th: model.loss(th) / scale
    x0 = model.make_feasible(x0)
    bounds = list(zip(model.lb, model.ub))
    if model.constrained:
        r = minimize(f, x0, method="SLSQP", bounds=bounds,
                     constraints=[{"type": "ineq", "fun": model.variance_slack}],
                     options={"maxiter": 2000, "ftol": 1e-15})
        x = np.clip(r.x, model.lb, model.ub)
        ok = model.variance_slack(x) >= -1e-9
        r_fun = f(x) if ok else np.inf
    else:
        r = minimize(f, x0, method="L-BFGS-B", bounds=bounds,
                     options={"maxiter": 5000, "ftol": 1e-15, "gtol": 1e-12})
        x, r_fun = r.x, r.fun
    fx0 = f(x0)
    return (x, r_fun) if np.isfinite(r_fun) and r_fun <= fx0 else (x0, fx0)


def _polish(model: Model, x0: np.ndarray, scale: float) -> tuple[np.ndarray, float]:
    """Derivative-free polish (robust to the kinks of the legacy log-scale row selection).
    Under the total-variance constraint, violations are penalised and the result is
    finished with the constrained local search."""
    def f(th):
        pen = 0.0
        if model.constrained:
            pen = 1e6 * min(model.variance_slack(th), 0.0) ** 2
        return model.loss(th) / scale + pen
    bounds = list(zip(model.lb, model.ub))
    best_x, best_f = x0, f(x0)
    for method, opts in (("Powell", {"maxiter": 20000, "xtol": 1e-10, "ftol": 1e-14}),
                         ("Nelder-Mead", {"maxiter": 20000, "xatol": 1e-10, "fatol": 1e-15})):
        r = minimize(f, best_x, method=method, bounds=bounds, options=opts)
        if np.isfinite(r.fun) and r.fun < best_f:
            best_x, best_f = r.x, r.fun
    return _local_search(model, best_x, scale)


def fit(data: pd.DataFrame, spec: FitSpec, *, n_starts: int = 12, seed: int = 0,
        starts: list[np.ndarray] | None = None, corr: np.ndarray | None = None,
        polish: bool = True) -> FitResult:
    """Fit `spec` to `data`, keeping the best of several starting points.

    The default start (2025 starting values) is always tried, plus `n_starts` random
    starts (and any explicit `starts`, e.g. a warm start for bootstrap replicates). The two
    best local solutions are polished with derivative-free searches. With n_starts=0 and
    polish=False only a local search from the default start is run.
    """
    if spec.sequential_base is not None:
        return _fit_sequential(data, spec, n_starts=n_starts, seed=seed, starts=starts,
                               corr=corr, polish=polish)
    model = Model(data, spec, corr=corr)
    if not model.free:
        theta = np.array([])
        p = model.unpack(theta)
        pred = model.predict(p)
        _, rows = model.residuals(pred)
        return FitResult(spec, model, theta, p, model.loss(theta), pred, rows, 0, 0,
                         fit_stats(model, pred))

    rng = np.random.default_rng(seed)
    x0s = [model.start()] + list(starts or []) + _random_starts(model, n_starts, rng)
    if not polish and n_starts == 0 and not starts:
        # purely local search from the default start (used to mimic the 2025 runs)
        theta, loss = _local_search(model, x0s[0], max(np.sum(model.w), 1e-300))
        theta = np.clip(theta, model.lb, model.ub)
        p = model.unpack(theta)
        pred = model.predict(p)
        _, rows = model.residuals(pred)
        return FitResult(spec, model, theta, p, model.loss(theta), pred, rows, 1, 1,
                         fit_stats(model, pred))
    scale = max(np.sum(model.w), 1e-300)   # optimise a normalised loss; same argmin
    outs = [_local_search(model, x0, scale) for x0 in x0s]
    if polish:
        for i in np.argsort([o[1] for o in outs])[:2]:
            outs.append(_polish(model, outs[i][0], scale))
    losses = np.array([o[1] for o in outs])
    best = int(np.argmin(losses))
    theta = np.clip(outs[best][0], model.lb, model.ub)
    n_best = int(np.sum(losses[:len(x0s)] <= losses[best] * (1 + 1e-6) + 1e-12))
    p = model.unpack(theta)
    pred = model.predict(p)
    _, rows = model.residuals(pred)
    return FitResult(spec, model, theta, p, model.loss(theta), pred, rows,
                     len(x0s), n_best, fit_stats(model, pred))


def _fit_sequential(data, spec: FitSpec, **kw) -> FitResult:
    """Stage 1: the base model; stage 2: the full model with the carried parameters fixed at
    their stage-1 values (sequential_carry=("all",): every parameter estimated in stage 1). Returns the stage-2 result, labelled with the original (sequential)
    spec so that refits (bootstrap, jackknife, CV) redo both stages. Its `model` is the stage-2
    model, so Wald/sandwich covariances are conditional on stage 1."""
    starts = kw.pop("starts", None)
    base = fit(data, spec.with_(components=spec.sequential_base, sequential_base=None), **kw)
    names = base.free if tuple(spec.sequential_carry) == ("all",) else spec.sequential_carry
    carried = tuple((p, float(base.params[p])) for p in names if p in base.params)
    stage2 = spec.with_(sequential_base=None, fixed=spec.fixed + carried)
    res = fit(data, stage2, starts=starts, **kw)
    res.spec = spec
    # degrees of freedom also spend the parameters estimated in stage 1
    n_carried = len(carried)
    res.stats["df"] -= n_carried
    res.stats["chi2_per_df"] = res.stats["chi2"] / res.stats["df"] if res.stats["df"] > 0 else np.nan
    res.stats["n_params_estimated"] = len(res.free) + n_carried
    res.stats.update({f"stage1_{p}": v for p, v in base.params.items()})
    res.stats["stage1_loss"] = base.loss
    return res


def fit_many(data: pd.DataFrame, specs: list[FitSpec], **kw) -> tuple[pd.DataFrame, pd.DataFrame, list[FitResult]]:
    """Fit a list of specs; return (summary table, long predictions table, results)."""
    results = [fit(data, s, **kw) for s in specs]
    summary = pd.DataFrame([r.summary() for r in results])
    summary.insert(0, "fit_id", range(len(results)))
    preds = []
    for i, r in enumerate(results):
        d = r.predictions()
        d.insert(0, "fit_id", i)
        preds.append(d)
    return summary, pd.concat(preds, ignore_index=True), results
