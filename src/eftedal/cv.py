"""Block cross-validation: refit without a block of relationship types, predict the block.

Scoring (fixed before looking at results):

* per held-out relationship type t: mse_t = sum_k w_k e_k^2 / sum_k w_k over its categories,
  with w = 1/SE^2 and e = observed - predicted (the within-type part of the hybrid weights);
* block score = mean of mse_t over the block's types (every type counts equally);
* overall score = mean of the block scores over the 10 single-block folds (every block counts
  equally); reported also as RMSE = sqrt(score);
* held-out chi2 = sum over held-out categories of (e / SE)^2 (SE-standardised alternative).

The group tests (all in-laws / all adoptive; 'all twins' is the twins block) are scored the same way but
reported separately. The 10 single-block refits also serve as the block jackknife.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import data as D
from .fit import derived_params, fit, variance_sum
from .model import FitSpec


def type_mse(df: pd.DataFrame) -> pd.DataFrame:
    """Within-type 1/SE^2-weighted MSE of held-out predictions (df has rel, corr, se, pred)."""
    d = df.assign(w=1 / df["se"] ** 2, e2=(df["corr"] - df["pred"]) ** 2)
    g = d.groupby("rel", sort=False)
    return pd.DataFrame({"mse": g.apply(lambda x: np.sum(x.w * x.e2) / np.sum(x.w), include_groups=False),
                         "chi2": g.apply(lambda x: np.sum(x.e2 * x.w), include_groups=False),
                         "n_cat": g.size()}).reset_index()


def cross_validate(data: pd.DataFrame, spec: FitSpec, *, n_starts: int = 12,
                   folds: dict[str, tuple[str, ...]] | None = None, seed: int = 0):
    """Run every fold for one spec.

    Returns (full_fit, fold_params, heldout_predictions, fold_scores):
      fold_params : one row per fold with every parameter of the refit
      heldout     : held-out categories with their predictions
      fold_scores : one row per fold: score (mean type mse), rmse, chi2, n types / categories
    """
    folds = D.cv_folds(data) if folds is None else folds
    # only folds that remove something from the spec's cohort
    in_cohort = (D.cohort_mask(data, spec.cohort, spec.drop_duplicate)
                 & ~data["rel"].isin(spec.exclude_types).to_numpy())
    present = set(data.loc[in_cohort, "rel"])
    folds = {k: v for k, v in folds.items() if present & set(v)}
    full = fit(data, spec, n_starts=n_starts, seed=seed)
    params, held, scores = [], [], []
    for name, rels in folds.items():
        r = fit(data, spec.with_(exclude_types=tuple(rels)), n_starts=n_starts, seed=seed,
                starts=[full.theta])
        p = dict(r.params)
        p["V_sum"] = variance_sum(r.params)
        p.update(derived_params(r.params))
        p.update(fold=name, loss=r.loss, n_converged_to_best=r.n_converged_to_best,
                 at_bound=",".join(r.at_bound()))
        params.append(p)
        out = data.loc[data["rel"].isin(rels), ["row", "rel", "side", "sex", "corr", "se", "cv_block"]].copy()
        out["pred"] = r.pred[out["row"].to_numpy()]
        out["pred_full"] = full.pred[out["row"].to_numpy()]
        out["fold"] = name
        held.append(out)
        tm = type_mse(out)
        scores.append({"fold": name, "group_test": name in D.CV_GROUP_TESTS,
                       "score": tm["mse"].mean(), "rmse": np.sqrt(tm["mse"].mean()),
                       "chi2": tm["chi2"].sum(), "n_types": len(tm), "n_cat": int(tm["n_cat"].sum())})
    return full, pd.DataFrame(params), pd.concat(held, ignore_index=True), pd.DataFrame(scores)


def summarize_scores(scores: pd.DataFrame) -> dict:
    """Overall CV summary from a fold_scores table (single-block folds only)."""
    s = scores[~scores["group_test"]]
    return {"cv_score": s["score"].mean(), "cv_rmse": np.sqrt(s["score"].mean()),
            "cv_chi2": s["chi2"].sum(), "n_folds": len(s)}


def jackknife(fold_params: pd.DataFrame, full_params: dict, params: list[str],
              blocks: list[str] | None = None) -> pd.DataFrame:
    """Block-jackknife SEs from the single-block refits:
    SE = sqrt((B - 1) / B * sum_b (theta_(b) - mean)^2); interval = estimate +/- 1.96 SE.
    'V_sum' (sum of the variance components) can be requested as a parameter."""
    blocks = D.CV_BLOCKS if blocks is None else blocks
    fp = fold_params[fold_params["fold"].isin(blocks)]     # single-block folds present
    full_params = dict(full_params)
    full_params.setdefault("V_sum", variance_sum(full_params))
    for k, v in derived_params(full_params).items():
        full_params.setdefault(k, v)
    B = len(fp)
    rows = []
    for p in params:
        th = fp[p].to_numpy(dtype=float)
        se = np.sqrt((B - 1) / B * np.sum((th - th.mean()) ** 2))
        est = full_params[p]
        rows.append({"param": p, "estimate": est, "se": se, "lo": est - 1.96 * se,
                     "hi": est + 1.96 * se, "jk_min": th.min(), "jk_max": th.max(),
                     "method": "block jackknife", "n_blocks": B})
    return pd.DataFrame(rows)


def paired_comparison(fold_scores: pd.DataFrame, by=("weights", "model")) -> pd.DataFrame:
    """Paired fold-by-fold comparison of every model with the best one (lowest mean score)
    within each weighting, over the single-block folds: mean difference in score and its
    standard error across the folds (the folds are the replicates)."""
    s = fold_scores[~fold_scores["group_test"]]
    rows = []
    for w, g in s.groupby(by[0]):
        wide = g.pivot_table(index="fold", columns=by[1], values="score")
        best = wide.mean().idxmin()
        for m in wide.columns:
            d = wide[m] - wide[best]
            rows.append({by[0]: w, by[1]: m, "best": best, "mean_diff": d.mean(),
                         "se_diff": d.std(ddof=1) / np.sqrt(len(d)),
                         "folds_worse": int((d > 0).sum()), "n_folds": len(d)})
    return pd.DataFrame(rows).sort_values([by[0], "mean_diff"])
