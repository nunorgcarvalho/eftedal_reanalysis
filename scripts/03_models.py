"""Main section 3 (core models) and the joint model set (supplement): fits, block
cross-validation and block-jackknife intervals, under both weightings. Writes results/models/
with one file group per model set ({set} = core | joint | core_am2):

{set}_fits.csv          one row per model x weighting: parameters, V_sum, loss, R2_lin, chi2/df
{set}_predictions.csv   predictions and per-component contributions for all 82 categories
{set}_jackknife.csv     estimate and block-jackknife 95% interval per model x weighting x parameter
                        (incl. V_sum, the sum of the variance components)
{set}_cv_params.csv     parameters of every leave-one-block-out refit (also the jackknife replicates)
{set}_cv_heldout.csv    held-out categories with their CV predictions
{set}_cv_scores.csv     CV score per model x weighting x fold (scoring rule: eftedal.cv)
{set}_cv_summary.csv    overall CV score per model x weighting
{set}_cv_paired.csv     paired fold-by-fold differences from the best model per weighting
"""

import argparse
from pathlib import Path

import pandas as pd

from eftedal import cv
from eftedal import data as D
from eftedal.model import CORE_MODELS, JOINT_MODELS
from eftedal.presets import MAIN_WEIGHTINGS, main_spec

# set name -> (models, FitSpec overrides)
SETS = {"core": (list(CORE_MODELS), {}),
        "joint": (list(JOINT_MODELS), {}),
        # supplement: the earlier GxG kernel, square of the assortative-mating kernel
        "core_am2": (list(CORE_MODELS), {"gxg_kernel": "am2"})}
PARAMS = ["V_A", "V_GxG", "V_house", "V_family", "V_sum", "V_A0", "pi", "r_f", "m", "r"]


def run_set(df, name, models, overrides, out, n_starts):
    fits, preds, jks, params, held, scores, summ = [], [], [], [], [], [], []
    for w in MAIN_WEIGHTINGS:
        for m in models:
            full, fp, ho, sc = cv.cross_validate(df, main_spec(m, w, **overrides), n_starts=n_starts)
            tag = dict(model=m, weights=w)
            fits.append({**tag, **full.summary()})
            preds.append(full.predictions().assign(**tag))
            fp = fp.assign(**tag)
            params.append(fp)
            ps = [p for p in PARAMS if p in full.params or p in ("V_sum", "V_A0")]
            jks.append(cv.jackknife(fp, full.params, ps,
                                    blocks=[f for f in fp["fold"] if f in D.CV_BLOCKS]).assign(**tag))
            held.append(ho.assign(**tag))
            scores.append(sc.assign(**tag))
            summ.append({**tag, **cv.summarize_scores(sc), "chi2_full": full.stats["chi2"],
                         "n_free": len(full.free)})
            print(name, w, m, f"cv_rmse={summ[-1]['cv_rmse']:.4f}", flush=True)
    write = lambda rows, f: (pd.concat(rows, ignore_index=True) if isinstance(rows[0], pd.DataFrame)
                             else pd.DataFrame(rows)).to_csv(out / f"{name}_{f}.csv", index=False)
    write(fits, "fits")
    write(preds, "predictions")
    write(jks, "jackknife")
    write(params, "cv_params")
    write(held, "cv_heldout")
    write(scores, "cv_scores")
    write(summ, "cv_summary")
    cv.paired_comparison(pd.concat(scores, ignore_index=True)).to_csv(out / f"{name}_cv_paired.csv", index=False)


def main(out: Path, n_starts: int, sets):
    out.mkdir(parents=True, exist_ok=True)
    df = D.load("MRE8")
    for name in sets:
        run_set(df, name, *SETS[name], out, n_starts)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=D.REPO_ROOT / "results" / "models")
    ap.add_argument("--n-starts", type=int, default=20)
    ap.add_argument("--sets", nargs="+", default=list(SETS), choices=list(SETS))
    a = ap.parse_args()
    main(a.out, a.n_starts, a.sets)
