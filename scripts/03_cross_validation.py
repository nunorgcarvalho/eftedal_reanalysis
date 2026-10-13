"""Block cross-validation of the 16 core models under both weightings (see `eftedal.cv` for
the folds and the scoring rule). Writes results/cv/:

cv_fold_params.csv    parameters of every refit (model x weighting x fold)
cv_heldout.csv        held-out categories with their CV predictions (and full-fit predictions)
cv_fold_scores.csv    score per model x weighting x fold
cv_summary.csv        overall score per model x weighting (single-block folds), plus the
                      in-sample score of the full fit on the same scale for comparison
"""

import argparse
from pathlib import Path

import pandas as pd

from eftedal import cv
from eftedal import data as D
from eftedal.model import CORE_MODELS
from eftedal.presets import MAIN_WEIGHTINGS, main_spec


def main(out: Path, n_starts: int):
    out.mkdir(parents=True, exist_ok=True)
    df = D.load("MRE8")
    params, held, scores, summ = [], [], [], []
    for w in MAIN_WEIGHTINGS:
        for m in CORE_MODELS:
            full, fp, ho, sc = cv.cross_validate(df, main_spec(m, w), n_starts=n_starts)
            tag = dict(model=m, weights=w)
            params.append(fp.assign(**tag))
            held.append(ho.assign(**tag))
            scores.append(sc.assign(**tag))
            insample = cv.type_mse(df.assign(pred=full.pred)[full.model.cohort])
            insample = insample.merge(df[["rel", "cv_block"]].drop_duplicates(), on="rel")
            in_score = insample.groupby("cv_block")["mse"].mean().mean()
            summ.append({**tag, **cv.summarize_scores(sc), "insample_score": in_score,
                         "n_free": len(full.free), "chi2_full": full.stats["chi2"]})
            print(w, m, f"{summ[-1]['cv_rmse']:.4f}", flush=True)
    pd.concat(params, ignore_index=True).to_csv(out / "cv_fold_params.csv", index=False)
    pd.concat(held, ignore_index=True).to_csv(out / "cv_heldout.csv", index=False)
    pd.concat(scores, ignore_index=True).to_csv(out / "cv_fold_scores.csv", index=False)
    pd.DataFrame(summ).to_csv(out / "cv_summary.csv", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=D.REPO_ROOT / "results" / "cv")
    ap.add_argument("--n-starts", type=int, default=20)
    a = ap.parse_args()
    main(a.out, a.n_starts)
