"""Supplement: sampling correlation between the parameters of each core model, from the
leave-one-type-out (block jackknife) refits. Writes results/models/:

core_param_corr.csv   per model x weighting x parameter pair: r = correlation of the 10 single-block
                      refit estimates (= the block-jackknife estimate of the sampling correlation of
                      the two estimates), and an approximate 95% interval from jackknifing r itself:
                      r_(i) is recomputed from the refits that leave out block i and one more block,
                      SE on Fisher's z = sqrt((B-1)/B sum_i (z_(i) - mean)^2), interval
                      tanh(atanh(r) +/- 1.96 SE). NaN when a parameter does not move across refits.
core_cv2_params.csv   parameters of the leave-two-blocks-out refits

The single-block refits are read from results/models/core_cv_params.csv (scripts/03_models.py).
"""

import argparse
import itertools
from pathlib import Path

import numpy as np
import pandas as pd

from eftedal import data as D
from eftedal import provenance
from eftedal.fit import fit
from eftedal.model import CORE_MODELS
from eftedal.presets import MAIN_WEIGHTINGS, main_spec

# the model's independent parameters as reported (m rather than r: r = m / V_A at equilibrium)
PARAMS = ["V_A", "m", "V_GxG", "V_family", "r_f"]


def corr(x, y):
    if np.std(x) < 1e-8 or np.std(y) < 1e-8:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def main(out: Path, n_starts: int):
    df = D.load("MRE8")
    folds = {b: D.rels_in_blocks(df, b) for b in D.CV_BLOCKS}
    cv1 = pd.read_csv(out / "core_cv_params.csv")
    rows, cv2 = [], []
    for w in MAIN_WEIGHTINGS:
        for m in CORE_MODELS:
            spec = main_spec(m, w)
            full = fit(df, spec, n_starts=20)
            ps = [p for p in PARAMS if p in full.params]
            one = cv1[(cv1.model == m) & (cv1.weights == w) & cv1.fold.isin(D.CV_BLOCKS)].set_index("fold")
            blocks = list(one.index)
            two = {}
            for a, b in itertools.combinations(blocks, 2):
                r = fit(df, spec.with_(exclude_types=tuple(folds[a]) + tuple(folds[b])),
                        n_starts=n_starts, starts=[full.theta])
                two[frozenset((a, b))] = r.params
                cv2.append({**r.params, "fold1": a, "fold2": b, "model": m, "weights": w})
            B = len(blocks)
            for p1, p2 in itertools.combinations(ps, 2):
                r_all = corr(one[p1].to_numpy(), one[p2].to_numpy())
                z = []
                for i in blocks:
                    th = [two[frozenset((i, j))] for j in blocks if j != i]
                    z.append(np.arctanh(np.clip(corr(np.array([t[p1] for t in th]),
                                                     np.array([t[p2] for t in th])), -0.999, 0.999)))
                z = np.array(z)
                se = np.sqrt((B - 1) / B * np.sum((z - z.mean()) ** 2)) if np.isfinite(z).all() else np.nan
                zr = np.arctanh(np.clip(r_all, -0.999, 0.999))
                rows.append({"model": m, "weights": w, "p1": p1, "p2": p2, "r": r_all,
                             "lo": np.tanh(zr - 1.96 * se), "hi": np.tanh(zr + 1.96 * se),
                             "se_z": se, "n_blocks": B})
            print(w, m, flush=True)
    pd.DataFrame(rows).to_csv(out / "core_param_corr.csv", index=False)
    pd.DataFrame(cv2).to_csv(out / "core_cv2_params.csv", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=D.REPO_ROOT / "results" / "models")
    ap.add_argument("--n-starts", type=int, default=6)
    a = ap.parse_args()
    main(a.out, a.n_starts)
    provenance.stamp(__file__)
