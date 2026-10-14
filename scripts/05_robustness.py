"""Supplement: robustness rows for the additive model and the in-law contrast for every
phenotype. Writes results/robustness/:

additive_rows.csv     (+ additive_predictions.csv) additive model under the main settings and one-at-a-time departures
                      (log scale with delta-method weights, separate m/r, Eftedal step-sib
                      kernel, block-balanced hybrid, duplicated row dropped, MZ fraction 0.42,
                      Eftedal's own settings), in the full and in Eftedal's bio+in-law cohort
phenotypes.csv        bio-only vs bio + in-law h2 for all 13 phenotypes (MRE8 is main section 2,
                      scripts/02_corrections.py)
"""

import argparse
from pathlib import Path

import pandas as pd

from eftedal import data as D
from eftedal import eftedal2025 as E
from eftedal.fit import fit
from eftedal.presets import main_spec

PHENOTYPES = ["MRE8", "MRE5", "MR5", "MR8", "MR9", "math5", "math8", "math9",
              "read5", "read8", "read9", "eng5", "eng8"]

# Eftedal et al.'s settings expressed in our framework (log scale, 1/SE^2 on the log scale,
# m and r separate, their step-sibling kernel).
EFTEDAL_SETTINGS = dict(scale="log", mr="free", stepsib="eftedal")   # with weights="invSqErr"

# Variants of the additive model: (name, weights, overrides)
VARIANTS = [
    ("main", "hybrid_block", {}),
    ("main", "invvar", {}),
    ("log scale, delta-method weights", "invvar", {"scale": "log"}),
    ("log scale, block hybrid", "hybrid_block", {"scale": "log"}),
    ("m, r separate", "hybrid_block", {"mr": "free"}),
    ("m, r separate", "invvar", {"mr": "free"}),
    ("Eftedal step-sib kernel", "hybrid_block", {"stepsib": "eftedal"}),
    ("Eftedal step-sib kernel", "invvar", {"stepsib": "eftedal"}),
    ("duplicated row dropped", "hybrid_block", {"drop_duplicate": True}),
    ("duplicated row dropped", "invvar", {"drop_duplicate": True}),
    ("MZ fraction 0.42", "hybrid_block", {"mz_frac": 0.42}),
    ("MZ fraction 0.42", "invvar", {"mz_frac": 0.42}),
    ("Eftedal settings", "invSqErr", EFTEDAL_SETTINGS),
]

# (name, weights, overrides, bio-only cohort, bio + in-law cohort)
CONTRAST_SETTINGS = [
    ("Eftedal settings", "invSqErr", EFTEDAL_SETTINGS, "bio_only", "ALL_but_adopted_UZ"),
    ("main", "hybrid_block", {}, "bio", "ALL_but_adopted"),
    ("main", "invvar", {}, "bio", "ALL_but_adopted"),
]


def _row(res, **tag):
    d = res.summary()
    d.update(tag)
    return d


def _preds(res, i, **tag):
    p = res.predictions()
    p.insert(0, "fit_id", i)
    return p.assign(**tag)


def additive_rows(df, n_starts):
    rows, preds = [], []
    for cohort in ("ALL", "ALL_but_adopted_UZ"):
        for name, w, kw in VARIANTS:
            res = fit(df, main_spec("A", w, cohort=cohort, **kw), n_starts=n_starts)
            tag = dict(variant=name, weighting=w, cohort_name=cohort)
            rows.append(_row(res, fit_id=len(rows), **tag))
            preds.append(_preds(res, len(rows) - 1, **tag))
    a = E.fit_bio_inlaw(df)
    rows.append({"variant": "Eftedal et al. code (exact)", "weighting": "invSqErr",
                 "cohort_name": "ALL_but_adopted_UZ", "V_A": a.h2, "m": a.m, "c": a.c,
                 "r": a.r_implied, "R2_lin": a.R2})
    return pd.DataFrame(rows), pd.concat(preds, ignore_index=True)


def contrast(df, n_starts, phenotype="MRE8"):
    rows = []
    for name, w, kw, c_bio, c_all in CONTRAST_SETTINGS:
        for part, cohort in (("bio only", c_bio), ("bio + in-law", c_all)):
            res = fit(df, main_spec("A", w, cohort=cohort, **kw), n_starts=n_starts)
            rows.append(_row(res, setting=name, weighting=w, part=part, cohort_name=cohort,
                             phenotype=phenotype))
    return pd.DataFrame(rows)


def main(out: Path, n_starts: int):
    out.mkdir(parents=True, exist_ok=True)
    df = D.load("MRE8")
    rows, preds = additive_rows(df, n_starts)
    rows.to_csv(out / "additive_rows.csv", index=False)
    preds.to_csv(out / "additive_predictions.csv", index=False)
    pd.concat([contrast(D.load(p), n_starts, p) for p in PHENOTYPES],
              ignore_index=True).to_csv(out / "phenotypes.csv", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=D.REPO_ROOT / "results" / "robustness")
    ap.add_argument("--n-starts", type=int, default=12)
    a = ap.parse_args()
    main(a.out, a.n_starts)
