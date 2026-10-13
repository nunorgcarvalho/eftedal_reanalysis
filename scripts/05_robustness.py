"""Robustness rows, the bio-only vs bio+in-law contrast, all phenotypes, and the 2025 model
set under the main settings. Writes results/robustness/:

additive_rows.csv     additive model under the main settings and one-at-a-time departures
                      (log scale with delta-method weights, separate m/r, Eftedal step-sib
                      kernel, block-balanced hybrid, duplicated row dropped, MZ fraction 0.42,
                      Eftedal's own settings), in the full and in Eftedal's bio+in-law cohort
contrast.csv          bio-only vs bio+in-law h2 under Eftedal's settings and ours (headline)
phenotypes.csv        the contrast for all 13 phenotypes
models_2025.csv       the 2025 model set (M01-M05b) under the main settings (supplement)
"""

import argparse
from pathlib import Path

import pandas as pd

from eftedal import data as D
from eftedal import eftedal2025 as E
from eftedal.fit import fit
from eftedal.model import MODELS_2025
from eftedal.presets import MAIN_WEIGHTINGS, main_spec

PHENOTYPES = ["MRE8", "MRE5", "MR5", "MR8", "MR9", "math5", "math8", "math9",
              "read5", "read8", "read9", "eng5", "eng8"]

# Eftedal et al.'s settings expressed in our framework (log scale, 1/SE^2 on the log scale,
# m and r separate, their step-sibling kernel).
EFTEDAL_SETTINGS = dict(scale="log", mr="free", stepsib="eftedal")   # with weights="invSqErr"

# Variants of the additive model: (name, weights, overrides)
VARIANTS = [
    ("main", "hybrid", {}),
    ("main", "invvar", {}),
    ("block-balanced hybrid", "hybrid_block", {}),
    ("log scale, delta-method weights", "invvar", {"scale": "log"}),
    ("log scale, hybrid", "hybrid", {"scale": "log"}),
    ("m, r separate", "hybrid", {"mr": "free"}),
    ("m, r separate", "invvar", {"mr": "free"}),
    ("Eftedal step-sib kernel", "hybrid", {"stepsib": "eftedal"}),
    ("Eftedal step-sib kernel", "invvar", {"stepsib": "eftedal"}),
    ("duplicated row dropped", "hybrid", {"drop_duplicate": True}),
    ("duplicated row dropped", "invvar", {"drop_duplicate": True}),
    ("MZ fraction 0.42", "hybrid", {"mz_frac": 0.42}),
    ("MZ fraction 0.42", "invvar", {"mz_frac": 0.42}),
    ("Eftedal settings", "invSqErr", EFTEDAL_SETTINGS),
]

CONTRAST_SETTINGS = [
    ("Eftedal settings", "invSqErr", EFTEDAL_SETTINGS),
    ("main", "hybrid", {}),
    ("main", "invvar", {}),
    ("block-balanced hybrid", "hybrid_block", {}),
    ("log scale, delta-method weights", "invvar", {"scale": "log"}),
]


def _row(res, **tag):
    d = res.summary()
    d.update(tag)
    return d


def additive_rows(df, n_starts):
    rows = []
    for cohort in ("ALL", "ALL_but_adopted_UZ"):
        for name, w, kw in VARIANTS:
            res = fit(df, main_spec("A", w, cohort=cohort, **kw), n_starts=n_starts)
            rows.append(_row(res, variant=name, weighting=w, cohort_name=cohort))
    a = E.fit_bio_inlaw(df)
    rows.append({"variant": "Eftedal et al. code (exact)", "weighting": "invSqErr",
                 "cohort_name": "ALL_but_adopted_UZ", "V_A": a.h2, "m": a.m, "c": a.c,
                 "r": a.r_implied, "R2_lin": a.R2})
    return pd.DataFrame(rows)


def contrast(df, n_starts, phenotype="MRE8"):
    rows = []
    for name, w, kw in CONTRAST_SETTINGS:
        for cohort in ("bio_only", "ALL_but_adopted_UZ"):
            res = fit(df, main_spec("A", w, cohort=cohort, **kw), n_starts=n_starts)
            rows.append(_row(res, setting=name, weighting=w, cohort_name=cohort, phenotype=phenotype))
    return pd.DataFrame(rows)


def models_2025(df, n_starts):
    rows = []
    for w in MAIN_WEIGHTINGS:
        for m in MODELS_2025:
            res = fit(df, main_spec(m, w), n_starts=n_starts)
            rows.append(_row(res, model=m, weighting=w))
    return pd.DataFrame(rows)


def main(out: Path, n_starts: int):
    out.mkdir(parents=True, exist_ok=True)
    df = D.load("MRE8")
    additive_rows(df, n_starts).to_csv(out / "additive_rows.csv", index=False)
    contrast(df, n_starts).to_csv(out / "contrast.csv", index=False)
    pd.concat([contrast(D.load(p), n_starts, p) for p in PHENOTYPES],
              ignore_index=True).to_csv(out / "phenotypes.csv", index=False)
    models_2025(df, n_starts).to_csv(out / "models_2025.csv", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=D.REPO_ROOT / "results" / "robustness")
    ap.add_argument("--n-starts", type=int, default=12)
    a = ap.parse_args()
    main(a.out, a.n_starts)
