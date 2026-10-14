"""Main sections 1 and 2: from Eftedal et al.'s estimator to ours, and the effect of in-laws.
Every estimate carries a block-jackknife 95% interval (leave one cv block of relationship types
out, refit; only blocks present in the fitted cohort). Writes results/corrections/:

steps.csv            additive model (A), Eftedal's estimator and our corrections applied
                     cumulatively (the cumulative order is the order of STEPS), plus one aside
                     (Eftedal's log scale with delta-method weights); one row per step x parameter
steps_jackknife.csv  the leave-one-block-out estimates behind steps.csv
inlaw.csv            h2 of the additive model fitted to biological relatives only vs. with
                     in-laws, and the difference, per estimator; the interval on the difference
                     is a paired jackknife (the same block left out of both fits)
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from eftedal import cv
from eftedal import data as D
from eftedal.presets import main_spec

REPORT = ["V_A", "m", "r"]

# Eftedal et al.'s estimator in our framework: log scale, 1/SE^2 weights on the log scale,
# m and r separate, their step-sibling kernel, bio + in-law cohort without UZ-linked categories.
EFTEDAL = dict(scale="log", weights="invSqErr", mr="free", stepsib="eftedal",
               cohort="ALL_but_adopted_UZ")

# (key, label, settings): cumulative corrections, then the change of weighting.
STEPS = [
    ("eftedal", "Eftedal et al. (log scale, 1/SE²)", EFTEDAL),
    ("linear", "+ linear scale (1/SE² = inverse variance)", {**EFTEDAL, "scale": "lin", "weights": "invvar"}),
    ("mr", "+ m = h²r", {**EFTEDAL, "scale": "lin", "weights": "invvar", "mr": "equilibrium"}),
    ("stepsib", "+ corrected step-sibling kernel",
     {**EFTEDAL, "scale": "lin", "weights": "invvar", "mr": "equilibrium", "stepsib": "corrected"}),
    ("all82", "+ all 82 categories (= ours, inverse variance)", dict(weights="invvar")),
    ("block", "ours, block-balanced hybrid weights", dict(weights="hybrid_block")),
]
ASIDE = [("eftedal_delta", "Eftedal, but delta-method weights on the log scale",
          {**EFTEDAL, "weights": "invvar"})]

# In-law contrast: (key, label, bio-only cohort, bio + in-law cohort, settings)
INLAW = [
    ("eftedal", "Eftedal et al.", "bio_only", "ALL_but_adopted_UZ",
     {k: v for k, v in EFTEDAL.items() if k != "cohort"}),
    ("invvar", "ours, inverse variance", "bio", "ALL_but_adopted", dict(weights="invvar")),
    ("hybrid_block", "ours, block-balanced hybrid", "bio", "ALL_but_adopted", dict(weights="hybrid_block")),
]


def _jk_rows(df, spec, n_starts, **tag):
    full, fp, _, sc = cv.cross_validate(df, spec, n_starts=n_starts,
                                       folds={b: D.rels_in_blocks(df, b) for b in D.CV_BLOCKS})
    jk = cv.jackknife(fp, full.params, REPORT, blocks=list(fp["fold"]))
    jk = jk.assign(**tag, R2_lin=full.stats["R2_lin"], chi2_per_df=full.stats["chi2_per_df"])
    return full, fp.assign(**tag), jk


def steps(df, n_starts):
    rows, jks = [], []
    for i, (key, label, kw) in enumerate(STEPS + ASIDE):
        _, fp, jk = _jk_rows(df, main_spec("A", **kw), n_starts, step=key, label=label,
                             order=i, aside=key in {a[0] for a in ASIDE})
        rows.append(jk)
        jks.append(fp)
        print("step", key, flush=True)
    return pd.concat(rows, ignore_index=True), pd.concat(jks, ignore_index=True)


def inlaw(df, n_starts):
    rows = []
    for key, label, c_bio, c_all, kw in INLAW:
        fb, fpb, jb = _jk_rows(df, main_spec("A", cohort=c_bio, **kw), n_starts)
        fa, fpa, ja = _jk_rows(df, main_spec("A", cohort=c_all, **kw), n_starts)
        # paired jackknife of the difference over the blocks of the larger cohort; a block absent
        # from the bio-only cohort leaves that fit unchanged
        b_bio = fpb.set_index("fold")["V_A"]
        reps = np.array([b_bio.get(b, fb.params["V_A"]) - v
                         for b, v in fpa.set_index("fold")["V_A"].items()])
        B = len(reps)
        se_gap = np.sqrt((B - 1) / B * np.sum((reps - reps.mean()) ** 2))
        gap = fb.params["V_A"] - fa.params["V_A"]
        for what, est, se in (("bio only", fb.params["V_A"], jb.set_index("param").loc["V_A", "se"]),
                              ("bio + in-law", fa.params["V_A"], ja.set_index("param").loc["V_A", "se"]),
                              ("difference", gap, se_gap)):
            rows.append({"estimator": key, "label": label, "quantity": what, "estimate": est,
                         "se": se, "lo": est - 1.96 * se, "hi": est + 1.96 * se})
        print("inlaw", key, flush=True)
    return pd.DataFrame(rows)


def main(out: Path, n_starts: int):
    out.mkdir(parents=True, exist_ok=True)
    df = D.load("MRE8")
    s, j = steps(df, n_starts)
    s.to_csv(out / "steps.csv", index=False)
    j.to_csv(out / "steps_jackknife.csv", index=False)
    inlaw(df, n_starts).to_csv(out / "inlaw.csv", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=D.REPO_ROOT / "results" / "corrections")
    ap.add_argument("--n-starts", type=int, default=12)
    a = ap.parse_args()
    main(a.out, a.n_starts)
