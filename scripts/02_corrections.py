"""Main sections 1 and 2: from Eftedal et al.'s estimator to ours, and the effect of in-laws.
Every estimate carries two 95% intervals (column `method`):
  block jackknife       ("marginal type" CI) leave one cv block of relationship types out, refit
                        (blocks present in the fitted cohort); sensitivity to which kinds of relatives exist
  parametric bootstrap  ("marginal pair" CI) correlations redrawn from N(rho, SE^2), refit; sampling
                        noise only. For the in-law difference the same draws are used for both fits (paired).
Writes results/corrections/:

steps.csv            additive model (A): Eftedal's estimator, Eftedal's with corrected (delta-method)
                     weights on the log scale, Eftedal's on the linear scale, and our revised
                     estimator with marginal-pair (1/SE^2) and marginal-type (block-balanced) weights;
                     one row per step x parameter x interval type
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
from eftedal import provenance
from eftedal import data as D
from eftedal.inference import parametric_bootstrap
from eftedal.presets import main_spec

REPORT = ["V_A", "V_A0", "m", "r"]

# Eftedal et al.'s estimator in our framework: log scale, 1/SE^2 weights on the log scale,
# m and r separate, their step-sibling kernel, bio + in-law cohort without UZ-linked categories.
EFTEDAL = dict(scale="log", weights="invSqErr", mr="free", stepsib="eftedal",
               cohort="ALL_but_adopted_UZ")

# (key, label, settings), in display order. 'revised_pair' applies every correction at once
# (linear scale with 1/SE^2 weights, m = h2 r, corrected step-sibling kernel, all 82 categories).
STEPS = [
    ("eftedal", "Eftedal et al.", EFTEDAL),
    ("eftedal_delta", "Eftedal et al. (corrected inverse variance weight)", {**EFTEDAL, "weights": "invvar"}),
    ("linear", "Eftedal et al. + linear scale", {**EFTEDAL, "scale": "lin", "weights": "invvar"}),
    ("revised_pair", "Revised (Marginal Pair)", dict(weights="invvar")),
    ("revised_type", "Revised (Marginal Type)", dict(weights="hybrid_block")),
]
ASIDE = []

# In-law contrast: (key, label, bio-only cohort, bio + in-law cohort, settings)
INLAW = [
    ("eftedal", "Eftedal et al.", "bio_only", "ALL_but_adopted_UZ",
     {k: v for k, v in EFTEDAL.items() if k != "cohort"}),
    ("invvar", "Revised (Marginal Pair)", "bio", "ALL_but_adopted", dict(weights="invvar")),
    ("hybrid_block", "Revised (Marginal Type)", "bio", "ALL_but_adopted", dict(weights="hybrid_block")),
]


def _jk_rows(df, spec, n_starts, n_boot, seed, **tag):
    """Full fit, jackknife replicates, and a table with both intervals; also the bootstrap
    replicates (for paired differences)."""
    full, fp, _, sc = cv.cross_validate(df, spec, n_starts=n_starts,
                                       folds={b: D.rels_in_blocks(df, b) for b in D.CV_BLOCKS})
    jk = cv.jackknife(fp, full.params, REPORT, blocks=list(fp["fold"]))
    reps = parametric_bootstrap(df, full, n_boot=n_boot, seed=seed, n_starts=2)
    est = {**full.params, "V_A0": full.summary()["V_A0"]}
    bt = pd.DataFrame([{"param": p, "estimate": est[p], "se": reps[p].std(),
                        "lo": reps[p].quantile(0.025), "hi": reps[p].quantile(0.975),
                        "method": "parametric bootstrap"} for p in REPORT])
    t = pd.concat([jk, bt], ignore_index=True).assign(
        **tag, R2_lin=full.stats["R2_lin"], chi2_per_df=full.stats["chi2_per_df"])
    return full, fp.assign(**tag), t, reps


def steps(df, n_starts, n_boot):
    rows, jks = [], []
    for i, (key, label, kw) in enumerate(STEPS + ASIDE):
        _, fp, jk, _ = _jk_rows(df, main_spec("A", **kw), n_starts, n_boot, 10 + i, step=key,
                                label=label, order=i, aside=key in {a[0] for a in ASIDE})
        rows.append(jk)
        jks.append(fp)
        print("step", key, flush=True)
    return pd.concat(rows, ignore_index=True), pd.concat(jks, ignore_index=True)


def inlaw(df, n_starts, n_boot):
    rows = []
    for i, (key, label, c_bio, c_all, kw) in enumerate(INLAW):
        # same seed -> same redrawn correlations for both fits (paired bootstrap)
        fb, fpb, jb, rb = _jk_rows(df, main_spec("A", cohort=c_bio, **kw), n_starts, n_boot, 50 + i)
        fa, fpa, ja, ra = _jk_rows(df, main_spec("A", cohort=c_all, **kw), n_starts, n_boot, 50 + i)
        # paired jackknife of the difference over the blocks of the larger cohort; a block absent
        # from the bio-only cohort leaves that fit unchanged
        b_bio = fpb.set_index("fold")["V_A"]
        reps = np.array([b_bio.get(b, fb.params["V_A"]) - v
                         for b, v in fpa.set_index("fold")["V_A"].items()])
        B = len(reps)
        se_gap = np.sqrt((B - 1) / B * np.sum((reps - reps.mean()) ** 2))
        gap = fb.params["V_A"] - fa.params["V_A"]
        pick = lambda t, m: t[(t.param == "V_A") & (t.method == m)].iloc[0]
        for what, est, jrow, brow in (("bio only", fb.params["V_A"], pick(jb, "block jackknife"), pick(jb, "parametric bootstrap")),
                                      ("bio + in-law", fa.params["V_A"], pick(ja, "block jackknife"), pick(ja, "parametric bootstrap"))):
            for r in (jrow, brow):
                rows.append({"estimator": key, "label": label, "quantity": what, "estimate": est,
                             "se": r["se"], "lo": r["lo"], "hi": r["hi"], "method": r["method"]})
        rows.append({"estimator": key, "label": label, "quantity": "difference", "estimate": gap,
                     "se": se_gap, "lo": gap - cv.jackknife_quantile(B) * se_gap, "hi": gap + cv.jackknife_quantile(B) * se_gap,
                     "method": "block jackknife"})
        d = rb["V_A"].to_numpy() - ra["V_A"].to_numpy()
        rows.append({"estimator": key, "label": label, "quantity": "difference", "estimate": gap,
                     "se": d.std(ddof=1), "lo": np.quantile(d, 0.025), "hi": np.quantile(d, 0.975),
                     "method": "parametric bootstrap"})
        print("inlaw", key, flush=True)
    return pd.DataFrame(rows)


def main(out: Path, n_starts: int, n_boot: int):
    out.mkdir(parents=True, exist_ok=True)
    df = D.load("MRE8")
    s, j = steps(df, n_starts, n_boot)
    s.to_csv(out / "steps.csv", index=False)
    j.to_csv(out / "steps_jackknife.csv", index=False)
    inlaw(df, n_starts, n_boot).to_csv(out / "inlaw.csv", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=D.REPO_ROOT / "results" / "corrections")
    ap.add_argument("--n-starts", type=int, default=12)
    ap.add_argument("--n-boot", type=int, default=300)
    a = ap.parse_args()
    main(a.out, a.n_starts, a.n_boot)
    provenance.stamp(__file__)
