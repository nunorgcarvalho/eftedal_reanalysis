"""Replicate the published Eftedal et al. (2025) estimates from their own analysis recipe.

Writes results/replication/eftedal2025_estimates.csv and eftedal2025_cis.csv.
"""

import argparse
from pathlib import Path

import pandas as pd

from eftedal import data as D
from eftedal import provenance
from eftedal import eftedal2025 as E

PUBLISHED = {
    # Table 1 of the paper
    ("bio_inlaw", "h2"): (0.503, 0.404, 0.572), ("bio_inlaw", "m"): (0.430, 0.365, 0.478),
    ("bio_inlaw", "c"): (0.972, 0.864, 1.127), ("bio_inlaw", "R2"): (0.873, None, None),
    ("bio_only", "h2"): (0.801, 0.749, 0.856), ("bio_only", "m"): (0.225, 0.204, 0.246),
    ("bio_only", "c"): (1.002, 0.896, 1.121), ("bio_only", "R2"): (0.988, None, None),
    # Discussion: classical twin design
    ("twin", "A"): (0.75, None, None), ("twin", "C"): (0.11, None, None),
}


def main(out: Path, n_boot: int):
    out.mkdir(parents=True, exist_ok=True)
    df = D.load("MRE8")

    bio = E.fit_bio_only(df)
    inl = E.fit_bio_inlaw(df)
    twin = E.twin_ace(df)
    est = [
        ("bio_only", "h2", bio.h2), ("bio_only", "m", bio.m), ("bio_only", "c", bio.c),
        ("bio_only", "r_implied", bio.r_implied), ("bio_only", "R2", bio.R2_log),
        ("bio_only", "R2_linear", bio.R2),
        ("bio_inlaw", "h2", inl.h2), ("bio_inlaw", "m", inl.m), ("bio_inlaw", "c", inl.c),
        ("bio_inlaw", "r_implied", inl.r_implied), ("bio_inlaw", "R2", inl.R2),
        ("bio_inlaw", "R2_log", inl.R2_log),
        ("twin", "A", twin["A"]), ("twin", "C", twin["C"]), ("twin", "E", twin["E"]),
    ]
    est = pd.DataFrame(est, columns=["model", "quantity", "ours"])
    est["published"] = [PUBLISHED.get((m, q), (None,))[0] for m, q in zip(est.model, est.quantity)]
    est.to_csv(out / "eftedal2025_estimates.csv", index=False)

    cis = pd.concat([
        bio.ci.assign(model="bio_only"),
        E.bootstrap_bio_inlaw(df, n_boot=n_boot).assign(model="bio_inlaw"),
        E.hessian_ci_bio_inlaw(df).assign(model="bio_inlaw"),
    ], ignore_index=True)
    pub = pd.DataFrame([(m, q, v[1], v[2]) for (m, q), v in PUBLISHED.items() if v[1] is not None],
                       columns=["model", "param", "published_lo", "published_hi"])
    cis = cis.merge(pub, on=["model", "param"], how="left")
    cis.to_csv(out / "eftedal2025_cis.csv", index=False)

    with pd.option_context("display.width", 200, "display.float_format", "{:.4f}".format):
        print(est.to_string(index=False))
        print(cis.to_string(index=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=D.REPO_ROOT / "results" / "replication")
    ap.add_argument("--n-boot", type=int, default=5000)
    a = ap.parse_args()
    main(a.out, a.n_boot)
    provenance.stamp(__file__)
