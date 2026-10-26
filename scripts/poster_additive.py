"""Poster figure input: the additive model under Eftedal et al.'s estimator and the revised
estimator (Marginal Pair, Marginal Type), each fitted to biological relatives only and to biological
relatives + in-laws, with h2 (V_A), m and r and both 95% intervals (Marginal Pair CI = parametric
bootstrap, Marginal Type CI = leave-one-type-out jackknife). Reuses scripts/02_corrections.py.
Writes results/poster/additive_by_cohort.csv.
"""

import argparse
import importlib.util
from pathlib import Path

import pandas as pd

from eftedal import data as D
from eftedal import provenance
from eftedal.presets import main_spec

spec = importlib.util.spec_from_file_location("corr", Path(__file__).with_name("02_corrections.py"))
corr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(corr)


def main(out: Path, n_starts: int, n_boot: int):
    out.mkdir(parents=True, exist_ok=True)
    df = D.load("MRE8")
    rows = []
    for i, (key, label, c_bio, c_all, kw) in enumerate(corr.INLAW):
        for cohort_name, cohort in (("bio only", c_bio), ("bio + in-laws", c_all)):
            _, _, t, _ = corr._jk_rows(df, main_spec("A", cohort=cohort, **kw), n_starts, n_boot, 50 + i,
                                       estimator=key, label=label, cohort_name=cohort_name, cohort=cohort)
            rows.append(t)
            print(key, cohort_name, flush=True)
    pd.concat(rows, ignore_index=True).to_csv(out / "additive_by_cohort.csv", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=D.REPO_ROOT / "results" / "poster")
    ap.add_argument("--n-starts", type=int, default=12)
    ap.add_argument("--n-boot", type=int, default=300)
    a = ap.parse_args()
    main(a.out, a.n_starts, a.n_boot)
    provenance.stamp(__file__)
