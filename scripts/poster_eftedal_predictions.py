"""Supplemental-slides input: per-category predictions of the additive model under Eftedal et al.'s
estimator (their settings, biological relatives + in-laws), with the additive contribution split
into its random-mating part and the excess under assortative mating (as in core_predictions.csv).
Intervals for the inset come from results/poster/additive_by_cohort.csv (scripts/poster_additive.py).
Writes results/poster/eftedal_predictions.csv.
"""

import argparse
import importlib.util
from pathlib import Path

from eftedal import data as D
from eftedal import provenance
from eftedal.fit import fit
from eftedal.presets import main_spec

spec = importlib.util.spec_from_file_location("corr", Path(__file__).with_name("02_corrections.py"))
corr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(corr)


def main(out: Path, n_starts: int):
    out.mkdir(parents=True, exist_ok=True)
    df = D.load("MRE8")
    key, label, _, cohort, kw = next(r for r in corr.INLAW if r[0] == "eftedal")
    res = fit(df, main_spec("A", modifiers=(), cohort=cohort, **kw), n_starts=n_starts)
    print(label, {k: round(v, 4) for k, v in res.params.items()})
    res.predictions().assign(model=label, cohort=cohort).to_csv(out / "eftedal_predictions.csv", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=D.REPO_ROOT / "results" / "poster")
    ap.add_argument("--n-starts", type=int, default=12)
    a = ap.parse_args()
    main(a.out, a.n_starts)
    provenance.stamp(__file__)
