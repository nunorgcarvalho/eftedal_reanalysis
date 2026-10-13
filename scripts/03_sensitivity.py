"""Sensitivity analyses under corrected settings. Writes results/sensitivity/:

fisher_{summary,predictions}.csv   additive-only model: cohort x fit setting x m-r mode x
                                   step-sibling formula
models_{summary,predictions}.csv   all named models on all 82 categories: fit setting x m-r
                                   mode x step-sibling formula
loo_summary.csv                    leave-one-relationship-type-out refits of key specs
extra_summary.csv                  one-at-a-time checks: duplicated row dropped, co-cousin
                                   family degree 3, MZ fraction 0.42, legacy UZ handling
"""

import argparse
from pathlib import Path

import pandas as pd

from eftedal import data as D
from eftedal.fit import fit, fit_many
from eftedal.presets import fisher_grid, key_specs, model_grid

KEY_SPECS = key_specs()


def write_grid(specs, out: Path, stem: str, n_starts: int):
    summary, preds, _ = fit_many(D.load("MRE8"), specs, n_starts=n_starts)
    summary.to_csv(out / f"{stem}_summary.csv", index=False)
    preds.to_csv(out / f"{stem}_predictions.csv", index=False)
    return summary


def leave_one_out(df, n_starts):
    rows = []
    for k, s in enumerate(KEY_SPECS):
        base = fit(df, s, n_starts=n_starts)
        types = df.loc[base.model.cohort, "rel"].unique()
        for t in [None, *types]:
            ss = s if t is None else s.with_(exclude_types=(t,))
            r = base if t is None else fit(df, ss, n_starts=n_starts, starts=[base.theta])
            d = r.summary()
            d.update(key_id=k, key_label=s.label, left_out=t or "(none)")
            rows.append(d)
    return pd.DataFrame(rows)


def one_at_a_time(df, n_starts):
    variants = {
        "baseline": {},
        "drop_duplicate_row": {"drop_duplicate": True},
        "cocousin_f=3": {"cocousin_family_degree": 3},
        "mz_frac=0.42": {"mz_frac": 0.42},
        "uz=exponent": {"uz": "exponent"},
    }
    rows = []
    for k, s in enumerate(KEY_SPECS):
        for name, kw in variants.items():
            if name == "cocousin_f=3" and "family" not in s.components:
                continue
            r = fit(df, s.with_(**kw), n_starts=n_starts)
            d = r.summary()
            d.update(key_id=k, key_label=s.label, variant=name)
            rows.append(d)
    return pd.DataFrame(rows)


def main(out: Path, n_starts: int):
    out.mkdir(parents=True, exist_ok=True)
    df = D.load("MRE8")
    write_grid(fisher_grid(), out, "fisher", n_starts)
    write_grid(model_grid(), out, "models", n_starts)
    leave_one_out(df, n_starts).to_csv(out / "loo_summary.csv", index=False)
    one_at_a_time(df, n_starts).to_csv(out / "extra_summary.csv", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=D.REPO_ROOT / "results" / "sensitivity")
    ap.add_argument("--n-starts", type=int, default=12)
    a = ap.parse_args()
    main(a.out, a.n_starts)
