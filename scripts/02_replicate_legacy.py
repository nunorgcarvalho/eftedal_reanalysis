"""Reproduce the 2025 analysis (NHGRI 2026 poster fits and the sensitivity write-up) with the
new code, and compare against the reference values from the 2025 R code.

For each 2025 fit this reports the 2025 loss evaluated by the new code at the 2025 parameters
(checks the model code) and the loss/parameters of a fresh multi-start fit (checks the
optimiser). Writes results/legacy/*.csv.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from eftedal import data as D
from eftedal.fit import fit, fit_many
from eftedal.model import Model, PARAM_ORDER
from eftedal.presets import nhgri_specs, sensitivity_doc_specs

REF = D.REPO_ROOT / "tests" / "reference"


def compare_nhgri(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    ref = pd.read_csv(REF / "nhgri_fits.csv")
    refpred = pd.read_csv(REF / "nhgri_predictions.csv")
    specs = nhgri_specs()
    summary, preds, results = fit_many(df, specs)
    rows = []
    for i, (spec, res) in enumerate(zip(specs, results)):
        f = ref.iloc[i]
        old = {p: f[p] for p in PARAM_ORDER if p in f and not pd.isna(f[p])}
        raw = Model(df, spec.with_(log_zero="legacy_raw"))
        old_pred = raw.predict(old)
        row = {"fit": i + 1, "label": spec.label,
               "pred_maxdiff_at_2025_params": float(np.abs(
                   old_pred - refpred.loc[refpred.fit == i + 1, "pred"].to_numpy()).max()),
               "loss_2025": f["objective"],
               "loss_new_code_at_2025_params": raw.loss([old[p] for p in raw.free]),
               "loss_refit": res.loss,
               "loss_refit_at_2025_rule": raw.loss(res.theta),
               "R2_2025": f["R2"], "R2_refit": res.stats["R2_fit"]}
        for p in PARAM_ORDER:
            if p in old:
                row[f"{p}_2025"] = old[p]
                row[f"{p}_refit"] = res.params[p]
        row["max_abs_param_diff"] = max(abs(res.params[p] - old[p]) for p in old)
        rows.append(row)
    return pd.DataFrame(rows), preds


def compare_sensitivity(df: pd.DataFrame) -> pd.DataFrame:
    ref = pd.read_csv(REF / "sensitivity_table.csv")
    ref = ref[ref["method"] == "non-linear optimization"].reset_index(drop=True)
    rows = []
    for spec, (_, f) in zip(sensitivity_doc_specs(), ref.iterrows()):
        res = fit(df, spec)
        old = np.clip([f["h2"], f["m"], f["r"]], 0, 1)
        rows.append({"label": spec.label, "cohort": spec.cohort, "scale": spec.scale,
                     "weights": spec.weights,
                     "h2_2025": f["h2"], "m_2025": f["m"], "r_2025": f["r"],
                     "h2_refit": res.params["V_A"], "m_refit": res.params["m"],
                     "r_refit": res.params["r"],
                     "loss_at_2025_params": res.model.loss(old), "loss_refit": res.loss})
    return pd.DataFrame(rows)


def main(out: Path):
    out.mkdir(parents=True, exist_ok=True)
    df = D.load("MRE8")
    nh, preds = compare_nhgri(df)
    nh.to_csv(out / "nhgri_comparison.csv", index=False)
    preds.to_csv(out / "nhgri_refit_predictions.csv", index=False)
    sens = compare_sensitivity(df)
    sens.to_csv(out / "sensitivity_comparison.csv", index=False)
    with pd.option_context("display.width", 250, "display.float_format", "{:.4g}".format):
        print(nh[["fit", "label", "pred_maxdiff_at_2025_params", "loss_2025",
                  "loss_new_code_at_2025_params", "loss_refit", "max_abs_param_diff"]].to_string(index=False))
        print(sens.to_string(index=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=D.REPO_ROOT / "results" / "legacy")
    main(ap.parse_args().out)
