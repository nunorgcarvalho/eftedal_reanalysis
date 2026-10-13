"""Fit the 16 core models (`model.CORE_MODELS`) under the main settings (`presets.MAIN`) with
both weightings, on all 82 categories. Writes results/main/:

fits_summary.csv       one row per model x weighting: parameters, loss, R2_lin, chi2, ...
fits_predictions.csv   predictions and per-component contributions for every category
"""

import argparse
from pathlib import Path

from eftedal import data as D
from eftedal.fit import fit_many
from eftedal.model import CORE_MODELS
from eftedal.presets import MAIN_WEIGHTINGS, main_spec


def main(out: Path, n_starts: int):
    out.mkdir(parents=True, exist_ok=True)
    df = D.load("MRE8")
    specs = [main_spec(m, w) for w in MAIN_WEIGHTINGS for m in CORE_MODELS]
    summary, preds, _ = fit_many(df, specs, n_starts=n_starts)
    summary.insert(1, "model", [s.label.split("|")[0] for s in specs])
    summary.to_csv(out / "fits_summary.csv", index=False)
    preds.to_csv(out / "fits_predictions.csv", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=D.REPO_ROOT / "results" / "main")
    ap.add_argument("--n-starts", type=int, default=20)
    a = ap.parse_args()
    main(a.out, a.n_starts)
