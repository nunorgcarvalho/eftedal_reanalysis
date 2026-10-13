"""Joint vs sequential estimation (base model A first, its V_A -- or V_A and r -- held fixed,
then the extra components estimated on top), for the models with GxG and/or the family
environment, under both weightings. Writes results/sequential/:

estimates.csv   joint and sequential fits side by side
intervals.csv   parametric-bootstrap intervals of the sequential estimator (both stages redone in
                every replicate) next to those of the joint fit, from results/uncertainty_main
"""

import argparse
from pathlib import Path

import pandas as pd

from eftedal import data as D
from eftedal.fit import fit
from eftedal.inference import parametric_bootstrap, percentile_table
from eftedal.presets import MAIN_WEIGHTINGS, main_spec

MODELS = ["A+GxG", "A+F", "A+GxG+F"]
VARIANTS = {"joint": None, "sequential (fix V_A)": ("V_A",), "sequential (fix V_A, r)": ("V_A", "r")}
REPORT = ["V_A", "V_GxG", "V_family", "r_f", "m", "r"]


def main(out: Path, n_boot: int, n_starts: int):
    out.mkdir(parents=True, exist_ok=True)
    df = D.load("MRE8")
    rows, ints = [], []
    for k, (w, m) in enumerate([(w, m) for w in MAIN_WEIGHTINGS for m in MODELS]):
        joint = None
        for name, carry in VARIANTS.items():
            s = main_spec(m, w)
            if carry is not None:
                s = s.with_(sequential_base=("A",), sequential_carry=carry)
            res = fit(df, s, n_starts=n_starts)
            joint = res if carry is None else joint
            d = res.summary()
            d.update(model=m, weighting=w, procedure=name, loss_vs_joint=res.loss / joint.loss)
            rows.append(d)
            if name == "sequential (fix V_A)":
                b = parametric_bootstrap(df, res, n_boot=n_boot, seed=500 + k, n_starts=3)
                est = {p: res.params[p] for p in REPORT if p in res.params}
                ints.append(percentile_table(b, est, method="param_boot x1").assign(
                    model=m, weights=w, procedure=name))
        print(w, m, flush=True)
    pd.DataFrame(rows).to_csv(out / "estimates.csv", index=False)
    main_iv = D.REPO_ROOT / "results" / "uncertainty_main" / "intervals.csv"
    iv = pd.concat(ints, ignore_index=True)
    if main_iv.exists():
        j = pd.read_csv(main_iv)
        j = j[(j["method"] == "param_boot x1") & j["model"].isin(MODELS) & j["param"].isin(REPORT)]
        iv = pd.concat([iv, j.assign(procedure="joint")[iv.columns.intersection(j.columns).tolist() + ["procedure"]]],
                       ignore_index=True)
    iv.to_csv(out / "intervals.csv", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=D.REPO_ROOT / "results" / "sequential")
    ap.add_argument("--n-boot", type=int, default=300)
    ap.add_argument("--n-starts", type=int, default=12)
    a = ap.parse_args()
    main(a.out, a.n_boot, a.n_starts)
