"""Component correlation decay across relationship types (main figure M4): for given parameters
(no fitting), each category's contribution of the additive (random-mating part and AM excess), GxG
and familial-environment components, scaled so a pair whose kernels are all 1 has correlation 1.
Writes results/decay/{name}.csv (one per example; plotted by R/make_figures.R).

    python scripts/07_decay.py                                     # the built-in examples
    python scripts/07_decay.py --name gen2 --model A+GxG --params V_A=0.4,m=0.2,V_GxG=0.3
"""

import argparse
from pathlib import Path

from eftedal import data as D
from eftedal.fit import component_decay
from eftedal.model import FitSpec
from eftedal.presets import MAIN

ABBR = {"A": "A", "GxG": "GxG", "H": "house", "F": "family"}

# name -> (components, as in model names: A, GxG, H, F; parameters)
EXAMPLES = {"genetics": ("A+GxG", dict(V_A=0.5, m=0.25, V_GxG=0.25)),
            "environment": ("F", dict(V_family=1.0, r_f=0.5))}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--name")
    ap.add_argument("--model", help="components, e.g. A+GxG or F (main settings)")
    ap.add_argument("--params", help="comma-separated key=value, e.g. V_A=0.5,m=0.25,V_GxG=0.25")
    ap.add_argument("--out", default=str(D.REPO_ROOT / "results" / "decay"))
    a = ap.parse_args()
    ex = EXAMPLES
    if a.name:
        ex = {a.name: (a.model, {k: float(v) for k, v in (kv.split("=") for kv in a.params.split(","))})}
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    df = D.load()
    for name, (model, params) in ex.items():
        comps = tuple(ABBR[c] for c in model.split("+"))
        spec = FitSpec(components=comps, **{**MAIN, **({} if "A" in comps else {"mr": "free"})})
        component_decay(df, spec, params).assign(model=model).to_csv(out / f"{name}.csv", index=False)
        print("wrote", out / f"{name}.csv")


if __name__ == "__main__":
    main()
