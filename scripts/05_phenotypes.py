"""Additive (Fisherian) model across every phenotype in the OSF files, under the main fit
settings, to check that the MRE8 pattern is general. Writes results/phenotypes/summary.csv."""

import argparse
from pathlib import Path

import pandas as pd

from eftedal import data as D
from eftedal.fit import fit
from eftedal.presets import spec

PHENOTYPES = ["MRE8", "MRE5", "MR5", "MR8", "MR9", "math5", "math8", "math9",
              "read5", "read8", "read9", "eng5", "eng8"]


def main(out: Path, n_starts: int):
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for ph in PHENOTYPES:
        df = D.load(ph)
        for cohort in ("ALL_but_adopted_UZ", "bio_only"):
            for setting in ("log-invSqErr", "log-invvar", "lin-invvar", "lin-hybrid"):
                for mr in ("free", "equilibrium"):
                    r = fit(df, spec("M01", cohort, setting, mr=mr), n_starts=n_starts)
                    d = r.summary()
                    d["phenotype"] = ph
                    d["setting"] = setting
                    rows.append(d)
        print(ph, flush=True)
    pd.DataFrame(rows).to_csv(out / "summary.csv", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=D.REPO_ROOT / "results" / "phenotypes")
    ap.add_argument("--n-starts", type=int, default=6)
    a = ap.parse_args()
    main(a.out, a.n_starts)
