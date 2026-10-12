"""Loading the relative-pair correlations and annotating each row.

One row of the loaded table is one category of relative (e.g. maternal same-sex-female
half siblings). Columns after `load()`:

From the authors' files (renamed):
    rel           relationship type (16 values, e.g. 'half sibs', 'cocousins')
    side          'maternal' | 'paternal' | 'other' (parent the pair is related through)
    sex           'ssf' | 'ssm' | 'os' (gender composition of the pair)
    corr, se      observed correlation and its standard error
    npairs        number of pairs

From `data/relationship_types.csv` (one row per `rel`):
    group              'bio' | 'inlaw' | 'adoptive'
    genetic_structure  how the genetic kernel is built (see `components.genetic_kernel`):
                       'full' (full-sib descent), 'full_uz' (through same-sex twins of
                       unknown zygosity), 'half', 'step', 'cocousin', 'none' (adoptive)
    n                  degree of relatedness as defined by Eftedal et al. (for in-laws, the
                       summed distance to the connecting mates; for t/uz classes, the mean
                       of the MZ and DZ alternatives)
    family_degree      our "family degree" f for the decaying family-environment component
    household          'full' | 'partial' | 'by_side' (full if maternal, partial if
                       paternal) | 'none' -- how much the pair shares a household
    cousin_class       1 = first-cousin-like, 2 = more distant cousin-like, 0 = neither
    broad_type         grouping used for the 'hybrid' weights (= rel)

Added here:
    row                 0..81, position in the authors' file
    suspect_duplicate   True for the 2nd cousins/other/ssf row that duplicates
                        2nd cousins/paternal/ssf (see data/README.md)
    uz_linked           True for classes whose relatedness depends on an unknown zygosity
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data"
CORR_DIR = DATA_DIR / "eftedal_osf" / "correlations"

# Order of relationship types as they appear in the authors' files.
REL_ORDER = [
    "mz", "uz", "dz", "sibs", "half sibs", "1st cousins t/uz", "1st cousins",
    "half cousins", "2nd cousins t/uz", "2nd cousins", "step sibs", "cocousins",
    "adoptive sibs", "adoptive 1st cousins", "adoptive 2nd cousins", "adoptive cocousins",
]


def load_relationship_types(path: Path | None = None) -> pd.DataFrame:
    return pd.read_csv(path or DATA_DIR / "relationship_types.csv")


def load(phenotype: str = "MRE8", path: Path | None = None,
         types: pd.DataFrame | None = None) -> pd.DataFrame:
    """Load one phenotype's 82 correlations, annotated with relationship-type attributes."""
    raw = pd.read_csv(path or CORR_DIR / f"{phenotype}.csv")
    df = raw.rename(columns={"type": "side", "pair": "sex", "SE": "se"})[
        ["rel", "side", "sex", "corr", "se", "npairs"]
    ].copy()
    df.insert(0, "row", np.arange(len(df)))

    types = load_relationship_types() if types is None else types
    df = df.merge(types, on="rel", how="left", validate="many_to_one")
    if df["genetic_structure"].isna().any():
        missing = sorted(df.loc[df["genetic_structure"].isna(), "rel"].unique())
        raise ValueError(f"relationship types missing from annotation: {missing}")
    df = df.sort_values("row").reset_index(drop=True)

    df["uz_linked"] = df["genetic_structure"].eq("full_uz")
    df["suspect_duplicate"] = _flag_duplicate(df)
    return df


def _flag_duplicate(df: pd.DataFrame) -> pd.Series:
    """Flag the 2nd-cousin 'other'/'ssf' row when it duplicates the 'paternal'/'ssf' row."""
    flag = pd.Series(False, index=df.index)
    sel = (df["rel"] == "2nd cousins") & (df["sex"] == "ssf")
    other = df.index[sel & (df["side"] == "other")]
    pat = df.index[sel & (df["side"] == "paternal")]
    if len(other) == 1 and len(pat) == 1:
        a, b = df.loc[other[0]], df.loc[pat[0]]
        if a["corr"] == b["corr"] and a["se"] == b["se"]:
            flag[other[0]] = True
    return flag


# --------------------------------------------------------------------------------------
# Cohorts: which rows a model is fitted to.
# --------------------------------------------------------------------------------------

def _all(df):
    return np.ones(len(df), dtype=bool)


def _no_adopted(df):
    return (df["group"] != "adoptive").to_numpy()


def _no_adopted_uz(df):
    return _no_adopted(df) & ~df["uz_linked"].to_numpy()


def _bio_only(df):
    return (df["group"] == "bio").to_numpy() & ~df["uz_linked"].to_numpy()


COHORTS = {
    # every category, including adoptive relatives and relatives through UZ twins
    "ALL": _all,
    "ALL_but_adopted": _no_adopted,
    # Eftedal et al. "Bio + In-law" (also excludes classes linked through UZ twins)
    "ALL_but_adopted_UZ": _no_adopted_uz,
    # Eftedal et al. "Bio only"
    "bio_only": _bio_only,
}


def cohort_mask(df: pd.DataFrame, cohort: str, drop_duplicate: bool = False) -> np.ndarray:
    mask = COHORTS[cohort](df)
    if drop_duplicate:
        mask = mask & ~df["suspect_duplicate"].to_numpy()
    return mask
