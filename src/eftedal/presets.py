"""Named collections of `FitSpec`s: the 2025 (legacy) analyses and the corrected defaults."""

from __future__ import annotations

import itertools

from .model import FitSpec, MODELS

# Options that reproduce the 2025 R code (archive/03_model_fitting.R, functions.R).
LEGACY = dict(uz="exponent", log_zero="legacy", adoptive="legacy_n99",
              cocousin_family_degree=3, mr="free", stepsib="eftedal")

# Corrected defaults for new analyses.
CORRECTED = dict(uz="mixture", log_zero="structural", adoptive="zero",
                 cocousin_family_degree=None, stepsib="eftedal")


def nhgri_specs() -> list[FitSpec]:
    """The 32 fits behind the NHGRI 2026 poster (results_ALL.rds), in the same order."""
    out = []
    for name in MODELS:
        for scale, weights in itertools.product(("log", "lin"), ("invSqErr", "hybrid")):
            out.append(FitSpec(components=MODELS[name], cohort="ALL", scale=scale,
                               weights=weights, label=f"{name} {scale},{weights}", **LEGACY))
    return out


def sensitivity_doc_specs() -> list[FitSpec]:
    """'non-linear optimization' rows of the 2025 sensitivity write-up (02_sensitivity_to_weights.R).

    That script used the older additive-only predictor, which did zero out adoptive pairs.
    The 'ALL' cohort was only fitted on the linear ('obs') scale.
    """
    out = []
    for cohort in ("ALL", "ALL_but_adopted", "ALL_but_adopted_UZ", "bio_only"):
        for scale in ("log", "lin"):
            if cohort == "ALL" and scale == "log":
                continue
            for weights in ("invSqErr", "uniform", "hybrid"):
                out.append(FitSpec(components=("A",), cohort=cohort, scale=scale, weights=weights,
                                   uz="exponent", log_zero="legacy", adoptive="zero",
                                   label=f"A {cohort} {scale},{weights}"))
    return out
