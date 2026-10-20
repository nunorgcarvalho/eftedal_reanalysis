"""Named collections of `FitSpec`s: the 2025 (legacy) analyses and the corrected defaults."""

from __future__ import annotations

import itertools

from .model import CORE_MODELS, JOINT_MODELS, PHASE2_MODELS, SEQUENTIAL_CARRY, FitSpec, MODELS

# Options that reproduce the 2025 R code (archive/03_model_fitting.R, functions.R).
LEGACY = dict(uz="exponent", log_zero="legacy", adoptive="legacy_n99",
              cocousin_family_degree=3, mr="free", stepsib="eftedal")

# Main analysis settings (agreed plan, 2026-10-12): linear scale, m = V_A r, corrected
# step-sibling kernel, co-cousin family degree 2 (annotation table), sum of V <= 1, UZ as an
# MZ/DZ mixture at 0.5. Weighting is chosen per fit: 'hybrid_block' (primary; every cv block
# of relationship types gets equal weight -- user decision 2026-10-13, replacing the 16-type
# 'hybrid') or 'invvar' (1/SE^2; the model taken literally).
MAIN = dict(scale="lin", mr="equilibrium", stepsib="corrected", uz="mixture", mz_frac=0.5,
            log_zero="structural", adoptive="zero", cocousin_family_degree=None,
            max_total_variance=1.0, gxg_kernel="a2")   # GxG kernel a^2 (user decision 2026-10-14)
MAIN_WEIGHTINGS = ("hybrid_block", "invvar")


def main_spec(model: str, weights: str = "hybrid_block", cohort: str = "ALL", **kw) -> FitSpec:
    """Spec under the main settings for a model of `CORE_MODELS` (incl. the two-step
    'A>GxG...' models), `JOINT_MODELS`, or the 2025 set `MODELS`."""
    extra = {}
    if model in CORE_MODELS:
        comps, base = CORE_MODELS[model]
        if base is not None:
            extra = dict(sequential_base=base, sequential_carry=SEQUENTIAL_CARRY)
    elif model in JOINT_MODELS:
        comps = JOINT_MODELS[model]
    elif model in PHASE2_MODELS:
        comps, mods = PHASE2_MODELS[model]
        extra = dict(modifiers=mods)
    else:
        comps = MODELS[model]
    return FitSpec(components=comps, cohort=cohort, weights=weights, **{**MAIN, **extra, **kw},
                   label=f"{model}|{cohort}|{weights}")


# Corrected defaults used by the exploratory 2026-10-12 grids.
CORRECTED = dict(uz="mixture", log_zero="structural", adoptive="zero",
                 cocousin_family_degree=None, stepsib="eftedal", max_total_variance=1.0)


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


# --------------------------------------------------------------------------------------
# Grids for the corrected analyses
# --------------------------------------------------------------------------------------

# (scale, weights) combinations and their short names
FIT_SETTINGS = {
    "log-invSqErr": ("log", "invSqErr"),   # Eftedal et al.
    "log-invvar": ("log", "invvar"),       # log scale, delta-method weights (corrected)
    "lin-invvar": ("lin", "invvar"),       # linear scale, 1/SE^2 (= Gaussian ML)
    "log-hybrid": ("log", "hybrid"),
    "lin-hybrid": ("lin", "hybrid"),       # NHGRI poster setting
    "lin-uniform": ("lin", "uniform"),
}


def fisher_grid() -> list[FitSpec]:
    """Additive-only (Fisherian) model: cohort x fit setting x m-r mode x step-sib formula."""
    out = []
    for cohort, setting, mr, stepsib in itertools.product(
            ("ALL_but_adopted_UZ", "ALL_but_adopted", "bio_only"), FIT_SETTINGS,
            ("free", "equilibrium"), ("eftedal", "corrected")):
        if cohort == "bio_only" and stepsib == "corrected":
            continue   # no step siblings in the bio-only cohort
        scale, weights = FIT_SETTINGS[setting]
        out.append(FitSpec(components=("A",), cohort=cohort, scale=scale, weights=weights,
                           mr=mr, stepsib=stepsib, **{k: v for k, v in CORRECTED.items() if k != "stepsib"},
                           label=f"M01|{cohort}|{setting}|{mr}|{stepsib}"))
    return out


def model_grid(settings=("log-invSqErr", "log-invvar", "lin-invvar", "lin-hybrid")) -> list[FitSpec]:
    """Every named model on all 82 categories: fit setting x m-r mode x step-sib formula."""
    out = []
    for name, setting, mr, stepsib in itertools.product(MODELS, settings, ("free", "equilibrium"),
                                                        ("eftedal", "corrected")):
        scale, weights = FIT_SETTINGS[setting]
        out.append(FitSpec(components=MODELS[name], cohort="ALL", scale=scale, weights=weights,
                           mr=mr, stepsib=stepsib, **{k: v for k, v in CORRECTED.items() if k != "stepsib"},
                           label=f"{name}|ALL|{setting}|{mr}|{stepsib}"))
    return out


def spec(model: str, cohort: str, setting: str, mr: str = "free", stepsib: str = "eftedal",
         **kw) -> FitSpec:
    """Shorthand for a corrected-defaults spec of a named model and fit setting."""
    scale, weights = FIT_SETTINGS[setting]
    base = {k: v for k, v in CORRECTED.items() if k != "stepsib"}
    return FitSpec(components=MODELS[model], cohort=cohort, scale=scale, weights=weights,
                   mr=mr, stepsib=stepsib, **{**base, **kw},
                   label=f"{model}|{cohort}|{setting}|{mr}|{stepsib}")


def key_specs() -> list[FitSpec]:
    """Specs probed by leave-one-out, one-at-a-time checks and uncertainty analyses."""
    return [
        spec("M01", "ALL_but_adopted_UZ", "log-invSqErr"),                    # Eftedal et al.
        spec("M01", "ALL_but_adopted_UZ", "log-invvar"),
        spec("M01", "ALL_but_adopted_UZ", "lin-invvar"),
        spec("M01", "ALL_but_adopted_UZ", "log-invvar", mr="equilibrium", stepsib="corrected"),
        spec("M01", "ALL", "lin-hybrid"),                                      # NHGRI fig 1
        spec("M02", "ALL", "lin-hybrid"),                                      # NHGRI fig 2
        spec("M04b", "ALL", "lin-hybrid"),                                     # NHGRI fig 3
        spec("M02", "ALL", "lin-invvar", mr="equilibrium", stepsib="corrected"),
        spec("M04b", "ALL", "lin-invvar", mr="equilibrium", stepsib="corrected"),
        spec("M05b", "ALL", "lin-invvar", mr="equilibrium", stepsib="corrected"),
    ]
