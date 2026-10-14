import numpy as np
import pandas as pd
import pytest

from eftedal import components as C
from eftedal import data as D
from eftedal import eftedal2025 as E
from eftedal.fit import fit
from eftedal.model import FitSpec, Model, PARAM_ORDER
from eftedal.presets import nhgri_specs

REF = D.REPO_ROOT / "tests" / "reference"


@pytest.fixture(scope="module")
def df():
    return D.load("MRE8")


# ---- data -------------------------------------------------------------------------------

def test_load(df):
    assert len(df) == 82
    assert df["rel"].nunique() == 16
    assert df["genetic_structure"].notna().all()
    assert df["suspect_duplicate"].sum() == 1
    dup = df[df["suspect_duplicate"]].iloc[0]
    assert (dup["rel"], dup["side"], dup["sex"]) == ("2nd cousins", "other", "ssf")


def test_cohort_sizes(df):
    sizes = {c: D.cohort_mask(df, c).sum() for c in D.COHORTS}
    assert sizes == {"ALL": 82, "ALL_but_adopted": 70, "ALL_but_adopted_UZ": 53, "bio": 58, "bio_only": 41}


# ---- kernels ------------------------------------------------------------------------------

def test_stepsib_kernels_match_writeup(df):
    """relative_covariance v1 / eftedal note: at h2=.503, m=.430 their step-sib prediction is
    0.111; the corrected one is 0.186 at r = m/h2 = 0.855 and 0.111 at r = m."""
    step = (df["rel"] == "step sibs").to_numpy()
    h2, m = 0.503, 0.430
    K, _ = C.genetic_kernels(df, m, 0.855, stepsib="eftedal")
    assert h2 * K[step][0] == pytest.approx(0.111, abs=5e-4)
    K, _ = C.genetic_kernels(df, m, m / h2, stepsib="corrected")
    assert h2 * K[step][0] == pytest.approx(0.186, abs=5e-4)
    K, _ = C.genetic_kernels(df, m, m, stepsib="corrected")
    assert h2 * K[step][0] == pytest.approx(0.111, abs=5e-4)


def test_uz_mixture_matches_eftedal_prediction(df):
    m = 0.43
    K, K2 = C.genetic_kernels(df, m, 0.3, uz="mixture", mz_frac=0.5)
    uz = (df["rel"] == "uz").to_numpy()
    assert K[uz][0] == pytest.approx(0.5 + 0.5 * (1 + m) / 2)
    assert K2[uz][0] == pytest.approx(0.5 + 0.5 * ((1 + m) / 2) ** 2)


def test_household_indicators(df):
    full, part = C.household_share(df)
    hs = df["rel"] == "half sibs"
    assert (full[hs & (df["side"] == "maternal")] == 1).all()
    assert (part[hs & (df["side"] == "paternal")] == 1).all()
    assert (part[df["rel"] == "step sibs"] == 1).all()
    assert full[df["rel"] == "adoptive sibs"].sum() == 3


# ---- model ---------------------------------------------------------------------------------

def test_log_invvar_weights_are_delta_method(df):
    m = Model(df, FitSpec(scale="log", weights="invvar", cohort="ALL_but_adopted"))
    sel = m.cohort & (df["corr"].to_numpy() > 0)
    np.testing.assert_allclose(m.w[sel], (df["corr"] ** 2 / df["se"] ** 2).to_numpy()[sel])


def test_equilibrium_mode(df):
    m = Model(df, FitSpec(components=("A",), mr="equilibrium"))
    assert m.free == ["V_A", "r"]
    p = m.unpack([0.6, 0.5])
    assert p["m"] == pytest.approx(0.3)


def test_fixed_params(df):
    m = Model(df, FitSpec(components=("A",), fixed=(("r", 0.47),)))
    assert m.free == ["V_A", "m"]
    assert m.unpack([0.5, 0.3])["r"] == 0.47


# ---- legacy reproduction -----------------------------------------------------------------

def test_legacy_predictions_and_losses(df):
    """At the archived 2025 parameters, predictions and losses match the 2025 R code."""
    ref = pd.read_csv(REF / "nhgri_fits.csv")
    refpred = pd.read_csv(REF / "nhgri_predictions.csv")
    for i, spec in enumerate(nhgri_specs()):
        f = ref.iloc[i]
        old = {p: f[p] for p in PARAM_ORDER if p in f and not pd.isna(f[p])}
        m = Model(df, spec.with_(log_zero="legacy_raw"))
        np.testing.assert_allclose(m.predict(old), refpred.loc[refpred.fit == i + 1, "pred"],
                                   atol=1e-9, err_msg=spec.label)
        assert m.loss([old[p] for p in m.free]) == pytest.approx(f["objective"], rel=1e-6)  # params stored to ~15 digits


@pytest.mark.parametrize("i", [0, 3, 7, 23, 31])
def test_legacy_refit(df, i):
    """Fresh fits land on the 2025 optimum (fits where the 2025 optimiser did not stop early
    and did not exploit the log-scale row-dropping rule)."""
    ref = pd.read_csv(REF / "nhgri_fits.csv").iloc[i]
    spec = nhgri_specs()[i]
    res = fit(df, spec, n_starts=6)
    assert res.loss <= ref["objective"] * (1 + 1e-6)
    for p in res.free:
        assert res.params[p] == pytest.approx(ref[p], abs=2e-3), p


# ---- Eftedal et al. replication ---------------------------------------------------------

def test_eftedal_bio_only(df):
    b = E.fit_bio_only(df)
    assert (b.h2, b.m, b.c) == pytest.approx((0.801, 0.225, 1.002), abs=1e-3)
    assert b.R2_log == pytest.approx(0.988, abs=1e-3)
    ci = b.ci.set_index("param")
    assert ci.loc["h2", ["lo", "hi"]].tolist() == pytest.approx([0.749, 0.856], abs=1e-3)
    assert ci.loc["m", ["lo", "hi"]].tolist() == pytest.approx([0.204, 0.246], abs=1e-3)


def test_eftedal_bio_inlaw(df):
    a = E.fit_bio_inlaw(df)
    # published 0.503 / 0.430 / 0.972 come from a slightly unconverged optimiser
    assert (a.h2, a.m, a.c) == pytest.approx((0.503, 0.430, 0.972), abs=5e-3)
    pred = E.predict(df[D.cohort_mask(df, "ALL_but_adopted_UZ")], 0.503, 0.430, 0.972)
    d = df[D.cohort_mask(df, "ALL_but_adopted_UZ")]
    assert E.weighted_r2(d["corr"].to_numpy(), pred, 1 / d["se"].to_numpy() ** 2) == pytest.approx(0.873, abs=1e-3)


def test_eftedal_matches_generic_model(df):
    """The generic model with Eftedal's settings reproduces their dedicated fit."""
    a = E.fit_bio_inlaw(df)
    res = fit(df, FitSpec(components=("A",), cohort="ALL_but_adopted_UZ", scale="log",
                          weights="invSqErr"), n_starts=4)
    assert res.params["V_A"] == pytest.approx(a.h2, abs=2e-3)
    assert res.params["m"] == pytest.approx(a.m, abs=2e-3)


def test_twin_ace(df):
    t = E.twin_ace(df)
    assert t["A"] == pytest.approx(0.75, abs=0.01)
    assert t["C"] == pytest.approx(0.11, abs=0.01)
