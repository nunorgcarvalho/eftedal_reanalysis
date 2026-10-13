import numpy as np
import pytest

from eftedal import data as D
from eftedal.fit import fit
from eftedal.inference import covariance, parametric_bootstrap
from eftedal.model import FitSpec, Model
from eftedal.presets import spec


@pytest.fixture(scope="module")
def df():
    return D.load("MRE8")


@pytest.fixture(scope="module")
def lin_fit(df):
    return fit(df, spec("M01", "ALL_but_adopted_UZ", "lin-invvar", mr="equilibrium"), n_starts=4)


def test_exclude_types(df):
    m = Model(df, FitSpec(cohort="ALL_but_adopted_UZ", exclude_types=("cocousins",)))
    assert m.cohort.sum() == 53 - 9
    assert not m.cohort[(df["rel"] == "cocousins").to_numpy()].any()


def test_hybrid_weights_sum_to_one_per_type(df):
    m = Model(df, FitSpec(cohort="ALL", weights="hybrid"))
    sums = df.assign(w=m.w).groupby("broad_type")["w"].sum()
    np.testing.assert_allclose(sums, 1.0)


def test_corrected_stepsib_equals_eftedal_when_r_equals_m(df):
    a = fit(df, spec("M01", "ALL_but_adopted_UZ", "lin-invvar", mr="equal", stepsib="eftedal"), n_starts=2)
    b = fit(df, spec("M01", "ALL_but_adopted_UZ", "lin-invvar", mr="equal", stepsib="corrected"), n_starts=2)
    assert a.loss == pytest.approx(b.loss, rel=1e-9)


def test_hessian_se_depends_on_weight_scale_sandwich_does_not(lin_fit):
    m = lin_fit.model
    w0 = m.w.copy()
    se = {}
    for k in (1.0, 1e-4):
        m.w = w0 * k
        se[k] = {meth: np.sqrt(np.diag(covariance(lin_fit, meth))) for meth in ("hessian_rss", "sandwich")}
    m.w = w0
    np.testing.assert_allclose(se[1e-4]["hessian_rss"], se[1.0]["hessian_rss"] * 100, rtol=1e-3)
    np.testing.assert_allclose(se[1e-4]["sandwich"], se[1.0]["sandwich"], rtol=1e-6)


def test_sandwich_matches_hessian_for_inverse_variance_weights(lin_fit):
    """With 1/SE^2 weights on the linear scale, the Gauss-Newton covariance is ~2 H^-1 of the RSS
    (equal up to the residual-curvature term of the full Hessian, a few % here because the
    residuals are large)."""
    s = np.sqrt(np.diag(covariance(lin_fit, "sandwich")))
    h = np.sqrt(np.diag(2 * covariance(lin_fit, "hessian_rss")))
    np.testing.assert_allclose(s, h, rtol=0.1)


def test_sandwich_agrees_with_parametric_bootstrap(df, lin_fit):
    reps = parametric_bootstrap(df, lin_fit, n_boot=150, seed=3)
    se_boot = reps[lin_fit.free].std().to_numpy()
    se_sand = np.sqrt(np.diag(covariance(lin_fit, "sandwich")))
    np.testing.assert_allclose(se_boot, se_sand, rtol=0.25)


def test_equilibrium_relation_holds_after_fit(lin_fit):
    p = lin_fit.params
    assert p["m"] == pytest.approx(p["V_A"] * p["r"])
