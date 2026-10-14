import numpy as np
import pytest

from eftedal import data as D
from eftedal.fit import fit
from eftedal.model import CORE_MODELS, CORE_MODELS_WITH_MATERNAL, FitSpec, Model
from eftedal.presets import MAIN, main_spec


@pytest.fixture(scope="module")
def df():
    return D.load("MRE8")


def test_core_model_set():
    assert len(CORE_MODELS) == 8
    assert not any("maternal" in c for c in CORE_MODELS.values())
    assert len(CORE_MODELS_WITH_MATERNAL) == 16
    assert CORE_MODELS["A"] == ("A",)
    assert CORE_MODELS_WITH_MATERNAL["A+GxG+H+F+mat"] == ("A", "GxG", "house", "family", "maternal")


def test_main_settings():
    s = main_spec("A+H")
    assert (s.scale, s.mr, s.stepsib, s.max_total_variance) == ("lin", "equilibrium", "corrected", 1.0)
    assert s.weights == "hybrid_block"


@pytest.mark.parametrize("pi, old", [(1.0, "house1"), (0.0, "house2")])
def test_household_pi_nests_house1_and_house2(df, pi, old):
    new = Model(df, FitSpec(components=("house",)))
    ref = Model(df, FitSpec(components=(old,)))
    np.testing.assert_allclose(new.predict({"V_house": 0.1, "pi": pi}),
                               ref.predict({f"V_{old}": 0.1}))


def test_household_part_time_rows(df):
    m = Model(df, FitSpec(components=("house",)))
    pred = m.predict({"V_house": 0.1, "pi": 0.3})
    part = ((df["rel"] == "step sibs") | ((df["rel"] == "half sibs") & (df["side"] == "paternal"))).to_numpy()
    full = df["rel"].isin(["mz", "uz", "dz", "sibs", "adoptive sibs"]).to_numpy() | \
        ((df["rel"] == "half sibs") & (df["side"] == "maternal")).to_numpy()
    np.testing.assert_allclose(pred[part], 0.03)
    np.testing.assert_allclose(pred[full], 0.1)
    assert (pred[~part & ~full] == 0).all()


def test_maternal_inflation_only_touches_maternal_rows(df):
    s = FitSpec(components=("A", "family", "maternal"), mr="equilibrium")
    m = Model(df, s)
    p = {"V_A": 0.6, "r": 0.4, "V_family": 0.1, "r_f": 0.5, "delta": 0.2}
    base = Model(df, FitSpec(components=("A", "family"), mr="equilibrium")).predict(m.unpack(
        [p[k] for k in m.free]))
    pred = m.predict(m.unpack([p[k] for k in m.free]))
    mat = (df["side"] == "maternal").to_numpy()
    np.testing.assert_allclose(pred[mat], 1.2 * base[mat])
    np.testing.assert_allclose(pred[~mat], base[~mat])
    comps = m.components(m.unpack([p[k] for k in m.free]))
    np.testing.assert_allclose(comps.sum(axis=1), pred)


def test_delta_bounds_allow_paternal_excess(df):
    m = Model(df, main_spec("A+mat"))
    i = m.free.index("delta")
    assert m.lb[i] < 0 < m.ub[i]


def test_cv_blocks(df):
    assert df["cv_block"].notna().all()
    assert set(df["cv_block"]) == set(D.CV_BLOCKS)
    folds = D.cv_folds(df)
    assert len(folds) == 12
    assert set(folds["twins"]) == {"mz", "uz", "dz"}
    assert set(folds["1st cousins"]) == {"1st cousins", "1st cousins t/uz"}


def test_all_core_models_fit(df):
    for name in CORE_MODELS:
        r = fit(df, main_spec(name), n_starts=1, polish=False)
        assert np.isfinite(r.loss)
        tot = sum(v for k, v in r.params.items() if k.startswith("V_"))
        assert tot <= 1 + 1e-8, name


def test_cross_validation_runs_and_jackknife(df):
    from eftedal import cv
    folds = {k: v for k, v in D.cv_folds(df).items() if k in ("half sibs", "cocousins", "all in-laws")}
    full, fp, held, sc = cv.cross_validate(df, main_spec("A"), n_starts=1, folds=folds)
    assert set(sc["fold"]) == set(folds)
    assert sc.loc[sc.fold == "all in-laws", "group_test"].item()
    # held-out rows are exactly the categories of the held-out types
    assert set(held.loc[held.fold == "cocousins", "rel"]) == {"cocousins"}
    jk = cv.jackknife(fp, full.params, full.free, blocks=["half sibs", "cocousins"])
    assert (jk["se"] >= 0).all()


def test_hybrid_block_weights_sum_to_one_per_block(df):
    m = Model(df, main_spec("A", "hybrid_block"))
    sums = df.assign(w=m.w).groupby("cv_block")["w"].sum()
    np.testing.assert_allclose(sums, 1.0)


def test_sequential_fit_equals_manual_two_stage(df):
    s = main_spec("A+GxG", "invvar")
    seq = fit(df, s.with_(sequential_base=("A",)), n_starts=4)
    a = fit(df, main_spec("A", "invvar"), n_starts=4)
    manual = fit(df, s.with_(fixed=(("V_A", a.params["V_A"]),)), n_starts=4)
    assert seq.params["V_A"] == pytest.approx(a.params["V_A"])
    assert seq.params["V_GxG"] == pytest.approx(manual.params["V_GxG"], abs=1e-6)
    assert seq.loss == pytest.approx(manual.loss, rel=1e-8)
    assert seq.spec.sequential_base == ("A",)        # refits redo both stages
    assert "V_A" not in seq.free


def test_sequential_bootstrap_redoes_stage_one(df):
    from eftedal.inference import parametric_bootstrap
    seq = fit(df, main_spec("A+GxG", "invvar").with_(sequential_base=("A",)), n_starts=2)
    reps = parametric_bootstrap(df, seq, n_boot=5, seed=1, n_starts=1)
    assert reps["V_A"].std() > 0                     # V_A re-estimated in every replicate
