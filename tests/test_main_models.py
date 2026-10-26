import numpy as np
import pytest

from eftedal import data as D
from eftedal.fit import fit
from eftedal.model import CORE_MODELS, JOINT_MODELS, FitSpec, Model
from eftedal.presets import MAIN, main_spec


@pytest.fixture(scope="module")
def df():
    return D.load("MRE8")


def test_model_sets():
    assert len(JOINT_MODELS) == 8 and JOINT_MODELS["A+GxG+H+F"] == ("A", "GxG", "house", "family")
    assert list(CORE_MODELS) == ["A", "A>GxG", "A+GxG", "A>F", "A+F", "A+GxG>F", "A+F>GxG", "A>GxG+F", "A+GxG+F"]
    s = main_spec("A>GxG+F")
    assert s.sequential_base == ("A",) and s.sequential_carry == ("all",)
    assert main_spec("A+F").sequential_base is None


def test_two_step_freezes_the_additive_part(df):
    a = fit(df, main_spec("A", "invvar"), n_starts=2)
    s = fit(df, main_spec("A>GxG+F", "invvar"), n_starts=2)  # a2 kernel (MAIN)
    for p in ("V_A", "r", "m"):
        assert s.params[p] == pytest.approx(a.params[p])
    # Phase 2: rho_g is estimated with A in stage 1 and carried; the familial modifiers enter with F
    assert s.params["rho_g"] == pytest.approx(a.params["rho_g"])
    assert set(s.free) == {"V_GxG", "V_family", "r_f", "rho_e", "k_side"}
    assert s.summary()["V_sum"] == pytest.approx(s.params["V_A"] + s.params["V_GxG"] + s.params["V_family"])
    # step 1 = A + F: V_A, V_family, r_f and r frozen; only V_GxG is fitted in step 2
    af = fit(df, main_spec("A+F", "invvar"), n_starts=2)
    t = fit(df, main_spec("A+F>GxG", "invvar"), n_starts=2)
    assert t.free == ["V_GxG"]
    for p in ("V_A", "V_family", "r_f", "r", "rho_g", "rho_e", "k_side"):
        assert t.params[p] == pytest.approx(af.params[p])


def test_jackknife_folds_restricted_to_cohort(df):
    from eftedal import cv
    full, fp, _, _ = cv.cross_validate(df, main_spec("A", "invvar", cohort="bio"), n_starts=1)
    assert set(fp["fold"]) == {"twins", "sibs", "half sibs", "1st cousins", "half cousins", "2nd cousins"}
    jk = cv.jackknife(fp, full.params, ["V_A", "V_sum"])
    assert jk.set_index("param").loc["V_sum", "estimate"] == pytest.approx(full.params["V_A"])


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
    m = Model(df, FitSpec(components=("A", "maternal")))
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
    for name in list(CORE_MODELS) + list(JOINT_MODELS):
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


def test_a2_gxg_kernel(df):
    from eftedal import components as C
    kd = C.KernelData(df)
    _, k_am = C.genetic_kernels(kd, 0.3, 0.4, stepsib="corrected")
    _, k_a2 = C.genetic_kernels(kd, 0.3, 0.4, stepsib="corrected", gxg="a2")
    _, k_a2b = C.genetic_kernels(kd, 0.0, 0.0, stepsib="corrected", gxg="a2")
    np.testing.assert_allclose(k_a2, k_a2b)                    # independent of m, r
    get = lambda k, rel: k[df.index[df.rel == rel][0]]
    assert get(k_a2, "mz") == 1 and get(k_a2, "sibs") == pytest.approx(0.25)
    assert get(k_a2, "half sibs") == pytest.approx(0.25 ** 2)
    assert get(k_a2, "1st cousins") == pytest.approx(0.25 ** 3)
    assert get(k_a2, "uz") == pytest.approx(0.5 * 1 + 0.5 * 0.25)
    assert get(k_a2, "step sibs") == 0 and get(k_a2, "cocousins") == 0
    assert get(k_am, "step sibs") > 0                           # the AM version is not


def test_random_mating_VA(df):
    r = fit(df, main_spec("A", "invvar"), n_starts=2)
    assert r.summary()["V_A0"] == pytest.approx(r.params["V_A"] * (1 - r.params["m"]))
    from eftedal import cv
    full, fp, _, _ = cv.cross_validate(df, main_spec("A", "invvar"), n_starts=1,
                                       folds={b: D.rels_in_blocks(df, b) for b in ("sibs", "twins")})
    np.testing.assert_allclose(fp["V_A0"], fp["V_A"] * (1 - fp["m"]))
    jk = cv.jackknife(fp, full.params, ["V_A0"], blocks=["sibs", "twins"]).iloc[0]
    assert jk["estimate"] == pytest.approx(r.params["V_A"] * (1 - r.params["m"]), rel=1e-4)


def test_additive_split_rm_and_am(df):
    r = fit(df, main_spec("A+GxG+F", "invvar"), n_starts=2)
    p = r.predictions()
    np.testing.assert_allclose(p.part_A0 + p.part_A_AM, p.comp_A)
    assert (p.part_A_AM >= -1e-12).all()                      # AM only adds sharing
    inlaw = p.rel.isin(["step sibs", "cocousins"])
    np.testing.assert_allclose(p.loc[inlaw, "part_A0"], 0)    # in-laws share nothing by descent
    sib = p.index[p.rel == "sibs"][0]
    assert p.part_A0[sib] == pytest.approx(r.summary()["V_A0"] * 0.5)


def test_component_decay_scaled_to_kernel_one(df):
    from eftedal.fit import component_decay
    g = component_decay(df, FitSpec(components=("A", "GxG"), **MAIN), dict(V_A=0.5, m=0.25, V_GxG=0.25))
    assert g.r.iloc[0] == pytest.approx(0.5)                  # r = m / V_A at equilibrium
    np.testing.assert_allclose(g.part_A0 + g.part_A_AM + g.comp_GxG, g.total)
    mz, sib = (g.rel == "mz").idxmax(), (g.rel == "sibs").idxmax()
    assert g.total[mz] == pytest.approx(1.0)
    assert g.part_A0[mz] == pytest.approx(0.375 / 0.75)          # V_A0 / (V_A + V_GxG)
    assert g.total[sib] == pytest.approx((0.5 * 1.25 / 2 + 0.25 / 4) / 0.75)
    e = component_decay(df, FitSpec(components=("family",), **{**MAIN, "mr": "free"}), dict(V_family=1.0, r_f=0.5))
    np.testing.assert_allclose(e.total, 0.5 ** df.family_degree.to_numpy())


def test_family_split(df):
    from eftedal.fit import variance_sum
    r = fit(df, main_spec("A+GxG+F", "invvar"), n_starts=2)
    s, p = r.summary(), r.predictions()
    assert s["V_fam_sib"] + s["V_fam_ext"] == pytest.approx(r.params["V_family"])
    assert s["V_sum"] == pytest.approx(variance_sum(r.params))   # derived parts not double-counted
    np.testing.assert_allclose(p.part_fam_sib + p.part_fam_ext, p.comp_family)
    f = Model(df, main_spec("A+GxG+F")).f
    np.testing.assert_allclose(p.part_fam_sib[f >= 1], 0)     # only pairs closer than degree 1
    sib = ((p.rel == "sibs") & (p.sex != "os")).idxmax()     # multiplier 1: same sex, side 'other'
    assert p.part_fam_sib[sib] == pytest.approx(s["V_fam_sib"])
    assert p.part_fam_ext[sib] == pytest.approx(s["V_fam_ext"])
    # r_f = 0: everything sib-specific
    e = Model(df, main_spec("A+F", modifiers=()))
    from eftedal.fit import family_split
    sib0, ext0 = family_split(e, {"V_family": 0.2, "r_f": 0.0})
    np.testing.assert_allclose(ext0, 0); np.testing.assert_allclose(sib0, 0.2 * (e.f == 0))


def test_phase2_reduces_to_phase1(df):
    """Phase-2 modifiers at their null values reproduce Phase 1 exactly."""
    from eftedal.model import PHASE2_MODELS
    p1 = fit(df, main_spec("A+GxG+F", "invvar", modifiers=()), n_starts=2)
    full = main_spec("A+GxG+F", "invvar")
    assert set(full.modifiers) == {"sex_g", "sex_e", "side"}
    null = dict(p1.params, rho_g=1.0, rho_e=1.0, k_side=0.0)
    np.testing.assert_allclose(Model(df, full).predict(null), p1.pred)
    # fixing the multipliers at 1 gives the Phase-1 fit
    r = fit(df, full.with_(fixed=(("rho_g", 1.0), ("rho_e", 1.0), ("k_side", 0.0))), n_starts=2)
    assert r.loss == pytest.approx(p1.loss, rel=1e-6)
    assert all(name.startswith("A+GxG+F | ") for name in PHASE2_MODELS)
    assert main_spec("A+GxG+F | none").modifiers == ()
    # Phase 2 is the default for the core set; the Additive model of the corrections stays Phase 1
    assert main_spec("A").modifiers == ("sex_g",) and main_spec("A", modifiers=()).modifiers == ()


def test_phase2_multipliers_act_where_they_should(df):
    spec = main_spec("A+GxG+F", "invvar")
    m = Model(df, spec)
    base = dict(V_A=0.4, V_GxG=0.2, V_family=0.15, r_f=0.5, r=0.5, m=0.2, rho_g=1.0, rho_e=1.0, k_side=0.0)
    c0 = m.components(base)
    c = m.components(dict(base, rho_g=0.9, rho_e=0.5, k_side=0.2))
    os_, mat, pat = m.os, m.side_mat, m.side_pat
    np.testing.assert_allclose(c.A[os_], 0.9 * c0.A[os_]); np.testing.assert_allclose(c.A[~os_], c0.A[~os_])
    np.testing.assert_allclose(c.GxG[os_], 0.81 * c0.GxG[os_])
    side = np.where(mat, 1.2, np.where(pat, 0.8, 1.0)) * np.where(os_, 0.5, 1.0)
    np.testing.assert_allclose(c.family, side * c0.family)
    np.testing.assert_allclose(m.predict(dict(base, rho_g=0.9, rho_e=0.5, k_side=0.2)), c.sum(axis=1))
    with pytest.raises(ValueError):
        FitSpec(components=("A",), modifiers=("side",), **{k: v for k, v in MAIN.items()})


def test_x_kernel(df):
    """The X grid (theory/x_kernel) covers every category, matches our autosomal kernel, and
    reproduces the random-mating X kernels of the report."""
    import pandas as pd
    from eftedal import components as C
    grid = C.x_kernel_grid(df)
    g = pd.read_csv(C.X_GRID_FILE)
    m = Model(df, main_spec("A", "invvar"))
    for mm in (0.0, 0.2, 0.4):
        K, _ = C.genetic_kernels(m.kd, mm, mm, stepsib="corrected", uz="mixture", mz_frac=0.5)
        theirs = g[g.m == mm].set_index(["rel", "side", "sex"]).Kauto
        keys = list(df[["rel", "side", "sex"]].itertuples(index=False, name=None))
        # kernels that also depend on r (Nagylaki half-sib types, step sibs, co-cousins) are compared
        # at m = 0 only: the grid ties r to its own h^2; D varies only weakly with it
        r_dep = ["half sibs", "1st cousins t/uz", "2nd cousins t/uz", "half cousins", "step sibs", "cocousins"]
        keep = np.ones(len(df), bool) if mm == 0 else ~df.rel.isin(r_dep).to_numpy()
        np.testing.assert_allclose(K[keep], theirs.loc[keys].to_numpy()[keep], atol=2e-3)
    D0 = C.x_kernel_adjustment(grid, 0.0)
    pick = lambda rel, side, sex: D0[((df.rel == rel) & (df.side == side) & (df.sex == sex)).to_numpy()][0]
    assert pick("sibs", "other", "os") == pytest.approx(1 / 3 - 1 / 2, abs=1e-5)       # K^X_FDC - a
    assert pick("half sibs", "paternal", "ssm") == pytest.approx(0 - 1 / 4, abs=1e-5)
    assert pick("1st cousins", "maternal", "ssm") == pytest.approx(0.5 - 1 / 8, abs=1e-5)
    # x_share = 0 is Phase 1; x_share > 0 moves only the additive part
    p = dict(V_A=0.5, r=0.4, m=0.2)
    a = Model(df, main_spec("A", "invvar", modifiers=())).predict(p)
    b = Model(df, main_spec("A", "invvar", modifiers=(), x_share=0.02)).predict(p)
    np.testing.assert_allclose(b - a, 0.5 * 0.02 * C.x_kernel_adjustment(grid, 0.2), atol=1e-12)


def test_phase2_splits_add_up(df):
    """With the multipliers away from 1, the additive and familial splits still add up."""
    r = fit(df, main_spec("A+GxG+F", "invvar"), n_starts=2)
    p = r.predictions()
    np.testing.assert_allclose(p.part_A0 + p.part_A_AM, p.comp_A)
    np.testing.assert_allclose(p.part_fam_sib + p.part_fam_ext, p.comp_family)
    np.testing.assert_allclose(p.comp_A + p.comp_GxG + p.comp_family, p.pred)


def test_scripts_compile():
    """Every pipeline script at least parses (they are not imported by the other tests)."""
    import py_compile
    from pathlib import Path
    root = Path(__file__).resolve().parents[1]
    for f in sorted((root / "scripts").glob("*.py")) + [root / "ASHG" / "QnA" / "make_tables.py"]:
        py_compile.compile(str(f), doraise=True)


def test_axc_kernel(df):
    """A x C contributes V_AxC K_A(m) c^f (with the Phase-2 multipliers), zero for adoptees."""
    from eftedal import components as C
    spec = main_spec("A+F+AxC", "invvar")
    assert spec.modifiers == ("sex_g", "sex_e", "side")
    m = Model(df, spec)
    p = dict(V_A=0.4, V_family=0.1, V_AxC=0.1, r_f=0.5, r=0.5, m=0.2, rho_g=1.0, rho_e=1.0, k_side=0.0)
    K, _ = C.genetic_kernels(m.kd, 0.2, 0.5, stepsib=spec.stepsib, uz=spec.uz, mz_frac=spec.mz_frac)
    comp = m.components(p)
    np.testing.assert_allclose(comp.AxC, 0.1 * K * 0.5 ** m.f)
    assert (comp.AxC[df.rel.str.startswith("adoptive").to_numpy()] == 0).all()
    np.testing.assert_allclose(m.predict(p), comp.sum(axis=1))
    q = dict(p, rho_g=0.9, rho_e=0.8, k_side=0.2)
    mult = np.where(m.os, 0.9 * 0.8, 1.0) * (1 + 0.2 * (m.side_mat.astype(float) - m.side_pat.astype(float)))
    np.testing.assert_allclose(m.components(q).AxC, comp.AxC * mult)
