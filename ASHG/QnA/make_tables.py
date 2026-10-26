"""Tables and in-text numbers for ASHG/QnA/QnA.tex, read from results/ (scripts 01-08 and
R/poster_figures.R must have run). Writes LaTeX snippets to results/qna/:
  numbers.tex      \\newcommand macros for every number quoted in the text
  t_weights.tex    share of the fit's weight by relative type, per estimator
  t_corrections.tex  Eftedal et al.'s estimator -> revised, step by step (additive model)
  t_inlaw.tex      h2 from biological relatives only vs. with in-laws, per estimator
  t_models.tex     fit statistics of the nine core models (Marginal Type)
  t_eftedal_ci.tex Eftedal et al.'s three CI methods vs. their published intervals
  t_jackknife.tex  jackknife bias vs. SE for selected parameters (interior vs. at a bound)
Run from anywhere: python ASHG/QnA/make_tables.py
"""

from pathlib import Path

import numpy as np
import pandas as pd

from eftedal import data as D
from eftedal import provenance

R = D.REPO_ROOT / "results"
OUT = R / "qna"
OUT.mkdir(parents=True, exist_ok=True)

CORE = ["A", "A>GxG", "A+GxG", "A>F", "A+F", "A+GxG>F", "A+F>GxG", "A>GxG+F", "A+GxG+F"]
BLOCK_LABELS = {"twins": "twins", "sibs": "full-sibs", "half sibs": "half-sibs", "1st cousins": "1st cousins",
                "half cousins": "half cousins", "2nd cousins": "2nd cousins", "step sibs": "step-sibs",
                "cocousins": "co-cousins", "adoptive sibs": "adoptive sibs", "adoptive cousins": "adoptive cousins"}
ADOPT = ["adoptive sibs", "adoptive cousins"]


def model_tex(m: str) -> str:
    """Model name with C for the shared environment and LaTeX arrows."""
    return m.replace("F", "C").replace(">", r"$\to$").replace("+", "+")


macros: dict[str, str] = {}


def num(name: str, value, fmt: str = "{:.2f}") -> None:
    assert name.isalpha(), name
    macros[name] = fmt.format(value) if not isinstance(value, str) else value


def write(name: str, text: str) -> None:
    (OUT / name).write_text(text)


def pct(x: float) -> str:
    if x < 0.001:
        return "$<$0.1\\%"
    return f"{100 * x:.1f}\\%" if x < 0.1 else f"{100 * x:.0f}\\%"


# ---- data ----------------------------------------------------------------------------------------
df = D.load("MRE8")
df["w"] = 1 / df["se"] ** 2
eft = df[D.cohort_mask(df, "ALL_but_adopted_UZ")].copy()
eft["w_eff"] = eft["w"] / eft["corr"] ** 2        # log-scale residual ~ raw residual / rho
num("NPairsMillions", df["npairs"].sum() / 1e6, "{:.1f}")
num("NCategories", len(df), "{:d}")
num("NCategoriesEftedal", len(eft), "{:d}")

# ---- T1: weight shares by relative type ------------------------------------------------------------
rows = []
for b in BLOCK_LABELS:
    a = df[df.cv_block == b]; e = eft[eft.cv_block == b]
    rows.append({"block": b, "pairs": a["npairs"].sum() / df["npairs"].sum(),
                 "eft_nominal": e["w"].sum() / eft["w"].sum(), "eft_eff": e["w_eff"].sum() / eft["w_eff"].sum(),
                 "mp": a["w"].sum() / df["w"].sum(), "mt": 0.1, "in_eft": len(e) > 0})
W = pd.DataFrame(rows)
lines = []
for r in W.itertuples():
    e1 = pct(r.eft_nominal) if r.in_eft else "--"
    e2 = pct(r.eft_eff) if r.in_eft else "--"
    lines.append(f"{BLOCK_LABELS[r.block]} & {pct(r.pairs)} & {e1} & {e2} & {pct(r.mp)} & 10\\% \\\\")
write("t_weights.tex", "\n".join(lines) + "\n")
dist = ["1st cousins", "2nd cousins", "cocousins"]
num("WeightDistantEftNominal", 100 * W[W.block.isin(dist)].eft_nominal.sum(), "{:.0f}")
num("WeightDistantEftEffective", 100 * W[W.block.isin(dist)].eft_eff.sum(), "{:.0f}")
num("WeightDistantPairs", 100 * W[W.block.isin(dist)].pairs.sum(), "{:.0f}")
num("WeightTwinsEftEffective", 100 * W[W.block == "twins"].eft_eff.sum(), "{:.2f}")
num("WeightTwinsEftNominal", 100 * W[W.block == "twins"].eft_nominal.sum(), "{:.1f}")

# ---- T2: corrections, step by step ---------------------------------------------------------------
st = pd.read_csv(R / "corrections" / "steps.csv")
st = st[~st["aside"].astype(bool)]
show = {0: "Eftedal et al.", 1: "Eftedal et al., corrected weights", 2: "Eftedal et al., raw scale",
        3: "Revised (Marginal Pair)", 4: "Revised (Marginal Type)"}
lines = []
for o, g in st.groupby("order"):
    jk = g[g.method == "block jackknife"].set_index("param")
    va = jk.loc["V_A"]
    lines.append(f"{show[o]} & {va.estimate:.2f} & [{max(va.lo, 0):.2f}, {min(va.hi, 1):.2f}] & "
                 f"{jk.loc['m'].estimate:.2f} & {jk.loc['r'].estimate:.2f} & {g.R2_lin.iloc[0]:.3f} & "
                 f"{g.chi2_per_df.iloc[0]:.0f} \\\\")
    key = ["Eft", "EftDelta", "EftLin", "RevPair", "RevType"][o]
    num(f"VA{key}", va.estimate)
    num(f"ChiDf{key}", g.chi2_per_df.iloc[0], "{:.0f}")
    num(f"Rsq{key}", g.R2_lin.iloc[0], "{:.3f}")
write("t_corrections.tex", "\n".join(lines) + "\n")
# Eftedal et al.'s estimator: h2 when one relative type is left out (its jackknife refits)
sj = pd.read_csv(R / "corrections" / "steps_jackknife.csv")
sj = sj[sj.step == "eftedal"].set_index("fold").V_A
num("VAEftNoCocousins", sj["cocousins"])
num("VAEftNoSibs", sj["sibs"])
jkA = st[(st.order == 4) & (st.method == "block jackknife")].set_index("param")
num("VAzeroRevType", jkA.loc["V_A0"].estimate)
num("mRevType", jkA.loc["m"].estimate)
num("rRevType", jkA.loc["r"].estimate)

# ---- T3: in-law contrast ---------------------------------------------------------------------------
il = pd.read_csv(R / "corrections" / "inlaw.csv")
lines = []
for key, g in il.groupby("estimator", sort=False):
    jk = g[g.method == "block jackknife"].set_index("quantity")
    bt = g[g.method == "parametric bootstrap"].set_index("quantity")
    d = jk.loc["difference"]; db = bt.loc["difference"]
    lab = g.label.iloc[0].replace("Eftedal et al.", "Eftedal et al.\\ (Marginal Pair)") if key == "eftedal" else g.label.iloc[0]
    lines.append(f"{lab} & {jk.loc['bio only'].estimate:.2f} & {jk.loc['bio + in-law'].estimate:.2f} & "
                 f"{d.estimate:.2f} & [{d.lo:.2f}, {d.hi:.2f}] & [{db.lo:.2f}, {db.hi:.2f}] \\\\")
    k = {"eftedal": "Eft", "invvar": "RevPair", "hybrid_block": "RevType"}[key]
    num(f"InlawDiff{k}", d.estimate)
    num(f"InlawDiffLo{k}", d.lo)
    num(f"InlawDiffHi{k}", d.hi)
    num(f"InlawBio{k}", jk.loc["bio only"].estimate)
    num(f"InlawAll{k}", jk.loc["bio + in-law"].estimate)
write("t_inlaw.tex", "\n".join(lines) + "\n")

# ---- T4: core models -------------------------------------------------------------------------------
fits = pd.read_csv(R / "models" / "core_fits.csv").query("weights == 'hybrid_block'").set_index("model")
cvs = pd.read_csv(R / "models" / "core_cv_summary.csv").query("weights == 'hybrid_block'").set_index("model")
pair = pd.read_csv(R / "models" / "core_cv_paired.csv").query("weights == 'hybrid_block'").set_index("model")
sc = pd.read_csv(R / "models" / "core_cv_scores.csv").query("weights == 'hybrid_block' and not group_test")
noad = sc[~sc.fold.isin(ADOPT)].groupby("model").score.mean().pow(0.5)
adshare = sc.assign(ad=sc.fold.isin(ADOPT)).groupby("model").apply(
    lambda g: g.loc[g.ad, "score"].sum() / g["score"].sum(), include_groups=False)
rank = cvs["cv_rmse"].rank().astype(int)
lines = []
for m in CORE:
    z = pair.loc[m, "mean_diff"] / pair.loc[m, "se_diff"] if pair.loc[m, "se_diff"] > 0 else np.nan
    zs = "--" if np.isnan(z) else f"{z:.1f}"
    bold = (lambda s: f"\\textbf{{{s}}}") if m == "A+GxG+F" else (lambda s: s)
    lines.append(f"{bold(model_tex(m))} & {fits.loc[m, 'R2_lin']:.3f} & {fits.loc[m, 'chi2_per_df']:.1f} & "
                 f"{cvs.loc[m, 'cv_rmse']:.3f} & {rank[m]} & {zs} & {noad[m]:.3f} & {pct(adshare[m])} \\\\")
write("t_models.tex", "\n".join(lines) + "\n")
for m, k in [("A", "A"), ("A+GxG", "AG"), ("A+F", "AC"), ("A+GxG+F", "AGC"), ("A>GxG", "AtoG"), ("A>F", "AtoC"),
             ("A+F>GxG", "ACtoG")]:
    num(f"CV{k}", cvs.loc[m, "cv_rmse"], "{:.3f}")
    num(f"CVnoad{k}", noad[m], "{:.3f}")
    num(f"ChiDfModel{k}", fits.loc[m, "chi2_per_df"], "{:.1f}")
    num(f"AdShare{k}", 100 * adshare[m], "{:.0f}")
num("ZA", pair.loc["A", "mean_diff"] / pair.loc["A", "se_diff"], "{:.1f}")
num("ZAC", pair.loc["A+F", "mean_diff"] / pair.loc["A+F", "se_diff"], "{:.1f}")
num("ZACtoG", pair.loc["A+F>GxG", "mean_diff"] / pair.loc["A+F>GxG", "se_diff"], "{:.1f}")
num("FoldsWorseA", int(pair.loc["A", "folds_worse"]), "{:d}")
num("ChiDfMin", fits["chi2_per_df"].min(), "{:.0f}")
num("ChiDfMax", fits["chi2_per_df"].max(), "{:.0f}")
jkm = pd.read_csv(R / "models" / "core_jackknife.csv").query("weights == 'hybrid_block'")
for m, k in [("A+GxG+F", "AGC"), ("A>GxG", "AtoG"), ("A+GxG", "AG")]:
    g = jkm[jkm.model == m].set_index("param")
    num(f"VA{k}", g.loc["V_A", "estimate"])
    num(f"VG{k}", g.loc["V_GxG", "estimate"])
    num(f"VGlo{k}", g.loc["V_GxG", "lo"])
    num(f"VGhi{k}", g.loc["V_GxG", "hi"])
    if "V_family" in g.index:
        num(f"VC{k}", g.loc["V_family", "estimate"])
        num(f"VClo{k}", g.loc["V_family", "lo"])
        num(f"VChi{k}", g.loc["V_family", "hi"])
        num(f"c{k}", g.loc["r_f", "estimate"])
    num(f"Vsum{k}", g.loc["V_sum", "estimate"])
num("VsumMin", jkm[jkm.param == "V_sum"].estimate.min())
num("VsumMax", jkm[jkm.param == "V_sum"].estimate.max())

# ---- pooled observed correlations used in the text ---------------------------------------------------
def pooled(rel):
    a = df[df.rel == rel]
    return (a["w"] * a["corr"]).sum() / a["w"].sum()
num("ObsAdoptSibs", pooled("adoptive sibs"))
num("ObsStep", pooled("step sibs"), "{:.3f}")
num("ObsMZ", pooled("mz"))

# step-sibling formula: predicted correlation and V_A under both kernels (additive model, all categories)
ar = pd.read_csv(R / "robustness" / "additive_rows.csv")
ap = pd.read_csv(R / "robustness" / "additive_predictions.csv")
for wt, k in [("hybrid_block", "Type"), ("invvar", "Pair")]:
    ids = ar[(ar.cohort_name == "ALL") & (ar.weighting == wt) & ar.variant.isin(["main", "Eftedal step-sib kernel"])]
    va = ids.set_index("variant").V_A
    num(f"StepVAdiff{k}", abs(va["main"] - va["Eftedal step-sib kernel"]), "{:.3f}")
    for var, kk in [("main", "Corr"), ("Eftedal step-sib kernel", "Eft")]:
        fid = ids[ids.variant == var].fit_id.iloc[0]
        num(f"StepPred{kk}{k}", ap[(ap.fit_id == fid) & (ap.rel == "step sibs")].pred.iloc[0], "{:.3f}")

# ---- T5: Eftedal et al.'s CI methods ---------------------------------------------------------------
ci = pd.read_csv(R / "replication" / "eftedal2025_cis.csv")
meth = {"lm t-interval": "$t$-interval of the log-linear regression (lm)",
        "case bootstrap (prob ~ 1/SE^2)": "bootstrap of categories, drawn $\\propto 1/\\mathrm{SE}^2$",
        "inverse Hessian of RSS": "inverse Hessian of the weighted RSS, $\\pm 2$ SE"}
modl = {"bio_only": "biological only", "bio_inlaw": "biological + in-laws"}
lines = []
for r in ci[ci.param == "h2"].itertuples():
    rep = "yes" if r.method != "inverse Hessian of RSS" else "no"
    lines.append(f"{modl[r.model]} & {meth[r.method]} & {r.estimate:.3f} & [{r.lo:.3f}, {r.hi:.3f}] & "
                 f"[{r.published_lo:.3f}, {r.published_hi:.3f}] & {rep} \\\\")
write("t_eftedal_ci.tex", "\n".join(lines) + "\n")
h = ci[(ci.param == "h2") & (ci.model == "bio_inlaw")].set_index("method")
ratio = (h.loc["case bootstrap (prob ~ 1/SE^2)", "hi"] - h.loc["case bootstrap (prob ~ 1/SE^2)", "lo"]) / \
        (h.loc["inverse Hessian of RSS", "hi"] - h.loc["inverse Hessian of RSS", "lo"])
num("HessianNarrower", ratio, "{:.0f}")

# ---- T6: jackknife bias at interior points vs. at a bound -----------------------------------------
cvp = pd.read_csv(R / "models" / "core_cv_params.csv")
cvp = cvp[cvp.fold.isin(D.CV_BLOCKS) & (cvp.weights == "hybrid_block")]
sym = {"V_A": "$V_A$", "m": "$m$", "r": "$r$", "V_GxG": "$V_{\\mathrm{GxG}}$", "V_family": "$V_C$", "r_f": "$c$"}
sel = [("A", "V_A"), ("A", "m"), ("A+GxG+F", "V_A"), ("A+GxG+F", "V_GxG"), ("A+GxG+F", "V_family"),
       ("A+GxG+F", "r_f"), ("A+F", "m"), ("A+F", "r"), ("A>F", "r_f"), ("A+GxG>F", "r_f")]
lines = []
for m, p in sel:
    th = cvp[cvp.model == m][p].to_numpy(float)
    g = jkm[(jkm.model == m) & (jkm.param == p)].iloc[0]
    B = len(th)
    bias = (B - 1) * (th.mean() - g.estimate)
    where = "at a bound" if g.estimate <= 1e-9 else "interior"
    lines.append(f"{model_tex(m)} & {sym[p]} & {g.estimate:.2f} & {where} & {th.mean():.2f} & {bias:+.3f} & "
                 f"{g.se:.3f} & {abs(bias) / g.se:.1f} \\\\")
write("t_jackknife.tex", "\n".join(lines) + "\n")

# ---- parameter correlations quoted in the text -------------------------------------------------------
pc = pd.read_csv(R / "models" / "core_param_corr.csv").query("weights == 'hybrid_block'")
def rpair(m, a, b):
    g = pc[(pc.model == m) & (((pc.p1 == a) & (pc.p2 == b)) | ((pc.p1 == b) & (pc.p2 == a)))]
    return g.r.iloc[0]
num("rAGCVAVG", rpair("A+GxG+F", "V_A", "V_GxG"))
num("rAGCVAVC", rpair("A+GxG+F", "V_A", "V_family"))
num("rAVAm", rpair("A", "V_A", "m"))

write("numbers.tex", "".join(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in macros.items()))
print(f"wrote {len(macros)} numbers and 6 tables to {OUT}")
provenance.stamp(__file__)
