# X-linked additive kernel for the Eftedal et al. (2025) categories

**Bottom line.**
- **Fixed share.** X carries about 2% of V_A (recommended fixed $s_X = 0.02$, plausible range 0.01–0.05).
- **Siblings and half siblings:** X predicts the right signs, but at most ~0.002–0.005 in correlation. That explains ≤ 17% of the sibling opposite-sex deficit and ≤ 18% of the maternal half-sib excess.
- **Cousins:** X could account for about a quarter of the small maternal > paternal excess at central values, and most of it only at the upper end ($s_X = 0.05$, $V_A = 0.7$). It cannot explain the cousins' same-sex excess, and it predicts the wrong sign for half cousins.
- **For the Phase-2 model:** add $s_X D_k(m)$, interpolated from `x_kernel_grid.csv`, to the autosomal kernel.

Full writeup with derivations, checks and figures: `x_kernel.pdf` (source `x_kernel.tex`). Files (all in this folder):

| file | content |
|---|---|
| `x_kinship.csv` | one row per (rel, side, sex): $\phi^X$, $a$, $K^X$ at random mating (FDC and NDC), mtDNA / Y line sharing, exact $K^X$ under AM at $m = 0.2, 0.43$, and two approximations |
| `x_configs.csv` | every pedigree configuration a category averages over, with its weight |
| `x_kernel_grid.csv` | exact $K^{auto}$, $K^X$ and $D = K^X - K^{auto}$ for $m = 0, 0.05, \dots, 0.6$ and $0.43$ (for the model) |
| `signatures.csv` | predicted vs observed contrasts |
| `pedigrees.py`, `x_kernel.py`, `compute.py`, `grid.py`, `signatures.py` | code (`PYTHONPATH=<worktree>/src python compute.py && python grid.py && python signatures.py`) |
| `make_tables.py`, `checks.py`, `figure.R`, `tex/`, `fig/` | tables, checks and figures for `x_kernel.pdf` |

## 1. What `side` means
**Rule** (Eftedal et al. 2025, Fig. 1 legend): a pair is maternal or paternal when the connecting pair one generation above them is two women or two men. It is `other` when that pair is a woman and a man, or when the pair is connected through both parents (full sibs).

| type | maternal | paternal | other |
|---|---|---|---|
| half sibs | shared mother | shared father | – |
| 1st cousins | parents are sisters | parents are brothers | brother and sister |
| 1st cousins t/uz | parents are twin sisters | twin brothers | – |
| half cousins | parents are half sisters | half brothers | one of each |
| 2nd cousins (and t/uz) | parents are two female 1st cousins | two male 1st cousins | one of each |
| co-cousins | the connecting mates are on both children's mothers' side | fathers' side | mixed |
| step sibs | – | – | always |

**What is averaged** (equal weights; listed per configuration in `x_configs.csv`):
- **Sexes the data do not record:** the sib pair at the top of 2nd cousins (FF, FM, MF, MM), the shared grandparent of half cousins, and the twins' sex for 2nd cousins t/uz.
- **Zygosity:** UZ and t/uz twins are ½ MZ + ½ DZ, the model's convention.
- **Orientation:** an `os` pair is averaged over which member is female, which matters only for asymmetric pedigrees.
- **Alternative weights:** using the observed maternal:paternal half-sib ratio (53:47) instead of 50:50 for half cousins changes their $K^X$ by ≤ 0.005.
- **Unresolved:** no side definition was ambiguous.

## 2. Derivation (random mating)
X kinship $\phi^X_{ij}$ is the probability that random X alleles of $i$ and $j$ are IBD, where a male has a single X. It is computed by the standard recursion on the younger member:

$$\phi^X_{\text{male }i,\,j} = \phi^X_{\text{mother}(i),\,j},\qquad \phi^X_{\text{female }i,\,j} = \tfrac12\big(\phi^X_{\text{mother}(i),\,j} + \phi^X_{\text{father}(i),\,j}\big),\qquad \phi^X_{ii} = 1\ (\text{male}),\ \tfrac12\ (\text{female}).$$

The kernel is the X genetic covariance divided by the **sex-averaged** X variance, so it sits on the autosomal kernel's scale (MZ averages 1 over the sexes):

$$\text{FDC: } K^X = \tfrac43\,\phi^X,\qquad \text{NDC: } K^X = \tfrac{4}{1.5}\,\phi^X\,2^{-\#\text{males in the pair}}.$$

Random-mating FDC kernels ($a$ = autosomal relatedness):

| type (a) | side | ssf | ssm | os |
|---|---|---|---|---|
| MZ (1) | | 0.667 | 1.333 | |
| sibs, DZ (1/2) | | 0.500 | 0.667 | 0.333 |
| half sibs (1/4) | maternal | 0.167 | 0.667 | 0.333 |
| | paternal | 0.333 | 0 | 0 |
| 1st cousins (1/8) | maternal | 0.125 | 0.500 | 0.250 |
| | paternal | 0.167 | 0 | 0 |
| | other | 0.083 | 0 | 0.083 |
| 2nd cousins (1/32) | maternal | 0.029 | 0.115 | 0.057 |
| | paternal | 0.031 | 0 | 0 |
| | other | 0.026 | 0 | 0.026 |

In-laws and adoptive pairs are 0 at random mating.

## 3. Dosage compensation and $s_X$
- **Convention: full dosage compensation (FDC).** Male X variance is 2× female. [Sidorenko et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6614401/) find a male/female ratio of X-linked $h^2_{SNP}$ of 1.9–2.2 across 20 UK Biobank traits; the X-relationship conventions follow Yang et al. 2011 and Fernando & Grossman 1990.
- **Alternative: no dosage compensation (NDC).** It reverses the male/female ordering (sisters get the largest X sharing). It changes magnitudes by less than a factor of 2.5 and does not change any conclusion below (`D_contrast_NDC_m0.0` in `signatures.csv`).
- **Expected share:** X is 5.0% of the haploid genome (156 of ~3,100 Mb; ~4% of protein-coding genes). Under FDC with per-active-copy effects like the autosomes, an X locus contributes ½ (males) and ¼ (females) of an autosomal locus's variance. That gives $s_X \approx 0.05 \times 0.375 / 0.97 \approx 0.02$.
- **Empirical check:** UK Biobank X $h^2_{SNP}$ is 0.63% (males) and 0.30% (females) of phenotypic variance, i.e. about 1–3% of typical total $h^2_{SNP}$.
- **Recommendation: fix $s_X = 0.02$**, with range 0.01–0.05. The upper end allows for brain or intellectual-disability gene enrichment on X and for escape from X inactivation.

## 4. How it enters the model, and assortative mating
- **Exact treatment under AM:** I built a pathMgr model with autosomal and X additive values per person and sex-specific X transmission. Under FDC a son receives his mother's X value with coefficient 1 and his father's with 0; a daughter receives ½ from each. Mates match on $Y$, and the equilibrium is solved separately for each sex.
- **Definition of $K^X(m)$:** at matched genotypic mate correlation $m$, $K^X$ is defined by $K^{total} = (1-s)K^{auto} + sK^X$. It does not depend on $s$ (sibs os: 0.5899, 0.5901, 0.5903 at $s$ = 0.005, 0.02, 0.05). It depends only weakly on $h^2$ at fixed $m$ (≤ 0.034 in $D$, i.e. ≤ 0.0005 in correlation).
- **Recommended implementation:** $K^{total}_k(m) = K^{auto}_k(m) + s_X\,D_k(m)$, with $D = K^X - K^{auto}$ linearly interpolated from `x_kernel_grid.csv`. $K^{auto}_k$ stays your own formula (Nagylaki half sibs, corrected step sibs, and so on). At $m = 0$, $D_k = K^X_{FDC} - a_k$ exactly.
- **Approximations and their errors:**
    - *Autosomal AM factor on the random-mating X kernel* ($K^X(0)\,K^{auto}(m)/K^{auto}(0)$): poor. Errors in $K^X$ reach 0.75 at $m = 0.43$ (median 0.14). AM's extra sharing comes from mate correlation on the *total* genetic value, and X takes part in it about additively, not proportionally.
    - *Constant X correction* ($D_k(m) \approx D_k(0)$): error in $K^X$ ≤ 0.19 (median 0.05), i.e. ≤ 0.002 in correlation at $s_X = 0.02$, $V_A = 0.5$. This is acceptable if interpolation is inconvenient.
- **In-laws:** step sibs and co-cousins get nonzero $D(m)$ under AM, because X transmission is sex-specific. For example, paternal co-cousins ssm have $D = -0.11$. These effects are ≤ 0.002 in correlation.
- **GxG:** no X version is needed. X-involved pairs are about $2s_X$ ≈ 4% of GxG variance, and the change to the $a^2$ kernel is ≤ 0.008 kernel units (≤ 0.002 in correlation).
- **mtDNA:** a maternal-line indicator (`mt_share`: 1 for sibs, maternal half sibs and maternal 1st cousins; ½ for maternal half cousins; ¼ for maternal 2nd cousins). Its expected variance for cognition is negligible, and it is confounded with any maternal-line environment.
- **Y:** ignored (`y_share` is given; it is nonzero only for paternal-line male pairs).

## 5. Predicted signatures vs the data
X adds $V_A\, s_X\, D_k(m)$ to each correlation. Contrasts are computed identically for data and predictions:
- ss − os: per side, mean(ssf, ssm) − os, then averaged over sides;
- maternal − paternal: per sex composition, then averaged.

Predictions are at $m = 0.43$. My aggregation reproduces the task's observed values to within 0.002.

| contrast | observed (SE) | X: $V_A$ 0.4–0.7, $s_X$ 0.02 | X: $s_X$ 0.05 |
|---|---|---|---|
| ss − os, sibs | 0.030 (0.003) | 0.001–0.002 | 0.003–0.005 |
| ss − os, half sibs | 0.007 (0.009) | 0.0006–0.001 | 0.001–0.003 |
| ss − os, 1st cousins | 0.010 (0.002) | 0.0002–0.0003 | ≤ 0.0007 |
| ss − os, 2nd cousins | 0.006 (0.002) | ≈ 0 | ≤ 0.0001 |
| mat − pat, half sibs | 0.050 (0.009) | 0.002–0.004 | 0.005–0.009 |
| mat − pat, 1st cousins | 0.014 (0.003) | 0.003–0.005 | 0.007–0.012 |
| mat − pat, 2nd cousins | 0.0095 (0.003) | 0.0014–0.0025 | 0.004–0.006 |
| mat − pat, half cousins | −0.013 (0.011) | +0.002–0.004 | +0.006–0.010 |
| mat − pat, co-cousins | 0.002 (0.003) | 0.0007–0.0013 | 0.002–0.003 |

**Verdict:**
- **Siblings and half siblings:** X is far too small, as expected (≤ 17–18%). The sex pattern also does not fit. X puts the maternal half-sib excess in brothers (+0.005) and opposite-sex pairs (+0.003), with sisters slightly *paternal* (−0.0005), whereas the data show it in all three compositions (ssf +0.035, ssm +0.062, os +0.052).
- **Cousins:** X can plausibly explain part of the small maternal excess: about 15–35% at central values, and most of it only at $s_X = 0.05$, $V_A = 0.7$. AM roughly triples the 2nd-cousin X contrast relative to random mating ($D$ 0.056 → 0.177), because mother-to-son X transmission carries AM-induced sharing with coefficient 1. X does not explain the cousins' same-sex excess (≤ 7%), and the half-cousin contrast has the opposite sign (within about 1.5 SE of X's prediction).
- **Overall:** the maternal and same-sex patterns need a non-X route, such as maternal-line environment or sex-specific environmental sharing.

## 6. Caveats
- **Assumptions:** primary phenotypic AM with Gaussian assortment at equilibrium; no X inactivation escape; unlinked loci; no sex-specific autosomal effects.
- **Duplicated row:** the 2nd cousins / other / ssf row is the suspect duplicate (`data/README.md`). It enters only the 2nd-cousin ss − os contrast.
- **Equal sex weights:** unrecorded sexes are averaged with equal weights. A 10-point change in any of those weights moves $K^X$ by ≤ 0.02, i.e. ≤ 0.0003 in correlation.
