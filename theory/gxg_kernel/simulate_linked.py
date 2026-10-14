"""Haplotype-level forward simulation of A + pairwise GxG under AM with a genetic map.

Extends simulate.py (same mating scheme, relative classes and measures) with phased haplotypes
and recombination. Map: --chrom chromosomes of 1 Morgan, each carrying genes of --per-gene loci
spread over --gene-span Morgans (gene starts uniform on the chromosome); Haldane map, no
interference, so adjacent-locus switches are independent with c_i = (1 - exp(-2 d_i))/2.

Interaction architectures, G = sum_{k<l} gamma_kl x_k x_l (scaled to Var[G] = V_G0 at LE):
    uniform      gamma_kl ~ N(0, 1)                                   (all pairs)
    beta         gamma_kl = beta_k beta_l delta_kl, delta ~ N(0, 1)
    beta_1mc     beta_k beta_l delta_kl (1 - c_kl)                    (the proposed proximity weight)
    beta_prox    beta_k beta_l delta_kl exp(-d_kl/lam), same chromosome only
    beta_prox_dir  beta_k beta_l exp(-d_kl/lam), same chromosome only (directional: delta = 1)

    python simulate_linked.py --arch beta_prox --lam 0.001 --reps 3 --out results/linked/strong

Writes <out>/<tag>_classes.csv, <tag>_population.csv, <tag>_bins.csv (per replicate: base variance
of G by recombination bin, for the path-model prediction); <tag> = arch[_lam].
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from simulate import classes, measure

C_EDGES = np.array([0, 1e-3, 5e-3, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.49999, 0.5 + 1e-9])


class LinkedSim:
    def __init__(self, arch, lam, N, n_chrom, per_gene, genes_per_chrom, gene_span,
                 V_A0, V_G0, V_E, noise, seed):
        self.rng = r = np.random.default_rng(seed)
        self.arch, self.N, self.V_E, self.noise = arch, N, V_E, noise
        chrom, pos = [], []
        for ch in range(n_chrom):
            starts = np.sort(r.uniform(0, 1 - gene_span, genes_per_chrom))
            for s in starts:
                pos.extend(np.sort(s + r.uniform(0, gene_span, per_gene)))
                chrom.extend([ch] * per_gene)
        self.chrom, self.pos = np.array(chrom), np.array(pos)
        M = self.M = len(pos)
        # switch probability before each locus (0.5 at chromosome starts = random first haplotype)
        d = np.diff(self.pos, prepend=0.0)
        self.c_int = (1 - np.exp(-2 * d)) / 2
        self.c_int[np.r_[True, self.chrom[1:] != self.chrom[:-1]]] = 0.5
        self.p = r.uniform(0.1, 0.9, M)
        self.scale = np.sqrt(2 * self.p * (1 - self.p))
        self.beta = r.normal(0, np.sqrt(V_A0 / M), M)
        same = self.chrom[:, None] == self.chrom[None, :]
        dist = np.where(same, np.abs(self.pos[:, None] - self.pos[None, :]), np.inf)
        self.c_pair = np.where(same, (1 - np.exp(-2 * np.where(same, dist, 0))) / 2, 0.5)
        bb = np.outer(self.beta, self.beta)
        delta = r.normal(0, 1, (M, M))
        delta = np.triu(delta, 1) + np.triu(delta, 1).T
        if arch == "uniform":
            Gam = delta
        elif arch == "beta":
            Gam = bb * delta
        elif arch == "beta_1mc":
            Gam = bb * delta * (1 - self.c_pair)
        elif arch == "beta_prox":
            Gam = bb * delta * np.exp(-dist / lam)
        elif arch == "beta_prox_dir":
            Gam = bb * np.exp(-dist / lam)
        else:
            raise ValueError(arch)
        np.fill_diagonal(Gam, 0.0)
        Gam *= np.sqrt(V_G0 / (np.sum(Gam ** 2) / 2))             # Var[G] = sum_{k<l} gamma^2 at LE
        self.Gam = Gam.astype(np.float32)
        # per-chromosome blocks when no cross-chromosome pairs (proximity archs)
        self.blockwise = arch in ("beta_prox", "beta_prox_dir")
        iu = np.triu_indices(M, 1)
        w = Gam[iu] ** 2
        cb = self.c_pair[iu]
        idx = np.digitize(cb, C_EDGES) - 1
        self.bins = pd.DataFrame(dict(bin=idx, w=w, cw=w * cb)).groupby("bin").sum().reset_index()
        self.bins["c"] = self.bins.cw / self.bins.w
        self.bins["V0"] = self.bins.w
        self.bins = self.bins[["bin", "c", "V0"]]

    def x(self, h):
        return ((h[:, 0] + h[:, 1]).astype(np.float32) - 2 * self.p.astype(np.float32)) / self.scale.astype(np.float32)

    def values(self, h):
        x = self.x(h)
        A = x @ self.beta.astype(np.float32)
        if self.blockwise:
            G = np.zeros(len(x), np.float32)
            for ch in np.unique(self.chrom):
                s = np.where(self.chrom == ch)[0]
                sl = slice(s[0], s[-1] + 1)
                G += np.einsum("ij,ij->i", x[:, sl] @ self.Gam[sl, sl], x[:, sl])
            G /= 2
        else:
            G = np.zeros(len(x), np.float32)
            for s in range(0, len(x), 4000):
                xs = x[s:s + 4000]
                G[s:s + 4000] = np.einsum("ij,ij->i", xs @ self.Gam, xs) / 2
        return A.astype(float), G.astype(float)

    def gamete(self, h):
        """h: (n, 2, M) -> (n, M) recombinant gamete."""
        n = len(h)
        out = np.empty((n, self.M), np.int8)
        for s in range(0, n, 4000):
            hs = h[s:s + 4000]
            sw = self.rng.random((len(hs), self.M), dtype=np.float32) < self.c_int.astype(np.float32)
            # parity within chromosome: cumulative switches (chromosome starts restart randomly,
            # which is the same as continuing the cumsum since the start switch is w.p. 1/2)
            par = (np.cumsum(sw, axis=1) & 1).astype(bool)
            out[s:s + 4000] = np.where(par, hs[:, 1], hs[:, 0])
        return out

    def child(self, hm, hf):
        return np.stack([self.gamete(hm), self.gamete(hf)], 1)

    def founders(self):
        return (self.rng.random((self.N, 2, self.M)) < self.p).astype(np.int8)

    def match(self, ym, yf):
        if self.noise is None:                                   # random mating
            return self.rng.permutation(len(yf))
        sd = yf.std()
        score_f = yf + self.rng.normal(0, self.noise * sd, len(yf))
        om, of = np.argsort(ym), np.argsort(score_f)
        partner = np.empty(len(ym), dtype=int)
        partner[om] = of
        return partner


def run(sim: LinkedSim, gens, second_frac):
    r = sim.rng
    h = sim.founders()
    N = len(h)
    sex = np.arange(N) % 2
    hist, pop_rows = [], []
    for t in range(gens + 1):
        A, G = sim.values(h)
        E = r.normal(0, np.sqrt(sim.V_E), len(h))
        Y = A + G + E
        hist.append(dict(A=A, G=G, Y=Y, sex=sex))
        if len(hist) > 4:                                        # classes() needs the last 4 gens
            hist[-5] = {k: v for k, v in hist[-5].items() if k == "kid_parents"}
        if t == gens:
            break
        fem, mal = np.where(sex == 0)[0], np.where(sex == 1)[0]
        n = min(len(fem), len(mal))
        fem, mal = r.permutation(fem)[:n], r.permutation(mal)[:n]
        partner = sim.match(Y[mal], Y[fem])
        moms, dads = fem[partner], mal
        pop_rows.append(dict(gen=t, V_A=A.var(), V_G=G.var(), c_AG=np.cov(A, G)[0, 1], V_Y=Y.var(),
                             rho_y=np.corrcoef(Y[moms], Y[dads])[0, 1],
                             skew_A=float(((A - A.mean()) ** 3).mean() / A.std() ** 3)))
        second = []
        if t == gens - 1 and second_frac > 0:
            k = int(second_frac * n)
            pm, pf = r.choice(n, k, replace=False), r.choice(n, k, replace=False)
            ok = pm != pf
            males, females = dads[pm[ok]], moms[pf[ok]]
            p2 = sim.match(Y[males], Y[females])
            second = [(females[p2[i]], males[i]) for i in range(len(males))]
        mom_id = np.r_[np.repeat(moms, 2), np.array([s[0] for s in second], int)]
        dad_id = np.r_[np.repeat(dads, 2), np.array([s[1] for s in second], int)]
        prim = np.r_[np.ones(2 * n, bool), np.zeros(len(second), bool)]
        h = sim.child(h[mom_id], h[dad_id])
        sex = r.integers(0, 2, len(h))
        hist[-1]["kid_parents"] = (mom_id, dad_id, prim)
    return hist, pd.DataFrame(pop_rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arch", required=True,
                    choices=["uniform", "beta", "beta_1mc", "beta_prox", "beta_prox_dir"])
    ap.add_argument("--lam", type=float, default=0.001, help="proximity scale (Morgans)")
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--N", type=int, default=20000)
    ap.add_argument("--chrom", type=int, default=20)
    ap.add_argument("--per-gene", type=int, default=4)
    ap.add_argument("--genes", type=int, default=25, help="genes per chromosome")
    ap.add_argument("--gene-span", type=float, default=0.001, help="Morgans (0.1 cM ~ 100 kb)")
    ap.add_argument("--VA", type=float, default=0.45)
    ap.add_argument("--VG", type=float, default=0.30)
    ap.add_argument("--VE", type=float, default=0.25)
    ap.add_argument("--noise", type=float, default=0.6, help="partner-score noise; <0 = random mating")
    ap.add_argument("--gens", type=int, default=20)
    ap.add_argument("--second", type=float, default=0.3)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    tag = a.arch + (f"_{a.lam:g}" if a.arch.startswith("beta_prox") else "")
    noise = None if a.noise < 0 else a.noise
    cls_rows, pop_rows, bin_rows = [], [], []
    for rep in range(a.reps):
        sim = LinkedSim(a.arch, a.lam, a.N, a.chrom, a.per_gene, a.genes, a.gene_span,
                        a.VA, a.VG, a.VE, noise, 2000 + rep)
        hist, pop = run(sim, a.gens, a.second)
        m = measure(hist, classes(hist)).assign(rep=rep)
        cls_rows.append(m)
        pop_rows.append(pop.assign(rep=rep))
        bin_rows.append(sim.bins.assign(rep=rep))
        print(tag, "rep", rep, flush=True)
        print(m.round(4).to_string(index=False), flush=True)
    pd.concat(cls_rows).to_csv(a.out / f"{tag}_classes.csv", index=False)
    pd.concat(pop_rows).to_csv(a.out / f"{tag}_population.csv", index=False)
    pd.concat(bin_rows).to_csv(a.out / f"{tag}_bins.csv", index=False)


if __name__ == "__main__":
    main()
