"""Allele-level forward simulation of pairwise (A x A) epistasis under primary phenotypic
assortative mating, to test which relative-covariance kernel GxG follows.

Model (per individual): standardized genotypes x_k = (g_k - 2p_k) / sqrt(2 p_k q_k), with the base
generation's allele frequencies (pinned base units), M unlinked loci, and
    A = sum_k beta_k x_k
    G = sum_{(k,l) in P} gamma_kl x_k x_l          (architecture "random": random pairs P, gamma
                                                     independent of beta)
      = lam/2 * (A^2 - sum_k beta_k^2 x_k^2)        (architecture "aligned": gamma_kl = lam beta_k
                                                     beta_l over all pairs, i.e. scale epistasis)
    Y = A + G + E
Mates are paired by rank-matching Y against a noisy copy of the partner's Y (realized mate
correlation is measured). Every primary couple has two children; in the last mating round a
fraction of individuals also has a second, assortatively chosen union with one child, which
creates half and step siblings. Relative classes are enumerated from the last three generations.

    python simulate.py --arch random --reps 3 --out results/

Writes <out>/<arch>_classes.csv (one row per replicate x class) and <arch>_population.csv.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


class Sim:
    def __init__(self, arch, N, M, P, V_A0, V_G0, V_E, noise, seed):
        self.rng = np.random.default_rng(seed)
        r = self.rng
        self.arch, self.N, self.M = arch, N, M
        self.p = r.uniform(0.1, 0.9, M)
        self.scale = np.sqrt(2 * self.p * (1 - self.p))
        self.beta = r.normal(0, np.sqrt(V_A0 / M), M)
        self.V_E, self.noise = V_E, noise
        if arch == "random":
            k = r.integers(0, M, 3 * P)
            l = r.integers(0, M, 3 * P)
            keep = k < l
            pairs = np.unique(np.stack([k[keep], l[keep]], 1), axis=0)[:P]
            self.pairs = pairs
            self.gamma = r.normal(0, np.sqrt(V_G0 / len(pairs)), len(pairs))
        elif arch == "aligned":
            sA2 = np.sum(self.beta ** 2)
            # Var[lam/2 (A^2 - sum b^2 x^2)] = lam^2/2 (sA2^2 - sum b^4) under linkage equilibrium
            self.lam = np.sqrt(2 * V_G0 / (sA2 ** 2 - np.sum(self.beta ** 4)))
        else:
            raise ValueError(arch)

    # ---- genetic values ---------------------------------------------------------------
    def x(self, g):
        return (g - 2 * self.p) / self.scale

    def values(self, g):
        x = self.x(g.astype(np.float32))
        A = x @ self.beta
        if self.arch == "random":
            G = np.zeros(len(g))
            for s in range(0, len(self.pairs), 2000):
                pk = self.pairs[s:s + 2000]
                G += (x[:, pk[:, 0]] * x[:, pk[:, 1]]) @ self.gamma[s:s + 2000]
        else:
            G = self.lam / 2 * (A ** 2 - (x ** 2) @ (self.beta ** 2))
        return A, G

    # ---- reproduction ------------------------------------------------------------------
    def gamete(self, g):
        return (self.rng.random(g.shape) < g / 2).astype(np.int8)

    def child(self, gm, gf):
        return self.gamete(gm) + self.gamete(gf)

    def match(self, ym, yf):
        """Rank-match: male i is paired with the female whose noisy score has his rank."""
        sd = yf.std()
        score_f = yf + self.rng.normal(0, self.noise * sd, len(yf))
        om, of = np.argsort(ym), np.argsort(score_f)
        partner = np.empty(len(ym), dtype=int)
        partner[om] = of
        return partner

    def founders(self):
        return (self.rng.random((self.N, self.M)) < self.p).astype(np.int8) + \
               (self.rng.random((self.N, self.M)) < self.p).astype(np.int8)


def run(arch, N, M, P, V_A0, V_G0, V_E, noise, gens, second_frac, seed):
    sim = Sim(arch, N, M, P, V_A0, V_G0, V_E, noise, seed)
    r = sim.rng
    g = sim.founders()
    sex = np.arange(N) % 2                                    # 0 female, 1 male
    hist = []                                                 # per generation: dict
    pop_rows = []
    for t in range(gens + 1):
        A, G = sim.values(g)
        E = r.normal(0, np.sqrt(V_E), len(g))
        Y = A + G + E
        hist.append(dict(g=g, A=A, G=G, Y=Y, sex=sex))
        if t == gens:
            break
        fem, mal = np.where(sex == 0)[0], np.where(sex == 1)[0]
        n = min(len(fem), len(mal))
        fem, mal = r.permutation(fem)[:n], r.permutation(mal)[:n]
        partner = sim.match(Y[mal], Y[fem])                    # mal[i] <-> fem[partner[i]]
        moms, dads = fem[partner], mal
        rho_y = np.corrcoef(Y[moms], Y[dads])[0, 1]
        pop_rows.append(dict(gen=t, V_A=A.var(), V_G=G.var(), c_AG=np.cov(A, G)[0, 1],
                             V_Y=Y.var(), rho_y=rho_y,
                             mate_AA=np.cov(A[moms], A[dads])[0, 1],
                             mate_GG=np.cov(G[moms], G[dads])[0, 1],
                             mate_AG=(np.cov(A[moms], G[dads])[0, 1] + np.cov(G[moms], A[dads])[0, 1]) / 2))
        couples = [(m, d) for m, d in zip(moms, dads)]
        last_round = t == gens - 1
        if last_round and second_frac > 0:
            # second unions: a male from one primary couple with a female from another, matched
            # assortatively on Y; one child each.
            k = int(second_frac * n)
            pm = r.choice(n, k, replace=False)                 # couples donating the male
            pf = r.choice(n, k, replace=False)                 # couples donating the female
            ok = pm != pf
            pm, pf = pm[ok], pf[ok]
            males, females = dads[pm], moms[pf]
            p2 = sim.match(Y[males], Y[females])
            second = [(females[p2[i]], males[i]) for i in range(len(males))]
        else:
            second = []
        kids, mom_id, dad_id, prim = [], [], [], []
        for (m, d) in couples:
            for _ in range(2):
                mom_id.append(m); dad_id.append(d); prim.append(True)
        for (m, d) in second:
            mom_id.append(m); dad_id.append(d); prim.append(False)
        mom_id, dad_id = np.array(mom_id), np.array(dad_id)
        g = sim.child(hist[-1]["g"][mom_id], hist[-1]["g"][dad_id])
        sex = r.integers(0, 2, len(g))
        hist[-1]["kid_parents"] = (mom_id, dad_id, np.array(prim))
    return sim, hist, pd.DataFrame(pop_rows)


# ---- relative classes -------------------------------------------------------------------
def classes(hist):
    """Pairs (generation, index) for each class, from the last three generations."""
    T = len(hist) - 1
    out = {}
    momT, dadT, primT = hist[T - 1]["kid_parents"]           # parents (in gen T-1) of gen T
    mom1, dad1, _ = hist[T - 2]["kid_parents"]               # parents (in gen T-2) of gen T-1
    mom0, dad0, _ = hist[T - 3]["kid_parents"]               # parents (gen T-3) of gen T-2

    def pairs_by(keys, idx):
        df = pd.DataFrame({"k": keys, "i": idx})
        res = []
        for _, grp in df.groupby("k"):
            ids = grp["i"].to_numpy()
            for a in range(len(ids)):
                for b in range(a + 1, len(ids)):
                    res.append((ids[a], ids[b]))
        return np.array(res, dtype=int).reshape(-1, 2)

    nT = len(momT)
    iT = np.arange(nT)
    # full sibs: same mother and father (primary couples have 2 kids)
    fs = pairs_by(momT * 10 ** 7 + dadT, iT)
    out["full sibs"] = (T, fs, T)
    # half sibs: share exactly one parent
    hs_m = pairs_by(momT, iT)
    hs_d = pairs_by(dadT, iT)
    hs = np.concatenate([hs_m, hs_d])
    hs = hs[(momT[hs[:, 0]] != momT[hs[:, 1]]) | (dadT[hs[:, 0]] != dadT[hs[:, 1]])]
    out["half sibs"] = (T, hs, T)
    # parent-offspring (gen T-1 mother, gen T child), primary kids
    po = np.stack([momT[primT], iT[primT]], 1)
    out["parent-offspring"] = (T - 1, po, T)
    # step sibs: child of (fA, mA) and child of (fB, mB) where mA--fB is a second union
    second = np.where(~primT)[0]
    prim_couple = {}
    for i in np.where(primT)[0]:
        prim_couple.setdefault(momT[i], []).append(i)
        prim_couple.setdefault(-1 - dadT[i], []).append(i)
    st = []
    for i in second:
        kidsA = prim_couple.get(-1 - dadT[i], [])            # kids of the male's primary couple
        kidsB = prim_couple.get(momT[i], [])                 # kids of the female's primary couple
        for a in kidsA:
            for b in kidsB:
                if momT[a] != momT[b] and dadT[a] != dadT[b]:
                    st.append((a, b))
    out["step sibs"] = (T, np.array(st, dtype=int).reshape(-1, 2), T)
    # first cousins: gen-T kids (primary) whose parents are full sibs in gen T-1
    sib_of = {}
    n1 = len(mom1)
    for a, b in pairs_by(mom1 * 10 ** 7 + dad1, np.arange(n1)):
        sib_of.setdefault(a, []).append(b)
        sib_of.setdefault(b, []).append(a)
    kids_of = {}
    for i in np.where(primT)[0]:
        kids_of.setdefault(momT[i], set()).add(i)
        kids_of.setdefault(dadT[i], set()).add(i)
    fc, cc = set(), set()
    for p, sibs in sib_of.items():
        for s in sibs:
            for a in kids_of.get(p, ()):
                for b in kids_of.get(s, ()):
                    if a < b:
                        fc.add((a, b))
    out["first cousins"] = (T, np.array(sorted(fc), dtype=int).reshape(-1, 2), T)
    # co-cousins: kids of p1 and p2, where p1's sib s1 is mated (primary) with t1, p2's sib
    spouse = {}
    for i in np.where(primT)[0]:
        spouse.setdefault(momT[i], set()).add(dadT[i])
        spouse.setdefault(dadT[i], set()).add(momT[i])
    for s1, t1s in spouse.items():
        for t1 in t1s:
            for p1 in sib_of.get(s1, ()):
                for p2 in sib_of.get(t1, ()):
                    for a in kids_of.get(p1, ()):
                        for b in kids_of.get(p2, ()):
                            if a != b and momT[a] != momT[b] and dadT[a] != dadT[b]:
                                cc.add((min(a, b), max(a, b)))
    out["co-cousins"] = (T, np.array(sorted(cc), dtype=int).reshape(-1, 2), T)
    # avuncular: gen T-1 individual and gen-T child of their full sib
    av = [(s, k) for p, sibs in sib_of.items() for s in sibs for k in kids_of.get(p, ())]
    out["avuncular"] = (T - 1, np.array(av, dtype=int).reshape(-1, 2), T)
    # grandparent: gen T-2 grandmother of gen T (via primary mother)
    gp = np.stack([mom1[momT[primT]], iT[primT]], 1)
    out["grandparent"] = (T - 2, gp, T)
    return out


def measure(hist, cls):
    rows = []
    for name, (ga, pairs, gb) in cls.items():
        if len(pairs) < 50:
            continue
        Ha, Hb = hist[ga], hist[gb]
        A_a, G_a = Ha["A"] - Ha["A"].mean(), Ha["G"] - Ha["G"].mean()
        A_b, G_b = Hb["A"] - Hb["A"].mean(), Hb["G"] - Hb["G"].mean()
        i, j = pairs[:, 0], pairs[:, 1]
        sA = np.sqrt(A_a.var() * A_b.var())
        sG = np.sqrt(G_a.var() * G_b.var())
        rows.append(dict(cls=name, n_pairs=len(pairs),
                         K_A=np.mean(A_a[i] * A_b[j]) / sA,
                         K_G=np.mean(G_a[i] * G_b[j]) / sG,
                         cross=(np.mean(A_a[i] * G_b[j]) + np.mean(G_a[i] * A_b[j])) / 2 / np.sqrt(sA * sG)))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--arch", choices=["random", "aligned"], required=True)
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--N", type=int, default=20000)
    ap.add_argument("--M", type=int, default=1000)
    ap.add_argument("--P", type=int, default=4000)
    ap.add_argument("--VA", type=float, default=0.45)
    ap.add_argument("--VG", type=float, default=0.30)
    ap.add_argument("--VE", type=float, default=0.25)
    ap.add_argument("--noise", type=float, default=0.6, help="partner-score noise (sets rho_y)")
    ap.add_argument("--gens", type=int, default=10)
    ap.add_argument("--second", type=float, default=0.3)
    ap.add_argument("--out", type=Path, default=Path(__file__).parent / "results")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    cls_rows, pop_rows = [], []
    for rep in range(a.reps):
        sim, hist, pop = run(a.arch, a.N, a.M, a.P, a.VA, a.VG, a.VE, a.noise, a.gens, a.second, 1000 + rep)
        m = measure(hist, classes(hist)).assign(rep=rep)
        cls_rows.append(m)
        pop_rows.append(pop.assign(rep=rep))
        print(a.arch, "rep", rep, "done", flush=True)
        print(m.round(4).to_string(index=False), flush=True)
    pd.concat(cls_rows).to_csv(a.out / f"{a.arch}_classes.csv", index=False)
    pd.concat(pop_rows).to_csv(a.out / f"{a.arch}_population.csv", index=False)


if __name__ == "__main__":
    main()
