"""Pedigree configurations for the 82 (rel, side, sex) categories of Eftedal et al. (2025).

Side rule (Eftedal et al. 2025, Fig. 1 legend): a pair is maternal / paternal when the connecting
pair one generation above them is two females / two males, and `other` when it is a female and a
male (or the pair is connected through both parents, as for full siblings). Concretely:

    half sibs        the shared parent is the mother (maternal) or the father (paternal)
    1st cousins      the cousins' parents are sisters (maternal), brothers (paternal), or a
                     brother and a sister (other)
    1st cousins t/uz the parents are same-sex twins: twin sisters (maternal) or twin brothers (paternal)
    half cousins     the parents are half sibs: half sisters / half brothers / one of each
    2nd cousins      the parents are 1st cousins: two women / two men / one of each
    2nd cousins t/uz as 2nd cousins, with same-sex twins as the grandparents' generation
    co-cousins       maternal: the mates connecting them are on both children's mother's side
    step sibs        always other (connected through opposite-sex mates)

Sexes the categories do not record (grandparents' generation, the shared grandparent of half
cousins, twins' sex for 2nd cousins t/uz) are averaged with equal weights; same-sex twins of
unknown zygosity are 1/2 MZ, 1/2 DZ (the model's UZ convention). An `os` pair is averaged over its
two orientations (which member is female), which matters only for asymmetric pedigrees.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import product

F, M = "F", "M"
OPP = {F: M, M: F}
PAIR_SEXES = {"ssf": [(F, F)], "ssm": [(M, M)], "os": [(F, M), (M, F)]}
SIDE_SEXES = {"maternal": [(F, F)], "paternal": [(M, M)], "other": [(F, M), (M, F)]}
ANY2 = [(F, F), (F, M), (M, F), (M, M)]


@dataclass
class Ped:
    ind: dict = field(default_factory=dict)        # key -> (sex, mother, father)
    unions: list = field(default_factory=list)     # (mother, father) incl. childless unions

    def add(self, k, sex, mother=None, father=None):
        if k in self.ind:
            return k
        self.ind[k] = (sex, mother, father)
        if mother is not None and (mother, father) not in self.unions:
            self.unions.append((mother, father))
        return k

    def child(self, k, sex, parent, mate):
        """Child k of `parent` and `mate` (a new founder of the opposite sex)."""
        ps = self.ind[parent][0]
        self.add(mate, OPP[ps])
        mo, fa = (parent, mate) if ps == F else (mate, parent)
        return self.add(k, sex, mo, fa)

    def mate(self, a, b):
        mo, fa = (a, b) if self.ind[a][0] == F else (b, a)
        if (mo, fa) not in self.unions:
            self.unions.append((mo, fa))

    def sex(self, k):
        return self.ind[k][0]


@dataclass
class Config:
    ped: Ped
    i: str
    j: str
    w: float
    label: str


def _sibs(si, sj):
    p = Ped(); p.add("mo", F); p.add("fa", M)
    p.add("i", si, "mo", "fa"); p.add("j", sj, "mo", "fa")
    return p


def _self(s):
    p = Ped(); p.add("i", s)
    return p


def _half_sibs(shared, si, sj):
    p = Ped(); p.add("g", shared)
    p.child("i", si, "g", "x1"); p.child("j", sj, "g", "x2")
    return p


def _first_cousins(sp1, sp2, si, sj, mz=False):
    p = Ped(); p.add("gm", F); p.add("gf", M)
    p.add("p1", sp1, "gm", "gf")
    p2 = "p1" if mz else p.add("p2", sp2, "gm", "gf")
    p.child("i", si, "p1", "o1"); p.child("j", sj, p2, "o2")
    return p


def _half_cousins(g, sp1, sp2, si, sj):
    p = Ped(); p.add("g", g)
    p.child("p1", sp1, "g", "h1"); p.child("p2", sp2, "g", "h2")
    p.child("i", si, "p1", "o1"); p.child("j", sj, "p2", "o2")
    return p


def _second_cousins(sp1, sp2, sq1, sq2, si, sj, mz=False):
    p = Ped(); p.add("gm", F); p.add("gf", M)
    p.add("p1", sp1, "gm", "gf")
    p2 = "p1" if mz else p.add("p2", sp2, "gm", "gf")
    p.child("q1", sq1, "p1", "r1"); p.child("q2", sq2, p2, "r2")
    p.child("i", si, "q1", "o1"); p.child("j", sj, "q2", "o2")
    return p


def _step_sibs(si, sj, via):
    """i child of (am, af), j child of (bm, bf); the connecting union is am x bf (via='am-bf')
    or bm x af. (They have a half sib in common; that child adds nothing genetic to i, j.)"""
    p = Ped(); p.add("am", F); p.add("af", M); p.add("bm", F); p.add("bf", M)
    p.add("i", si, "am", "af"); p.add("j", sj, "bm", "bf")
    if via == "am-bf":
        p.mate("am", "bf")
    else:
        p.mate("bm", "af")
    return p


def _cocousins(sp1, sp2, ss1, si, sj):
    """p1 and its sib s1 (parents g1, g2); p2 and its sib t1 (parents g3, g4); s1 x t1 mated."""
    p = Ped(); p.add("g1", F); p.add("g2", M); p.add("g3", F); p.add("g4", M)
    p.add("p1", sp1, "g1", "g2"); p.add("s1", ss1, "g1", "g2")
    p.add("p2", sp2, "g3", "g4"); p.add("t1", OPP[ss1], "g3", "g4")
    p.mate("s1", "t1")
    p.child("i", si, "p1", "o1"); p.child("j", sj, "p2", "o2")
    return p


def configs(rel: str, side: str, sex: str) -> list[Config]:
    """All pedigree configurations for one category, with weights summing to 1."""
    out = []
    pairs = PAIR_SEXES[sex]

    def emit(builder_list):
        n = len(builder_list)
        for ped, lab, w in builder_list:
            out.append(Config(ped, "i", "j" if "j" in ped.ind else "i", w, lab))

    L = []
    if rel == "mz":
        (si, _), = pairs
        L = [(_self(si), "MZ", 1.0)]
    elif rel == "uz":
        (si, sj), = pairs
        L = [(_self(si), "MZ", 0.5), (_sibs(si, sj), "DZ", 0.5)]
    elif rel in ("dz", "sibs"):
        L = [(_sibs(si, sj), f"i{si} j{sj}", 1 / len(pairs)) for si, sj in pairs]
    elif rel == "half sibs":
        g = F if side == "maternal" else M
        L = [(_half_sibs(g, si, sj), f"shared {g}; i{si} j{sj}", 1 / len(pairs)) for si, sj in pairs]
    elif rel == "1st cousins":
        cs = [(a, b, si, sj) for (a, b) in SIDE_SEXES[side] for (si, sj) in pairs]
        L = [(_first_cousins(a, b, si, sj), f"p {a}{b}; i{si} j{sj}", 1 / len(cs)) for a, b, si, sj in cs]
    elif rel == "1st cousins t/uz":
        (a, b), = SIDE_SEXES[side]
        for zyg, mz in (("MZ", True), ("DZ", False)):
            for si, sj in pairs:
                L.append((_first_cousins(a, b, si, sj, mz=mz), f"{zyg} twins {a}{b}; i{si} j{sj}", 0.5 / len(pairs)))
    elif rel == "half cousins":
        cs = [(g, a, b, si, sj) for g in (F, M) for (a, b) in SIDE_SEXES[side] for (si, sj) in pairs]
        L = [(_half_cousins(g, a, b, si, sj), f"shared gp {g}; p {a}{b}; i{si} j{sj}", 1 / len(cs))
             for g, a, b, si, sj in cs]
    elif rel == "2nd cousins":
        cs = [(a, b, c, d, si, sj) for (a, b) in ANY2 for (c, d) in SIDE_SEXES[side] for (si, sj) in pairs]
        L = [(_second_cousins(a, b, c, d, si, sj), f"gp-sibs {a}{b}; p {c}{d}; i{si} j{sj}", 1 / len(cs))
             for a, b, c, d, si, sj in cs]
    elif rel == "2nd cousins t/uz":
        cs = [(t, c, d, si, sj) for t in (F, M) for (c, d) in SIDE_SEXES[side] for (si, sj) in pairs]
        for zyg, mz in (("MZ", True), ("DZ", False)):
            for t, c, d, si, sj in cs:
                L.append((_second_cousins(t, t, c, d, si, sj, mz=mz),
                          f"{zyg} twins {t}{t}; p {c}{d}; i{si} j{sj}", 0.5 / len(cs)))
    elif rel == "step sibs":
        cs = [(si, sj, v) for (si, sj) in pairs for v in ("am-bf", "bm-af")]
        L = [(_step_sibs(si, sj, v), f"union {v}; i{si} j{sj}", 1 / len(cs)) for si, sj, v in cs]
    elif rel == "cocousins":
        cs = [(a, b, s1, si, sj) for (a, b) in SIDE_SEXES[side] for s1 in (F, M) for (si, sj) in pairs]
        L = [(_cocousins(a, b, s1, si, sj), f"p {a}{b}; s1 {s1}; i{si} j{sj}", 1 / len(cs))
             for a, b, s1, si, sj in cs]
    elif rel.startswith("adoptive"):
        L = []                                   # genetically unrelated: kernel 0
    else:
        raise ValueError(rel)
    emit(L)
    assert not out or abs(sum(c.w for c in out) - 1) < 1e-12, rel
    return out
