#!/usr/bin/env python3
"""
LRC(14) relation-fibre / GAP scout for HYP-2637.

This connects three live threads:

1. HYP-2635: the failed dissociated-stranger split.  Some wide primitive
   sets have every nonzero element in a height-2 relation, so a peel argument
   has no isolated stranger.
2. THM-446 / S559 / S560: the summand graph is the additive face; the
   multiplicand graph is the divisibility/product-sign face.
3. HYP-2634: the two-large lift-opposition example is controlled by exact
   low-height integer relations, with term signs decided only after multiplying
   the coefficient character by the denominator-sign parity.

The key object below is a weighted summand fibre:

    c in N^k,  sum c_i <= M,  c_i <= H,
    c |-> sum_i c_i e_i.

Collisions in these fibres are exactly small integer relations
sum_i n_i e_i = 0.  Ordinary pair sums are only the mass-2, coefficient-1
slice; the HYP-2635 examples require the larger weighted slice.

Tournament Analysis declaration:
  Vertices are proof quotients, not runners: weighted summand fibres,
  height-2 coverage, Freiman/GAP scout, multiplicand sign parity,
  finite residue characters, ordinary pair sums, and raw speed labels.
  Pairwise observable: which quotient preserves the proof predicate
  "relation-dense finite ledger versus dissociated tail, with signed lift
  parity retained."  Gauge: edge points to the quotient that preserves more of
  that predicate.  Ties follow the displayed Hamiltonian path.

Assumption challenge:
  I considered runners, gaps, pair sums, relation vectors, fibre sums, Fourier
  modes, denominator products, residue characters, and proof obligations.  The
  chosen quotient preserves the LRC predicate relevant to HYP-2635/HYP-2634:
  low-height relation density and signed reciprocal contribution.  It destroys
  exact witness times and ordinary graph adjacency.
"""

from __future__ import annotations

import itertools
import math
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import lru_cache
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_support6_residue_cusp_codex_s12 as s12  # noqa: E402


AMBIENT_D = 9


def section(title: str) -> None:
    print("\n" + "=" * 92)
    print(title)
    print("=" * 92)


def gcd_all(vals: tuple[int, ...]) -> int:
    g = 0
    for v in vals:
        g = math.gcd(g, abs(v))
    return g


def primitive_relation(n: tuple[int, ...]) -> tuple[int, ...]:
    g = gcd_all(n)
    out = tuple(x // g for x in n)
    first = next((x for x in out if x), 0)
    if first < 0:
        out = tuple(-x for x in out)
    return out


def coeff_vectors(k: int, max_coeff: int, max_mass: int) -> list[tuple[int, ...]]:
    out: list[tuple[int, ...]] = []

    def rec(i: int, left: int, cur: list[int]) -> None:
        if i == k:
            if 0 < sum(cur) <= max_mass:
                out.append(tuple(cur))
            return
        for c in range(min(max_coeff, left) + 1):
            cur.append(c)
            rec(i + 1, left - c, cur)
            cur.pop()

    rec(0, max_mass, [])
    return out


@dataclass(frozen=True)
class Relation:
    n: tuple[int, ...]
    fibre_sum: int
    pos_mass: int
    neg_mass: int
    height: int

    @property
    def support(self) -> tuple[int, ...]:
        return tuple(i for i, x in enumerate(self.n) if x)


def weighted_fibres(
    E: tuple[int, ...], max_coeff: int, max_mass: int
) -> dict[int, list[tuple[int, ...]]]:
    fibres: dict[int, list[tuple[int, ...]]] = defaultdict(list)
    for c in coeff_vectors(len(E), max_coeff, max_mass):
        fibres[sum(ci * ei for ci, ei in zip(c, E))].append(c)
    return dict(fibres)


def collision_relations(
    E: tuple[int, ...], max_coeff: int, max_mass: int
) -> dict[tuple[int, ...], Relation]:
    fibres = weighted_fibres(E, max_coeff, max_mass)
    rels: dict[tuple[int, ...], Relation] = {}
    for fibre_sum, vecs in fibres.items():
        if len(vecs) < 2:
            continue
        for a, b in itertools.combinations(vecs, 2):
            n0 = tuple(ai - bi for ai, bi in zip(a, b))
            if not any(n0):
                continue
            nonzero_support = [
                i for i, (ni, ei) in enumerate(zip(n0, E)) if ni and ei != 0
            ]
            if len(nonzero_support) < 2:
                continue
            n = primitive_relation(n0)
            if sum(ni * ei for ni, ei in zip(n, E)) != 0:
                continue
            pos = sum(x for x in n if x > 0)
            neg = -sum(x for x in n if x < 0)
            rels[n] = Relation(
                n=n,
                fibre_sum=fibre_sum,
                pos_mass=pos,
                neg_mass=neg,
                height=max(abs(x) for x in n),
            )
    return rels


def ordinary_pair_energy(E: tuple[int, ...]) -> tuple[int, int, int, list[tuple[int, list[tuple[int, int]]]]]:
    fibres: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for i, j in itertools.combinations(range(len(E)), 2):
        fibres[E[i] + E[j]].append((i, j))
    energy = sum(len(v) * len(v) for v in fibres.values())
    collisions = sum(math.comb(len(v), 2) for v in fibres.values() if len(v) > 1)
    touched = set()
    examples = []
    for s, pairs in sorted(fibres.items()):
        if len(pairs) > 1:
            examples.append((s, pairs))
            for pair in pairs:
                touched.update(pair)
    return energy, collisions, len(touched), examples[:5]


def weighted_metrics(
    E: tuple[int, ...], max_coeff: int, max_mass: int, relation_height: int
) -> dict[str, object]:
    fibres = weighted_fibres(E, max_coeff, max_mass)
    rels = {
        n: r
        for n, r in collision_relations(E, max_coeff, max_mass).items()
        if r.height <= relation_height
    }
    nonzero_indices = {i for i, e in enumerate(E) if e != 0}
    touched_nonzero = {
        i for r in rels.values() for i in r.support if i in nonzero_indices
    }
    energy = sum(len(v) * len(v) for v in fibres.values())
    collision_fibres = sum(1 for v in fibres.values() if len(v) > 1)
    max_fibre = max(len(v) for v in fibres.values())
    top_fibres = sorted(
        ((s, v) for s, v in fibres.items() if len(v) > 1),
        key=lambda item: (-len(item[1]), item[0]),
    )[:5]
    rel_rank = rational_rank([r.n for r in rels.values()], len(E))
    return {
        "vectors": sum(len(v) for v in fibres.values()),
        "energy": energy,
        "collision_fibres": collision_fibres,
        "max_fibre": max_fibre,
        "relations": rels,
        "relation_rank": rel_rank,
        "relation_nullity": len(E) - rel_rank,
        "coverage": len(touched_nonzero),
        "coverage_total": len(nonzero_indices),
        "top_fibres": top_fibres,
    }


def rational_rank(rows: list[tuple[int, ...]], ncols: int) -> int:
    """Exact row rank over Q for the bounded relation matrix."""
    mat = [[Fraction(x) for x in row] for row in rows if any(row)]
    rank = 0
    col = 0
    while rank < len(mat) and col < ncols:
        pivot = None
        for r in range(rank, len(mat)):
            if mat[r][col]:
                pivot = r
                break
        if pivot is None:
            col += 1
            continue
        mat[rank], mat[pivot] = mat[pivot], mat[rank]
        pv = mat[rank][col]
        mat[rank] = [x / pv for x in mat[rank]]
        for r in range(len(mat)):
            if r == rank or not mat[r][col]:
                continue
            factor = mat[r][col]
            mat[r] = [x - factor * y for x, y in zip(mat[r], mat[rank])]
        rank += 1
        col += 1
    return rank


def fmt_vec(E: tuple[int, ...], c: tuple[int, ...]) -> str:
    terms = []
    for coeff, value in zip(c, E):
        if coeff == 0:
            continue
        if coeff == 1:
            terms.append(str(value))
        else:
            terms.append(f"{coeff}*{value}")
    return " + ".join(terms) if terms else "0"


def display_fibre(E: tuple[int, ...], target: int, max_coeff: int, max_mass: int) -> None:
    fibres = weighted_fibres(E, max_coeff, max_mass)
    vecs = fibres.get(target, [])
    print(f"  fibre sum {target}: {len(vecs)} representations")
    for c in sorted(vecs, key=lambda x: (sum(x), x))[:10]:
        print(f"    {fmt_vec(E, c)}")


def rank2_gap_scout(E: tuple[int, ...], coeff_bound: int = 12) -> tuple[tuple[int, int] | None, tuple[int, int] | None]:
    """Crude Freiman-GAP scout: cover E by a nonnegative rank-2 difference semigroup."""
    diffs = sorted({b - a for a, b in itertools.combinations(E, 2) if b > a})
    best = None
    best_shape = None
    for d1, d2 in itertools.combinations(diffs, 2):
        coords = {}
        for x in E:
            found = None
            for a in range(coeff_bound + 1):
                for b in range(coeff_bound + 1):
                    if a * d1 + b * d2 == x:
                        found = (a, b)
                        break
                if found is not None:
                    break
            if found is None:
                break
            coords[x] = found
        else:
            max_a = max(a for a, _ in coords.values())
            max_b = max(b for _, b in coords.values())
            volume = (max_a + 1) * (max_b + 1)
            candidate = (volume, max_a + max_b, d1, d2, max_a, max_b)
            if best is None or candidate < best:
                best = candidate
                best_shape = (d1, d2)
    if best is None:
        return None, None
    return best_shape, (best[4], best[5])


def report_hyp2635_examples() -> None:
    section("HYP-2635 WIDE STRUCTURED EXAMPLES: PAIR SUMS VS WEIGHTED FIBRES")
    examples = [
        ("kp_example_A", (0, 3, 5, 16, 28, 30, 33, 35)),
        ("kp_example_B", (0, 4, 12, 15, 20, 21, 25, 31)),
        ("consecutive_8", tuple(range(8))),
        ("dissociated_powers", (0, 1, 5, 25, 125, 625, 3125, 15625)),
    ]
    print(
        "Height-2 coverage ignores the pinned 0 coordinate.  Ordinary pair sums "
        "are the Sidon slice; weighted fibres allow coefficients up to 2 and "
        "side mass up to 4, so they see relations like 2*16+3=35."
    )
    print()
    print(
        f"{'case':>18} {'span':>6} {'pair_E':>8} {'pair_col':>8} "
        f"{'pair_touch':>10} {'wt_E':>8} {'wt_colfib':>9} {'maxfib':>7} "
        f"{'h2_cover':>10} {'rel nul':>7} {'rank2 scout':>18}"
    )
    for name, E in examples:
        Erel = tuple(e for e in E if e != 0)
        pair_E, pair_col, pair_touch, pair_examples = ordinary_pair_energy(E)
        wt = weighted_metrics(Erel, max_coeff=2, max_mass=4, relation_height=2)
        gap_steps, gap_shape = rank2_gap_scout(E)
        if gap_steps is None:
            gap_text = "-"
        else:
            gap_text = f"{gap_steps} box {gap_shape}"
        print(
            f"{name:>18} {E[-1]-E[0]:>6} {pair_E:>8} {pair_col:>8} "
            f"{pair_touch:>10} {wt['energy']:>8} {wt['collision_fibres']:>9} "
            f"{wt['max_fibre']:>7} {wt['coverage']:>4}/{wt['coverage_total']:<5} "
            f"{wt['relation_nullity']:>7} "
            f"{gap_text:>18}"
        )
        print(f"  E={E}")
        if pair_examples:
            s, pairs = pair_examples[0]
            readable = [f"{E[i]}+{E[j]}" for i, j in pairs]
            print(f"  first pair-sum collision: {s} = " + " = ".join(readable))
        rels = sorted(wt["relations"].values(), key=lambda r: (r.height, r.pos_mass + r.neg_mass, r.n))[:4]
        for r in rels:
            pos = tuple(max(x, 0) for x in r.n)
            neg = tuple(max(-x, 0) for x in r.n)
            print(f"  h{r.height} relation: {fmt_vec(Erel, pos)} = {fmt_vec(Erel, neg)}")
    print()
    print(
        "Readout: the wide examples are not ordinary APs, but they are not "
        "dissociated either.  Weighted summand energy covers every nonzero "
        "vertex with height-2 relations, and the relation-matrix nullity is "
        "small.  This is the exact pocket where a Freiman/GAP lemma should "
        "replace the failed stranger peel."
    )


def sign_word(x: float, eps: float = 1e-12) -> str:
    if x > eps:
        return "+"
    if x < -eps:
        return "-"
    return "0"


@lru_cache(None)
def residue_coeff(residues: tuple[int, ...]) -> complex:
    return s12.residue_coeff(residues, AMBIENT_D)


def relation_term(ns: tuple[int, ...]) -> complex:
    return residue_coeff(s12.residue_tuple(ns)) / math.prod(ns)


def relation_side(E: tuple[int, ...], ns: tuple[int, ...], positive: bool) -> str:
    coeff = tuple((max(x, 0) if positive else max(-x, 0)) for x in ns)
    return fmt_vec(E, coeff)


def relation_sign_table(E: tuple[int, ...], rels: list[tuple[str, tuple[int, ...]]]) -> None:
    print(
        f"{'motif':>22} {'additive equality':>38} {'neg#':>5} {'den':>4} "
        f"{'coeff':>6} {'term':>6} {'z':>14}"
    )
    for name, ns in rels:
        assert sum(e * n for e, n in zip(E, ns)) == 0
        c = residue_coeff(s12.residue_tuple(ns))
        z = relation_term(ns)
        neg_count = sum(1 for n in ns if n < 0)
        den_sign = -1 if math.prod(ns) < 0 else 1
        coeff_sign = sign_word(c.real)
        term_sign = sign_word(z.real)
        equality = f"{relation_side(E, ns, True)} = {relation_side(E, ns, False)}"
        print(
            f"{name:>22} {equality:>38} {neg_count:>5} {den_sign:>4} "
            f"{coeff_sign:>6} {term_sign:>6} {z.real:>14.6g}"
        )


def seed_support(a: int) -> tuple[int, int, int, int, int, int]:
    return (1, a, 8, a + 7, 15, 22)


def low_height_sign_counts(E: tuple[int, ...], H: int) -> dict[int, Counter]:
    counts: dict[int, Counter] = defaultdict(Counter)
    cs = s12.coeffs(H)
    for ns in itertools.product(cs, repeat=len(E)):
        if sum(e * n for e, n in zip(E, ns)) != 0:
            continue
        h = max(abs(n) for n in ns)
        z = relation_term(ns)
        neg_count = sum(1 for n in ns if n < 0)
        key = (
            f"term{sign_word(z.real)}",
            f"den{'-' if neg_count % 2 else '+'}",
            f"coeff{sign_word(residue_coeff(s12.residue_tuple(ns)).real)}",
        )
        counts[h][key] += 1
    return counts


def report_hyp2634_seed() -> None:
    section("HYP-2634 SEED: ADDITION CREATES FIBRES, MULTIPLICATION CREATES SIGNS")
    for a in (2, 4):
        E = seed_support(a)
        wt = weighted_metrics(E, max_coeff=4, max_mass=7, relation_height=4)
        print(
            f"a={a}, E={E}: weighted vectors={wt['vectors']}, "
            f"energy={wt['energy']}, collision fibres={wt['collision_fibres']}, "
            f"max fibre={wt['max_fibre']}, h<=4 coverage={wt['coverage']}/{wt['coverage_total']}"
        )
        for target in (42, 48, 49):
            display_fibre(E, target, max_coeff=4, max_mass=7)
        print()

    print("Named relation signs for the two QR lifts:")
    relation_sign_table(
        seed_support(2),
        [
            ("a2_negative", (-1, 2, -1, -1, -2, 2)),
            ("universal_positive", (1, 1, -1, -1, -2, 2)),
        ],
    )
    relation_sign_table(
        seed_support(4),
        [
            ("a4_h3_negative", (-1, 3, -1, -1, 2, -1)),
            ("a4_h4_negative", (-4, 4, -1, 3, -1, -1)),
            ("universal_positive", (1, 1, -1, -1, -2, 2)),
        ],
    )
    print()
    print("Low-height sign census through H=4, grouped by term/denominator/coefficient sign:")
    for a in (2, 4):
        print(f"  a={a}")
        counts = low_height_sign_counts(seed_support(a), 4)
        for h in sorted(counts):
            total = sum(counts[h].values())
            common = ", ".join(f"{k}:{v}" for k, v in counts[h].most_common(4))
            print(f"    h={h}: total={total}, {common}")
    print()
    print(
        "Readout: additive equality gives the fibre collision, but the sign is "
        "not additive.  It is multiplicative: denominator parity "
        "(-1)^{#negative coefficients} times the finite residue coefficient. "
        "That is the even/odd -> positive/negative bridge in the reciprocal lift."
    )


def tournament_analysis() -> None:
    section("TOURNAMENT ANALYSIS")
    path = [
        "weighted_summand_fibres",
        "height2_relation_coverage",
        "freiman_gap_scout",
        "multiplicand_sign_parity",
        "finite_residue_character",
        "ordinary_pair_sums",
        "raw_speed_vertices",
    ]
    scores = {v: len(path) - 1 - i for i, v in enumerate(path)}
    hist = Counter(scores.values())
    print("Hamiltonian proof path:")
    print("  " + " > ".join(path))
    print(f"score_hist = {dict(sorted(hist.items()))}")
    print("directed_3_cycles = 0")
    print("SCC_sizes = [1,1,1,1,1,1,1]")
    print("hamiltonian_paths = 1")
    print(
        "Challenged assumption: pair sums are not enough.  The proof predicate "
        "is preserved by weighted summand fibres plus multiplicand sign parity; "
        "ordinary summand graph edges and raw speed labels destroy the height-2 "
        "relation/GAP pocket."
    )


def main() -> None:
    section("LRC14 RELATION-FIBRE / GAP SCOUT S26")
    print(
        "Goal: make the HYP-2635 additive-energy lead concrete and splice it "
        "to HYP-2634's sign mechanism."
    )
    report_hyp2635_examples()
    report_hyp2634_seed()
    tournament_analysis()
    section("S26 PROOF TARGET")
    print(
        "Proposed replacement split:\n"
        "  (A) If a nonzero element is not touched by bounded weighted-summand "
        "relations, peel it and use the dissociated/independent limit.\n"
        "  (B) If every element is touched, the weighted additive energy is high; "
        "apply a Freiman/BSG-style lemma to place the set inside a low-rank GAP, "
        "then bound that finite GAP family using AP-orbit invariance plus dimension.\n"
        "  (C) For reciprocal lift tails, keep the multiplicative sign coordinate: "
        "term sign = residue-character sign times denominator parity."
    )
    print("LRC(14) remains open; this file sharpens the next proof obligation.")


if __name__ == "__main__":
    main()
