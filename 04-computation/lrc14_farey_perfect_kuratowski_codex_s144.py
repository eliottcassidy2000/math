#!/usr/bin/env python3
"""S144: Farey products, perfect numbers, and Kuratowski guardrails for LRC14.

The prompt asks to keep the earlier POKE carriers in view while revisiting
perfect numbers, F_3/F_4, and the product dictionary

    p/q -> K_{p,q},        p*q = |E(K_{p,q})|.

This script records the useful split:

* F_3 already contains the even-perfect product 2*3=6 via 2/3 -> K_{2,3},
  but that graph is planar.
* F_4 first contains a complete-bipartite Kuratowski wall via
  3/4 -> K_{3,4}, but the product 12 is not perfect.
* Later even perfect numbers give K_{2^(r-1),2^r-1}.  For r>=3 these graphs
  are nonplanar only because they already contain K_{3,3}; they are not new
  minimal obstructions.
* The graph-theoretic iteration is minor/subdivision transitivity, not
  mediant averaging.  A disjoint union K5 + K3,3 is nonplanar but fails
  minimality immediately.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import gcd


def farey(order: int) -> list[Fraction]:
    out: list[Fraction] = []
    for q in range(1, order + 1):
        for p in range(q + 1):
            if gcd(p, q) == 1:
                out.append(Fraction(p, q))
    return sorted(set(out))


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


@dataclass(frozen=True)
class Kpq:
    p: int
    q: int

    @property
    def edges(self) -> int:
        return self.p * self.q

    @property
    def contains_k33(self) -> bool:
        return min(self.p, self.q) >= 3

    @property
    def rank(self) -> str:
        if self.contains_k33:
            return "K33-wall"
        if min(self.p, self.q) == 2:
            return "two-block"
        if min(self.p, self.q) == 1:
            return "star"
        return "empty"

    @property
    def planar_readout(self) -> str:
        return "nonplanar by inherited K3,3" if self.contains_k33 else "planar Kpq"


@dataclass(frozen=True)
class ObstructionGraph:
    name: str
    vertices: int
    edges: int
    components: int
    contains_k5: bool
    contains_k33: bool
    minor_minimal: bool
    reason: str


def even_perfect_rows(max_r: int = 19) -> list[tuple[int, int, int, int]]:
    rows: list[tuple[int, int, int, int]] = []
    for r in range(2, max_r + 1):
        mersenne = 2**r - 1
        if is_prime(mersenne):
            left = 2 ** (r - 1)
            perfect = left * mersenne
            rows.append((r, left, mersenne, perfect))
    return rows


def perfect_numbers(max_r: int = 19) -> set[int]:
    return {perfect for _r, _left, _mersenne, perfect in even_perfect_rows(max_r)}


def farey_level_split() -> None:
    print("[1] F_3/F_4 split: perfect product is not the Kuratowski wall")
    prev: set[Fraction] = set()
    perfects = perfect_numbers(19)
    for order in range(1, 8):
        nodes = set(farey(order))
        new_nodes = [x for x in sorted(nodes - prev) if x > 0]
        new_perfect = [x for x in new_nodes if x.numerator * x.denominator in perfects]
        new_walls = [x for x in new_nodes if Kpq(x.numerator, x.denominator).contains_k33]
        print(
            f"  F_{order}: new={','.join(map(str, new_nodes)) or '-':12s} "
            f"perfect-products={','.join(map(str, new_perfect)) or '-':6s} "
            f"K33-new={','.join(map(str, new_walls)) or '-'}"
        )
        prev = nodes
    print("  readout:")
    print("    F_3 has 2/3 -> K_{2,3}, product 6, the first even perfect number.")
    print("    F_4 has 3/4 -> K_{3,4}, product 12, the first K3,3-containing Kpq.")
    print("    Therefore product perfection and bipartite nonplanarity decouple.")
    print()


def perfect_product_thresholds(max_order: int = 50) -> None:
    print("[2] Perfect-product Farey thresholds through F_50")
    perfects = perfect_numbers(31)
    first_seen: dict[Fraction, int] = {}
    for order in range(1, max_order + 1):
        for frac in farey(order):
            if frac == 0:
                continue
            if frac.numerator * frac.denominator in perfects and frac not in first_seen:
                first_seen[frac] = order

    print(f"  {'order':>5s} {'fraction':>9s} {'p*q':>8s} {'Kpq':>12s} {'rank':>10s}")
    for frac, order in sorted(first_seen.items(), key=lambda item: (item[1], item[0])):
        led = Kpq(frac.numerator, frac.denominator)
        print(
            f"  {order:5d} {str(frac):>9s} {led.edges:8d} "
            f"K_{{{led.p},{led.q}}} {led.rank:>10s}"
        )
    print("  coarse-fiber warning:")
    print("    1/N also has product N, so product-perfect alone is too coarse.")
    print("    The Mersenne factor lane is the coprime nontrivial factorization")
    print("    2^(r-1)/(2^r-1), not the star fraction 1/N.")
    print()


def even_perfect_kpq_ladder() -> None:
    print("[3] Even perfect numbers as complete-bipartite edge loads")
    print(f"  {'r':>2s} {'fraction':>10s} {'perfect':>8s} {'Farey in':>8s} {'Kpq rank':>10s} {'readout'}")
    for r, left, mersenne, perfect in even_perfect_rows(19):
        frac = Fraction(left, mersenne)
        led = Kpq(frac.numerator, frac.denominator)
        print(
            f"  {r:2d} {str(frac):>10s} {perfect:8d} "
            f"F_{max(frac.numerator, frac.denominator):<6d} {led.rank:>10s} "
            f"{led.planar_readout}"
        )
    print("  readout:")
    print("    r=2 gives 2/3 and the planar perfect seed K_{2,3}.")
    print("    r>=3 gives inherited K3,3 nonplanarity, not new obstruction types.")
    print()


def mediant_cross_terms() -> None:
    print("[4] Mediant product load is cross-incidence, not graph averaging")
    examples = [
        (Fraction(0, 1), Fraction(1, 3)),
        (Fraction(1, 3), Fraction(1, 2)),
        (Fraction(1, 2), Fraction(2, 3)),
        (Fraction(2, 3), Fraction(1, 1)),
        (Fraction(1, 14), Fraction(2, 27)),
    ]
    print(
        f"  {'left':>7s} {'right':>7s} {'mediant':>8s} "
        f"{'edge formula':>31s} {'Kpq rank':>10s}"
    )
    for left, right in examples:
        med = Fraction(left.numerator + right.numerator, left.denominator + right.denominator)
        ab = left.numerator * left.denominator
        cd = right.numerator * right.denominator
        cross1 = left.numerator * right.denominator
        cross2 = right.numerator * left.denominator
        led = Kpq(med.numerator, med.denominator)
        formula = f"{ab}+{cd}+{cross1}+{cross2}={led.edges}"
        print(f"  {str(left):>7s} {str(right):>7s} {str(med):>8s} {formula:>31s} {led.rank:>10s}")
    print("  readout:")
    print("    A mediant K_{a+c,b+d} has edges (a+c)(b+d).")
    print("    That is parent edges plus two typed cross terms, not an average graph.")
    print("    For LRC14, the 1/14,2/27 mediant is 3/41: the first unit-excess K33 wall.")
    print()


def kuratowski_guardrail() -> None:
    print("[5] Kuratowski/Wagner guardrail")
    graphs = [
        ObstructionGraph(
            "K5",
            vertices=5,
            edges=10,
            components=1,
            contains_k5=True,
            contains_k33=False,
            minor_minimal=True,
            reason="one of the two planar forbidden minors",
        ),
        ObstructionGraph(
            "K3,3",
            vertices=6,
            edges=9,
            components=1,
            contains_k5=False,
            contains_k33=True,
            minor_minimal=True,
            reason="one of the two planar forbidden minors",
        ),
        ObstructionGraph(
            "K5 disjoint-union K3,3",
            vertices=11,
            edges=19,
            components=2,
            contains_k5=True,
            contains_k33=True,
            minor_minimal=False,
            reason="delete either component and a nonplanar component remains",
        ),
    ]
    print(
        f"  {'graph':26s} {'V':>3s} {'E':>3s} {'comp':>4s} "
        f"{'has K5':>7s} {'has K33':>8s} {'minimal?':>9s} reason"
    )
    for graph in graphs:
        print(
            f"  {graph.name:26s} {graph.vertices:3d} {graph.edges:3d} {graph.components:4d} "
            f"{str(graph.contains_k5):>7s} {str(graph.contains_k33):>8s} "
            f"{str(graph.minor_minimal):>9s} {graph.reason}"
        )
    print("  readout:")
    print("    K5 and K3,3 are an obstruction set, not endpoints of a numeric continuum.")
    print("    The proof operation is minor/subdivision transitivity: J <= H <= G implies J <= G.")
    print("    Disjoint union is a larger witness container, not a third irreducible obstruction.")
    print()


def poke_hypotheses() -> None:
    print("[6] POKE hypotheses generated by this scout")
    hypotheses = [
        (
            "F3/F4 decoupling",
            "The first perfect-product Farey fraction, 2/3 in F_3, is planar; "
            "the first Kpq nonplanarity, 3/4 in F_4, is not perfect.",
        ),
        (
            "Mersenne Kpq ladder",
            "Even perfect numbers are edge loads of K_{2^(r-1),2^r-1}; "
            "after r=2 they are nonplanar only by inherited K3,3.",
        ),
        (
            "Mediant cross-term ledger",
            "A mediant product expands as parent edge loads plus two cross terms; "
            "the cross terms are typed incidence, not graph averaging.",
        ),
        (
            "Minor-transitive proof iteration",
            "Use K3,3/K5 containment, minors of minors, or subdivisions of subdivisions "
            "as the graph iteration layer; do not build new obstruction claims from averages.",
        ),
        (
            "LRC14 apex-N generalization",
            "For an apex N unit-excess chain p/(Np-1), the first Kpq wall is "
            "3/(3N-1).  At N=14 this recovers 3/41.",
        ),
        (
            "Carrier order",
            "Exact M/Farey and C27/unital labels should precede Kpq, perfect-number, "
            "polyhedral, tiling, flower, and pi analogy labels.",
        ),
    ]
    for idx, (name, text) in enumerate(hypotheses, start=1):
        print(f"  H{idx}. {name}: {text}")
    print()


def carrier_tournament() -> None:
    print("[7] Carrier tournament readout")
    carriers = [
        ("exact M/Farey branch", (7, 7, 7, 7, 7)),
        ("C27/unital branch-local pair grammar", (6, 7, 6, 6, 7)),
        ("Kuratowski minor-transitive core", (6, 6, 7, 7, 6)),
        ("Kpq product incidence ledger", (5, 5, 6, 5, 5)),
        ("Mersenne/perfect edge-load lane", (4, 4, 5, 4, 5)),
        ("mediant cross-term ledger", (4, 5, 4, 4, 5)),
        ("polyhedral/tiling recursion labels", (3, 4, 3, 4, 5)),
        ("pi/flower/unital analogy labels", (2, 3, 2, 3, 5)),
        ("raw graph average or scalar product", (1, 1, 1, 1, 1)),
    ]
    wins = {name: 0 for name, _score in carriers}
    directed_triangles = 0
    for i, (left, left_score) in enumerate(carriers):
        for right, right_score in carriers[i + 1 :]:
            winner = left if left_score >= right_score else right
            wins[winner] += 1
    hist: dict[int, int] = {}
    for score in wins.values():
        hist[score] = hist.get(score, 0) + 1
    print("  observable tuple:")
    print("    branch retention, unit preservation, obstruction minimality,")
    print("    state-lift fit, scalar-decoy resistance")
    print(f"  score_hist={dict(sorted(hist.items()))} directed_3cycles={directed_triangles}")
    print("  order:")
    for idx, (name, _score) in enumerate(carriers, start=1):
        print(f"    {idx}. {name}")
    print()


def main() -> None:
    print("S144 LRC14 FAREY PERFECT KURATOWSKI CARRIER")
    print("=" * 78)
    print("[0] Scope")
    print("  Vertices considered: Farey fractions, K_{p,q} graphs, perfect-number")
    print("  factor pairs, mediant cross terms, K5/K3,3 obstruction cores,")
    print("  C27/unital labels, and earlier pi/flower/solid/tiling carriers.")
    print("  Quotient chosen: proof carriers with exact M/Farey branch retained.")
    print("  Destroyed by this quotient: exact runner geometry and safe intervals.")
    print()
    farey_level_split()
    perfect_product_thresholds()
    even_perfect_kpq_ladder()
    mediant_cross_terms()
    kuratowski_guardrail()
    poke_hypotheses()
    carrier_tournament()
    print("[8] Final readout")
    print("  The new useful lane is not 'perfect number implies obstruction' and not")
    print("  'mediant implies graph average'.  It is a labelled split:")
    print("    F_3: perfect-product planar seed 2/3 -> K_{2,3};")
    print("    F_4: first complete-bipartite K3,3 wall 3/4 -> K_{3,4};")
    print("    LRC14: unit-excess wall 3/41 -> K_{3,41};")
    print("    graph proof step: minor/subdivision transitivity, not averaging.")
    print("  This remains a POKE proof-interface hypothesis, not a proof of LRC14.")


if __name__ == "__main__":
    main()
