#!/usr/bin/env python3
"""HYP-2705 / T940: projective metric and synchronizer scout.

This is a small exact-fraction bridge from the new external analogies back to
the current LRC14 true-wide proof route.

Objects:
* Fubini-Study lens: for the death-chain transition law, the squared
  projective overlap with a subspace is the probability mass in that subspace.
  We keep cos^2(theta) as a rational, avoiding numerical angles.
* Road-coloring lens: a missed-sector set is synchronized by a word when all
  of its colors have appeared.  The iid synchronization probability is the
  exact coupon-collector transition probability to depth 0.
* Gibbs lens: negative survival currency is the low-temperature obstruction
  state.  Two decorrelated hits remove all negative survival coefficients in
  the depth quotient; true rows only remain hard through signed resonant
  deviation away from this quotient.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from itertools import permutations
from math import comb


MIDDLE = {1, 2, 3, 4}


def fmt(q: F) -> str:
    return str(q)


def transition_prob(t: int, s: int, r: int) -> F:
    """Probability t missed sectors become s missed after r iid color hits."""

    if not 0 <= s <= t <= 6:
        return F(0)
    need_hit = t - s
    total = F(0)
    for j in range(need_hit + 1):
        total += ((-1) ** j) * comb(need_hit, j) * F(7 - s - j, 7) ** r
    return F(comb(t, s)) * total


def currency_coeff(t: int) -> F:
    """Coefficient in C=p1+p2+p3+p4-4p6."""

    if 1 <= t <= 4:
        return F(1)
    if t == 6:
        return F(-4)
    return F(0)


def expected_currency(t: int, r: int) -> F:
    return sum(transition_prob(t, s, r) * currency_coeff(s) for s in range(t + 1))


def middle_projection_cos2(t: int, r: int) -> F:
    """Squared FS overlap with the middle-survival subspace span{s=1..4}."""

    return sum(transition_prob(t, s, r) for s in MIDDLE if s <= t)


def sync_probability(t: int, r: int) -> F:
    """Probability all initially missed colors have appeared after r hits."""

    return transition_prob(t, 0, r)


def directed_3cycles(edges: set[tuple[int, int]], n: int) -> int:
    total = 0
    for a in range(n):
        for b in range(a + 1, n):
            for c in range(b + 1, n):
                if (
                    (a, b) in edges and (b, c) in edges and (c, a) in edges
                ) or (
                    (a, c) in edges and (c, b) in edges and (b, a) in edges
                ):
                    total += 1
    return total


def hamiltonian_path_count(edges: set[tuple[int, int]], n: int) -> int:
    count = 0
    for perm in permutations(range(n)):
        if all((perm[i], perm[i + 1]) in edges for i in range(n - 1)):
            count += 1
    return count


def print_fs_table() -> None:
    print("FUBINI-STUDY PROJECTION LEDGER")
    print(
        "  For |psi_{t,r}> with coordinates sqrt(Pr[t->s after r hits]),"
    )
    print(
        "  cos^2(theta_middle)=||Proj_{s=1..4} psi||^2 is the exact middle mass."
    )
    print("  Start from full missed state t=6:")
    print("       r   cos2_middle   sync_to_0   E[C after r hits]")
    for r in range(0, 9):
        print(
            f"      {r:>2}   {fmt(middle_projection_cos2(6, r)):>12}"
            f"   {fmt(sync_probability(6, r)):>9}"
            f"   {fmt(expected_currency(6, r)):>16}"
        )
    print()
    print("  Two-hit depth profile:")
    print("       t   cos2_middle   sync_to_0   E[C after 2 hits]")
    for t in range(7):
        print(
            f"      {t:>2}   {fmt(middle_projection_cos2(t, 2)):>12}"
            f"   {fmt(sync_probability(t, 2)):>9}"
            f"   {fmt(expected_currency(t, 2)):>16}"
        )
    print()


def print_synchronizer_table() -> None:
    print("ROAD-COLORING / SYNCHRONIZER LEDGER")
    print("  Known missed set of size t: shortest synchronizing word length is t.")
    print("  Unknown missed subset of the six inner colors: length 6 is necessary and sufficient.")
    print("  Iid hit synchronization probability from t=6:")
    print("       r   Pr(depth 0)        Pr(middle)")
    for r in range(0, 13):
        print(
            f"      {r:>2}   {fmt(sync_probability(6, r)):>11}"
            f"   {fmt(middle_projection_cos2(6, r)):>11}"
        )
    print()


def print_gibbs_obstruction() -> None:
    print("GIBBS / LOW-TEMPERATURE OBSTRUCTION LEDGER")
    print("  Energy proxy: negative survival currency depth after decorrelated hits.")
    for r in range(0, 5):
        values = [(t, expected_currency(t, r)) for t in range(7)]
        negative = [(t, v) for t, v in values if v < 0]
        min_t, min_v = min(values, key=lambda item: (item[1], -item[0]))
        if negative:
            neg_text = ", ".join(f"t={t}:{fmt(v)}" for t, v in negative)
        else:
            neg_text = "none"
        print(
            f"    r={r}: min=t{min_t}:{fmt(min_v)}, negative_depths={neg_text}"
        )
    print(
        "  Consequence: in the missed-count quotient, two decorrelated far hits"
    )
    print(
        "  eliminate negative currency.  The LRC14 difficulty is the signed"
    )
    print("  resonant deviation from this quotient, not the quotient itself.")
    print()


def print_tournament_analysis() -> None:
    print("TOURNAMENT ANALYSIS")
    vertices = [
        "fs_middle_projection",
        "road_coloring_sync",
        "gibbs_currency_gap",
        "hebbian_density_matrix",
        "propagator_kernel",
        "clifford_stabilizer_cutspace",
        "cat_map_markov_partition",
        "crossing_squarefree_carrier",
        "raw_analogy",
    ]
    # Rank by how directly the lens preserves the current LRC cap predicate.
    rank = {
        "fs_middle_projection": 0,
        "road_coloring_sync": 1,
        "gibbs_currency_gap": 2,
        "hebbian_density_matrix": 3,
        "propagator_kernel": 4,
        "clifford_stabilizer_cutspace": 5,
        "cat_map_markov_partition": 6,
        "crossing_squarefree_carrier": 7,
        "raw_analogy": 8,
    }
    edges: set[tuple[int, int]] = set()
    for i, a in enumerate(vertices):
        for j, b in enumerate(vertices):
            if i == j:
                continue
            if rank[a] < rank[b]:
                edges.add((i, j))
    scores = Counter()
    for i in range(len(vertices)):
        scores[sum(1 for edge in edges if edge[0] == i)] += 1
    print("  vertices: proof lenses, not runners or arcs")
    print("  pairwise observable: how much LRC cap/survival data the lens preserves")
    print("  switch/gauge: orient toward less scalarization before margin comparison")
    print(f"  Hamiltonian path: {' > '.join(vertices)}")
    print(f"  score_hist={dict(sorted(scores.items()))}")
    print(f"  directed_3cycles={directed_3cycles(edges, len(vertices))}")
    print(
        f"  hamiltonian_path_count={hamiltonian_path_count(edges, len(vertices))}"
    )
    print(
        "  challenged assumption: the new physics/neural analogies are not"
    )
    print(
        "  separate routes until they preserve the missed-sector cap predicate."
    )
    print()


def main() -> None:
    print("LRC14 PROJECTIVE/SYNCHRONIZER/GIBBS SCOUT (exact Fraction)")
    print("=" * 78)
    print_fs_table()
    print_synchronizer_table()
    print_gibbs_obstruction()
    print_tournament_analysis()


if __name__ == "__main__":
    main()
