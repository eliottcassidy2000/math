#!/usr/bin/env python3
"""Boundary-lift analogy atlas across irreducibility, LRC, unit distance, and towers.

codex-2026-06-13

Recent agents have independently landed several instances of the same pattern:

    visible scalar / boundary total  ->  hidden lift / support certificate

This script records a small, reproducible atlas and a Tournament Analysis over
candidate transfer routes.  It is not a proof of any frontier theorem; it is a
hypothesis generator meant to prevent the next session from collapsing distinct
fibers into a single numerology.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations


@dataclass(frozen=True)
class Carrier:
    name: str
    scalar_shadow: str
    hidden_lift: str
    irreducible_face: str
    local_anchor: str
    live_gap: str
    scores: dict[str, int]


CRITERIA = (
    "exact_certificate",
    "hidden_lift_visibility",
    "irreducibility_signal",
    "lrc14_transfer",
    "unit_distance_transfer",
    "code72_transfer",
    "computable_now",
    "scalar_warning",
)


CARRIERS = (
    Carrier(
        name="polynomial_convolution_lift",
        scalar_shadow="coefficient row a_k",
        hidden_lift="integer factor grid with diagonal sums a_k=sum b_i c_j",
        irreducible_face="irreducible iff no nontrivial integral convolution lift",
        local_anchor="HYP-2452: degree-4/5 primitive scans have zero mismatches",
        live_gap="bounded ILP/SAT lifts beyond degree 5 plus Newton-boundary layers",
        scores={
            "exact_certificate": 5,
            "hidden_lift_visibility": 5,
            "irreducibility_signal": 5,
            "lrc14_transfer": 4,
            "unit_distance_transfer": 2,
            "code72_transfer": 5,
            "computable_now": 5,
            "scalar_warning": 4,
        },
    ),
    Carrier(
        name="lrc_q27_support_lift",
        scalar_shadow="denominator q is blocked / plain shell horizon",
        hidden_lift="runner support, Pisano class, divisor fiber, Bprime owner",
        irreducible_face="hard row survives only if no AP/Vstar/2AP/descent support lift appears",
        local_anchor="HYP-2444: Q27 covers all one-stranger rows; residual is +-10 mod 27 plus 13-clock",
        live_gap="multi-stranger rows keeping low-clock coverage while spending shell-27 resources",
        scores={
            "exact_certificate": 4,
            "hidden_lift_visibility": 5,
            "irreducibility_signal": 4,
            "lrc14_transfer": 5,
            "unit_distance_transfer": 3,
            "code72_transfer": 4,
            "computable_now": 5,
            "scalar_warning": 5,
        },
    ),
    Carrier(
        name="unit_distance_moser_irreducible_fiber",
        scalar_shadow="edge count u(N) versus 3N and product average degree",
        hidden_lift="rank-4 Moser-ring incidence, direction packets, ears, exact embedding",
        irreducible_face="N* crosser is non-product; products tie at 27 and first beat only at 32",
        local_anchor="OPEN-Q-057/THM-433/THM-437/HYP-2301: H(3,3) rigid at 81; Moser row gives 85 at 28",
        live_gap="prove u(27)=81 or find a genuinely irreducible non-product 82-edge row",
        scores={
            "exact_certificate": 4,
            "hidden_lift_visibility": 4,
            "irreducibility_signal": 5,
            "lrc14_transfer": 3,
            "unit_distance_transfer": 5,
            "code72_transfer": 2,
            "computable_now": 4,
            "scalar_warning": 5,
        },
    ),
    Carrier(
        name="triangular_moment_fractional_address",
        scalar_shadow="equal interval sums or equal square sums",
        hidden_lift="moment layer, endpoint Pell address, fractional Bernoulli offset",
        irreducible_face="integer power centers appear only for p=1,2 in the checked family",
        local_anchor="HYP-2453/HYP-2454: unique Q_L(3)=F_R(4), 78/90 support shadow",
        live_gap="prove p>=3 bracket and attach 78/90 to an actual support ledger",
        scores={
            "exact_certificate": 4,
            "hidden_lift_visibility": 4,
            "irreducibility_signal": 3,
            "lrc14_transfer": 4,
            "unit_distance_transfer": 2,
            "code72_transfer": 4,
            "computable_now": 5,
            "scalar_warning": 4,
        },
    ),
    Carrier(
        name="code72_support_design_lift",
        scalar_shadow="Type II weight enumerator feasibility and lambda_5=78",
        hidden_lift="binary support/design/matroid incidence for weight-16 words",
        irreducible_face="a feasible enumerator is reducible noise unless a support incidence lift exists",
        local_anchor="HYP-2452/HYP-2454/HYP-2445 transfer: 78 appears as lambda_5 and D7 index",
        live_gap="turn 78/90 and convolution-lift grammar into support-incidence constraints",
        scores={
            "exact_certificate": 2,
            "hidden_lift_visibility": 5,
            "irreducibility_signal": 4,
            "lrc14_transfer": 4,
            "unit_distance_transfer": 2,
            "code72_transfer": 5,
            "computable_now": 3,
            "scalar_warning": 5,
        },
    ),
    Carrier(
        name="p_curvature_operator_ledger",
        scalar_shadow="naive mod-p scalar rank / local residue pass",
        hidden_lift="operator p-curvature plus carry/local-section compatibility",
        irreducible_face="local scalar tests are reducible shadows unless operator sections cohere",
        local_anchor="HYP-2446: 1/(1-z) and z/(1-z) give opposite scalar/operator false lessons",
        live_gap="define K_q(S) operator-style for LRC14 blockers and descent failures",
        scores={
            "exact_certificate": 4,
            "hidden_lift_visibility": 5,
            "irreducibility_signal": 3,
            "lrc14_transfer": 5,
            "unit_distance_transfer": 1,
            "code72_transfer": 3,
            "computable_now": 4,
            "scalar_warning": 5,
        },
    ),
    Carrier(
        name="product_quotient_diagonal_support_gate",
        scalar_shadow="Shioda supersingularity / group index arithmetic",
        hidden_lift="diagonal symmetric forms across asymmetric Frobenius twists",
        irreducible_face="scalar supersingularity is too reducible; diagonal support forces descent",
        local_anchor="HYP-2445: 91=C(14,2) and 78=C(13,2) match LRC/code support gates",
        live_gap="translate projection-degree descent into finite support-gate constraints",
        scores={
            "exact_certificate": 4,
            "hidden_lift_visibility": 5,
            "irreducibility_signal": 4,
            "lrc14_transfer": 4,
            "unit_distance_transfer": 2,
            "code72_transfer": 5,
            "computable_now": 3,
            "scalar_warning": 5,
        },
    ),
)


def beats(a: Carrier, b: Carrier) -> bool:
    """Majority tournament over criteria, with declaration order as tie path."""
    votes_a = 0
    votes_b = 0
    for criterion in CRITERIA:
        av = a.scores[criterion]
        bv = b.scores[criterion]
        if av > bv:
            votes_a += 1
        elif bv > av:
            votes_b += 1
    if votes_a != votes_b:
        return votes_a > votes_b
    return CARRIERS.index(a) < CARRIERS.index(b)


def scc_sizes(edges: dict[str, set[str]]) -> list[int]:
    names = [carrier.name for carrier in CARRIERS]
    rev: dict[str, set[str]] = {name: set() for name in names}
    for u, outs in edges.items():
        for v in outs:
            rev[v].add(u)

    seen: set[str] = set()
    order: list[str] = []

    def dfs(u: str) -> None:
        seen.add(u)
        for v in edges[u]:
            if v not in seen:
                dfs(v)
        order.append(u)

    for name in names:
        if name not in seen:
            dfs(name)

    seen.clear()
    sizes: list[int] = []

    def rdfs(u: str) -> int:
        seen.add(u)
        total = 1
        for v in rev[u]:
            if v not in seen:
                total += rdfs(v)
        return total

    for name in reversed(order):
        if name not in seen:
            sizes.append(rdfs(name))
    return sorted(sizes, reverse=True)


def hamiltonian_paths(edges: dict[str, set[str]]) -> int:
    names = [carrier.name for carrier in CARRIERS]
    n = len(names)
    index = {name: i for i, name in enumerate(names)}
    dp: dict[tuple[int, int], int] = {(1 << i, i): 1 for i in range(n)}
    for mask in range(1 << n):
        for last in range(n):
            count = dp.get((mask, last), 0)
            if not count:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if names[nxt] in edges[names[last]]:
                    dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + count
    full = (1 << n) - 1
    return sum(dp.get((full, i), 0) for i in range(n))


def directed_triangles(edges: dict[str, set[str]]) -> int:
    total = 0
    names = [carrier.name for carrier in CARRIERS]
    for a, b, c in combinations(names, 3):
        wins = Counter()
        for u, v in ((a, b), (a, c), (b, c)):
            wins[u if v in edges[u] else v] += 1
        if sorted(wins.values()) == [1, 1, 1]:
            total += 1
    return total


def tournament_fingerprint() -> tuple[dict[str, set[str]], dict[str, int]]:
    names = [carrier.name for carrier in CARRIERS]
    edges = {name: set() for name in names}
    for a, b in combinations(CARRIERS, 2):
        if beats(a, b):
            edges[a.name].add(b.name)
        else:
            edges[b.name].add(a.name)
    scores = {name: len(outs) for name, outs in edges.items()}
    return edges, scores


def print_carrier_atlas() -> None:
    print("Boundary-lift carrier atlas")
    print("===========================")
    print("Thesis: recent agents are seeing the same proof grammar in different clothes.")
    print("A scalar shadow can pass while the hidden lift fails, or vice versa.")
    print()
    for carrier in CARRIERS:
        print(f"* {carrier.name}")
        print(f"  scalar:      {carrier.scalar_shadow}")
        print(f"  hidden lift: {carrier.hidden_lift}")
        print(f"  irreducible: {carrier.irreducible_face}")
        print(f"  anchor:      {carrier.local_anchor}")
        print(f"  live gap:    {carrier.live_gap}")
        print()


def print_cross_analogies() -> None:
    rows = (
        (
            "Polynomial reducibility",
            "nontrivial factor grid exists",
            "absence of every integral convolution lift",
            "LRC/code support ledgers should be solved as lift feasibility, not scalar counts",
        ),
        (
            "Unit-distance product family",
            "Cartesian/Minkowski product attains edge count",
            "N* crosser below 32 must be non-product Moser fiber",
            "treat product constructions as reducible rows; search the irreducible fiber",
        ),
        (
            "LRC plain q-blocking",
            "scalar denominator appears blocked",
            "Q27/Pisano/owner lift opens or certifies the row",
            "multi-stranger frontier should track support allocations like factor tokens",
        ),
        (
            "Triangular equalities",
            "p=1 and p=2 interval sums balance",
            "p>=3 needs fractional/moment address in checked range",
            "78/90 is a support address, not a standalone numerological proof",
        ),
        (
            "[72,36,16] code",
            "Gleason/enumerator totals pass",
            "binary support-design incidence lift remains unresolved",
            "adapt convolution/UD non-product tests to rule in/out support realizability",
        ),
    )
    print("Irreducible analogy table")
    print("=========================")
    print("domain | reducible shadow | irreducible/lift face | transfer instruction")
    print("--- | --- | --- | ---")
    for row in rows:
        print(" | ".join(row))
    print()


def print_tournament() -> None:
    edges, scores = tournament_fingerprint()
    score_hist = Counter(scores.values())
    ranking = sorted(CARRIERS, key=lambda c: (-scores[c.name], c.name))
    print("Tournament Analysis")
    print("===================")
    print("Vertices: carriers/proof obligations, not runners, points, primes, or polynomial rows.")
    print("Pairwise observable: majority comparison across criteria:")
    print("  " + ", ".join(CRITERIA))
    print("Switch/gauge: A -> B when A wins more criteria; tie Hamiltonian path is declaration order.")
    print(f"score_hist={dict(sorted(score_hist.items()))}")
    print(f"directed_3cycles={directed_triangles(edges)}")
    print(f"scc_sizes={scc_sizes(edges)}")
    print(f"hamiltonian_paths={hamiltonian_paths(edges)}")
    print("ranking:")
    for carrier in ranking:
        print(f"  score={scores[carrier.name]} {carrier.name}")
    print()
    print("Nontrivial edge flips against scalar-warning-only order:")
    scalar_order = sorted(CARRIERS, key=lambda c: (-c.scores["scalar_warning"], c.name))
    rank = {carrier.name: i for i, carrier in enumerate(scalar_order)}
    flips = []
    for a, b in combinations(CARRIERS, 2):
        scalar_pref = a.name if rank[a.name] < rank[b.name] else b.name
        tour_pref = a.name if b.name in edges[a.name] else b.name
        if scalar_pref != tour_pref:
            flips.append((scalar_pref, tour_pref))
    print(f"edge_flips_vs_scalar_warning={len(flips)}/{len(CARRIERS) * (len(CARRIERS) - 1) // 2}")
    for scalar_pref, tour_pref in flips[:8]:
        print(f"  scalar would choose {scalar_pref}; majority carrier chooses {tour_pref}")
    print()


def print_next_moves() -> None:
    print("Procedural next moves")
    print("=====================")
    print("1. Build one common 'lift feasibility' data model:")
    print("   boundary totals, candidate hidden cells, local gates, surviving allocations, and proof owner.")
    print("2. Polynomial side: extend HYP-2452 with bounded ILP/SAT degree-6 lifts.")
    print("3. LRC side: turn Q27/Pisano one-stranger rows into multi-stranger allocation ledgers.")
    print("4. Unit-distance side: classify product-reducible versus Moser-irreducible extension fibers at N=27/28.")
    print("5. Code side: encode [72,36,16] minimum-word supports as binary incidence lifts over the 78/90 address.")
    print()
    print("Assumption challenge")
    print("====================")
    print("Alternate vertices considered: polynomial rows, factor coefficients, prime tokens,")
    print("LRC runners, residues, Q27 divisors, unit-distance points, unit directions,")
    print("Moser ears, triangular endpoints, code supports, Fourier modes, and proof obligations.")
    print("Chosen vertices are carriers/proof obligations because they preserve the predicate")
    print("'does a hidden lift exist?' across all four domains.  The quotient destroys")
    print("domain-specific geometry, so each transfer must reattach its local side channels")
    print("before being treated as evidence.")


def main() -> None:
    print_carrier_atlas()
    print_cross_analogies()
    print_tournament()
    print_next_moves()


if __name__ == "__main__":
    main()
