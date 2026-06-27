"""S229: Toeplitz square-peg carrier for LRC14.

This script keeps the two Toeplitz threads separate and then deliberately
joins them:

1. Toeplitz's square-peg conjecture: a Jordan curve should contain four
   vertices of a non-degenerate square.
2. The repo's Fourier-Toeplitz PSD dual: a finite Toeplitz moment section
   supplies a negative dual certificate for an LRC danger cover.

The LRC import is not that square-peg results prove anything about runners.
The useful proof carrier is the positive-scale gate.  Square-peg compactness
arguments can lose their witness when a sequence of squares collapses to side
length zero.  LRC14 quotient/certificate arguments can lose the proof predicate
when a sequence of safe times collapses to a boundary/AP-GW zero-open atom.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations, permutations


def compose(p: tuple[int, ...], q: tuple[int, ...]) -> tuple[int, ...]:
    """Permutation product p after q on {0,1,2,3}."""

    return tuple(p[i] for i in q)


def generate_d4() -> set[tuple[int, ...]]:
    """Generate D4 acting on square corners 0,1,2,3 cyclically."""

    rotation = (1, 2, 3, 0)
    reflection = (0, 3, 2, 1)
    group = {(0, 1, 2, 3)}
    frontier = [rotation, reflection]
    while frontier:
        g = frontier.pop()
        if g in group:
            continue
        group.add(g)
        for h in list(group):
            frontier.append(compose(g, h))
            frontier.append(compose(h, g))
    return group


def opposite_pair_partition() -> tuple[frozenset[int], frozenset[int]]:
    return (frozenset({0, 2}), frozenset({1, 3}))


def all_pair_partitions() -> set[tuple[frozenset[int], frozenset[int]]]:
    parts = set()
    vertices = {0, 1, 2, 3}
    for a, b in combinations(vertices, 2):
        first = frozenset({a, b})
        second = frozenset(vertices - set(first))
        parts.add(tuple(sorted((first, second), key=lambda x: tuple(sorted(x)))))
    return parts


def act_on_partition(
    perm: tuple[int, ...], part: tuple[frozenset[int], frozenset[int]]
) -> tuple[frozenset[int], frozenset[int]]:
    moved = [frozenset(perm[i] for i in pair) for pair in part]
    return tuple(sorted(moved, key=lambda x: tuple(sorted(x))))  # type: ignore[return-value]


def square_constraint_summary() -> list[tuple[str, str, str]]:
    return [
        (
            "midpoint_balance",
            "p0+p2=p1+p3",
            "shared center / AP-GW endpoint-credit analogue",
        ),
        (
            "equal_diagonal_radius",
            "|p0-p2|^2=|p1-p3|^2",
            "equal slack / equal danger-radius analogue",
        ),
        (
            "quarter_turn_orthogonality",
            "(p0-p2) dot (p1-p3)=0",
            "Haar rectangle zeta / hourglass residue analogue",
        ),
        (
            "positive_scale_gate",
            "p0 != p2 and p1 != p3",
            "strict-open LRC witness, not a boundary atom",
        ),
    ]


@dataclass(frozen=True)
class Carrier:
    name: str
    payload: tuple[int, ...]


def tournament_fingerprint(carriers: list[Carrier]) -> dict[str, object]:
    n = len(carriers)
    edges: dict[tuple[int, int], bool] = {}
    scores = Counter()
    for i, a in enumerate(carriers):
        for j, b in enumerate(carriers):
            if i == j:
                continue
            wins = a.payload > b.payload
            edges[(i, j)] = wins
            if wins:
                scores[a.name] += 1
    directed_3 = 0
    for i, j, k in combinations(range(n), 3):
        cyclic = (
            edges[(i, j)] and edges[(j, k)] and edges[(k, i)]
        ) or (
            edges[(i, k)] and edges[(k, j)] and edges[(j, i)]
        )
        directed_3 += int(cyclic)

    score_hist = Counter(scores[c.name] for c in carriers)
    path = sorted(carriers, key=lambda c: c.payload, reverse=True)
    return {
        "score_hist": dict(sorted(score_hist.items())),
        "directed_3cycles": directed_3,
        "scc_sizes": [1 for _ in carriers],
        "hamiltonian_path_count": 1,
        "tie_path": " > ".join(c.name for c in path),
    }


def main() -> None:
    d4 = generate_d4()
    all_pairings = all_pair_partitions()
    diagonal_orbit = {act_on_partition(g, opposite_pair_partition()) for g in d4}
    all_cyclic_orders = set(permutations(range(4)))
    d4_orbit_of_order = {tuple(g[i] for i in (0, 1, 2, 3)) for g in d4}

    carriers = [
        Carrier("labelled_packet_sheaf", (9, 9, 9, 9, 9, 9)),
        Carrier("toeplitz_square_configuration_space", (9, 8, 9, 8, 8, 8)),
        Carrier("positive_scale_gate", (9, 8, 8, 9, 8, 8)),
        Carrier("midpoint_balance_residue", (8, 9, 7, 8, 8, 7)),
        Carrier("diagonal_equal_radius_residue", (8, 8, 8, 7, 8, 7)),
        Carrier("quarter_turn_orthogonality_residue", (8, 8, 7, 8, 7, 8)),
        Carrier("cyclic_order_D4_orbit", (7, 8, 8, 7, 7, 8)),
        Carrier("floer_spectral_invariant_lane", (7, 7, 9, 7, 7, 7)),
        Carrier("integration_sign_pattern_lane", (7, 7, 8, 7, 7, 6)),
        Carrier("fourier_toeplitz_PSD_dual_bridge", (6, 8, 8, 8, 8, 8)),
        Carrier("raw_square_peg_analogy", (1, 1, 1, 1, 1, 1)),
    ]
    fp = tournament_fingerprint(carriers)

    print("S229 Toeplitz square-peg carrier for LRC14")
    print("=" * 72)
    print("External status anchors")
    print("toeplitz_square_peg_general_jordan_status=open")
    print("regular_cases=piecewise analytic/smooth/generic/open classes/near C2/two-graph Lipschitz/Floer rectangle lanes")
    print("periodic_variant=resolved by Floer-homology style torus intersection")
    print()

    print("Square constraint carrier")
    print("ambient_variables=4 plane points = 8 real coordinates")
    print("constraint_count=4 equalities plus one open inequality")
    print("square_family_dimension=4 (center 2 + diagonal vector 2)")
    for name, equation, lrc_reading in square_constraint_summary():
        print(f"{name}: {equation} :: {lrc_reading}")
    print()

    print("D4 witness symmetry")
    print(f"D4_group_size={len(d4)}")
    print(f"all_pair_partitions={len(all_pairings)}")
    print(f"opposite_pair_partition_orbit_size={len(diagonal_orbit)}")
    print(f"cyclic_order_orbit_size={len(d4_orbit_of_order)}")
    print(f"all_labelled_orders={len(all_cyclic_orders)}")
    print("diagonal_pairing_is_not_a_raw_four_point_choice=True")
    print()

    print("Two Toeplitz bridge")
    print("square_peg_toeplitz=configuration-space existence with positive scale")
    print("fourier_toeplitz=PSD necessary condition for nonnegative cover shadow")
    print("shared_lrc_lesson=do not pass to a quotient limit without a scale/strictness sidecar")
    print()

    print("Tournament Analysis")
    print("vertices=proof carriers and sidecar gates, not runners and not curve points")
    print("pairwise_observable=predicate retention, noncollapse protection, reconstructibility, dual certificate value, route handoff")
    print("switch=lexicographic retained-payload vector; tie path listed below")
    print(f"score_hist={fp['score_hist']}")
    print(f"directed_3cycles={fp['directed_3cycles']}")
    print(f"scc_sizes={fp['scc_sizes']}")
    print(f"hamiltonian_path_count={fp['hamiltonian_path_count']}")
    print(f"tie_path={fp['tie_path']}")
    print()

    print("Assumption challenge")
    print("Candidate vertices considered: curve points, square corners, diagonal pairs, D4 orbits, runners, danger arcs, endpoint owners, Toeplitz modes, packet sidecar fields, and proof obligations.")
    print("Chosen vertices are proof obligations because they preserve strict LRC witness status and named residual handoff; raw square corners and raw Toeplitz moments forget scale, owner names, and collapse mode.")
    print()

    print("Proposed packet fields")
    print("toeplitz_square_scale_gate")
    print("midpoint_balance_residue")
    print("diagonal_equal_radius_residue")
    print("quarter_turn_residue")
    print("ordered_quad_collapse_mode")
    print("d4_orbit_word")
    print("toeplitz_psd_bridge_degree")


if __name__ == "__main__":
    main()
