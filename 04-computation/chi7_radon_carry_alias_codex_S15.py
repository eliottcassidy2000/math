#!/usr/bin/env python3
"""Exact Radon/carry audit for THM-850's depth-11 chi-seven plane.

The scalar affine Radon transform on F_7^2 is injective when all eight
directions are retained.  The 78 degree-11 address monomials do not inject
into that plane: reduction mod 7 has 33 singleton fibres, 15 triple fibres,
and one empty residue.  This script certifies the resulting rank-48 pushforward
and identifies the missing A/B/C carry role.

Tournament Analysis uses affine needle directions, not runners or tournament
arcs, as vertices.  The pairwise observable is line-load concentration.  The
gauge switch reverses whether localization or equidistribution is preferred;
the declared direction order is the tie Hamiltonian path.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from itertools import combinations, permutations
from pathlib import Path


P = 7
MODULUS = 101
R = 11
PLANE_SUM = R % P
DIRECTIONS = (
    (1, 0), (0, 1), (1, 1),
    (1, 2), (1, 4), (1, 6),
    (1, 3), (1, 5),
)


def modular_rank(rows: list[list[int]], modulus: int = MODULUS) -> int:
    """Rank of an integer matrix over F_modulus."""
    if not rows:
        return 0
    matrix = [[value % modulus for value in row] for row in rows]
    height, width = len(matrix), len(matrix[0])
    pivot_row = 0
    for column in range(width):
        pivot = next(
            (row for row in range(pivot_row, height) if matrix[row][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        inverse = pow(matrix[pivot_row][column], -1, modulus)
        matrix[pivot_row] = [value * inverse % modulus for value in matrix[pivot_row]]
        for row in range(height):
            if row == pivot_row or not matrix[row][column]:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                (left - factor * right) % modulus
                for left, right in zip(matrix[row], matrix[pivot_row], strict=True)
            ]
        pivot_row += 1
        if pivot_row == height:
            break
    return pivot_row


def normalize_direction(triple: tuple[int, int, int]) -> tuple[int, int]:
    """Pass from a coefficient triple modulo constants to P^1(F_7)."""
    left = (triple[0] - triple[2]) % P
    right = (triple[1] - triple[2]) % P
    assert left or right
    scale = pow(left if left else right, -1, P)
    return left * scale % P, right * scale % P


def direction_orbit(seed: tuple[int, int, int]) -> frozenset[tuple[int, int]]:
    return frozenset(
        normalize_direction(tuple(seed[index] for index in order))
        for order in permutations(range(3))
    )


def incidence_rows(
    directions: tuple[tuple[int, int], ...],
    domain: tuple[tuple[int, ...], ...],
) -> list[list[int]]:
    rows = []
    for left, right in directions:
        for level in range(P):
            rows.append(
                [int((left * point[0] + right * point[1]) % P == level) for point in domain]
            )
    return rows


def tournament_fingerprint(values: list[int], reverse: bool = False) -> dict[str, object]:
    """Transitive comparison tournament with input order as the tie path."""
    count = len(values)
    arcs: set[tuple[int, int]] = set()
    scores = [0] * count
    for left in range(count):
        for right in range(left + 1, count):
            lhs, rhs = values[left], values[right]
            if reverse:
                lhs, rhs = -lhs, -rhs
            winner, loser = (left, right) if lhs >= rhs else (right, left)
            arcs.add((winner, loser))
            scores[winner] += 1

    paths = [[0] * count for _ in range(1 << count)]
    for vertex in range(count):
        paths[1 << vertex][vertex] = 1
    for support in range(1 << count):
        for last, number in enumerate(paths[support]):
            if not number:
                continue
            for nxt in range(count):
                if not (support >> nxt) & 1 and (last, nxt) in arcs:
                    paths[support | (1 << nxt)][nxt] += number
    directed_triangles = sum(
        (a, b) in arcs and (b, c) in arcs and (c, a) in arcs
        for a, b, c in combinations(range(count), 3)
    )
    return {
        "arcs": arcs,
        "score_histogram": dict(sorted(Counter(scores).items())),
        "directed_triangles": directed_triangles,
        "scc_sizes": [1] * count if directed_triangles == 0 else None,
        "hamiltonian_paths": sum(paths[-1]),
    }


def run() -> str:
    addresses = tuple(
        (alpha, beta, R - alpha - beta)
        for alpha in range(R + 1)
        for beta in range(R - alpha + 1)
    )
    plane = tuple(
        (left, right, (PLANE_SUM - left - right) % P)
        for left in range(P)
        for right in range(P)
    )
    assert len(addresses) == 78 and len(plane) == 49

    buckets: dict[tuple[int, int, int], list[tuple[int, int, int]]] = defaultdict(list)
    for address in addresses:
        residue = tuple(value % P for value in address)
        buckets[residue].append(address)
    fibre_histogram = Counter(len(buckets.get(residue, ())) for residue in plane)
    assert fibre_histogram == {0: 1, 1: 33, 3: 15}
    assert [residue for residue in plane if residue not in buckets] == [(6, 6, 6)]

    carry_labels = []
    carry_roles = []
    for address in addresses:
        residue = tuple(value % P for value in address)
        quotient = tuple((value - base) // P for value, base in zip(address, residue, strict=True))
        level = sum(quotient)
        assert level in (0, 1)
        role = -1 if level == 0 else quotient.index(1)
        carry_labels.append((*residue, role))
        carry_roles.append(role)
    assert len(set(carry_labels)) == 78

    # A nonzero determinant modulo 101 certifies the same lower rank over Q;
    # the universal 1+6d row relation supplies the matching characteristic-zero upper bound.
    subset_rank_rows = []
    plane_xy = tuple(point[:2] for point in plane)
    address_xy = tuple(point[:2] for point in addresses)
    for size in range(1, len(DIRECTIONS) + 1):
        plane_ranks = Counter()
        address_ranks = Counter()
        for subset in combinations(DIRECTIONS, size):
            plane_rank = modular_rank(incidence_rows(subset, plane_xy))
            address_rank = modular_rank(incidence_rows(subset, address_xy))
            expected_plane = 1 + 6 * size
            expected_address = expected_plane if size < 8 else 48
            assert plane_rank == expected_plane
            assert address_rank == expected_address
            plane_ranks[plane_rank] += 1
            address_ranks[address_rank] += 1
        subset_rank_rows.append((size, plane_ranks, address_ranks))

    expected_profiles = {
        (1, 0): (17, 15, 13, 11, 9, 7, 6),
        (0, 1): (17, 15, 13, 11, 9, 7, 6),
        (1, 1): (9, 11, 13, 15, 17, 6, 7),
        (1, 2): (11, 11, 11, 11, 12, 11, 11),
        (1, 4): (11, 11, 12, 11, 11, 11, 11),
        (1, 6): (12, 11, 11, 11, 11, 11, 11),
        (1, 3): (11, 11, 11, 12, 11, 11, 11),
        (1, 5): (11, 12, 11, 11, 11, 11, 11),
    }
    profiles = {}
    for direction in DIRECTIONS:
        profile = tuple(
            sum(
                (direction[0] * alpha + direction[1] * beta) % P == level
                for alpha, beta, _ in addresses
            )
            for level in range(P)
        )
        profiles[direction] = profile
    assert profiles == expected_profiles

    coordinate_orbit = direction_orbit((1, 0, 0))
    mixed_orbit = direction_orbit((1, 2, 0))
    chi_orbit = direction_orbit((1, 3, 0))
    assert coordinate_orbit == frozenset({(1, 0), (0, 1), (1, 1)})
    assert mixed_orbit == frozenset({(1, 2), (1, 4), (1, 6)})
    assert chi_orbit == frozenset({(1, 3), (1, 5)})
    assert coordinate_orbit | mixed_orbit | chi_orbit == frozenset(DIRECTIONS)
    assert not (coordinate_orbit & mixed_orbit or coordinate_orbit & chi_orbit or mixed_orbit & chi_orbit)

    # q=(1,2,4): subtract the gamma coefficient and projectively normalize.
    assert normalize_direction((1, 2, 4)) == (1, 3)
    charge_profile = tuple(
        sum((alpha + 2 * beta + 4 * gamma) % P == level for alpha, beta, gamma in addresses)
        for level in range(P)
    )
    assert charge_profile == (12, 11, 11, 11, 11, 11, 11)

    # Once carry is retained as four channels, six common directions suffice
    # and are necessary. True endpoint reflection is weaker than formal S3.
    carry_domains = {
        role: tuple(addresses[index][:2] for index, value in enumerate(carry_roles) if value == role)
        for role in (-1, 0, 1, 2)
    }
    assert {role: len(domain) for role, domain in carry_domains.items()} == {
        -1: 33, 0: 15, 1: 15, 2: 15
    }
    reflection_six = frozenset({(1, 0), (0, 1), (1, 2), (1, 4), (1, 3), (1, 5)})
    triality_six = coordinate_orbit | mixed_orbit

    def channel_ranks(subset: frozenset[tuple[int, int]]) -> tuple[int, ...]:
        ordered = tuple(direction for direction in DIRECTIONS if direction in subset)
        return tuple(modular_rank(incidence_rows(ordered, carry_domains[role])) for role in (-1, 0, 1, 2))

    def reflect_direction(direction: tuple[int, int]) -> tuple[int, int]:
        return normalize_direction((direction[1], direction[0], 0))

    assert frozenset(map(reflect_direction, reflection_six)) == reflection_six
    assert channel_ranks(reflection_six) == (33, 15, 15, 15)
    assert sum(channel_ranks(reflection_six)) == 78
    assert channel_ranks(triality_six) == (32, 15, 15, 15)
    assert sum(channel_ranks(triality_six)) == 77

    five_max_no_carry = max(
        modular_rank(incidence_rows(subset, carry_domains[-1]))
        for subset in combinations(DIRECTIONS, 5)
    )
    assert five_max_no_carry <= 31 < 33
    invariant_full_six = []
    for subset_tuple in combinations(DIRECTIONS, 6):
        subset = frozenset(subset_tuple)
        if frozenset(map(reflect_direction, subset)) != subset:
            continue
        if channel_ranks(subset) == (33, 15, 15, 15):
            invariant_full_six.append(subset)
    assert invariant_full_six == [reflection_six]

    concentration = [
        sum((P * load - len(addresses)) ** 2 for load in profiles[direction])
        for direction in DIRECTIONS
    ]
    assert concentration[:3] == [4942, 4942, 4942]
    assert concentration[3:] == [42, 42, 42, 42, 42]
    localized = tournament_fingerprint(concentration)
    uniform = tournament_fingerprint(concentration, reverse=True)
    edge_flips = len(localized["arcs"].symmetric_difference(uniform["arcs"])) // 2
    assert edge_flips == 15
    assert localized["directed_triangles"] == uniform["directed_triangles"] == 0
    assert localized["hamiltonian_paths"] == uniform["hamiltonian_paths"] == 1

    lines = [
        "THM-850: CHI7 RADON/CARRY ALIAS EXACT AUDIT",
        "=" * 72,
        "depth address simplex: r=11, addresses=78",
        "residue plane H_4: cells=49, occupied=48",
        "residue fibres: empty=1, singleton=33, triple=15",
        "empty canonical residue: (6,6,6)",
        "pushforward rank/kernel over Q: 48/30 = 48/(15*(3-1))",
        "carry repair: singleton role=-1; triple roles=A/B/C; repaired labels=78",
        "",
        "RADON SUBSET RANKS (plane/address)",
    ]
    for size, plane_ranks, address_ranks in subset_rank_rows:
        lines.append(
            f"d={size}: plane={dict(sorted(plane_ranks.items()))}, "
            f"address={dict(sorted(address_ranks.items()))}"
        )
    lines.extend(
        [
            "",
            "S3 DIRECTION ORBITS",
            f"coordinate={sorted(coordinate_orbit)} profile multiset=(17,15,13,11,9,7,6)",
            f"mixed={sorted(mixed_orbit)} profile multiset=(12,11,11,11,11,11,11)",
            f"chi={sorted(chi_orbit)} profile multiset=(12,11,11,11,11,11,11)",
            f"q=(1,2,4) normalizes to (1,3); actual charge profile={charge_profile}",
            "",
            "CARRY-CHANNEL TOMOGRAPHY",
            f"channel sizes none/A/B/C={tuple(len(carry_domains[role]) for role in (-1,0,1,2))}",
            f"reflection-invariant six={sorted(reflection_six)} ranks={channel_ranks(reflection_six)} total=78",
            f"triality-invariant six={sorted(triality_six)} ranks={channel_ranks(triality_six)} total=77",
            f"maximum no-carry rank with five directions={five_max_no_carry}<33",
            "unique reflection-invariant full-rank six-subset: PASS",
            "",
            "TOURNAMENT ANALYSIS ON NEEDLE DIRECTIONS",
            "pairwise observable=line-load concentration; switch=localization vs uniformity",
            f"localized score/C3/SCC/HP={localized['score_histogram']}/0/{localized['scc_sizes']}/1",
            f"uniform score/C3/SCC/HP={uniform['score_histogram']}/0/{uniform['scc_sizes']}/1",
            f"edge flips between gauges={edge_flips}",
            "PRESERVES: residue-aggregated address values and all eight affine line pencils",
            "DESTROYS: 30 independent within-residue A/B/C carry contrasts",
            "CHALLENGED ASSUMPTION: needle directions, not runners or arcs, are the vertices",
            "ALL ASSERTIONS PASSED",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    if not __debug__:
        raise RuntimeError("Radon/carry verification requires assertions; do not use python -O")
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("05-knowledge/results/chi7_radon_carry_alias_codex_S15.out"),
    )
    args = parser.parse_args()
    text = run()
    args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
