#!/usr/bin/env python3
"""Exact THM-850 addendum: Fano-flag homology and chi7 Radon repair.

The 21 j=4 flood labels are K7 edges.  Sending {a,b} to the incident
Fano flag (a xor b, {a,b,a xor b}) identifies their missing point/Fano
marginal sector with the cycle space of the Heawood incidence graph.
Oriented Heawood hexagons then give the smallest natural curl probes.

On the depth-11 no-carry address channel, the script also identifies the
one triality-six Radon alias as an alternating chi7 Vandermonde character
and certifies the unique reflection-invariant full-rank six-direction deck.

Tournament Analysis uses eight oriented hexagon probes, not runners or arcs,
as vertices.  Its pairwise observable is absolute V1 cycle charge, the
larger-charge gauge orients comparisons, and basis order is the tie path.

Scope: this is an exact carrier-reconstruction theorem.  It neither defines
the missing flood-to-operation role map nor closes a THM-741 metric sweep.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product
from pathlib import Path


P = 7
R = 11
HERE = Path(__file__).resolve().parent
THM741_PATH = HERE / "lrc14_thm741_2002_body_j4_tree_kps_S128c5.py"

POINTS = tuple(range(1, 8))
EDGES = tuple(combinations(POINTS, 2))
EDGE_INDEX = {frozenset(edge): index for index, edge in enumerate(EDGES)}
LINES = tuple(
    sorted(
        {
            tuple(sorted((left, right, left ^ right)))
            for left, right in EDGES
        }
    )
)

DIRECTIONS = (
    (1, 0),
    (0, 1),
    (1, 1),
    (1, 2),
    (1, 4),
    (1, 6),
    (1, 3),
    (1, 5),
)
TRIALITY_SIX = (
    (1, 0),
    (0, 1),
    (1, 1),
    (1, 2),
    (1, 4),
    (1, 6),
)
REFLECTION_SIX = (
    (1, 0),
    (0, 1),
    (1, 2),
    (1, 4),
    (1, 3),
    (1, 5),
)

HEXAGON_BASIS = (
    (1, 2, 5),
    (1, 2, 6),
    (1, 2, 7),
    (1, 3, 4),
    (1, 3, 5),
    (1, 3, 6),
    (1, 3, 7),
    (1, 4, 6),
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rational_rank(rows: list[list[int | Fraction]]) -> int:
    """Exact row rank over Q."""
    if not rows:
        return 0
    matrix = [[Fraction(value) for value in row] for row in rows]
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
        scale = matrix[pivot_row][column]
        matrix[pivot_row] = [value / scale for value in matrix[pivot_row]]
        for row in range(height):
            if row == pivot_row or not matrix[row][column]:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                left - factor * right
                for left, right in zip(matrix[row], matrix[pivot_row], strict=True)
            ]
        pivot_row += 1
        if pivot_row == height:
            break
    return pivot_row


def determinant(rows: list[list[int | Fraction]]) -> Fraction:
    """Exact determinant of a square matrix over Q."""
    matrix = [[Fraction(value) for value in row] for row in rows]
    size = len(matrix)
    assert all(len(row) == size for row in matrix)
    value = Fraction(1)
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if matrix[row][column]),
            None,
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
            value = -value
        scale = matrix[column][column]
        value *= scale
        for row in range(column + 1, size):
            if not matrix[row][column]:
                continue
            factor = matrix[row][column] / scale
            matrix[row] = [
                left - factor * right
                for left, right in zip(matrix[row], matrix[column], strict=True)
            ]
    return value


def load_thm741():
    spec = importlib.util.spec_from_file_location("thm741_heawood_dependency", THM741_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def gl32_maps() -> tuple[tuple[int, ...], ...]:
    maps = []
    for x, y, z in product(POINTS, repeat=3):
        span = {x, y, z, x ^ y, x ^ z, y ^ z, x ^ y ^ z}
        if len(span) != 7:
            continue
        image = [0]
        for point in POINTS:
            image.append(
                (x if point & 1 else 0)
                ^ (y if point & 2 else 0)
                ^ (z if point & 4 else 0)
            )
        maps.append(tuple(image))
    return tuple(maps)


def flag(edge: tuple[int, int]) -> tuple[int, tuple[int, int, int]]:
    left, right = edge
    point = left ^ right
    line = tuple(sorted((left, right, point)))
    assert point in line and line in LINES
    return point, line


def hexagon_current(triple: tuple[int, int, int]) -> list[int]:
    """Alternating current around the Heawood hexagon from a Fano basis."""
    a, b, c = triple
    assert a ^ b ^ c != 0
    point_line_steps = (
        (a, a, b),
        (b, a, b),
        (b, b, c),
        (c, b, c),
        (c, c, a),
        (a, c, a),
    )
    current = [0] * len(EDGES)
    for step, (point, left, right) in enumerate(point_line_steps):
        line = {left, right, left ^ right}
        edge = frozenset(line - {point})
        current[EDGE_INDEX[edge]] += 1 if step % 2 == 0 else -1
    assert sum(value != 0 for value in current) == 6
    return current


def dot(left: list[int], right: list[int | Fraction]) -> Fraction:
    return sum(
        (Fraction(a) * Fraction(b) for a, b in zip(left, right, strict=True)),
        Fraction(0),
    )


def normalize_direction(triple: tuple[int, int, int]) -> tuple[int, int]:
    left = (triple[0] - triple[2]) % P
    right = (triple[1] - triple[2]) % P
    assert left or right
    inverse = pow(left if left else right, -1, P)
    return left * inverse % P, right * inverse % P


def reflect_direction(direction: tuple[int, int]) -> tuple[int, int]:
    return normalize_direction((direction[1], direction[0], 0))


def radon_rows(
    directions: tuple[tuple[int, int], ...],
    domain: tuple[tuple[int, int], ...],
) -> list[list[int]]:
    return [
        [int((left * a + right * b) % P == level) for a, b in domain]
        for left, right in directions
        for level in range(P)
    ]


def chi7(value: int) -> int:
    value %= P
    if value == 0:
        return 0
    return 1 if value in (1, 2, 4) else -1


def format_fraction_tuple(values: tuple[Fraction | int, ...]) -> str:
    return "(" + ",".join(str(value) for value in values) + ")"


def run() -> str:
    # K7 edges are Fano flags, i.e. the 21 edges of the Heawood graph.
    flags = tuple(flag(edge) for edge in EDGES)
    assert len(set(flags)) == 21
    point_rows = [[int(vertex in edge) for edge in EDGES] for vertex in POINTS]
    fano_rows = [
        [int(set(edge).issubset(set(line))) for edge in EDGES]
        for line in LINES
    ]
    flag_point_rows = [
        [int(flag(edge)[0] == vertex) for edge in EDGES]
        for vertex in POINTS
    ]
    for vertex_index, vertex in enumerate(POINTS):
        reconstructed = [
            sum(
                fano_rows[line_index][edge_index]
                for line_index, line in enumerate(LINES)
                if vertex in line
            )
            - point_rows[vertex_index][edge_index]
            for edge_index in range(len(EDGES))
        ]
        assert reconstructed == flag_point_rows[vertex_index]
    assert (
        rational_rank(point_rows),
        rational_rank(fano_rows),
        rational_rank(point_rows + fano_rows),
        rational_rank(flag_point_rows + fano_rows),
    ) == (7, 7, 13, 13)

    maps = gl32_maps()
    assert len(maps) == 168
    line_sets = {frozenset(line) for line in LINES}
    assert all(
        frozenset(point_map[point] for point in line) in line_sets
        for point_map in maps
        for line in LINES
    )

    # Connectedness and absence of four-cycles certify Heawood girth six.
    graph: dict[tuple[str, object], set[tuple[str, object]]] = {}
    for point, line in flags:
        point_node = ("p", point)
        line_node = ("l", line)
        graph.setdefault(point_node, set()).add(line_node)
        graph.setdefault(line_node, set()).add(point_node)
    seen = {next(iter(graph))}
    todo = list(seen)
    while todo:
        node = todo.pop()
        for neighbor in graph[node]:
            if neighbor not in seen:
                seen.add(neighbor)
                todo.append(neighbor)
    assert len(seen) == 14
    assert all(
        len(set(left).intersection(right)) == 1
        for left, right in combinations(LINES, 2)
    )

    noncollinear = tuple(
        triple for triple in combinations(POINTS, 3) if triple[0] ^ triple[1] ^ triple[2]
    )
    assert len(noncollinear) == 28
    all_hexagons = [hexagon_current(triple) for triple in noncollinear]
    basis_currents = [hexagon_current(triple) for triple in HEXAGON_BASIS]
    assert rational_rank(all_hexagons) == rational_rank(basis_currents) == 8
    assert all(
        dot(current, marginal) == 0
        for current in all_hexagons
        for marginal in point_rows + fano_rows
    )

    pivot_edges = ((1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (2, 4), (2, 5), (2, 6))
    basis_minor = [
        [current[EDGE_INDEX[frozenset(edge)]] for edge in pivot_edges]
        for current in basis_currents
    ]
    basis_det = determinant(basis_minor)
    assert basis_det == 1
    recovery_rows = point_rows + fano_rows[:6] + basis_currents
    recovery_det = determinant(recovery_rows)
    assert recovery_det == -100842  # -2*3*7^5

    # Exact THM-741 observables on the 21 flood bodies.
    thm741 = load_thm741()
    observables: dict[str, list[int | Fraction]] = {"r": [], "m": [], "V1": []}
    for edge in EDGES:
        body = tuple(sorted((*range(8, 15), *edge)))
        _, r_edge, m_edge = thm741.good_norm(body)
        threshold = 3 * m_edge / (thm741.S2 * r_edge)
        v1_edge = thm741.minV(4, threshold.numerator, threshold.denominator)
        observables["r"].append(r_edge)
        observables["m"].append(m_edge)
        observables["V1"].append(v1_edge)

    basis_charges = {
        name: tuple(dot(current, values) for current in basis_currents)
        for name, values in observables.items()
    }
    assert basis_charges["r"] == (Fraction(0),) * 8
    assert basis_charges["m"] == tuple(
        map(
            Fraction,
            ("1/140", "-1/42", "-1/210", "-1/42", "-41/1260", "-1/35", "-1/28", "-1/42"),
        )
    )
    assert basis_charges["V1"] == tuple(map(Fraction, (44, -20, 15, 15, 65, -6, 25, 1)))
    all_charges = {
        name: tuple(dot(current, values) for current in all_hexagons)
        for name, values in observables.items()
    }
    assert all(value == 0 for value in all_charges["r"])
    m_zero_triples = tuple(
        triple for triple, value in zip(noncollinear, all_charges["m"], strict=True) if value == 0
    )
    assert m_zero_triples == ((1, 2, 4), (4, 5, 7))
    assert all(value != 0 for value in all_charges["V1"])

    # Depth-11 address/carry channels.
    addresses = tuple(
        (alpha, beta, R - alpha - beta)
        for alpha in range(R + 1)
        for beta in range(R - alpha + 1)
    )
    carry_roles = []
    for address in addresses:
        residue = tuple(value % P for value in address)
        quotient = tuple(
            (value - base) // P for value, base in zip(address, residue, strict=True)
        )
        carry_roles.append(-1 if sum(quotient) == 0 else quotient.index(1))
    carry_domains = {
        role: tuple(
            address[:2]
            for address, address_role in zip(addresses, carry_roles, strict=True)
            if address_role == role
        )
        for role in (-1, 0, 1, 2)
    }
    no_carry = tuple(
        (a, b) for a in range(P) for b in range(P) if 5 <= a + b <= 11
    )
    assert carry_domains[-1] == no_carry
    assert tuple(len(carry_domains[role]) for role in (-1, 0, 1, 2)) == (33, 15, 15, 15)

    def channel_ranks(subset: tuple[tuple[int, int], ...]) -> tuple[int, ...]:
        return tuple(
            rational_rank(radon_rows(subset, carry_domains[role]))
            for role in (-1, 0, 1, 2)
        )

    reflection_orbits = (
        frozenset(((1, 0), (0, 1))),
        frozenset(((1, 1),)),
        frozenset(((1, 2), (1, 4))),
        frozenset(((1, 6),)),
        frozenset(((1, 3), (1, 5))),
    )
    assert all(
        frozenset(map(reflect_direction, orbit)) == orbit for orbit in reflection_orbits
    )
    invariant_five = tuple(
        subset
        for subset in combinations(DIRECTIONS, 5)
        if frozenset(map(reflect_direction, subset)) == frozenset(subset)
    )
    invariant_six = tuple(
        subset
        for subset in combinations(DIRECTIONS, 6)
        if frozenset(map(reflect_direction, subset)) == frozenset(subset)
    )
    assert len(invariant_five) == 6 and len(invariant_six) == 4
    five_rank_rows = tuple(channel_ranks(subset) for subset in invariant_five)
    six_rank_rows = tuple(channel_ranks(subset) for subset in invariant_six)
    assert Counter(row[0] for row in five_rank_rows) == {28: 1, 29: 2, 30: 3}
    assert all(row[1:] == (15, 15, 15) for row in five_rank_rows)
    assert Counter(row[0] for row in six_rank_rows) == {31: 1, 32: 2, 33: 1}
    assert all(row[1:] == (15, 15, 15) for row in six_rank_rows)
    assert channel_ranks(TRIALITY_SIX) == (32, 15, 15, 15)
    assert channel_ranks(REFLECTION_SIX) == (33, 15, 15, 15)
    assert [
        subset for subset in invariant_six if channel_ranks(subset) == (33, 15, 15, 15)
    ] == [REFLECTION_SIX]

    # Five arbitrary pencils have universal rank at most 1+6*5=31; this one attains it.
    rank_31_five = ((1, 0), (0, 1), (1, 2), (1, 4), (1, 3))
    assert rational_rank(radon_rows(rank_31_five, no_carry)) == 31

    # The unique triality-six alias is an alternating chi7 Vandermonde character.
    vandermonde = []
    for a, b in no_carry:
        c = R - a - b
        assert 0 <= c < P
        vandermonde.append(-chi7((a - b) * (b - c) * (c - a)))
    assert Counter(vandermonde) == {-1: 12, 0: 9, 1: 12}
    triality_profiles = tuple(
        tuple(
            sum(
                value
                for (a, b), value in zip(no_carry, vandermonde, strict=True)
                if (left * a + right * b) % P == level
            )
            for level in range(P)
        )
        for left, right in TRIALITY_SIX
    )
    assert triality_profiles == ((0,) * P,) * len(TRIALITY_SIX)
    chi_profiles = tuple(
        tuple(
            sum(
                value
                for (a, b), value in zip(no_carry, vandermonde, strict=True)
                if (left * a + right * b) % P == level
            )
            for level in range(P)
        )
        for left, right in ((1, 3), (1, 5))
    )
    assert chi_profiles == (
        (4, -3, -3, -3, 4, 4, -3),
        (-4, 3, 3, 3, -4, 3, -4),
    )

    full_rows = radon_rows(REFLECTION_SIX, no_carry)
    full_minor_indices = (
        *range(0, 13),
        *range(14, 20),
        *range(21, 27),
        *range(28, 34),
        35,
        36,
    )
    assert len(full_minor_indices) == 33
    full_det = determinant([full_rows[index] for index in full_minor_indices])
    assert full_det == -352947  # -3*7^6

    # Tournament fingerprint on the eight minimal cycle probes.
    magnitudes = [abs(value) for value in basis_charges["V1"]]
    scores = [0] * len(HEXAGON_BASIS)
    edge_flips_vs_basis = 0
    for left, right in combinations(range(len(HEXAGON_BASIS)), 2):
        winner = left if magnitudes[left] >= magnitudes[right] else right
        scores[winner] += 1
        edge_flips_vs_basis += int(winner != left)
    hamiltonian_path = tuple(
        HEXAGON_BASIS[index]
        for index in sorted(
            range(len(HEXAGON_BASIS)),
            key=lambda index: (-magnitudes[index], index),
        )
    )
    assert sorted(scores) == list(range(8))
    assert edge_flips_vs_basis == 8
    assert hamiltonian_path == (
        (1, 3, 5),
        (1, 2, 5),
        (1, 3, 7),
        (1, 2, 6),
        (1, 2, 7),
        (1, 3, 4),
        (1, 3, 6),
        (1, 4, 6),
    )

    lines = [
        "THM-850 ADDENDUM: FANO-FLAG/HEAWOOD AND CHI7 RADON REPAIR",
        "=" * 76,
        "method: exact Fraction row reduction, determinants, and finite incidence",
        "",
        "FANO FLAGS AND THE INVISIBLE EDGE SECTOR",
        f"K7 edges / Fano flags / Heawood edges={len(EDGES)}; vertices=14; connected=yes; girth=6",
        "point-star/Fano/joint marginal ranks=7/7/13",
        "flag identity: row_v=sum_(L contains v) column_L-star_v",
        "invisible sector=ker Heawood boundary=H1(Heawood;Q), dimension=21-14+1=8",
        f"GL(3,2) maps audited={len(maps)}; marginal orbit span remains rank 13",
        "",
        "ORIENTED HEXAGON REPAIR",
        f"noncollinear triples / Heawood hexagons={len(noncollinear)}; orbit rank={rational_rank(all_hexagons)}",
        f"basis={HEXAGON_BASIS}",
        f"basis pivot edges={pivot_edges}; minor determinant={basis_det}",
        f"13 marginals + 8 cycles recovery determinant={recovery_det}=-2*3*7^5",
        "minimum invisible support=6 (Heawood girth); minimum scalar repair count=8",
        f"r basis charges={format_fraction_tuple(basis_charges['r'])}",
        f"m basis charges={format_fraction_tuple(basis_charges['m'])}",
        f"V1 basis charges={format_fraction_tuple(basis_charges['V1'])}",
        f"all-hexagon nonzero counts r/m/V1=0/{28-len(m_zero_triples)}/28; m zeros={m_zero_triples}",
        "",
        "DEPTH-11 CARRY RADON REPAIR",
        "channel sizes no-carry/A/B/C=(33,15,15,15)",
        f"reflection direction orbits={reflection_orbits}",
        "invariant five no-carry rank census={28:1,29:2,30:3}; carry ranks always 15,15,15",
        "general five-direction upper/attained rank=31/31<33",
        "invariant six no-carry rank census={31:1,32:2,33:1}; carry ranks always 15,15,15",
        f"triality six ranks={channel_ranks(TRIALITY_SIX)}; nullity=1",
        "triality kernel f(a,b,c)=-chi7((a-b)(b-c)(c-a)); values -/0/+=12/9/12",
        f"chi-pair detection profiles={chi_profiles}",
        f"unique reflection-full six={REFLECTION_SIX}; ranks={channel_ranks(REFLECTION_SIX)}",
        f"full-rank 33x33 minor determinant={full_det}=-3*7^6",
        "",
        "TOURNAMENT ANALYSIS ON EIGHT HEXAGON PROBES (NOT RUNNERS/ARCS)",
        "pair observable=absolute V1 cycle charge; switch=larger wins; tie=basis order",
        f"score_hist={dict(Counter(scores))}; C3=0; SCC=8x1; edge_flips_vs_basis={edge_flips_vs_basis}; HP=1",
        f"Hamiltonian path={hamiltonian_path}",
        "PRESERVES: all point/Fano marginals plus the full eight-dimensional flag curl",
        "DESTROYS: interval locations, lower-tree chronology, owner data, and the LRC metric predicate",
        "CHALLENGED ASSUMPTION: oriented Fano flags/hexagons, not runners or arcs, are the probe vertices",
        "",
        "SCOPE: exact carrier recovery only; no THM-741 sweep or j=4 metric closure",
        f"source_sha256={sha256(Path(__file__).resolve())}",
        f"thm741_dependency_sha256={sha256(THM741_PATH)}",
        "ALL EXACT ASSERTIONS PASSED",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    if not __debug__:
        raise RuntimeError("verification requires assertions; do not use python -O")
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "05-knowledge/results/lrc14_fano_flag_heawood_radon_repair_codex_S15.out"
        ),
    )
    args = parser.parse_args()
    text = run()
    args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
