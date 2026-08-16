#!/usr/bin/env python3
"""Exact role-chart probe from a THM-3479 relation tuple to the 7x13 carrier.

This is a finite-exact structural component of promoted THM-3479.  It maps
relation-coordinate *potentials* to exact graph coboundaries.  It does not map
THM-2334 endpoint currents, phases, or physical LRC words to carrier edges.
"""

from __future__ import annotations

import ast
from fractions import Fraction
from functools import reduce
from hashlib import sha256
from itertools import permutations
from json import dumps
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

SOURCE_PINS = (
    (
        "THM-3473",
        ROOT / "01-canon/theorems/THM-3473-three-times-p-eight-owner-private-sheet-partition-and-irredundancy.md",
        "05236a7338a4c92b443b2798508d407e2e4f82d942111d022064f6ec9fe86ca2",
        "PROVED",
    ),
    (
        "THM-3479",
        ROOT / "01-canon/theorems/THM-3479-literal-half-twist-relation-current-two-transplant-certificate.md",
        "025998551e3cdf3c6e4db5c0a4f208dd32f6845970fd4729d4a276035e0fdfeb",
        "PROVED-STRUCTURAL-AUDITED",
    ),
)

EXPECTED_SEMANTIC_SHA256 = (
    "1992d46a4df3a3c862e9c36ecfc4f992eea654d9373b852836c9652962ad49e7"
)

LABELS = ("c1", "c2", "c3", "H", "q1", "q2", "q3", "q4", "q5")
RELATION = (-27, -27, -27, 20110798, -41, -27, -27, -27, 38)
REL = dict(zip(LABELS, RELATION))

U_FULL = {
    "c1": 13,
    "c2": 2197,
    "c3": 742586,
    "H": 1,
    "q1": 183,
    "q2": 27,
    "q3": 131,
    "q4": 53,
    "q5": 313,
}
U_CLOCK = {
    "c1": 65,
    "c2": 2197,
    "c3": 742586,
    "H": 5,
    "q1": 661549,
    "q2": 655231,
    "q3": 658533,
    "q4": 661445,
    "q5": 291,
}

# Vertices are u1,...,u8, encoded as 0,...,7.  Edges inherit the owner-order
# orientation from the smaller endpoint to the larger endpoint.
EDGES = (
    (0, 3), (0, 4), (0, 5),
    (1, 2), (1, 4), (1, 7),
    (2, 4), (2, 7),
    (3, 4), (3, 5),
    (4, 5), (4, 6), (4, 7),
)
EDGE_INDEX = {edge: index for index, edge in enumerate(EDGES)}
HUB = 4
LEAF = 6
WINGS = ((0, 3, 5), (1, 2, 7))
BLOCKERS = ("c1", "c2", "c3")
MIDDLE_UNITS = ("q2", "q3", "q4")

# Three fundamental triangles in each K4.
CYCLES = (
    (0, 3, 4), (0, 3, 5), (0, 4, 5),
    (1, 2, 4), (1, 2, 7), (1, 4, 7),
)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def source_hashes() -> tuple[tuple[str, str, str], ...]:
    rows = tuple(
        (label, status, lf_sha256(path))
        for label, path, _expected, status in SOURCE_PINS
    )
    for row, pin in zip(rows, SOURCE_PINS):
        expected = pin[2]
        if expected != "UNFROZEN":
            require(row[2] == expected, ("source hash drift", row, expected))
    return rows


def security_certificate() -> tuple[object, ...]:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
            "assert statement forbidden")
    forbidden = {"eval", "exec", "compile", "__import__"}
    calls = {
        node.func.id
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }
    require(not (calls & forbidden), ("dynamic execution", calls & forbidden))
    imports = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            imports.add(node.module or "")
    allowed = {
        "__future__", "ast", "fractions", "functools", "hashlib",
        "itertools", "json", "math", "pathlib",
    }
    require(imports <= allowed, ("unexpected imports", imports - allowed))
    return ("NO_ASSERT_AST", tuple(sorted(imports)), tuple(sorted(forbidden)))


def dot_relation(values: dict[str, int | Fraction]) -> Fraction:
    return sum((Fraction(REL[label]) * values[label] for label in LABELS),
               Fraction(0))


def recover_q1(other_values: dict[str, int | Fraction]) -> Fraction:
    numerator = sum(
        (Fraction(REL[label]) * other_values[label]
         for label in LABELS if label != "q1"),
        Fraction(0),
    )
    return -numerator / REL["q1"]


def matrix_rank(matrix: list[list[int | Fraction]]) -> int:
    rows = [[Fraction(value) for value in row] for row in matrix]
    if not rows:
        return 0
    height = len(rows)
    width = len(rows[0])
    pivot_row = 0
    for column in range(width):
        pivot = next(
            (row for row in range(pivot_row, height) if rows[row][column]),
            None,
        )
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        scale = rows[pivot_row][column]
        rows[pivot_row] = [value / scale for value in rows[pivot_row]]
        for row in range(height):
            if row == pivot_row or rows[row][column] == 0:
                continue
            factor = rows[row][column]
            rows[row] = [
                left - factor * right
                for left, right in zip(rows[row], rows[pivot_row])
            ]
        pivot_row += 1
        if pivot_row == height:
            break
    return pivot_row


def determinant(matrix: list[list[int]]) -> int:
    rows = [row[:] for row in matrix]
    size = len(rows)
    require(all(len(row) == size for row in rows), "nonsquare determinant")
    if size == 0:
        return 1
    sign = 1
    previous = 1
    for pivot_index in range(size - 1):
        if rows[pivot_index][pivot_index] == 0:
            swap = next(
                (row for row in range(pivot_index + 1, size)
                 if rows[row][pivot_index]),
                None,
            )
            if swap is None:
                return 0
            rows[pivot_index], rows[swap] = rows[swap], rows[pivot_index]
            sign = -sign
        pivot = rows[pivot_index][pivot_index]
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                numerator = (
                    rows[row][column] * pivot
                    - rows[row][pivot_index] * rows[pivot_index][column]
                )
                require(numerator % previous == 0,
                        ("Bareiss exact division", numerator, previous))
                rows[row][column] = numerator // previous
        previous = pivot
        for row in range(pivot_index + 1, size):
            rows[row][pivot_index] = 0
    return sign * rows[-1][-1]


def permutation_order(permutation: tuple[int, ...]) -> int:
    seen = set()
    answer = 1
    for start in range(len(permutation)):
        if start in seen:
            continue
        cursor = start
        length = 0
        while cursor not in seen:
            seen.add(cursor)
            cursor = permutation[cursor]
            length += 1
        answer = answer * length // gcd(answer, length)
    return answer


def graph_automorphisms() -> tuple[tuple[int, ...], ...]:
    edge_set = frozenset(frozenset(edge) for edge in EDGES)
    autos = []
    for permutation in permutations(range(8)):
        image = frozenset(
            frozenset((permutation[left], permutation[right]))
            for left, right in EDGES
        )
        if image == edge_set:
            autos.append(tuple(permutation))
    return tuple(autos)


def role_charts() -> tuple[tuple[str, ...], ...]:
    charts = []
    for swap in (0, 1):
        blocker_wing = WINGS[swap]
        unit_wing = WINGS[1 - swap]
        for blocker_order in permutations(BLOCKERS):
            for unit_order in permutations(MIDDLE_UNITS):
                chart = {HUB: "H", LEAF: "q5"}
                chart.update(zip(blocker_wing, blocker_order))
                chart.update(zip(unit_wing, unit_order))
                charts.append(tuple(chart[vertex] for vertex in range(8)))
    return tuple(sorted(set(charts)))


def automorphism_orbit(
    chart: tuple[str, ...], autos: tuple[tuple[int, ...], ...]
) -> tuple[tuple[str, ...], ...]:
    orbit = set()
    for permutation in autos:
        image = [""] * 8
        for vertex, label in enumerate(chart):
            image[permutation[vertex]] = label
        orbit.add(tuple(image))
    return tuple(sorted(orbit))


def edge_weights(
    values: dict[str, int], chart: tuple[str, ...]
) -> tuple[int, ...]:
    potentials = tuple(values[label] for label in chart)
    return tuple(potentials[left] - potentials[right] for left, right in EDGES)


def cycle_vector(cycle: tuple[int, ...]) -> tuple[int, ...]:
    vector = [0] * len(EDGES)
    closed = cycle + (cycle[0],)
    for left, right in zip(closed, closed[1:]):
        edge = tuple(sorted((left, right)))
        require(edge in EDGE_INDEX, ("nonedge in cycle", cycle, edge))
        vector[EDGE_INDEX[edge]] += 1 if left < right else -1
    return tuple(vector)


def weighted_tree_determinant(weights: tuple[int, ...]) -> int:
    laplacian = [[0] * 8 for _ in range(8)]
    for (left, right), weight in zip(EDGES, weights):
        laplacian[left][left] += weight
        laplacian[right][right] += weight
        laplacian[left][right] -= weight
        laplacian[right][left] -= weight
    retained = tuple(vertex for vertex in range(8) if vertex != HUB)
    minor = [[laplacian[row][column] for column in retained] for row in retained]
    return determinant(minor)


def k4_tree_sum(
    values: dict[str, int], chart: tuple[str, ...], vertices: tuple[int, ...]
) -> int:
    positions = {vertex: index for index, vertex in enumerate(vertices)}
    laplacian = [[0] * 4 for _ in range(4)]
    for left, right in EDGES:
        if left not in positions or right not in positions:
            continue
        weight = values[chart[left]] - values[chart[right]]
        i = positions[left]
        j = positions[right]
        laplacian[i][i] += weight
        laplacian[j][j] += weight
        laplacian[i][j] -= weight
        laplacian[j][i] -= weight
    minor = [row[:-1] for row in laplacian[:-1]]
    return determinant(minor)


def valuation(value: int, prime: int) -> int:
    require(value != 0, "valuation of zero")
    exponent = 0
    value = abs(value)
    while value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


def chart_census(
    name: str,
    values: dict[str, int],
    charts: tuple[tuple[str, ...], ...],
    cycle_basis: tuple[tuple[int, ...], ...],
) -> tuple[object, ...]:
    rows = []
    for chart in charts:
        weights = edge_weights(values, chart)
        pairings = tuple(
            sum(coefficient * weight
                for coefficient, weight in zip(cycle, weights))
            for cycle in cycle_basis
        )
        require(pairings == (0,) * 6, (name, "nonzero H1 pairing", pairings))
        determinant_value = weighted_tree_determinant(weights)
        wing_a = k4_tree_sum(values, chart, WINGS[0] + (HUB,))
        wing_b = k4_tree_sum(values, chart, WINGS[1] + (HUB,))
        bridge = weights[EDGE_INDEX[(HUB, LEAF)]]
        require(determinant_value == wing_a * wing_b * bridge,
                (name, "tree factorization", determinant_value,
                 wing_a, wing_b, bridge))
        rows.append((chart, weights, pairings, wing_a, wing_b, bridge,
                     determinant_value))

    determinants = tuple(row[-1] for row in rows)
    require(all(determinants), (name, "zero role-chart determinant"))
    valuations = tuple(valuation(value, 13) for value in determinants)
    require(set(valuations) == {4}, (name, "13-adic valuation", valuations))
    determinant_gcd = reduce(gcd, (abs(value) for value in determinants))
    signs = (
        sum(value > 0 for value in determinants),
        sum(value < 0 for value in determinants),
    )
    canonical = (
        "c1", "q2", "q3", "c2", "H", "c3", "q5", "q4",
    )
    canonical_row = next(row for row in rows if row[0] == canonical)
    return (
        name,
        len(rows),
        sum(value != 0 for value in determinants),
        len(set(determinants)),
        signs,
        min(abs(value) for value in determinants),
        max(abs(value) for value in determinants),
        determinant_gcd,
        valuations[0],
        canonical_row[3],
        canonical_row[4],
        canonical_row[5],
        canonical_row[6],
        tuple(values[label] for label in canonical),
    )


def main() -> None:
    hashes = source_hashes()
    security = security_certificate()

    require(tuple(label for label, coefficient in zip(LABELS, RELATION)
                  if coefficient == -27)
            == ("c1", "c2", "c3", "q2", "q3", "q4"),
            "six-slot coefficient packet drift")
    require(REL["q1"] == -41 and REL["H"] == 20110798
            and REL["q5"] == 38, "distinguished coefficient drift")

    for name, values in (("U_full", U_FULL), ("U_clock", U_CLOCK)):
        require(dot_relation(values) == 0, (name, "relation failure"))
        require(recover_q1(values) == values["q1"],
                (name, "q1 chart recovery"))

    generic_one = {
        label: Fraction(1) for label in LABELS if label != "q1"
    }
    gauge_q1 = recover_q1(generic_one)
    require(gauge_q1 == Fraction(20110674, 41), gauge_q1)
    require(gcd(abs(REL["q1"]), abs(REL["q5"])) == 1,
            "integral chart index witness lost")

    charts = role_charts()
    require(len(charts) == 72, len(charts))
    autos = graph_automorphisms()
    require(len(autos) == 72, len(autos))
    require(tuple(sorted({permutation_order(auto) for auto in autos}))
            == (1, 2, 3, 4, 6), "automorphism order drift")
    canonical = (
        "c1", "q2", "q3", "c2", "H", "c3", "q5", "q4",
    )
    require(automorphism_orbit(canonical, autos) == charts,
            "role charts are not one automorphism orbit")

    # The potential-to-edge matrix has rank seven and kernel the diagonal.
    columns = []
    for source_vertex in range(8):
        potential = [0] * 8
        potential[source_vertex] = 1
        columns.append(tuple(
            potential[left] - potential[right] for left, right in EDGES
        ))
    coboundary_matrix = [
        [columns[column][row] for column in range(8)]
        for row in range(13)
    ]
    require(matrix_rank(coboundary_matrix) == 7,
            "coboundary rank drift")
    require(all(sum(row) == 0 for row in coboundary_matrix),
            "diagonal is not the kernel")

    cycle_basis = tuple(cycle_vector(cycle) for cycle in CYCLES)
    require(matrix_rank([list(cycle) for cycle in cycle_basis]) == 6,
            "cycle basis rank drift")
    require(13 - 8 + 1 == 6, "H1 dimension drift")

    full_census = chart_census(
        "U_full", U_FULL, charts, cycle_basis
    )
    clock_census = chart_census(
        "U_clock", U_CLOCK, charts, cycle_basis
    )
    require(full_census[7] == 10967424, full_census[7])
    require(clock_census[7] == 5026736, clock_census[7])
    require(full_census[9:13] == (
        -408230882586462649,
        140608,
        -312,
        17908964716879810126984704,
    ), full_census[9:13])
    require(clock_census[9:13] == (
        -408094818615760385,
        -282469436237523288,
        -286,
        -32968453616912573701504956325888921680,
    ), clock_census[9:13])

    flat_weights = edge_weights({label: 1 for label in LABELS}, canonical)
    require(flat_weights == (0,) * 13, flat_weights)
    require(weighted_tree_determinant(flat_weights) == 0,
            "flat hostile survives")

    semantic_payload = {
        "source_hashes": hashes,
        "security": security,
        "labels": LABELS,
        "relation": RELATION,
        "equal_packet": ("c1", "c2", "c3", "q2", "q3", "q4"),
        "gauge": ("q1", -41, str(gauge_q1), 41),
        "graph": (EDGES, HUB, LEAF, WINGS),
        "automorphisms": (len(autos), (1, 2, 3, 4, 6)),
        "role_charts": (len(charts), 1),
        "cohomology": (7, 6, tuple((0,) * 6 for _ in range(2))),
        "full_census": full_census,
        "clock_census": clock_census,
        "flat_hostile": 0,
    }
    semantic_hash = sha256(
        dumps(semantic_payload, sort_keys=True, separators=(",", ":"),
              default=str).encode("ascii")
    ).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "UNFROZEN":
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256,
                ("semantic digest", semantic_hash,
                 EXPECTED_SEMANTIC_SHA256))

    print("LRC RELATION ROLE-CHART WEIGHTED CLOSURE FINITE-EXACT SIDECAR")
    print("STATUS: FINITE-EXACT STRUCTURAL COMPONENT OF PROVED THM-3479; INDEPENDENTLY AUDITED")
    print(f"SOURCE_HASHES: {hashes}")
    print(f"SECURITY: {security}")
    print(f"RELATION_LABELS: {LABELS}")
    print(f"RELATION: {RELATION}")
    print("SIX_EQUAL_MINUS27_SLOTS: (c1,c2,c3,q2,q3,q4)")
    print("ROLE_CONTRACT: delete q1 as rational relation gauge; H->u5 hub; q5->u7 leaf; blocker triple and q2/q3/q4 triple -> the two K4 wings")
    print(f"Q1_GAUGE_LIFT_OF_CONSTANT_POTENTIAL: {gauge_q1}")
    print("INTEGRAL_BOUNDARY: q1 projection has index 41 in Z^8; the chart is an isomorphism only over Q")
    print(f"GRAPH_AUTOMORPHISMS_AND_ORDERS: {(len(autos), (1, 2, 3, 4, 6))}")
    print(f"ROLE_CHARTS_AND_AUT_ORBITS: {(len(charts), 1)}")
    print("COCHAIN_MAP: ker(a:Q^9->Q)/<diagonal gauge> ~= B^1(G;Q), rank 7")
    print("ABSOLUTE_H1: dim H^1(G;Q)=6; all six cycle pairings vanish in every chart")
    print(f"U_FULL_CENSUS: {full_census}")
    print(f"U_CLOCK_CENSUS: {clock_census}")
    print("UNIFORM_13_ADIC_LAW: every one of 72 determinants has valuation_13 exactly 4 for each tuple")
    print("CANONICAL_FACTORIZATION: determinant=(u5-u7 bridge weight)*(first K4 tree sum)*(second K4 tree sum)")
    print("HOSTILE_FLAT_POTENTIAL: determinant=0")
    print(f"SEMANTIC_SHA256: {semantic_hash}")
    print("VERDICT: the role-labelled relation-potential transplant is unique up to graph automorphism and has robust nonzero weighted closure on both tuples, but it is an exact coboundary with zero absolute H1 and is not an endpoint current, phase map, bispectrum, physical row, or LRC(14) proof")


if __name__ == "__main__":
    main()
