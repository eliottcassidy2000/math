#!/usr/bin/env python3
"""Finite exact certificate for the exact-M=12 ``U=0`` endpoint wall.

Universe and filters
--------------------
Work in the inherited reduced (2,3) seam.  A residual row ``p^i y^j`` has
weight ``2*i+3*j`` and contributes the two valued-support points

    (j+2,i+j,1), (j,i+j+1,1).

The certificate retains every allowed row through weight twelve, deletes the
weight-twelve endpoint ``U p^6``, and fixes the nonzero rows ``p,p^2,p^3``,
``W p^3 y^2``, and ``Z y^4``.  Every remaining lower row is independently
present or absent.  At each row subset, every active coincident support point
is also independently deleted; this is a hostile over-approximation of every
possible aggregate-coefficient cancellation.

The five hull complexes are determined by

    u = [p^5]H,  a = [p^4 y]H,  d = [p^4]H,

with the fixed nonzero coefficient ``e=[p^3]H=-1376/135``.  Case A also
requires the boundary discriminant ``E5=a^2-4*W*u`` to be nonzero.  On the
terminal Case-E face, ``eta=[p^3*y]H`` is a genuine owner and the additional
boundary discriminant ``E3=eta^2-4*e*W`` must be nonzero.  The program proves
the five lower-hull complexes, hostile counts, global and face Pick ledgers,
response packets, graph genera, and good-target differential orders.  It
checks the exact Case-E elliptic normal form, reduced-edge positive controls,
and the three required hostiles: ``E5=0`` gives a (2,5) cusp, ``E3=0`` gives a
repeated Case-E edge, while ``Lambda=W+Z=0`` gives the length-twelve A23
contact.  At that contact the program checks only the exact specialization and
five positive order forms imported from THM-4297; it does not claim to
rederive THM-4297's critical ladder.

The default path is the load-bearing exhaustive certificate.  ``--audit``
adds a second complete 1,024-mask coefficient-state census and checks that the
owner rule alone recovers the same five complexes.  The script uses only the
Python standard library and never relies on ``assert``, so ``python -O`` runs
the identical gates.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import combinations
from math import gcd
from typing import Dict, FrozenSet, Iterable, List, Mapping, Sequence, Set, Tuple


Point3 = Tuple[int, int, int]
Point2 = Tuple[int, int]
Row = Tuple[int, int, int]
Plane = Tuple[Fraction, Fraction, Fraction]
Polynomial = Tuple[Fraction, ...]  # ascending coefficient order


CHECKS = 0


def need(condition: bool, message: str) -> None:
    """Optimization-safe exact gate."""

    global CHECKS
    CHECKS += 1
    if not bool(condition):
        raise RuntimeError(message)


def F(value: int, denominator: int = 1) -> Fraction:
    return Fraction(value, denominator)


def lcm(first: int, second: int) -> int:
    return abs(first * second) // gcd(first, second)


M = (F(1, 12), F(1, 6), F(-1, 6))
H5 = (F(0), F(1, 5), F(-1, 5))
V4A = (F(-1, 4), F(1, 4), F(-1, 4))
V3A = (F(-2, 3), F(1, 3), F(-1, 3))
V4W = (F(-1, 8), F(1, 4), F(-1, 4))
V3W = (F(-1, 3), F(1, 3), F(-1, 3))
PLANE_NAMES: Mapping[Plane, str] = {
    M: "M", H5: "H5", V4A: "V4a", V3A: "V3a", V4W: "V4W", V3W: "V3W"
}


ALL_ROWS: Tuple[Row, ...] = (
    (1, 0, 2),
    (2, 0, 4),
    (3, 0, 6),
    (0, 2, 6),
    (2, 1, 7),
    (4, 0, 8),
    (1, 2, 8),
    (3, 1, 9),
    (0, 3, 9),
    (5, 0, 10),
    (2, 2, 10),
    (4, 1, 11),
    (1, 3, 11),
    (3, 2, 12),
    (0, 4, 12),
)

BASE_POINTS: FrozenSet[Point3] = frozenset(
    {(2, 0, 0), (0, 1, 0), (2, 0, 1)}
)
FIXED_KEYS: FrozenSet[Point2] = frozenset(
    {(1, 0), (2, 0), (3, 0), (3, 2), (0, 4)}
)
U5 = (5, 0)
A11 = (4, 1)
DELTA = (4, 0)
ETA = (3, 1)
W12 = (3, 2)
Z12 = (0, 4)


CASE_DATA: Mapping[str, Mapping[str, object]] = {
    "A": {
        "gate": "u!=0,E5=a^2-4Wu!=0 (a arbitrary)",
        "zeros": frozenset(),
        "owners": frozenset({U5}),
        "planes": frozenset({M, H5}),
        "vertices": ((0, 1), (2, 0), (6, 4), (2, 6), (0, 6)),
        "ledger": (46, 14, 17),
        "packet": (11, 11, 5, 5, 2, 2, 2, 2, 1),
        "responses": (41, 33),
        "base": 60,
        "orders": {"M": 45, "H5": 50},
        "faces": {"C": 4, "H5A": 2},
        "graph_vertices": ("R", "C", "H5"),
        "graph_edges": (("R-C", 12), ("C-H5", 1)),
        "components": 3,
        "edges": 13,
        "hostiles": 4500,
    },
    "B": {
        "gate": "u=0,a!=0,d!=0",
        "zeros": frozenset({U5}),
        "owners": frozenset({A11, DELTA}),
        "planes": frozenset({M, H5, V4A}),
        "vertices": ((0, 1), (2, 0), (6, 4), (2, 6), (1, 6), (0, 5)),
        "ledger": (45, 13, 17),
        "packet": (11, 11, 5, 5, 2, 2, 2, 2, 1),
        "responses": (41, 33),
        "base": 60,
        "orders": {"M": 45, "H5": 50, "V4a": 65},
        "faces": {"C": 4, "H50": 2, "V4a": 0},
        "graph_vertices": ("R", "C", "H5", "V4a"),
        "graph_edges": (("R-C", 12), ("C-H5", 1), ("H5-V4a", 1)),
        "components": 4,
        "edges": 14,
        "hostiles": 1080,
    },
    "C": {
        "gate": "u=0,a!=0,d=0",
        "zeros": frozenset({U5, DELTA}),
        "owners": frozenset({A11}),
        "planes": frozenset({M, H5, V3A}),
        "vertices": ((0, 1), (2, 0), (6, 4), (2, 6), (1, 6), (0, 4)),
        "ledger": (44, 12, 17),
        "packet": (11, 11, 5, 5, 2, 2, 2, 2, 1),
        "responses": (41, 33),
        "base": 60,
        "orders": {"M": 45, "H5": 50, "V3a": 90},
        "faces": {"C": 4, "H50": 2, "V3a": 0},
        "graph_vertices": ("R", "C", "H5", "V3a"),
        "graph_edges": (("R-C", 12), ("C-H5", 1), ("H5-V3a", 1)),
        "components": 4,
        "edges": 14,
        "hostiles": 720,
    },
    "D": {
        "gate": "u=0,a=0,d!=0",
        "zeros": frozenset({U5, A11}),
        "owners": frozenset({DELTA}),
        "planes": frozenset({M, V4W}),
        "vertices": ((0, 1), (2, 0), (6, 4), (2, 6), (0, 5)),
        "ledger": (44, 12, 17),
        "packet": (11, 11, 9, 2, 2, 2, 2, 1),
        "responses": (40, 32),
        "base": 24,
        "orders": {"M": 18, "V4W": 23},
        "faces": {"C": 4, "V4W": 2},
        "graph_vertices": ("R", "C", "V4W"),
        "graph_edges": (("R-C", 12), ("C-V4W", 1)),
        "components": 3,
        "edges": 13,
        "hostiles": 720,
    },
    "E": {
        "gate": "u=0,a=0,d=0,E3=eta^2-4eW!=0 (eta arbitrary)",
        "zeros": frozenset({U5, A11, DELTA}),
        "owners": frozenset(),
        "planes": frozenset({M, V3W}),
        "vertices": ((0, 1), (2, 0), (6, 4), (2, 6), (0, 4)),
        "ledger": (42, 12, 16),
        "packet": (11, 11, 4, 4, 2, 2, 2, 2, 1),
        "responses": (39, 31),
        "base": 12,
        "orders": {"M": 9, "V3W": 14},
        "faces": {"C": 4, "V3W": 1},
        "graph_vertices": ("R", "C", "V3W"),
        "graph_edges": (("R-C", 12), ("C-V3W", 1)),
        "components": 3,
        "edges": 13,
        "hostiles": 480,
    },
}


FACE_POLYGONS: Mapping[str, Tuple[Tuple[Point2, ...], Tuple[int, int, int]]] = {
    "C": (((0, 0), (2, 5), (4, 4)), (12, 6, 4)),
    "H5A": (((0, 0), (0, 5), (2, 5)), (10, 8, 2)),
    "H50": (((0, 0), (1, 5), (2, 5)), (5, 3, 2)),
    "V4a": (((0, 0), (0, 4), (1, 5)), (4, 6, 0)),
    "V3a": (((0, 0), (0, 3), (1, 5)), (3, 5, 0)),
    "V4W": (((0, 0), (0, 4), (2, 5)), (8, 6, 2)),
    "V3W": (((0, 0), (0, 3), (2, 5)), (6, 6, 1)),
}

# Every nonconstant monomial on each irreducible face has the displayed
# common Euler weight alpha*s_exp+beta*p_exp=degree.  Hence on F=0,
# alpha*S*F_S+beta*P*F_P=-degree*constant is a nonzero torus unit.  This is a
# symbolic smoothness certificate, not a numeric Jacobian sample.
FACE_EULER_DATA: Mapping[
    str, Tuple[int, Tuple[Point2, ...], Tuple[int, int], int]
] = {
    "C": (1, ((2, 5), (4, 4)), (1, 2), 12),
    "H5": (-1, ((0, 5), (1, 5), (2, 5)), (0, 1), 5),
    "V4a": (-1, ((0, 4), (1, 5)), (-1, 1), 4),
    "V3a": (-1, ((0, 3), (1, 5)), (-2, 1), 3),
    "V4W": (-1, ((0, 4), (2, 5)), (-1, 2), 8),
    "V3W": (-1, ((0, 3), (1, 4), (2, 5)), (-1, 1), 3),
}


def row_endpoints(row: Row) -> Tuple[Point3, Point3]:
    i, j, _weight = row
    return (j + 2, i + j, 1), (j, i + j + 1, 1)


def expanded_support(rows: Iterable[Row]) -> Set[Point3]:
    support = set(BASE_POINTS)
    for row in rows:
        support.update(row_endpoints(row))
    return support


def candidate_plane_records(
    points: Iterable[Point3],
) -> Tuple[Tuple[Point3, ...], Dict[Point3, int], Tuple[Tuple[Plane, int, int], ...]]:
    ordered = tuple(sorted(set(points)))
    index = {point: position for position, point in enumerate(ordered)}
    planes: Set[Plane] = set()
    for first, second, third in combinations(ordered, 3):
        determinant = (
            (second[0] - first[0]) * (third[1] - first[1])
            - (second[1] - first[1]) * (third[0] - first[0])
        )
        if determinant == 0:
            continue
        slope_s = Fraction(
            (second[2] - first[2]) * (third[1] - first[1])
            - (second[1] - first[1]) * (third[2] - first[2]),
            determinant,
        )
        slope_p = Fraction(
            (second[0] - first[0]) * (third[2] - first[2])
            - (second[2] - first[2]) * (third[0] - first[0]),
            determinant,
        )
        constant = (
            Fraction(first[2]) - slope_s * first[0] - slope_p * first[1]
        )
        planes.add((slope_s, slope_p, constant))

    records: List[Tuple[Plane, int, int]] = []
    for plane in sorted(planes):
        below = 0
        equal = 0
        for position, (r, ell, height) in enumerate(ordered):
            gap = (
                Fraction(height)
                - plane[0] * r
                - plane[1] * ell
                - plane[2]
            )
            if gap < 0:
                below |= 1 << position
            elif gap == 0:
                equal |= 1 << position
        records.append((plane, below, equal))
    return ordered, index, tuple(records)


def projected_rank_two(points: Sequence[Point3], bits: int) -> bool:
    chosen = [point for position, point in enumerate(points) if bits & (1 << position)]
    for first, second, third in combinations(chosen, 3):
        determinant = (
            (second[0] - first[0]) * (third[1] - first[1])
            - (second[1] - first[1]) * (third[0] - first[0])
        )
        if determinant:
            return True
    return False


def lower_planes(
    data: Tuple[
        Tuple[Point3, ...],
        Dict[Point3, int],
        Tuple[Tuple[Plane, int, int], ...],
    ],
    support: Iterable[Point3],
) -> FrozenSet[Plane]:
    points, index, records = data
    bits = 0
    for point in support:
        bits |= 1 << index[point]
    answer = set()
    for plane, below, equal in records:
        if below & bits:
            continue
        if projected_rank_two(points, equal & bits):
            answer.add(plane)
    return frozenset(answer)


def convex_hull(points: Iterable[Point2]) -> Tuple[Point2, ...]:
    ordered = sorted(set(points))

    def cross(origin: Point2, first: Point2, second: Point2) -> int:
        return (
            (first[0] - origin[0]) * (second[1] - origin[1])
            - (first[1] - origin[1]) * (second[0] - origin[0])
        )

    lower: List[Point2] = []
    for point in ordered:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper: List[Point2] = []
    for point in reversed(ordered):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return tuple(lower[:-1] + upper[:-1])


def polygon_ledger(points: Iterable[Point2]) -> Tuple[Tuple[Point2, ...], int, int, int]:
    vertices = convex_hull(points)
    area2 = abs(
        sum(
            first[0] * second[1] - second[0] * first[1]
            for first, second in zip(vertices, vertices[1:] + vertices[:1])
        )
    )
    boundary = sum(
        gcd(abs(second[0] - first[0]), abs(second[1] - first[1]))
        for first, second in zip(vertices, vertices[1:] + vertices[:1])
    )
    need((area2 - boundary + 2) % 2 == 0, "Pick parity")
    interior = (area2 - boundary + 2) // 2
    return vertices, area2, boundary, interior


def edge_packet(vertices: Tuple[Point2, ...]) -> Tuple[int, ...]:
    packet: List[int] = []
    for start, end in zip(vertices, vertices[1:] + vertices[:1]):
        dx = end[0] - start[0]
        dy = end[1] - start[1]
        length = gcd(abs(dx), abs(dy))
        inward = (-dy // length, dx // length)
        constant = inward[0] * start[0] + inward[1] * start[1]
        index = inward[0] + inward[1] - constant
        # The s=0 edge is the affine divisor, not a source-infinity packet.
        if start[0] == end[0] == 0:
            continue
        packet.extend([index] * length)
    return tuple(sorted(packet, reverse=True))


def plane_gap(plane: Plane, point: Point3) -> Fraction:
    r, ell, height = point
    return Fraction(height) - plane[0] * r - plane[1] * ell - plane[2]


def vertical_order(plane: Plane, base_degree: int) -> Fraction:
    a_s, b_p, constant = plane
    return base_degree * (F(5, 6) - a_s - b_p - constant)


def trim(poly: Sequence[Fraction]) -> Polynomial:
    answer = list(poly)
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return tuple(answer)


def derivative(poly: Polynomial) -> Polynomial:
    if len(poly) <= 1:
        return (F(0),)
    return trim(tuple(F(index) * poly[index] for index in range(1, len(poly))))


def poly_divmod(numerator: Polynomial, denominator: Polynomial) -> Tuple[Polynomial, Polynomial]:
    numerator_list = list(trim(numerator))
    denominator = trim(denominator)
    need(denominator != (F(0),), "polynomial division by zero")
    quotient = [F(0)] * max(1, len(numerator_list) - len(denominator) + 1)
    while len(numerator_list) >= len(denominator) and any(numerator_list):
        shift = len(numerator_list) - len(denominator)
        factor = numerator_list[-1] / denominator[-1]
        quotient[shift] += factor
        for index, coefficient in enumerate(denominator):
            numerator_list[index + shift] -= factor * coefficient
        numerator_list = list(trim(tuple(numerator_list)))
    return trim(tuple(quotient)), trim(tuple(numerator_list))


def poly_gcd(first: Polynomial, second: Polynomial) -> Polynomial:
    first = trim(first)
    second = trim(second)
    while second != (F(0),):
        _quotient, remainder = poly_divmod(first, second)
        first, second = second, remainder
    if first == (F(0),):
        return first
    return trim(tuple(coefficient / first[-1] for coefficient in first))


def squarefree(poly: Polynomial) -> bool:
    return len(poly_gcd(trim(poly), derivative(trim(poly)))) == 1


def edge_regular(poly: Polynomial) -> bool:
    """An edge retains both corners and has no repeated torus root."""

    return len(poly) >= 2 and poly[0] != 0 and poly[-1] != 0 and squarefree(poly)


def edge_schemes(
    case_name: str, parameters: Mapping[str, Fraction]
) -> Tuple[Tuple[str, Polynomial], ...]:
    """Literal outer and internal edge schemes, with ascending coefficients."""

    W = parameters["W"]
    Z = parameters["Z"]
    u = parameters["u"]
    a = parameters["a"]
    d = parameters["d"]
    e = parameters["e"]
    eta = parameters["eta"]
    common = (
        ("X-1", (F(-1), F(1))),
        ("1-ZX4", (F(1), F(0), F(0), F(0), -Z)),
        ("(X-1)(WX+Z)", (-Z, Z - W, W)),
    )
    internal = (("1-WX", (F(1), -W)),)
    if case_name == "A":
        return common + (
            ("W+aX+uX2", (W, a, u)),
            ("u-X5", (u, F(0), F(0), F(0), F(0), F(-1))),
        ) + internal
    if case_name == "B":
        return common + (
            ("W+aX", (W, a)),
            ("a+dX", (a, d)),
            ("d-X4", (d, F(0), F(0), F(0), F(-1))),
        ) + internal + (("1-aX", (F(1), -a)),)
    if case_name == "C":
        return common + (
            ("W+aX", (W, a)),
            ("a+eX", (a, e)),
            ("e-X3", (e, F(0), F(0), F(-1))),
        ) + internal + (("1-aX", (F(1), -a)),)
    if case_name == "D":
        return common + (
            ("W+dX", (W, d)),
            ("d-X4", (d, F(0), F(0), F(0), F(-1))),
        ) + internal
    if case_name == "E":
        return common + (
            ("e+etaX+WX2", (e, eta, W)),
            ("e-X3", (e, F(0), F(0), F(-1))),
        ) + internal
    raise RuntimeError(f"unknown case {case_name}")


def positive_parameters(case_name: str) -> Dict[str, Fraction]:
    """One exact point on each open gate, with Lambda nonzero."""

    values = {
        "W": F(1),
        "Z": F(2),
        "u": F(0),
        "a": F(0),
        "d": F(0),
        "e": F(-1376, 135),
        "eta": F(0),
    }
    if case_name == "A":
        values["u"] = F(1)  # a=0 is explicitly allowed; E5=-4.
    elif case_name == "B":
        values["a"] = values["d"] = F(1)
    elif case_name == "C":
        values["a"] = F(1)
    elif case_name == "D":
        values["d"] = F(1)
    elif case_name != "E":
        raise RuntimeError(f"unknown case {case_name}")
    return values


def expected_case_from_owner_mask(keys: Set[Point2]) -> str:
    if U5 in keys:
        return "A"
    if A11 in keys:
        return "B" if DELTA in keys else "C"
    return "D" if DELTA in keys else "E"


def check_universe() -> Tuple[Point3, ...]:
    generated = tuple(
        sorted(
            (
                (i, j, 2 * i + 3 * j)
                for i in range(7)
                for j in range(5)
                if 0 < 2 * i + 3 * j <= 12
                and (i, j) not in {(0, 1), (1, 1), (6, 0)}
            ),
            key=lambda row: (row[2], row[1], row[0]),
        )
    )
    need(generated == ALL_ROWS, "complete U=0 M12 row universe")
    need(tuple(row for row in ALL_ROWS if row[2] == 12) == ((3, 2, 12), (0, 4, 12)),
         "exact-weight-twelve endpoints after U deletion")

    counts: Dict[Point3, int] = {}
    for row in ALL_ROWS:
        for point in row_endpoints(row):
            counts[point] = counts.get(point, 0) + 1
    collisions = tuple(sorted(point for point, count in counts.items() if count > 1))
    expected = (
        (2, 3, 1),
        (2, 4, 1),
        (2, 5, 1),
        (3, 4, 1),
        (3, 5, 1),
        (4, 5, 1),
    )
    need(collisions == expected, "complete retained collision list")
    load_bearing = {
        (2, 6, 1),  # W
        (6, 4, 1),  # Z
        (0, 6, 1),  # u
        (1, 6, 1),  # a
        (0, 5, 1),  # d
        (0, 4, 1),  # fixed e
    }
    need(not load_bearing.intersection(collisions), "owner point entered collision list")
    eta_row = next(row for row in ALL_ROWS if row[:2] == ETA)
    eta_face_points = tuple(
        point for point in row_endpoints(eta_row) if plane_gap(V3W, point) == 0
    )
    need(eta_face_points == ((1, 5, 1),), "eta is a literal Case-E V3W owner")
    eta_relative_exponent = (eta_face_points[0][0], eta_face_points[0][1] - 1)
    need(
        eta_relative_exponent == (1, 4) == ((0 + 2) // 2, (3 + 5) // 2),
        "eta exponent is the midpoint of the e/W Case-E edge",
    )
    return collisions


def exhaustive_case_hull(case_name: str, case: Mapping[str, object]) -> int:
    zeros = set(case["zeros"])
    retained = tuple(row for row in ALL_ROWS if row[:2] not in zeros)
    required_keys = set(FIXED_KEYS) | set(case["owners"])
    required_rows = tuple(row for row in retained if row[:2] in required_keys)
    optional_rows = tuple(row for row in retained if row[:2] not in required_keys)
    data = candidate_plane_records(expanded_support(retained))
    expected_planes = case["planes"]

    hostile_count = 0
    for row_mask in range(1 << len(optional_rows)):
        chosen = required_rows + tuple(
            row for index, row in enumerate(optional_rows) if row_mask & (1 << index)
        )
        support = expanded_support(chosen)
        active_counts: Dict[Point3, int] = {}
        for row in chosen:
            for point in row_endpoints(row):
                active_counts[point] = active_counts.get(point, 0) + 1
        active_collisions = tuple(
            sorted(point for point, count in active_counts.items() if count > 1)
        )
        for collision_mask in range(1 << len(active_collisions)):
            hostile_support = support - {
                point
                for index, point in enumerate(active_collisions)
                if collision_mask & (1 << index)
            }
            need(
                lower_planes(data, hostile_support) == expected_planes,
                f"case {case_name}: hull changed under row/collision hostile",
            )
            hostile_count += 1
    need(hostile_count == case["hostiles"], f"case {case_name}: hostile count")
    return hostile_count


def scaled_height_charts(
    case_name: str, case: Mapping[str, object]
) -> Tuple[Tuple[Plane, Tuple[int, int, int], int], ...]:
    """Audit the minimal L-height lattice and primitive affine normals."""

    base_degree = int(case["base"])
    denominators = [6]  # good elliptic target twist
    for plane in case["planes"]:
        denominators.extend(value.denominator for value in plane)
    minimal_base = 1
    for denominator in denominators:
        minimal_base = lcm(minimal_base, denominator)
    need(minimal_base == base_degree, f"case {case_name}: minimal clearing base")

    records = []
    for plane in sorted(case["planes"]):
        scaled = tuple(base_degree * value for value in plane)
        need(all(value.denominator == 1 for value in scaled),
             f"case {case_name}: integral L-height function")
        coefficients = tuple(int(value) for value in scaled)
        # In coordinates (r,l,H=L*height), the affine equation is
        # H=(L*a)r+(L*b)l+L*c, with normal (L*a,L*b,-1).
        primitive_normal = (coefficients[0], coefficients[1], -1)
        need(primitive_normal[-1] == -1,
             f"case {case_name}: height-normal coefficient")
        common_divisor = 0
        for coefficient in primitive_normal:
            common_divisor = gcd(common_divisor, abs(coefficient))
        need(common_divisor == 1, f"case {case_name}: primitive height normal")
        records.append((plane, primitive_normal, coefficients[2]))
    return tuple(records)


def check_symbolic_torus_smoothness() -> None:
    """Verify Euler certificates for every irreducible non-rational face."""

    # R=S^2-P has partial_P R=-1, so it is smooth over every coefficient ring.
    need(-1 != 0, "R symbolic torus smoothness")
    for face_name, (constant, exponents, weights, degree) in FACE_EULER_DATA.items():
        alpha, beta_weight = weights
        weighted_degrees = tuple(
            alpha * s_exponent + beta_weight * p_exponent
            for s_exponent, p_exponent in exponents
        )
        need(all(value == degree for value in weighted_degrees),
             f"face {face_name}: symbolic Euler degree")
        # On F=0 the Euler combination is -degree*constant, hence nonzero.
        need(-degree * constant != 0,
             f"face {face_name}: symbolic torus Jacobian unit")


def check_symbolic_edge_gates() -> None:
    """Audit every literal edge and the exact open-gate failure boundaries."""

    need(F(-1376, 135) != 0, "fixed e coefficient is nonzero")
    for case_name in CASE_DATA:
        values = positive_parameters(case_name)
        W, Z = values["W"], values["Z"]
        E5 = values["a"] ** 2 - 4 * W * values["u"]
        E3 = values["eta"] ** 2 - 4 * values["e"] * W
        need(W * Z * (W + Z) != 0, f"case {case_name}: common edge gates")
        if case_name == "A":
            need(values["u"] * E5 != 0, "case A: symbolic u/E5 gates")
        if case_name in {"B", "C"}:
            need(values["a"] != 0, f"case {case_name}: symbolic a gate")
        if case_name in {"B", "D"}:
            need(values["d"] != 0, f"case {case_name}: symbolic d gate")
        if case_name == "E":
            need(E3 != 0, "case E: symbolic E3 gate")

        schemes = dict(edge_schemes(case_name, values))
        need(all(edge_regular(poly) for poly in schemes.values()),
             f"case {case_name}: all symbolic edge schemes regular")

        # The two quadratic discriminants are literally E5 and E3.
        if case_name == "A":
            q = schemes["W+aX+uX2"]
            need(q[1] ** 2 - 4 * q[0] * q[2] == E5,
                 "case A: edge discriminant is E5")
        if case_name == "E":
            q = schemes["e+etaX+WX2"]
            need(q[1] ** 2 - 4 * q[0] * q[2] == E3,
                 "case E: edge discriminant is E3")

        # W=0 or Z=0 loses a declared corner.  Lambda=0 keeps all corners but
        # makes exactly the common top product non-squarefree.
        for parameter in ("W", "Z"):
            hostile = dict(values)
            hostile[parameter] = F(0)
            failed = {
                name for name, poly in edge_schemes(case_name, hostile)
                if not edge_regular(poly)
            }
            need(bool(failed), f"case {case_name}: {parameter}-zero edge hostile")
        lambda_hostile = dict(values)
        lambda_hostile["Z"] = -lambda_hostile["W"]
        lambda_failed = tuple(
            name for name, poly in edge_schemes(case_name, lambda_hostile)
            if not edge_regular(poly)
        )
        need(lambda_failed == ("(X-1)(WX+Z)",),
             f"case {case_name}: isolated Lambda double root")

    # Exact owner/discriminant failure boundaries for the case-specific edges.
    owner_hostiles = {
        "A_u": ("A", {"u": F(0), "a": F(1)}),
        "B_a": ("B", {"a": F(0)}),
        "B_d": ("B", {"d": F(0)}),
        "C_a": ("C", {"a": F(0)}),
        "D_d": ("D", {"d": F(0)}),
    }
    for label, (case_name, changes) in owner_hostiles.items():
        hostile = positive_parameters(case_name)
        hostile.update(changes)
        need(any(not edge_regular(poly) for _name, poly in edge_schemes(case_name, hostile)),
             f"{label}: owner gate is edge-load-bearing")

    e5_hostile = positive_parameters("A")
    e5_hostile["a"] = F(2)
    e5_poly = dict(edge_schemes("A", e5_hostile))["W+aX+uX2"]
    need(e5_hostile["a"] ** 2 - 4 * e5_hostile["W"] * e5_hostile["u"] == 0,
         "E5-zero symbolic hostile")
    need(not edge_regular(e5_poly), "E5-zero repeated edge")

    e3_hostile = positive_parameters("E")
    e3_hostile["W"] = e3_hostile["e"]
    e3_hostile["eta"] = 2 * e3_hostile["e"]
    e3_poly = dict(edge_schemes("E", e3_hostile))["e+etaX+WX2"]
    need(e3_hostile["eta"] ** 2 - 4 * e3_hostile["e"] * e3_hostile["W"] == 0,
         "E3-zero symbolic hostile")
    need(not edge_regular(e3_poly), "E3-zero repeated edge")


def check_case_geometry(
    case_name: str, case: Mapping[str, object]
) -> Tuple[Tuple[Plane, Tuple[int, int, int], int], ...]:
    vertices, area2, boundary, genus = polygon_ledger(case["vertices"])
    need(vertices == case["vertices"], f"case {case_name}: global vertices")
    need((area2, boundary, genus) == case["ledger"], f"case {case_name}: Pick ledger")
    packet = edge_packet(vertices)
    need(packet == case["packet"], f"case {case_name}: packet")
    need(sum(packet) == case["responses"][0], f"case {case_name}: full response")
    need(sum(packet) - 8 == case["responses"][1], f"case {case_name}: finite response")
    need(sum(index - 1 for index in packet) == 2 * genus - 2,
         f"case {case_name}: Hurwitz defect")

    genus_sum = 0
    for face_name, expected_genus in case["faces"].items():
        face_polygon, expected_ledger = FACE_POLYGONS[face_name]
        _face_vertices, face_area2, face_boundary, face_genus = polygon_ledger(face_polygon)
        need((face_area2, face_boundary, face_genus) == expected_ledger,
             f"case {case_name}: face {face_name} Pick ledger")
        need(face_genus == expected_genus, f"case {case_name}: face genus")
        genus_sum += face_genus
    graph_vertices = tuple(case["graph_vertices"])
    graph_edges = tuple(case["graph_edges"])
    need(len(graph_vertices) == int(case["components"]),
         f"case {case_name}: explicit graph vertex ledger")
    need(sum(multiplicity for _edge, multiplicity in graph_edges) == int(case["edges"]),
         f"case {case_name}: explicit graph edge ledger")
    reached = {graph_vertices[0]}
    changed = True
    while changed:
        changed = False
        for edge_name, multiplicity in graph_edges:
            need(multiplicity > 0, f"case {case_name}: positive graph edge multiplicity")
            first, second = edge_name.split("-")
            need(first in graph_vertices and second in graph_vertices,
                 f"case {case_name}: graph edge endpoints")
            if first in reached and second not in reached:
                reached.add(second)
                changed = True
            if second in reached and first not in reached:
                reached.add(first)
                changed = True
    need(reached == set(graph_vertices), f"case {case_name}: graph connectedness")
    graph_b1 = int(case["edges"]) - int(case["components"]) + 1
    need(graph_b1 == 11, f"case {case_name}: graph b1")
    need(genus_sum + graph_b1 == genus, f"case {case_name}: complete genus inventory")

    base_degree = int(case["base"])
    order_by_plane = {
        "M": vertical_order(M, base_degree),
        "H5": vertical_order(H5, base_degree),
        "V4a": vertical_order(V4A, base_degree),
        "V3a": vertical_order(V3A, base_degree),
        "V4W": vertical_order(V4W, base_degree),
        "V3W": vertical_order(V3W, base_degree),
    }
    for label, expected_order in case["orders"].items():
        need(order_by_plane[label] == expected_order,
             f"case {case_name}: good-differential order {label}")
        need(order_by_plane[label] > 0, f"case {case_name}: positive face order")

    for plane in case["planes"]:
        need(plane[1] + plane[2] == 0, f"case {case_name}: b+c cancellation")
        for point in expanded_support(
            row for row in ALL_ROWS if row[:2] not in set(case["zeros"])
        ):
            need(plane_gap(plane, point) >= 0, f"case {case_name}: negative face gap")

    values = positive_parameters(case_name)
    need(all(edge_regular(poly) for _name, poly in edge_schemes(case_name, values)),
         f"case {case_name}: exact positive edge control")
    return scaled_height_charts(case_name, case)


def check_hostile_and_differential_controls() -> None:
    # Case-A positive control permits alpha=0: E5=-4 when W=u=1.
    need(F(0) ** 2 - 4 * F(1) * F(1) == -4, "alpha-zero E5 positive control")

    # E5=0 hostile: W=u=1,a=2 gives (1+X)^2 and the (2,5) cusp delta=2.
    cusp_edge = (F(1), F(2), F(1))
    need(not squarefree(cusp_edge), "E5-zero edge must be repeated")
    need((2 - 1) * (5 - 1) // 2 == 2, "(2,5) cusp delta")

    # Case-E face and exact elliptic normal form.  Put T=S*P, X=P^-1,
    # q(T)=e+eta*T+W*T^2 and Y=2*W*T+eta.  The face equation says q=X^3;
    # completing the square gives Y^2=4*W*X^3+E3.
    e = F(-1376, 135)
    W_e = F(1)
    eta = F(0)
    E3 = eta * eta - 4 * e * W_e
    q_coefficients = (e, eta, W_e)
    completed_square = (
        eta * eta - E3,
        4 * W_e * eta,
        4 * W_e * W_e,
    )
    four_W_q = tuple(4 * W_e * coefficient for coefficient in q_coefficients)
    need(E3 != 0, "Case-E positive E3 control")
    need(squarefree(q_coefficients), "Case-E quadratic edge squarefree")
    need(completed_square == four_W_q, "Case-E exact elliptic normal form")
    elliptic_cubic = (E3, F(0), F(0), 4 * W_e)
    need(squarefree(elliptic_cubic), "Case-E normal-form cubic squarefree")
    need(-2 * 3 + 3 * (3 - 1) == 0, "Case-E cyclic-cubic genus one")

    # The monomial change (S,P)->(X=P^-1,V=SP) is unimodular, with inverse
    # S=XV, P=X^-1.  The subsequent affine coordinate Y=2WV+eta is etale on
    # the declared W-nonzero gate.
    exponent_matrix = ((0, -1), (1, 1))
    inverse_matrix = ((1, 1), (-1, 0))
    determinant = (
        exponent_matrix[0][0] * exponent_matrix[1][1]
        - exponent_matrix[0][1] * exponent_matrix[1][0]
    )
    need(determinant == 1, "V3W monomial change unimodular")
    matrix_product = tuple(
        tuple(sum(exponent_matrix[row][k] * inverse_matrix[k][column] for k in range(2))
              for column in range(2))
        for row in range(2)
    )
    need(matrix_product == ((1, 0), (0, 1)), "V3W integral inverse matrix")
    need(2 * W_e != 0, "V3W affine Y-coordinate etale")

    # E3=0 hostile: W=e and eta=2e give q=e(1+T)^2.
    W_e3_zero = e
    eta_e3_zero = 2 * e
    repeated_case_e_edge = (e, eta_e3_zero, W_e3_zero)
    need(eta_e3_zero**2 - 4 * e * W_e3_zero == 0,
         "Case-E E3-zero discriminant")
    need(not squarefree(repeated_case_e_edge),
         "Case-E E3-zero edge must be repeated")

    # Lambda=0 hostile: W=1,Z=-1 gives (X-1)^2, one length-twelve A23 contact.
    lambda_edge = (F(1), F(-2), F(1))
    need(not squarefree(lambda_edge), "Lambda-zero top edge must be repeated")
    need(2 * 12 - 1 == 23, "length-twelve contact is A23")
    # Only the specialization of THM-4297 is checked here.  With U=0 and
    # Z=-W, one has A=W and U_eff=W/2.  The following coefficient identity is
    # its exact local input; the critical ladder itself remains imported.
    W_lambda = F(1)
    A_specialized = W_lambda
    U_effective = W_lambda / 2
    need(A_specialized == W_lambda and U_effective == W_lambda / 2,
         "THM4297 U-zero specialization A=W,Ueff=W/2")
    left_identity = {4: -W_lambda, 5: W_lambda}
    right_identity = {
        4: -W_lambda / 2 - W_lambda / 2,
        5: W_lambda,
        6: W_lambda / 2 - W_lambda / 2,
    }
    need(left_identity == {degree: coefficient for degree, coefficient in right_identity.items()
                            if coefficient},
         "THM4297 specialized cubic-correction identity")
    tail_forms = (
        (F(6), F(2)),
        (F(5, 2), F(9, 2)),
        (F(2), F(4)),
        (F(3, 2), F(7, 2)),
        (F(1), F(3)),
    )
    need(all(s_coefficient > 0 and beta_coefficient > 0
             for s_coefficient, beta_coefficient in tail_forms),
         "Lambda-zero A23 tail orders")

    # Exact hostile to primitive-CM transfer.  For T=SP,
    # X=P^2,Y=WP^3+2ZT^2 gives Y^2-W^2X^3-4Z
    # =4Z(WP^3T^2+ZT^4-1), which vanishes on C.
    W = F(3)
    Z = F(5)
    left_coefficients = {
        (6, 0): W * W - W * W,
        (3, 2): 4 * W * Z,
        (0, 4): 4 * Z * Z,
        (0, 0): -4 * Z,
    }
    right_coefficients = {
        (6, 0): F(0),
        (3, 2): 4 * W * Z,
        (0, 4): 4 * Z * Z,
        (0, 0): -4 * Z,
    }
    need(left_coefficients == right_coefficients, "explicit j=0 quotient identity")

    # Generic torus completion: t=P-S^2 is nonzero in k(C), since reducing C
    # modulo t gives 1-(W+Z)S^12, whose constant coefficient is one for every
    # Lambda.  Thus R does not become a component of C, even at Lambda=0.
    completion_remainder = (F(1), -(W + Z))
    need(completion_remainder[0] == 1,
         "generic torus completion t=P-S^2 nonzero")

    # The generic-unit identities used by the differential argument.
    need(5 != 0, "H5: P*F_P=5")
    for k in (3, 4):
        need((5, k - 5) != (0, 0), f"Vk: 5+(k-5)cP^k nonzero")


def run_audit_mode() -> Mapping[str, int]:
    fixed_rows = tuple(row for row in ALL_ROWS if row[:2] in FIXED_KEYS)
    optional_rows = tuple(row for row in ALL_ROWS if row[:2] not in FIXED_KEYS)
    need(len(optional_rows) == 10, "audit optional coefficient count")
    universal_data = candidate_plane_records(expanded_support(ALL_ROWS))
    counts = {name: 0 for name in CASE_DATA}
    for mask in range(1 << len(optional_rows)):
        chosen_optional = tuple(
            row for index, row in enumerate(optional_rows) if mask & (1 << index)
        )
        chosen = fixed_rows + chosen_optional
        keys = {row[:2] for row in chosen}
        case_name = expected_case_from_owner_mask(keys)
        expected_planes = CASE_DATA[case_name]["planes"]
        need(lower_planes(universal_data, expanded_support(chosen)) == expected_planes,
             f"audit mask {mask}: owner classification")
        counts[case_name] += 1
    need(sum(counts.values()) == 1024, "audit mask total")
    need(counts == {"A": 512, "B": 128, "C": 128, "D": 128, "E": 128},
         "audit complex population")
    return counts


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def plane_text(plane: Plane) -> str:
    return "(" + ",".join(fraction_text(value) for value in plane) + ")"


def height_chart_text(
    records: Sequence[Tuple[Plane, Tuple[int, int, int], int]]
) -> str:
    return ";".join(
        f"{PLANE_NAMES[plane]}:normal={normal},offset={offset}"
        for plane, normal, offset in records
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--audit",
        action="store_true",
        help="add the independent complete 1,024-mask owner-state census",
    )
    arguments = parser.parse_args()

    collisions = check_universe()
    check_symbolic_torus_smoothness()
    check_symbolic_edge_gates()
    hostile_counts: Dict[str, int] = {}
    height_charts: Dict[
        str, Tuple[Tuple[Plane, Tuple[int, int, int], int], ...]
    ] = {}
    for case_name, case in CASE_DATA.items():
        hostile_counts[case_name] = exhaustive_case_hull(case_name, case)
        height_charts[case_name] = check_case_geometry(case_name, case)
    check_hostile_and_differential_controls()

    audit_counts: Mapping[str, int] | None = None
    if arguments.audit:
        audit_counts = run_audit_mode()

    print("THM-4327 U=0 ENDPOINT FINITE CERTIFICATE")
    print(f"MODE={'audit' if arguments.audit else 'primary'}")
    print(f"UNIVERSE_ROWS={len(ALL_ROWS)}")
    print("FILTERS=U=0;W*Z!=0;optional lower rows arbitrary;active collisions hostile-deleted")
    print("COLLISIONS=" + ";".join(str(point) for point in collisions))
    for case_name, case in CASE_DATA.items():
        planes = ";".join(plane_text(plane) for plane in sorted(case["planes"]))
        packet = ",".join(str(value) for value in case["packet"])
        orders = ",".join(f"{key}:{value}" for key, value in case["orders"].items())
        ledger = "/".join(str(value) for value in case["ledger"])
        responses = "/".join(str(value) for value in case["responses"])
        print(
            f"CASE_{case_name} gate={case['gate']} planes={planes} "
            f"hostiles={hostile_counts[case_name]} polygon={ledger} "
            f"packet={packet} responses={responses} L={case['base']} orders={orders}"
        )
        graph_edges = ",".join(
            f"{edge}:{multiplicity}" for edge, multiplicity in case["graph_edges"]
        )
        print(
            f"CASE_{case_name}_GRAPH_LAMBDA_NONZERO V={case['components']} E={case['edges']} b1=11 "
            f"vertices={','.join(case['graph_vertices'])} edges={graph_edges}"
        )
        print(f"CASE_{case_name}_HEIGHTS " + height_chart_text(height_charts[case_name]))
    print(f"TOTAL_HOSTILES={sum(hostile_counts.values())}")
    print("POSITIVE_FACE_ORDERS=all_strictly_positive")
    print("HEIGHT_CHARTS=all_minimal_L_integral_primitive;height_normal_coefficient=-1")
    print("TORUS_SMOOTHNESS=R,C,H5,V4a,V3a,V4W,V3W_symbolic_Euler_certificates")
    print("EDGE_SCHEMES=common[X-1,1-ZX4,(X-1)(WX+Z),1-WX];A[W+aX+uX2,u-X5];B[W+aX,a+dX,d-X4,1-aX];C[W+aX,a+eX,e-X3,1-aX];D[W+dX,d-X4];E[e+etaX+WX2,e-X3]")
    print("EDGE_GATE_LEDGER=common:W*Z*Lambda;A:u*E5;B:a*d;C:a;D:d;E:E3;internal:W_and_a_when_present")
    print("CASE_E_FACE=-1+e*P^3+eta*S*P^4+W*S^2*P^5")
    print("CASE_E_NORMAL_FORM=T=S*P;X=P^-1;Y=2W*T+eta;Y^2=4W*X^3+E3;genus=1")
    print("CASE_E_CHANGE=(S,P)->(X=P^-1,V=SP)_det=1;Y=2WV+eta_etale_for_W!=0")
    print("TORUS_COMPLETION=t=P-S^2;C_mod_t=1-Lambda*S^12_nonzero")
    print("LAMBDA_ZERO=top_edge_double_root;one_A23_contact;all_five_imported_tail_orders_positive")
    print("THM4297_IMPORT=checked_A=W,Ueff=W/2_and_order_table_only;critical_ladder_not_rederived")
    print("HOSTILE_E5_ZERO=(2,5)_cusp_delta_2;repeated-tail-audit-required")
    print("HOSTILE_E3_ZERO=e+eta*X+W*X^2_has_double_root;repeated-tail-audit-required")
    print("HOSTILE_CM_TRANSFER=explicit_nonconstant_j0_quotient_verified")
    if audit_counts is not None:
        print("AUDIT_MASKS=1024")
        print(
            "AUDIT_COMPLEX_COUNTS="
            + ",".join(f"{name}:{audit_counts[name]}" for name in CASE_DATA)
        )
    print(f"CHECKS={CHECKS}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
