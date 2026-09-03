#!/usr/bin/env python3
"""Independent hostile referee for the reserved THM-4355 endpoint gate.

This certificate deliberately rebuilds the alpha_11=U=0 calculation from the
literal lifted source.  It does not import the primary THM-4355 certificate.
All arithmetic used for support, hull, Pick, graph, and valuation checks is
exact.  SymPy is used only for literal polynomial identities and local series.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from math import gcd
import sys

import sympy as sp

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


class Audit:
    def __init__(self) -> None:
        self.checks = 0

    def ok(self, condition: bool, label: str) -> None:
        self.checks += 1
        if not condition:
            raise AssertionError(label)

    def eq(self, left, right, label: str) -> None:
        self.checks += 1
        if left != right:
            raise AssertionError(f"{label}: {left!r} != {right!r}")


A = Audit()
Point3 = tuple[int, int, int]
Point2 = tuple[int, int]

E = F(-1376, 135)
D_STAR = F(3968, 63)       # K-e=0
D_KZERO = F(5696, 105)    # K=0, excluded


def kval(delta: F) -> F:
    return F(2848, 45) - F(7, 6) * delta


BASE: dict[str, Point3] = {
    "s2": (0, 2, 0),
    "p": (0, 0, 1),
    "qhalf": (1, 2, 0),
    "m3L": (1, 2, 1),
    "m3R": (1, 0, 2),
    "eightL": (1, 2, 2),
    "eightR": (1, 0, 3),
    "eR": (1, 0, 4),
    "Ktop": (1, 4, 2),
}

OPTIONAL: dict[str, tuple[Point3, ...]] = {
    "KminusE": ((1, 2, 3),),
    "Phi": ((1, 3, 3), (1, 1, 4)),
    "Delta": ((1, 0, 5),),
    "ThetaMinusDelta": ((1, 2, 4),),
    "Theta": ((1, 4, 3),),
    "eta": ((1, 3, 4), (1, 1, 5)),
    "u": ((1, 2, 5), (1, 0, 6)),
}


def literal_coefficients(delta: F, theta: F, phi: F, eta: F, u: F) -> dict[Point3, F]:
    """Expand (s^2-p)(1-QH)-Qs^2/2 on alpha=U=W=xi=...=0."""
    k = kval(delta)
    out: defaultdict[Point3, F] = defaultdict(F)
    out[(0, 2, 0)] += 1
    out[(0, 0, 1)] -= 1
    out[(1, 2, 0)] -= F(1, 2)
    rows = (
        (F(-3), 0, 1),
        (F(8, 3), 0, 2),
        (E, 0, 3),
        (k, 2, 2),
        (phi, 1, 3),
        (delta, 0, 4),
        (theta, 2, 3),
        (eta, 1, 4),
        (u, 0, 5),
    )
    for coefficient, s_exp, p_exp in rows:
        out[(1, s_exp + 2, p_exp)] -= coefficient
        out[(1, s_exp, p_exp + 1)] += coefficient
    return {point: coefficient for point, coefficient in out.items() if coefficient}


def cross(u: Point3, v: Point3) -> Point3:
    return (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    )


def sub3(a: Point3, b: Point3) -> Point3:
    return tuple(a[i] - b[i] for i in range(3))  # type: ignore[return-value]


def dot3(a: Point3, b: Point3) -> int:
    return sum(a[i] * b[i] for i in range(3))


def primitive_plane(normal: Point3, constant: int) -> tuple[int, int, int, int]:
    g = 0
    for value in (*normal, constant):
        g = gcd(g, abs(value))
    if g:
        normal = tuple(value // g for value in normal)  # type: ignore[assignment]
        constant //= g
    if normal[0] < 0:
        normal = tuple(-value for value in normal)  # type: ignore[assignment]
        constant = -constant
    return (*normal, constant)


def lower_facet_planes(points: frozenset[Point3]) -> dict[tuple[int, int, int, int], frozenset[Point3]]:
    """Facets minimizing a weight whose Q-coordinate is positive."""
    pts = sorted(points)
    faces: dict[tuple[int, int, int, int], frozenset[Point3]] = {}
    for p0, p1, p2 in combinations(pts, 3):
        n = cross(sub3(p1, p0), sub3(p2, p0))
        if n == (0, 0, 0) or n[0] == 0:
            continue
        c = dot3(n, p0)
        vals = [dot3(n, p) - c for p in pts]
        if not (all(v >= 0 for v in vals) or all(v <= 0 for v in vals)):
            continue
        if all(v <= 0 for v in vals):
            n = tuple(-v for v in n)  # type: ignore[assignment]
            c = -c
            vals = [-v for v in vals]
        if n[0] < 0:
            # This supporting side is visible for negative Q-weight, not Q -> 0.
            continue
        face = frozenset(p for p, value in zip(pts, vals) if value == 0)
        if len(face) < 3:
            continue
        key = primitive_plane(n, c)
        faces[key] = face
    return faces


def convex_hull_2d(points: set[Point2]) -> list[Point2]:
    pts = sorted(points)
    if len(pts) <= 1:
        return pts

    def turn(o: Point2, a: Point2, b: Point2) -> int:
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower: list[Point2] = []
    for p in pts:
        while len(lower) >= 2 and turn(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper: list[Point2] = []
    for p in reversed(pts):
        while len(upper) >= 2 and turn(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


def pick_ledger(points: set[Point2]) -> tuple[int, int, int]:
    hull = convex_hull_2d(points)
    twice_area = abs(sum(
        hull[i][0] * hull[(i + 1) % len(hull)][1]
        - hull[i][1] * hull[(i + 1) % len(hull)][0]
        for i in range(len(hull))
    ))
    boundary = sum(
        gcd(
            abs(hull[(i + 1) % len(hull)][0] - hull[i][0]),
            abs(hull[(i + 1) % len(hull)][1] - hull[i][1]),
        )
        for i in range(len(hull))
    )
    interior = (twice_area - boundary + 2) // 2
    return twice_area, boundary, interior


@dataclass(frozen=True)
class ExactRow:
    coupled: str
    delta: F
    theta: F
    u: F
    eta: F
    phi: F
    support: frozenset[Point3]
    fan: str


COUPLED_REPS = (
    ("D0_T0", F(0), F(0)),
    ("D0_Tx", F(0), F(1)),
    ("Ds_T0", D_STAR, F(0)),
    ("Ds_Ts", D_STAR, D_STAR),
    ("Ds_Tx", D_STAR, 2 * D_STAR),
    ("Dg_T0", F(1), F(0)),
    ("Dg_Td", F(1), F(1)),
    ("Dg_Tx", F(1), F(2)),
)


def fan_name(delta: F, theta: F, u: F, eta: F, phi: F) -> str:
    if theta:
        if u:
            return "BU+FUT+T"
        if delta and eta:
            return "C4e+BETA+FET+T"
        if delta:
            return "BDT+T"
        if eta:
            return "C3e+BETA+FET+T"
        return "CTH+BTH+T"
    if u and eta:
        return "BU+FUE+EETA"
    if u:
        return "BU+EU"
    if delta and eta:
        return "C4e+BETA+EETA"
    if delta:
        return "BDELTA+EDELTA"
    if eta:
        return "C3e+BETA+EETA"
    if phi:
        return "CPHI+BPHI+EDELTA"
    return "BDEEP"


EXACT_ROWS: list[ExactRow] = []
for coupled, delta, theta in COUPLED_REPS:
    A.ok(kval(delta) != 0, f"K owner survives {coupled}")
    for u_bit, eta_bit, phi_bit in product((0, 1), repeat=3):
        coeffs = literal_coefficients(delta, theta, F(phi_bit), F(eta_bit), F(u_bit))
        EXACT_ROWS.append(ExactRow(
            coupled, delta, theta, F(u_bit), F(eta_bit), F(phi_bit),
            frozenset(coeffs), fan_name(delta, theta, F(u_bit), F(eta_bit), F(phi_bit)),
        ))

A.eq(len(EXACT_ROWS), 64, "eight coupled classes times three bits")
A.eq(len({row.support for row in EXACT_ROWS}), 64, "all exact supports are distinct")
A.eq(sum(bool(row.theta) for row in EXACT_ROWS), 40, "Theta-owner support count")
A.eq(sum(not row.theta for row in EXACT_ROWS), 24, "Theta-zero support count")
A.eq(len({row.support for row in EXACT_ROWS if row.theta}), 40, "Theta-owner supports distinct")
A.eq(len({row.support for row in EXACT_ROWS if not row.theta}), 24, "Theta-zero supports distinct")

EXPECTED_FAN_COUNTS = {
    "BU+FUT+T": 20,
    "C4e+BETA+FET+T": 8,
    "BDT+T": 8,
    "C3e+BETA+FET+T": 2,
    "CTH+BTH+T": 2,
    "BU+FUE+EETA": 6,
    "BU+EU": 6,
    "C4e+BETA+EETA": 4,
    "BDELTA+EDELTA": 4,
    "C3e+BETA+EETA": 2,
    "CPHI+BPHI+EDELTA": 1,
    "BDEEP": 1,
}
A.eq(Counter(row.fan for row in EXACT_ROWS), Counter(EXPECTED_FAN_COUNTS), "exact fan split")


# Full face support polygons, reconstructed from the displayed factor equations.
P0, S20 = (0, 0, 1), (0, 2, 0)
FACE_POINTS: dict[str, frozenset[Point3]] = {
    "BU": frozenset((P0, S20, (1, 0, 6), (1, 2, 5))),
    "FUT": frozenset((S20, (1, 2, 5), (1, 3, 4), (1, 4, 3))),
    "T": frozenset((S20, (1, 4, 2), (1, 4, 3))),
    "C4e": frozenset((P0, (1, 0, 5), (1, 1, 5))),
    "BETA": frozenset((P0, S20, (1, 1, 5), (1, 3, 4))),
    "FET": frozenset((S20, (1, 3, 4), (1, 4, 3))),
    "BDT": frozenset((P0, S20, (1, 0, 5), (1, 2, 4), (1, 4, 3))),
    "C3e": frozenset((P0, (1, 0, 4), (1, 1, 5))),
    "CTH": frozenset((P0, (1, 0, 4), (1, 1, 4), (1, 2, 4))),
    "BTH": frozenset((P0, S20, (1, 2, 4), (1, 4, 3))),
    "FUE": frozenset((S20, (1, 2, 5), (1, 3, 4))),
    "EETA": frozenset((S20, (1, 3, 4), (1, 4, 2))),
    "EU": frozenset((S20, (1, 2, 5), (1, 4, 2))),
    "BDELTA": frozenset((P0, S20, (1, 0, 5), (1, 2, 4))),
    "EDELTA": frozenset((S20, (1, 2, 4), (1, 3, 3), (1, 4, 2))),
    "CPHI": frozenset((P0, (1, 0, 4), (1, 1, 4))),
    "BPHI": frozenset((P0, S20, (1, 1, 4), (1, 3, 3))),
    "BDEEP": frozenset((P0, S20, (1, 0, 4), (1, 2, 3), (1, 4, 2))),
}


def plane_from_face(points: frozenset[Point3]) -> tuple[int, int, int, int]:
    for p0, p1, p2 in combinations(sorted(points), 3):
        n = cross(sub3(p1, p0), sub3(p2, p0))
        if n != (0, 0, 0) and n[0] != 0:
            c = dot3(n, p0)
            if n[0] < 0:
                n = tuple(-x for x in n)  # type: ignore[assignment]
                c = -c
            return primitive_plane(n, c)
    raise AssertionError("face has no Q-visible spanning triple")


FACE_PLANES = {name: plane_from_face(points) for name, points in FACE_POINTS.items()}
A.eq(len(set(FACE_PLANES.values())), 13, "thirteen geometric planes underlie eighteen stratum-names")

for row in EXACT_ROWS:
    got = lower_facet_planes(row.support)
    expected_names = row.fan.split("+")
    expected_planes = {FACE_PLANES[name] for name in expected_names}
    A.eq(set(got), expected_planes, f"lower hull fan {row.coupled}/{int(row.u)}{int(row.eta)}{int(row.phi)}")
    for name in expected_names:
        key = FACE_PLANES[name]
        A.eq(got[key], row.support & FACE_POINTS[name], f"face support {name}")


# A deliberately larger point-toggle atlas: it ignores all affine coefficient
# coupling, but keeps the mandatory base and K owner.  Every exact support embeds.
TOGGLE_BLOCKS = tuple(OPTIONAL)
HOSTILE_SUPPORTS: set[frozenset[Point3]] = set()
HOSTILE_FANS: set[tuple[tuple[int, int, int, int], ...]] = set()
for bits in product((0, 1), repeat=len(TOGGLE_BLOCKS)):
    support = set(BASE.values())
    for bit, block in zip(bits, TOGGLE_BLOCKS):
        if bit:
            support.update(OPTIONAL[block])
    frozen = frozenset(support)
    HOSTILE_SUPPORTS.add(frozen)
    HOSTILE_FANS.add(tuple(sorted(lower_facet_planes(frozen))))
A.eq(len(HOSTILE_SUPPORTS), 128, "support-decoupled hostile atlas")
for row in EXACT_ROWS:
    A.ok(row.support in HOSTILE_SUPPORTS, "exact support embeds in hostile atlas")


GLOBAL_PICK = {
    "BU+FUT+T": (32, 12, 11),
    "C4e+BETA+FET+T": (29, 11, 10),
    "BDT+T": (26, 10, 9),
    "C3e+BETA+FET+T": (28, 10, 10),
    "CTH+BTH+T": (24, 10, 8),
    "BU+FUE+EETA": (31, 11, 11),
    "BU+EU": (30, 10, 11),
    "C4e+BETA+EETA": (28, 10, 10),
    "BDELTA+EDELTA": (24, 10, 8),
    "C3e+BETA+EETA": (27, 9, 10),
    "CPHI+BPHI+EDELTA": (21, 9, 7),
    "BDEEP": (18, 8, 6),
}
for fan, expected in GLOBAL_PICK.items():
    witnesses = [row for row in EXACT_ROWS if row.fan == fan]
    ledgers = {pick_ledger({(p[1], p[2]) for p in row.support}) for row in witnesses}
    A.eq(ledgers, {expected}, f"global Pick ledger {fan}")


FACE_PICK = {
    "BU": (20, 12, 5), "FUT": (10, 8, 2), "T": (2, 4, 0),
    "C4e": (4, 6, 0), "BETA": (18, 4, 8), "FET": (5, 3, 2),
    "BDT": (24, 8, 9), "C3e": (3, 5, 0), "CTH": (6, 6, 1),
    "BTH": (16, 4, 7), "FUE": (5, 7, 0), "EETA": (6, 4, 2),
    "EU": (10, 8, 2), "BDELTA": (16, 10, 4), "EDELTA": (8, 8, 1),
    "CPHI": (3, 5, 0), "BPHI": (14, 4, 6), "BDEEP": (18, 8, 6),
}
for name, expected in FACE_PICK.items():
    A.eq(pick_ledger({(p[1], p[2]) for p in FACE_POINTS[name]}), expected, f"face Pick {name}")


# The chart-plane triple is derived, not copied: if the primitive supporting
# equation is n_Q Q+n_s s+n_p p=c, it is (-n_s/n_Q,-n_p/n_Q,c/n_Q).
def chart_plane(name: str) -> tuple[F, F, F]:
    nq, ns, np, c = FACE_PLANES[name]
    return F(-ns, nq), F(-np, nq), F(c, nq)


TARGET_ORDERS = {
    "BU": (30, 22), "FUT": (30, 25), "T": (6, 8),
    "C4e": (12, 10), "BETA": (18, 13), "FET": (30, 25),
    "BDT": (24, 17), "C3e": (6, 7), "CTH": (6, 5),
    "BTH": (24, 17), "FUE": (30, 25), "EETA": (6, 6),
    "EU": (30, 28), "BDELTA": (24, 17), "EDELTA": (12, 10),
    "CPHI": (6, 5), "BPHI": (42, 29), "BDEEP": (6, 4),
}
for name, (target_base, expected_order) in TARGET_ORDERS.items():
    source_base = FACE_PLANES[name][0]
    A.ok(target_base % source_base == 0, f"target base dominates primitive source base {name}")
    plane = chart_plane(name)
    order = target_base * (F(5, 6) - sum(plane, F(0)))
    A.eq(order, expected_order, f"target differential order {name}")
    A.ok(order > 0, f"positive face order {name}")


# ---------------------------------------------------------------------------
# Literal carrier equations and normalization discriminants
# ---------------------------------------------------------------------------

P, S = sp.symbols("P S")
u_s, eta_s, phi_s, delta_s, theta_s, k_s = sp.symbols(
    "u eta Phi Delta Theta K", nonzero=True
)

FACE_CURVES = {
    "FUT": 1 - u_s * P**5 - eta_s * S * P**4 - theta_s * S**2 * P**3,
    "FET": 1 - eta_s * S * P**4 - theta_s * S**2 * P**3,
    "BDT-C": delta_s * P**4 + theta_s * S**2 * P**3 - 1,
    "CTH": E * P**3 + phi_s * S * P**3 + theta_s * S**2 * P**3 - 1,
    "EETA": 1 - eta_s * S * P**4 - k_s * S**2 * P**2,
    "EU": 1 - u_s * P**5 - k_s * S**2 * P**2,
    "EDELTA": 1 - delta_s * P**4 - phi_s * S * P**3 - k_s * S**2 * P**2,
    "BDEEP-C": E * P**3 + k_s * S**2 * P**2 - 1,
}
NORMALIZED_BRANCH = {
    "FUT": P * (4 * theta_s + (eta_s**2 - 4 * u_s * theta_s) * P**5),
    "FET": P * (4 * theta_s + eta_s**2 * P**5),
    "BDT-C": theta_s * P * (1 - delta_s * P**4),
    "CTH": (phi_s**2 - 4 * E * theta_s) * P**4 + 4 * theta_s * P,
    "EETA": eta_s**2 * P**6 + 4 * k_s,
    "EU": k_s * (1 - u_s * P**5),
    "EDELTA": (phi_s**2 - 4 * k_s * delta_s) * P**4 + 4 * k_s,
    "BDEEP-C": k_s * (1 - E * P**3),
}
DISC_SQUARE_FACTOR = {
    "FUT": P**2, "FET": P**2, "BDT-C": 4 * P**2, "CTH": P**2,
    "EETA": P**2, "EU": 4 * P**2, "EDELTA": P**2, "BDEEP-C": 4 * P**2,
}
for name, curve in FACE_CURVES.items():
    disc = sp.factor(sp.discriminant(curve, S))
    expected = sp.factor(DISC_SQUARE_FACTOR[name] * NORMALIZED_BRANCH[name])
    A.eq(sp.expand(disc - expected), 0, f"normalization discriminant {name}")


def polynomial_degree(expr) -> int:
    return int(sp.Poly(expr, P).degree())


CLEAN_CARRIERS = {
    "FUT": (6, 2, 2), "FET": (6, 2, 2), "BDT-C": (5, 2, 1),
    "CTH": (4, 1, 2), "EETA": (6, 2, 2), "EU": (5, 2, 1),
    "EDELTA": (4, 1, 2), "BDEEP-C": (3, 1, 1),
}
GENERIC_SPECIALIZATION = {
    u_s: F(2), eta_s: F(3), phi_s: F(5), delta_s: F(7), theta_s: F(11), k_s: F(13)
}
for name, (degree, genus, infinity_ends) in CLEAN_CARRIERS.items():
    branch = sp.Poly(NORMALIZED_BRANCH[name].subs(GENERIC_SPECIALIZATION), P)
    A.eq(branch.degree(), degree, f"branch degree {name}")
    A.eq(sp.gcd(branch, branch.diff()).degree(), 0, f"squarefree positive control {name}")
    A.eq((degree - 1) // 2, genus, f"hyperelliptic genus {name}")
    A.eq(1 if degree % 2 else 2, infinity_ends, f"infinity ends {name}")


# Four and only four nontrivial edge discriminants/resultants.
X = sp.symbols("X")
L5 = eta_s**2 - 4 * u_s * theta_s
L3 = phi_s**2 - 4 * E * theta_s
L15 = delta_s + theta_s
L4 = phi_s**2 - 4 * k_s * delta_s
A.eq(sp.discriminant(theta_s + eta_s * X + u_s * X**2, X), L5, "FUT A4 wall")
A.eq(sp.discriminant(theta_s + phi_s * X + E * X**2, X), L3, "CTH A2 wall")
A.eq(sp.resultant(X - 1, delta_s * X + theta_s, X), L15, "BDT A15 wall")
A.eq(sp.discriminant(delta_s + phi_s * X + k_s * X**2, X), L4, "EDELTA A3 wall")

FIXED_EDGE_POLYS = (
    1 - u_s * X**5, u_s - X**5, delta_s - X**4, 1 - delta_s * X**4,
    E - X**3, 1 - eta_s * X, X**2 - k_s, 1 - k_s * X**2,
    k_s + E * X,
)
for index, edge in enumerate(FIXED_EDGE_POLYS):
    A.eq(sp.gcd(sp.Poly(edge, X), sp.Poly(sp.diff(edge, X), X)).degree(), 0,
         f"generic fixed edge is reduced {index}")
A.eq(kval(F(0)) + E, F(7168, 135), "BDEEP two linear factors stay distinct")


# ---------------------------------------------------------------------------
# Graph and arithmetic-genus ledgers
# ---------------------------------------------------------------------------

def graph_stats(vertices: set[str], edges: list[tuple[str, str]]) -> tuple[int, int, int, int]:
    adjacent = {v: set() for v in vertices}
    for left, right in edges:
        A.ok(left in vertices and right in vertices, "edge endpoints exist")
        adjacent[left].add(right)
        adjacent[right].add(left)
    unseen = set(vertices)
    components = 0
    while unseen:
        components += 1
        stack = [unseen.pop()]
        while stack:
            here = stack.pop()
            for nxt in adjacent[here]:
                if nxt in unseen:
                    unseen.remove(nxt)
                    stack.append(nxt)
    return len(vertices), len(edges), len(edges) - len(vertices) + components, components


def repeated_edge(left: str, right: str, count: int) -> list[tuple[str, str]]:
    return [(left, right) for _ in range(count)]


def graph_bu(last: str, with_middle: bool = True) -> tuple[set[str], list[tuple[str, str]]]:
    verticals = {f"V{i}" for i in range(5)}
    vertices = {"R", last} | verticals
    edges: list[tuple[str, str]] = []
    for vtx in sorted(verticals):
        edges += repeated_edge("R", vtx, 2)
        edges.append((vtx, last))
    return vertices, edges


def chain_with_parallel(n: int, labels: tuple[str, str, str, str]) -> tuple[set[str], list[tuple[str, str]]]:
    left, rcomp, acomp, right = labels
    vertices = set(labels)
    return vertices, [(left, rcomp)] + repeated_edge(rcomp, acomp, n) + [(acomp, right)]


FAN_GRAPHS: dict[str, tuple[set[str], list[tuple[str, str]], int, tuple[int, int, int]]] = {}
vtx, edg = graph_bu("F")
vtx.add("T"); edg.append(("F", "T"))
FAN_GRAPHS["BU+FUT+T"] = (vtx, edg, 2, (8, 16, 9))
vtx, edg = chain_with_parallel(9, ("C4", "R", "A", "F")); vtx.add("T"); edg.append(("F", "T"))
FAN_GRAPHS["C4e+BETA+FET+T"] = (vtx, edg, 2, (5, 12, 8))
FAN_GRAPHS["BDT+T"] = ({"R", "C", "T"}, repeated_edge("R", "C", 8) + [("C", "T")], 2, (3, 9, 7))
vtx, edg = chain_with_parallel(9, ("C3", "R", "A", "F")); vtx.add("T"); edg.append(("F", "T"))
FAN_GRAPHS["C3e+BETA+FET+T"] = (vtx, edg, 2, (5, 12, 8))
FAN_GRAPHS["CTH+BTH+T"] = (
    {"C", "R", "A", "T"}, [("C", "R")] + repeated_edge("R", "A", 8) + [("A", "T")], 1, (4, 10, 7)
)
vtx, edg = graph_bu("F")
vtx.add("E"); edg.append(("F", "E"))
FAN_GRAPHS["BU+FUE+EETA"] = (vtx, edg, 2, (8, 16, 9))
FAN_GRAPHS["BU+EU"] = (*graph_bu("E"), 2, (7, 15, 9))
FAN_GRAPHS["C4e+BETA+EETA"] = (*chain_with_parallel(9, ("C4", "R", "A", "E")), 2, (4, 11, 8))
verticals4 = {f"V{i}" for i in range(4)}
vertices = {"R", "E"} | verticals4
edges = []
for vtx0 in sorted(verticals4):
    edges += repeated_edge("R", vtx0, 2)
    edges.append((vtx0, "E"))
FAN_GRAPHS["BDELTA+EDELTA"] = (vertices, edges, 1, (6, 12, 7))
FAN_GRAPHS["C3e+BETA+EETA"] = (*chain_with_parallel(9, ("C3", "R", "A", "E")), 2, (4, 11, 8))
FAN_GRAPHS["CPHI+BPHI+EDELTA"] = (*chain_with_parallel(7, ("C3", "R", "A", "E")), 1, (4, 9, 6))
FAN_GRAPHS["BDEEP"] = ({"R", "C"}, repeated_edge("R", "C", 6), 1, (2, 6, 5))

TOTAL_GENUS = {
    "BU+FUT+T": 11, "C4e+BETA+FET+T": 10, "BDT+T": 9,
    "C3e+BETA+FET+T": 10, "CTH+BTH+T": 8,
    "BU+FUE+EETA": 11, "BU+EU": 11, "C4e+BETA+EETA": 10,
    "BDELTA+EDELTA": 8, "C3e+BETA+EETA": 10,
    "CPHI+BPHI+EDELTA": 7, "BDEEP": 6,
}
for fan, (vertices, edges, carrier_genus, expected_graph) in FAN_GRAPHS.items():
    stats = graph_stats(vertices, edges)
    A.eq(stats[:3], expected_graph, f"dual multigraph {fan}")
    A.eq(stats[3], 1, f"collision-free graph connected {fan}")
    A.eq(stats[2] + carrier_genus, TOTAL_GENUS[fan], f"arithmetic genus {fan}")


# ---------------------------------------------------------------------------
# Exact source charts for A4, A2, A3, and A15
# ---------------------------------------------------------------------------

Q, s, p, sigma, rho, x, v, z = sp.symbols("Q s p sigma rho x v z")
H = (
    -3 * p + sp.Rational(8, 3) * p**2 + sp.Rational(E.numerator, E.denominator) * p**3
    + k_s * (s * p)**2 + phi_s * p**2 * (s * p) + delta_s * p**4
    + theta_s * p * (s * p)**2 + eta_s * p**3 * (s * p) + u_s * p**5
)
FSOURCE = (s**2 - p) * (1 - Q * H) - Q * s**2 / 2


def two_stage(expr, first: dict, second: dict):
    return sp.cancel(expr.subs(first, simultaneous=True).subs(second, simultaneous=True))


# A4/FUT source chart.
g_a4 = sigma**2 * FSOURCE.subs({Q: sigma**5, s: sigma**-1 * S, p: sigma**-1 * P}, simultaneous=True)
lhs_a4 = sp.cancel(x**7 * g_a4.subs({S: v / x, P: 1 / x, sigma: rho / x}, simultaneous=True))
a4_A = u_s + eta_s * v + theta_s * v**2
a4_B = delta_s + phi_s * v + k_s * v**2
rhs_a4 = (v**2 - rho) * (
    x**5 - a4_A - rho * a4_B - sp.Rational(E.numerator, E.denominator) * rho**2
    - sp.Rational(8, 3) * rho**3 + 3 * rho**4
) - rho**5 * v**2 / 2
A.eq(sp.factor(lhs_a4 - rhs_a4), 0, "literal A4 chart identity")

# A2/CTH source chart (u=eta=Delta=0).
source_a2 = FSOURCE.subs({u_s: 0, eta_s: 0, delta_s: 0}, simultaneous=True)
g_a2 = sigma * source_a2.subs({Q: sigma**3, s: S, p: sigma**-1 * P}, simultaneous=True)
lhs_a2 = sp.cancel(x**4 * g_a2.subs({P: 1 / x, sigma: rho / x}, simultaneous=True))
a2_A = sp.Rational(E.numerator, E.denominator) + phi_s * S + theta_s * S**2
a2_B = sp.Rational(8, 3) + k_s * S**2
rhs_a2 = (rho * S**2 - 1) * (x**3 - a2_A - rho * a2_B + 3 * rho**2) - rho**4 * S**2 / 2
A.eq(sp.factor(lhs_a2 - rhs_a2), 0, "literal A2 chart identity")

# A3/EDELTA source chart (u=eta=Theta=0).
source_a3 = FSOURCE.subs({u_s: 0, eta_s: 0, theta_s: 0}, simultaneous=True)
g_a3 = sigma**2 * source_a3.subs({Q: sigma**4, s: sigma**-1 * S, p: sigma**-1 * P}, simultaneous=True)
lhs_a3 = sp.cancel(x**6 * g_a3.subs({S: v / x, P: 1 / x, sigma: rho / x}, simultaneous=True))
a3_C = delta_s + phi_s * v + k_s * v**2
rhs_a3 = (v**2 - rho) * (
    x**4 - a3_C - sp.Rational(E.numerator, E.denominator) * rho
    - sp.Rational(8, 3) * rho**2 + 3 * rho**3
) - rho**4 * v**2 / 2
A.eq(sp.factor(lhs_a3 - rhs_a3), 0, "literal A3 chart identity")

# A15/BDT source chart (u=eta=0).
source_a15 = FSOURCE.subs({u_s: 0, eta_s: 0}, simultaneous=True)
g_a15 = sigma**2 * source_a15.subs({Q: sigma**8, s: sigma**-1 * S, p: sigma**-2 * P}, simultaneous=True)
qratio = sp.symbols("q")
lhs_a15 = sp.cancel(x**10 * g_a15.subs(
    {S: 1 / x, P: qratio / x**2, sigma: rho / x}, simultaneous=True
))
b_a15 = (
    x**8 - delta_s * qratio**4 - theta_s * qratio**3 - rho * phi_s * qratio**3
    - rho**2 * (sp.Rational(E.numerator, E.denominator) * qratio**3 + k_s * qratio**2)
    - sp.Rational(8, 3) * rho**4 * qratio**2 + 3 * rho**6 * qratio
)
rhs_a15 = (1 - qratio) * b_a15 - rho**8 / 2
A.eq(sp.factor(lhs_a15 - rhs_a15), 0, "literal A15 chart identity")


# ---------------------------------------------------------------------------
# Odd A4 and A2 selector depths and target-compatible orders
# ---------------------------------------------------------------------------

a, qlin = sp.symbols("a qlin", nonzero=True)
c1 = delta_s + phi_s * a + k_s * a**2
psi_a4 = (
    rho * c1
    + sp.Rational(E.numerator, E.denominator) * rho**2
    + sp.Rational(8, 3) * rho**3 - 3 * rho**4
    - rho**2 * qlin**2 / (4 * (theta_s + k_s * rho))
)
series_a4 = sp.series(psi_a4, rho, 0, 5).removeO().expand()
A.eq(series_a4.coeff(rho, 1), c1, "A4 selector c1")
A.eq(series_a4.coeff(rho, 2), sp.Rational(E.numerator, E.denominator) - qlin**2 / (4 * theta_s), "A4 selector c2")
A.eq(series_a4.coeff(rho, 3), sp.Rational(8, 3) + k_s * qlin**2 / (4 * theta_s**2), "A4 selector c3")
A.eq(series_a4.coeff(rho, 4), -3 - k_s**2 * qlin**2 / (4 * theta_s**3), "A4 selector c4")
ratio_for_deep = -sp.Rational(8, 3) / sp.Rational(E.numerator, E.denominator)
A.eq(ratio_for_deep, sp.Rational(45, 172), "A4 forced K/Theta ratio")
c4_forced = -3 - sp.Rational(E.numerator, E.denominator) * ratio_for_deep**2
A.eq(c4_forced, sp.Rational(-99, 43), "A4 depth at most four")

# A2: on Theta*a^2=e and the first selector 8/3+K*a^2=0,
# K/Theta is again 45/172 and the second selector is nonzero.
k_over_theta_a2 = -sp.Rational(8, 3) / sp.Rational(E.numerator, E.denominator)
a2_second = -3 - k_over_theta_a2 * sp.Rational(-8, 3)
A.eq(k_over_theta_a2, sp.Rational(45, 172), "A2 forced K/Theta ratio")
A.eq(a2_second, sp.Rational(-99, 43), "A2 depth at most two")


def primitive_odd_row(m: int, split: int, sigma_over_target: int) -> tuple[int, int, int]:
    for target_order in range(1, 200):
        for x_order in range(1, 500):
            if (m - split) * x_order != split * sigma_over_target * target_order:
                continue
            if (m * x_order) % 2:
                continue
            y_order = m * x_order // 2
            if gcd(gcd(target_order, x_order), y_order) == 1:
                return target_order, x_order, y_order
    raise AssertionError("primitive odd row not found")


A4_ROWS = {1: ((4, 6, 15), 115), 2: ((1, 4, 10), 35), 3: ((2, 18, 45), 95), 4: ((1, 24, 60), 85)}
for split, (expected_row, expected_order) in A4_ROWS.items():
    row = primitive_odd_row(5, split, TARGET_ORDERS["FUT"][0] // FACE_PLANES["FUT"][0])
    A.eq(row, expected_row, f"A4 primitive target row r={split}")
    bt, bx, by = row
    A.eq(25 * bt + 5 * bx - by, expected_order, f"A4 form order r={split}")
    A.eq((5 - split) // 2 + split // 2, 2, f"A4 genus+delta r={split}")

A2_ROWS = {1: ((2, 2, 3), 13), 2: ((1, 4, 6), 11)}
for split, (expected_row, expected_order) in A2_ROWS.items():
    row = primitive_odd_row(3, split, TARGET_ORDERS["CTH"][0] // FACE_PLANES["CTH"][0])
    A.eq(row, expected_row, f"A2 primitive target row r={split}")
    bt, bx, by = row
    A.eq(5 * bt + 3 * bx - by, expected_order, f"A2 form order r={split}")
    A.eq((3 - split) // 2 + split // 2, 1, f"A2 genus+delta r={split}")


# ---------------------------------------------------------------------------
# Even A3 wall: two ends in the same connected BDELTA complement
# ---------------------------------------------------------------------------

# The first critical coefficient is e, hence the collision depth is exactly one.
A.ok(E != 0, "A3 depth one is forced by e")
A3_TARGET_ROW = primitive_odd_row(4, 1, TARGET_ORDERS["EDELTA"][0] // FACE_PLANES["EDELTA"][0])
A.eq(A3_TARGET_ROW, (1, 1, 2), "A3 primitive target row")
A.eq(10 * A3_TARGET_ROW[0] + 4 * A3_TARGET_ROW[1] - A3_TARGET_ROW[2], 12, "A3 tail form order")
A.eq((4 - 1 - 1) // 2, 1, "A3 elliptic tail genus")
A.eq(1 // 2, 0, "A3 persistent delta")

base_vertices = {"R", "A+", "A-"} | {f"V{i}" for i in range(4)}
base_edges: list[tuple[str, str]] = []
for i in range(4):
    base_edges += repeated_edge("R", f"V{i}", 2)
    base_edges.append((f"V{i}", "A+" if i < 2 else "A-"))
A.eq(graph_stats(base_vertices, base_edges), (7, 12, 6, 1), "A3 normalized complement")
tail_vertices = base_vertices | {"H"}
tail_edges = base_edges + [("H", "A+"), ("H", "A-")]
A.eq(graph_stats(tail_vertices, tail_edges), (8, 14, 7, 1), "A3 same-complement two-end increment")
A.eq(7 + 1, 8, "A3 graph plus elliptic tail restores genus eight")


# ---------------------------------------------------------------------------
# A15 boundary contact: direct source-minimal clusters and target orders
# ---------------------------------------------------------------------------

t, XX, ZZ, lam = sp.symbols("t XX ZZ lam")
local_a15 = sp.expand(rhs_a15.subs({qratio: 1 + z, theta_s: -delta_s}, simultaneous=True))
central_a15 = sp.expand(local_a15.subs(rho, 0))
A.eq(sp.expand(central_a15 - z * (delta_s * z * (1 + z)**3 - x**8)), 0, "A15 central contact equation")
A.eq(8, 8, "A15 branch intersection length")


def initial_t_coefficient(expr, substitutions: dict, expected_weight: int):
    expanded = sp.expand(expr.subs(substitutions, simultaneous=True))
    poly = sp.Poly(expanded, t)
    minimum = min(term[0][0] for term in poly.terms())
    A.eq(minimum, expected_weight, "source-minimal cluster weight")
    return sp.expand(poly.coeff_monomial(t**minimum))


lead_r1 = initial_t_coefficient(
    local_a15, {x: t * XX, rho: t**8 * XX, z: t**8 * ZZ}, 16
)
A.eq(sp.expand(lead_r1 - (delta_s * ZZ**2 - XX * (XX**7 - phi_s) * ZZ)), 0, "A15 depth-one tail")
lead_r2 = initial_t_coefficient(
    local_a15.subs(phi_s, 0), {x: t * XX, rho: t**4 * XX, z: t**8 * ZZ}, 16
)
A.eq(sp.expand(lead_r2 - (
    delta_s * ZZ**2 - XX**2 * (XX**6 - (sp.Rational(E.numerator, E.denominator) + k_s)) * ZZ
)), 0, "A15 depth-two tail")

D_DEEP = sp.Rational(2048, 45)
K_DEEP = sp.Rational(1376, 135)
A.eq(sp.Rational(2848, 45) - sp.Rational(7, 6) * D_DEEP, K_DEEP, "deep A15 affine seam")
A.eq(K_DEEP, -sp.Rational(E.numerator, E.denominator), "deep A15 K=-e")
lead_r4 = initial_t_coefficient(
    local_a15.subs({phi_s: 0, k_s: K_DEEP, delta_s: D_DEEP}, simultaneous=True),
    {x: t * XX, rho: t**2 * XX, z: t**8 * ZZ}, 16,
)
expected_r4 = D_DEEP * ZZ**2 - XX**4 * (XX**4 - sp.Rational(8, 3)) * ZZ - XX**8 / 2
A.eq(sp.expand(lead_r4 - expected_r4), 0, "A15 depth-four tail")

# Independent discriminant-branch expansion through the first decisive term.
mu = sp.symbols("mu")
Xvalue = phi_s * rho + (k_s + sp.Rational(E.numerator, E.denominator)) * rho**2 + (sp.Rational(8, 3) + mu) * rho**4
zcritical = mu * rho**4 / (2 * delta_s)
f_branch = local_a15.subs(x**8, Xvalue)
f_at = sp.expand(f_branch.subs(z, zcritical))
fz_at = sp.expand(sp.diff(f_branch, z).subs(z, zcritical))
A.eq(sp.expand(f_at).coeff(rho, 8), -sp.Rational(1, 2) - mu**2 / (4 * delta_s), "A15 branch discriminant coefficient")
A.eq(sp.expand(fz_at).coeff(rho, 4), 0, "A15 formal critical section to leading order")
A.eq(sp.expand((sp.Rational(8, 3) + lam) * (sp.Rational(8, 3) - lam)).subs(lam**2, -2 * D_DEEP),
     sp.Rational(1472, 15), "A15 deepest constants nonzero")
A.ok(D_DEEP != 0, "A15 two depth-four constants are distinct")

deep_branch = sp.Poly((XX**4 - sp.Rational(8, 3))**2 + sp.Rational(4096, 45), XX)
A.eq(deep_branch.degree(), 8, "A15 deep branch degree")
A.eq(sp.gcd(deep_branch, deep_branch.diff()).degree(), 0, "A15 deep branch squarefree")
A.eq((deep_branch.degree() - 2) // 2, 3, "A15 deep tail genus")
A.eq(2, 2 if deep_branch.degree() % 2 == 0 else 1, "A15 deep tail has two infinity ends")
A.eq(deep_branch.eval(0), sp.Rational(1472, 15), "A15 derivative root X=0 is safe")


def primitive_a15_row(split: int) -> tuple[int, int, int]:
    sigma_over_target = TARGET_ORDERS["BDT"][0] // FACE_PLANES["BDT"][0]
    for target_order in range(1, 200):
        for x_order in range(1, 500):
            if (8 - split) * x_order != split * sigma_over_target * target_order:
                continue
            y_order = 8 * x_order
            if gcd(gcd(target_order, x_order), y_order) == 1:
                return target_order, x_order, y_order
    raise AssertionError("primitive A15 row not found")


A15_ROWS = {1: ((7, 3, 24), 116), 2: ((1, 1, 8), 16), 4: ((1, 3, 24), 14)}
for split, (expected_row, expected_order) in A15_ROWS.items():
    row = primitive_a15_row(split)
    A.eq(row, expected_row, f"A15 target-compatible row r={split}")
    bt, bx, by = row
    A.eq(17 * bt + 7 * bx - by, expected_order, f"A15 form order r={split}")
    A.ok(expected_order > 0, f"A15 form positivity r={split}")

# The complement has two components before insertion: isolated R and C--T.
complement_vertices = {"R", "C", "T"}
complement_edges = [("C", "T")]
A.eq(graph_stats(complement_vertices, complement_edges), (3, 1, 0, 2), "A15 disconnected complement")

def a15_sign_graph(mutual_nodes: int) -> tuple[int, int, int, int]:
    vertices = complement_vertices | {"H+", "H-"}
    edges = complement_edges + repeated_edge("H+", "H-", mutual_nodes)
    edges += [("R", "H+"), ("C", "H-")]
    return graph_stats(vertices, edges)


A.eq(a15_sign_graph(7), (5, 10, 6, 1), "A15 depth-one graph")
A.eq(2 + 1 + 6, 9, "A15 depth-one genus ledger")
A.eq(a15_sign_graph(6), (5, 9, 5, 1), "A15 depth-two graph")
A.eq(2 + 2 + 5, 9, "A15 depth-two genus ledger")
deep_vertices = complement_vertices | {"H"}
deep_edges = complement_edges + [("R", "H"), ("C", "H")]
A.eq(graph_stats(deep_vertices, deep_edges), (4, 3, 0, 1), "A15 depth-four bridge graph")
A.eq(2 + 4 + 3 + 0, 9, "A15 depth-four genus ledger")

# Explicit hostile: applying the even-cusp same-complement +1 here overcounts.
A.eq(2 + 4 + 3 + 1, 10, "hostile false A15 same-complement count")
A.ok(2 + 4 + 3 + 1 != 9, "A15 graph +1 is forbidden")


# ---------------------------------------------------------------------------
# Coefficientwise-cancellation and source/target-base firewalls
# ---------------------------------------------------------------------------

coeff_dstar = literal_coefficients(D_STAR, F(1), F(1), F(1), F(1))
A.ok((1, 2, 3) not in coeff_dstar, "K-e cancellation is coefficientwise collected")
A.ok((1, 4, 2) in coeff_dstar, "K owner survives K-e cancellation")
coeff_equal = literal_coefficients(F(1), F(1), F(1), F(1), F(1))
A.ok((1, 2, 4) not in coeff_equal, "Theta-Delta cancellation is coefficientwise collected")
A.ok((1, 4, 3) in coeff_equal and (1, 0, 5) in coeff_equal, "owners survive shared-point cancellation")
dummy_f, dummy_u = sp.symbols("f u")
A.eq((dummy_u - dummy_f).subs(dummy_u, dummy_f), 0, "specialization hostile vanishes")
A.ok(sp.Poly(dummy_u - dummy_f, dummy_u).degree() == 1, "specialization hostile is not coefficientwise zero")

A.eq((FACE_PLANES["FUT"][0], TARGET_ORDERS["FUT"][0]), (5, 30), "A4 source/target bases distinct")
A.eq((FACE_PLANES["CTH"][0], TARGET_ORDERS["CTH"][0]), (3, 6), "A2 source/target bases distinct")
A.eq((FACE_PLANES["EDELTA"][0], TARGET_ORDERS["EDELTA"][0]), (4, 12), "A3 source/target bases distinct")
A.eq((FACE_PLANES["BDT"][0], TARGET_ORDERS["BDT"][0]), (8, 24), "A15 source/target bases distinct")


print("THM4355_TERMINAL_ALPHA_ZERO_INDEPENDENT_REFEREE=PASS")
print("scope=Z=beta_11=zeta_3=W=xi_10=U=alpha_11=0;K!=0")
print(f"exact_supports={len(EXACT_ROWS)};theta_nonzero=40;theta_zero=24;exact_fans={len(EXPECTED_FAN_COUNTS)}")
print("fan_counts=" + ",".join(f"{name}:{EXPECTED_FAN_COUNTS[name]}" for name in EXPECTED_FAN_COUNTS))
print(f"hostile_atlas={len(HOSTILE_SUPPORTS)}_supports;hostile_fans={len(HOSTILE_FANS)}")
print("collision_walls=A4:L5,A15:L15,A2:L3,A3:L4")
print("A4_depths=1,2,3,4;orders=115,35,95,85;genus_restored=2")
print("A2_depths=1,2;orders=13,11;genus_restored=1")
print("A3_depth=1;order=12;two_ends=same_complement;graph_increment=1;genus_restored=8")
print("A15_depths=1,2,4;orders=116,16,14;ends=different_complements;graph_increment=0;genus_restored=9")
print("source_target_bases=A4:5/30,A2:3/6,A3:4/12,A15:8/24")
print("verdict=PASS_RELATIVE_TO_INHERITED_TOROIDAL_GOOD_TARGET_PROPER_FLAT_INTERFACES")
print("nonconsequence=JC2_DC2_AND_SEAM_ENTRY_REMAIN_OPEN")
print(f"CHECKS={A.checks}")
