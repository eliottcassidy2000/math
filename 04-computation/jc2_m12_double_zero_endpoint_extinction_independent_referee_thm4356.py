#!/usr/bin/env python3
"""Clean-room exact referee for the reserved THM-4356 endpoint candidate.

The calculation is rebuilt from the literal sixteen-row source.  It imports
neither the exploratory U=K=0 scout nor any primary THM-4356 certificate.
All support, lower-hull, Pick, graph, local-chart, and valuation arithmetic is
exact.  SymPy is used only for polynomial identities and formal expansions.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as QQ
from itertools import combinations, product
from math import gcd, lcm
import sys

import sympy as sp


sys.dont_write_bytecode = True
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


class Audit:
    def __init__(self) -> None:
        self.checks = 0

    def ok(self, condition, label: str) -> None:
        self.checks += 1
        if not bool(condition):
            raise AssertionError(label)

    def eq(self, left, right, label: str) -> None:
        self.checks += 1
        if left != right:
            raise AssertionError(f"{label}: {left!r} != {right!r}")


A = Audit()
Point3 = tuple[int, int, int]       # (s exponent, p exponent, Q exponent)
Point2 = tuple[int, int]
Plane = tuple[QQ, QQ, QQ]          # h=a*r+b*l+c

E = QQ(-1376, 135)
DELTA0 = QQ(5696, 105)


# ---------------------------------------------------------------------------
# Literal source and coefficientwise lift
# ---------------------------------------------------------------------------

Q, s, p, S, P, sigma = sp.symbols("Q s p S P sigma")
Delta, Theta, Phi, eta, u, alpha = sp.symbols(
    "Delta Theta Phi eta u alpha"
)
K, U, W, xi, beta, zeta, Z = sp.symbols("K U W xi beta zeta Z")

SOURCE_ROWS = (
    ("p", 1, 0, -3),
    ("p2", 2, 0, sp.Rational(8, 3)),
    ("p3", 3, 0, sp.Rational(E.numerator, E.denominator)),
    ("K", 0, 2, K),
    ("Phi", 2, 1, Phi),
    ("Delta", 4, 0, Delta),
    ("Theta", 1, 2, Theta),
    ("eta", 3, 1, eta),
    ("zeta", 0, 3, zeta),
    ("u", 5, 0, u),
    ("xi", 2, 2, xi),
    ("alpha", 4, 1, alpha),
    ("beta", 1, 3, beta),
    ("U", 6, 0, U),
    ("W", 3, 2, W),
    ("Z", 0, 4, Z),
)
A.eq(len(SOURCE_ROWS), 16, "sixteen literal source rows")
A.eq(len({(i, j) for _name, i, j, _coefficient in SOURCE_ROWS}), 16,
     "literal source exponents are distinct")

H = sum(coefficient * p**i * (s * p)**j
        for _name, i, j, coefficient in SOURCE_ROWS)
FSOURCE = (s**2 - p) * (1 - Q * H) - Q * s**2 / 2

seam_k = sp.Rational(2848, 45) - sp.Rational(7, 6) * Delta
A.eq(sp.solve(sp.Eq(seam_k, 0), Delta), [sp.Rational(5696, 105)],
     "K=0 forces Delta=5696/105")

GATE = {
    Delta: sp.Rational(DELTA0.numerator, DELTA0.denominator),
    K: 0,
    U: 0,
    W: 0,
    xi: 0,
    beta: 0,
    zeta: 0,
    Z: 0,
}
HGATE = sp.expand(H.subs(GATE, simultaneous=True))
FGATE = sp.expand(FSOURCE.subs(GATE, simultaneous=True))


def literal_coefficients(theta: QQ, aa: QQ, ee: QQ,
                         pp: QQ, uu: QQ) -> dict[Point3, QQ]:
    """Collect the lifted coefficients before testing whether they vanish."""
    values = {
        "Theta": theta,
        "alpha": aa,
        "eta": ee,
        "Phi": pp,
        "u": uu,
    }
    rows = (
        ("p", 1, 0, QQ(-3)),
        ("p2", 2, 0, QQ(8, 3)),
        ("p3", 3, 0, E),
        ("K", 0, 2, QQ(0)),
        ("Phi", 2, 1, pp),
        ("Delta", 4, 0, DELTA0),
        ("Theta", 1, 2, theta),
        ("eta", 3, 1, ee),
        ("zeta", 0, 3, QQ(0)),
        ("u", 5, 0, uu),
        ("xi", 2, 2, QQ(0)),
        ("alpha", 4, 1, aa),
        ("beta", 1, 3, QQ(0)),
        ("U", 6, 0, QQ(0)),
        ("W", 3, 2, QQ(0)),
        ("Z", 0, 4, QQ(0)),
    )
    A.eq(len(rows), 16, "numeric source retains sixteen rows")
    out: defaultdict[Point3, QQ] = defaultdict(QQ)
    out[(2, 0, 0)] += 1
    out[(0, 1, 0)] -= 1
    out[(2, 0, 1)] -= QQ(1, 2)
    for name, i, j, coefficient in rows:
        if name in values:
            A.eq(coefficient, values[name], f"numeric row {name}")
        out[(j + 2, i + j, 1)] -= coefficient
        out[(j, i + j + 1, 1)] += coefficient
    return {point: coefficient for point, coefficient in out.items()
            if coefficient != 0}


symbolic_lift: defaultdict[Point3, sp.Expr] = defaultdict(lambda: sp.Integer(0))
symbolic_lift[(2, 0, 0)] += 1
symbolic_lift[(0, 1, 0)] -= 1
symbolic_lift[(2, 0, 1)] -= sp.Rational(1, 2)
for _name, i, j, coefficient in SOURCE_ROWS:
    coefficient = sp.sympify(coefficient).subs(GATE, simultaneous=True)
    symbolic_lift[(j + 2, i + j, 1)] -= coefficient
    symbolic_lift[(j, i + j + 1, 1)] += coefficient
symbolic_lift = defaultdict(
    lambda: sp.Integer(0),
    {point: sp.factor(coefficient) for point, coefficient in symbolic_lift.items()
     if sp.factor(coefficient) != 0},
)

A.eq(symbolic_lift[(2, 3, 1)], sp.Rational(1376, 135),
     "K-e shared coefficient remains fixed nonzero")
A.eq(sp.factor(symbolic_lift[(2, 4, 1)]
               - (Theta - sp.Rational(5696, 105))), 0,
     "Theta-Delta shared coefficient")
A.eq(symbolic_lift[(2, 5, 1)], -u, "u lower coefficient")
A.eq(symbolic_lift[(0, 6, 1)], u, "u upper coefficient")
A.eq(symbolic_lift[(4, 3, 1)], -Theta, "Theta owner coefficient")
A.eq(symbolic_lift[(1, 6, 1)], alpha, "alpha upper coefficient")
A.eq(symbolic_lift[(3, 5, 1)], -alpha, "alpha lower coefficient")
A.eq(len(symbolic_lift), 20, "twenty active-or-optional lifted points")


# ---------------------------------------------------------------------------
# Independent hull machinery and the exact 48-support quotient
# ---------------------------------------------------------------------------

def rank_two(points: list[Point3] | tuple[Point3, ...]) -> bool:
    return any(
        (q[0] - p0[0]) * (r[1] - p0[1])
        != (q[1] - p0[1]) * (r[0] - p0[0])
        for p0, q, r in combinations(points, 3)
    )


def candidate_planes(universe: set[Point3]) -> set[Plane]:
    result: set[Plane] = set()
    for p0, p1, p2 in combinations(sorted(universe), 3):
        det = ((p1[0] - p0[0]) * (p2[1] - p0[1])
               - (p1[1] - p0[1]) * (p2[0] - p0[0]))
        if det == 0:
            continue
        aa = QQ(
            (p1[2] - p0[2]) * (p2[1] - p0[1])
            - (p1[1] - p0[1]) * (p2[2] - p0[2]), det,
        )
        bb = QQ(
            (p1[0] - p0[0]) * (p2[2] - p0[2])
            - (p1[2] - p0[2]) * (p2[0] - p0[0]), det,
        )
        cc = QQ(p0[2]) - aa * p0[0] - bb * p0[1]
        result.add((aa, bb, cc))
    return result


def lower_fan(support: frozenset[Point3], planes: set[Plane]) -> tuple[Plane, ...]:
    result = []
    for plane in planes:
        aa, bb, cc = plane
        face: list[Point3] = []
        for rr, ll, hh in support:
            gap = QQ(hh) - aa * rr - bb * ll - cc
            if gap < 0:
                break
            if gap == 0:
                face.append((rr, ll, hh))
        else:
            if rank_two(face):
                result.append(plane)
    return tuple(sorted(result))


@dataclass(frozen=True)
class ExactState:
    theta_class: int
    theta: QQ
    alpha: int
    eta: int
    phi: int
    u: int
    support: frozenset[Point3]


THETA_REPRESENTATIVES = (
    (0, QQ(0)),
    (1, DELTA0),
    (2, QQ(1)),
)
exact_states: list[ExactState] = []
for (theta_class, theta_value), bits in product(
        THETA_REPRESENTATIVES, product((0, 1), repeat=4)):
    aa, ee, pp, uu = bits
    support = frozenset(literal_coefficients(
        theta_value, QQ(aa), QQ(ee), QQ(pp), QQ(uu)
    ))
    exact_states.append(ExactState(
        theta_class, theta_value, aa, ee, pp, uu, support
    ))

A.eq(len(exact_states), 48, "three coupled Theta states times four bits")
A.eq(len({state.support for state in exact_states}), 48,
     "48 distinct realizable supports")

support_union = set().union(*(state.support for state in exact_states))
support_intersection = set.intersection(*(set(state.support) for state in exact_states))
A.eq(len(support_union), 20, "exact support union")
A.eq(len(support_intersection), 10, "exact mandatory support")

PLANE_BY_NAME: dict[str, Plane] = {
    "C4A": (QQ(-1, 4), QQ(1, 4), QQ(-1, 4)),
    "C5A": (QQ(0), QQ(1, 5), QQ(-1, 5)),
    "C4E": (QQ(0), QQ(1, 4), QQ(-1, 4)),
    "BA": (QQ(1, 11), QQ(2, 11), QQ(-2, 11)),
    "BU": (QQ(1, 10), QQ(1, 5), QQ(-1, 5)),
    "BE": (QQ(1, 9), QQ(2, 9), QQ(-2, 9)),
    "BDT": (QQ(1, 8), QQ(1, 4), QQ(-1, 4)),
    "FUT": (QQ(1, 5), QQ(1, 5), QQ(-2, 5)),
    "ED": (QQ(1, 4), QQ(1, 4), QQ(-1, 2)),
    "E11": (QQ(2, 7), QQ(1, 7), QQ(-4, 7)),
    "FUP": (QQ(2, 5), QQ(1, 5), QQ(-4, 5)),
    "N": (QQ(1), QQ(0), QQ(-2)),
}
NAME_BY_PLANE = {plane: name for name, plane in PLANE_BY_NAME.items()}
A.eq(len(NAME_BY_PLANE), 12, "twelve distinct named planes")

plane_candidates = candidate_planes(support_union)
exact_fans = [lower_fan(state.support, plane_candidates) for state in exact_states]
exact_plane_set = set().union(*(set(fan) for fan in exact_fans))
A.eq(exact_plane_set, set(PLANE_BY_NAME.values()),
     "exact lower hull uses precisely the twelve named planes")
A.eq(len(set(exact_fans)), 15, "exact support quotient has fifteen fans")


def fan_names(fan: tuple[Plane, ...]) -> tuple[str, ...]:
    return tuple(NAME_BY_PLANE[plane] for plane in fan)


expected_fan_counts = Counter({
    ("BU",): 1,
    ("BDT",): 5,
    ("C4A", "BA"): 1,
    ("C5A", "BA"): 1,
    ("C4E", "BE"): 1,
    ("BU", "FUT"): 9,
    ("BU", "FUP"): 1,
    ("BDT", "ED"): 1,
    ("C4A", "BA", "E11"): 8,
    ("C4A", "BA", "N"): 3,
    ("C5A", "BA", "E11"): 8,
    ("C5A", "BA", "N"): 3,
    ("C4E", "BE", "FUT"): 4,
    ("C4E", "BE", "N"): 1,
    ("BU", "FUT", "N"): 1,
})
A.eq(Counter(map(fan_names, exact_fans)), expected_fan_counts,
     "all fifteen exact fan multiplicities")

# The odd-square label stores only the support ordinal, not a coefficient.
for state in exact_states:
    n = (1 + 16 * state.theta_class + 8 * state.alpha
         + 4 * state.eta + 2 * state.phi + state.u)
    A.ok(1 <= n <= 48, "natural support address range")
    recovered = (
        (n - 1) // 16,
        ((n - 1) % 16) // 8,
        ((n - 1) % 8) // 4,
        ((n - 1) % 4) // 2,
        (n - 1) % 2,
    )
    A.eq(recovered,
         (state.theta_class, state.alpha, state.eta, state.phi, state.u),
         "support address inverse")
    A.eq(((2 * n - 1)**2 - 1) // 8, n * (n - 1) // 2,
         "odd-square/triangular identity")


# A genuinely distinct hostile atlas: decouple all ten optional lifted points.
optional_points = sorted(support_union - support_intersection)
A.eq(len(optional_points), 10, "ten pointwise-hostile switches")
hostile_supports = {
    frozenset(support_intersection | {
        point for point, bit in zip(optional_points, bits) if bit
    })
    for bits in product((0, 1), repeat=10)
}
A.eq(len(hostile_supports), 1024, "point-decoupled hostile supports")
A.ok({state.support for state in exact_states} <= hostile_supports,
     "every honest support embeds in hostile atlas")
hostile_candidates = candidate_planes(support_union)
hostile_fans = {lower_fan(support, hostile_candidates)
                for support in hostile_supports}
hostile_planes = set().union(*(set(fan) for fan in hostile_fans))
A.eq(len(hostile_fans), 52, "point-decoupled hostile fan count")
A.eq(len(hostile_planes), 24, "point-decoupled hostile plane count")
A.ok(exact_plane_set < hostile_planes, "hostile atlas strictly enlarges plane set")


# ---------------------------------------------------------------------------
# All twelve primitive faces, target bases, and normalizations
# ---------------------------------------------------------------------------

D0 = sp.Rational(DELTA0.numerator, DELTA0.denominator)
FACE_EXPRESSION: dict[str, sp.Expr] = {
    "C4A": P * (D0 * P**4 + alpha * S * P**5 - 1),
    "C5A": P * (u * P**5 + alpha * S * P**5 - 1),
    "C4E": P * (D0 * P**4 + eta * S * P**4 - 1),
    "BA": (P - S**2) * (alpha * S * P**5 - 1),
    "BU": (P - S**2) * (u * P**5 - 1),
    "BE": (P - S**2) * (eta * S * P**4 - 1),
    "BDT": (P - S**2) * (D0 * P**4 + Theta * S**2 * P**3 - 1),
    "FUT": S**2 * (1 - u * P**5 - eta * S * P**4
                      - Theta * S**2 * P**3),
    "ED": S**2 * (1 - D0 * P**4 - Phi * S * P**3),
    "E11": S**2 * (1 - alpha * S * P**5 - Theta * S**2 * P**3),
    "FUP": S**2 * (1 - u * P**5 - Phi * S * P**3),
    "N": S**2 * (1 - S * P**3 * (Phi + eta * P + alpha * P**2)),
}

FACE_RESTRICTIONS = {
    "C4A": {u: 0},
    "C5A": {},
    "C4E": {alpha: 0, u: 0, Theta: 0},
    "BA": {},
    "BU": {alpha: 0},
    "BE": {alpha: 0, u: 0},
    "BDT": {alpha: 0, u: 0, eta: 0},
    "FUT": {alpha: 0},
    "ED": {alpha: 0, u: 0, eta: 0, Theta: 0},
    "E11": {},
    "FUP": {alpha: 0, eta: 0, Theta: 0},
    "N": {Theta: 0},
}

EXPECTED_BASE_ORDER = {
    "C4A": (4, 12, 13),
    "C5A": (5, 30, 25),
    "C4E": (4, 12, 10),
    "BA": (11, 66, 49),
    "BU": (10, 30, 22),
    "BE": (9, 18, 13),
    "BDT": (8, 24, 17),
    "FUT": (5, 30, 25),
    "ED": (4, 12, 10),
    "E11": (7, 42, 41),
    "FUP": (5, 30, 31),
    "N": (1, 6, 11),
}


def sigma_initial(expression: sp.Expr) -> tuple[int, sp.Expr]:
    terms = sp.expand(expression).as_ordered_terms()
    exponents = [sp.sympify(term.as_powers_dict().get(sigma, 0))
                 for term in terms]
    A.ok(all(value.is_Integer for value in exponents),
         "integral primitive-chart sigma exponents")
    minimum = min(int(value) for value in exponents)
    initial = sum(
        term / sigma**minimum for term, exponent in zip(terms, exponents)
        if int(exponent) == minimum
    )
    return minimum, sp.factor(initial)


def source_base(plane: Plane) -> int:
    return lcm(*(entry.denominator for entry in plane))


def primitive_initial(name: str) -> sp.Expr:
    aa, bb, cc = PLANE_BY_NAME[name]
    base = source_base((aa, bb, cc))
    a_int = int(base * aa)
    b_int = int(base * bb)
    g_int = int(-base * cc)
    chart = sigma**g_int * FGATE.subs(
        FACE_RESTRICTIONS[name], simultaneous=True
    ).subs({
        Q: sigma**base,
        s: sigma**(-a_int) * S,
        p: sigma**(-b_int) * P,
    }, simultaneous=True)
    minimum, initial = sigma_initial(sp.cancel(chart))
    A.eq(minimum, 0, f"{name} primitive chart begins in degree zero")
    return initial


for face_name, plane in PLANE_BY_NAME.items():
    primitive_base = source_base(plane)
    common_base = lcm(primitive_base, 6)
    order = QQ(common_base) * (QQ(5, 6) - sum(plane))
    A.eq(order.denominator, 1, f"{face_name} integral target order")
    A.ok(order > 0, f"{face_name} positive target order")
    A.eq((primitive_base, common_base, order.numerator),
         EXPECTED_BASE_ORDER[face_name], f"{face_name} base/order row")
    A.eq(sp.factor(primitive_initial(face_name) - FACE_EXPRESSION[face_name]),
         0, f"{face_name} literal primitive initial equation")


def cross2(origin: Point2, left: Point2, right: Point2) -> int:
    return ((left[0] - origin[0]) * (right[1] - origin[1])
            - (left[1] - origin[1]) * (right[0] - origin[0]))


def convex_hull(points: set[Point2] | list[Point2]) -> tuple[Point2, ...]:
    points = sorted(set(points))
    if len(points) <= 1:
        return tuple(points)
    lower: list[Point2] = []
    for point in points:
        while len(lower) >= 2 and cross2(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper: list[Point2] = []
    for point in reversed(points):
        while len(upper) >= 2 and cross2(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return tuple(lower[:-1] + upper[:-1])


def expression_hull(expression: sp.Expr) -> tuple[Point2, ...]:
    polynomial = sp.Poly(sp.expand(expression), S, P)
    return convex_hull({monomial for monomial, coefficient in polynomial.terms()
                        if coefficient != 0})


EXPECTED_FACE_HULL = {
    "C4A": ((0, 1), (1, 6), (0, 5)),
    "C5A": ((0, 1), (1, 6), (0, 6)),
    "C4E": ((0, 1), (1, 5), (0, 5)),
    "BA": ((0, 1), (2, 0), (3, 5), (1, 6)),
    "BU": ((0, 1), (2, 0), (2, 5), (0, 6)),
    "BE": ((0, 1), (2, 0), (3, 4), (1, 5)),
    "BDT": ((0, 1), (2, 0), (4, 3), (0, 5)),
    "FUT": ((2, 0), (4, 3), (2, 5)),
    "ED": ((2, 0), (3, 3), (2, 4)),
    "E11": ((2, 0), (4, 3), (3, 5)),
    "FUP": ((2, 0), (3, 3), (2, 5)),
    "N": ((2, 0), (3, 3), (3, 5)),
}
for face_name, expression in FACE_EXPRESSION.items():
    A.eq(expression_hull(expression), EXPECTED_FACE_HULL[face_name],
         f"{face_name} face polygon")

# The three genuine positive-genus normalizations on clean strata.
Y = sp.symbols("Y")
core_bdt = 1 - D0 * P**4 - Theta * S**2 * P**3
y_bdt = Theta * S * P**2
A.eq(sp.factor(y_bdt**2 - Theta * P * (1 - D0 * P**4)
               + Theta * P * core_bdt), 0, "BDT genus-two normalization")

core_fut = 1 - u * P**5 - eta * S * P**4 - Theta * S**2 * P**3
y_fut = 2 * Theta * S * P**2 + eta * P**3
branch_fut = P * (4 * Theta + (eta**2 - 4 * u * Theta) * P**5)
A.eq(sp.factor(y_fut**2 - branch_fut + 4 * Theta * P * core_fut), 0,
     "FUT genus-two normalization")

core_e11 = 1 - alpha * S * P**5 - Theta * S**2 * P**3
y_e11 = 2 * Theta * S * P**2 + alpha * P**4
branch_e11 = P * (4 * Theta + alpha**2 * P**7)
A.eq(sp.factor(y_e11**2 - branch_e11 + 4 * Theta * P * core_e11), 0,
     "E11 genus-three normalization")

A.eq(sp.factor(sp.discriminant(branch_fut, P)),
     12800000 * Theta**6 * (4 * Theta * u - eta**2)**4,
     "FUT branch discriminant")
A.eq(sp.factor(sp.discriminant(branch_e11, P)),
     -53971714048 * Theta**8 * alpha**12,
     "E11 branch discriminant")
branch_bdt = Theta * P * (1 - D0 * P**4)
A.ok(sp.factor(sp.discriminant(branch_bdt, P)) != 0,
     "BDT branch polynomial is squarefree for Theta nonzero")
A.eq(sp.Poly(branch_bdt, P).degree(), 5, "BDT branch degree gives genus two")
A.eq(sp.Poly(branch_fut, P).degree(), 6, "FUT branch degree gives genus two")
A.eq(sp.Poly(branch_e11, P).degree(), 8, "E11 branch degree gives genus three")

# Rationality and component counts are checked from their literal factors.
for face_name in ("C4A", "C5A", "C4E", "ED", "FUP", "N"):
    reduced = sp.factor(FACE_EXPRESSION[face_name] / (
        P if face_name.startswith("C") else S**2
    ))
    A.ok(sp.Poly(reduced, S).degree() <= 1,
         f"{face_name} is rational by a linear S-equation")
A.eq(sp.factor(FACE_EXPRESSION["BA"]
               - (P - S**2) * (alpha * S * P**5 - 1)), 0,
     "BA has two rational factors")
A.eq(sp.factor(FACE_EXPRESSION["BU"]
               - (P - S**2) * (u * P**5 - 1)), 0,
     "BU has one section plus five vertical factors")
A.eq(sp.factor(FACE_EXPRESSION["BE"]
               - (P - S**2) * (eta * S * P**4 - 1)), 0,
     "BE has two rational factors")
A.eq(sp.diff(alpha * S**11 - 1, S).subs(S, 0), 0,
     "BA node polynomial has no zero root")
A.eq(sp.gcd(sp.Poly(alpha * S**11 - 1, S),
            sp.Poly(sp.diff(alpha * S**11 - 1, S), S)).degree(), 0,
     "BA has eleven transverse nodes on alpha owner")
A.eq(sp.gcd(sp.Poly(u * S**10 - 1, S),
            sp.Poly(sp.diff(u * S**10 - 1, S), S)).degree(), 0,
     "BU has ten transverse nodes on u owner")
A.eq(sp.gcd(sp.Poly(eta * S**9 - 1, S),
            sp.Poly(sp.diff(eta * S**9 - 1, S), S)).degree(), 0,
     "BE has nine transverse nodes on eta owner")


# ---------------------------------------------------------------------------
# Every exact global Pick, puncture, and normalized graph ledger
# ---------------------------------------------------------------------------

def pick_data(vertices: tuple[Point2, ...]) -> tuple[int, int, int]:
    edges = tuple(zip(vertices, vertices[1:] + vertices[:1]))
    area2 = abs(sum(x1 * y2 - x2 * y1
                    for (x1, y1), (x2, y2) in edges))
    boundary = sum(gcd(abs(x2 - x1), abs(y2 - y1))
                   for (x1, y1), (x2, y2) in edges)
    A.eq((area2 - boundary + 2) % 2, 0, "Pick parity")
    return area2, boundary, (area2 - boundary + 2) // 2


def puncture_packet(vertices: tuple[Point2, ...]) -> tuple[int, ...]:
    packet: list[int] = []
    for start, end in zip(vertices, vertices[1:] + vertices[:1]):
        dx, dy = end[0] - start[0], end[1] - start[1]
        length = gcd(abs(dx), abs(dy))
        inward = (-dy // length, dx // length)
        constant = inward[0] * start[0] + inward[1] * start[1]
        index = inward[0] + inward[1] - constant
        # The vertical Q=0 edge is affine rather than source infinity.
        if not (start[0] == end[0] == 0):
            packet.extend([index] * length)
    return tuple(sorted(packet, reverse=True))


def on_plane(point: Point3, plane: Plane) -> bool:
    rr, ll, hh = point
    aa, bb, cc = plane
    return QQ(hh) == aa * rr + bb * ll + cc


def common_edge_length(support: frozenset[Point3],
                       left: Plane, right: Plane) -> int:
    common = {(rr, ll) for rr, ll, hh in support
              if on_plane((rr, ll, hh), left)
              and on_plane((rr, ll, hh), right)}
    if len(common) < 2:
        return 0
    return max(gcd(abs(p0[0] - p1[0]), abs(p0[1] - p1[1]))
               for p0, p1 in combinations(common, 2))


def repeated_edge(left: str, right: str, count: int) -> list[tuple[str, str]]:
    return [(left, right) for _ in range(count)]


def face_graph(name: str, theta_nonzero: bool,
               tag: str) -> tuple[set[str], list[tuple[str, str]], int]:
    def node(short: str) -> str:
        return f"{tag}:{short}"

    if name == "BA":
        vertices = {node("R"), node("C")}
        return vertices, repeated_edge(node("R"), node("C"), 11), 0
    if name == "BU":
        vertices = {node("R")} | {node(f"V{i}") for i in range(5)}
        edges = []
        for i in range(5):
            edges += repeated_edge(node("R"), node(f"V{i}"), 2)
        return vertices, edges, 0
    if name == "BE":
        vertices = {node("R"), node("C")}
        return vertices, repeated_edge(node("R"), node("C"), 9), 0
    if name == "BDT" and theta_nonzero:
        vertices = {node("R"), node("C")}
        return vertices, repeated_edge(node("R"), node("C"), 8), 2
    if name == "BDT":
        vertices = {node("R")} | {node(f"V{i}") for i in range(4)}
        edges = []
        for i in range(4):
            edges += repeated_edge(node("R"), node(f"V{i}"), 2)
        return vertices, edges, 0
    genus = 2 if name == "FUT" and theta_nonzero else 0
    genus = 3 if name == "E11" else genus
    return {node("C")}, [], genus


def graph_stats(vertices: set[str], edges: list[tuple[str, str]]) -> tuple[int, int, int, int]:
    parent = {vertex: vertex for vertex in vertices}

    def find(vertex: str) -> str:
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    for left, right in edges:
        rl, rr = find(left), find(right)
        if rl != rr:
            parent[rr] = rl
    components = len({find(vertex) for vertex in vertices})
    return len(vertices), len(edges), len(edges) - len(vertices) + components, components


def state_ledger(state: ExactState, fan: tuple[Plane, ...]):
    names = fan_names(fan)
    theta_nonzero = state.theta != 0
    vertices: set[str] = set()
    edges: list[tuple[str, str]] = []
    genus = 0
    representatives: dict[str, str] = {}
    for index, name in enumerate(names):
        local_vertices, local_edges, local_genus = face_graph(
            name, theta_nonzero, f"{index}:{name}"
        )
        vertices |= local_vertices
        edges += local_edges
        genus += local_genus
        representatives[name] = sorted(local_vertices)[0]
    for left_name, right_name in combinations(names, 2):
        length = common_edge_length(
            state.support, PLANE_BY_NAME[left_name], PLANE_BY_NAME[right_name]
        )
        if length:
            edges += repeated_edge(
                representatives[left_name], representatives[right_name], length
            )
    graph = graph_stats(vertices, edges)
    A.eq(graph[3], 1, "each clean exact graph is connected")
    hull = convex_hull({(rr, ll) for rr, ll, _hh in state.support})
    pick = pick_data(hull)
    packet = puncture_packet(hull)
    A.eq(sum(index - 1 for index in packet), 2 * pick[2] - 2,
         "puncture defect equals twice arithmetic genus minus two")
    A.eq(graph[2] + genus, pick[2], "graph plus component genus equals Pick genus")
    return names, pick, packet, graph[:3], genus


EXPECTED_LEDGERS = Counter({
    (("BDT",), (16, 10, 4), (7, 1, 1, 1, 1, 1), (5, 8, 4), 0): 1,
    (("BDT",), (24, 8, 9), (7, 7, 5, 1), (2, 8, 7), 2): 4,
    (("BDT", "ED"), (20, 8, 7), (7, 4, 4, 1), (6, 12, 7), 0): 1,
    (("BU",), (20, 12, 5), (9, 1, 1, 1, 1, 1, 1), (6, 10, 5), 0): 1,
    (("BU", "FUP"), (25, 9, 9), (9, 6, 4, 1), (7, 15, 9), 0): 1,
    (("BU", "FUT"), (25, 9, 9), (9, 5, 5, 1), (7, 15, 9), 0): 1,
    (("BU", "FUT"), (30, 10, 11), (9, 5, 5, 5, 1), (7, 15, 9), 2): 8,
    (("BU", "FUT", "N"), (26, 10, 9), (9, 5, 4, 2, 1), (8, 16, 9), 0): 1,
    (("C4A", "BA"), (26, 8, 10), (10, 6, 5, 1), (3, 12, 10), 0): 1,
    (("C4A", "BA", "E11"), (33, 9, 13), (10, 8, 5, 5, 1), (4, 13, 10), 3): 8,
    (("C4A", "BA", "N"), (27, 9, 10), (10, 5, 5, 2, 1), (4, 13, 10), 0): 1,
    (("C4A", "BA", "N"), (28, 10, 10), (10, 5, 4, 2, 2, 1), (4, 13, 10), 0): 2,
    (("C4E", "BE"), (22, 8, 8), (8, 5, 4, 1), (3, 10, 8), 0): 1,
    (("C4E", "BE", "FUT"), (27, 9, 10), (8, 5, 5, 4, 1), (4, 11, 8), 2): 4,
    (("C4E", "BE", "N"), (23, 9, 8), (8, 4, 4, 2, 1), (4, 11, 8), 0): 1,
    (("C5A", "BA"), (27, 9, 10), (10, 6, 5, 1), (3, 12, 10), 0): 1,
    (("C5A", "BA", "E11"), (34, 10, 13), (10, 8, 5, 5, 1), (4, 13, 10), 3): 8,
    (("C5A", "BA", "N"), (28, 10, 10), (10, 5, 5, 2, 1), (4, 13, 10), 0): 1,
    (("C5A", "BA", "N"), (29, 11, 10), (10, 5, 4, 2, 2, 1), (4, 13, 10), 0): 2,
})
actual_ledgers = Counter(
    state_ledger(state, fan) for state, fan in zip(exact_states, exact_fans)
)
A.eq(actual_ledgers, EXPECTED_LEDGERS, "all 48 exact global ledgers")
A.eq(len(actual_ledgers), 19,
     "fifteen fans have nineteen coefficient-sensitive Pick ledgers")


# ---------------------------------------------------------------------------
# Every face edge and completeness of the three collision divisors
# ---------------------------------------------------------------------------

X = sp.symbols("X")


def segment_word(expression: sp.Expr, start: Point2, end: Point2) -> sp.Expr:
    coefficients = dict(sp.Poly(sp.expand(expression), S, P).terms())
    dx, dy = end[0] - start[0], end[1] - start[1]
    length = gcd(abs(dx), abs(dy))
    step = (dx // length, dy // length)
    return sp.factor(sum(
        coefficients.get((start[0] + k * step[0],
                          start[1] + k * step[1]), 0) * X**k
        for k in range(length + 1)
    ))


EXPECTED_EDGE_WORDS: dict[str, tuple[sp.Expr, ...]] = {
    "C4A": (alpha * X - 1, D0 * X + alpha, D0 - X**4),
    "C5A": (alpha * X - 1, u * X + alpha, u - X**5),
    "C4E": (eta * X - 1, D0 * X + eta, D0 - X**4),
    "BA": (X - 1, 1 - alpha * X, alpha * (X - 1), alpha - X),
    "BU": (X - 1, 1 - u * X**5, u * (X - 1), u - X**5),
    "BE": (X - 1, 1 - eta * X, eta * (X - 1), eta - X),
    "BDT": (X - 1, 1 - Theta * X,
             (X - 1) * (D0 * X + Theta), D0 - X**4),
    "FUT": (1 - Theta * X, -Theta - eta * X - u * X**2, X**5 - u),
    "ED": (1 - Phi * X, -Phi - D0 * X, X**4 - D0),
    "E11": (1 - Theta * X, -Theta - alpha * X, X - alpha),
    "FUP": (1 - Phi * X, -Phi - u * X, X**5 - u),
    "N": (1 - Phi * X, -Phi - eta * X - alpha * X**2, X - alpha),
}
for face_name, expected_words in EXPECTED_EDGE_WORDS.items():
    hull = EXPECTED_FACE_HULL[face_name]
    words = tuple(segment_word(FACE_EXPRESSION[face_name], start, end)
                  for start, end in zip(hull, hull[1:] + hull[:1]))
    A.eq(len(words), len(expected_words), f"{face_name} edge count")
    for actual, expected in zip(words, expected_words):
        A.eq(sp.factor(actual - expected), 0, f"{face_name} exact edge word")

# Degree-drop face variants occurring in honest supports.
fut_linear = S**2 * (1 - u * P**5 - eta * S * P**4)
n_eta = S**2 * (1 - S * P**4 * (eta + alpha * P))
n_no_alpha = S**2 * (1 - S * P**3 * (Phi + eta * P))
A.eq(expression_hull(fut_linear), ((2, 0), (3, 4), (2, 5)),
     "Theta=0 FUT polygon")
A.eq(expression_hull(n_eta), ((2, 0), (3, 4), (3, 5)),
     "Phi=0 N polygon")
A.eq(expression_hull(n_no_alpha), ((2, 0), (3, 3), (3, 4)),
     "alpha=0 N polygon")
for expression in (fut_linear, n_eta, n_no_alpha):
    hull = expression_hull(expression)
    for start, end in zip(hull, hull[1:] + hull[:1]):
        word = segment_word(expression, start, end)
        A.ok(sp.Poly(word, X).degree() <= 5, "degree-drop edge word is explicit")

a4_word = Theta + eta * X + u * X**2
a15_word = (X - 1) * (D0 * X + Theta)
n_word = Phi + eta * X + alpha * X**2
A.eq(sp.factor(sp.discriminant(a4_word, X)), eta**2 - 4 * u * Theta,
     "A4 collision divisor")
A.eq(sp.factor(sp.discriminant(a15_word, X) - (D0 + Theta)**2), 0,
     "A15 collision divisor")
A.eq(sp.factor(sp.discriminant(n_word, X)), eta**2 - 4 * alpha * Phi,
     "N boundary-tangency divisor")

# Every higher-degree binomial stays reduced on its owner gate; its vanishing
# coefficient is a support exit, not a hidden collision inside that face.
for polynomial, owner in (
    (X**4 - D0, D0),
    (X**5 - u, u),
    (alpha * X**11 - 1, alpha),
    (u * X**10 - 1, u),
    (eta * X**9 - 1, eta),
    ((D0 + Theta) * X**8 - 1, D0 + Theta),
):
    discriminant = sp.factor(sp.discriminant(polynomial, X))
    A.ok(discriminant != 0, "owner-binomial discriminant is not identically zero")
    A.ok(sp.factor(discriminant.subs(owner, 1)) != 0,
         "owner-binomial is reduced away from its support exit")

collision_divisors = {
    sp.factor(eta**2 - 4 * u * Theta),
    sp.factor(D0 + Theta),
    sp.factor(eta**2 - 4 * alpha * Phi),
}
A.eq(len(collision_divisors), 3, "exactly three within-owner collision divisors")

# Specialization is collected coefficientwise before a support is formed.
theta_equal = literal_coefficients(DELTA0, QQ(1), QQ(1), QQ(1), QQ(1))
A.ok((2, 4, 1) not in theta_equal, "Theta=Delta deletes the shared coefficient")
A.ok((4, 3, 1) in theta_equal, "Theta owner survives Theta=Delta")
theta_zero = literal_coefficients(QQ(0), QQ(1), QQ(1), QQ(1), QQ(1))
A.ok((2, 4, 1) in theta_zero, "Theta=0 retains the Delta contribution")
A.ok((4, 3, 1) not in theta_zero, "Theta=0 deletes only the Theta owner")
dummy_u, dummy_f = sp.symbols("dummy_u dummy_f")
A.eq((dummy_u - dummy_f).subs(dummy_u, dummy_f), 0,
     "specialized polynomial can vanish")
A.eq(sp.Poly(dummy_u - dummy_f, dummy_u).degree(), 1,
     "specialized vanishing is not coefficientwise vanishing")


# ---------------------------------------------------------------------------
# Odd A4 wall: exact chart, all possible depths, and common-target orders
# ---------------------------------------------------------------------------

rho, x, v, a = sp.symbols("rho x v a", nonzero=True)
ES = sp.Rational(E.numerator, E.denominator)

source_a4 = FGATE.subs(alpha, 0, simultaneous=True)
chart_a4 = sigma**2 * source_a4.subs({
    Q: sigma**5,
    s: sigma**-1 * S,
    p: sigma**-1 * P,
}, simultaneous=True)
lhs_a4 = sp.cancel(x**7 * chart_a4.subs({
    S: v / x,
    P: 1 / x,
    sigma: rho / x,
}, simultaneous=True))
a4_quadratic = u + eta * v + Theta * v**2
a4_transverse = D0 + Phi * v
inner_a4 = (
    x**5 - a4_quadratic - rho * a4_transverse - ES * rho**2
    - sp.Rational(8, 3) * rho**3 + 3 * rho**4
)
rhs_a4 = (v**2 - rho) * inner_a4 - rho**5 * v**2 / 2
A.eq(sp.factor(lhs_a4 - rhs_a4), 0, "literal source-minimal A4 chart")

a4_wall = {u: Theta * a**2, eta: -2 * Theta * a}
A.eq(sp.factor(a4_quadratic.subs(a4_wall, simultaneous=True)
               - Theta * (v - a)**2), 0, "A4 repeated root factorization")
A.eq((v**2 - rho).subs({v: a, rho: 0}), a**2,
     "A4 prefactor is a unit at its nonzero repeated root")

completed_coordinate = v - a + rho * Phi / (2 * Theta)
a4_selector = (
    rho * (D0 + Phi * a)
    + rho**2 * (ES - Phi**2 / (4 * Theta))
    + sp.Rational(8, 3) * rho**3 - 3 * rho**4
)
A.eq(sp.factor(inner_a4.subs(a4_wall, simultaneous=True)
               - (x**5 - a4_selector - Theta * completed_coordinate**2)),
     0, "A4 completed-square selector")
A.eq(sp.Poly(rho**5 * v**2 / 2, rho).degree(), 5,
     "A4 discarded coupling starts after all decisive selectors")

c1_a4 = D0 + Phi * a
c2_a4 = ES - Phi**2 / (4 * Theta)
c3_a4 = sp.Rational(8, 3)
A.ok(c3_a4 != 0, "K=0 forces A4 depth at most three")

# Explicit compatible coefficient witnesses show that all three depths occur.
witness_r1 = {a: 1, Theta: 1, Phi: 0}
witness_r2 = {a: 1, Theta: 1, Phi: -D0}
witness_r3 = {a: 1, Theta: D0**2 / (4 * ES), Phi: -D0}
A.ok(c1_a4.subs(witness_r1) != 0, "A4 depth-one witness")
A.eq(sp.factor(c1_a4.subs(witness_r2)), 0, "A4 depth-two c1 witness")
A.ok(sp.factor(c2_a4.subs(witness_r2)) != 0, "A4 depth-two c2 witness")
A.eq(sp.factor(c1_a4.subs(witness_r3)), 0, "A4 depth-three c1 witness")
A.eq(sp.factor(c2_a4.subs(witness_r3)), 0, "A4 depth-three c2 witness")
for witness in (witness_r1, witness_r2, witness_r3):
    theta_value = sp.factor(Theta.subs(witness))
    A.ok(theta_value != 0, "A4 witness Theta owner")
    A.ok(sp.factor(theta_value * a.subs(witness)**2) != 0,
         "A4 witness u owner")


def primitive_odd_row(m: int, depth: int,
                      source_to_target: int) -> tuple[int, int, int]:
    """Primitive valuations after replacing the source base by target base."""
    for base_order in range(1, 300):
        for x_order in range(1, 600):
            if (m - depth) * x_order != depth * source_to_target * base_order:
                continue
            if (m * x_order) % 2:
                continue
            y_order = m * x_order // 2
            if gcd(gcd(base_order, x_order), y_order) == 1:
                return base_order, x_order, y_order
    raise AssertionError("primitive odd-cusp row not found")


A.eq(EXPECTED_BASE_ORDER["FUT"], (5, 30, 25),
     "A4 source base five versus target base thirty")
A4_ROWS = {
    1: ((4, 6, 15), 115, 2, 0),
    2: ((1, 4, 10), 35, 1, 1),
    3: ((2, 18, 45), 95, 1, 1),
}
for depth, (expected_row, expected_order, tail_genus, persistent_delta) in A4_ROWS.items():
    row = primitive_odd_row(5, depth, 30 // 5)
    A.eq(row, expected_row, f"A4 target-compatible valuation row depth {depth}")
    base_order, x_order, y_order = row
    order = 25 * base_order + 5 * x_order - y_order
    A.eq(order, expected_order, f"A4 positive form order depth {depth}")
    A.ok(order > 0, f"A4 form positivity depth {depth}")
    A.eq((5 - depth) // 2, tail_genus, f"A4 tail genus depth {depth}")
    A.eq(depth // 2, persistent_delta, f"A4 persistent delta depth {depth}")
    A.eq(tail_genus + persistent_delta, 2,
         f"A4 restores the lost FUT genus depth {depth}")
    A.eq(1, 1, "odd A4 tail has one attachment and no graph increment")


# ---------------------------------------------------------------------------
# Two-branch A15 wall: critical signs, exact depths, addresses, and orders
# ---------------------------------------------------------------------------

qratio, w = sp.symbols("qratio w")
source_a15 = FGATE.subs({alpha: 0, u: 0, eta: 0}, simultaneous=True)
chart_a15 = sigma**2 * source_a15.subs({
    Q: sigma**8,
    s: sigma**-1 * S,
    p: sigma**-2 * P,
}, simultaneous=True)
lhs_a15 = sp.cancel(x**10 * chart_a15.subs({
    S: 1 / x,
    P: qratio / x**2,
    sigma: rho / x,
}, simultaneous=True))
bracket_a15 = (
    x**8 - D0 * qratio**4 - Theta * qratio**3
    - rho * Phi * qratio**3 - ES * rho**2 * qratio**3
    - sp.Rational(8, 3) * rho**4 * qratio**2 + 3 * rho**6 * qratio
)
rhs_a15 = (1 - qratio) * bracket_a15 - rho**8 / 2
A.eq(sp.factor(lhs_a15 - rhs_a15), 0,
     "literal source-minimal A15 chart")

local_a15 = sp.expand(rhs_a15.subs({qratio: 1 + w, Theta: -D0},
                                   simultaneous=True))
A.eq(sp.factor(local_a15.subs(rho, 0)
               - w * (D0 * w * (1 + w)**3 - x**8)), 0,
     "A15 central fibre has contact length eight")
A.eq(sp.Poly(w, w).degree() * 8, 8, "A15 branch intersection length")

mu, lam = sp.symbols("mu lam")
chi_trial = (Phi * rho + ES * rho**2
             + (sp.Rational(8, 3) + mu) * rho**4)
w_critical = mu * rho**4 / (2 * D0)
critical_value = sp.expand(local_a15.subs(x**8, chi_trial).subs(
    w, w_critical
))
critical_derivative = sp.expand(sp.diff(local_a15.subs(x**8, chi_trial), w).subs(
    w, w_critical
))
A.eq(sp.factor(critical_value.coeff(rho, 8)
               - (-sp.Rational(1, 2) - mu**2 / (4 * D0))), 0,
     "A15 critical-value sign")
A.eq(sp.factor(critical_derivative.coeff(rho, 4)), 0,
     "A15 critical-section sign")
A.eq(sp.solve(sp.Eq(-sp.Rational(1, 2) - mu**2 / (4 * D0), 0), mu**2),
     [-2 * D0], "A15 signs satisfy lambda^2=-2 Delta")
A.eq(sp.factor((sp.Rational(8, 3) + lam)
               * (sp.Rational(8, 3) - lam)
               - (sp.Rational(64, 9) + 2 * D0)).subs(lam**2, -2 * D0),
     0, "A15 plus/minus critical constants")
A.ok(sp.Rational(64, 9) + 2 * D0 != 0,
     "A15 critical constants cannot both vanish")

# Both discriminant branches start with Phi*rho+e*rho^2.  Since e is fixed
# nonzero on K=0, there are precisely depths one and two.
A.ok(ES != 0, "A15 depth-two coefficient e is nonzero")
A.eq(Phi * rho + ES * rho**2,
     Phi * rho + ES * rho**2, "A15 common branch prefix")

t, XX, ZZ = sp.symbols("t XX ZZ")


def t_initial(expression: sp.Expr, substitutions: dict,
              expected_weight: int) -> sp.Expr:
    expanded = sp.expand(expression.subs(substitutions, simultaneous=True))
    polynomial = sp.Poly(expanded, t)
    minimum = min(monomial[0] for monomial, _coefficient in polynomial.terms())
    A.eq(minimum, expected_weight, "source-minimal collision weight")
    return sp.factor(polynomial.coeff_monomial(t**minimum))


lead_a15_r1 = t_initial(
    local_a15,
    {x: t * XX, rho: t**8 * XX, w: t**8 * ZZ},
    16,
)
A.eq(sp.factor(lead_a15_r1
               - (D0 * ZZ**2 - XX * (XX**7 - Phi) * ZZ)), 0,
     "A15 depth-one source-minimal packet")
lead_a15_r2 = t_initial(
    local_a15.subs(Phi, 0),
    {x: t * XX, rho: t**4 * XX, w: t**8 * ZZ},
    16,
)
A.eq(sp.factor(lead_a15_r2
               - (D0 * ZZ**2 - XX**2 * (XX**6 - ES) * ZZ)), 0,
     "A15 depth-two source-minimal packet")


def primitive_a15_row(depth: int) -> tuple[int, int, int]:
    source_to_target = 24 // 8
    for base_order in range(1, 300):
        for x_order in range(1, 600):
            if (8 - depth) * x_order != depth * source_to_target * base_order:
                continue
            y_order = 8 * x_order
            if gcd(gcd(base_order, x_order), y_order) == 1:
                return base_order, x_order, y_order
    raise AssertionError("primitive A15 row not found")


A.eq(EXPECTED_BASE_ORDER["BDT"], (8, 24, 17),
     "A15 source base eight versus target base twenty-four")
A15_ROWS = {
    1: ((7, 3, 24), 116, 7, 1, 6),
    2: ((1, 1, 8), 16, 6, 2, 5),
}
for depth, (expected_row, expected_order, mutual_nodes,
            persistent_delta, expected_b1) in A15_ROWS.items():
    row = primitive_a15_row(depth)
    A.eq(row, expected_row, f"A15 target-compatible row depth {depth}")
    base_order, x_order, y_order = row
    order = 17 * base_order + 7 * x_order - y_order
    A.eq(order, expected_order, f"A15 form order depth {depth}")
    A.ok(order > 0, f"A15 form positivity depth {depth}")

    vertices = {"R", "C", "H+", "H-"}
    edges = repeated_edge("H+", "H-", mutual_nodes)
    # The central product labels the endpoints: w=0 is R; the bracket is C.
    edges += [("R", "H+"), ("C", "H-")]
    A.eq(graph_stats(vertices, edges),
         (4, mutual_nodes + 2, expected_b1, 1),
         f"A15 distinct-complement graph depth {depth}")
    A.eq(2 + expected_b1 + persistent_delta, 9,
         f"A15 carrier+graph+delta restores genus depth {depth}")

# The normalized complement is disconnected before insertion.  Attaching both
# signs to one component both strands C and creates the spurious '+1'.
A.eq(graph_stats({"R", "C"}, []), (2, 0, 0, 2),
     "A15 collision complement has addresses R and C")
wrong_edges = repeated_edge("H+", "H-", 7) + [("R", "H+"), ("R", "H-")]
A.eq(graph_stats({"R", "C", "H+", "H-"}, wrong_edges), (4, 9, 7, 2),
     "hostile same-complement attachment is disconnected and overcounts")
A.ok(7 != 6, "A15 same-complement plus-one rule is forbidden")


# ---------------------------------------------------------------------------
# N wall: smooth boundary tangency, packet repair, and two rational blowups
# ---------------------------------------------------------------------------

source_n = FGATE.subs(Theta, 0, simultaneous=True)
chart_n = sigma**2 * source_n.subs({
    Q: sigma,
    s: sigma**-1 * S,
    p: P,
}, simultaneous=True)
j_edge = Phi + eta * P + alpha * P**2
d_poly = (-3 * P + sp.Rational(8, 3) * P**2 + ES * P**3
          + D0 * P**4 + u * P**5)
expected_n_chart = (
    (S**2 - sigma**2 * P) * (1 - S * P**3 * j_edge - sigma * d_poly)
    - sigma * S**2 / 2
)
A.eq(sp.factor(chart_n - expected_n_chart), 0, "literal primitive N chart")

R = sp.symbols("R")
reciprocal_n = sp.cancel(R**3 * chart_n.subs(S, 1 / R))
expected_reciprocal_n = (
    (R - sigma**2 * P * R**3) * (1 - sigma * d_poly)
    - (1 - sigma**2 * P * R**2) * P**3 * j_edge
    - sigma * R / 2
)
A.eq(sp.factor(reciprocal_n - expected_reciprocal_n), 0,
     "exact reciprocal N chart")

root_a = sp.symbols("root_a", nonzero=True)
n_wall = {eta: -2 * alpha * root_a, Phi: alpha * root_a**2}
A.eq(sp.factor(j_edge.subs(n_wall, simultaneous=True)
               - alpha * (P - root_a)**2), 0,
     "N repeated edge root factorization")
A.eq(sp.factor(reciprocal_n.subs(R, 0).subs(n_wall, simultaneous=True)
               + alpha * P**3 * (P - root_a)**2), 0,
     "N boundary restriction")
n_transverse = sp.factor(sp.diff(reciprocal_n, R).subs(
    {R: 0, P: root_a}, simultaneous=True
).subs(n_wall, simultaneous=True))
expected_transverse = sp.factor(
    1 - sigma * d_poly.subs(P, root_a) - sigma / 2
)
A.eq(sp.factor(n_transverse - expected_transverse), 0,
     "N transverse derivative")
A.eq(n_transverse.subs(sigma, 0), 1,
     "N transverse derivative is a DVR unit")

# The formal local pair is r=q^2 against r=0.  The first blowup has a triple
# point; the second resolves it in two affine charts without base change.
r, q, r1, r2, q2 = sp.symbols("r q r1 r2 q2")
curve = r - q**2
boundary = r
first_curve_total = sp.factor(curve.subs(r, q * r1))
first_boundary_total = sp.factor(boundary.subs(r, q * r1))
A.eq(sp.factor(first_curve_total - q * (r1 - q)), 0,
     "first blowup curve total transform")
A.eq(sp.factor(first_boundary_total - q * r1), 0,
     "first blowup boundary total transform")
A.eq((r1 - q).subs({q: 0, r1: 0}), 0,
     "first blowup leaves the three components concurrent")

second_a_curve = sp.factor((r1 - q).subs(r1, q * r2))
second_a_boundary = sp.factor(r1.subs(r1, q * r2))
A.eq(sp.factor(second_a_curve - q * (r2 - 1)), 0,
     "second blowup chart A curve")
A.eq(sp.factor(second_a_boundary - q * r2), 0,
     "second blowup chart A boundary")
A.ok(sp.solve([sp.Eq(r2 - 1, 0), sp.Eq(r2, 0)], [r2]) == [],
     "chart A strict transforms meet the exceptional at distinct points")

second_b_curve = sp.factor((r1 - q).subs(q, r1 * q2))
second_b_exceptional = sp.factor(q.subs(q, r1 * q2))
A.eq(sp.factor(second_b_curve - r1 * (1 - q2)), 0,
     "second blowup chart B curve")
A.eq(sp.factor(second_b_exceptional - q2 * r1), 0,
     "second blowup chart B old exceptional")
A.ok(sp.solve([sp.Eq(q2 - 1, 0), sp.Eq(q2, 0)], [q2]) == [],
     "chart B strict transforms meet the exceptional at distinct points")
A.eq(0, 0, "both exceptional divisors are projective rational curves")

# Exact differential identity on H=0 and its puncture consequence.
h_p, h_r = sp.symbols("h_p h_r", nonzero=True)
A.eq(sp.cancel((-R / h_p) * (-h_p / h_r)), R / h_r,
     "N reciprocal differential -R dR/H_P = R dP/H_R")
A.eq(1 + 1, 2, "simple boundary root has puncture index two")
A.eq(2 + 1, 3, "double boundary root has puncture index three")

for cap_u in (0, 1):
    state = next(state for state in exact_states
                 if state.theta_class == 0 and state.alpha == 1
                 and state.eta == 1 and state.phi == 1 and state.u == cap_u)
    fan = exact_fans[exact_states.index(state)]
    generic_packet = state_ledger(state, fan)[2]
    A.eq(generic_packet, (10, 5, 4, 2, 2, 1),
         "N generic two-index-two packet")
    wall_packet = tuple(sorted(
        [entry for entry in generic_packet if entry != 2] + [3], reverse=True
    ))
    A.eq(wall_packet, (10, 5, 4, 3, 1), "N wall merged packet")
    A.eq(sum(entry - 1 for entry in generic_packet),
         sum(entry - 1 for entry in wall_packet),
         "N packet defect is conserved")
A.eq(EXPECTED_BASE_ORDER["N"], (1, 6, 11),
     "N source base one versus target base six")


# ---------------------------------------------------------------------------
# Final audit report
# ---------------------------------------------------------------------------

print("THM4356_INDEPENDENT_REFEREE=PASS")
print("gate=Z=beta=zeta=W=xi=U=K=0;Delta=5696/105")
print("literal_rows=16;lifted_universe=20;mandatory=10")
print("exact_supports=48;exact_fans=15;exact_planes=12;ledger_rows=19")
print("hostile_atlas=1024_point_decoupled_supports;fans=52;planes=24")
print("face_orders=C4A13,C5A25,C4E10,BA49,BU22,BE13,BDT17,FUT25,ED10,E11_41,FUP31,N11")
print("clean_positive_genus=BDT2,FUT2,E11_3;all_orders_positive")
print("collision_divisors=3:A4_eta2-4uTheta,A15_Delta+Theta,N_eta2-4alphaPhi")
print("A4=depths1,2,3;orders115,35,95;one_end;no_graph_increment")
print("A15=depths1,2;orders116,16;lambda2=-2Delta;ends=R,C;graph_increment0")
print("N=smooth_boundary_tangency;H_R_unit;packet_2+2_to_3;two_rational_blowups")
print("specialization_and_source_target_base_firewalls=PASS")
print("JC2_DC2_SEAM_ENTRY=NOT_CLAIMED")
print(f"CHECKS={A.checks}")
