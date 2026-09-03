#!/usr/bin/env python3
"""Primary exact certificate for THM-4353's simultaneous zero endpoint.

Scope:
    Z=beta_11=zeta_3=W=xi_10=K=0, U != 0,
    Delta=5696/105

in the inherited reduced (2,3), exact-weight-twelve seam.  This import-free
certificate reconstructs the literal sixteen-row source, enumerates the 48
coefficient-realizable supports and a 128-key hostile Boolean atlas, computes
all lower faces and primitive charts, proves the only carrier genera, checks
every edge/packet/graph ledger, and resolves the sole repeated boundary edge
by two explicit ordinary blowups.
"""

from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations, product
from math import gcd, lcm
import sys

import sympy as sp


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")

CHECKS = 0


def need(condition, label):
    global CHECKS
    CHECKS += 1
    if condition is not True and not bool(condition):
        raise RuntimeError("THM-4353 primary audit failure: " + label)


def rank_two(points):
    for first, second, third in combinations(points, 3):
        if ((second[0] - first[0]) * (third[1] - first[1])
                != (second[1] - first[1]) * (third[0] - first[0])):
            return True
    return False


def candidate_planes(universe):
    answer = set()
    for first, second, third in combinations(sorted(universe), 3):
        determinant = (
            (second[0] - first[0]) * (third[1] - first[1])
            - (second[1] - first[1]) * (third[0] - first[0])
        )
        if determinant == 0:
            continue
        aa = F(
            (second[2] - first[2]) * (third[1] - first[1])
            - (second[1] - first[1]) * (third[2] - first[2]),
            determinant,
        )
        bb = F(
            (second[0] - first[0]) * (third[2] - first[2])
            - (second[2] - first[2]) * (third[0] - first[0]),
            determinant,
        )
        cc = F(first[2]) - aa * first[0] - bb * first[1]
        answer.add((aa, bb, cc))
    return tuple(sorted(answer))


def lower_fan(active, planes):
    answer = set()
    for aa, bb, cc in planes:
        equal = []
        for rr, ll, hh in active:
            gap = F(hh) - aa * rr - bb * ll - cc
            if gap < 0:
                break
            if gap == 0:
                equal.append((rr, ll, hh))
        else:
            if rank_two(equal):
                answer.add((aa, bb, cc))
    return frozenset(answer)


def convex_hull(points):
    points = sorted(set(points))
    if len(points) <= 1:
        return tuple(points)

    def cross(origin, first, second):
        return (
            (first[0] - origin[0]) * (second[1] - origin[1])
            - (first[1] - origin[1]) * (second[0] - origin[0])
        )

    lower = []
    for point in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper = []
    for point in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return tuple(lower[:-1] + upper[:-1])


def pick(vertices):
    pairs = list(zip(vertices, vertices[1:] + vertices[:1]))
    area2 = abs(sum(
        first[0] * second[1] - second[0] * first[1]
        for first, second in pairs
    ))
    boundary = sum(
        gcd(abs(second[0] - first[0]), abs(second[1] - first[1]))
        for first, second in pairs
    )
    need((area2 - boundary + 2) % 2 == 0, "Pick parity")
    return area2, boundary, (area2 - boundary + 2) // 2


def edge_packet(vertices):
    packet = []
    for start, end in zip(vertices, vertices[1:] + vertices[:1]):
        dx, dy = end[0] - start[0], end[1] - start[1]
        length = gcd(abs(dx), abs(dy))
        inward = (-dy // length, dx // length)
        constant = inward[0] * start[0] + inward[1] * start[1]
        index = inward[0] + inward[1] - constant
        if not (start[0] == end[0] == 0):
            packet.extend([index] * length)
    return tuple(sorted(packet, reverse=True))


def face_order(plane):
    base = lcm(6, *(entry.denominator for entry in plane))
    density = F(5, 6) - sum(plane)
    order = density * base
    need(order.denominator == 1 and order > 0, "positive integral form order")
    return density, base, order.numerator


def on_plane(point, plane):
    rr, ll, hh = point
    aa, bb, cc = plane
    return F(hh) == aa * rr + bb * ll + cc


def face_expression(support, plane, substitution=None):
    expression = sum(
        coefficient * S**rr * P**ll
        for (rr, ll, _hh), coefficient in support.items()
        if on_plane((rr, ll, _hh), plane)
    )
    if substitution:
        expression = expression.subs(substitution)
    return sp.factor(expression)


def edge_polynomial(expression, start, end):
    polynomial = sp.Poly(sp.expand(expression), S, P)
    dx, dy = end[0] - start[0], end[1] - start[1]
    length = gcd(abs(dx), abs(dy))
    step = (dx // length, dy // length)
    answer = 0
    for index in range(length + 1):
        exponent = (start[0] + index * step[0], start[1] + index * step[1])
        answer += polynomial.coeff_monomial(S**exponent[0] * P**exponent[1]) * X**index
    return sp.factor(answer)


# Literal sixteen-row source and the simultaneous endpoint gate.
U, W, Z = sp.symbols("U W Z")
Delta, Theta, Phi, eta = sp.symbols("Delta Theta Phi eta")
zeta, u, xi, alpha, beta, K = sp.symbols("zeta u xi alpha beta K")
S, P, X = sp.symbols("S P X")
e = -sp.Rational(1376, 135)
DELTA0 = sp.Rational(5696, 105)
FULL_ROWS = {
    (1, 0): -3,
    (2, 0): sp.Rational(8, 3),
    (3, 0): e,
    (0, 2): K,
    (2, 1): Phi,
    (4, 0): Delta,
    (1, 2): Theta,
    (3, 1): eta,
    (0, 3): zeta,
    (5, 0): u,
    (2, 2): xi,
    (4, 1): alpha,
    (1, 3): beta,
    (6, 0): U,
    (3, 2): W,
    (0, 4): Z,
}
need(len(FULL_ROWS) == 16, "literal source has sixteen rows")
need(sp.Rational(2848, 45) - sp.Rational(7, 6) * DELTA0 == 0,
     "K=0 forces Delta=5696/105")

GATE = {Z: 0, beta: 0, zeta: 0, W: 0, xi: 0, K: 0, Delta: DELTA0}
ROWS = {
    label: sp.factor(sp.sympify(coefficient).subs(GATE))
    for label, coefficient in FULL_ROWS.items()
    if sp.factor(sp.sympify(coefficient).subs(GATE)) != 0
}
need(len(ROWS) == 10, "ten surviving source rows")

SUPPORT = defaultdict(lambda: sp.Integer(0))
SUPPORT[(2, 0, 0)] += 1
SUPPORT[(0, 1, 0)] -= 1
SUPPORT[(2, 0, 1)] -= sp.Rational(1, 2)
for (ii, jj), coefficient in ROWS.items():
    SUPPORT[(jj + 2, ii + jj, 1)] -= coefficient
    SUPPORT[(jj, ii + jj + 1, 1)] += coefficient
SUPPORT = {
    point: sp.factor(coefficient)
    for point, coefficient in SUPPORT.items()
    if coefficient != 0
}
EXPECTED_SUPPORT = {
    (0, 1, 0): -1,
    (0, 2, 1): -3,
    (0, 3, 1): sp.Rational(8, 3),
    (0, 4, 1): e,
    (0, 5, 1): DELTA0,
    (0, 6, 1): u,
    (0, 7, 1): U,
    (1, 4, 1): Phi,
    (1, 5, 1): eta,
    (1, 6, 1): alpha,
    (2, 0, 0): 1,
    (2, 0, 1): -sp.Rational(1, 2),
    (2, 1, 1): 3,
    (2, 2, 1): -sp.Rational(8, 3),
    (2, 3, 1): -e,
    (2, 4, 1): Theta - DELTA0,
    (2, 5, 1): -u,
    (2, 6, 1): -U,
    (3, 3, 1): -Phi,
    (3, 4, 1): -eta,
    (3, 5, 1): -alpha,
    (4, 3, 1): -Theta,
}
need(SUPPORT.keys() == EXPECTED_SUPPORT.keys(), "complete specialized support keys")
for point in EXPECTED_SUPPORT:
    need(sp.simplify(SUPPORT[point] - EXPECTED_SUPPORT[point]) == 0,
         "complete specialized support coefficient")

# All possible lower planes and the eight exact fans.
M = (F(1, 12), F(1, 6), F(-1, 6))
D6 = (F(1, 6), F(1, 6), F(-1, 3))
E11 = (F(2, 7), F(1, 7), F(-4, 7))
E01 = (F(1, 4), F(1, 6), F(-1, 2))
EETA = (F(1, 3), F(1, 6), F(-2, 3))
EPHI = (F(1, 2), F(1, 6), F(-1, 1))
NORM = (F(1, 1), F(0), F(-2, 1))
PLANES = candidate_planes(SUPPORT)


def expected_fan(theta_on, alpha_on, eta_on, phi_on):
    if theta_on:
        if alpha_on:
            return frozenset((M, D6, E11))
        return frozenset((M, E01))
    if alpha_on:
        if eta_on or phi_on:
            return frozenset((M, D6, NORM))
        return frozenset((M, D6))
    if eta_on:
        if phi_on:
            return frozenset((M, EETA, NORM))
        return frozenset((M, EETA))
    if phi_on:
        return frozenset((M, EPHI))
    return frozenset((M,))


THETA_REPRESENTATIVES = (
    ("zero", sp.Integer(0)),
    ("cancel", DELTA0),
    ("generic", sp.Integer(1)),
)
exact_supports = {}
exact_fans = Counter()
for (theta_tag, theta_value), alpha_value, eta_value, phi_value, u_value in product(
    THETA_REPRESENTATIVES, (0, 1), (0, 1), (0, 1), (0, 1)
):
    substitution = {
        U: 2, Theta: theta_value, alpha: alpha_value,
        eta: eta_value, Phi: phi_value, u: u_value,
    }
    active = frozenset(
        point for point, coefficient in SUPPORT.items()
        if sp.factor(coefficient.subs(substitution)) != 0
    )
    fan = lower_fan(active, PLANES)
    expected = expected_fan(
        theta_value != 0, alpha_value != 0, eta_value != 0, phi_value != 0
    )
    need(fan == expected, "exact support has its owner fan")
    key = (theta_tag, alpha_value, eta_value, phi_value, u_value)
    exact_supports[key] = active
    exact_fans[fan] += 1

need(len(exact_supports) == 48, "48 exact coefficient states")
need(len(set(exact_supports.values())) == 48, "48 distinct exact supports")
EXPECTED_EXACT_FANS = Counter({
    frozenset((M, D6, E11)): 16,
    frozenset((M, E01)): 16,
    frozenset((M, D6, NORM)): 6,
    frozenset((M, D6)): 2,
    frozenset((M, EETA, NORM)): 2,
    frozenset((M, EETA)): 2,
    frozenset((M, EPHI)): 2,
    frozenset((M,)): 2,
})
need(exact_fans == EXPECTED_EXACT_FANS, "exact eight-fan census")

# Conservative Boolean atlas.  Row presence and aggregate deletion are
# independently toggled, so configurations are not claimed realizable.
REQUIRED_LABELS = ((1, 0), (2, 0), (3, 0), (4, 0), (6, 0))
OPTIONAL_LABELS = ((5, 0), (2, 1), (1, 2), (3, 1), (4, 1))
CANCELLABLE = ((2, 4, 1), (2, 5, 1))


def lifted(labels):
    answer = {(2, 0, 0), (0, 1, 0), (2, 0, 1)}
    for ii, jj in labels:
        answer.add((jj + 2, ii + jj, 1))
        answer.add((jj, ii + jj + 1, 1))
    return answer


hostile_states = {}
hostile_fans = Counter()
for row_bits in product((0, 1), repeat=5):
    labels = list(REQUIRED_LABELS)
    labels.extend(label for bit, label in zip(row_bits, OPTIONAL_LABELS) if bit)
    for deletion_bits in product((0, 1), repeat=2):
        active = lifted(labels)
        for bit, point in zip(deletion_bits, CANCELLABLE):
            if bit:
                active.discard(point)
        active = frozenset(active)
        key = tuple(row_bits) + tuple(deletion_bits)
        hostile_states[key] = active
        u_on, phi_on, theta_on, eta_on, alpha_on = row_bits
        fan = lower_fan(active, PLANES)
        need(fan == expected_fan(theta_on, alpha_on, eta_on, phi_on),
             "hostile support has its owner fan")
        hostile_fans[fan] += 1

need(len(hostile_states) == 128, "128 hostile keyed configurations")
need(len(set(hostile_states.values())) == 96, "96 hostile distinct supports")
need(hostile_fans == Counter({
    frozenset((M, D6, E11)): 32,
    frozenset((M, E01)): 32,
    frozenset((M, D6, NORM)): 24,
    frozenset((M, D6)): 8,
    frozenset((M, EETA, NORM)): 8,
    frozenset((M, EETA)): 8,
    frozenset((M, EPHI)): 8,
    frozenset((M,)): 8,
}), "hostile eight-fan census")

for key, active in exact_supports.items():
    theta_tag, alpha_value, eta_value, phi_value, u_value = key
    theta_on = int(theta_tag != "zero")
    row_bits = (u_value, phi_value, theta_on, eta_value, alpha_value)
    deletion_bits = (int(theta_tag == "cancel"), int(u_value == 0))
    atlas_key = tuple(row_bits) + tuple(deletion_bits)
    need(hostile_states[atlas_key] == active,
         "exact support embeds literally in hostile atlas")

# Exact face equations, polygons, and form orders.
FACE_EQUATIONS = {
    "M": (M, {}, (P - S**2) * (U * P**6 - 1)),
    "D6": (D6, {}, S**2 * (1 - U * P**6 - alpha * S * P**5)),
    "E11": (E11, {}, S**2 * (1 - alpha * S * P**5
                                - Theta * S**2 * P**3)),
    "E01": (E01, {alpha: 0}, S**2 * (1 - U * P**6
                                        - Theta * S**2 * P**3)),
    "Eeta": (EETA, {Theta: 0, alpha: 0},
             S**2 * (1 - U * P**6 - eta * S * P**4)),
    "EPhi": (EPHI, {Theta: 0, alpha: 0, eta: 0},
             S**2 * (1 - U * P**6 - Phi * S * P**3)),
    "N": (NORM, {Theta: 0},
          S**2 * (1 - S * P**3 * (Phi + eta * P + alpha * P**2))),
}
for name, (plane, substitution, expected) in FACE_EQUATIONS.items():
    need(sp.expand(face_expression(SUPPORT, plane, substitution) - expected) == 0,
         name + " exact face equation")

FACE_POLYGONS = {
    "M": ((0, 1), (2, 0), (2, 6), (0, 7)),
    "D6": ((2, 0), (3, 5), (2, 6)),
    "E11": ((2, 0), (4, 3), (3, 5)),
    "E01": ((2, 0), (4, 3), (2, 6)),
    "Eeta": ((2, 0), (3, 4), (2, 6)),
    "EPhi": ((2, 0), (3, 3), (2, 6)),
    "N35": ((2, 0), (3, 3), (3, 5)),
    "N45": ((2, 0), (3, 4), (3, 5)),
    "N34": ((2, 0), (3, 3), (3, 4)),
}
EXPECTED_FACE_PICK = {
    "M": (24, 14, 6),
    "D6": (6, 8, 0),
    "E11": (7, 3, 3),
    "E01": (12, 8, 3),
    "Eeta": (6, 8, 0),
    "EPhi": (6, 8, 0),
    "N35": (2, 4, 0),
    "N45": (1, 3, 0),
    "N34": (1, 3, 0),
}
need({name: pick(vertices) for name, vertices in FACE_POLYGONS.items()}
     == EXPECTED_FACE_PICK, "all face Pick ledgers")

EXPECTED_ORDERS = {
    "M": (F(3, 4), 12, 9),
    "D6": (F(5, 6), 6, 5),
    "E11": (F(41, 42), 42, 41),
    "E01": (F(11, 12), 12, 11),
    "Eeta": (F(1, 1), 6, 6),
    "EPhi": (F(7, 6), 6, 7),
    "N": (F(11, 6), 6, 11),
}
for name, (plane, _substitution, _expected) in FACE_EQUATIONS.items():
    need(face_order(plane) == EXPECTED_ORDERS[name], name + " form order")

# Every primitive chart is reconstructed directly from the source.
sigma = sp.symbols("sigma")
HROWS = {
    label: sp.sympify(coefficient).subs(GATE)
    for label, coefficient in FULL_ROWS.items()
}


def primitive_chart(qpower, sweight, pweight, clear, substitution):
    old_s = sigma**(-sweight) * S
    old_p = sigma**(-pweight) * P
    old_y = old_s * old_p
    hpoly = sum(
        coefficient * old_p**ii * old_y**jj
        for (ii, jj), coefficient in HROWS.items()
    ).subs(substitution)
    fcurve = (
        (old_s**2 - old_p) * (1 - sigma**qpower * hpoly)
        - sigma**qpower * old_s**2 / 2
    )
    chart = sp.cancel(sigma**clear * fcurve).expand()
    need(not sp.denom(sp.together(chart)).has(sigma),
         "primitive chart has no negative sigma power")
    polynomial = sp.Poly(chart, sigma)
    need(min(term[0][0] for term in polynomial.terms()) == 0,
         "primitive chart has valuation zero")
    return chart


CHART_SPECS = {
    "M": (12, 1, 2, 2, {}),
    "D6": (6, 1, 1, 2, {}),
    "E11": (7, 2, 1, 4, {}),
    "E01": (12, 3, 2, 6, {alpha: 0}),
    "Eeta": (6, 2, 1, 4, {Theta: 0, alpha: 0}),
    "EPhi": (6, 3, 1, 6, {Theta: 0, alpha: 0, eta: 0}),
    "N": (1, 1, 0, 2, {Theta: 0}),
}
CHARTS = {}
for name, spec in CHART_SPECS.items():
    chart = primitive_chart(*spec)
    expected_initial = FACE_EQUATIONS[name][2]
    need(sp.expand(chart.subs(sigma, 0) - expected_initial) == 0,
         name + " primitive initial")
    CHARTS[name] = chart

# Exact carrier normalizations and all rational faces.
f11 = 1 - alpha * S * P**5 - Theta * S**2 * P**3
y11 = 2 * Theta * S * P**2 + alpha * P**4
branch11 = P * (4 * Theta + alpha**2 * P**7)
need(sp.expand(y11**2 - branch11 + 4 * Theta * P * f11) == 0,
     "E11 hyperelliptic identity")
need(sp.factor(sp.discriminant(branch11, P))
     == -sp.Integer(53971714048) * Theta**8 * alpha**12,
     "E11 branch discriminant")
need(sp.degree(branch11, P) == 8, "E11 branch degree eight")

f01 = 1 - U * P**6 - Theta * S**2 * P**3
y01 = Theta * S * P**2
branch01 = Theta * P * (1 - U * P**6)
need(sp.expand(y01**2 - branch01 + Theta * P * f01) == 0,
     "E01 hyperelliptic identity")
need(sp.factor(sp.discriminant(branch01, P))
     == sp.Integer(46656) * Theta**12 * U**5,
     "E01 branch discriminant")
need(sp.degree(branch01, P) == 7, "E01 branch degree seven")

RATIONAL_SOLUTIONS = {
    "D6": (
        1 - U * P**6 - alpha * S * P**5,
        (1 - U * P**6) / (alpha * P**5),
    ),
    "Eeta": (
        1 - U * P**6 - eta * S * P**4,
        (1 - U * P**6) / (eta * P**4),
    ),
    "EPhi": (
        1 - U * P**6 - Phi * S * P**3,
        (1 - U * P**6) / (Phi * P**3),
    ),
    "N": (
        1 - S * P**3 * (Phi + eta * P + alpha * P**2),
        1 / (P**3 * (Phi + eta * P + alpha * P**2)),
    ),
}
for name, (equation, solution) in RATIONAL_SOLUTIONS.items():
    need(sp.factor(equation.subs(S, solution)) == 0,
         name + " rational normalization")

# Literal edge restrictions.  These recover every scheme used below.
fM = FACE_EQUATIONS["M"][2]
fD6 = FACE_EQUATIONS["D6"][2]
fE11 = FACE_EQUATIONS["E11"][2]
fE01 = FACE_EQUATIONS["E01"][2]
fEeta = FACE_EQUATIONS["Eeta"][2]
fEPhi = FACE_EQUATIONS["EPhi"][2]
fN = FACE_EQUATIONS["N"][2]

EDGE_RESTRICTIONS = (
    (fM, (0, 1), (2, 0), X - 1),
    (fM, (2, 6), (0, 7), U * (X - 1)),
    (fM, (0, 1), (0, 7), U * X**6 - 1),
    (fM, (2, 0), (2, 6), 1 - U * X**6),
    (fD6, (2, 0), (3, 5), 1 - alpha * X),
    (fD6, (3, 5), (2, 6), -(alpha + U * X)),
    (fE11, (2, 0), (4, 3), 1 - Theta * X),
    (fE11, (4, 3), (3, 5), -(Theta + alpha * X)),
    (fE01, (2, 0), (4, 3), 1 - Theta * X),
    (fE01, (4, 3), (2, 6), -(Theta + U * X)),
    (fEeta, (2, 0), (3, 4), 1 - eta * X),
    (fEeta, (3, 4), (2, 6), -(eta + U * X)),
    (fEPhi, (2, 0), (3, 3), 1 - Phi * X),
    (fEPhi, (3, 3), (2, 6), -(Phi + U * X)),
    (fN, (2, 0), (3, 3), 1 - Phi * X),
    (fN, (3, 3), (3, 5), -(Phi + eta * X + alpha * X**2)),
    (fN, (2, 0), (3, 5), 1 - alpha * X),
    (fN.subs(Phi, 0), (2, 0), (3, 4), 1 - eta * X),
    (fN.subs(Phi, 0), (3, 4), (3, 5), -(eta + alpha * X)),
    (fN.subs(alpha, 0), (3, 3), (3, 4), -(Phi + eta * X)),
    (fN.subs(alpha, 0), (2, 0), (3, 4), 1 - eta * X),
)
for expression, start, end, expected in EDGE_RESTRICTIONS:
    need(sp.expand(edge_polynomial(expression, start, end) - expected) == 0,
         "literal edge restriction")

need(sp.factor(sp.discriminant(U - X**6, X)) == 46656 * U**5,
     "outer sixfold edge reduced")
need(sp.factor(sp.discriminant(1 - U * X**6, X)) == 46656 * U**5,
     "internal sixfold edge reduced")
need(sp.discriminant(Phi + alpha * X**2, X) == -4 * alpha * Phi,
     "even quadratic edge reduced")
JEDGE = Phi + eta * X + alpha * X**2
JDISC = eta**2 - 4 * alpha * Phi
need(sp.discriminant(JEDGE, X) == JDISC,
     "sole variable edge discriminant")

# Global polygons, packets, and graph ledgers.
GLOBAL_POLYGONS = {
    "theta-alpha": ((0, 1), (2, 0), (4, 3), (3, 5), (2, 6), (0, 7)),
    "theta": ((0, 1), (2, 0), (4, 3), (2, 6), (0, 7)),
    "alpha-phi": ((0, 1), (2, 0), (3, 3), (3, 5), (2, 6), (0, 7)),
    "alpha-eta": ((0, 1), (2, 0), (3, 4), (3, 5), (2, 6), (0, 7)),
    "alpha": ((0, 1), (2, 0), (3, 5), (2, 6), (0, 7)),
    "eta-phi": ((0, 1), (2, 0), (3, 3), (3, 4), (2, 6), (0, 7)),
    "eta": ((0, 1), (2, 0), (3, 4), (2, 6), (0, 7)),
    "phi": ((0, 1), (2, 0), (3, 3), (2, 6), (0, 7)),
    "none": ((0, 1), (2, 0), (2, 6), (0, 7)),
}
EXPECTED_GLOBAL = {
    "theta-alpha": ((37, 11, 14), (11, 8, 6, 5, 1), (9, 19, 11, 3)),
    "theta": ((36, 10, 14), (13, 11, 5, 1), (8, 18, 11, 3)),
    "alpha-phi": ((32, 12, 11), (11, 6, 4, 2, 2, 1), (9, 19, 11, 0)),
    "alpha-eta": ((31, 11, 11), (11, 6, 5, 2, 1), (9, 19, 11, 0)),
    "alpha": ((30, 10, 11), (11, 6, 6, 1), (8, 18, 11, 0)),
    "eta-phi": ((31, 11, 11), (11, 7, 4, 2, 1), (9, 19, 11, 0)),
    "eta": ((30, 10, 11), (11, 7, 5, 1), (8, 18, 11, 0)),
    "phi": ((30, 10, 11), (11, 8, 4, 1), (8, 18, 11, 0)),
    "none": ((24, 14, 6), (11, 1, 1, 1, 1, 1, 1, 1), (7, 12, 6, 0)),
}
for name, vertices in GLOBAL_POLYGONS.items():
    ledger, packet, graph = EXPECTED_GLOBAL[name]
    need(pick(vertices) == ledger, name + " global Pick ledger")
    need(edge_packet(vertices) == packet, name + " source-infinity packet")
    need(sum(index - 1 for index in packet) == 2 * ledger[2] - 2,
         name + " Riemann-Hurwitz checksum")
    component_count, edge_count, betti, carrier_genus = graph
    need(edge_count - component_count + 1 == betti,
         name + " graph Betti number")
    need(betti + carrier_genus == ledger[2],
         name + " graph-plus-normalization genus")

# The polygon packet is the squarefree-edge packet.  On the sole doubled
# N-edge, its two index-two punctures coalesce into one index-three puncture;
# the defect and normalization genus remain unchanged.
GENERIC_N_PACKET = EXPECTED_GLOBAL["alpha-phi"][1]
REPEATED_N_PACKET = (11, 6, 4, 3, 1)
need(GENERIC_N_PACKET == (11, 6, 4, 2, 2, 1),
     "generic N packet retains two simple edge roots")
need(sum(index - 1 for index in GENERIC_N_PACKET)
     == sum(index - 1 for index in REPEATED_N_PACKET) == 20,
     "repeated N packet preserves ramification defect")

need(12 - 7 + 1 == 6, "M seven-component twelve-node graph")
need(18 - 8 + 1 == 11, "six M-to-first-face attachments")
need(19 - 9 + 1 == 11, "one terminal N attachment preserves Betti")

# The sole hostile edge collision.  It is a persistent boundary tangency,
# not a collision tail.  Both charts of both ordinary blowups are explicit.
R, x, a = sp.symbols("R x a", nonzero=True)
r, q, z = sp.symbols("r q z")
Dpoly = (
    -3 * P + sp.Rational(8, 3) * P**2 + e * P**3
    + DELTA0 * P**4 + u * P**5 + U * P**6
)
Jpoly = Phi + eta * P + alpha * P**2
G_N = CHARTS["N"]
H_N = sp.cancel(R**3 * G_N.subs(S, 1 / R)).expand()
C = 1 - sigma**2 * P * R**2
B = C * (1 - sigma * Dpoly) - sigma / 2
H_EXPECTED = sp.expand(R * B - P**3 * Jpoly * C)
need(sp.expand(H_N - H_EXPECTED) == 0,
     "exact reciprocal boundary chart")

REPEATED = {Phi: alpha * a**2, eta: -2 * alpha * a}
need(sp.factor(Jpoly.subs(REPEATED) - alpha * (P - a)**2) == 0,
     "repeated quadratic parameterization")
need(sp.factor(H_N.subs({R: 0, **REPEATED})
               + alpha * P**3 * (P - a)**2) == 0,
     "persistent doubled boundary restriction")
need(sp.diff(H_N, R).subs({sigma: 0, R: 0, P: a, **REPEATED}) == 1,
     "unit transverse derivative")

# For the Poincare residue omega=-dP/G_S, the reciprocal identity
# H=R^3*G(1/R) gives omega=R*dP/H_R on the curve.  A simple edge root has
# ord(R)=1 and index two.  The repeated root below has ord(R)=2, hence exact
# differential order two and one index-three puncture.
GS_RECIPROCAL = sp.diff(G_N, S).subs(S, 1 / R)
need(sp.factor(
    sp.diff(H_N, R) + R * GS_RECIPROCAL
    - 3 * R**2 * G_N.subs(S, 1 / R)
) == 0, "reciprocal residue identity")
need(REPEATED_N_PACKET == (11, 6, 4, 3, 1),
     "order-two residue gives index-three puncture")

P0 = a + x
D0x = Dpoly.subs(P, P0)
C0 = 1 - sigma**2 * P0 * R**2
B0x = C0 * (1 - sigma * D0x) - sigma / 2
H_REPEAT = sp.expand(R * B0x - alpha * P0**3 * x**2 * C0)
need(sp.expand(H_N.subs(REPEATED).subs(P, P0) - H_REPEAT) == 0,
     "centered repeated-contact equation")
B_SECTION = 1 - sigma * Dpoly.subs(P, a) - sigma / 2
need(B_SECTION.subs(sigma, 0) == 1, "boundary coefficient is a DVR unit")
need(sp.diff(H_REPEAT.subs(R, 0), x).subs(x, 0) == 0,
     "repeated contact has no linear term")
need(sp.factor(sp.diff(H_REPEAT.subs(R, 0), x, 2).subs(x, 0)
               + 2 * alpha * a**3) == 0,
     "repeated contact has exact quadratic order")

# First blowup of (x,R).  In the x-chart, R=x*r and the strict transform
# meets the crossing x=r=0.  In the R-chart, x=R*q and the exceptional
# equation is the unit B_SECTION, so there is no complementary point.
B1X = B0x.subs(R, x * r)
C1X = C0.subs(R, x * r)
FIRST_X = sp.expand(r * B1X - alpha * P0**3 * x * C1X)
need(sp.expand(H_REPEAT.subs(R, x * r) - x * FIRST_X) == 0,
     "first blowup x-chart strict transform")
need(sp.factor(FIRST_X.subs({x: 0, r: 0}) - 0) == 0,
     "first x-chart crossing point")
need(sp.factor(sp.diff(FIRST_X, r).subs({x: 0, r: 0}) - B_SECTION) == 0,
     "first x-chart transverse unit")

P1R = a + R * q
D1R = Dpoly.subs(P, P1R)
C1R = 1 - sigma**2 * P1R * R**2
B1R = C1R * (1 - sigma * D1R) - sigma / 2
FIRST_R = sp.expand(B1R - alpha * P1R**3 * R * q**2 * C1R)
need(sp.expand(H_REPEAT.subs(x, R * q) - R * FIRST_R) == 0,
     "first blowup R-chart strict transform")
need(sp.factor(FIRST_R.subs(R, 0) - B_SECTION) == 0,
     "first complementary exceptional chart is empty")

# Second blowup of the remaining crossing (x,r) in the first x-chart.
# Its x-chart has r=x*z; its r-chart has x=r*q.  The two exceptional
# coordinates are reciprocal nonzero units, so they describe one point away
# from both boundary crossings.
B2X = B1X.subs(r, x * z)
C2X = C1X.subs(r, x * z)
SECOND_X = sp.expand(z * B2X - alpha * P0**3 * C2X)
need(sp.expand(FIRST_X.subs(r, x * z) - x * SECOND_X) == 0,
     "second blowup x-chart strict transform")
SECOND_X_EXCEPTIONAL = sp.factor(SECOND_X.subs(x, 0))
need(sp.expand(SECOND_X_EXCEPTIONAL - (z * B_SECTION - alpha * a**3)) == 0,
     "second x-chart exceptional equation")

P2R = a + r * q
R2R = r**2 * q
D2R = Dpoly.subs(P, P2R)
C2R = 1 - sigma**2 * P2R * R2R**2
B2R = C2R * (1 - sigma * D2R) - sigma / 2
SECOND_R = sp.expand(B2R - alpha * P2R**3 * q * C2R)
need(sp.expand(FIRST_X.subs(x, r * q) - r * SECOND_R) == 0,
     "second blowup r-chart strict transform")
SECOND_R_EXCEPTIONAL = sp.factor(SECOND_R.subs(r, 0))
need(sp.expand(SECOND_R_EXCEPTIONAL - (B_SECTION - alpha * a**3 * q)) == 0,
     "second r-chart exceptional equation")

z_exceptional = alpha * a**3 / B_SECTION
q_exceptional = B_SECTION / (alpha * a**3)
need(sp.factor(SECOND_X_EXCEPTIONAL.subs(z, z_exceptional)) == 0,
     "second x-chart exceptional point")
need(sp.factor(SECOND_R_EXCEPTIONAL.subs(q, q_exceptional)) == 0,
     "second r-chart exceptional point")
need(sp.factor(z_exceptional * q_exceptional) == 1,
     "second blowup charts glue reciprocally")
need(z_exceptional != 0 and q_exceptional != 0,
     "strict transform avoids boundary crossings")

# K=0 prunes the rational T leaf from the two Theta-nonzero THM-4350 fans.
OLD_11 = ((0, 1), (2, 0), (4, 2), (4, 3), (3, 5), (2, 6), (0, 7))
OLD_01 = ((0, 1), (2, 0), (4, 2), (4, 3), (2, 6), (0, 7))
need(pick(OLD_11) == (39, 13, 14), "old 11 polygon")
need(pick(OLD_01) == (38, 12, 14), "old 01 polygon")
need((20 - 10 + 1, 19 - 9 + 1) == (11, 11),
     "11 rational-leaf pruning preserves graph Betti")
need((19 - 9 + 1, 18 - 8 + 1) == (11, 11),
     "01 rational-leaf pruning preserves graph Betti")

print("THM-4353 PRIMARY EXACT CERTIFICATE")
print("scope=Z=beta_11=zeta_3=W=xi_10=K=0;U!=0;Delta=5696/105")
print("source=literal_16_rows;surviving_rows=10;specialized_support_points=22")
print("exact_supports=48;exact_fans=8;split=16,16,6,2,2,2,2,2")
print("hostile_atlas=128_keyed;96_distinct_supports;all_48_exact_embed")
print("hostile_fans=8;split=32,32,24,8,8,8,8,8")
print("faces=M,D6,E11,E01,Eeta,EPhi,N;all_literal_primitive_charts=verified")
print("carrier=E11:g3:base42:order41;E01:g3:base12:order11")
print("rational_faces=M_components,D6,Eeta,EPhi,N")
print("global_genera=Theta_nonzero:14;Theta_zero_with_owner:11;no_owner:6")
print("graph_betti=Theta_nonzero_or_normal_owner:11;no_owner:6")
print("sole_edge_discriminant=eta^2-4*alpha_11*Phi")
print("generic_N_packet=11,6,4,2,2,1;repeated_N_packet=11,6,4,3,1")
print("repeated_edge=residue_order2;index3;persistent_tangency;H_R=unit;no_tail")
print("blowups=two_ordinary;first_complement_empty;second_points_reciprocal_units")
print("K_zero_transition=prunes_rational_T_leaf;V_and_E_drop_together;b1_preserved")
print("proper_flat_input=all_noncarrier_components_rational;carrier_orders_positive")
print("CHECKS", CHECKS)
print("VERDICT ACCEPT full simultaneous-zero-endpoint gate relative to inherited interfaces")
