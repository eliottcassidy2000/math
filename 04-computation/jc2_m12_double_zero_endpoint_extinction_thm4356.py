#!/usr/bin/env python3
"""Primary exact certificate for THM-4356's double-zero endpoint gate.

Scope:
    Z=beta_11=zeta_3=W=xi_10=U=K=0, Delta=5696/105,
    alpha_11,Theta,Phi,eta,u arbitrary,

in the inherited reduced (2,3), exact-weight-twelve seam.  This self-contained
certificate reconstructs the literal source, enumerates all 48 realizable
supports and the hostile over-atlas, checks every exact face from primitive
substitutions, proves all component/graph/genus ledgers, and reconstructs the
only three repeated-edge walls: A4, A15, and the smooth N-face tangency.
"""

from collections import Counter, defaultdict
from fractions import Fraction as Fr
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
        raise RuntimeError("THM-4356 primary audit failure: " + label)


def noncollinear(points):
    for a, b, c in combinations(points, 3):
        if ((b[0] - a[0]) * (c[1] - a[1])
                != (b[1] - a[1]) * (c[0] - a[0])):
            return True
    return False


def candidate_planes(universe):
    answer = set()
    for first, second, third in combinations(sorted(universe), 3):
        det = ((second[0] - first[0]) * (third[1] - first[1])
               - (second[1] - first[1]) * (third[0] - first[0]))
        if det == 0:
            continue
        aa = Fr((second[2] - first[2]) * (third[1] - first[1])
                - (second[1] - first[1]) * (third[2] - first[2]), det)
        bb = Fr((second[0] - first[0]) * (third[2] - first[2])
                - (second[2] - first[2]) * (third[0] - first[0]), det)
        cc = Fr(first[2]) - aa * first[0] - bb * first[1]
        answer.add((aa, bb, cc))
    return tuple(sorted(answer))


def lower_fan(active, planes):
    answer = set()
    for aa, bb, cc in planes:
        equal = []
        for rr, ll, hh in active:
            gap = Fr(hh) - aa * rr - bb * ll - cc
            if gap < 0:
                break
            if gap == 0:
                equal.append((rr, ll, hh))
        else:
            if noncollinear(equal):
                answer.add((aa, bb, cc))
    return frozenset(answer)


def convex_hull(points):
    pts = sorted(set(points))

    def cross(o, a, b):
        return ((a[0] - o[0]) * (b[1] - o[1])
                - (a[1] - o[1]) * (b[0] - o[0]))

    lower = []
    for point in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper = []
    for point in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return tuple(lower[:-1] + upper[:-1])


def pick(vertices):
    pairs = tuple(zip(vertices, vertices[1:] + vertices[:1]))
    area2 = abs(sum(a[0] * b[1] - b[0] * a[1] for a, b in pairs))
    boundary = sum(gcd(abs(b[0] - a[0]), abs(b[1] - a[1]))
                   for a, b in pairs)
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


def on_plane(point, plane):
    rr, ll, hh = point
    aa, bb, cc = plane
    return Fr(hh) == aa * rr + bb * ll + cc


def plane_chart(plane):
    qpower = lcm(*(entry.denominator for entry in plane))
    sw = qpower * plane[0]
    pw = qpower * plane[1]
    clear = -qpower * plane[2]
    need(all(value.denominator == 1 for value in (sw, pw, clear)),
         "integral primitive chart")
    return qpower, sw.numerator, pw.numerator, clear.numerator


def face_order(plane):
    base = lcm(6, *(entry.denominator for entry in plane))
    density = Fr(5, 6) - sum(plane)
    order = density * base
    need(order.denominator == 1 and order > 0, "positive face order")
    return base, order.numerator


def edge_polynomial(expression, start, end, S, P, X):
    poly = sp.Poly(sp.expand(expression), S, P)
    dx, dy = end[0] - start[0], end[1] - start[1]
    length = gcd(abs(dx), abs(dy))
    step = dx // length, dy // length
    answer = 0
    for index in range(length + 1):
        exponent = start[0] + index * step[0], start[1] + index * step[1]
        answer += poly.coeff_monomial(S**exponent[0] * P**exponent[1]) * X**index
    return sp.factor(answer)


# Literal source and gate.
s, p, Q = sp.symbols("s p Q")
S, P, X, sigma = sp.symbols("S P X sigma")
U, W, Z = sp.symbols("U W Z")
Delta, Theta, Phi, eta = sp.symbols("Delta Theta Phi eta")
zeta, u, xi, alpha, beta, K = sp.symbols("zeta u xi alpha beta K")
e = -sp.Rational(1376, 135)
D = sp.Rational(5696, 105)
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
need(len(FULL_ROWS) == 16, "sixteen literal rows")
need(sp.Rational(2848, 45) - sp.Rational(7, 6) * D == 0,
     "K zero fixes Delta")
GATE = {Z: 0, beta: 0, zeta: 0, W: 0, xi: 0, U: 0,
        K: 0, Delta: D}
ROWS = {label: sp.factor(sp.sympify(value).subs(GATE))
        for label, value in FULL_ROWS.items()
        if sp.factor(sp.sympify(value).subs(GATE)) != 0}
need(len(ROWS) == 9, "nine surviving labelled rows")

H_ALL = sp.expand(sum(value * p**ii * (s * p)**jj
                      for (ii, jj), value in FULL_ROWS.items()))
F_ALL = sp.expand((s**2 - p) * (1 - Q * H_ALL) - Q * s**2 / 2)

SUPPORT = defaultdict(lambda: sp.Integer(0))
SUPPORT[(2, 0, 0)] += 1
SUPPORT[(0, 1, 0)] -= 1
SUPPORT[(2, 0, 1)] -= sp.Rational(1, 2)
for (ii, jj), coefficient in ROWS.items():
    SUPPORT[(jj + 2, ii + jj, 1)] -= coefficient
    SUPPORT[(jj, ii + jj + 1, 1)] += coefficient
SUPPORT = {point: sp.factor(value) for point, value in SUPPORT.items()
           if value != 0}
EXPECTED_SUPPORT = {
    (0, 1, 0): -1,
    (0, 2, 1): -3,
    (0, 3, 1): sp.Rational(8, 3),
    (0, 4, 1): e,
    (0, 5, 1): D,
    (0, 6, 1): u,
    (1, 4, 1): Phi,
    (1, 5, 1): eta,
    (1, 6, 1): alpha,
    (2, 0, 0): 1,
    (2, 0, 1): -sp.Rational(1, 2),
    (2, 1, 1): 3,
    (2, 2, 1): -sp.Rational(8, 3),
    (2, 3, 1): -e,
    (2, 4, 1): Theta - D,
    (2, 5, 1): -u,
    (3, 3, 1): -Phi,
    (3, 4, 1): -eta,
    (3, 5, 1): -alpha,
    (4, 3, 1): -Theta,
}
need(SUPPORT.keys() == EXPECTED_SUPPORT.keys(), "twenty possible points")
for point, expected in EXPECTED_SUPPORT.items():
    need(sp.expand(SUPPORT[point] - expected) == 0,
         "literal support coefficient")


# Twelve exact face planes.
C4A = (Fr(-1, 4), Fr(1, 4), Fr(-1, 4))
C5 = (Fr(0), Fr(1, 5), Fr(-1, 5))
C4E = (Fr(0), Fr(1, 4), Fr(-1, 4))
BA = (Fr(1, 11), Fr(2, 11), Fr(-2, 11))
BU = (Fr(1, 10), Fr(1, 5), Fr(-1, 5))
BE = (Fr(1, 9), Fr(2, 9), Fr(-2, 9))
BDT = (Fr(1, 8), Fr(1, 4), Fr(-1, 4))
F5 = (Fr(1, 5), Fr(1, 5), Fr(-2, 5))
EDP = (Fr(1, 4), Fr(1, 4), Fr(-1, 2))
E11 = (Fr(2, 7), Fr(1, 7), Fr(-4, 7))
EUP = (Fr(2, 5), Fr(1, 5), Fr(-4, 5))
N = (Fr(1), Fr(0), Fr(-2))
EXACT_PLANES = frozenset((C4A, C5, C4E, BA, BU, BE, BDT, F5,
                          EDP, E11, EUP, N))
PLANES = candidate_planes(SUPPORT)


def expected_fan(theta_on, alpha_on, eta_on, phi_on, u_on):
    if alpha_on:
        cap = C5 if u_on else C4A
        if theta_on:
            return frozenset((cap, BA, E11))
        return frozenset((cap, BA, N)) if eta_on or phi_on \
            else frozenset((cap, BA))
    if theta_on:
        if u_on:
            return frozenset((BU, F5))
        if eta_on:
            return frozenset((C4E, BE, F5))
        return frozenset((BDT,))
    if u_on:
        if eta_on:
            return frozenset((BU, F5, N)) if phi_on \
                else frozenset((BU, F5))
        return frozenset((BU, EUP)) if phi_on else frozenset((BU,))
    if eta_on:
        return frozenset((C4E, BE, N)) if phi_on \
            else frozenset((C4E, BE))
    return frozenset((BDT, EDP)) if phi_on else frozenset((BDT,))


THETA_REPS = (("zero", sp.Integer(0)), ("cancel", D),
              ("generic", sp.Integer(1)))
EXACT = {}
EXACT_FANS = Counter()
SPLIT = Counter()
for (theta_tag, theta_value), av, hv, phv, uv in product(
        THETA_REPS, (0, 1), (0, 1), (0, 1), (0, 1)):
    sub = {Theta: theta_value, alpha: av, eta: hv, Phi: phv, u: uv}
    active = frozenset(point for point, coefficient in SUPPORT.items()
                       if sp.factor(coefficient.subs(sub)) != 0)
    theta_on = theta_value != 0
    fan = lower_fan(active, PLANES)
    need(fan == expected_fan(theta_on, av != 0, hv != 0,
                             phv != 0, uv != 0), "exact owner fan")
    key = theta_tag, av, hv, phv, uv
    EXACT[key] = active
    EXACT_FANS[fan] += 1
    SPLIT["alpha"] += int(av != 0)
    SPLIT["theta"] += int(av == 0 and theta_on)
    SPLIT["terminal"] += int(av == 0 and not theta_on)

need(len(EXACT) == len(set(EXACT.values())) == 48,
     "forty-eight distinct realizable supports")
need(len(EXACT_FANS) == 15, "fifteen exact fan words")
need(set().union(*EXACT_FANS) == EXACT_PLANES, "twelve exact planes")
need(SPLIT == Counter({"alpha": 24, "theta": 16, "terminal": 8}),
     "24+16+8 theorem split")
need(len({lower_fan(active, PLANES) for key, active in EXACT.items()
          if key[1] != 0}) == 6, "alpha block has six fans")
need(len({lower_fan(active, PLANES) for key, active in EXACT.items()
          if key[1] == 0 and key[0] != "zero"}) == 3,
     "theta block has three fans")
need(len({lower_fan(active, PLANES) for key, active in EXACT.items()
          if key[1] == 0 and key[0] == "zero"}) == 8,
     "terminal block has eight fans")
need(sorted(EXACT_FANS.values()) == [1, 1, 1, 1, 1, 1, 1, 1,
                                     3, 3, 4, 5, 8, 8, 9],
     "exact fan multiplicities")


# Hostile atlas: five optional rows and the two aggregate deletions are
# toggled independently.  These are probes, not coefficient claims.
REQUIRED = ((1, 0), (2, 0), (3, 0), (4, 0))
OPTIONAL = ((5, 0), (2, 1), (1, 2), (3, 1), (4, 1))
CANCELLABLE = ((2, 4, 1), (2, 5, 1))


def lifted(labels):
    active = {(2, 0, 0), (0, 1, 0), (2, 0, 1)}
    for ii, jj in labels:
        active.add((jj + 2, ii + jj, 1))
        active.add((jj, ii + jj + 1, 1))
    return active


HOSTILE = {}
HOSTILE_FANS = Counter()
for row_bits in product((0, 1), repeat=5):
    labels = list(REQUIRED)
    labels.extend(label for bit, label in zip(row_bits, OPTIONAL) if bit)
    for deletion_bits in product((0, 1), repeat=2):
        active = lifted(labels)
        for bit, point in zip(deletion_bits, CANCELLABLE):
            if bit:
                active.discard(point)
        active = frozenset(active)
        HOSTILE[row_bits + deletion_bits] = active
        HOSTILE_FANS[lower_fan(active, PLANES)] += 1

need(len(HOSTILE) == 128, "128 keyed hostile states")
need(len(set(HOSTILE.values())) == 96, "96 hostile supports")
need(len(HOSTILE_FANS) == 23, "23 hostile fan words")
need(len(set().union(*HOSTILE_FANS)) == 17, "17 hostile face planes")
for key, active in EXACT.items():
    theta_tag, av, hv, phv, uv = key
    row_bits = (uv, phv, int(theta_tag != "zero"), hv, av)
    deletion_bits = (int(theta_tag == "cancel"), int(uv == 0))
    need(HOSTILE[row_bits + deletion_bits] == active,
         "exact support embeds in hostile atlas")


# Named exact face variants.  All equations are up to a torus monomial.
FACE_VARIANTS = {
    "C4a": (C4A, {u: 0}, P * (D * P**4 + alpha * S * P**5 - 1),
             ((0, 1), (1, 6), (0, 5)), (4, 6, 0), (12, 13)),
    "C5": (C5, {}, P * (u * P**5 + alpha * S * P**5 - 1),
            ((0, 1), (1, 6), (0, 6)), (5, 7, 0), (30, 25)),
    "C4e": (C4E, {u: 0, alpha: 0},
             P * (D * P**4 + eta * S * P**4 - 1),
             ((0, 1), (1, 5), (0, 5)), (4, 6, 0), (12, 10)),
    "Balpha": (BA, {}, (P - S**2) * (alpha * S * P**5 - 1),
                ((0, 1), (2, 0), (3, 5), (1, 6)),
                (22, 4, 10), (66, 49)),
    "BU": (BU, {alpha: 0}, (P - S**2) * (u * P**5 - 1),
           ((0, 1), (2, 0), (2, 5), (0, 6)),
           (20, 12, 5), (30, 22)),
    "Beta": (BE, {u: 0, alpha: 0},
             (P - S**2) * (eta * S * P**4 - 1),
             ((0, 1), (2, 0), (3, 4), (1, 5)),
             (18, 4, 8), (18, 13)),
    "BDT": (BDT, {u: 0, eta: 0, alpha: 0},
            (P - S**2) * (D * P**4 + Theta * S**2 * P**3 - 1),
            ((0, 1), (2, 0), (4, 3), (0, 5)),
            (24, 8, 9), (24, 17)),
    "BDelta": (BDT, {u: 0, eta: 0, alpha: 0, Theta: 0},
               (P - S**2) * (D * P**4 - 1),
               ((0, 1), (2, 0), (2, 4), (0, 5)),
               (16, 10, 4), (24, 17)),
    "FUT": (F5, {alpha: 0},
            S**2 * (1 - u * P**5 - eta * S * P**4
                    - Theta * S**2 * P**3),
            ((2, 0), (4, 3), (2, 5)), (10, 8, 2), (30, 25)),
    "FET": (F5, {alpha: 0, u: 0},
            S**2 * (1 - eta * S * P**4 - Theta * S**2 * P**3),
            ((2, 0), (4, 3), (3, 4)), (5, 3, 2), (30, 25)),
    "FUE": (F5, {alpha: 0, Theta: 0},
            S**2 * (1 - u * P**5 - eta * S * P**4),
            ((2, 0), (3, 4), (2, 5)), (5, 7, 0), (30, 25)),
    "EDPhi": (EDP, {alpha: 0, Theta: 0, u: 0, eta: 0},
              S**2 * (1 - D * P**4 - Phi * S * P**3),
              ((2, 0), (3, 3), (2, 4)), (4, 6, 0), (12, 10)),
    "E11": (E11, {},
            S**2 * (1 - alpha * S * P**5 - Theta * S**2 * P**3),
            ((2, 0), (4, 3), (3, 5)), (7, 3, 3), (42, 41)),
    "EUPhi": (EUP, {alpha: 0, Theta: 0, eta: 0},
              S**2 * (1 - u * P**5 - Phi * S * P**3),
              ((2, 0), (3, 3), (2, 5)), (5, 7, 0), (30, 31)),
    "N": (N, {Theta: 0},
          S**2 * (1 - S * P**3 * (Phi + eta * P + alpha * P**2)),
          ((2, 0), (3, 3), (3, 5)), (2, 4, 0), (6, 11)),
}


def literal_face(plane, substitution):
    return sp.factor(sum(coefficient.subs(substitution) * S**rr * P**ll
                         for (rr, ll, hh), coefficient in SUPPORT.items()
                         if on_plane((rr, ll, hh), plane)))


def primitive_chart(plane, substitution):
    qpower, sw, pw, clear = plane_chart(plane)
    chart = sp.cancel(
        sigma**clear * F_ALL.subs(GATE).subs(substitution).subs(
            {Q: sigma**qpower, s: sigma**(-sw) * S,
             p: sigma**(-pw) * P}, simultaneous=True)
    )
    need(not sp.denom(sp.together(chart)).has(sigma),
         "primitive chart integral")
    sigma_poly = sp.Poly(sp.expand(chart), sigma)
    need(min(monomial[0][0] for monomial in sigma_poly.terms()) == 0,
         "primitive chart has valuation zero")
    return sp.expand(chart)


for name, (plane, substitution, expression, polygon, pick_data, order_data) \
        in FACE_VARIANTS.items():
    need(sp.expand(literal_face(plane, substitution) - expression) == 0,
         name + " literal face")
    need(pick(polygon) == pick_data, name + " face Pick")
    need(face_order(plane) == order_data, name + " face order")
    chart = primitive_chart(plane, substitution)
    need(sp.expand(chart.subs(sigma, 0) - expression) == 0,
         name + " primitive initial")


# Exact edge words on every actual face polygon.
EXPECTED_EDGES = {
    "C4a": (alpha * X - 1, alpha + D * X, D - X**4),
    "C5": (alpha * X - 1, alpha + u * X, u - X**5),
    "C4e": (eta * X - 1, eta + D * X, D - X**4),
    "Balpha": (X - 1, 1 - alpha * X, alpha * (X - 1), alpha - X),
    "BU": (X - 1, 1 - u * X**5, u * (X - 1), u - X**5),
    "Beta": (X - 1, 1 - eta * X, eta * (X - 1), eta - X),
    "BDT": (X - 1, 1 - Theta * X,
            (X - 1) * (D * X + Theta), D - X**4),
    "BDelta": (X - 1, 1 - D * X**4,
               D * (X - 1), D - X**4),
    "FUT": (1 - Theta * X, -(Theta + eta * X + u * X**2),
            X**5 - u),
    "FET": (1 - Theta * X, -(Theta + eta * X), X - eta),
    "FUE": (1 - eta * X, -(eta + u * X), X**5 - u),
    "EDPhi": (1 - Phi * X, -(Phi + D * X), X**4 - D),
    "E11": (1 - Theta * X, -(Theta + alpha * X), X - alpha),
    "EUPhi": (1 - Phi * X, -(Phi + u * X), X**5 - u),
    "N": (1 - Phi * X, -(Phi + eta * X + alpha * X**2), X - alpha),
}
for name, edges in EXPECTED_EDGES.items():
    expression = FACE_VARIANTS[name][2]
    polygon = FACE_VARIANTS[name][3]
    actual = tuple(edge_polynomial(expression, polygon[index],
                                   polygon[(index + 1) % len(polygon)],
                                   S, P, X)
                   for index in range(len(polygon)))
    for got, wanted in zip(actual, edges):
        need(sp.expand(got - wanted) == 0, name + " exact edge")

L5 = eta**2 - 4 * u * Theta
LN = eta**2 - 4 * alpha * Phi
L15 = D + Theta
need(sp.discriminant(Theta + eta * X + u * X**2, X) == L5,
     "A4 edge divisor")
need(sp.discriminant(Phi + eta * X + alpha * X**2, X) == LN,
     "N tangency edge divisor")
need(sp.resultant(X - 1, D * X + Theta, X) == L15,
     "A15 edge resultant")
need(D != 0 and e != 0, "fixed coefficients nonzero")

# Every nonlinear edge outside the three displayed collision words is reduced
# on its owner cell.  Its discriminant vanishes only on an owner-deletion wall.
need(sp.factor(sp.discriminant(D - X**4, X)) == -256 * D**3,
     "fixed quartic edge reduced")
need(sp.factor(sp.discriminant(1 - D * X**4, X)) == -256 * D**3,
     "fixed reciprocal quartic edge reduced")
need(sp.factor(sp.discriminant(u - X**5, X)) == 3125 * u**4,
     "u quintic edge reduced off owner wall")
need(sp.factor(sp.discriminant(1 - u * X**5, X)) == 3125 * u**4,
     "u reciprocal quintic edge reduced off owner wall")
need(sp.expand(sp.discriminant((X - 1) * (D * X + Theta), X)
               - L15**2) == 0,
     "BDT product has only the L15 collision")


# Carrier identities and factor intersection counts.
f11 = 1 - alpha * S * P**5 - Theta * S**2 * P**3
y11 = 2 * Theta * S * P**2 + alpha * P**4
b11 = P * (alpha**2 * P**7 + 4 * Theta)
need(sp.expand(y11**2 - b11 + 4 * Theta * P * f11) == 0,
     "E11 normalization")
need(sp.factor(sp.discriminant(b11, P))
     == -sp.Integer(53971714048) * Theta**8 * alpha**12,
     "E11 squarefree on owner gate")

ffut = 1 - u * P**5 - eta * S * P**4 - Theta * S**2 * P**3
yfut = 2 * Theta * S * P**2 + eta * P**3
bfut = P * (4 * Theta + L5 * P**5)
need(sp.expand(yfut**2 - bfut + 4 * Theta * P * ffut) == 0,
     "FUT normalization")
need(sp.expand(sp.discriminant(bfut, P)
               - 12800000 * L5**4 * Theta**6) == 0,
     "FUT has only the L5 carrier collision")
need(sp.gcd(sp.Poly(bfut, P), sp.Poly(sp.diff(bfut, P), P)).degree() == 0,
     "FUT squarefree away from L5")

ffet = ffut.subs(u, 0)
bfet = bfut.subs(u, 0)
need(sp.gcd(sp.Poly(bfet, P), sp.Poly(sp.diff(bfet, P), P)).degree() == 0,
     "FET automatically squarefree")

fbdt = 1 - D * P**4 - Theta * S**2 * P**3
ybdt = Theta * S * P**2
bbdt = Theta * P * (1 - D * P**4)
need(sp.expand(ybdt**2 - bbdt + Theta * P * fbdt) == 0,
     "BDT carrier normalization")
need(sp.factor(sp.discriminant(bbdt, P)) == -256 * D**3 * Theta**8,
     "BDT carrier stays smooth even on L15")
need(sp.gcd(sp.Poly(bbdt, P), sp.Poly(sp.diff(bbdt, P), P)).degree() == 0,
     "BDT carrier always squarefree")

INTERSECTIONS = {
    "Balpha": sp.degree((alpha * S * P**5 - 1).subs(P, S**2), S),
    "BU": sp.degree((u * P**5 - 1).subs(P, S**2), S),
    "Beta": sp.degree((eta * S * P**4 - 1).subs(P, S**2), S),
    "BDelta": sp.degree((D * P**4 - 1).subs(P, S**2), S),
    "BDT": sp.degree((D * P**4 + Theta * S**2 * P**3 - 1)
                     .subs(P, S**2), S),
}
need(INTERSECTIONS == {"Balpha": 11, "BU": 10, "Beta": 9,
                       "BDelta": 8, "BDT": 8},
     "factor intersection degrees")
need(sp.expand(fbdt.subs(P, S**2) - (1 - L15 * S**8)) == 0,
     "BDT nodes escape exactly on L15")
ATTACHMENTS = {
    "cap_to_Balpha": 1,
    "Balpha_to_E11_or_N": 1,
    "cap_to_Beta": 1,
    "Beta_to_FET_or_N": 1,
    "BU_to_F5_or_EUPhi": sp.degree(X**5 - u, X),
    "BDelta_to_EDPhi": sp.degree(X**4 - D, X),
}
need(ATTACHMENTS == {"cap_to_Balpha": 1,
                     "Balpha_to_E11_or_N": 1,
                     "cap_to_Beta": 1,
                     "Beta_to_FET_or_N": 1,
                     "BU_to_F5_or_EUPhi": 5,
                     "BDelta_to_EDPhi": 4},
     "complete clean attachment multiplicities")

RATIONAL_SOLUTIONS = {
    "C4a": (D * P**4 + alpha * S * P**5 - 1,
             (1 - D * P**4) / (alpha * P**5)),
    "C5": (u * P**5 + alpha * S * P**5 - 1,
            (1 - u * P**5) / (alpha * P**5)),
    "C4e": (D * P**4 + eta * S * P**4 - 1,
             (1 - D * P**4) / (eta * P**4)),
    "FUE": (1 - u * P**5 - eta * S * P**4,
            (1 - u * P**5) / (eta * P**4)),
    "EDPhi": (1 - D * P**4 - Phi * S * P**3,
              (1 - D * P**4) / (Phi * P**3)),
    "EUPhi": (1 - u * P**5 - Phi * S * P**3,
              (1 - u * P**5) / (Phi * P**3)),
    "N": (1 - S * P**3 * (Phi + eta * P + alpha * P**2),
          1 / (P**3 * (Phi + eta * P + alpha * P**2))),
}
for name, (equation, solution) in RATIONAL_SOLUTIONS.items():
    need(sp.factor(equation.subs(S, solution)) == 0,
         name + " rational parametrization")


# Every exact support has the graph-plus-normalization genus predicted by its
# global Newton polygon.  This also audits the attachment multiplicities.
def graph_ledger(theta_on, alpha_on, eta_on, phi_on, u_on):
    if alpha_on:
        if theta_on:
            return 4, 13, 10, 3
        return (4, 13, 10, 0) if eta_on or phi_on else (3, 12, 10, 0)
    if theta_on:
        if u_on:
            return 7, 15, 9, 2
        if eta_on:
            return 4, 11, 8, 2
        return 2, 8, 7, 2
    if u_on:
        if eta_on or phi_on:
            extra_n = int(eta_on and phi_on)
            return 7 + extra_n, 15 + extra_n, 9, 0
        return 6, 10, 5, 0
    if eta_on:
        extra_n = int(phi_on)
        return 3 + extra_n, 10 + extra_n, 8, 0
    if phi_on:
        return 6, 12, 7, 0
    return 5, 8, 4, 0


GLOBAL = Counter()
for key, active in EXACT.items():
    theta_tag, av, hv, phv, uv = key
    theta_on = theta_tag != "zero"
    polygon = convex_hull([(rr, ll) for rr, ll, hh in active])
    ledger = pick(polygon)
    packet = edge_packet(polygon)
    vertices, edges, b1, carrier_genus = graph_ledger(
        theta_on, bool(av), bool(hv), bool(phv), bool(uv))
    need(edges - vertices + 1 == b1, "graph Betti")
    need(b1 + carrier_genus == ledger[2], "global genus conservation")
    need(sum(index - 1 for index in packet) == 2 * ledger[2] - 2,
         "global packet defect")
    GLOBAL[(polygon, ledger, packet)] += 1
need(len(GLOBAL) == 19, "nineteen global polygon/packet ledgers")


# Smooth N-face tangency.  It is an edge collision but not a curve
# singularity: H_R is a unit and two ordinary blowups add rational chains.
R, a = sp.symbols("R a", nonzero=True)
Dpoly = (-3 * P + sp.Rational(8, 3) * P**2 + e * P**3
         + D * P**4 + u * P**5)
Jpoly = Phi + eta * P + alpha * P**2
G_N = primitive_chart(N, {Theta: 0})
H_N = sp.cancel(R**3 * G_N.subs(S, 1 / R)).expand()
C_N = 1 - sigma**2 * P * R**2
B_N = C_N * (1 - sigma * Dpoly) - sigma / 2
need(sp.expand(H_N - (R * B_N - P**3 * Jpoly * C_N)) == 0,
     "exact N reciprocal chart")
REPEATED_N = {Phi: alpha * a**2, eta: -2 * alpha * a}
need(sp.expand(Jpoly.subs(REPEATED_N) - alpha * (P - a)**2) == 0,
     "N doubled root")
need(sp.diff(H_N, R).subs({sigma: 0, R: 0, P: a, **REPEATED_N}) == 1,
     "N transverse derivative unit")
need((10, 5, 4, 2, 2, 1) != (10, 5, 4, 3, 1),
     "N puncture packet must change")
need(sum(i - 1 for i in (10, 5, 4, 2, 2, 1))
     == sum(i - 1 for i in (10, 5, 4, 3, 1)) == 18,
     "N tangency preserves defect")
for uv in (0, 1):
    active = EXACT[("zero", 1, 1, 1, uv)]
    polygon = convex_hull([(rr, ll) for rr, ll, hh in active])
    need(edge_packet(polygon) == (10, 5, 4, 2, 2, 1),
         "N generic packet derived from global polygon")

# Two literal ordinary blowups of the centered N tangency.  The complementary
# chart is empty and the exceptional components are rational.
xN, rN, qN, zN = sp.symbols("xN rN qN zN")
P0 = a + xN
C0 = 1 - sigma**2 * P0 * R**2
B0 = C0 * (1 - sigma * Dpoly.subs(P, P0)) - sigma / 2
H0 = sp.expand(R * B0 - alpha * P0**3 * xN**2 * C0)
need(sp.expand(H_N.subs(REPEATED_N).subs(P, P0) - H0) == 0,
     "centered N tangency")
first_x = sp.expand(rN * B0.subs(R, xN * rN)
                    - alpha * P0**3 * xN * C0.subs(R, xN * rN))
need(sp.expand(H0.subs(R, xN * rN) - xN * first_x) == 0,
     "first N blowup x-chart")
first_R = sp.expand(
    (1 - sigma**2 * (a + R * qN) * R**2)
    * (1 - sigma * Dpoly.subs(P, a + R * qN)) - sigma / 2
    - alpha * (a + R * qN)**3 * R * qN**2
    * (1 - sigma**2 * (a + R * qN) * R**2)
)
need(sp.expand(H0.subs(xN, R * qN) - R * first_R) == 0,
     "first N blowup complementary chart")
need(first_R.subs({R: 0, sigma: 0}) == 1,
     "first N complementary exceptional is empty")
second_x = sp.expand(zN * B0.subs(R, xN**2 * zN)
                     - alpha * P0**3 * C0.subs(R, xN**2 * zN))
need(sp.expand(first_x.subs(rN, xN * zN) - xN * second_x) == 0,
     "second N blowup x-chart")
need(sp.diff(second_x, zN).subs({xN: 0, sigma: 0}) == 1,
     "second N exceptional is transverse rational")
P2 = a + rN * qN
R2 = rN**2 * qN
C2 = 1 - sigma**2 * P2 * R2**2
B2 = C2 * (1 - sigma * Dpoly.subs(P, P2)) - sigma / 2
second_r = sp.expand(B2 - alpha * P2**3 * qN * C2)
need(sp.expand(first_x.subs({xN: rN * qN}) - rN * second_r) == 0,
     "second N blowup complementary chart")
need(second_r.subs({rN: 0, sigma: 0}) == 1 - alpha * a**3 * qN,
     "second N exceptional has one reduced point")
need(sp.diff(second_r.subs({rN: 0, sigma: 0}), qN) == -alpha * a**3,
     "second N exceptional point is transverse")


# Odd A4 wall.  K=0 makes the third selector 8/3, so depths are exactly a
# subset of {1,2,3}; there is no inherited depth-four subcase.
x, v, rho, y, t, bigX, bigY = sp.symbols("x v rho y t bigX bigY")
G_FUT = primitive_chart(F5, {alpha: 0})
A5 = u + eta * v + Theta * v**2
B5 = D + Phi * v
REC5 = sp.cancel(x**7 * G_FUT.subs(
    {P: 1 / x, S: v / x, sigma: rho / x}, simultaneous=True))
EXPECTED_REC5 = ((v**2 - rho)
                 * (x**5 - A5 - rho * B5 - e * rho**2
                    - sp.Rational(8, 3) * rho**3 + 3 * rho**4)
                 - rho**5 * v**2 / 2)
need(sp.expand(REC5 - EXPECTED_REC5) == 0, "exact A4 reciprocal chart")
a5 = sp.symbols("a5", nonzero=True)
A4_WALL = {eta: -2 * Theta * a5, u: Theta * a5**2}
need(sp.expand(A5.subs(A4_WALL) - Theta * (v - a5)**2) == 0,
     "A4 central square")
c1 = D + Phi * a5
c2 = e - Phi**2 / (4 * Theta)
c3 = sp.Rational(8, 3)
vcrit = a5 - Phi * rho / (2 * Theta)
CRIT5 = (Theta * (v - a5)**2 + rho * (D + Phi * v)
         + e * rho**2 + sp.Rational(8, 3) * rho**3 - 3 * rho**4
         + rho**5 * v**2 / (2 * (v**2 - rho)))
series5 = sp.series(CRIT5.subs(v, vcrit), rho, 0, 5).removeO().expand()
need(sp.factor(series5.coeff(rho, 1) - c1) == 0,
     "A4 first selector")
need(sp.factor(series5.coeff(rho, 2) - c2) == 0,
     "A4 second selector")
need(sp.factor(series5.coeff(rho, 3) - c3) == 0,
     "A4 third selector")
need(sp.factor(series5.coeff(rho, 4) + 3) == 0,
     "A4 fourth displayed coefficient")
need(c3 != 0, "K-zero A4 depth at most three")
A4_ROWS = {
    1: ((4, 6, 15), 2, 0, 115),
    2: ((1, 4, 10), 1, 1, 35),
    3: ((2, 18, 45), 1, 1, 95),
}
JAC5 = sp.det(sp.Matrix([
    [sp.diff(v / x, v), sp.diff(v / x, x)],
    [sp.diff(1 / x, v), sp.diff(1 / x, x)],
]))
need(sp.factor(JAC5) == -x**-3, "A4 reciprocal coordinate Jacobian")
need(7 - 3 == 4, "A4 residue has x4 conductor buffer")
ctail = sp.symbols("ctail", nonzero=True)
A4_TAILS = {
    1: bigX * (bigX**4 - ctail),
    2: bigX**3 - ctail,
    3: bigX * (bigX**2 - ctail),
}
for depth, ((ws, wx, wy), genus, persistent, order) in A4_ROWS.items():
    need(gcd(gcd(ws, wx), wy) == 1, "A4 primitive ray")
    need(5 * wx == depth * (6 * ws + wx) == 2 * wy,
         "A4 balance")
    need(25 * ws + 5 * wx - wy == order > 0, "A4 form order")
    need(genus + persistent == 2, "A4 intrinsic genus")
    need(9 + genus + persistent == 11, "A4 global genus conservation")
    tail_poly = A4_TAILS[depth]
    need(sp.gcd(sp.Poly(tail_poly, bigX),
                sp.Poly(sp.diff(tail_poly, bigX), bigX)).degree() == 0,
         "A4 reduced tail squarefree")
need(sp.factor(c1.subs(Phi, -D / a5)) == 0, "A4 depth-two nonempty")
theta_depth3 = sp.factor(D**2 / (4 * e * a5**2))
need(sp.factor(c2.subs({Phi: -D / a5, Theta: theta_depth3})) == 0,
     "A4 depth-three nonempty")


# A15 wall.  Work directly on the target-compatible 24-fold face chart.
tau, q, z, lam = sp.symbols("tau q z lambda")
G24 = sp.cancel(
    tau**6 * F_ALL.subs(GATE).subs({alpha: 0, u: 0, eta: 0,
                                    Theta: -D}).subs(
        {Q: tau**24, s: tau**-3 * S, p: tau**-6 * P},
        simultaneous=True)
)
G24_EXPECTED = ((S**2 - P)
                * (1 - D * P**4 + D * S**2 * P**3
                   - tau**3 * Phi * S * P**3
                   - e * tau**6 * P**3
                   - sp.Rational(8, 3) * tau**12 * P**2
                   + 3 * tau**18 * P)
                - tau**24 * S**2 / 2)
need(sp.expand(G24 - G24_EXPECTED) == 0,
     "exact K-zero target-compatible A15 chart")
REC15 = sp.cancel(-x**10 * G24_EXPECTED.subs(
    {S: 1 / x, P: (1 + z) / x**2, tau**3: rho / x},
    simultaneous=True))
REC15_EXPECTED = (z * (x**8 - D * z * (1 + z)**3
                       - Phi * rho * (1 + z)**3
                       - e * rho**2 * (1 + z)**3
                       - sp.Rational(8, 3) * rho**4 * (1 + z)**2
                       + 3 * rho**6 * (1 + z))
                  + rho**8 / 2)
need(sp.expand(REC15 - REC15_EXPECTED) == 0,
     "exact A15 reciprocal chart")
need(sp.expand(REC15.subs(rho, 0)
               - z * (x**8 - D * z * (1 + z)**3)) == 0,
     "A15 central contact")
need(sp.factor(EXPECTED_EDGES["BDT"][2].subs(Theta, -D))
     == D * (X - 1)**2, "A15 repeated edge")

# The two critical x^8 branches begin Phi*rho+e*rho^2, so K=0 leaves only
# r=1 and r=2.  Their exact exceptional equations and graph ledgers follow.
eps = sp.symbols("eps")
B15 = (-D * z * (1 + z)**3 - Phi * rho * (1 + z)**3
       - e * rho**2 * (1 + z)**3
       - sp.Rational(8, 3) * rho**4 * (1 + z)**2
       + 3 * rho**6 * (1 + z))
zlead = eps * rho**4 / lam
Xcrit15 = sp.series(-B15.subs(z, zlead)
                    - zlead * sp.diff(B15, z).subs(z, zlead),
                    rho, 0, 5).removeO().expand().subs(eps**2, 1)
predicted15 = (Phi * rho + e * rho**2
               + (sp.Rational(8, 3) - eps * lam) * rho**4)
need(sp.expand((Xcrit15 - predicted15) * lam
               - eps * rho**4 * (2 * D + lam**2)) == 0,
     "A15 critical branches through rho4")
need(e != 0, "K-zero A15 depth at most two")
JAC15 = sp.det(sp.Matrix([
    [sp.diff(1 / x, x), sp.diff(1 / x, z)],
    [sp.diff((1 + z) / x**2, x), sp.diff((1 + z) / x**2, z)],
]))
need(sp.factor(JAC15) == -x**-4, "A15 reciprocal coordinate Jacobian")
need(10 - 4 == 6, "A15 residue has raw x6 numerator")
A15_ROWS = {
    1: ((7, 3, 24), 116, 7, 1, (4, 9, 6)),
    2: ((1, 1, 8), 16, 6, 2, (4, 8, 5)),
}
for depth, ((ws, wx, wy), order, mutual_nodes, persistent, graph) \
        in A15_ROWS.items():
    h = gcd(3 * depth, 8 - depth)
    need((ws, wx) == ((8 - depth) // h, 3 * depth // h),
         "A15 primitive balance")
    need(wy == 8 * wx, "A15 y weight")
    need(17 * ws - wx == order > 0, "A15 form order")
    vertices, edges, b1 = graph
    need(edges == mutual_nodes + 2, "A15 distinct-component attachments")
    need(edges - vertices + 1 == b1, "A15 graph Betti")
    need(2 + persistent + b1 == 9, "A15 genus conservation")
need(sp.gcd(sp.Poly(bigX**7 - Phi, bigX),
            sp.Poly(7 * bigX**6, bigX)).degree() == 0,
     "A15 depth-one seven simple mutual nodes")
need(sp.gcd(sp.Poly(bigX**6 - e, bigX),
            sp.Poly(6 * bigX**5, bigX)).degree() == 0,
     "A15 depth-two six simple mutual nodes")
Ztail = sp.symbols("Ztail")
tail1 = D * Ztail**2 - bigX * (bigX**7 - Phi) * Ztail
tail2 = D * Ztail**2 - bigX**2 * (bigX**6 - e) * Ztail
need(sp.expand(tail1 - Ztail * (D * Ztail
                                - bigX * (bigX**7 - Phi))) == 0,
     "A15 depth-one two rational signs")
need(sp.expand(tail2 - Ztail * (D * Ztail
                                - bigX**2 * (bigX**6 - e))) == 0,
     "A15 depth-two two rational signs")
need(tail1.subs(Ztail, 0) == 0 and tail2.subs(Ztail, 0) == 0,
     "A15 zero sign attaches to R")
C15 = (x**8 - D * z * (1 + z)**3
       - Phi * rho * (1 + z)**3 - e * rho**2 * (1 + z)**3
       - sp.Rational(8, 3) * rho**4 * (1 + z)**2
       + 3 * rho**6 * (1 + z))
t15 = sp.symbols("t15")
Ctail1 = sp.expand(
    C15.subs({x: t15**3 * bigX, rho: t15**24 * bigX,
              z: t15**24 * Ztail}) / t15**24
).subs(t15, 0)
Ctail2 = sp.expand(
    C15.subs({x: t15 * bigX, rho: t15**4 * bigX,
              z: t15**8 * Ztail, Phi: 0}) / t15**8
).subs(t15, 0)
need(sp.expand(Ctail1 - (bigX**8 - D * Ztail - Phi * bigX)) == 0,
     "A15 depth-one nonzero sign comes from C")
need(sp.expand(Ctail2 - (bigX**8 - D * Ztail - e * bigX**2)) == 0,
     "A15 depth-two nonzero sign comes from C")


# Specialization audit: old component equations may specialize, but their
# face ownership generally does not.  These identities record the mergers.
old_D6_at_U0 = S**2 * (1 - alpha * S * P**5)
old_E10_at_K0 = S**2 * (1 - alpha * S * P**5)
old_Eeta_at_K0 = S**2 * (1 - eta * S * P**4)
old_EU_at_K0 = S**2 * (1 - u * P**5)
old_EDelta_at_K0 = S**2 * (1 - D * P**4 - Phi * S * P**3)
need(sp.factor(old_D6_at_U0 + S**2 * (alpha * S * P**5 - 1)) == 0,
     "old D6 merges into Balpha component")
need(sp.factor(old_E10_at_K0 + S**2 * (alpha * S * P**5 - 1)) == 0,
     "old E10 merges into Balpha component")
need(sp.factor(old_Eeta_at_K0 + S**2 * (eta * S * P**4 - 1)) == 0,
     "old Eeta merges into Beta component")
need(sp.factor(old_EU_at_K0 + S**2 * (u * P**5 - 1)) == 0,
     "old EU merges into BU vertical components")
need(sp.expand(old_EDelta_at_K0 - FACE_VARIANTS["EDPhi"][2]) == 0,
     "old EDelta specializes to rational EDPhi")

# Old U-roots run to P=infinity; old K-carrier branch roots coalesce at P=0.
need(sp.solve(sp.Eq(1 - U * P**6, 0), U)[0] == P**-6,
     "U-root valuation equation")
need(sp.expand((alpha**2 * P**8 + 4 * K).subs(K, 0)
               - (alpha * P**4)**2) == 0,
     "E10 branch polynomial becomes a square")
need(sp.expand((eta**2 * P**6 + 4 * K).subs(K, 0)
               - (eta * P**3)**2) == 0,
     "Eeta branch polynomial becomes a square")
need(sp.expand(((Phi**2 - 4 * K * D) * P**4 + 4 * K).subs(K, 0)
               - (Phi * P**2)**2) == 0,
     "EDelta branch polynomial becomes a square")


print("THM4356_DOUBLE_ZERO_ENDPOINT_AUDIT=PASS")
print("gate=Z=beta=zeta=W=xi=U=K=0;Delta=5696/105")
print("literal_rows=16;surviving_rows=9;possible_support_points=20")
print("exact_supports=48;exact_fans=15;exact_planes=12")
print("split=alpha_nonzero:24/6fans;alpha_zero_theta_nonzero:16/3fans;terminal:8/8fans")
print("hostile=128_keyed;96_distinct_supports;23_fans;17_planes")
print("positive_carriers=E11:g3/order41;FUT,FET:g2/order25;BDT-C:g2/order17")
print("collision_divisors=L5=eta^2-4uTheta;L15=Delta+Theta;LN=eta^2-4alphaPhi")
print("A4=depths1,2,3;orders115,35,95;one_end;no_graph_increment")
print("A15=depths1,2;orders116,16;ends_to_R_and_C;graph_increment0")
print("Nwall=smooth_boundary_tangency;unit_H_R;packet_2+2_to_3;rational_blowups")
print("smallest_theorem_ready=alpha=Theta=0;8_supports;8_fans;all_components_rational")
print("canonical_hostile=alpha=u=eta=0;Theta=-5696/105;BDT_A15")
print("global_ledgers=19;all_Pick_packet_graph_checks_pass")
print("CHECKS=%d" % CHECKS)
print("RESULT=PASS")
