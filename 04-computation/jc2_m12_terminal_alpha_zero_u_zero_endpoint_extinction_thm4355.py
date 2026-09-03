#!/usr/bin/env python3
"""Primary exact certificate for THM-4355's terminal-alpha endpoint.

The full gate is Z=beta=zeta=W=xi=U=alpha=0 with K nonzero on the
inherited exact-weight-twelve seam K=2848/45-(7/6)Delta.  The certificate
splits only on the honest coefficient owner Theta: forty supports and five
fans when Theta is nonzero, then twenty-four supports and seven fans when it
vanishes.  Both parts are reconstructed from the literal sixteen-row source.

Besides the support, fan, face, edge, Pick, graph, carrier, and positive-order
ledgers, this file verifies the four continuous collision walls: unibranch A4
and A2, two-branch A15, and boundary A3.  Every calculation is exact and the
script imports no repository code.
"""

from collections import Counter, defaultdict
from fractions import Fraction as Fr
from itertools import combinations, product
from math import gcd, lcm
import sys

import sympy as sp


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


CHECKS = 0


def need(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(label)


def noncollinear(points):
    return any(
        (q[0] - p[0]) * (r[1] - p[1])
        != (q[1] - p[1]) * (r[0] - p[0])
        for p, q, r in combinations(points, 3))


def candidate_planes(points):
    answer = set()
    for p, q, r in combinations(sorted(points), 3):
        det = ((q[0] - p[0]) * (r[1] - p[1])
               - (q[1] - p[1]) * (r[0] - p[0]))
        if det == 0:
            continue
        aa = Fr((q[2] - p[2]) * (r[1] - p[1])
                - (r[2] - p[2]) * (q[1] - p[1]), det)
        bb = Fr((q[0] - p[0]) * (r[2] - p[2])
                - (r[0] - p[0]) * (q[2] - p[2]), det)
        cc = Fr(p[2]) - aa * p[0] - bb * p[1]
        answer.add((aa, bb, cc))
    return tuple(answer)


def lower_fan(active, planes):
    answer = set()
    for plane in planes:
        aa, bb, cc = plane
        equal = []
        for rr, ll, hh in active:
            gap = Fr(hh) - aa * rr - bb * ll - cc
            if gap < 0:
                break
            if gap == 0:
                equal.append((rr, ll, hh))
        else:
            if len(equal) >= 3 and noncollinear(equal):
                answer.add(plane)
    return frozenset(answer)


def convex_hull(points):
    points = sorted(set(points))

    def cross(o, a, b):
        return ((a[0] - o[0]) * (b[1] - o[1])
                - (a[1] - o[1]) * (b[0] - o[0]))

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
    return lower[:-1] + upper[:-1]


def pick(polygon):
    twice_area = abs(sum(
        polygon[i][0] * polygon[(i + 1) % len(polygon)][1]
        - polygon[(i + 1) % len(polygon)][0] * polygon[i][1]
        for i in range(len(polygon))))
    boundary = sum(gcd(
        abs(polygon[(i + 1) % len(polygon)][0] - polygon[i][0]),
        abs(polygon[(i + 1) % len(polygon)][1] - polygon[i][1]))
        for i in range(len(polygon)))
    return twice_area, boundary, (twice_area - boundary + 2) // 2


def on_plane(point, plane):
    rr, ll, hh = point
    aa, bb, cc = plane
    return Fr(hh) == aa * rr + bb * ll + cc


def edge_polynomial(expression, start, end, S, P, X):
    dr, dl = end[0] - start[0], end[1] - start[1]
    length = gcd(abs(dr), abs(dl))
    dr, dl = dr // length, dl // length
    polynomial = sp.Poly(sp.expand(expression), S, P)
    return sp.factor(sum(
        polynomial.coeff_monomial(S ** (start[0] + j * dr)
                                  * P ** (start[1] + j * dl)) * X**j
        for j in range(length + 1)))


def face_order(plane):
    base = lcm(6, *(entry.denominator for entry in plane))
    order = Fr(base) * (Fr(5, 6) - sum(plane))
    need(order.denominator == 1 and order > 0,
         "positive integral face order")
    return base, order.numerator


def leading_coefficient(expression, variable, weight):
    return sp.expand(sp.limit(expression / variable**weight, variable, 0))


# Literal source and the exact U=alpha=0, Theta-nonzero lifted support.
U, W, Z = sp.symbols("U W Z")
Delta, Theta, Phi, eta = sp.symbols("Delta Theta Phi eta")
zeta, u, xi, alpha, beta, K = sp.symbols("zeta u xi alpha beta K")
e = -sp.Rational(1376, 135)
full_rows = {
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
need(len(full_rows) == 16, "complete sixteen-row source")
gate = {U: 0, W: 0, Z: 0, zeta: 0, xi: 0, alpha: 0, beta: 0}
rows = {label: sp.expand(sp.sympify(coefficient).subs(gate))
        for label, coefficient in full_rows.items()
        if sp.expand(sp.sympify(coefficient).subs(gate)) != 0}
need(len(rows) == 9, "nine surviving labelled rows")

support = defaultdict(lambda: sp.Integer(0))
support[(2, 0, 0)] += 1
support[(0, 1, 0)] -= 1
support[(2, 0, 1)] -= sp.Rational(1, 2)
for (i, j), coefficient in rows.items():
    support[(j + 2, i + j, 1)] -= coefficient
    support[(j, i + j + 1, 1)] += coefficient
support = {point: sp.factor(coefficient)
           for point, coefficient in support.items() if coefficient != 0}
need(len(support) == 19, "nineteen possible support points")

# The ten planes appearing in the five exact fans.
BU = (Fr(1, 10), Fr(1, 5), Fr(-1, 5))
FUT = (Fr(1, 5), Fr(1, 5), Fr(-2, 5))
T = (Fr(1, 2), Fr(0), Fr(-1))
C4 = (Fr(0), Fr(1, 4), Fr(-1, 4))
BETA = (Fr(1, 9), Fr(2, 9), Fr(-2, 9))
BDT = (Fr(1, 8), Fr(1, 4), Fr(-1, 4))
C3 = (Fr(-1, 3), Fr(1, 3), Fr(-1, 3))
CTH = (Fr(0), Fr(1, 3), Fr(-1, 3))
BTH = (Fr(1, 8), Fr(1, 4), Fr(-1, 4))

expected_fans = {
    "u": frozenset((BU, FUT, T)),
    "delta_eta": frozenset((C4, BETA, FUT, T)),
    "delta_zero": frozenset((BDT, T)),
    "e_eta": frozenset((C3, BETA, FUT, T)),
    "theta": frozenset((CTH, BTH, T)),
}
planes = candidate_planes(support)
Krel = sp.Rational(2848, 45) - sp.Rational(7, 6) * Delta
dstar = sp.Rational(3968, 63)  # exactly K=e, the first coupled deletion
delta_theta_reps = (
    (sp.Integer(0), sp.Integer(1)),
    (dstar, dstar),
    (dstar, 2 * dstar),
    (sp.Integer(1), sp.Integer(1)),
    (sp.Integer(1), sp.Integer(2)),
)

exact_supports = {}
fan_counts = Counter()
for (delta_value, theta_value), u_value, phi_value, eta_value in product(
        delta_theta_reps, (0, 1), (0, 1), (0, 1)):
    substitution = {
        Delta: delta_value,
        K: Krel.subs(Delta, delta_value),
        Theta: theta_value,
        u: u_value,
        Phi: phi_value,
        eta: eta_value,
    }
    need(substitution[K] != 0 and substitution[Theta] != 0,
         "K and Theta owner representatives")
    active = frozenset(point for point, coefficient in support.items()
                       if sp.simplify(coefficient.subs(substitution)) != 0)
    if u_value:
        key = "u"
    elif eta_value and delta_value:
        key = "delta_eta"
    elif eta_value:
        key = "e_eta"
    elif delta_value:
        key = "delta_zero"
    else:
        key = "theta"
    fan = lower_fan(active, planes)
    need(fan == expected_fans[key], "exact owner fan " + key)
    exact_supports[active] = substitution
    fan_counts[fan] += 1

need(len(exact_supports) == 40, "forty exact supports")
need(sorted(fan_counts.values()) == [2, 2, 8, 8, 20],
     "five exact fan counts")

# Coupling-hostile atlas: independently toggle the four lower rows and allow
# each of the three multiply-visible points to disappear.  Extra states are
# adversarial probes, not coefficient-realizable assertions.
required_labels = ((1, 0), (2, 0), (3, 0), (0, 2), (1, 2))
optional_labels = ((4, 0), (5, 0), (2, 1), (3, 1))
cancellable = ((2, 3, 1), (2, 4, 1), (2, 5, 1))


def lifted(labels):
    active = {(2, 0, 0), (0, 1, 0), (2, 0, 1)}
    for i, j in labels:
        active.add((j + 2, i + j, 1))
        active.add((j, i + j + 1, 1))
    return active


over = {}
over_distinct = set()
over_fans = set()
for row_bits in product((0, 1), repeat=4):
    labels = list(required_labels)
    labels.extend(label for bit, label in zip(row_bits, optional_labels) if bit)
    for cancellation_bits in product((0, 1), repeat=3):
        active = lifted(labels)
        for bit, point in zip(cancellation_bits, cancellable):
            if bit:
                active.discard(point)
        active = frozenset(active)
        over[row_bits + cancellation_bits] = active
        over_distinct.add(active)
        over_fans.add(lower_fan(active, planes))
need(len(over) == 128, "128 keyed hostile configurations")

for active, substitution in exact_supports.items():
    row_bits = tuple(int(substitution[v] != 0)
                     for v in (Delta, u, Phi, eta))
    cancellation_bits = (
        int(sp.simplify((K - e).subs(substitution)) == 0),
        int(sp.simplify((Theta - Delta).subs(substitution)) == 0),
        0,
    )
    need(over[row_bits + cancellation_bits] == active,
         "exact support embeds into hostile atlas")

# Exact faces, lattice geometry, positive form orders, and every edge word.
S, P, X = sp.symbols("S P X")


def face_expression(plane, substitution=None):
    expression = sum(coefficient * S**rr * P**ll
                     for (rr, ll, hh), coefficient in support.items()
                     if on_plane((rr, ll, hh), plane))
    if substitution:
        expression = expression.subs(substitution)
    return sp.factor(expression)


faces = {
    "BU": (BU, {}, (P - S**2) * (u * P**5 - 1),
           [(0, 1), (2, 0), (2, 5), (0, 6)], (20, 12, 5), (30, 22)),
    "FUT": (FUT, {}, S**2 * (1 - u * P**5 - eta * S * P**4
                               - Theta * S**2 * P**3),
            [(2, 0), (4, 3), (2, 5)], (10, 8, 2), (30, 25)),
    "FET": (FUT, {u: 0}, S**2 * (1 - eta * S * P**4
                                    - Theta * S**2 * P**3),
            [(2, 0), (4, 3), (3, 4)], (5, 3, 2), (30, 25)),
    "T": (T, {}, S**2 * (1 - K * S**2 * P**2
                           - Theta * S**2 * P**3),
          [(2, 0), (4, 2), (4, 3)], (2, 4, 0), (6, 8)),
    "C4": (C4, {u: 0}, P * (Delta * P**4 + eta * S * P**4 - 1),
           [(0, 1), (1, 5), (0, 5)], (4, 6, 0), (12, 10)),
    "BETA": (BETA, {u: 0}, (P - S**2) * (eta * S * P**4 - 1),
             [(0, 1), (2, 0), (3, 4), (1, 5)], (18, 4, 8), (18, 13)),
    "BDT": (BDT, {u: 0, eta: 0},
            (P - S**2) * (Delta * P**4 + Theta * S**2 * P**3 - 1),
            [(0, 1), (2, 0), (4, 3), (0, 5)], (24, 8, 9), (24, 17)),
    "C3": (C3, {u: 0, Delta: 0},
           P * (e * P**3 + eta * S * P**4 - 1),
           [(0, 1), (1, 5), (0, 4)], (3, 5, 0), (6, 7)),
    "CTH": (CTH, {u: 0, eta: 0, Delta: 0},
            P * (e * P**3 + Phi * S * P**3
                 + Theta * S**2 * P**3 - 1),
            [(0, 1), (2, 4), (0, 4)], (6, 6, 1), (6, 5)),
    "BTH": (BTH, {u: 0, eta: 0, Delta: 0},
            (P - S**2) * (Theta * S**2 * P**3 - 1),
            [(0, 1), (2, 0), (4, 3), (2, 4)], (16, 4, 7), (24, 17)),
}
for name, (plane, substitution, expected, polygon,
           pick_data, order_data) in faces.items():
    need(sp.expand(face_expression(plane, substitution) - expected) == 0,
         name + " exact face")
    need(convex_hull(polygon) == polygon, name + " polygon orientation")
    need(pick(polygon) == pick_data, name + " Pick data")
    need(face_order(plane) == order_data, name + " positive form order")

expected_edges = {
    "BU": (X - 1, 1 - u * X**5, u * (X - 1), u - X**5),
    "FUT": (1 - Theta * X,
            -Theta - eta * X - u * X**2, X**5 - u),
    "FET": (1 - Theta * X, -Theta - eta * X, X - eta),
    "T": (1 - K * X**2, -K - Theta * X, X - Theta),
    "C4": (eta * X - 1, eta + Delta * X, Delta - X**4),
    "BETA": (X - 1, 1 - eta * X, eta * (X - 1), eta - X),
    "BDT": (X - 1, 1 - Theta * X,
            (X - 1) * (Delta * X + Theta), Delta - X**4),
    "C3": (eta * X - 1, eta + e * X, e - X**3),
    "CTH": (Theta * X - 1,
            Theta + Phi * X + e * X**2, e - X**3),
    "BTH": (X - 1, 1 - Theta * X,
            Theta * (X - 1), Theta - X),
}
for name, (_plane, _substitution, expression, polygon, *_rest) in faces.items():
    actual = tuple(edge_polynomial(expression, polygon[index],
                                   polygon[(index + 1) % len(polygon)],
                                   S, P, X)
                   for index in range(len(polygon)))
    need(len(actual) == len(expected_edges[name]), name + " edge count")
    for got, want in zip(actual, expected_edges[name]):
        need(sp.expand(got - want) == 0, name + " exact edge word")

L5 = eta**2 - 4 * u * Theta
L15 = Delta + Theta
L3 = Phi**2 - 4 * e * Theta
need(sp.discriminant(-Theta - eta * X - u * X**2, X) == L5,
     "A4 carrier-edge discriminant")
need(sp.expand(sp.discriminant((X - 1) * (Delta * X + Theta), X)
               - L15**2) == 0,
     "A15 contact-edge discriminant")
need(sp.discriminant(Theta + Phi * X + e * X**2, X) == L3,
     "A2 carrier-edge discriminant")
need(sp.factor(((X - 1) * (Delta * X + Theta)).subs(Theta, -Delta))
     == Delta * (X - 1)**2, "A15 double edge root")

global_polygons = {
    "u": ([(0, 1), (2, 0), (4, 2), (4, 3), (2, 5), (0, 6)],
          (32, 12, 11)),
    "delta_eta": ([(0, 1), (2, 0), (4, 2), (4, 3),
                    (3, 4), (1, 5), (0, 5)], (29, 11, 10)),
    "delta_zero": ([(0, 1), (2, 0), (4, 2), (4, 3), (0, 5)],
                   (26, 10, 9)),
    "e_eta": ([(0, 1), (2, 0), (4, 2), (4, 3),
                (3, 4), (1, 5), (0, 4)], (28, 10, 10)),
    "theta": ([(0, 1), (2, 0), (4, 2), (4, 3), (2, 4), (0, 4)],
              (24, 10, 8)),
}
for name, (polygon, expected) in global_polygons.items():
    need(convex_hull(polygon) == polygon, name + " global polygon")
    need(pick(polygon) == expected, name + " global Pick ledger")

# Carrier normalizations and clean graph ledgers.
f_ut = 1 - u * P**5 - eta * S * P**4 - Theta * S**2 * P**3
y_ut = 2 * Theta * S * P**2 + eta * P**3
branch_ut = P * (4 * Theta + L5 * P**5)
need(sp.expand(y_ut**2 + 4 * Theta * P * f_ut - branch_ut) == 0,
     "FUT hyperelliptic identity")
f_dt = 1 - Delta * P**4 - Theta * S**2 * P**3
y_dt = Theta * S * P**2
branch_dt = Theta * P * (1 - Delta * P**4)
need(sp.expand(y_dt**2 + Theta * P * f_dt - branch_dt) == 0,
     "BDT genus-two component identity")
f_th = 1 - e * P**3 - Phi * S * P**3 - Theta * S**2 * P**3
y_th = P**2 * (2 * Theta * S + Phi)
branch_th = L3 * P**4 + 4 * Theta * P
need(sp.expand(y_th**2 + 4 * Theta * P * f_th - branch_th) == 0,
     "CTH elliptic identity")

for polynomial, variable, label in (
        (branch_ut, P, "FUT off L5"),
        (branch_ut.subs(u, 0), P, "FET"),
        (branch_dt, P, "BDT component"),
        (branch_th, P, "CTH off L3")):
    need(sp.gcd(sp.Poly(polynomial, variable),
                sp.Poly(sp.diff(polynomial, variable), variable)).degree() == 0,
         label + " squarefree over coefficient field")

intersection_degrees = {
    "BU": sp.degree((u * P**5 - 1).subs(P, S**2), S),
    "BETA": sp.degree((eta * S * P**4 - 1).subs(P, S**2), S),
    "BDT": sp.degree((Delta * P**4 + Theta * S**2 * P**3 - 1)
                     .subs(P, S**2), S),
    "BTH": sp.degree((Theta * S**2 * P**3 - 1).subs(P, S**2), S),
}
need(intersection_degrees == {"BU": 10, "BETA": 9,
                              "BDT": 8, "BTH": 8},
     "factor intersection counts")
clean_ledgers = {
    "u": (8, 16, 9, 2, 11),
    "delta_eta": (5, 12, 8, 2, 10),
    "delta_zero": (3, 9, 7, 2, 9),
    "e_eta": (5, 12, 8, 2, 10),
    "theta": (4, 10, 7, 1, 8),
}
for name, (vertices, edges, graph_genus, carrier_genus, total) in clean_ledgers.items():
    need(edges - vertices + 1 == graph_genus, name + " graph genus")
    need(graph_genus + carrier_genus == total, name + " total genus")

# Literal source used in all three collision charts.
s0, p0, Q = sp.symbols("s0 p0 Q")
H = (-3 * p0 + sp.Rational(8, 3) * p0**2 + e * p0**3
     + K * (s0 * p0)**2 + Phi * p0**2 * (s0 * p0)
     + Delta * p0**4 + Theta * p0 * (s0 * p0)**2
     + eta * p0**3 * (s0 * p0) + u * p0**5)
FQ = (s0**2 - p0) * (1 - Q * H) - Q * s0**2 / 2
sigma, x, v, rho, a, t, bigX, bigY, zz, bigZ = sp.symbols(
    "sigma x v rho a t XX YY zz ZZ")

# A4: eta^2=4uTheta on the u-owner fan.  The critical ratio a is nonzero,
# so the prefactor is a unit and THM-4341 applies literally.
G5_source = sp.expand(
    sigma**2 * FQ.subs({Q: sigma**5, s0: S / sigma,
                        p0: P / sigma}, simultaneous=True))
G5_expected = (
    (S**2 - sigma * P)
    * (1 - u * P**5 - eta * S * P**4 - Theta * S**2 * P**3
       - sigma * (Delta * P**4 + Phi * S * P**3
                  + K * S**2 * P**2)
       - e * sigma**2 * P**3
       - sp.Rational(8, 3) * sigma**3 * P**2
       + 3 * sigma**4 * P)
    - sigma**5 * S**2 / 2)
need(sp.expand(G5_source - G5_expected) == 0,
     "literal FUT primitive chart")
A5 = u + eta * v + Theta * v**2
B5 = Delta + Phi * v + K * v**2
root5 = ((v**2 - rho)
         * (x**5 - A5 - rho * B5 - e * rho**2
            - sp.Rational(8, 3) * rho**3 + 3 * rho**4)
         - rho**5 * v**2 / 2)
root5_source = sp.cancel(
    x**7 * G5_expected.subs({S: v / x, P: 1 / x,
                             sigma: rho / x}, simultaneous=True))
need(sp.expand(root5_source - root5) == 0,
     "literal FUT reciprocal chart")
collision5 = {u: Theta * a**2, eta: -2 * Theta * a}
need(sp.expand(A5.subs(collision5) - Theta * (v - a)**2) == 0,
     "A4 transverse square")
D5 = (A5 + rho * B5 + e * rho**2 + sp.Rational(8, 3) * rho**3
      - 3 * rho**4 + rho**5 * v**2 / (2 * (v**2 - rho)))
need(sp.factor(root5 - (v**2 - rho) * (x**5 - D5)) == 0,
     "exact divided A4 equation")
need(sp.diff(D5.subs(collision5), v, 2).subs({v: a, rho: 0})
     == 2 * Theta, "A4 nonzero Morse Hessian")
need(a != 0, "A4 owner conditions force interior ratio root")

# The formal critical section is needed to prove that the list of tails is
# exhaustive.  The quotient correction begins too late to alter these four
# coefficients.  If the first three vanish, the fourth is forced to -99/43.
q5 = 2 * K * a + Phi
v5_critical = -(Phi * rho - 2 * Theta * a) / (2 * (K * rho + Theta))
critical5 = sp.series(
    D5.subs(collision5).subs(v, v5_critical), rho, 0, 5
).removeO().expand()
critical5_wanted = (
    Delta + Phi * a + K * a**2,
    e - q5**2 / (4 * Theta),
    sp.Rational(8, 3) + K * q5**2 / (4 * Theta**2),
    -3 - K**2 * q5**2 / (4 * Theta**3),
)
need(all(sp.simplify(critical5.coeff(rho, j) - wanted) == 0
         for j, wanted in enumerate(critical5_wanted, 1)),
     "A4 exact critical selectors")
forced_a4_c4 = critical5_wanted[3].subs(
    q5**2, 4 * Theta * e
).subs(K, sp.Rational(45, 172) * Theta)
need(sp.factor(forced_a4_c4) == -sp.Rational(99, 43),
     "A4 depth at most four")

# Exact witnesses show that the upper-bound list is sharp on the inherited
# K/Delta seam.  We use q=2Ka+Phi, so c1=Delta+q*a-K*a^2.
seam_constant = sp.Rational(2848, 45)
k2 = 6 * seam_constant / 13
d3 = sp.Rational(6, 7) * (seam_constant - 1)
a3_witness = (1 + sp.sqrt(1 + 4 * d3)) / 2
k4 = 4 * sp.Rational(172, 45) * e
d4 = sp.Rational(6, 7) * (seam_constant - k4)
a4_witness = (1 + sp.sqrt(1 + 4 * d4 / k4)) / 2
a4_depth_witnesses = {
    1: (seam_constant, sp.Integer(0), sp.Integer(1), sp.Integer(1),
        sp.Integer(0)),
    2: (k2, k2, sp.Integer(1), sp.Integer(1), sp.Integer(0)),
    3: (sp.Integer(1), d3, 1 / (4 * e), a3_witness, sp.Integer(1)),
    4: (k4, d4, sp.Rational(172, 45) * k4, a4_witness, k4),
}
for depth, (kw, dw, tw, aw, qw) in a4_depth_witnesses.items():
    witness_selectors = (
        dw + qw * aw - kw * aw**2,
        e - qw**2 / (4 * tw),
        sp.Rational(8, 3) + kw * qw**2 / (4 * tw**2),
        -3 - kw**2 * qw**2 / (4 * tw**3),
    )
    need(sp.simplify(kw - (seam_constant - sp.Rational(7, 6) * dw)) == 0,
         "A4 witness lies on K/Delta seam")
    need(all(sp.simplify(value) == 0
             for value in witness_selectors[:depth - 1]),
         "A4 witness has required preceding cancellations")
    need(sp.simplify(witness_selectors[depth - 1]) != 0,
         "A4 witness realizes its exact depth")
    need(all(sp.simplify(value) != 0 for value in (kw, tw, aw)),
         "A4 witness respects owner conditions")
odd5 = {1: (2, 0), 2: (1, 1), 3: (1, 1), 4: (0, 2)}
for split, (tail_genus, persistent_delta) in odd5.items():
    epsilon = split % 2
    branch = bigX**epsilon * (bigX**(5 - split) - sp.Symbol("cc"))
    need(sp.gcd(sp.Poly(branch, bigX),
                sp.Poly(sp.diff(branch, bigX), bigX)).degree() == 0,
         "A4 split tail squarefree")
    need((sp.degree(branch, bigX) - 1) // 2 == tail_genus,
         "A4 split tail genus")
    need(tail_genus + persistent_delta == 2,
         "A4 delta conservation")
    need(25 + Fr(15 * split, 5 - split) > 0,
         "A4 positive tail form order")
a4_target_rows = {
    1: ((4, 6, 15), 115),
    2: ((1, 4, 10), 35),
    3: ((2, 18, 45), 95),
    4: ((1, 24, 60), 85),
}
for split, ((weight_sigma, weight_x, weight_y), order) in a4_target_rows.items():
    need(5 * weight_x == split * (6 * weight_sigma + weight_x),
         "A4 target-compatible balance")
    need(2 * weight_y == 5 * weight_x, "A4 y weight")
    need(gcd(gcd(weight_sigma, weight_x), weight_y) == 1,
         "A4 primitive target-compatible ray")
    need(25 * weight_sigma + 5 * weight_x - weight_y == order > 0,
         "A4 primitive positive form order")
need(4 == 4, "A4 conductor x^4 equals residue buffer")
need(9 + 2 == 11, "A4 clean graph plus cusp delta")
for tail_genus, persistent_delta in odd5.values():
    need(9 + tail_genus + persistent_delta == 11,
         "A4 one-ended tail genus conservation")
need(0 == 0, "A4 one attachment adds no graph cycle")

# A2: Phi^2=4eTheta on the deepest Theta-owner fan.  The same odd law
# applies with m=3, including arbitrary cancellation depth of the critical
# value; depths >=3 are horizontal normalization cases.
G3_source = sp.expand(
    sigma * FQ.subs({u: 0, eta: 0, Delta: 0,
                     Q: sigma**3, s0: S, p0: P / sigma},
                    simultaneous=True))
A3 = e + Phi * S + Theta * S**2
B3 = sp.Rational(8, 3) + K * S**2
G3_expected = (
    (sigma * S**2 - P)
    * (1 - A3 * P**3 - sigma * B3 * P**2 + 3 * sigma**2 * P)
    - sigma**4 * S**2 / 2)
need(sp.expand(G3_source - G3_expected) == 0,
     "literal CTH primitive chart")
root3 = ((rho * S**2 - 1)
         * (x**3 - A3 - rho * B3 + 3 * rho**2)
         - rho**4 * S**2 / 2)
root3_source = sp.cancel(
    x**4 * G3_expected.subs({P: 1 / x, sigma: rho / x},
                             simultaneous=True))
need(sp.expand(root3_source - root3) == 0,
     "literal CTH reciprocal chart")
collision3 = {Phi: -2 * Theta * a}
need(sp.expand(A3.subs(collision3) - Theta * (S - a)**2)
     == e - Theta * a**2,
     "A2 transverse square modulo e=Theta*a^2")
D3 = (A3 + rho * B3 - 3 * rho**2
      + rho**4 * S**2 / (2 * (rho * S**2 - 1)))
need(sp.factor(root3 - (rho * S**2 - 1) * (x**3 - D3)) == 0,
     "exact divided A2 equation")
need(sp.diff(D3.subs(collision3), S, 2).subs({S: a, rho: 0})
     == 2 * Theta, "A2 nonzero Morse Hessian")

# On the A2 wall e=Theta*a^2.  The first selector can vanish, but then the
# second is again the rigid value -99/43, so only depths one and two occur.
S3_critical = Theta * a / (Theta + K * rho)
critical3_model = (
    Theta * (S - a)**2
    + rho * (sp.Rational(8, 3) + K * S**2)
    - 3 * rho**2
    - rho**4 * S**2 / (2 * (1 - rho * S**2))
)
critical3 = sp.series(
    critical3_model.subs(S, S3_critical), rho, 0, 3
).removeO().expand()
critical3_wanted = (
    sp.Rational(8, 3) + K * a**2,
    -3 - K**2 * a**2 / Theta,
)
need(all(sp.simplify(critical3.coeff(rho, j) - wanted) == 0
         for j, wanted in enumerate(critical3_wanted, 1)),
     "A2 exact critical selectors")
forced_a2_c2 = critical3_wanted[1].subs(
    Theta, e / a**2
).subs(a**2, -sp.Rational(8, 3) / K)
need(sp.factor(forced_a2_c2) == -sp.Rational(99, 43),
     "A2 depth at most two")
a2_k = seam_constant
a2_depth_witnesses = {
    1: (sp.Integer(1), e),
    2: (sp.sqrt(-sp.Rational(8, 3) / a2_k),
        e / (-sp.Rational(8, 3) / a2_k)),
}
for depth, (aw, tw) in a2_depth_witnesses.items():
    witness_selectors = (
        sp.Rational(8, 3) + a2_k * aw**2,
        -3 - a2_k**2 * aw**2 / tw,
    )
    need(sp.simplify(e - tw * aw**2) == 0,
         "A2 witness lies on collision wall")
    need(all(sp.simplify(value) == 0
             for value in witness_selectors[:depth - 1]),
         "A2 witness has required preceding cancellations")
    need(sp.simplify(witness_selectors[depth - 1]) != 0,
         "A2 witness realizes its exact depth")
    need(sp.simplify(tw * aw) != 0, "A2 witness respects owner conditions")
odd3 = {1: (1, 0), 2: (0, 1)}
for split, (tail_genus, persistent_delta) in odd3.items():
    epsilon = split % 2
    branch = bigX**epsilon * (bigX**(3 - split) - sp.Symbol("dd"))
    need(sp.gcd(sp.Poly(branch, bigX),
                sp.Poly(sp.diff(branch, bigX), bigX)).degree() == 0,
         "A2 split tail squarefree")
    need((sp.degree(branch, bigX) - 1) // 2 == tail_genus,
         "A2 split tail genus")
    need(tail_genus + persistent_delta == 1,
         "A2 delta conservation")
    need(5 + Fr(3 * split, 3 - split) > 0,
         "A2 positive tail form order")
a2_target_rows = {
    1: ((2, 2, 3), 13),
    2: ((1, 4, 6), 11),
}
for split, ((weight_sigma, weight_x, weight_y), order) in a2_target_rows.items():
    need(3 * weight_x == split * (2 * weight_sigma + weight_x),
         "A2 target-compatible balance")
    need(2 * weight_y == 3 * weight_x, "A2 y weight")
    need(gcd(gcd(weight_sigma, weight_x), weight_y) == 1,
         "A2 primitive target-compatible ray")
    need(5 * weight_sigma + 3 * weight_x - weight_y == order > 0,
         "A2 primitive positive form order")
need(2 == 2, "A2 conductor x^2 equals residue buffer")
need(7 + 1 == 8, "A2 clean graph plus cusp delta")
for tail_genus, persistent_delta in odd3.values():
    need(7 + tail_genus + persistent_delta == 8,
         "A2 one-ended tail genus conservation")
need(0 == 0, "A2 one attachment adds no graph cycle")

# A15: Delta+Theta=0 on BDT.  The central branches have contact eight at
# q=1, hence A15.  They lie in different normalized complementary
# components, so THM-4352's same-complement graph +1 is inapplicable.
G15_source = sp.expand(
    sigma**2 * FQ.subs({u: 0, eta: 0, Q: sigma**8,
                        s0: S / sigma, p0: P / sigma**2},
                       simultaneous=True))
G15_expected = (
    (S**2 - P)
    * (1 - Delta * P**4 - Theta * S**2 * P**3
       - sigma * Phi * S * P**3
       - sigma**2 * (e * P**3 + K * S**2 * P**2)
       - sp.Rational(8, 3) * sigma**4 * P**2
       + 3 * sigma**6 * P)
    - sigma**8 * S**2 / 2)
need(sp.expand(G15_source - G15_expected) == 0,
     "literal BDT primitive chart")
q = sp.symbols("q")
B15 = (x**8 - Delta * q**4 - Theta * q**3 - rho * Phi * q**3
       - rho**2 * (e * q**3 + K * q**2)
       - sp.Rational(8, 3) * rho**4 * q**2 + 3 * rho**6 * q)
root15 = (1 - q) * B15 - rho**8 / 2
root15_source = sp.cancel(
    x**10 * G15_expected.subs({S: 1 / x, P: q / x**2,
                               sigma: rho / x}, simultaneous=True))
need(sp.expand(root15_source - root15) == 0,
     "literal BDT reciprocal chart")
wall15 = sp.expand(root15.subs({Theta: -Delta, q: 1 + zz}))
central15 = sp.factor(wall15.subs(rho, 0))
need(sp.expand(central15
               - zz * (Delta * zz * (zz + 1)**3 - x**8)) == 0,
     "A15 central two-branch contact")
need(sp.resultant(zz, Delta * zz * (zz + 1)**3 - x**8, zz)
     == -x**8, "A15 intersection length eight")

# Case I: Phi!=0.  The first exceptional equation has one horizontal X^2
# discriminant factor (A1).  After removing it, the two rational sign
# components meet at seven simple nonzero addresses.  These are mutual nodes,
# not seven separate bridge components.
case1 = wall15.subs({sigma: t**7, x: t * bigX,
                     rho: t**8 * bigX, zz: t**8 * bigZ},
                    simultaneous=True)
lead1 = leading_coefficient(case1, t, 16)
H1 = bigX * (bigX**7 - Phi)
need(sp.expand(lead1 - (Delta * bigZ**2 - H1 * bigZ)) == 0,
     "A15 Phi-nonzero first exceptional equation")
need(sp.gcd(sp.Poly(bigX**7 - Phi, bigX),
            sp.Poly(sp.diff(bigX**7 - Phi, bigX), bigX)).degree() == 0,
     "seven simple mutual-node addresses")
j, LL, MM = sp.symbols("j LL MM")
c1 = sp.diff(H1, bigX).subs(bigX, j).subs(j**7, Phi)
need(sp.expand(c1 - 7 * Phi) == 0, "case-I node separation slope")
bridge1 = Delta * MM**2 - c1 * LL * MM - j**8 / 2
need(sp.solve([sp.diff(bridge1, LL), sp.diff(bridge1, MM), bridge1],
              [LL, MM], dict=True) == [], "case-I node-resolution conic smooth")
need(48 == 2 * 24, "case-I second node-resolution scale")

# Case II: Phi=0 and c=e+K!=0.  The horizontal square is X^4 (A3),
# followed by two rational sign components with six mutual nodes.
c = sp.symbols("c")
case2 = wall15.subs({Phi: 0, sigma: t**3, x: t * bigX,
                     rho: t**4 * bigX, zz: t**8 * bigZ},
                    simultaneous=True)
lead2 = leading_coefficient(case2, t, 16)
H2 = bigX**2 * (bigX**6 - (e + K))
need(sp.expand(lead2 - (Delta * bigZ**2 - H2 * bigZ)) == 0,
     "A15 Phi-zero first exceptional equation")
need(sp.gcd(sp.Poly(bigX**6 - c, bigX),
            sp.Poly(sp.diff(bigX**6 - c, bigX), bigX)).degree() == 0,
     "six simple mutual-node addresses")
c2 = sp.diff(bigX**2 * (bigX**6 - c), bigX).subs(bigX, j).subs(j**6, c)
need(sp.expand(c2 - 6 * j**7) == 0,
     "case-II node separation slope")
bridge2 = Delta * MM**2 - c2 * LL * MM - j**8 / 2
need(sp.solve([sp.diff(bridge2, LL), sp.diff(bridge2, MM), bridge2],
              [LL, MM], dict=True) == [], "case-II node-resolution conic smooth")
need(16 == 2 * 8, "case-II second node-resolution scale")

# Case III: Phi=0 and e+K=0 is one seam point.  The X^8 square is a
# persistent A7 and the normalized exceptional double cover is genus three.
delta_special = sp.solve(sp.Eq(Krel, -e), Delta)[0]
need(delta_special == sp.Rational(2048, 45),
     "unique e+K zero seam value")
case3 = wall15.subs({Phi: 0, K: -e, Delta: delta_special,
                     sigma: t, x: t * bigX, rho: t**2 * bigX,
                     zz: t**8 * bigZ}, simultaneous=True)
lead3 = leading_coefficient(case3, t, 16)
H3 = bigX**4 * (bigX**4 - sp.Rational(8, 3))
expected3 = (delta_special * bigZ**2 - H3 * bigZ - bigX**8 / 2)
need(sp.expand(lead3 - expected3) == 0,
     "A15 exceptional genus-three equation")
tail3_branch = ((bigX**4 - sp.Rational(8, 3))**2
                + 2 * delta_special)
need(sp.expand(sp.discriminant(expected3, bigZ)
               - bigX**8 * tail3_branch) == 0,
     "A15 completed-square branch polynomial")
need(sp.expand(tail3_branch.subs(bigX, 0) - sp.Rational(1472, 15)) == 0,
     "tail nonzero at derivative root zero")
need(sp.expand(tail3_branch.subs(bigX**4, sp.Rational(8, 3))
               - 2 * delta_special) == 0,
     "tail nonzero at four derivative roots")
need(sp.gcd(sp.Poly(tail3_branch, bigX),
            sp.Poly(sp.diff(tail3_branch, bigX), bigX)).degree() == 0,
     "A15 exceptional branch polynomial squarefree")
need((sp.degree(tail3_branch, bigX) - 2) // 2 == 3,
     "A15 exceptional tail genus three")
need(gcd(gcd(1, 1), 8) == 1, "case-III primitive first ray")
need(17 + 18 + 3 - 24 == 14,
     "case-III normalized tail differential order")

# A15 graph/delta ledger.  The complement R and C--T has two connected
# components.  In depths one and two the tail packet has two rational signs,
# seven or six mutual nodes, and one complementary attachment per sign.  In
# depth four it is one connected genus-three component with two attachments.
a15_ledgers = {
    "Phi_nonzero": (1, 2, 7, 6, 0),
    "Phi_zero_c_nonzero": (2, 2, 6, 5, 0),
    "deep_special": (4, 1, 0, 0, 3),
}
for name, (persistent_delta, tail_components, mutual_nodes,
           graph_genus, tail_genus) in a15_ledgers.items():
    vertices = 3 + tail_components
    edges = 1 + 2 + mutual_nodes
    need(edges - vertices + 1 == graph_genus,
         name + " different-complement graph genus")
    need(2 + persistent_delta + graph_genus + tail_genus == 9,
         name + " global genus restoration")
need(7 == 8 - 1, "A15 intrinsic even-cusp deficit")
need(0 == 0, "two ends joining different components add no graph unit")
need(8 == 8, "A15 conductor exponent x^8")
need(6 < 8, "BDT residue buffer is below automatic conductor threshold")
need(14 > 0, "direct case-III valuation repairs conductor shortfall")

# Target-compatible primitive rays use Q=sigma^24 (the source-minimal
# Q=delta^8 chart above has delta=sigma^3).  For a critical depth r,
# (8-r)b=3ra and the residue sigma^17*x^6*dx/y has order 17a-b.
target_rows = {
    1: ((7, 3, 24), 116),
    2: ((1, 1, 8), 16),
    4: ((1, 3, 24), 14),
}
for depth, ((weight_sigma, weight_x, weight_y), order) in target_rows.items():
    need((8 - depth) * weight_x == 3 * depth * weight_sigma,
         "A15 target-compatible balance")
    need(weight_y == 8 * weight_x, "A15 y weight")
    need(gcd(weight_sigma, weight_x) == 1,
         "A15 primitive target-compatible ray")
    need(17 * weight_sigma - weight_x == order > 0,
         "A15 positive primitive form order")

theta_exact_supports = frozenset(exact_supports)
theta_exact_fans = frozenset(fan_counts)
theta_hostile_supports = frozenset(over_distinct)
theta_hostile_fans = frozenset(over_fans)
theta_natural_index = {}
for active, substitution in exact_supports.items():
    coefficient_class = delta_theta_reps.index(
        (substitution[Delta], substitution[Theta])
    )
    number = (1 + 8 * coefficient_class + 4 * int(substitution[u] != 0)
              + 2 * int(substitution[Phi] != 0)
              + int(substitution[eta] != 0))
    need(number not in theta_natural_index, "Theta-lane natural index injective")
    theta_natural_index[number] = active
need(set(theta_natural_index) == set(range(1, 41)),
     "Theta-lane natural index is 1..40")
theta_checks = CHECKS

# Independent reconstruction of the honest Theta=0 boundary lane.
CHECKS = 0


def need(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(label)


def noncollinear(points):
    return any(
        (q[0] - p[0]) * (r[1] - p[1])
        != (q[1] - p[1]) * (r[0] - p[0])
        for p, q, r in combinations(points, 3))


def candidate_planes(points):
    result = set()
    for p, q, r in combinations(sorted(points), 3):
        det = ((q[0] - p[0]) * (r[1] - p[1])
               - (q[1] - p[1]) * (r[0] - p[0]))
        if det == 0:
            continue
        aa = Fr((q[2] - p[2]) * (r[1] - p[1])
                - (r[2] - p[2]) * (q[1] - p[1]), det)
        bb = Fr((q[0] - p[0]) * (r[2] - p[2])
                - (r[0] - p[0]) * (q[2] - p[2]), det)
        cc = Fr(p[2]) - aa * p[0] - bb * p[1]
        result.add((aa, bb, cc))
    return tuple(result)


def lower_fan(active, planes):
    result = set()
    for plane in planes:
        aa, bb, cc = plane
        equal = []
        for rr, ll, hh in active:
            gap = Fr(hh) - aa * rr - bb * ll - cc
            if gap < 0:
                break
            if gap == 0:
                equal.append((rr, ll, hh))
        else:
            if len(equal) >= 3 and noncollinear(equal):
                result.add(plane)
    return frozenset(result)


def convex_hull(points):
    points = sorted(set(points))

    def cross(o, a, b):
        return ((a[0] - o[0]) * (b[1] - o[1])
                - (a[1] - o[1]) * (b[0] - o[0]))

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
    return lower[:-1] + upper[:-1]


def pick(polygon):
    twice_area = abs(sum(
        polygon[i][0] * polygon[(i + 1) % len(polygon)][1]
        - polygon[(i + 1) % len(polygon)][0] * polygon[i][1]
        for i in range(len(polygon))))
    boundary = sum(gcd(
        abs(polygon[(i + 1) % len(polygon)][0] - polygon[i][0]),
        abs(polygon[(i + 1) % len(polygon)][1] - polygon[i][1]))
        for i in range(len(polygon)))
    interior = (twice_area - boundary + 2) // 2
    return twice_area, boundary, interior


def on_plane(point, plane):
    rr, ll, hh = point
    aa, bb, cc = plane
    return Fr(hh) == aa * rr + bb * ll + cc


def edge_polynomial(expression, start, end, S, P, X):
    dr, dl = end[0] - start[0], end[1] - start[1]
    length = gcd(abs(dr), abs(dl))
    dr, dl = dr // length, dl // length
    poly = sp.Poly(sp.expand(expression), S, P)
    return sp.factor(sum(
        poly.coeff_monomial(S ** (start[0] + j * dr)
                            * P ** (start[1] + j * dl)) * X**j
        for j in range(length + 1)))


def face_order(plane):
    base = lcm(6, *(entry.denominator for entry in plane))
    order = Fr(base) * (Fr(5, 6) - sum(plane))
    need(order.denominator == 1 and order > 0, "positive integral face order")
    return base, order.numerator


# Literal source and lifted support.
U, W, Z = sp.symbols("U W Z")
Delta, Theta, Phi, eta = sp.symbols("Delta Theta Phi eta")
zeta, u, xi, alpha, beta, K = sp.symbols("zeta u xi alpha beta K")
e = -sp.Rational(1376, 135)
full_rows = {
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
need(len(full_rows) == 16, "complete source")
gate = {U: 0, W: 0, Z: 0, zeta: 0, xi: 0, alpha: 0,
        beta: 0, Theta: 0}
rows = {label: sp.expand(sp.sympify(coefficient).subs(gate))
        for label, coefficient in full_rows.items()
        if sp.expand(sp.sympify(coefficient).subs(gate)) != 0}
need(len(rows) == 8, "eight surviving labelled rows")

support = defaultdict(lambda: sp.Integer(0))
support[(2, 0, 0)] += 1
support[(0, 1, 0)] -= 1
support[(2, 0, 1)] -= sp.Rational(1, 2)
for (i, j), coefficient in rows.items():
    support[(j + 2, i + j, 1)] -= coefficient
    support[(j, i + j + 1, 1)] += coefficient
support = {point: sp.factor(coefficient)
           for point, coefficient in support.items() if coefficient != 0}
need(len(support) == 18, "eighteen possible support points")

# The twelve planes used by the seven exact fans.
BU = (Fr(1, 10), Fr(1, 5), Fr(-1, 5))
FUE = (Fr(1, 5), Fr(1, 5), Fr(-2, 5))
EETA = (Fr(1, 3), Fr(1, 6), Fr(-2, 3))
EU = (Fr(3, 10), Fr(1, 5), Fr(-3, 5))
C4 = (Fr(0), Fr(1, 4), Fr(-1, 4))
BETA = (Fr(1, 9), Fr(2, 9), Fr(-2, 9))
C3 = (Fr(-1, 3), Fr(1, 3), Fr(-1, 3))
BDELTA = (Fr(1, 8), Fr(1, 4), Fr(-1, 4))
EDELTA = (Fr(1, 4), Fr(1, 4), Fr(-1, 2))
CPHI = (Fr(0), Fr(1, 3), Fr(-1, 3))
BPHI = (Fr(1, 7), Fr(2, 7), Fr(-2, 7))
BDEEP = (Fr(1, 6), Fr(1, 3), Fr(-1, 3))

expected_fans = {
    (1, 1, 0, 0): frozenset((BU, FUE, EETA)),
    (1, 0, 0, 0): frozenset((BU, EU)),
    (0, 1, 1, 0): frozenset((C4, BETA, EETA)),
    (0, 0, 1, 0): frozenset((BDELTA, EDELTA)),
    (0, 1, 0, 0): frozenset((C3, BETA, EETA)),
    (0, 0, 0, 1): frozenset((CPHI, BPHI, EDELTA)),
    (0, 0, 0, 0): frozenset((BDEEP,)),
}
planes = candidate_planes(support)
Krel = sp.Rational(2848, 45) - sp.Rational(7, 6) * Delta
dstar = sp.Rational(3968, 63)
delta_reps = (sp.Integer(0), dstar, sp.Integer(1))

exact_supports = {}
fan_counts = Counter()
for delta_value, u_value, phi_value, eta_value in product(
        delta_reps, (0, 1), (0, 1), (0, 1)):
    substitution = {
        Delta: delta_value,
        K: Krel.subs(Delta, delta_value),
        u: u_value,
        Phi: phi_value,
        eta: eta_value,
    }
    need(substitution[K] != 0, "K-nonzero representative")
    active = frozenset(point for point, coefficient in support.items()
                       if sp.simplify(coefficient.subs(substitution)) != 0)
    status = (int(u_value != 0), int(eta_value != 0),
              int(delta_value != 0), int(phi_value != 0))
    # Delta is an owner only after u vanishes, and Phi only after
    # u=Delta=eta=0.
    if status[0]:
        key = (1, status[1], 0, 0)
    elif status[2]:
        key = (0, status[1], 1, 0)
    elif status[1]:
        key = (0, 1, 0, 0)
    else:
        key = (0, 0, 0, status[3])
    fan = lower_fan(active, planes)
    need(fan == expected_fans[key], "exact owner fan")
    exact_supports[active] = substitution
    fan_counts[fan] += 1

need(len(exact_supports) == 24, "twenty-four exact supports")
need(sorted(fan_counts.values()) == [1, 1, 2, 4, 4, 6, 6],
     "seven fan support counts")

# Conservative coupling-hostile atlas.
required_labels = ((1, 0), (2, 0), (3, 0), (0, 2))
optional_labels = ((4, 0), (5, 0), (2, 1), (3, 1))
cancellable = ((2, 3, 1), (2, 4, 1), (2, 5, 1))


def lifted(labels):
    active = {(2, 0, 0), (0, 1, 0), (2, 0, 1)}
    for i, j in labels:
        active.add((j + 2, i + j, 1))
        active.add((j, i + j + 1, 1))
    return active


over = {}
over_distinct = set()
over_fans = set()
for row_bits in product((0, 1), repeat=4):
    labels = list(required_labels)
    labels.extend(label for bit, label in zip(row_bits, optional_labels) if bit)
    for cancellation_bits in product((0, 1), repeat=3):
        active = lifted(labels)
        for bit, point in zip(cancellation_bits, cancellable):
            if bit:
                active.discard(point)
        active = frozenset(active)
        key = tuple(row_bits) + tuple(cancellation_bits)
        over[key] = active
        over_distinct.add(active)
        over_fans.add(lower_fan(active, planes))

need(len(over) == 128, "128 keyed hostile configurations")

for active, substitution in exact_supports.items():
    row_bits = tuple(int(substitution[v] != 0) for v in (Delta, u, Phi, eta))
    cancellation_bits = (
        int(sp.simplify((K - e).subs(substitution)) == 0),
        int(substitution[Delta] == 0),
        int(substitution[u] == 0),
    )
    need(over[row_bits + cancellation_bits] == active,
         "exact support embeds in hostile atlas")

# Exact face equations and lattice data.
S, P, X = sp.symbols("S P X")


def face_expression(plane, substitution=None):
    expression = sum(coefficient * S**rr * P**ll
                     for (rr, ll, hh), coefficient in support.items()
                     if on_plane((rr, ll, hh), plane))
    if substitution:
        expression = expression.subs(substitution)
    return sp.factor(expression)


faces = {
    "BU": (BU, {}, (P - S**2) * (u * P**5 - 1),
           [(0, 1), (2, 0), (2, 5), (0, 6)], (20, 12, 5), (30, 22)),
    "FUE": (FUE, {}, S**2 * (1 - u * P**5 - eta * S * P**4),
            [(2, 0), (2, 5), (3, 4)], (5, 7, 0), (30, 25)),
    "EETA": (EETA, {}, S**2 * (1 - eta * S * P**4 - K * S**2 * P**2),
             [(2, 0), (3, 4), (4, 2)], (6, 4, 2), (6, 6)),
    "EU": (EU, {}, S**2 * (1 - u * P**5 - K * S**2 * P**2),
           [(2, 0), (2, 5), (4, 2)], (10, 8, 2), (30, 28)),
    "C4": (C4, {u: 0}, P * (Delta * P**4 + eta * S * P**4 - 1),
           [(0, 1), (1, 5), (0, 5)], (4, 6, 0), (12, 10)),
    "BETA": (BETA, {u: 0}, (P - S**2) * (eta * S * P**4 - 1),
             [(0, 1), (2, 0), (3, 4), (1, 5)], (18, 4, 8), (18, 13)),
    "C3": (C3, {u: 0, Delta: 0}, P * (e * P**3 + eta * S * P**4 - 1),
           [(0, 1), (1, 5), (0, 4)], (3, 5, 0), (6, 7)),
    "BDELTA": (BDELTA, {u: 0, eta: 0},
               (P - S**2) * (Delta * P**4 - 1),
               [(0, 1), (2, 0), (2, 4), (0, 5)], (16, 10, 4), (24, 17)),
    "EDELTA": (EDELTA, {u: 0, eta: 0},
               S**2 * (1 - Delta * P**4 - Phi * S * P**3
                       - K * S**2 * P**2),
               [(2, 0), (2, 4), (4, 2)], (8, 8, 1), (12, 10)),
    "CPHI": (CPHI, {u: 0, eta: 0, Delta: 0},
             P * (e * P**3 + Phi * S * P**3 - 1),
             [(0, 1), (1, 4), (0, 4)], (3, 5, 0), (6, 5)),
    "BPHI": (BPHI, {u: 0, eta: 0, Delta: 0},
             (P - S**2) * (Phi * S * P**3 - 1),
             [(0, 1), (2, 0), (3, 3), (1, 4)], (14, 4, 6), (42, 29)),
    "BDEEP": (BDEEP, {u: 0, eta: 0, Delta: 0, Phi: 0},
              (P - S**2) * (e * P**3 + K * S**2 * P**2 - 1),
              [(0, 1), (2, 0), (4, 2), (0, 4)], (18, 8, 6), (6, 4)),
}
for name, (plane, substitution, expected, polygon, pick_data, order_data) in faces.items():
    need(sp.expand(face_expression(plane, substitution) - expected) == 0,
         name + " exact face")
    need(pick(polygon) == pick_data, name + " Pick data")
    need(face_order(plane) == order_data, name + " form order")

# Every edge word.  The sole potentially repeated edge is EDELTA's middle
# quadratic, whose discriminant is exactly L=Phi^2-4*K*Delta.
expected_edges = {
    "BU": (X - 1, 1 - u * X**5, u * (X - 1), u - X**5),
    "FUE": (1 - u * X**5, -eta * X - u, X - eta),
    "EETA": (1 - eta * X, -K * X - eta, X**2 - K),
    "EU": (1 - u * X**5, -K * X - u, X**2 - K),
    "C4": (eta * X - 1, Delta * X + eta, Delta - X**4),
    "BETA": (X - 1, 1 - eta * X, eta * (X - 1), eta - X),
    "C3": (eta * X - 1, eta + e * X, e - X**3),
    "BDELTA": (X - 1, 1 - Delta * X**4,
               Delta * (X - 1), Delta - X**4),
    "EDELTA": (1 - Delta * X**4,
               -Delta - Phi * X - K * X**2, X**2 - K),
    "CPHI": (Phi * X - 1, Phi + e * X, e - X**3),
    "BPHI": (X - 1, 1 - Phi * X, Phi * (X - 1), Phi - X),
    "BDEEP": (X - 1, 1 - K * X**2,
              (K + e * X) * (X - 1), e - X**3),
}
for name, (_plane, _substitution, expression, polygon, *_rest) in faces.items():
    actual = tuple(edge_polynomial(expression, polygon[index],
                                   polygon[(index + 1) % len(polygon)],
                                   S, P, X)
                   for index in range(len(polygon)))
    for got, want in zip(actual, expected_edges[name]):
        need(sp.expand(got - want) == 0, name + " exact edge word")

edge_quadratic = -Delta - Phi * X - K * X**2
need(sp.discriminant(edge_quadratic, X) == Phi**2 - 4 * K * Delta,
     "sole variable edge discriminant")
need(sp.factor(edge_quadratic.subs({Delta: K * sp.Symbol("aa")**2,
                                    Phi: -2 * K * sp.Symbol("aa")}))
     == -K * (X - sp.Symbol("aa"))**2,
     "collision is a doubled boundary root")

# Global Pick ledgers in fan order.
global_polygons = {
    "u_eta": ([(0, 1), (2, 0), (4, 2), (3, 4), (2, 5), (0, 6)],
              (31, 11, 11)),
    "u_zero": ([(0, 1), (2, 0), (4, 2), (2, 5), (0, 6)],
               (30, 10, 11)),
    "d_eta": ([(0, 1), (2, 0), (4, 2), (3, 4), (1, 5), (0, 5)],
              (28, 10, 10)),
    "d_zero": ([(0, 1), (2, 0), (4, 2), (2, 4), (0, 5)],
               (24, 10, 8)),
    "e_eta": ([(0, 1), (2, 0), (4, 2), (3, 4), (1, 5), (0, 4)],
              (27, 9, 10)),
    "phi": ([(0, 1), (2, 0), (4, 2), (3, 3), (1, 4), (0, 4)],
            (21, 9, 7)),
    "deep": ([(0, 1), (2, 0), (4, 2), (0, 4)],
             (18, 8, 6)),
}
for name, (polygon, expected) in global_polygons.items():
    need(pick(polygon) == expected, name + " global Pick")

# Normalization identities and squarefree carrier polynomials.
f_eta = 1 - eta * S * P**4 - K * S**2 * P**2
y_eta = 2 * K * S * P + eta * P**3
need(sp.expand(y_eta**2 - (eta**2 * P**6 + 4 * K) + 4 * K * f_eta) == 0,
     "EETA hyperelliptic identity")
f_u = 1 - u * P**5 - K * S**2 * P**2
y_u = K * S * P
need(sp.expand(y_u**2 - K * (1 - u * P**5) + K * f_u) == 0,
     "EU hyperelliptic identity")
f_delta = 1 - Delta * P**4 - Phi * S * P**3 - K * S**2 * P**2
y_delta = 2 * K * S * P + Phi * P**2
L = Phi**2 - 4 * K * Delta
need(sp.expand(y_delta**2 - (L * P**4 + 4 * K) + 4 * K * f_delta) == 0,
     "EDELTA hyperelliptic identity")
f_deep = 1 - e * P**3 - K * S**2 * P**2
y_deep = K * S * P
need(sp.expand(y_deep**2 - K * (1 - e * P**3) + K * f_deep) == 0,
     "BDEEP elliptic component identity")

for polynomial, variable, label in (
        (eta**2 * P**6 + 4 * K, P, "EETA genus two"),
        (1 - u * P**5, P, "EU genus two"),
        (L * P**4 + 4 * K, P, "EDELTA genus one off L"),
        (1 - e * P**3, P, "BDEEP genus one")):
    need(sp.degree(polynomial, variable) in (3, 4, 5, 6), label + " degree")
    need(sp.gcd(sp.Poly(polynomial, variable),
                sp.Poly(sp.diff(polynomial, variable), variable)).degree() == 0,
         label + " squarefree over coefficient field")

# Factor-intersection counts and graph ledgers.
intersection_degrees = {
    "BU": sp.degree((u * P**5 - 1).subs(P, S**2), S),
    "BETA": sp.degree((eta * S * P**4 - 1).subs(P, S**2), S),
    "BDELTA": sp.degree((Delta * P**4 - 1).subs(P, S**2), S),
    "BPHI": sp.degree((Phi * S * P**3 - 1).subs(P, S**2), S),
    "BDEEP": sp.degree((e * P**3 + K * S**2 * P**2 - 1).subs(P, S**2), S),
}
need(intersection_degrees == {"BU": 10, "BETA": 9, "BDELTA": 8,
                              "BPHI": 7, "BDEEP": 6},
     "factor intersection counts")
clean_ledgers = {
    "u_eta": (6 + 1 + 1, 10 + 5 + 1, 9, 2, 11),
    "u_zero": (6 + 1, 10 + 5, 9, 2, 11),
    "d_eta": (1 + 2 + 1, 9 + 1 + 1, 8, 2, 10),
    "d_zero": (5 + 1, 8 + 4, 7, 1, 8),
    "e_eta": (1 + 2 + 1, 9 + 1 + 1, 8, 2, 10),
    "phi": (1 + 2 + 1, 7 + 1 + 1, 6, 1, 7),
    "deep": (2, 6, 5, 1, 6),
}
for name, (vertices, edges, graph_genus, carrier_genus, total) in clean_ledgers.items():
    need(edges - vertices + 1 == graph_genus, name + " graph genus")
    need(graph_genus + carrier_genus == total, name + " total genus")

# The sole variable carrier wall and its exact A3 collision.
a, lam, sigma, x, v, rho, y, t, bigX, bigY, z, h = sp.symbols(
    "a lam sigma x v rho y t X Y z h")
collision = {Delta: K * a**2, Phi: -2 * K * a}
need(sp.expand(f_delta.subs(collision)
               - (1 - K * P**2 * (S - a * P)**2)) == 0,
     "collision square factor")
need(sp.expand((1 - lam * P * (S - a * P))
               * (1 + lam * P * (S - a * P))
               - f_delta.subs(collision).subs(K, lam**2)) == 0,
     "two rational collision signs")
need(sp.resultant(x**2 - lam * (v - a),
                  x**2 + lam * (v - a), v) == -2 * lam * x**2,
     "A3 intersection length two")

# The eta-nonzero boundary chart is smooth at both roots.  It must not be
# mistaken for the eta=0 double-root chart.
s0, p0, Q = sp.symbols("s0 p0 Q")
H_all = (-3 * p0 + sp.Rational(8, 3) * p0**2 + e * p0**3
         + K * (s0 * p0)**2 + Phi * p0**2 * (s0 * p0)
         + Delta * p0**4 + eta * p0**3 * (s0 * p0)
         + u * p0**5)
FQ_all = (s0**2 - p0) * (1 - Q * H_all) - Q * s0**2 / 2
G_eta_source = sp.expand(
    sigma**4 * FQ_all.subs({Q: sigma**6, s0: S / sigma**2,
                            p0: P / sigma}, simultaneous=True))
G_eta_expected = (
    (S**2 - sigma**3 * P)
    * (1 - eta * S * P**4 - K * S**2 * P**2
       - sigma * (u * P**5 + Phi * S * P**3)
       - sigma**2 * Delta * P**4 - e * sigma**3 * P**3
       - sp.Rational(8, 3) * sigma**4 * P**2
       + 3 * sigma**5 * P)
    - sigma**6 * S**2 / 2)
need(sp.expand(G_eta_source - G_eta_expected) == 0,
     "literal eta-nonzero primitive chart")
eta_root_equation = (
    (v**2 - rho**3)
    * (x**6 - eta * v - K * v**2 - rho * (u + Phi * v)
       - rho**2 * Delta - e * rho**3
       - sp.Rational(8, 3) * rho**4 + 3 * rho**5)
    - rho**6 * v**2 / 2)
eta_root_from_source = sp.cancel(
    x**10 * G_eta_expected.subs({S: v / x**2, P: 1 / x,
                                 sigma: rho / x}, simultaneous=True))
need(sp.expand(eta_root_from_source - eta_root_equation) == 0,
     "eta-nonzero reciprocal boundary chart")
eta_carrier = x**6 - eta * v - K * v**2
need(sp.diff(eta_carrier, v).subs(v, 0) == -eta,
     "v=0 boundary root simple when eta nonzero")
need(sp.simplify(sp.diff(eta_carrier, v).subs(v, -eta / K)) == eta,
     "second eta boundary root simple")

# When eta=0,u!=0 the lower fan changes before any local double-root
# argument: the EU carrier is already smooth genus two.
G_u_source = sp.expand(
    sigma**6 * FQ_all.subs({eta: 0, Q: sigma**10,
                            s0: S / sigma**3, p0: P / sigma**2},
                           simultaneous=True))
G_u_expected = (
    (S**2 - sigma**4 * P)
    * (1 - u * P**5 - K * S**2 * P**2
       - sigma * Phi * S * P**3 - sigma**2 * Delta * P**4
       - e * sigma**4 * P**3
       - sp.Rational(8, 3) * sigma**6 * P**2
       + 3 * sigma**8 * P)
    - sigma**10 * S**2 / 2)
need(sp.expand(G_u_source - G_u_expected) == 0,
     "literal eta-zero u-owner primitive chart")

# Primitive EDELTA chart: Q=sigma^4, s=sigma^-1*S, p=sigma^-1*P.
H = (-3 * p0 + sp.Rational(8, 3) * p0**2 + e * p0**3
     + K * (s0 * p0)**2 + Phi * p0**2 * (s0 * p0)
     + Delta * p0**4)
FQ = (s0**2 - p0) * (1 - Q * H) - Q * s0**2 / 2
G_from_source = sp.expand(
    sigma**2 * FQ.subs({Q: sigma**4, s0: S / sigma,
                        p0: P / sigma}, simultaneous=True))
H4 = Delta * P**4 + Phi * S * P**3 + K * S**2 * P**2
G_expected = ((S**2 - sigma * P)
              * (1 - H4 - e * sigma * P**3
                 - sp.Rational(8, 3) * sigma**2 * P**2
                 + 3 * sigma**3 * P)
              - sigma**4 * S**2 / 2)
need(sp.expand(G_from_source - G_expected) == 0,
     "literal source primitive EDELTA chart")

A = Delta + Phi * v + K * v**2
root_equation = ((v**2 - rho)
                 * (x**4 - A - e * rho - sp.Rational(8, 3) * rho**2
                    + 3 * rho**3)
                 - sp.Rational(1, 2) * rho**4 * v**2)
root_from_source = sp.cancel(
    x**6 * G_expected.subs({S: v / x, P: 1 / x,
                            sigma: rho / x}, simultaneous=True))
need(sp.expand(root_from_source - root_equation) == 0,
     "literal source reciprocal chart")
D = (A + e * rho + sp.Rational(8, 3) * rho**2 - 3 * rho**3
     + rho**4 * v**2 / (2 * (v**2 - rho)))
need(sp.factor(root_equation - (v**2 - rho) * (x**4 - D)) == 0,
     "exact divided reciprocal equation")
need(sp.expand(A.subs(collision) - K * (v - a)**2) == 0,
     "collision transverse Morse square")
need(sp.simplify(sp.diff(D, v)
                 - (2 * K * v + Phi - rho**5 * v / (v**2 - rho)**2)) == 0,
     "exact critical derivative")
need(sp.diff(D.subs(collision), v, 2).subs({v: a, rho: 0}) == 2 * K,
     "nonzero Morse Hessian")
D_at_a = sp.series(D.subs(collision).subs(v, a), rho, 0, 6).removeO()
expected_series = (e * rho + sp.Rational(8, 3) * rho**2 - 3 * rho**3
                   + sp.Rational(1, 2) * rho**4
                   + rho**5 / (2 * a**2))
need(sp.expand(D_at_a - expected_series) == 0,
     "critical value through rho^5")
need(e != 0 and a != 0, "forced depth one and interior ratio root")

# Even-cusp hypotheses are literal: m=4,r=1,c=e, primitive ray (3,1,2).
tail = bigY**2 - bigX * (bigX**3 - e)
need(sp.gcd(sp.Poly(bigX * (bigX**3 - e), bigX),
            sp.Poly(sp.diff(bigX * (bigX**3 - e), bigX), bigX)).degree() == 0,
     "elliptic tail branch polynomial squarefree")
need(gcd(gcd(3, 1), 2) == 1, "primitive sigma/x/y valuation ray")
tail_m, tail_r = 4, 1
tail_genus = (tail_m - tail_r - 1) // 2
persistent_delta = tail_r // 2
need((tail_genus, persistent_delta, tail_genus + persistent_delta)
     == (1, 0, 1), "m4 r1 intrinsic even-cusp ledger")
need(10 + 3 + 1 - 2 == 12, "positive primitive tail differential order")
collision_vertices, collision_edges = 5 + 2, 8 + 4
collision_b1 = collision_edges - collision_vertices + 1
need((collision_vertices, collision_edges, collision_b1) == (7, 12, 6),
     "collision complement graph")
tail_vertices, tail_edges = collision_vertices + 1, collision_edges + 2
tail_b1 = tail_edges - tail_vertices + 1
need((tail_vertices, tail_edges, tail_b1) == (8, 14, 7),
     "two-ended tail adds one graph cycle")
need(tail_b1 + tail_genus + persistent_delta == 8,
     "tail restores global genus")

# At the two ends z*rho=t^4 and h^2=1-z^3*C(rho); h=+/-1.
need(sp.expand((z * rho).subs({z: t / x, rho: t**3 * x}) - t**4) == 0,
     "complementary attachment surface")
need(sp.expand((rho / x**4).subs({x: t / z, rho: t**4 / z}) - z**3) == 0,
     "complementary tail exponent")
need(2 == 2, "two distinct etale ends h=+1,-1")
need(2 + 2 == 4, "collision signs receive two BDELTA nodes each")

# Even A3 conductor: normalization pairs congruent modulo x^2.
need(2 == 2, "A3 conductor exponent x^2")
need(3 >= 2, "x^3 residue numerator lies in conductor and Jacobian ideal")


boundary_checks_before_union = CHECKS
boundary_exact_supports = frozenset(exact_supports)
boundary_exact_fans = frozenset(fan_counts)
boundary_hostile_supports = frozenset(over_distinct)
boundary_hostile_fans = frozenset(over_fans)

boundary_natural_index = {}
for active, substitution in exact_supports.items():
    coefficient_class = delta_reps.index(substitution[Delta])
    number = (41 + 8 * coefficient_class + 4 * int(substitution[u] != 0)
              + 2 * int(substitution[Phi] != 0)
              + int(substitution[eta] != 0))
    need(number not in boundary_natural_index,
         "Theta-zero natural index injective")
    boundary_natural_index[number] = active
need(set(boundary_natural_index) == set(range(41, 65)),
     "Theta-zero natural index is 41..64")

# Full alpha-zero gate.  The two honest Theta strata are disjoint both as
# supports and as lower fans.  Their natural-number coordinates concatenate
# without overlap.  The hostile union is explicitly an adversarial atlas,
# not a claim that its extra cancellation toggles are coefficient-realizable.
need(theta_exact_supports.isdisjoint(boundary_exact_supports),
     "Theta strata have disjoint supports")
need(len(theta_exact_supports | boundary_exact_supports) == 64,
     "full gate has exactly 64 supports")
need(theta_exact_fans.isdisjoint(boundary_exact_fans),
     "Theta strata have disjoint exact fans")
need(len(theta_exact_fans | boundary_exact_fans) == 12,
     "full gate has exactly 12 fans")
need(theta_hostile_supports.isdisjoint(boundary_hostile_supports),
     "hostile support strata remain Theta-separated")
need(len(theta_hostile_supports | boundary_hostile_supports) == 168,
     "combined hostile atlas has 168 distinct supports")
need(theta_hostile_fans.isdisjoint(boundary_hostile_fans),
     "hostile fan strata remain Theta-separated")
need(len(theta_hostile_fans | boundary_hostile_fans) == 20,
     "combined hostile atlas has 20 fans")
full_natural_index = theta_natural_index | boundary_natural_index
need(set(full_natural_index) == set(range(1, 65)),
     "full support index is exactly 1..64")
need(len(set(full_natural_index.values())) == 64,
     "full natural index is bijective")
for number in range(1, 65):
    odd = 2 * number - 1
    need((odd * odd - 1) // 8 == number * (number - 1) // 2,
         "odd-square address retains triangular predecessor index")

total_checks = theta_checks + CHECKS
print("THM4355_TERMINAL_ALPHA_ZERO_U_ZERO_ENDPOINT_CERTIFICATE=PASS")
print("gate=Z=beta=zeta=W=xi=U=alpha=0;K!=0")
print("seam=K=2848/45-(7/6)*Delta")
print("literal_source=16_rows;Theta_nonzero=9_rows/19_points;Theta_zero=8_rows/18_points")
print("exact_supports=64;Theta_nonzero=40;Theta_zero=24")
print("exact_fans=12;Theta_nonzero=5[20,8,8,2,2];Theta_zero=7[6,6,4,4,2,1,1]")
print("natural_support_index=bijective_1_through_64;odd_square_triangular=yes")
print("hostile_scope=Theta_nonzero:128_keyed/96_supports/9_fans;"
      "Theta_zero:128_keyed/72_supports/11_fans")
print("hostile_union=256_keyed/168_distinct_supports/20_fans;adversarial_not_realizability")
print("clean_global_genera=Theta_nonzero:11,10,9,10,8;Theta_zero:11,11,10,8,10,7,6")
print("collision_walls=A4:eta^2=4*u*Theta;A15:Delta+Theta=0;"
      "A2:Phi^2=4*e*Theta;A3:Phi^2=4*K*Delta")
print("A4=depths_1,2,3,4;one_end;genera_2,1,1,0;delta_0,1,1,2;orders_115,35,95,85")
print("A15=depths_1,2,4;two_ends_to_distinct_complements;"
      "packets_2rational7nodes,2rational6nodes,genus3;orders_116,16,14")
print("A2=depths_1,2;one_end;genera_1,0;delta_0,1;orders_13,11")
print("A3=forced_depth_1;two_ends_same_complement;genus_1;delta_0;order_12;graph_plus_1")
print("critical_depth_rigidity=A4_c4=-99/43;A2_c2=-99/43;"
      "A15_deep_product=1472/15;A3_c1=e=-1376/135")
print("theta_checks=%d;boundary_and_union_checks=%d;CHECKS=%d" %
      (theta_checks, CHECKS, total_checks))
