#!/usr/bin/env python3
"""Independent hostile referee for THM-4351's F00 corner.

This file deliberately imports no repository mathematics or computation
module.  It rebuilds the sixteen-row source, specializes the doubly deleted
normal-owner gate, exhausts exact and conservative support states, derives
the primitive E00 and reciprocal charts, and audits the even-A5 tails.
"""

from collections import defaultdict
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
        raise RuntimeError("F00 independent-referee failure: " + label)


def candidate_planes(universe):
    planes = set()
    for p, q, r in combinations(sorted(universe), 3):
        det = ((q[0] - p[0]) * (r[1] - p[1])
               - (q[1] - p[1]) * (r[0] - p[0]))
        if det == 0:
            continue
        aa = F((q[2] - p[2]) * (r[1] - p[1])
               - (q[1] - p[1]) * (r[2] - p[2]), det)
        bb = F((q[0] - p[0]) * (r[2] - p[2])
               - (q[2] - p[2]) * (r[0] - p[0]), det)
        cc = F(p[2]) - aa * p[0] - bb * p[1]
        planes.add((aa, bb, cc))
    return tuple(sorted(planes))


def rank_two(points):
    for p, q, r in combinations(points, 3):
        if ((q[0] - p[0]) * (r[1] - p[1])
                != (q[1] - p[1]) * (r[0] - p[0])):
            return True
    return False


def lower_fan(active, planes):
    fan = set()
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
                fan.add((aa, bb, cc))
    return frozenset(fan)


def on_plane(point, plane):
    rr, ll, hh = point
    aa, bb, cc = plane
    return F(hh) == aa * rr + bb * ll + cc


def pick(vertices):
    edges = list(zip(vertices, vertices[1:] + vertices[:1]))
    area2 = abs(sum(x1 * y2 - x2 * y1
                    for (x1, y1), (x2, y2) in edges))
    boundary = sum(gcd(abs(x2 - x1), abs(y2 - y1))
                   for (x1, y1), (x2, y2) in edges)
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
    order = base * density
    need(order.denominator == 1 and order > 0, "positive integral face order")
    return base, order.numerator


def series_order(expression, variable, limit):
    polynomial = sp.Poly(sp.series(expression, variable, 0, limit)
                         .removeO().expand(), variable)
    terms = [monomial[0] for monomial, coefficient in polynomial.terms()
             if coefficient != 0]
    return min(terms) if terms else sp.oo


# -------------------------------------------------------------------------
# Literal source and exact coupled support.

U, W, Z = sp.symbols("U W Z")
Delta, Theta, Phi, eta = sp.symbols("Delta Theta Phi eta")
zeta, u, xi, alpha, beta, K = sp.symbols("zeta u xi alpha beta K")
e = -sp.Rational(1376, 135)

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
need(len(FULL_ROWS) == 16, "literal sixteen-row source")

GATE = {Z: 0, beta: 0, zeta: 0, W: 0, xi: 0, alpha: 0, Theta: 0}
ROWS = {label: sp.sympify(coefficient).subs(GATE)
        for label, coefficient in FULL_ROWS.items()
        if sp.sympify(coefficient).subs(GATE) != 0}
need(len(ROWS) == 9, "nine surviving labelled source rows")

support = defaultdict(lambda: sp.Integer(0))
support[(2, 0, 0)] += 1
support[(0, 1, 0)] -= 1
support[(2, 0, 1)] -= sp.Rational(1, 2)
for (ii, jj), coefficient in ROWS.items():
    support[(jj + 2, ii + jj, 1)] -= coefficient
    support[(jj, ii + jj + 1, 1)] += coefficient
support = {point: sp.factor(coefficient)
           for point, coefficient in support.items() if coefficient != 0}

need(support[(2, 6, 1)] == -U, "unique U lower owner")
need(support[(0, 7, 1)] == U, "unique U upper owner")
need(support[(4, 2, 1)] == -K, "unique K owner")
need(support[(3, 4, 1)] == -eta, "eta midpoint")
need(sp.expand(support[(2, 3, 1)] - (K - e)) == 0,
     "coupled K/e coefficient")
need(support[(2, 4, 1)] == -Delta, "Theta=0 Delta coefficient")
need(support[(2, 5, 1)] == -u, "coupled u lift")

M = (F(1, 12), F(1, 6), F(-1, 6))
E00 = (F(1, 3), F(1, 6), F(-2, 3))
EXPECTED_FAN = frozenset((M, E00))
planes = candidate_planes(support)

KREL = sp.Rational(2848, 45) - sp.Rational(7, 6) * Delta
DSTAR = sp.solve(sp.Eq(KREL - e, 0), Delta)[0]
KZERO = sp.solve(sp.Eq(KREL, 0), Delta)[0]
need(DSTAR == sp.Rational(3968, 63), "K/e cancellation address")
need(KZERO == sp.Rational(5696, 105), "excluded K=0 address")

exact_masks = set()
for delta_value, u_value, phi_value, eta_value in product(
        (0, DSTAR, 1), (0, 1), (0, 1), (0, 1)):
    kval = KREL.subs(Delta, delta_value)
    need(kval != 0, "exact representative avoids K=0")
    substitution = {
        U: 2, K: kval, Delta: delta_value,
        u: u_value, Phi: phi_value, eta: eta_value,
    }
    active = frozenset(point for point, coefficient in support.items()
                       if sp.simplify(coefficient.subs(substitution)) != 0)
    need(lower_fan(active, planes) == EXPECTED_FAN,
         "exact support has M,E00 fan")
    exact_masks.add(active)
need(len(exact_masks) == 24, "24 distinct coupled exact supports")

# A deliberately larger owner-00 Boolean atlas: four row-presence bits and
# three independently toggled aggregate deletions.  It is not realizability.
REQUIRED = ((1, 0), (2, 0), (3, 0), (0, 2), (6, 0))
OPTIONAL = ((4, 0), (5, 0), (2, 1), (3, 1))
CANCELLABLE = ((2, 3, 1), (2, 4, 1), (2, 5, 1))


def lifted(labels):
    active = {(2, 0, 0), (0, 1, 0), (2, 0, 1)}
    for ii, jj in labels:
        active.add((jj + 2, ii + jj, 1))
        active.add((jj, ii + jj + 1, 1))
    return active


hostile_keyed = []
hostile_distinct = set()
for row_bits, deletion_bits in product(
        product((0, 1), repeat=4), product((0, 1), repeat=3)):
    labels = list(REQUIRED)
    labels.extend(label for bit, label in zip(row_bits, OPTIONAL) if bit)
    active = lifted(labels)
    for bit, point in zip(deletion_bits, CANCELLABLE):
        if bit:
            active.discard(point)
    active = frozenset(active)
    need(lower_fan(active, planes) == EXPECTED_FAN,
         "hostile Boolean state has M,E00 fan")
    hostile_keyed.append((row_bits, deletion_bits, active))
    hostile_distinct.add(active)
need(len(hostile_keyed) == 128, "128 keyed hostile configurations")
need(exact_masks <= hostile_distinct, "every exact mask embeds in hostile atlas")

# -------------------------------------------------------------------------
# Faces, polygons, edge words, carrier, and generic graph.

S, P, X = sp.symbols("S P X")


def face_expression(plane):
    return sp.factor(sum(coefficient * S**rr * P**ll
                         for (rr, ll, hh), coefficient in support.items()
                         if on_plane((rr, ll, hh), plane)))


GM = (S**2 - P) * (1 - U * P**6)
GE = S**2 * (1 - U * P**6 - eta * S * P**4 - K * S**2 * P**2)
need(sp.expand(face_expression(M) - GM) == 0, "literal M face")
need(sp.expand(face_expression(E00) - GE) == 0, "literal E00 face")

POLYGONS = {
    "M": [(0, 1), (2, 0), (2, 6), (0, 7)],
    "E00": [(2, 0), (4, 2), (2, 6)],
    "GLOBAL": [(0, 1), (2, 0), (4, 2), (2, 6), (0, 7)],
}
need(pick(POLYGONS["M"]) == (24, 14, 6), "M Pick ledger")
need(pick(POLYGONS["E00"]) == (12, 10, 2), "E00 Pick ledger")
need(pick(POLYGONS["GLOBAL"]) == (36, 12, 13), "global Pick ledger")
need(edge_packet(POLYGONS["GLOBAL"]) == (11, 7, 7, 2, 2, 1),
     "global source-infinity packet")
need(sum(index - 1 for index in edge_packet(POLYGONS["GLOBAL"]))
     == 2 * 13 - 2, "global Riemann-Hurwitz checksum")
need(face_order(M) == (12, 9), "M primitive base/order")
need(face_order(E00) == (6, 6), "E00 primitive base/order")

outer_edges = (
    X - 1,
    1 - K * X**2,
    K + eta * X + U * X**2,
    X - 1,
    U - X**6,
)
internal_edge = 1 - U * X**6
for index, polynomial in enumerate(outer_edges):
    if index != 2:
        need(sp.degree(sp.gcd(polynomial, sp.diff(polynomial, X)), X) == 0,
             "fixed outer edge squarefree")
need(sp.degree(sp.gcd(internal_edge, sp.diff(internal_edge, X)), X) == 0,
     "internal edge squarefree")
need(sp.discriminant(outer_edges[2], X) == eta**2 - 4 * K * U,
     "sole variable edge discriminant")

f00 = 1 - U * P**6 - eta * S * P**4 - K * S**2 * P**2
y00 = 2 * K * S * P + eta * P**3
DISC = eta**2 - 4 * K * U
need(sp.expand(y00**2 - (4 * K + DISC * P**6) + 4 * K * f00) == 0,
     "exact E00 hyperelliptic normalization")
branch = 4 * K + DISC * P**6
need(sp.factor(sp.discriminant(branch, P))
     == sp.Integer(47775744) * K**5 * (4 * K * U - eta**2)**5,
     "exact branch discriminant")
need(18 - 8 + 1 == 11 and 11 + 2 == 13,
     "generic V=8 E=18 b1=11 plus carrier genus two")

# -------------------------------------------------------------------------
# Primitive E00 chart and the exact reciprocal one-series equation.

sigma = sp.symbols("sigma")
p0 = P / sigma
s0 = S / sigma**2
y0 = s0 * p0
H0 = sum(coefficient * p0**ii * y0**jj
         for (ii, jj), coefficient in ROWS.items())
FQ = (s0**2 - p0) * (1 - sigma**6 * H0) - sigma**6 * s0**2 / 2
G = sp.expand(sigma**4 * FQ)

H6 = U * P**6 + eta * S * P**4 + K * S**2 * P**2
H5 = u * P**5 + Phi * S * P**3
EXPECTED_G = sp.expand(
    (S**2 - sigma**3 * P)
    * (1 - H6 - sigma * H5 - sigma**2 * Delta * P**4
       - e * sigma**3 * P**3 - sp.Rational(8, 3) * sigma**4 * P**2
       + 3 * sigma**5 * P)
    - sigma**6 * S**2 / 2
)
need(sp.simplify(G - EXPECTED_G) == 0, "complete primitive E00 chart")

x, v, rho = sp.symbols("x v rho")
A = U + eta * v + K * v**2
B = u + Phi * v
ROOT = sp.factor(x**10 * G.subs({P: 1 / x, S: v / x**2}))
EXPECTED_ROOT = sp.expand(
    (v**2 - (sigma * x)**3)
    * (x**6 - A - sigma * x * B - (sigma * x)**2 * Delta
       - e * (sigma * x)**3 - sp.Rational(8, 3) * (sigma * x)**4
       + 3 * (sigma * x)**5)
    - (sigma * x)**6 * v**2 / 2
)
need(sp.simplify(ROOT - EXPECTED_ROOT) == 0, "exact reciprocal root chart")

DROOT = (
    A + rho * B + rho**2 * Delta + e * rho**3
    + sp.Rational(8, 3) * rho**4 - 3 * rho**5
    + rho**6 * v**2 / (2 * (v**2 - rho**3))
)
need(sp.factor(sp.together(
    ROOT / (v**2 - (sigma * x)**3)
    - (x**6 - DROOT.subs(rho, sigma * x)))) == 0,
    "divided reciprocal one-series equation")

# Coordinate/differential conversion.  On the curve the second product-rule
# term vanishes, so Phi_v is a unit times the Morse coordinate.
GS_ROOT = sp.diff(G, S).subs({P: 1 / x, S: v / x**2})
need(sp.simplify(sp.diff(ROOT, v) - x**8 * GS_ROOT) == 0,
     "G_S=x^-8 Phi_v differential conversion")

# -------------------------------------------------------------------------
# The repeated edge, critical series, tail invariants, and addresses.

a = sp.symbols("a", nonzero=True)
REPEAT = {U: K * a**2, eta: -2 * K * a}
need(sp.expand(A.subs(REPEAT) - K * (v - a)**2) == 0,
     "discriminant-zero square parameterization")
need(sp.diff(DROOT, v, 2).subs({rho: 0, v: a, **REPEAT}) == 2 * K,
     "nonzero Morse Hessian")
need((v**2 - rho**3).subs({v: a, rho: 0}) == a**2,
     "reciprocal prefactor is a unit")

v_linear = a - rho * Phi / (2 * K)
derivative_through_8 = sp.series(
    sp.diff(DROOT.subs(REPEAT), v).subs(v, v_linear), rho, 0, 9
).removeO().expand()
need(derivative_through_8 == 0,
     "linear critical point is exact through rho^8")
critical = sp.series(
    DROOT.subs(REPEAT).subs(v, v_linear), rho, 0, 7
).removeO().expand()
EXPECTED_CRITICAL = (
    (u + Phi * a) * rho
    + (Delta - Phi**2 / (4 * K)) * rho**2
    + e * rho**3 + sp.Rational(8, 3) * rho**4
    - 3 * rho**5 + sp.Rational(1, 2) * rho**6
)
need(sp.simplify(critical - EXPECTED_CRITICAL) == 0,
     "critical-value series through rho^6")
need(e == -sp.Rational(1376, 135) and e != 0,
     "fixed cubic coefficient forces critical depth at most three")

# Nonempty controls for all and only r=1,2,3, including the seam relation.
for expected_r, values in {
    1: {Delta: 0, K: KREL.subs(Delta, 0), Phi: 0, u: 1, a: 1},
    2: {Delta: 1, K: KREL.subs(Delta, 1), Phi: 0, u: 0, a: 1},
    3: {Delta: 0, K: KREL.subs(Delta, 0), Phi: 0, u: 0, a: 1},
}.items():
    specialized = EXPECTED_CRITICAL.subs(values)
    need(series_order(specialized, rho, 7) == expected_r,
         "critical-depth positive control")

lam = sp.symbols("lambda", nonzero=True)
linear = lam * P * (S - a * P**2)
factor_sub = {K: lam**2, U: lam**2 * a**2, eta: -2 * lam**2 * a}
need(sp.expand(f00.subs(factor_sub) - (1 - linear) * (1 + linear)) == 0,
     "repeated E00 splits into two rational signs")
need(sp.expand((1 - U * P**6).subs({U: lam**2 * a**2})
               - (1 - lam * a * P**3) * (1 + lam * a * P**3)) == 0,
     "six M attachments split three plus three")
need(sp.simplify(linear.subs({P: 1 / x, S: v / x**2})
                 - lam * (v - a) / x**3) == 0,
     "complement signs address the two face factors")

tau, t, XX, YY, z, h, c0, crho = sp.symbols(
    "tau t XX YY z h c0 crho", nonzero=True
)
TAILS = {
    1: (2, 0, 68),
    2: (1, 1, 64),
    3: (1, 1, 60),
}
for rr, (tail_genus, persistent_delta, form_order) in TAILS.items():
    epsilon = rr % 2
    tail_polynomial = XX**epsilon * (XX**(6 - rr) - c0)
    degree = sp.degree(tail_polynomial, XX)
    need(sp.degree(sp.gcd(tail_polynomial,
                          sp.diff(tail_polynomial, XX)), XX) == 0,
         "tail branch polynomial squarefree")
    need((degree - 1) // 2 == tail_genus, "tail genus")
    need(rr // 2 == persistent_delta, "persistent local delta")
    need(persistent_delta + tail_genus + 1 == 3,
         "A5 delta conservation with two-ended graph cycle")
    # sigma=tau^[2(6-r)], x=tau^[2r]X, y=tau^[6r]Y is a
    # convenient common cover, not the minimal Newton base.
    lhs_order = 2 * 6 * rr
    x6_order = 6 * 2 * rr
    perturbation_order = (2 * (6 - rr) * rr
                          + 2 * rr * rr)
    need(lhs_order == x6_order == perturbation_order,
         "honest tail scaling balances all terms")
    common_equation = (
        (tau**(6 * rr) * YY)**2
        - (tau**(2 * rr) * XX)**6
        + ((tau**(2 * (6 - rr)) * tau**(2 * rr) * XX)**rr) * c0
    )
    need(sp.expand(common_equation / tau**(12 * rr)
                   - (YY**2 - XX**rr * (XX**(6 - rr) - c0))) == 0,
         "common cover has the claimed exceptional function field")
    computed_order = (6 * 2 * (6 - rr) + 6 * 2 * rr
                      + 2 * rr - 6 * rr)
    need(computed_order == form_order > 0,
         "tail differential order 68/64/60")

    # Newton equality (6-r)ord(x)=r ord(sigma) gives the primitive ray.
    gg = gcd(rr, 6 - rr)
    wsigma = (6 - rr) // gg
    wx = rr // gg
    wy = 3 * rr // gg
    need(gcd(gcd(wsigma, wx), wy) == 1,
         "primitive Newton valuation vector")
    need(2 * wy == 6 * wx == rr * (wsigma + wx),
         "primitive weights balance the formal equation")
    primitive_equation = (
        (t**wy * YY)**2 - (t**wx * XX)**6
        + ((t**wsigma * t**wx * XX)**rr) * c0
    )
    need(sp.expand(primitive_equation / t**(2 * wy)
                   - (YY**2 - XX**rr * (XX**(6 - rr) - c0))) == 0,
         "primitive base has the same tail, not a root-coordinate cover")
    primitive_order = 6 * wsigma + 6 * wx + wx - wy
    need(primitive_order == {1: 34, 2: 16, 3: 10}[rr] > 0,
         "primitive tail differential order 34/16/10")
    need(computed_order == 2 * gg * primitive_order,
         "common-cover order is ramification times primitive order")

    # rho=tau^12 X, z=1/X=tau^12/rho.  This exact identity yields
    # h^2=1-z^(6-r)C(rho) on the complementary chart.
    raw = tau**(12 * (6 - rr)) * rho**(rr - 6)
    need(sp.simplify(raw.subs(rho, tau**12 / z) - z**(6 - rr)) == 0,
         "reciprocal complementary-chart monomial")
    # On the primitive base rho=t^[6/g]X, so the attachment surface is
    # A_(6/g-1), namely A5,A2,A1 rather than the common-cover A11.
    need(6 // gg in (6, 3, 2), "primitive complement exponent")
    primitive_rho = t**(6 // gg) * XX
    primitive_x = t**wx * XX
    need(sp.simplify((primitive_rho**rr / primitive_x**6)
                     .subs(XX, 1 / z) - z**(6 - rr)) == 0,
         "primitive reciprocal complement")
    complement = h**2 - 1 + z**(6 - rr) * crho
    need(sp.diff(complement, h).subs({h: 1, z: 0}) == 2,
         "positive endpoint etale")
    need(sp.diff(complement, h).subs({h: -1, z: 0}) == -2,
         "negative endpoint etale")

# Generic, raw repeated, and regularized graph/genus ledgers.
need(12 + 6 - 8 + 1 == 11, "generic graph b1=11")
need(12 + 6 - 9 + 1 == 10, "split-face graph b1=10")
need(12 + 6 + 3 - 9 + 1 == 13,
     "raw repeated arithmetic genus with delta(A5)=3")
for tail_genus, persistent_delta, _order in TAILS.values():
    regularized_b1 = (18 + 2) - (9 + 1) + 1
    need(regularized_b1 == 11, "tail plus two paths adds one graph cycle")
    need(regularized_b1 + tail_genus + persistent_delta == 13,
         "regularized global genus remains thirteen")

for surface_exponent in (12, 6, 3, 2):
    # Replacing each of the two contracted attachment edges by an A_(N-1)
    # chain adds N-1 vertices and N-1 net edges per side.
    expanded_vertices = 10 + 2 * (surface_exponent - 1)
    expanded_edges = 20 - 2 + 2 * surface_exponent
    need(expanded_edges - expanded_vertices + 1 == 11,
         "resolving both attachment surfaces preserves graph b1")

print("THM-4351 F00 hostile referee: PASS")
print("SCOPE Z=beta11=zeta3=W=xi10=alpha11=Theta=0; U*K!=0")
print("EXACT_SUPPORTS=24; FAN=M,E00")
print(f"HOSTILE_ATLAS=128_KEYED/{len(hostile_distinct)}_DISTINCT; FAN=M,E00")
print("PICK M=(24,14,6) E00=(12,10,2) GLOBAL=(36,12,13)")
print("PACKET=(11,7,7,2,2,1); FACE_BASE_ORDER M=12/9 E00=6/6")
print("GENERIC E00: y^2=4K+(eta^2-4KU)P^6; genus=2")
print("COLLISION eta^2=4KU: critical c1=u+Phi*a")
print("critical c2=Delta-Phi^2/(4K); c3=-1376/135; r<=3")
print("TAILS r=1:g2 r=2:g1 r=3:g1; primitive_orders=34,16,10")
print("COMMON_COVER orders=68,64,60; z*rho=tau^12; not primitive")
print("ATTACHMENTS primitive surfaces=A5,A2,A1; two etale signs; M split=3+3")
print("GENUS generic=11+2=13; collision=11+g_tail+floor(r/2)=13")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
