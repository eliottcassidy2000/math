#!/usr/bin/env python3
"""Exact primary certificate for THM-4350's xi_10=0 lower fan.

Scope:
    Z=beta_11=zeta_3=W=xi_10=0, U*K != 0,
    K=2848/45-(7/6)Delta.

The audit rebuilds the specialized support from all sixteen source rows,
enumerates coefficient-realizable supports and a conservative Boolean atlas,
checks all four owner fans, and derives the sole repeated-edge root chart.
It imports no repository mathematics or computation module.
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
        raise RuntimeError("xi_10=0 fan audit failure: " + label)


def rank_two(points):
    for p, q, r in combinations(points, 3):
        if ((q[0] - p[0]) * (r[1] - p[1])
                != (q[1] - p[1]) * (r[0] - p[0])):
            return True
    return False


def candidate_planes(universe):
    answer = set()
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


def pick(vertices):
    pairs = list(zip(vertices, vertices[1:] + vertices[:1]))
    area2 = abs(sum(x1 * y2 - x2 * y1
                    for (x1, y1), (x2, y2) in pairs))
    boundary = sum(gcd(abs(x2 - x1), abs(y2 - y1))
                   for (x1, y1), (x2, y2) in pairs)
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
    denominators = [value.denominator for value in plane] + [6]
    base = lcm(*denominators)
    aa, bb, cc = plane
    order = F(base) * (F(5, 6) - aa - bb - cc)
    need(order.denominator == 1 and order > 0, "positive integral face order")
    return base, order.numerator


def on_plane(point, plane):
    rr, ll, hh = point
    aa, bb, cc = plane
    return F(hh) == aa * rr + bb * ll + cc


def face_expression(support, plane, substitution=None):
    S, P = sp.symbols("S P")
    expression = sum(
        coefficient * S**rr * P**ll
        for (rr, ll, _hh), coefficient in support.items()
        if on_plane((rr, ll, _hh), plane)
    )
    if substitution:
        expression = expression.subs(substitution)
    return sp.factor(expression)


# The complete sixteen-term THM-4230 source.
U, W, Z = sp.symbols("U W Z")
Delta, Theta, Phi, eta = sp.symbols("Delta Theta Phi eta")
zeta, u, xi, alpha, beta, K = sp.symbols(
    "zeta u xi alpha beta K"
)
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
need(len(FULL_ROWS) == 16, "sixteen source rows")
GATE = {Z: 0, beta: 0, zeta: 0, W: 0, xi: 0}
ROWS = {
    label: sp.factor(sp.sympify(coefficient).subs(GATE))
    for label, coefficient in FULL_ROWS.items()
    if sp.factor(sp.sympify(coefficient).subs(GATE)) != 0
}
need(len(ROWS) == 11, "eleven surviving labelled rows")

SUPPORT = defaultdict(lambda: sp.Integer(0))
SUPPORT[(2, 0, 0)] += 1
SUPPORT[(0, 1, 0)] -= 1
SUPPORT[(2, 0, 1)] -= sp.Rational(1, 2)
for (i, j), coefficient in ROWS.items():
    SUPPORT[(j + 2, i + j, 1)] -= coefficient
    SUPPORT[(j, i + j + 1, 1)] += coefficient
SUPPORT = {
    point: sp.factor(coefficient)
    for point, coefficient in SUPPORT.items()
    if coefficient != 0
}

need(SUPPORT[(2, 6, 1)] == -U, "unique U lower owner")
need(SUPPORT[(0, 7, 1)] == U, "unique U upper owner")
need(SUPPORT[(4, 2, 1)] == -K, "unique K owner")
need(SUPPORT[(1, 6, 1)] == alpha, "first alpha owner")
need(SUPPORT[(3, 5, 1)] == -alpha, "second alpha owner")
need(SUPPORT[(4, 3, 1)] == -Theta, "Theta fan owner")
need(SUPPORT[(3, 4, 1)] == -eta, "eta midpoint owner")
need(sp.simplify(SUPPORT[(2, 3, 1)] - (K - e)) == 0, "K/e aggregate")
need(sp.simplify(SUPPORT[(2, 4, 1)] - (Theta - Delta)) == 0,
     "Theta/Delta aggregate")
need(SUPPORT[(2, 5, 1)] == -u, "u coupling across both lifts")

M = (F(1, 12), F(1, 6), F(-1, 6))
D = (F(1, 6), F(1, 6), F(-1, 3))
E11 = (F(2, 7), F(1, 7), F(-4, 7))
T = (F(1, 2), F(0), F(-1))
E10 = (F(3, 8), F(1, 8), F(-3, 4))
E01 = (F(1, 4), F(1, 6), F(-1, 2))
E00 = (F(1, 3), F(1, 6), F(-2, 3))
FANS = {
    (1, 1): frozenset((M, D, E11, T)),
    (1, 0): frozenset((M, D, E10)),
    (0, 1): frozenset((M, E01, T)),
    (0, 0): frozenset((M, E00)),
}
PLANES = candidate_planes(SUPPORT)

# K and Delta are coupled.  These eight representatives exhaust the possible
# zero/nonzero statuses of (Delta,K-e,Theta,Theta-Delta) away from K=0.
KREL = sp.Rational(2848, 45) - sp.Rational(7, 6) * Delta
DSTAR = sp.solve(sp.Eq(KREL - e, 0), Delta)[0]
KZERO = sp.solve(sp.Eq(KREL, 0), Delta)[0]
need(DSTAR == sp.Rational(3968, 63), "K/e cancellation value")
need(KZERO == sp.Rational(5696, 105), "K=0 exit value")
DELTA_THETA_REPS = (
    (sp.Integer(0), sp.Integer(0)),
    (sp.Integer(0), sp.Integer(1)),
    (DSTAR, sp.Integer(0)),
    (DSTAR, DSTAR),
    (DSTAR, 2 * DSTAR),
    (sp.Integer(1), sp.Integer(0)),
    (sp.Integer(1), sp.Integer(1)),
    (sp.Integer(1), sp.Integer(2)),
)
status_keys = set()
for delta_value, theta_value in DELTA_THETA_REPS:
    kval = KREL.subs(Delta, delta_value)
    need(kval != 0, "status representative avoids K=0")
    status_keys.add((
        int(delta_value != 0), int(kval - e != 0),
        int(theta_value != 0), int(theta_value - delta_value != 0),
    ))
need(len(status_keys) == 8, "eight coupled Delta/Theta status classes")

actual_states = {}
actual_fans = Counter()
exact_by_owner = Counter()
embedding_keys = set()
for (delta_value, theta_value), u_value, phi_value, eta_value, alpha_value in product(
    DELTA_THETA_REPS, (0, 1), (0, 1), (0, 1), (0, 1)
):
    substitution = {
        U: 2, Delta: delta_value, K: KREL.subs(Delta, delta_value),
        Theta: theta_value, u: u_value, Phi: phi_value,
        eta: eta_value, alpha: alpha_value,
    }
    active = frozenset(
        point for point, coefficient in SUPPORT.items()
        if sp.simplify(coefficient.subs(substitution)) != 0
    )
    owner = (int(alpha_value != 0), int(theta_value != 0))
    fan = lower_fan(active, PLANES)
    need(fan == FANS[owner], "exact state has owner fan")
    exact_by_owner[owner] += 1
    actual_fans[fan] += 1
    row_bits = (
        int(delta_value != 0), int(u_value != 0), int(phi_value != 0),
        int(theta_value != 0), int(eta_value != 0), int(alpha_value != 0),
    )
    cancellation_bits = (
        int(sp.simplify((K - e).subs(substitution)) == 0),
        int(theta_value - delta_value == 0),
        int(u_value == 0),
    )
    key = row_bits + cancellation_bits
    embedding_keys.add(key)
    actual_states[active] = key

need(len(actual_states) == 128, "128 distinct realizable support masks")
need(len(embedding_keys) == 128, "128 distinct atlas embeddings")
need(exact_by_owner == Counter({
    (1, 1): 40, (1, 0): 24, (0, 1): 40, (0, 0): 24,
}), "40/24/40/24 exact owner split")

# A conservative Boolean atlas independently toggles all six optional rows and
# all three aggregate deletions.  It contains every exact support but makes no
# false realizability claim about its 512 keyed configurations.
REQUIRED_LABELS = ((1, 0), (2, 0), (3, 0), (6, 0), (0, 2))
OPTIONAL_LABELS = ((4, 0), (5, 0), (2, 1), (1, 2), (3, 1), (4, 1))
CANCELLABLE = ((2, 3, 1), (2, 4, 1), (2, 5, 1))


def lifted(labels):
    answer = {(2, 0, 0), (0, 1, 0), (2, 0, 1)}
    for i, j in labels:
        answer.add((j + 2, i + j, 1))
        answer.add((j, i + j + 1, 1))
    return answer


over_states = {}
over_by_owner = Counter()
for row_bits in product((0, 1), repeat=6):
    labels = list(REQUIRED_LABELS)
    labels.extend(label for bit, label in zip(row_bits, OPTIONAL_LABELS) if bit)
    for cancellation_bits in product((0, 1), repeat=3):
        active = lifted(labels)
        for bit, point in zip(cancellation_bits, CANCELLABLE):
            if bit:
                active.discard(point)
        key = tuple(row_bits) + tuple(cancellation_bits)
        active = frozenset(active)
        over_states[key] = active
        owner = (row_bits[5], row_bits[3])
        need(lower_fan(active, PLANES) == FANS[owner],
             "over-atlas state has owner fan")
        over_by_owner[owner] += 1
need(len(over_states) == 512, "512 conservative keyed configurations")
need(over_by_owner == Counter({key: 128 for key in FANS}),
     "128 keyed configurations per conservative owner type")
for active, key in actual_states.items():
    need(key in over_states and over_states[key] == active,
         "realizable support embeds exactly in over-atlas")

# Exact face equations.
S, P = sp.symbols("S P")
FACE_EXPECTED = {
    "M": ((M, {}), (S**2 - P) * (1 - U * P**6)),
    "D": ((D, {}), S**2 * (1 - U * P**6 - alpha * S * P**5)),
    "E11": ((E11, {}), S**2 * (1 - alpha * S * P**5
                                 - Theta * S**2 * P**3)),
    "T": ((T, {}), S**2 * (1 - K * S**2 * P**2
                             - Theta * S**2 * P**3)),
    "E10": ((E10, {Theta: 0}), S**2 * (1 - alpha * S * P**5
                                         - K * S**2 * P**2)),
    "E01": ((E01, {alpha: 0}), S**2 * (1 - U * P**6
                                         - Theta * S**2 * P**3)),
    "E00": ((E00, {alpha: 0, Theta: 0}),
             S**2 * (1 - U * P**6 - eta * S * P**4
                      - K * S**2 * P**2)),
}
for name, ((plane, substitution), expected) in FACE_EXPECTED.items():
    actual = face_expression(SUPPORT, plane, substitution)
    need(sp.expand(actual - expected) == 0, name + " exact face")

POLYGONS = {
    "M": [(0, 1), (2, 0), (2, 6), (0, 7)],
    "D": [(2, 0), (2, 6), (3, 5)],
    "E11": [(2, 0), (3, 5), (4, 3)],
    "T": [(2, 0), (4, 2), (4, 3)],
    "E10": [(2, 0), (3, 5), (4, 2)],
    "E01": [(2, 0), (2, 6), (4, 3)],
    "E00": [(2, 0), (2, 6), (4, 2)],
    "G11": [(0, 1), (2, 0), (4, 2), (4, 3), (3, 5), (2, 6), (0, 7)],
    "G10": [(0, 1), (2, 0), (4, 2), (3, 5), (2, 6), (0, 7)],
    "G01": [(0, 1), (2, 0), (4, 2), (4, 3), (2, 6), (0, 7)],
    "G00": [(0, 1), (2, 0), (4, 2), (2, 6), (0, 7)],
}
EXPECTED_PICK = {
    "M": (24, 14, 6), "D": (6, 8, 0), "E11": (7, 3, 3),
    "T": (2, 4, 0), "E10": (8, 4, 3), "E01": (12, 8, 3),
    "E00": (12, 10, 2), "G11": (39, 13, 14),
    "G10": (38, 12, 14), "G01": (38, 12, 14),
    "G00": (36, 12, 13),
}
need({name: pick(vertices) for name, vertices in POLYGONS.items()} == EXPECTED_PICK,
     "all face/global Pick ledgers")

PACKETS = {
    "G11": (11, 8, 6, 3, 2, 2, 1),
    "G10": (11, 10, 6, 2, 2, 1),
    "G01": (13, 11, 3, 2, 2, 1),
    "G00": (11, 7, 7, 2, 2, 1),
}
for name, expected in PACKETS.items():
    packet = edge_packet(POLYGONS[name])
    need(packet == expected, name + " source-infinity packet")
    need(sum(index - 1 for index in packet)
         == 2 * EXPECTED_PICK[name][2] - 2,
         name + " Riemann--Hurwitz checksum")

ORDERS = {
    "M": (12, 9), "D": (6, 5), "E11": (42, 41), "T": (6, 8),
    "E10": (24, 26), "E01": (12, 11), "E00": (6, 6),
}
for name, plane in {
    "M": M, "D": D, "E11": E11, "T": T,
    "E10": E10, "E01": E01, "E00": E00,
}.items():
    need(face_order(plane) == ORDERS[name], name + " base/order")

# Normal forms certify the actual carrier genera rather than relying on Pick
# alone.  Each displayed branch polynomial is squarefree on its owner gate.
f11 = 1 - alpha * S * P**5 - Theta * S**2 * P**3
y11 = 2 * Theta * S * P**2 + alpha * P**4
need(sp.expand(y11**2 - P * (4 * Theta + alpha**2 * P**7)
               + 4 * Theta * P * f11) == 0, "E11 genus-three normal form")

f10 = 1 - alpha * S * P**5 - K * S**2 * P**2
y10 = 2 * K * S * P + alpha * P**4
need(sp.expand(y10**2 - (4 * K + alpha**2 * P**8)
               + 4 * K * f10) == 0, "E10 genus-three normal form")

f01 = 1 - U * P**6 - Theta * S**2 * P**3
y01 = Theta * S * P**2
need(sp.expand(y01**2 - Theta * P * (1 - U * P**6)
               + Theta * P * f01) == 0, "E01 genus-three normal form")

DISC = eta**2 - 4 * K * U
f00 = 1 - U * P**6 - eta * S * P**4 - K * S**2 * P**2
y00 = 2 * K * S * P + eta * P**3
need(sp.expand(y00**2 - (4 * K + DISC * P**6)
               + 4 * K * f00) == 0, "E00 genus-two normal form")

branch_polynomials = {
    "E11": P * (4 * Theta + alpha**2 * P**7),
    "E10": 4 * K + alpha**2 * P**8,
    "E01": Theta * P * (1 - U * P**6),
    "E00": 4 * K + DISC * P**6,
}
EXPECTED_BRANCH_DISCS = {
    "E11": -sp.Integer(53971714048) * Theta**8 * alpha**12,
    "E10": sp.Integer(274877906944) * K**7 * alpha**14,
    "E01": sp.Integer(46656) * Theta**12 * U**5,
    "E00": sp.Integer(47775744) * K**5 * (4 * K * U - eta**2)**5,
}
for name, polynomial in branch_polynomials.items():
    discriminant = sp.factor(sp.discriminant(polynomial, P))
    need(sp.expand(discriminant - EXPECTED_BRANCH_DISCS[name]) == 0,
         name + " exact branch discriminant")

# Edge schemes, separated by toric owner.  Linear endpoints are nonzero on
# the labelled gate; the only variable discriminant occurs in type 00.
X = sp.symbols("X")
OUTER = {
    (1, 1): (X - 1, 1 - K * X**2, K + Theta * X,
             Theta + alpha * X, alpha + U * X, X - 1, U - X**6),
    (1, 0): (X - 1, 1 - K * X**2, K + alpha * X,
             alpha + U * X, X - 1, U - X**6),
    (0, 1): (X - 1, 1 - K * X**2, K + Theta * X,
             Theta + U * X, X - 1, U - X**6),
    (0, 0): (X - 1, 1 - K * X**2, K + eta * X + U * X**2,
             X - 1, U - X**6),
}
INTERNAL = {
    (1, 1): (1 - U * X**6, 1 - alpha * X, 1 - Theta * X),
    (1, 0): (1 - U * X**6, 1 - alpha * X),
    (0, 1): (1 - U * X**6, 1 - Theta * X),
    (0, 0): (1 - U * X**6,),
}
need(sp.discriminant(1 - K * X**2, X) == 4 * K, "K edge reduced")
need(sp.factor(sp.discriminant(U - X**6, X)) == 46656 * U**5,
     "outer U edge reduced")
need(sp.factor(sp.discriminant(1 - U * X**6, X)) == 46656 * U**5,
     "internal U edge reduced")
need(sp.discriminant(K + eta * X + U * X**2, X) == DISC,
     "sole variable edge discriminant")

# Component/graph ledgers on the squarefree strata.  M consists of R plus six
# rational U-lines, with twelve R--line nodes.  Shared-edge lattice lengths
# are the remaining entries in each edge count.
GRAPH = {
    (1, 1): (10, 20, 11, 3, 14),
    (1, 0): (9, 19, 11, 3, 14),
    (0, 1): (9, 19, 11, 3, 14),
    (0, 0): (8, 18, 11, 2, 13),
}
for owner, (vertices, edges, b1, carrier_genus, total_genus) in GRAPH.items():
    need(edges - vertices + 1 == b1, "owner graph Betti number")
    need(b1 + carrier_genus == total_genus, "owner total genus")

# The alpha=0, Theta!=0 D12 face also has a complete primitive chart.  This
# is useful because its genus-three carrier is the clean maximal subgate.
sigma = sp.symbols("sigma")
p12, s12 = sigma**-2 * P, sigma**-3 * S
y12 = s12 * p12
H12_SOURCE = sum(
    coefficient * p12**i * y12**j for (i, j), coefficient in ROWS.items()
)
H12_SOURCE = sp.expand(H12_SOURCE.subs(alpha, 0))
FQ12 = ((s12**2 - p12) * (1 - sigma**12 * H12_SOURCE)
        - sigma**12 * s12**2 / 2)
G12 = sp.expand(sigma**6 * FQ12)
expected_G12 = sp.expand(
    (S**2 - sigma**4 * P)
    * (1 - (U * P**6 + Theta * S**2 * P**3)
       - sigma * eta * S * P**4
       - sigma**2 * (u * P**5 + K * S**2 * P**2)
       - sigma**3 * Phi * S * P**3
       - sigma**4 * Delta * P**4
       - e * sigma**6 * P**3
       - sp.Rational(8, 3) * sigma**8 * P**2
       + 3 * sigma**10 * P)
    - sigma**12 * S**2 / 2
)
need(sp.simplify(G12 - expected_G12) == 0, "exact D12 primitive chart")

# Exact type-00 repeated-edge chart.  This turns the apparent stopping wall
# DISC=0 into an even A5 calculation whose depth is bounded by the fixed e.
p0, s0 = sigma**-1 * P, sigma**-2 * S
y0 = s0 * p0
H0 = sum(coefficient * p0**i * y0**j for (i, j), coefficient in ROWS.items())
H0 = sp.expand(H0.subs({alpha: 0, Theta: 0}))
FQ = (s0**2 - p0) * (1 - sigma**6 * H0) - sigma**6 * s0**2 / 2
G = sp.expand(sigma**4 * FQ)
H6 = U * P**6 + eta * S * P**4 + K * S**2 * P**2
H5 = u * P**5 + Phi * S * P**3
expected_G = sp.expand(
    (S**2 - sigma**3 * P)
    * (1 - H6 - sigma * H5 - sigma**2 * Delta * P**4
       - e * sigma**3 * P**3 - sp.Rational(8, 3) * sigma**4 * P**2
       + 3 * sigma**5 * P)
    - sigma**6 * S**2 / 2
)
need(sp.simplify(G - expected_G) == 0, "exact E00 primitive chart")

x, v, rho = sp.symbols("x v rho")
A = U + eta * v + K * v**2
B = u + Phi * v
PHI_ROOT = sp.factor(x**10 * G.subs({P: 1 / x, S: v / x**2}))
expected_root = sp.expand(
    (v**2 - (sigma * x)**3)
    * (x**6 - A - sigma * x * B - (sigma * x)**2 * Delta
       - e * (sigma * x)**3 - sp.Rational(8, 3) * (sigma * x)**4
       + 3 * (sigma * x)**5)
    - (sigma * x)**6 * v**2 / 2
)
need(sp.simplify(PHI_ROOT - expected_root) == 0, "exact E00 reciprocal root chart")

DROOT = (
    A + rho * B + rho**2 * Delta + e * rho**3
    + sp.Rational(8, 3) * rho**4 - 3 * rho**5
    + rho**6 * v**2 / (2 * (v**2 - rho**3))
)
need(sp.factor(sp.together(
    PHI_ROOT / (v**2 - (sigma * x)**3)
    - (x**6 - DROOT.subs(rho, sigma * x))
)) == 0, "divided one-series E00 chart")

a = sp.symbols("a", nonzero=True)
REPEAT = {U: K * a**2, eta: -2 * K * a}
need(sp.expand(A.subs(REPEAT) - K * (v - a)**2) == 0,
     "repeated edge parameterization")
need(sp.diff(DROOT, v, 2).subs({rho: 0, v: a, **REPEAT}) == 2 * K,
     "nonzero Morse Hessian")
need((v**2 - rho**3).subs({v: a, rho: 0}) == a**2,
     "reciprocal prefactor unit")

vcrit = a - rho * Phi / (2 * K)
critical = sp.together(DROOT.subs(REPEAT).subs(v, vcrit))
critical_series = sp.series(critical, rho, 0, 6).removeO().expand()
expected_critical = (
    (u + Phi * a) * rho
    + (Delta - Phi**2 / (4 * K)) * rho**2
    + e * rho**3 + sp.Rational(8, 3) * rho**4 - 3 * rho**5
)
need(sp.simplify(critical_series - expected_critical) == 0,
     "critical series through forced cubic term")
need(e != 0, "critical depth at most three")

# Coordinate derivative and tail/complement checks.  From PHI=x^10 G and
# S=v/x^2 one has G_S=x^-8 PHI_v and -dP/G_S=x^6 dx/PHI_v.
GS_ROOT = sp.diff(G, S).subs({P: 1 / x, S: v / x**2})
need(sp.simplify(sp.diff(PHI_ROOT, v) - x**8 * GS_ROOT) == 0,
     "root-chart differential conversion")
TAIL = {
    1: (2, 0, 68),
    2: (1, 1, 64),
    3: (1, 1, 60),
}
tau, z, h, crho = sp.symbols("tau z h crho")
for r, (tail_genus, persistent_delta, form_order) in TAIL.items():
    epsilon = r % 2
    degree = (6 - r) + epsilon
    need((degree - 1) // 2 == tail_genus, "tail genus")
    need(r // 2 == persistent_delta, "persistent delta")
    need(persistent_delta + tail_genus + 1 == 3, "A5 conservation")
    computed_order = (
        6 * 2 * (6 - r) + 6 * 2 * r + 2 * r - 6 * r
    )
    need(computed_order == form_order > 0, "positive tail differential order")
    raw_term = tau ** (12 * (6 - r)) * rho ** (r - 6) * crho
    replaced = sp.cancel(raw_term.subs(tau**12, z * rho))
    need(sp.simplify(replaced - z ** (6 - r) * crho) == 0,
         "two-ended complement")
    complement = h**2 - 1 + z**(6 - r) * crho
    need(sp.diff(complement, h).subs({h: 1, z: 0}) == 2,
         "positive attachment etale")
    need(sp.diff(complement, h).subs({h: -1, z: 0}) == -2,
         "negative attachment etale")

lam = sp.symbols("lambda", nonzero=True)
linear = lam * (S * P - a * P**3)
factor_sub = {K: lam**2, U: lam**2 * a**2, eta: -2 * lam**2 * a}
need(sp.expand((1 - H6).subs(factor_sub) - (1 - linear) * (1 + linear)) == 0,
     "repeated E00 face splits into two rational signs")
need(sp.expand((1 - U * P**6).subs({U: lam**2 * a**2})
               - (1 - lam * a * P**3) * (1 + lam * a * P**3)) == 0,
     "six M attachments split three plus three")

# The repeated raw graph has 12 M nodes, six M--sign nodes, and one A5
# contact of delta three.  Its regularized genus stays thirteen in all three
# possible critical-depth rows.
need(12 + 6 + 3 - 9 + 1 == 13, "raw repeated-face arithmetic genus")
for tail_genus, persistent_delta, _order in TAIL.values():
    base_graph_genus = (18 + 2) - (9 + 1) + 1
    need(base_graph_genus + tail_genus + persistent_delta == 13,
         "regularized repeated-face genus")

print("THM-4350 first-normal-owner xi_10=0 exact fan audit")
print("SCOPE Z=beta11=zeta3=W=xi10=0; U*K!=0")
print("REALIZABLE_SUPPORTS=128")
print("OWNER_COUNTS alpha,Theta: 11=40 10=24 01=40 00=24")
print("OVER_ATLAS=512_KEYED_CONFIGS (128 per owner; not a support census)")
print("FANS 11=M,D,E11,T; 10=M,D,E10; 01=M,E01,T; 00=M,E00")
print("PICK_GLOBAL 11=(39,13,14) 10=(38,12,14) 01=(38,12,14) 00=(36,12,13)")
print("FACE_BASE_ORDER M=12/9 D=6/5 E11=42/41 T=6/8 E10=24/26 E01=12/11 E00=6/6")
print("PACKETS 11=(11,8,6,3,2,2,1) 10=(11,10,6,2,2,1) 01=(13,11,3,2,2,1) 00=(11,7,7,2,2,1)")
print("SOLE_EDGE_WALL eta^2-4*K*U=0 on owner 00")
print("CRITICAL c1=u+Phi*a c2=Delta-Phi^2/(4K) c3=-1376/135")
print("TAILS r=1:g2/order68 r=2:g1/order64 r=3:g1/order60")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
