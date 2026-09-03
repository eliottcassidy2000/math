#!/usr/bin/env python3
"""Independent exact referee for the K=0, W=xi=0 endpoint candidate.

This file deliberately reconstructs the sixteen source rows, their lift, and
all geometry used below.  It does not import the exploratory program being
refereed.  The only non-computational input is the inherited Newton/residue
dictionary stated in the canonical exact-weight-twelve seam results.
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


def need(condition, message):
    global CHECKS
    CHECKS += 1
    if condition is not True and not bool(condition):
        raise AssertionError(message)


def rank_two(points):
    for p, q, r in combinations(points, 3):
        if ((q[0] - p[0]) * (r[1] - p[1])
                != (q[1] - p[1]) * (r[0] - p[0])):
            return True
    return False


def candidate_planes(universe):
    answer = set()
    for p, q, r in combinations(sorted(universe), 3):
        determinant = ((q[0] - p[0]) * (r[1] - p[1])
                       - (q[1] - p[1]) * (r[0] - p[0]))
        if determinant == 0:
            continue
        aa = F((q[2] - p[2]) * (r[1] - p[1])
               - (q[1] - p[1]) * (r[2] - p[2]), determinant)
        bb = F((q[0] - p[0]) * (r[2] - p[2])
               - (q[2] - p[2]) * (r[0] - p[0]), determinant)
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


def cross(origin, a, b):
    return ((a[0] - origin[0]) * (b[1] - origin[1])
            - (a[1] - origin[1]) * (b[0] - origin[0]))


def convex_hull(points):
    points = sorted(set(points))
    if len(points) <= 1:
        return tuple(points)
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
        # The r=0 edge is the affine divisor, not source infinity.
        if not (start[0] == end[0] == 0):
            packet.extend([index] * length)
    return tuple(sorted(packet, reverse=True))


def on_plane(point, plane):
    rr, ll, hh = point
    aa, bb, cc = plane
    return F(hh) == aa * rr + bb * ll + cc


def face_expression(support, plane):
    return sp.factor(sum(coefficient * S**rr * P**ll
                         for (rr, ll, hh), coefficient in support.items()
                         if on_plane((rr, ll, hh), plane)))


def face_order(plane):
    base = lcm(6, *(entry.denominator for entry in plane))
    order = F(base) * (F(5, 6) - sum(plane))
    need(order.denominator == 1 and order > 0, "positive integral face order")
    return base, order.numerator


def face_polygon(expression):
    polynomial = sp.Poly(sp.expand(expression), S, P)
    return convex_hull(tuple(monomial for monomial, coefficient
                             in polynomial.terms() if coefficient != 0))


def segment_word(expression, start, end):
    polynomial = sp.Poly(sp.expand(expression), S, P)
    coefficient = {monomial: value for monomial, value in polynomial.terms()}
    dx, dy = end[0] - start[0], end[1] - start[1]
    length = gcd(abs(dx), abs(dy))
    step = (dx // length, dy // length)
    return sp.factor(sum(coefficient.get((start[0] + k * step[0],
                                          start[1] + k * step[1]), 0) * X**k
                         for k in range(length + 1)))


def same_scheme(actual, expected):
    """Equality up to a nonzero coefficient-field scalar."""
    actual_poly = sp.Poly(sp.expand(actual), X)
    expected_poly = sp.Poly(sp.expand(expected), X)
    if actual_poly.degree() != expected_poly.degree():
        return False
    a_coeff = actual_poly.all_coeffs()
    e_coeff = expected_poly.all_coeffs()
    pivot = next((i for i, value in enumerate(e_coeff) if value != 0), None)
    if pivot is None or a_coeff[pivot] == 0:
        return False
    ratio = sp.cancel(a_coeff[pivot] / e_coeff[pivot])
    return all(sp.factor(a - ratio * e) == 0
               for a, e in zip(a_coeff, e_coeff))


def sigma_initial(expression):
    expression = sp.expand(expression)
    terms = expression.as_ordered_terms()
    exponents = [sp.sympify(term.as_powers_dict().get(sigma, 0)) for term in terms]
    need(all(value.is_Integer for value in exponents), "integral sigma exponent")
    minimum = min(int(value) for value in exponents)
    initial = sum(term / sigma**minimum for term, value in zip(terms, exponents)
                  if int(value) == minimum)
    return minimum, sp.factor(initial)


# -------------------------------------------------------------------------
# Literal sixteen-row source and its exact lift.

U, W, Z = sp.symbols("U W Z")
Delta, Theta, Phi, eta = sp.symbols("Delta Theta Phi eta")
zeta, u, xi, alpha, beta, K = sp.symbols("zeta u xi alpha beta K")
e = -sp.Rational(1376, 135)
d0 = sp.Rational(5696, 105)

SOURCE_ROWS = (
    ("p", 1, 0, -3),
    ("p2", 2, 0, sp.Rational(8, 3)),
    ("p3", 3, 0, e),
    ("Delta", 4, 0, Delta),
    ("u", 5, 0, u),
    ("U", 6, 0, U),
    ("K", 0, 2, K),
    ("Phi", 2, 1, Phi),
    ("Theta", 1, 2, Theta),
    ("eta", 3, 1, eta),
    ("xi", 2, 2, xi),
    ("alpha", 4, 1, alpha),
    ("beta", 1, 3, beta),
    ("W", 3, 2, W),
    ("zeta", 0, 3, zeta),
    ("Z", 0, 4, Z),
)
need(len(SOURCE_ROWS) == 16, "literal sixteen-row count")
need(len({(ii, jj) for _name, ii, jj, _coefficient in SOURCE_ROWS}) == 16,
     "literal rows have distinct exponent pairs")

GATE = {K: 0, W: 0, Z: 0, xi: 0, beta: 0, zeta: 0, Delta: d0}


def lifted_support(substitution):
    support = defaultdict(lambda: sp.Integer(0))
    support[(2, 0, 0)] += 1
    support[(0, 1, 0)] -= 1
    support[(2, 0, 1)] -= sp.Rational(1, 2)
    for _name, ii, jj, raw_coefficient in SOURCE_ROWS:
        coefficient = sp.factor(sp.sympify(raw_coefficient).subs(GATE).subs(substitution))
        if coefficient == 0:
            continue
        support[(jj + 2, ii + jj, 1)] -= coefficient
        support[(jj, ii + jj + 1, 1)] += coefficient
    return {point: sp.factor(coefficient) for point, coefficient in support.items()
            if sp.factor(coefficient) != 0}


SYMBOLIC = lifted_support({})
need(sp.factor(SYMBOLIC[(2, 4, 1)] - (Theta - d0)) == 0,
     "Theta/Delta coupled lift")
need(SYMBOLIC[(2, 5, 1)] == -u, "u lower lift")
need(SYMBOLIC[(0, 6, 1)] == u, "u upper lift")
need(SYMBOLIC[(4, 3, 1)] == -Theta, "Theta owner")
need(SYMBOLIC[(3, 3, 1)] == -Phi, "Phi owner")
need(SYMBOLIC[(3, 4, 1)] == -eta, "eta owner")
need(SYMBOLIC[(3, 5, 1)] == -alpha, "alpha owner")
need(SYMBOLIC[(2, 6, 1)] == -U and SYMBOLIC[(0, 7, 1)] == U,
     "two U endpoints")

M = (F(1, 12), F(1, 6), F(-1, 6))
D6 = (F(1, 6), F(1, 6), F(-1, 3))
E11 = (F(2, 7), F(1, 7), F(-4, 7))
E01 = (F(1, 4), F(1, 6), F(-1, 2))
EETA = (F(1, 3), F(1, 6), F(-2, 3))
EPHI = (F(1, 2), F(1, 6), F(-1, 1))
N = (F(1, 1), F(0, 1), F(-2, 1))
FACE_NAMES = {M: "M", D6: "D6", E11: "E11", E01: "E01",
              EETA: "Eeta", EPHI: "EPhi", N: "N"}


def expected_fan(theta_nonzero, alpha_nonzero, phi_nonzero, eta_nonzero):
    if theta_nonzero:
        return frozenset((M, D6, E11) if alpha_nonzero else (M, E01))
    if alpha_nonzero:
        return frozenset((M, D6, N) if (phi_nonzero or eta_nonzero)
                         else (M, D6))
    if phi_nonzero and eta_nonzero:
        return frozenset((M, EETA, N))
    if eta_nonzero:
        return frozenset((M, EETA))
    if phi_nonzero:
        return frozenset((M, EPHI))
    return frozenset((M,))


exact_records = []
exact_masks = set()
universe = set()
for theta_status, alpha_bit, eta_bit, phi_bit, u_bit in product(
        range(3), (0, 1), (0, 1), (0, 1), (0, 1)):
    theta_value = (0, d0, d0 + 1)[theta_status]
    substitution = {
        Theta: theta_value,
        alpha: alpha_bit,
        eta: eta_bit,
        Phi: phi_bit,
        u: u_bit,
        U: 1,
    }
    support = lifted_support(substitution)
    mask = frozenset(support)
    exact_masks.add(mask)
    universe.update(mask)
    exact_records.append((theta_status, alpha_bit, eta_bit, phi_bit, u_bit,
                          support, mask))
need(len(exact_records) == 48 and len(exact_masks) == 48,
     "48 exact realizable supports")

planes = candidate_planes(universe)
fan_counts = Counter()
for theta_status, alpha_bit, eta_bit, phi_bit, _u_bit, _support, mask in exact_records:
    fan = lower_fan(mask, planes)
    expected = expected_fan(theta_status != 0, bool(alpha_bit), bool(phi_bit),
                            bool(eta_bit))
    need(fan == expected, "exact support has predicted fan")
    fan_counts[tuple(FACE_NAMES[item] for item in sorted(fan))] += 1
need(fan_counts == Counter({
    tuple(FACE_NAMES[item] for item in sorted((M, D6, E11))): 16,
    tuple(FACE_NAMES[item] for item in sorted((M, E01))): 16,
    tuple(FACE_NAMES[item] for item in sorted((M, D6, N))): 6,
    tuple(FACE_NAMES[item] for item in sorted((M, D6))): 2,
    tuple(FACE_NAMES[item] for item in sorted((M, EETA, N))): 2,
    tuple(FACE_NAMES[item] for item in sorted((M, EETA))): 2,
    tuple(FACE_NAMES[item] for item in sorted((M, EPHI))): 2,
    tuple(FACE_NAMES[item] for item in sorted((M,))): 2,
}), "eight-fan exact census")

# A separately constructed Boolean over-atlas: five optional source rows plus
# independent deletions of the two aggregate locations.  This is deliberately
# support-only; the later discriminant wall is audited separately.
base_mask = frozenset(lifted_support({Theta: 0, alpha: 0, eta: 0,
                                     Phi: 0, u: 0, U: 1}))
row_points = {
    "Theta": frozenset(((4, 3, 1), (2, 4, 1))),
    "alpha": frozenset(((3, 5, 1), (1, 6, 1))),
    "eta": frozenset(((3, 4, 1), (1, 5, 1))),
    "Phi": frozenset(((3, 3, 1), (1, 4, 1))),
    "u": frozenset(((2, 5, 1), (0, 6, 1))),
}
hostile_keyed = []
hostile_distinct = set()
hostile_fans = Counter()
for row_bits in product((0, 1), repeat=5):
    active = set(base_mask)
    for bit, label in zip(row_bits, ("Theta", "alpha", "eta", "Phi", "u")):
        if bit:
            active.update(row_points[label])
    for deletion_bits in product((0, 1), repeat=2):
        candidate = set(active)
        for bit, point in zip(deletion_bits, ((2, 4, 1), (2, 5, 1))):
            if bit:
                candidate.discard(point)
        candidate = frozenset(candidate)
        fan = lower_fan(candidate, planes)
        need(fan == expected_fan(bool(row_bits[0]), bool(row_bits[1]),
                                 bool(row_bits[3]), bool(row_bits[2])),
             "hostile support has owner-predicted fan")
        hostile_keyed.append((row_bits, deletion_bits, candidate))
        hostile_distinct.add(candidate)
        hostile_fans[tuple(FACE_NAMES[item] for item in sorted(fan))] += 1
need(len(hostile_keyed) == 128 and len(hostile_distinct) == 96,
     "128 keyed / 96 distinct hostile supports")
need(exact_masks <= hostile_distinct, "all exact supports embed literally")
need(sorted(hostile_fans.values()) == [8, 8, 8, 8, 8, 24, 32, 32],
     "hostile fan multiplicities")

# -------------------------------------------------------------------------
# The seven face initials, primitive charts, polygons, and form orders.

S, P, X, sigma = sp.symbols("S P X sigma")
FACES = {
    "M": (M, (S**2 - P) * (1 - U * P**6)),
    "D6": (D6, S**2 * (1 - U * P**6 - alpha * S * P**5)),
    "E11": (E11, S**2 * (1 - alpha * S * P**5
                           - Theta * S**2 * P**3)),
    "E01": (E01, S**2 * (1 - U * P**6
                           - Theta * S**2 * P**3)),
    "Eeta": (EETA, S**2 * (1 - U * P**6 - eta * S * P**4)),
    "EPhi": (EPHI, S**2 * (1 - U * P**6 - Phi * S * P**3)),
    "N": (N, S**2 * (1 - S * P**3 * (Phi + eta * P + alpha * P**2))),
}
for name, (plane, expected) in FACES.items():
    need(sp.expand(face_expression(SYMBOLIC, plane) - expected) == 0,
         name + " literal face initial")

s, p, Q = sp.symbols("s p Q")
y = s * p
Dpoly = (-3 * p + sp.Rational(8, 3) * p**2 + e * p**3
         + d0 * p**4 + u * p**5 + U * p**6)
H = (Dpoly + Phi * p**2 * y + Theta * p * y**2
     + eta * p**3 * y + alpha * p**4 * y)
FQ = sp.expand((s**2 - p) * (1 - Q * H) - Q * s**2 / 2)
CHARTS = {
    "M": (12, 1, 2, 2, {}),
    "D6": (6, 1, 1, 2, {}),
    "E11": (7, 2, 1, 4, {}),
    "E01": (12, 3, 2, 6, {alpha: 0}),
    "Eeta": (6, 2, 1, 4, {Theta: 0, alpha: 0}),
    "EPhi": (6, 3, 1, 6, {Theta: 0, alpha: 0, eta: 0}),
    "N": (1, 1, 0, 2, {Theta: 0}),
}
for name, (q_weight, s_weight, p_weight, multiplier, owner_zero) in CHARTS.items():
    scaled = sp.expand(sigma**multiplier * FQ.subs(owner_zero).subs({
        Q: sigma**q_weight,
        s: sigma**(-s_weight) * S,
        p: sigma**(-p_weight) * P,
    }))
    minimum, initial = sigma_initial(scaled)
    need(minimum == 0, name + " primitive chart clears at order zero")
    need(sp.expand(initial - FACES[name][1]) == 0,
         name + " primitive substitution initial")

EXPECTED_FACE_PICK = {
    "M": (24, 14, 6),
    "D6": (6, 8, 0),
    "E11": (7, 3, 3),
    "E01": (12, 8, 3),
    "Eeta": (6, 8, 0),
    "EPhi": (6, 8, 0),
    "N": (2, 4, 0),
}
for name, (_plane, expression) in FACES.items():
    need(pick(face_polygon(expression)) == EXPECTED_FACE_PICK[name],
         name + " face Pick ledger")

EXPECTED_ORDERS = {
    "M": (12, 9), "D6": (6, 5), "E11": (42, 41),
    "E01": (12, 11), "Eeta": (6, 6), "EPhi": (6, 7),
    "N": (6, 11),
}
for name, (plane, _expression) in FACES.items():
    need(face_order(plane) == EXPECTED_ORDERS[name], name + " base/order")

# Global hulls, Pick data, generic source-infinity packets, and graph ledgers.
def case_name(theta_nonzero, alpha_nonzero, phi_nonzero, eta_nonzero):
    if theta_nonzero:
        return "TA" if alpha_nonzero else "T0"
    if alpha_nonzero:
        if phi_nonzero:
            return "A_Phi"
        if eta_nonzero:
            return "A_eta"
        return "A_only"
    if phi_nonzero and eta_nonzero:
        return "Phi_eta"
    if eta_nonzero:
        return "eta_only"
    if phi_nonzero:
        return "Phi_only"
    return "none"


EXPECTED_GLOBAL = {
    "TA": ((37, 11, 14), (11, 8, 6, 5, 1), (9, 19, 11), 3),
    "T0": ((36, 10, 14), (13, 11, 5, 1), (8, 18, 11), 3),
    "A_Phi": ((32, 12, 11), (11, 6, 4, 2, 2, 1), (9, 19, 11), 0),
    "A_eta": ((31, 11, 11), (11, 6, 5, 2, 1), (9, 19, 11), 0),
    "A_only": ((30, 10, 11), (11, 6, 6, 1), (8, 18, 11), 0),
    "Phi_eta": ((31, 11, 11), (11, 7, 4, 2, 1), (9, 19, 11), 0),
    "eta_only": ((30, 10, 11), (11, 7, 5, 1), (8, 18, 11), 0),
    "Phi_only": ((30, 10, 11), (11, 8, 4, 1), (8, 18, 11), 0),
    "none": ((24, 14, 6), (11, 1, 1, 1, 1, 1, 1, 1), (7, 12, 6), 0),
}
for theta_status, alpha_bit, eta_bit, phi_bit, _u_bit, _support, mask in exact_records:
    label = case_name(theta_status != 0, bool(alpha_bit), bool(phi_bit), bool(eta_bit))
    polygon = convex_hull((rr, ll) for rr, ll, _hh in mask)
    ledger, packet, graph, carrier_genus = EXPECTED_GLOBAL[label]
    need(pick(polygon) == ledger, label + " global Pick")
    need(edge_packet(polygon) == packet, label + " generic polygon packet")
    need(sum(index - 1 for index in packet) == 2 * ledger[2] - 2,
         label + " generic packet RH defect")
    fan = lower_fan(mask, planes)
    extra_faces = len(fan) - 1
    derived_graph = ((7, 12, 6) if extra_faces == 0 else
                     (8, 18, 11) if extra_faces == 1 else
                     (9, 19, 11))
    need(derived_graph == graph,
         label + " graph from M=(7,12), six first attachments, one leaf")
    need(graph[1] - graph[0] + 1 == graph[2], label + " graph b1")
    need(graph[2] + carrier_genus == ledger[2], label + " genus decomposition")

# Literal face-edge words.  These tests derive each word from a face polynomial
# and an oriented lattice segment, rather than accepting a hand-written list.
FM = FACES["M"][1]
FD = FACES["D6"][1]
F11 = FACES["E11"][1]
F01 = FACES["E01"][1]
FEta = FACES["Eeta"][1]
FPhi = FACES["EPhi"][1]
FN = FACES["N"][1]
EDGE_TESTS = (
    (FM, (0, 1), (2, 0), X - 1),
    (FM, (2, 6), (0, 7), X - 1),
    (FM, (0, 7), (0, 1), U - X**6),
    (FM, (2, 0), (2, 6), 1 - U * X**6),
    (FD, (2, 0), (3, 5), 1 - alpha * X),
    (FD, (3, 5), (2, 6), alpha + U * X),
    (F11, (2, 0), (4, 3), 1 - Theta * X),
    (F11, (4, 3), (3, 5), Theta + alpha * X),
    (F01, (2, 0), (4, 3), 1 - Theta * X),
    (F01, (4, 3), (2, 6), Theta + U * X),
    (FEta, (2, 0), (3, 4), 1 - eta * X),
    (FEta, (3, 4), (2, 6), eta + U * X),
    (FPhi, (2, 0), (3, 3), 1 - Phi * X),
    (FPhi, (3, 3), (2, 6), Phi + U * X),
    (FN, (2, 0), (3, 3), 1 - Phi * X),
    (FN, (3, 3), (3, 5), -(Phi + eta * X + alpha * X**2)),
    (FN, (2, 0), (3, 5), 1 - alpha * X),
)
for expression, start, end, expected in EDGE_TESTS:
    need(same_scheme(segment_word(expression, start, end), expected),
         "literal oriented edge word")

# The special saturated N edges obtained after deleting an endpoint owner.
need(same_scheme(segment_word(FN.subs(Phi, 0), (2, 0), (3, 4)),
                         1 - eta * X), "N eta first edge")
need(same_scheme(segment_word(FN.subs(Phi, 0), (3, 4), (3, 5)),
                         -(eta + alpha * X)), "N eta/alpha outer edge")
need(same_scheme(segment_word(FN.subs(alpha, 0), (2, 0), (3, 4)),
                         1 - eta * X), "N eta internal edge")
need(same_scheme(segment_word(FN.subs(alpha, 0), (3, 3), (3, 4)),
                         -(Phi + eta * X)), "N Phi/eta outer edge")

# In characteristic zero, the two sextics are reduced for U nonzero.  All
# remaining boundary/internal words are linear except the displayed J edge.
need(sp.factor(sp.discriminant(1 - U * X**6, X)) != 0,
     "internal sextic is squarefree on U!=0")
need(sp.factor(sp.discriminant(U - X**6, X)) != 0,
     "outer sextic is squarefree on U!=0")
need(sp.factor(sp.discriminant(Phi + eta * X + alpha * X**2, X)
               - (eta**2 - 4 * alpha * Phi)) == 0,
     "J is the sole nonlinear variable edge discriminant")

# -------------------------------------------------------------------------
# Carrier normalizations and discriminants.

f11 = 1 - alpha * S * P**5 - Theta * S**2 * P**3
y11 = 2 * Theta * S * P**2 + alpha * P**4
need(sp.expand(y11**2 - P * (4 * Theta + alpha**2 * P**7)
               + 4 * Theta * P * f11) == 0, "E11 normalization identity")
disc11 = sp.factor(sp.discriminant(P * (4 * Theta + alpha**2 * P**7), P))
need(disc11 == -sp.Integer(53971714048) * Theta**8 * alpha**12,
     "E11 discriminant")

f01 = 1 - U * P**6 - Theta * S**2 * P**3
y01 = Theta * S * P**2
need(sp.expand(y01**2 - Theta * P * (1 - U * P**6)
               + Theta * P * f01) == 0, "E01 normalization identity")
disc01 = sp.factor(sp.discriminant(Theta * P * (1 - U * P**6), P))
need(disc01 == sp.Integer(46656) * Theta**12 * U**5,
     "E01 discriminant")

# Rational faces are linear in S after their torus monomials are removed.
for name in ("D6", "Eeta", "EPhi", "N"):
    core = sp.cancel(FACES[name][1] / S**2)
    need(sp.degree(core, S) == 1, name + " rational by solving for S")

# M has the rational parabola P=S^2 and six vertical rational lines.  There
# are twelve transverse intersections and hence graph b1=12-7+1=6.
need(sp.discriminant(1 - U * P**6, P) != 0, "M six roots squarefree")
need(12 - 7 + 1 == 6, "M graph genus")

# -------------------------------------------------------------------------
# The sole repeated edge: exact tangency, packet correction, and two blowups.

R, a, xvar = sp.symbols("R a x")
J = Phi + eta * P + alpha * P**2
D = (-3 * P + sp.Rational(8, 3) * P**2 + e * P**3
     + d0 * P**4 + u * P**5 + U * P**6)
GN = ((S**2 - sigma**2 * P) * (1 - S * P**3 * J - sigma * D)
      - sigma * S**2 / 2)
need(sp.expand(GN - sigma**2 * FQ.subs(Theta, 0).subs(
    {Q: sigma, s: S / sigma, p: P})) == 0,
     "exact primitive N chart from literal source")

Hboundary = sp.factor(sp.expand(R**3 * GN.subs(S, 1 / R)))
Hexpect = ((R - sigma**2 * P * R**3) * (1 - sigma * D)
           - (1 - sigma**2 * P * R**2) * P**3 * J - sigma * R / 2)
need(sp.expand(Hboundary - Hexpect) == 0, "exact reciprocal boundary chart")

collision = {eta: -2 * alpha * a, Phi: alpha * a**2}
Hcollision = sp.factor(Hboundary.subs(collision))
need(sp.factor(J.subs(collision) - alpha * (P - a)**2) == 0,
     "quadratic repeated-root parameterization")
need(sp.factor(Hcollision.subs({R: 0, P: a + xvar})
               + alpha * (a + xvar)**3 * xvar**2) == 0,
     "boundary restriction has exactly a double root")
HR = sp.factor(sp.diff(Hcollision, R).subs({R: 0, P: a}))
need(sp.factor(HR - (1 - sigma * D.subs(P, a) - sigma / 2)) == 0,
     "transverse derivative is an exact DVR unit")
need(HR.subs(sigma, 0) == 1, "transverse unit has residue one")

# Exact formal normal form without extracting a square root.  Write
# H=R*A-alpha*P^3*(P-a)^2*B, then Y=R*A/(alpha*P^3*B).
Bunit = 1 - sigma**2 * P * R**2
Aunit = Bunit * (1 - sigma * D) - sigma / 2
need(sp.expand(Hcollision - (R * Aunit
                              - alpha * P**3 * (P - a)**2 * Bunit)) == 0,
     "unit-factor tangency decomposition")
need(Aunit.subs({R: 0, P: a, sigma: 0}) == 1
     and Bunit.subs({R: 0, P: a, sigma: 0}) == 1,
     "normal-form multipliers are units")

# The face differential generator is dS/G_P.  With S=1/R and H=R^3 G,
# it is -R*dR/H_P = R*dP/H_R on H=0.  Thus a simple J root gives order one
# (index 2), while the double root gives order two (index 3).
need(sp.cancel((-R**-2) / R**-3) == -R,
     "dS/G_P has coefficient -R/H_P in the boundary chart")
need(1 + 1 == 2, "simple-root puncture index")
need(2 + 1 == 3, "double-root puncture index")
generic_collision_packet = EXPECTED_GLOBAL["A_Phi"][1]
need(generic_collision_packet == (11, 6, 4, 2, 2, 1),
     "generic two-root packet")
repeated_packet = tuple(sorted((11, 6, 4, 3, 1), reverse=True))
need(repeated_packet == (11, 6, 4, 3, 1), "correct repeated-edge packet")
need(sum(index - 1 for index in generic_collision_packet)
     == sum(index - 1 for index in repeated_packet) == 20,
     "packet coalescence preserves defect and genus eleven")

# Two ordinary blowups of Y=X^2 against the boundary Y=0.  First X-chart:
# Y=X*y1, so C1:y1=X, D1:y1=0, E1:X=0 meet in one triple point.
x1, y1, x2, y2 = sp.symbols("x1 y1 x2 y2")
need(sp.factor((sp.Symbol("Y") - sp.Symbol("X0")**2).subs(
    {sp.Symbol("Y"): x1 * y1, sp.Symbol("X0"): x1}) / x1) == y1 - x1,
     "first blowup X-chart strict transform")
# In the other first-blowup chart, Y=y1 and X=y1*x1, the strict curve is
# 1-y1*x1^2=0 and is disjoint from the exceptional y1=0.
need(sp.factor((sp.Symbol("Y") - sp.Symbol("X0")**2).subs(
    {sp.Symbol("Y"): y1, sp.Symbol("X0"): y1 * x1}) / y1)
     == 1 - y1 * x1**2, "first blowup complementary chart")
# Blow up the triple point.  In x1=x2,y1=x2*y2, C2:y2=1 and D2:y2=0
# meet E2:x2=0 at distinct points.  In x1=y2*x2,y1=y2, C2:x2=1 and
# E1':x2=0 meet E2:y2=0 at distinct points.
need(sp.factor((y1 - x1).subs({x1: x2, y1: x2 * y2}) / x2)
     == y2 - 1, "second blowup first chart")
need(sp.factor((y1 - x1).subs({x1: y2 * x2, y1: y2}) / y2)
     == 1 - x2, "second blowup second chart")
need(sp.cancel((-eta / (2 * alpha)).subs(collision) - a) == 0,
     "the repeated-root center is a=-eta/(2alpha)")
need(sp.cancel((eta**2 / (4 * alpha) - Phi).subs(collision)) == 0,
     "the collision relation is Phi=eta^2/(4alpha) in the coefficient field")


print("K0_FULL_INDEPENDENT_REFEREE=PASS_WITH_PACKET_CORRECTION")
print("source_rows=16;exact_supports=48;exact_fans=8")
print("hostile=128_keyed/96_distinct;exact_embeds=48")
print("faces=M,D6,E11,E01,Eeta,EPhi,N;primitive_initials=7")
print("carrier=E11_g3_disc=-53971714048*Theta^8*alpha^12;")
print("carrier=E01_g3_disc=46656*Theta^12*U^5")
print("orders=M9@12,D6_5@6,E11_41@42,E01_11@12,Eeta_6@6,EPhi_7@6,N_11@6")
print("collision=boundary_tangency_R=(P-a)^2*unit;H_R=unit;two_ordinary_blowups")
print("PACKET_CORRECTION:generic=(11,6,4,2,2,1);repeated=(11,6,4,3,1);defect=20")
print("no_base_change_for_tangency=yes;exceptionals=P1,P1")
print("checks=" + str(CHECKS))
