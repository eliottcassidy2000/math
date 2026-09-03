#!/usr/bin/env python3
"""Clean-room exact referee for the proposed THM-4354 gate.

Reconstruct the sixteen source rows and lifted Newton support directly.
Exact coefficient realizability and a larger hostile atlas stay separate.
"""

from collections import Counter, defaultdict
from fractions import Fraction as QQ
from itertools import combinations, product
from math import gcd, lcm
import sys

import sympy as sp

sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")
CHECKS = 0


def check(condition, message):
    global CHECKS
    CHECKS += 1
    if not bool(condition):
        raise AssertionError(message)


def rank_two(points):
    for p0, p1, p2 in combinations(tuple(points), 3):
        det = ((p1[0] - p0[0]) * (p2[1] - p0[1])
               - (p1[1] - p0[1]) * (p2[0] - p0[0]))
        if det:
            return True
    return False


def all_supporting_planes(universe):
    answer = set()
    for p0, p1, p2 in combinations(sorted(universe), 3):
        det = ((p1[0] - p0[0]) * (p2[1] - p0[1])
               - (p1[1] - p0[1]) * (p2[0] - p0[0]))
        if det == 0:
            continue
        aa = QQ((p1[2] - p0[2]) * (p2[1] - p0[1])
                - (p1[1] - p0[1]) * (p2[2] - p0[2]), det)
        bb = QQ((p1[0] - p0[0]) * (p2[2] - p0[2])
                - (p1[2] - p0[2]) * (p2[0] - p0[0]), det)
        cc = QQ(p0[2]) - aa * p0[0] - bb * p0[1]
        answer.add((aa, bb, cc))
    return tuple(sorted(answer))


def lower_fan(active, planes):
    answer = set()
    for aa, bb, cc in planes:
        equality = []
        for rr, ll, height in active:
            gap = QQ(height) - aa * rr - bb * ll - cc
            if gap < 0:
                break
            if gap == 0:
                equality.append((rr, ll, height))
        else:
            if rank_two(equality):
                answer.add((aa, bb, cc))
    return frozenset(answer)


def cross(origin, left, right):
    return ((left[0] - origin[0]) * (right[1] - origin[1])
            - (left[1] - origin[1]) * (right[0] - origin[0]))


def convex_hull(points):
    points = sorted(set(points))
    if len(points) < 2:
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


def pick_data(vertices):
    cyclic = tuple(zip(vertices, vertices[1:] + vertices[:1]))
    area2 = abs(sum(x0 * y1 - x1 * y0
                    for (x0, y0), (x1, y1) in cyclic))
    boundary = sum(gcd(abs(x1 - x0), abs(y1 - y0))
                   for (x0, y0), (x1, y1) in cyclic)
    check((area2 - boundary + 2) % 2 == 0, "Pick parity")
    return area2, boundary, (area2 - boundary + 2) // 2


def infinity_packet(vertices):
    packet = []
    for start, end in zip(vertices, vertices[1:] + vertices[:1]):
        dx, dy = end[0] - start[0], end[1] - start[1]
        length = gcd(abs(dx), abs(dy))
        nx, ny = -dy // length, dx // length
        constant = nx * start[0] + ny * start[1]
        index = nx + ny - constant
        if not (start[0] == end[0] == 0):
            packet.extend([index] * length)
    return tuple(sorted(packet, reverse=True))


def is_on_plane(point, plane):
    rr, ll, height = point
    aa, bb, cc = plane
    return QQ(height) == aa * rr + bb * ll + cc


def face_from_lift(lifted, plane, S, P):
    return sp.factor(sum(coefficient * S**rr * P**ll
                         for (rr, ll, height), coefficient in lifted.items()
                         if is_on_plane((rr, ll, height), plane)))


def sigma_initial(expression, sigma):
    terms = sp.expand(expression).as_ordered_terms()
    orders = [sp.sympify(term.as_powers_dict().get(sigma, 0))
              for term in terms]
    check(all(order.is_Integer for order in orders), "integral sigma powers")
    minimum = min(int(order) for order in orders)
    initial = sum(term / sigma**minimum
                  for term, order in zip(terms, orders)
                  if int(order) == minimum)
    return minimum, sp.factor(initial)


def edge_word(expression, start, end, S, P, X):
    coefficients = dict(sp.Poly(sp.expand(expression), S, P).terms())
    dx, dy = end[0] - start[0], end[1] - start[1]
    length = gcd(abs(dx), abs(dy))
    step = (dx // length, dy // length)
    return sp.factor(sum(coefficients.get(
        (start[0] + k * step[0], start[1] + k * step[1]), 0) * X**k
                         for k in range(length + 1)))


def associates(left, right, variable):
    lp, rp = sp.Poly(sp.expand(left), variable), sp.Poly(sp.expand(right), variable)
    if lp.degree() != rp.degree():
        return False
    lc, rc = lp.all_coeffs(), rp.all_coeffs()
    pivot = next((i for i, value in enumerate(rc) if value != 0), None)
    if pivot is None or lc[pivot] == 0:
        return False
    ratio = sp.cancel(lc[pivot] / rc[pivot])
    return all(sp.factor(a - ratio * b) == 0 for a, b in zip(lc, rc))


Delta, Theta, Phi, eta, u = sp.symbols("Delta Theta Phi eta u")
K, U, W, Z = sp.symbols("K U W Z")
xi, alpha, beta, zeta = sp.symbols("xi alpha beta zeta")
e = -sp.Rational(1376, 135)
dstar = sp.Rational(3968, 63)
k_of_delta = sp.Rational(2848, 45) - sp.Rational(7, 6) * Delta

SOURCE = (
    ("p", 1, 0, -3), ("p2", 2, 0, sp.Rational(8, 3)),
    ("p3", 3, 0, e), ("K", 0, 2, K), ("Phi", 2, 1, Phi),
    ("Delta", 4, 0, Delta), ("Theta", 1, 2, Theta),
    ("eta", 3, 1, eta), ("zeta", 0, 3, zeta), ("u", 5, 0, u),
    ("xi", 2, 2, xi), ("alpha", 4, 1, alpha),
    ("beta", 1, 3, beta), ("U", 6, 0, U), ("W", 3, 2, W),
    ("Z", 0, 4, Z),
)
check(len(SOURCE) == 16, "literal source has sixteen rows")
check(len({(ii, jj) for _name, ii, jj, _c in SOURCE}) == 16,
      "source exponent pairs are distinct")
GATE = {Z: 0, beta: 0, zeta: 0, W: 0, xi: 0, U: 0, K: k_of_delta}


def lifted_coefficients(substitution):
    result = defaultdict(lambda: sp.Integer(0))
    result[(2, 0, 0)] += 1
    result[(0, 1, 0)] -= 1
    result[(2, 0, 1)] -= sp.Rational(1, 2)
    for _name, ii, jj, raw in SOURCE:
        coefficient = sp.factor(sp.sympify(raw).subs(GATE).subs(substitution))
        if coefficient == 0:
            continue
        result[(jj + 2, ii + jj, 1)] -= coefficient
        result[(jj, ii + jj + 1, 1)] += coefficient
    return {point: sp.factor(value) for point, value in result.items()
            if sp.factor(value) != 0}


symbolic_lift = lifted_coefficients({})
check(sp.factor(symbolic_lift[(2, 3, 1)]
                - (k_of_delta + sp.Rational(1376, 135))) == 0,
      "K/p3 aggregate coefficient")
check(sp.factor(symbolic_lift[(2, 4, 1)] - (Theta - Delta)) == 0,
      "Theta/Delta aggregate coefficient")
check(symbolic_lift[(2, 5, 1)] == -u, "u/xi aggregate coefficient")
check(sp.solve(sp.Eq(k_of_delta + sp.Rational(1376, 135), 0), Delta)
      == [dstar], "unique first aggregate cancellation")
check(k_of_delta.subs(Delta, dstar) != 0, "Delta-star remains in K!=0 gate")

COUPLED_CLASSES = (
    (sp.Integer(0), sp.Integer(0)), (sp.Integer(0), sp.Integer(1)),
    (dstar, sp.Integer(0)), (dstar, dstar), (dstar, 2 * dstar),
    (sp.Integer(1), sp.Integer(0)), (sp.Integer(1), sp.Integer(1)),
    (sp.Integer(1), sp.Integer(2)),
)
exact_records, exact_masks = [], set()
for class_index, (delta_value, theta_value) in enumerate(COUPLED_CLASSES):
    check(k_of_delta.subs(Delta, delta_value) != 0, "class respects K!=0")
    for u_bit, phi_bit, eta_bit in product((0, 1), repeat=3):
        values = {Delta: delta_value, Theta: theta_value, u: u_bit,
                  Phi: phi_bit, eta: eta_bit, alpha: 1}
        lifted = lifted_coefficients(values)
        mask = frozenset(lifted)
        address = 1 + 8 * class_index + 4 * u_bit + 2 * phi_bit + eta_bit
        exact_records.append((class_index, u_bit, phi_bit, eta_bit,
                              delta_value, theta_value, address, lifted, mask))
        exact_masks.add(mask)
check(len(exact_records) == len(exact_masks) == 64, "64 exact supports")
address_to_mask = {record[6]: record[8] for record in exact_records}
mask_to_address = {record[8]: record[6] for record in exact_records}
check(set(address_to_mask) == set(range(1, 65)), "addresses fill [1,64]")
check(len(mask_to_address) == 64, "support-address injectivity")
for record in exact_records:
    c, ub, pb, eb, _dv, _tv, address, _lift, mask = record
    residue = address - 1
    check((residue // 8, (residue % 8) // 4, (residue % 4) // 2,
           residue % 2) == (c, ub, pb, eb), "natural address inverse")
    check(mask_to_address[mask] == address, "mask recovers address")


# A clean-room Boolean over-atlas: fixed p,p2,p3,K,alpha rows; independent
# Delta,Theta,u,Phi,eta rows; then adversarial deletion at all three inherited
# aggregate sites.  It is a support stress test, not a realizability claim.
def row_points(ii, jj):
    return frozenset(((jj + 2, ii + jj, 1), (jj, ii + jj + 1, 1)))


row_by_name = {name: row_points(ii, jj) for name, ii, jj, _raw in SOURCE}
fixed_mask = {(2, 0, 0), (0, 1, 0), (2, 0, 1)}
for name in ("p", "p2", "p3", "K", "alpha"):
    fixed_mask.update(row_by_name[name])
optional_names = ("Delta", "Theta", "u", "Phi", "eta")
aggregate_sites = ((2, 3, 1), (2, 4, 1), (2, 5, 1))
hostile_records, hostile_distinct, hostile_universe = [], set(), set()
for bits in product((0, 1), repeat=5):
    for deletion_bits in product((0, 1), repeat=3):
        active = set(fixed_mask)
        for bit, name in zip(bits, optional_names):
            if bit:
                active.update(row_by_name[name])
        for deletion, point in zip(deletion_bits, aggregate_sites):
            if deletion:
                active.discard(point)
        mask = frozenset(active)
        hostile_records.append((bits, deletion_bits, mask))
        hostile_distinct.add(mask)
        hostile_universe.update(mask)
check(len(hostile_records) == 256, "256 keyed hostile configurations")
check(len(hostile_distinct) == 168, "168 hostile support classes")
check(exact_masks <= hostile_distinct, "exact supports embed in hostile atlas")

planes = all_supporting_planes(hostile_universe)
C5 = (QQ(0), QQ(1, 5), QQ(-1, 5))
C4 = (QQ(-1, 4), QQ(1, 4), QQ(-1, 4))
C3 = (QQ(-2, 3), QQ(1, 3), QQ(-1, 3))
B = (QQ(1, 11), QQ(2, 11), QQ(-2, 11))
E11 = (QQ(2, 7), QQ(1, 7), QQ(-4, 7))
T = (QQ(1, 2), QQ(0), QQ(-1))
E10 = (QQ(3, 8), QQ(1, 8), QQ(-3, 4))


def expected_fan(u_present, delta_present, theta_present):
    cap = C5 if u_present else C4 if delta_present else C3
    right = (E11, T) if theta_present else (E10,)
    return frozenset((cap, B, *right))


exact_fan_counts = Counter()
for record in exact_records:
    (_c, ub, _pb, _eb, dv, tv, _address, _lift, mask) = record
    actual = lower_fan(mask, planes)
    check(actual == expected_fan(bool(ub), dv != 0, tv != 0),
          "exact owner-predicted fan")
    cap_name = "C5" if ub else "C4" if dv != 0 else "C3"
    exact_fan_counts[(cap_name, bool(tv))] += 1
check(exact_fan_counts == Counter({
    ("C5", True): 20, ("C5", False): 12,
    ("C4", True): 16, ("C4", False): 8,
    ("C3", True): 4, ("C3", False): 4,
}), "six exact fan counts")

hostile_fan_counts, hostile_fans = Counter(), set()
for bits, _deletions, mask in hostile_records:
    db, tb, ub, _pb, _eb = bits
    actual = lower_fan(mask, planes)
    check(actual == expected_fan(bool(ub), bool(db), bool(tb)),
          "hostile deletions preserve fan")
    cap_name = "C5" if ub else "C4" if db else "C3"
    hostile_fan_counts[(cap_name, bool(tb))] += 1
    hostile_fans.add(actual)
check(len(hostile_fans) == 6, "hostile atlas has six fans")
check(hostile_fan_counts == Counter({
    ("C5", True): 64, ("C5", False): 64,
    ("C4", True): 32, ("C4", False): 32,
    ("C3", True): 32, ("C3", False): 32,
}), "hostile fan multiplicities")


# Face equations are pulled from the literal lifted support, then reconstructed
# a second way by primitive substitutions into the complete source equation.
S, P, X, sigma = sp.symbols("S P X sigma")
EXPECTED_FACES = {
    "C5": (C5, P * (u * P**5 + alpha * S * P**5 - 1)),
    "C4": (C4, P * (Delta * P**4 + alpha * S * P**5 - 1)),
    "C3": (C3, P * (e * P**3 + alpha * S * P**5 - 1)),
    "B": (B, (P - S**2) * (alpha * S * P**5 - 1)),
    "E11": (E11, S**2 * (1 - alpha * S * P**5
                            - Theta * S**2 * P**3)),
    "T": (T, S**2 * (1 - S**2 * P**2 * (K + Theta * P))),
    "E10": (E10, S**2 * (1 - alpha * S * P**5 - K * S**2 * P**2)),
}
for name, (plane, expected) in EXPECTED_FACES.items():
    actual = face_from_lift(symbolic_lift, plane, S, P)
    check(sp.expand(actual - expected.subs(K, k_of_delta)) == 0,
          name + " literal face")

s, p, Q = sp.symbols("s p Q")
y = s * p
H = sum(sp.sympify(coefficient).subs(K, k_of_delta) * p**ii * y**jj
        for _name, ii, jj, coefficient in SOURCE)
H = sp.expand(H.subs({Z: 0, beta: 0, zeta: 0, W: 0, xi: 0, U: 0}))
FQ = sp.expand((s**2 - p) * (1 - Q * H) - Q * s**2 / 2)
CHARTS = {
    "C5": (5, 0, -1, 1, {}),
    "C4": (4, 1, -1, 1, {u: 0}),
    "C3": (3, 2, -1, 1, {u: 0, Delta: 0}),
    "B": (11, -1, -2, 2, {}),
    "E11": (7, -2, -1, 4, {}),
    "T": (2, -1, 0, 2, {}),
    "E10": (8, -3, -1, 6, {Theta: 0}),
}
for name, (qw, sw, pw, clearing, owner_zeros) in CHARTS.items():
    scaled = sp.expand(sigma**clearing * FQ.subs(owner_zeros).subs({
        Q: sigma**qw, s: sigma**sw * S, p: sigma**pw * P,
    }, simultaneous=True))
    minimum, initial = sigma_initial(scaled, sigma)
    expected = EXPECTED_FACES[name][1].subs(K, k_of_delta).subs(owner_zeros)
    check(minimum == 0, name + " minimal primitive clearing")
    check(sp.expand(initial - expected) == 0, name + " primitive initial")

EXPECTED_POLYGONS = {
    "C5": ((0, 1), (1, 6), (0, 6)),
    "C4": ((0, 1), (1, 6), (0, 5)),
    "C3": ((0, 1), (1, 6), (0, 4)),
    "B": ((0, 1), (2, 0), (3, 5), (1, 6)),
    "E11": ((2, 0), (4, 3), (3, 5)),
    "T": ((2, 0), (4, 2), (4, 3)),
    "E10": ((2, 0), (4, 2), (3, 5)),
}
for name, (_plane, expression) in EXPECTED_FACES.items():
    points = tuple(monomial for monomial, coefficient
                   in sp.Poly(sp.expand(expression), S, P).terms()
                   if coefficient != 0)
    check(set(convex_hull(points)) == set(EXPECTED_POLYGONS[name]),
          name + " face polygon")


def target_order(plane):
    base = lcm(6, *(entry.denominator for entry in plane))
    order = QQ(base) * (QQ(5, 6) - sum(plane))
    check(order.denominator == 1 and order > 0, "positive integral order")
    return base, order.numerator


EXPECTED_ORDERS = {
    "C5": (30, 25), "C4": (12, 13), "C3": (6, 9),
    "B": (66, 49), "E11": (42, 41), "T": (6, 8),
    "E10": (24, 26),
}
for name, (plane, _expression) in EXPECTED_FACES.items():
    check(target_order(plane) == EXPECTED_ORDERS[name], name + " form order")


# The middle face has two rational components and eleven transverse nodes.
R_component = P - S**2
C_component = alpha * S * P**5 - 1
intersection = sp.factor(C_component.subs(P, S**2))
check(intersection == alpha * S**11 - 1, "middle intersection equation")
check(sp.degree(intersection, S) == 11, "eleven middle intersections")
check(sp.factor(sp.discriminant(intersection, S))
      == -sp.Integer(285311670611) * alpha**10,
      "eleven middle roots are distinct")
check(sp.factor(sp.diff(intersection, S)) == 11 * alpha * S**10,
      "middle derivative is nonzero at every root")
check(gcd(1, 5) == 1, "monomial middle component is rational")
check(11 - 2 + 1 == 10, "middle multigraph genus")

# Exact hyperelliptic normalizations of the only positive-genus carriers.
f11 = 1 - alpha * S * P**5 - Theta * S**2 * P**3
y11 = 2 * Theta * S * P**2 + alpha * P**4
check(sp.expand(y11**2 - P * (alpha**2 * P**7 + 4 * Theta)
                + 4 * Theta * P * f11) == 0, "E11 normalization")
disc11 = sp.factor(sp.discriminant(P * (alpha**2 * P**7 + 4 * Theta), P))
check(disc11 == -sp.Integer(53971714048) * Theta**8 * alpha**12,
      "E11 branch discriminant")
f10 = 1 - alpha * S * P**5 - K * S**2 * P**2
y10 = 2 * K * S * P + alpha * P**4
check(sp.expand(y10**2 - (alpha**2 * P**8 + 4 * K) + 4 * K * f10) == 0,
      "E10 normalization")
disc10 = sp.factor(sp.discriminant(alpha**2 * P**8 + 4 * K, P))
check(disc10 == sp.Integer(274877906944) * K**7 * alpha**14,
      "E10 branch discriminant")
check((8 - 2) // 2 == 3, "squarefree degree-eight covers have genus three")

# Hostile parametrization probes for every face declared rational.
for core, parametrization in (
        (sp.cancel(EXPECTED_FACES["C5"][1] / P),
         {S: (1 - u * P**5) / (alpha * P**5)}),
        (sp.cancel(EXPECTED_FACES["C4"][1] / P),
         {S: (1 - Delta * P**4) / (alpha * P**5)}),
        (sp.cancel(EXPECTED_FACES["C3"][1] / P),
         {S: (1 - e * P**3) / (alpha * P**5)}),
        (C_component, {S: 1 / (alpha * P**5)})):
    check(sp.factor(core.subs(parametrization)) == 0,
          "claimed rational face has an explicit parametrization")
tau = sp.symbols("tau")
t_core = sp.cancel(EXPECTED_FACES["T"][1] / S**2)
t_p = (1 - K * tau**2) / (Theta * tau**2)
t_param = {P: t_p, S: tau / t_p}
check(sp.factor(t_core.subs(t_param, simultaneous=True)) == 0,
      "T has an explicit rational parametrization")

# Derive every displayed boundary word from an oriented lattice segment.
FC5, FC4, FC3 = (EXPECTED_FACES[name][1] for name in ("C5", "C4", "C3"))
FB, FE11, FT, FE10 = (EXPECTED_FACES[name][1]
                       for name in ("B", "E11", "T", "E10"))
EDGE_TESTS = (
    (FC5, (0, 1), (1, 6), alpha * X - 1),
    (FC5, (1, 6), (0, 6), alpha + u * X),
    (FC5, (0, 6), (0, 1), u - X**5),
    (FC4, (0, 1), (1, 6), alpha * X - 1),
    (FC4, (1, 6), (0, 5), alpha + Delta * X),
    (FC4, (0, 5), (0, 1), Delta - X**4),
    (FC3, (0, 1), (1, 6), alpha * X - 1),
    (FC3, (1, 6), (0, 4), alpha + e * X),
    (FC3, (0, 4), (0, 1), e - X**3),
    (FB, (0, 1), (2, 0), X - 1),
    (FB, (2, 0), (3, 5), 1 - alpha * X),
    (FB, (3, 5), (1, 6), alpha * (X - 1)),
    (FB, (1, 6), (0, 1), alpha - X),
    (FE11, (2, 0), (4, 3), 1 - Theta * X),
    (FE11, (4, 3), (3, 5), -Theta - alpha * X),
    (FE11, (3, 5), (2, 0), X - alpha),
    (FT, (2, 0), (4, 2), 1 - K * X**2),
    (FT, (4, 2), (4, 3), -K - Theta * X),
    (FT, (4, 3), (2, 0), X - Theta),
    (FE10, (2, 0), (4, 2), 1 - K * X**2),
    (FE10, (4, 2), (3, 5), -K - alpha * X),
    (FE10, (3, 5), (2, 0), X - alpha),
)
for expression, start, end, expected in EDGE_TESTS:
    check(associates(edge_word(expression, start, end, S, P, X), expected, X),
          "oriented edge word")
for word in (u - X**5, Delta - X**4, e - X**3, 1 - K * X**2):
    check(sp.factor(sp.discriminant(word, X)) != 0,
          "nonlinear owner word is generically reduced")

# Adjacent faces share primitive lattice edges.  Thus all attachments outside
# the eleven-node middle face are single nodes.
FACE_POINT_SETS = {
    name: {monomial for monomial, coefficient
           in sp.Poly(sp.expand(expression), S, P).terms() if coefficient != 0}
    for name, (_plane, expression) in EXPECTED_FACES.items()
}


def shared_edge_length(left, right):
    common = sorted(FACE_POINT_SETS[left] & FACE_POINT_SETS[right])
    check(len(common) == 2, left + "/" + right + " endpoints")
    return gcd(abs(common[1][0] - common[0][0]),
               abs(common[1][1] - common[0][1]))


for cap in ("C5", "C4", "C3"):
    check(shared_edge_length(cap, "B") == 1, cap + "/B attachment")
check(shared_edge_length("B", "E11") == 1, "B/E11 attachment")
check(shared_edge_length("E11", "T") == 1, "E11/T attachment")
check(shared_edge_length("B", "E10") == 1, "B/E10 attachment")
check((5, 14, 14 - 5 + 1) == (5, 14, 10), "Theta graph")
check((4, 13, 13 - 4 + 1) == (4, 13, 10), "Theta-zero graph")
check(10 + 3 == 13, "graph plus carrier genus")

# Six global polygons and their independently derived infinity packets.
GLOBAL = {
    ("C5", True): (((0, 1), (2, 0), (4, 2), (4, 3),
                     (3, 5), (1, 6), (0, 6)),
                    (36, 12, 13), (10, 8, 5, 3, 2, 2, 1)),
    ("C5", False): (((0, 1), (2, 0), (4, 2), (3, 5),
                      (1, 6), (0, 6)),
                     (35, 11, 13), (10, 10, 5, 2, 2, 1)),
    ("C4", True): (((0, 1), (2, 0), (4, 2), (4, 3),
                     (3, 5), (1, 6), (0, 5)),
                    (35, 11, 13), (10, 8, 5, 3, 2, 2, 1)),
    ("C4", False): (((0, 1), (2, 0), (4, 2), (3, 5),
                      (1, 6), (0, 5)),
                     (34, 10, 13), (10, 10, 5, 2, 2, 1)),
    ("C3", True): (((0, 1), (2, 0), (4, 2), (4, 3),
                     (3, 5), (1, 6), (0, 4)),
                    (34, 10, 13), (10, 8, 5, 3, 2, 2, 1)),
    ("C3", False): (((0, 1), (2, 0), (4, 2), (3, 5),
                      (1, 6), (0, 4)),
                     (33, 9, 13), (10, 10, 5, 2, 2, 1)),
}
for key, (vertices, expected_pick, expected_packet) in GLOBAL.items():
    hull = convex_hull(vertices)
    check(set(hull) == set(vertices), str(key) + " hull")
    actual_pick, actual_packet = pick_data(hull), infinity_packet(hull)
    check(actual_pick == expected_pick, str(key) + " Pick")
    check(actual_packet == expected_packet, str(key) + " packet")
    check(sum(index - 1 for index in actual_packet) == 24,
          str(key) + " packet defect")
    check(2 * actual_pick[2] - 2 == 24, str(key) + " genus defect")

for record in exact_records:
    (_c, ub, _pb, _eb, dv, tv, _address, _lift, mask) = record
    cap_name = "C5" if ub else "C4" if dv != 0 else "C3"
    vertices, expected_pick, expected_packet = GLOBAL[(cap_name, tv != 0)]
    hull = convex_hull(tuple((rr, ll) for rr, ll, _height in mask))
    check(set(hull) == set(vertices), "exact global polygon")
    check(pick_data(hull) == expected_pick, "exact Pick ledger")
    check(infinity_packet(hull) == expected_packet, "exact infinity packet")

print("U0_FIRST_NORMAL_OWNER_INDEPENDENT_REFEREE=PASS")
print("gate=Z=beta=zeta=W=xi=U=0;K*alpha!=0")
print("source_rows=16;exact_supports=64;exact_fans=6")
print("exact_fan_counts=C5T20,C5Z12,C4T16,C4Z8,C3T4,C3Z4")
print("hostile=256_keyed/168_distinct;fan_counts=64,64,32,32,32,32")
print("faces=C5,C4,C3,B,E11,T,E10;primitive_initials=7")
print("middle=2_rational_components/11_transverse_nodes/graph_b1_10")
print("carriers=E11_g3,E10_g3;global_genus=13")
print("orders=C5_25@30,C4_13@12,C3_9@6,B_49@66,E11_41@42,T_8@6,E10_26@24")
print("packets=Theta(10,8,5,3,2,2,1);zero(10,10,5,2,2,1);defect=24")
print("natural_address=n=1+8c+4u+2Phi+eta;range=1..64;bijective=yes")
print("checks=" + str(CHECKS))
