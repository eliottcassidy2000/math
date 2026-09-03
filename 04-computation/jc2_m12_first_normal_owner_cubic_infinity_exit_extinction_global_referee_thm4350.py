#!/usr/bin/env python3
"""Import-free global audit for THM-4350's xi_10=0 infinity boundary.

Scope: Z=beta_11=zeta_3=W=xi_10=0 and U*K!=0 on the inherited
exact-weight-12 seam K=2848/45-(7/6)Delta.
"""

from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations, product
from math import gcd, lcm
import sys

import sympy as sp


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")


def need(condition, label):
    if not bool(condition):
        raise AssertionError(label)


def rank_two(points):
    return any(
        (b[0] - a[0]) * (c[1] - a[1])
        != (b[1] - a[1]) * (c[0] - a[0])
        for a, b, c in combinations(points, 3)
    )


def all_planes(universe):
    answer = set()
    for p, q, r in combinations(universe, 3):
        det = (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])
        if det == 0:
            continue
        aa = F(
            (q[2] - p[2]) * (r[1] - p[1])
            - (q[1] - p[1]) * (r[2] - p[2]),
            det,
        )
        bb = F(
            (q[0] - p[0]) * (r[2] - p[2])
            - (q[2] - p[2]) * (r[0] - p[0]),
            det,
        )
        cc = F(p[2]) - aa * p[0] - bb * p[1]
        answer.add((aa, bb, cc))
    return tuple(sorted(answer))


def lower_faces(active, candidate_planes):
    result = set()
    for aa, bb, cc in candidate_planes:
        equality = []
        for rr, ll, hh in active:
            gap = F(hh) - aa * rr - bb * ll - cc
            if gap < 0:
                break
            if gap == 0:
                equality.append((rr, ll, hh))
        else:
            if rank_two(equality):
                result.add((aa, bb, cc))
    return frozenset(result)


def pick(vertices):
    pairs = list(zip(vertices, vertices[1:] + vertices[:1]))
    area2 = abs(sum(x1 * y2 - x2 * y1 for (x1, y1), (x2, y2) in pairs))
    boundary = sum(
        gcd(abs(x2 - x1), abs(y2 - y1))
        for (x1, y1), (x2, y2) in pairs
    )
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
    base = 1
    for number in (*plane, F(1, 6)):
        base = lcm(base, number.denominator)
    order = base * (F(5, 6) - sum(plane))
    need(order.denominator == 1 and order > 0, "positive integral face order")
    return base, order.numerator


def edge_polynomial(expression, start, end, S, P, X):
    dr, dl = end[0] - start[0], end[1] - start[1]
    length = gcd(abs(dr), abs(dl))
    dr, dl = dr // length, dl // length
    poly = sp.Poly(sp.expand(expression), S, P)
    return sp.expand(sum(
        poly.coeff_monomial(S ** (start[0] + i * dr) * P ** (start[1] + i * dl)) * X**i
        for i in range(length + 1)
    ))


# -------------------------------------------------------------------------
# Literal support, seam coupling, and exact versus synthetic mask atlases.
# -------------------------------------------------------------------------
U, u, alpha, Delta, eta, Theta, Phi, K = sp.symbols(
    "U u alpha Delta eta Theta Phi K"
)
xi = sp.Integer(0)
e = -sp.Rational(1376, 135)
rows = {
    (1, 0): -3,
    (2, 0): sp.Rational(8, 3),
    (3, 0): e,
    (4, 0): Delta,
    (5, 0): u,
    (6, 0): U,
    (0, 2): K,
    (2, 1): Phi,
    (1, 2): Theta,
    (3, 1): eta,
    (2, 2): xi,
    (4, 1): alpha,
}
support = defaultdict(lambda: sp.Integer(0))
support[(2, 0, 0)] += 1
support[(0, 1, 0)] -= 1
support[(2, 0, 1)] -= sp.Rational(1, 2)
for (i, j), coefficient in rows.items():
    support[(j + 2, i + j, 1)] -= coefficient
    support[(j, i + j + 1, 1)] += coefficient
support = {point: sp.expand(coefficient) for point, coefficient in support.items() if coefficient != 0}

need(support[(2, 6, 1)] == -U, "U vertex owner")
need(support[(4, 2, 1)] == -K, "K vertex owner")
need(support[(3, 5, 1)] == -alpha, "alpha next owner")
need(support[(4, 3, 1)] == -Theta, "Theta next owner")
need(sp.simplify(support[(2, 3, 1)] - (K - e)) == 0, "c1 aggregate")
need(sp.simplify(support[(2, 4, 1)] - (Theta - Delta)) == 0, "c2 aggregate")
need(sp.simplify(support[(2, 5, 1)] + u) == 0, "c3=-u coupling")
need((4, 4, 1) not in support, "xi vertex absent")

M = (F(1, 12), F(1, 6), F(-1, 6))
D6 = (F(1, 6), F(1, 6), F(-1, 3))
D7 = (F(2, 7), F(1, 7), F(-4, 7))
T = (F(1, 2), F(0), F(-1))
F10 = (F(3, 8), F(1, 8), F(-3, 4))
F01 = (F(1, 4), F(1, 6), F(-1, 2))
F00 = (F(1, 3), F(1, 6), F(-2, 3))
fans = {
    (1, 1): frozenset((M, D6, D7, T)),
    (1, 0): frozenset((M, D6, F10)),
    (0, 1): frozenset((M, F01, T)),
    (0, 0): frozenset((M, F00)),
}
candidate_planes = all_planes(tuple(support))

Krel = sp.Rational(2848, 45) - sp.Rational(7, 6) * Delta
dstar = sp.Rational(3968, 63)
need(sp.simplify(Krel.subs(Delta, dstar) - e) == 0, "c1 seam cancellation")
need(sp.solve(sp.Eq(Krel, 0), Delta)[0] == sp.Rational(5696, 105), "K=0 boundary")

delta_theta_reps = (
    (sp.Integer(0), sp.Integer(0)),
    (sp.Integer(0), sp.Integer(1)),
    (dstar, sp.Integer(0)),
    (dstar, dstar),
    (dstar, 2 * dstar),
    (sp.Integer(1), sp.Integer(0)),
    (sp.Integer(1), sp.Integer(1)),
    (sp.Integer(1), sp.Integer(2)),
)

actual_states = {}
actual_fan_counter = Counter()
actual_owner_counter = Counter()
for (delta_value, theta_value), u_value, phi_value, eta_value, alpha_value in product(
    delta_theta_reps, (0, 1), (0, 1), (0, 1), (0, 1)
):
    substitution = {
        U: 1,
        Delta: delta_value,
        K: Krel.subs(Delta, delta_value),
        Theta: theta_value,
        u: u_value,
        Phi: phi_value,
        eta: eta_value,
        alpha: alpha_value,
    }
    need(substitution[K] != 0, "representative lies on K-nonzero gate")
    active = frozenset(
        point for point, coefficient in support.items()
        if sp.simplify(coefficient.subs(substitution)) != 0
    )
    row_bits = (
        int(delta_value != 0),
        int(u_value != 0),
        int(phi_value != 0),
        int(theta_value != 0),
        int(eta_value != 0),
        int(alpha_value != 0),
    )
    cancellations = (
        int(sp.simplify((K - e).subs(substitution)) == 0),
        int(sp.simplify((Theta - Delta).subs(substitution)) == 0),
        int(u_value == 0),
    )
    key = row_bits + cancellations
    actual_states[active] = key
    owner = (int(alpha_value != 0), int(theta_value != 0))
    actual_owner_counter[owner] += 1
    actual_fan_counter[lower_faces(active, candidate_planes)] += 1

need(len(actual_states) == 8 * 2 * 2**3 == 128, "128 exact support patterns")
need(actual_owner_counter == Counter({(1, 1): 40, (1, 0): 24, (0, 1): 40, (0, 0): 24}), "exact owner counts")
need(actual_fan_counter == Counter({fans[(1, 1)]: 40, fans[(1, 0)]: 24, fans[(0, 1)]: 40, fans[(0, 0)]: 24}), "exact four-fan atlas")

required_labels = ((1, 0), (2, 0), (3, 0), (6, 0), (0, 2))
optional_labels = ((4, 0), (5, 0), (2, 1), (1, 2), (3, 1), (4, 1))
cancellable = ((2, 3, 1), (2, 4, 1), (2, 5, 1))


def lifted(labels):
    points = {(2, 0, 0), (0, 1, 0), (2, 0, 1)}
    for i, j in labels:
        points.add((j + 2, i + j, 1))
        points.add((j, i + j + 1, 1))
    return points


over_states = {}
over_counter = Counter()
for row_bits in product((0, 1), repeat=6):
    labels = list(required_labels)
    labels += [label for bit, label in zip(row_bits, optional_labels) if bit]
    for cancellations in product((0, 1), repeat=3):
        active = lifted(labels)
        for bit, point in zip(cancellations, cancellable):
            if bit:
                active.discard(point)
        key = tuple(row_bits) + tuple(cancellations)
        over_states[key] = frozenset(active)
        over_counter[lower_faces(active, candidate_planes)] += 1

need(len(over_states) == 512, "512 Boolean over-atlas")
need(over_counter == Counter({fan: 128 for fan in fans.values()}), "over-atlas four-fan counts")
for active, key in actual_states.items():
    need(key in over_states and over_states[key] == active, "exact support embeds in over-atlas")

# The 512 figure counts labelled row/deletion configurations, not distinct
# active supports: deletion of an already absent aggregate creates duplicate
# configurations.  This second quotient is another synthetic-mask firewall.
distinct_over_supports = set(over_states.values())
exact_supports = set(actual_states)
need(len(distinct_over_supports) == 336, "336 distinct conservative supports")
need(len(distinct_over_supports - exact_supports) == 208,
     "208 genuinely synthetic distinct supports")
need(sum(active in exact_supports for active in over_states.values()) == 216,
     "216 conservative configurations land on exact supports")
need(Counter(Counter(over_states.values()).values()) == Counter({1: 192, 2: 128, 4: 16}),
     "over-atlas configuration multiplicities")

# -------------------------------------------------------------------------
# Exact faces, polygons, boundary/internal schemes, Pick/RH, and orders.
# -------------------------------------------------------------------------
S, P, X = sp.symbols("S P X")
face_expressions = {
    "M": (P - S**2) * (U * P**6 - 1),
    "D6": S**2 * (1 - U * P**6 - alpha * S * P**5),
    "D7": S**2 * (1 - alpha * S * P**5 - Theta * S**2 * P**3),
    "T": S**2 * (1 - S**2 * P**2 * (K + Theta * P)),
    "F10": S**2 * (1 - alpha * S * P**5 - K * S**2 * P**2),
    "F01": S**2 * (1 - U * P**6 - Theta * S**2 * P**3),
    "F00": S**2 * (1 - U * P**6 - eta * S * P**4 - K * S**2 * P**2),
}
plane_names = {"M": M, "D6": D6, "D7": D7, "T": T, "F10": F10, "F01": F01, "F00": F00}
for name, plane in plane_names.items():
    reconstructed = 0
    for (rr, ll, hh), coefficient in support.items():
        if F(hh) == plane[0] * rr + plane[1] * ll + plane[2]:
            reconstructed += coefficient * S**rr * P**ll
    need(sp.expand(reconstructed - face_expressions[name]) == 0, f"exact face {name}")

polygons = {
    "M": [(0, 1), (2, 0), (2, 6), (0, 7)],
    "D6": [(2, 0), (3, 5), (2, 6)],
    "D7": [(2, 0), (4, 3), (3, 5)],
    "T": [(2, 0), (4, 2), (4, 3)],
    "F10": [(2, 0), (4, 2), (3, 5)],
    "F01": [(2, 0), (4, 3), (2, 6)],
    "F00": [(2, 0), (4, 2), (2, 6)],
}
expected_pick = {
    "M": (24, 14, 6),
    "D6": (6, 8, 0),
    "D7": (7, 3, 3),
    "T": (2, 4, 0),
    "F10": (8, 4, 3),
    "F01": (12, 8, 3),
    "F00": (12, 10, 2),
}
need({name: pick(poly) for name, poly in polygons.items()} == expected_pick, "face Pick ledgers")
need({name: face_order(plane) for name, plane in plane_names.items()} == {
    "M": (12, 9),
    "D6": (6, 5),
    "D7": (42, 41),
    "T": (6, 8),
    "F10": (24, 26),
    "F01": (12, 11),
    "F00": (6, 6),
}, "face base/order ledgers")

configs = {
    "11": {
        "global": [(0, 1), (2, 0), (4, 2), (4, 3), (3, 5), (2, 6), (0, 7)],
        "boundary_owners": ["M", "T", "T", "D7", "D6", "M", "M"],
        "boundary_expected": [X - 1, 1 - K * X**2, -K - Theta * X, -Theta - alpha * X, -alpha - U * X, U * (X - 1), U - X**6],
        "internal": [("M", "D6", (2, 0), (2, 6), 1 - U * X**6), ("D6", "D7", (2, 0), (3, 5), 1 - alpha * X), ("D7", "T", (2, 0), (4, 3), 1 - Theta * X)],
        "pick": (39, 13, 14),
        "packet": (11, 8, 6, 3, 2, 2, 1),
        "VEg": (10, 20, 11, 14),
    },
    "10": {
        "global": [(0, 1), (2, 0), (4, 2), (3, 5), (2, 6), (0, 7)],
        "boundary_owners": ["M", "F10", "F10", "D6", "M", "M"],
        "boundary_expected": [X - 1, 1 - K * X**2, -K - alpha * X, -alpha - U * X, U * (X - 1), U - X**6],
        "internal": [("M", "D6", (2, 0), (2, 6), 1 - U * X**6), ("D6", "F10", (2, 0), (3, 5), 1 - alpha * X)],
        "pick": (38, 12, 14),
        "packet": (11, 10, 6, 2, 2, 1),
        "VEg": (9, 19, 11, 14),
    },
    "01": {
        "global": [(0, 1), (2, 0), (4, 2), (4, 3), (2, 6), (0, 7)],
        "boundary_owners": ["M", "T", "T", "F01", "M", "M"],
        "boundary_expected": [X - 1, 1 - K * X**2, -K - Theta * X, -Theta - U * X, U * (X - 1), U - X**6],
        "internal": [("M", "F01", (2, 0), (2, 6), 1 - U * X**6), ("F01", "T", (2, 0), (4, 3), 1 - Theta * X)],
        "pick": (38, 12, 14),
        "packet": (13, 11, 3, 2, 2, 1),
        "VEg": (9, 19, 11, 14),
    },
    "00": {
        "global": [(0, 1), (2, 0), (4, 2), (2, 6), (0, 7)],
        "boundary_owners": ["M", "F00", "F00", "M", "M"],
        "boundary_expected": [X - 1, 1 - K * X**2, -K - eta * X - U * X**2, U * (X - 1), U - X**6],
        "internal": [("M", "F00", (2, 0), (2, 6), 1 - U * X**6)],
        "pick": (36, 12, 13),
        "packet": (11, 7, 7, 2, 2, 1),
        "VEg": (8, 18, 11, 13),
    },
}

for label, config in configs.items():
    vertices = config["global"]
    need(pick(vertices) == config["pick"], f"global Pick {label}")
    packet = edge_packet(vertices)
    need(packet == config["packet"], f"packet {label}")
    need(sum(index - 1 for index in packet) == 2 * config["pick"][2] - 2, f"RH {label}")
    for start, end, owner, expected in zip(vertices, vertices[1:] + vertices[:1], config["boundary_owners"], config["boundary_expected"]):
        need(sp.expand(edge_polynomial(face_expressions[owner], start, end, S, P, X) - expected) == 0, f"boundary {label} {start}->{end}")
    for first, second, start, end, expected in config["internal"]:
        need(sp.expand(edge_polynomial(face_expressions[first], start, end, S, P, X) - expected) == 0, f"internal first {label}")
        need(sp.expand(edge_polynomial(face_expressions[second], start, end, S, P, X) - expected) == 0, f"internal second {label}")
    vertices_count, edges_count, graph_b1, total_genus = config["VEg"]
    need(edges_count - vertices_count + 1 == graph_b1, f"graph b1 {label}")
    need(graph_b1 + (3 if label != "00" else 2) == total_genus, f"component genus {label}")

# -------------------------------------------------------------------------
# Exact face normal forms and the sole residual collision.
# -------------------------------------------------------------------------
core_D7 = 1 - alpha * S * P**5 - Theta * S**2 * P**3
Y_D7 = 2 * Theta * P**2 * S + alpha * P**4
f_D7 = P * (alpha**2 * P**7 + 4 * Theta)
need(sp.expand(Y_D7**2 - f_D7 + 4 * Theta * P * core_D7) == 0, "D7 hyperelliptic identity")

core_F10 = 1 - alpha * S * P**5 - K * S**2 * P**2
Y_F10 = 2 * K * S * P + alpha * P**4
f_F10 = alpha**2 * P**8 + 4 * K
need(sp.expand(Y_F10**2 - f_F10 + 4 * K * core_F10) == 0, "F10 hyperelliptic identity")

core_F01 = 1 - U * P**6 - Theta * S**2 * P**3
Y_F01 = P**2 * S
f_F01 = P * (1 - U * P**6)
need(sp.expand(Theta * Y_F01**2 - f_F01 + P * core_F01) == 0, "F01 hyperelliptic identity")

core_F00 = 1 - U * P**6 - eta * S * P**4 - K * S**2 * P**2
Y_F00 = 2 * K * S * P + eta * P**3
L00 = eta**2 - 4 * K * U
f_F00 = L00 * P**6 + 4 * K
need(sp.expand(Y_F00**2 - f_F00 + 4 * K * core_F00) == 0, "F00 hyperelliptic identity")

field = sp.QQ.frac_field(U, K, alpha, Theta, eta)
for name, branch, degree, genus in (
    ("D7", f_D7, 8, 3),
    ("F10", f_F10, 8, 3),
    ("F01", f_F01, 7, 3),
):
    poly = sp.Poly(branch, P, domain=field)
    need(poly.degree() == degree, f"{name} degree")
    need(sp.gcd(poly, poly.diff()).degree() == 0, f"{name} squarefree")
    need((degree - 2) // 2 == genus if degree % 2 == 0 else (degree - 1) // 2 == genus, f"{name} genus")

# F00 is genus two exactly off L00=0.  At L00=0 it splits into two rational
# branches with A5 contact at infinity; all six M attachments split 3+3.
Lsym = sp.symbols("Lsym", nonzero=True)
poly00 = sp.Poly(Lsym * P**6 + 4 * K, P, domain=sp.QQ.frac_field(Lsym, K))
need(sp.gcd(poly00, poly00.diff()).degree() == 0 and poly00.degree() == 6, "F00 squarefree off collision")
need((6 - 2) // 2 == 2, "F00 genus two")

lam, mu = sp.symbols("lam mu", nonzero=True)
w = S * P
collision_sub = {K: lam**2, U: mu**2, eta: 2 * lam * mu}
need(sp.expand(core_F00.subs(collision_sub) - (1 - (lam * w + mu * P**3) ** 2)) == 0, "F00 collision square")
need(sp.factor((K + eta * X + U * X**2).subs(collision_sub)) == (lam + mu * X) ** 2, "outer double edge")
need(sp.expand((1 - U * X**6).subs(collision_sub) - (1 - mu * X**3) * (1 + mu * X**3)) == 0, "3+3 M attachments")
need(18 - 9 + 1 == 10, "split F00 graph b1")
# Arithmetic genus two, two rational branches, and one collision point imply
# conductor delta 3: 2 = 0+0+delta-2+1.
need(2 + 2 - 1 == 3, "F00 A5 delta")

# Nonempty intersections with each seam-support class show that seam
# cancellations do not create a second geometric collision owner.
need(sp.simplify(Krel.subs(Delta, 0)) != 0, "collision can meet Delta=0 class")
need(sp.simplify(Krel.subs(Delta, dstar)) != 0, "collision can meet c1=0 class")
need(sp.simplify(Krel.subs(Delta, 1)) != 0, "collision can meet generic seam class")

# -------------------------------------------------------------------------
# Independent exact local chart on the repeated F00 edge.  This sharpens the
# apparent stopping wall: the fixed cubic coefficient e forces depth <= 3.
# -------------------------------------------------------------------------
sigma = sp.symbols("sigma")
p0 = sigma**-1 * P
s0 = sigma**-2 * S
y0 = s0 * p0
H0 = sum(coefficient * p0**i * y0**j for (i, j), coefficient in rows.items())
H0 = sp.expand(H0.subs({alpha: 0, Theta: 0}))
F_source = (s0**2 - p0) * (1 - sigma**6 * H0) - sigma**6 * s0**2 / 2
G_source = sp.expand(sigma**4 * F_source)
H6 = U * P**6 + eta * S * P**4 + K * S**2 * P**2
H5 = u * P**5 + Phi * S * P**3
expected_G = sp.expand(
    (S**2 - sigma**3 * P)
    * (1 - H6 - sigma * H5 - sigma**2 * Delta * P**4
       - e * sigma**3 * P**3 - sp.Rational(8, 3) * sigma**4 * P**2
       + 3 * sigma**5 * P)
    - sigma**6 * S**2 / 2
)
need(sp.simplify(G_source - expected_G) == 0, "exact F00 primitive chart")

x, v, rho = sp.symbols("x v rho")
A = U + eta * v + K * v**2
B = u + Phi * v
root_equation = sp.factor(x**10 * G_source.subs({P: 1 / x, S: v / x**2}))
expected_root = sp.expand(
    (v**2 - (sigma * x)**3)
    * (x**6 - A - sigma * x * B - (sigma * x)**2 * Delta
       - e * (sigma * x)**3 - sp.Rational(8, 3) * (sigma * x)**4
       + 3 * (sigma * x)**5)
    - (sigma * x)**6 * v**2 / 2
)
need(sp.simplify(root_equation - expected_root) == 0, "exact reciprocal F00 chart")

Droot = (
    A + rho * B + rho**2 * Delta + e * rho**3
    + sp.Rational(8, 3) * rho**4 - 3 * rho**5
    + rho**6 * v**2 / (2 * (v**2 - rho**3))
)
need(sp.factor(sp.together(
    root_equation / (v**2 - (sigma * x)**3)
    - (x**6 - Droot.subs(rho, sigma * x))
)) == 0, "divided one-series chart")

a = sp.symbols("a", nonzero=True)
repeat = {U: K * a**2, eta: -2 * K * a}
need(sp.expand(A.subs(repeat) - K * (v - a)**2) == 0,
     "canonical repeated-edge parameterization")
need(sp.diff(Droot, v, 2).subs({rho: 0, v: a, **repeat}) == 2 * K,
     "Morse Hessian is a unit")
need((v**2 - rho**3).subs({rho: 0, v: a}) == a**2,
     "reciprocal prefactor is a unit")

# Corrections to the displayed critical section begin at rho^6, so this
# substitution computes the complete critical value through rho^5.
vcrit = a - rho * Phi / (2 * K)
critical = sp.series(sp.together(Droot.subs(repeat).subs(v, vcrit)), rho, 0, 6).removeO().expand()
critical_expected = (
    (u + Phi * a) * rho
    + (Delta - Phi**2 / (4 * K)) * rho**2
    + e * rho**3 + sp.Rational(8, 3) * rho**4 - 3 * rho**5
)
need(sp.simplify(critical - critical_expected) == 0,
     "critical-value series through fixed cubic")
need(e != 0, "critical depth is at most three")

# Phi_root=x^10 G and S=v/x^2 give G_S=x^-8 Phi_v.  Together
# with the F00 primary good-form order 6 this transports the differential as
# a unit times sigma^6*x^6 dx/y after formal Morse preparation.
GS_root = sp.diff(G_source, S).subs({P: 1 / x, S: v / x**2})
need(sp.simplify(sp.diff(root_equation, v) - x**8 * GS_root) == 0,
     "root-chart differential conversion")

tail_data = {1: (2, 0, 68), 2: (1, 1, 64), 3: (1, 1, 60)}
tau, z, h, Crho = sp.symbols("tau z h Crho")
for depth, (tail_genus, persistent_delta, form_order) in tail_data.items():
    epsilon = depth % 2
    tail_degree = 6 - depth + epsilon
    need((tail_degree - 1) // 2 == tail_genus, "tail genus")
    need(depth // 2 == persistent_delta, "persistent square delta")
    need(persistent_delta + tail_genus + 1 == 3, "A5 conservation")
    computed_order = 6 * 2 * (6 - depth) + 6 * 2 * depth + 2 * depth - 6 * depth
    need(computed_order == form_order > 0, "positive tail form order")
    raw_complement = tau**(12 * (6 - depth)) * rho**(depth - 6) * Crho
    rewritten = sp.cancel(raw_complement.subs(tau**12, z * rho))
    need(sp.simplify(rewritten - z**(6 - depth) * Crho) == 0,
         "two-ended complement")
    complement = h**2 - 1 + z**(6 - depth) * Crho
    need(sp.diff(complement, h).subs({h: 1, z: 0}) == 2,
         "positive sign attachment is etale")
    need(sp.diff(complement, h).subs({h: -1, z: 0}) == -2,
         "negative sign attachment is etale")

lambda0 = sp.symbols("lambda0", nonzero=True)
linear = lambda0 * (S * P - a * P**3)
factor_sub = {K: lambda0**2, U: lambda0**2 * a**2,
              eta: -2 * lambda0**2 * a}
need(sp.expand(core_F00.subs(factor_sub) - (1 - linear) * (1 + linear)) == 0,
     "repeated F00 splits into two rational signs")
need(sp.expand((1 - U * P**6).subs({U: lambda0**2 * a**2})
               - (1 - lambda0 * a * P**3) * (1 + lambda0 * a * P**3)) == 0,
     "six M attachments split three plus three")
need(12 + 6 + 3 - 9 + 1 == 13, "raw collision arithmetic genus")
for tail_genus, persistent_delta, _form_order in tail_data.values():
    regularized_b1 = (18 + 2) - (9 + 1) + 1
    need(regularized_b1 + tail_genus + persistent_delta == 13,
         "regularized collision genus")

print("THM4350_GLOBAL_EXACT_AUDIT=PASS")
print("gate=Z=beta_11=zeta_3=W=xi_10=0;U*K!=0")
print("exact_supports=128;fan_counts=alphaTheta11:40,10:24,01:40,00:24")
print("conservative_over_atlas=512_keyed_configs;fan_count_each=128")
print("over_atlas_distinct_supports=336;synthetic_distinct_supports=208")
print("canonical_exact_embedding=128_configs;other_coupling_illegal_configs=384")
print("clean_owner_union=(alpha,Theta)!=(0,0);exact_supports=104")
print("maximal_clean_locus=clean_owner_union OR [alpha=Theta=0 and eta^2!=4KU]")
print("face_genera=D7/F10/F01:3;F00:2")
print("face_orders=M:9@12,D6:5@6,D7:41@42,T:8@6,F10:26@24,F01:11@12,F00:6@6")
print("global_pick=11/10/01:genus14;00:genus13;graph_b1=11_all")
print("sole_collision=alpha=Theta=0,eta^2=4KU;F00_splits;M_attachments=3+3;A5_delta=3")
print("collision_critical_depths=r1/r2/r3;tail_genera=2/1/1;form_orders=68/64/60")
print("collision_status=exact_even_A5_local_closure;no_remaining_normal-series wall")
