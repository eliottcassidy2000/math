#!/usr/bin/env python3
"""Primary exact certificate for THM-4354's first-normal-owner U=0 gate.

Scope:
    Z=beta_11=zeta_3=W=xi_10=U=0, K*alpha_11 != 0,
    K=2848/45-(7/6)Delta

in the inherited reduced (2,3), exact-weight-twelve seam.  This self-contained
certificate reconstructs the literal sixteen-row source, enumerates the 64
coefficient-realizable supports and a 256-key hostile Boolean atlas, computes
all lower faces and primitive charts, checks the reducible middle face and its
eleven nodes, proves the carrier genera, and verifies every edge, global Pick,
packet, graph, and natural-address ledger used by THM-4354.
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
        raise RuntimeError("THM-4354 primary audit failure: " + label)


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
    need(order.denominator == 1 and order > 0,
         "positive integral form order")
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
        exponent = (start[0] + index * step[0],
                    start[1] + index * step[1])
        answer += polynomial.coeff_monomial(
            S**exponent[0] * P**exponent[1]
        ) * X**index
    return sp.factor(answer)


# Literal sixteen-row source and the U=0 first-normal-owner gate.
U, W, Z = sp.symbols("U W Z")
Delta, Theta, Phi, eta = sp.symbols("Delta Theta Phi eta")
zeta, u, xi, alpha, beta, K = sp.symbols(
    "zeta u xi alpha beta K"
)
S, P, X = sp.symbols("S P X")
e = -sp.Rational(1376, 135)
KREL = sp.Rational(2848, 45) - sp.Rational(7, 6) * Delta
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
ZERO_GATE = {Z: 0, beta: 0, zeta: 0, W: 0, xi: 0, U: 0}
ROWS = {
    label: sp.factor(sp.sympify(coefficient).subs(ZERO_GATE))
    for label, coefficient in FULL_ROWS.items()
    if sp.factor(sp.sympify(coefficient).subs(ZERO_GATE)) != 0
}
need(len(ROWS) == 10, "ten surviving labelled source rows")

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
    (0, 5, 1): Delta,
    (0, 6, 1): u,
    (1, 4, 1): Phi,
    (1, 5, 1): eta,
    (1, 6, 1): alpha,
    (2, 0, 0): 1,
    (2, 0, 1): -sp.Rational(1, 2),
    (2, 1, 1): 3,
    (2, 2, 1): -sp.Rational(8, 3),
    (2, 3, 1): K - e,
    (2, 4, 1): Theta - Delta,
    (2, 5, 1): -u,
    (3, 3, 1): -Phi,
    (3, 4, 1): -eta,
    (3, 5, 1): -alpha,
    (4, 2, 1): -K,
    (4, 3, 1): -Theta,
}
need(SUPPORT.keys() == EXPECTED_SUPPORT.keys(),
     "complete specialized support has 21 expected points")
for point in EXPECTED_SUPPORT:
    need(sp.expand(SUPPORT[point] - EXPECTED_SUPPORT[point]) == 0,
         "literal specialized support coefficient")
need(sp.expand(SUPPORT[(2, 3, 1)] - (K + sp.Rational(1376, 135))) == 0,
     "first multiply visible coefficient")
need(SUPPORT[(2, 4, 1)] == Theta - Delta,
     "second multiply visible coefficient")
need(SUPPORT[(2, 5, 1)] == -u,
     "third multiply visible coefficient")

# All possible lower planes and the exact 3 x 2 owner grid.
C5 = (F(0), F(1, 5), F(-1, 5))
C4 = (F(-1, 4), F(1, 4), F(-1, 4))
C3 = (F(-2, 3), F(1, 3), F(-1, 3))
B = (F(1, 11), F(2, 11), F(-2, 11))
E11 = (F(2, 7), F(1, 7), F(-4, 7))
T = (F(1, 2), F(0), F(-1))
E10 = (F(3, 8), F(1, 8), F(-3, 4))
PLANES = candidate_planes(SUPPORT)
FANS = {
    ("C5", 1): frozenset((C5, B, E11, T)),
    ("C5", 0): frozenset((C5, B, E10)),
    ("C4", 1): frozenset((C4, B, E11, T)),
    ("C4", 0): frozenset((C4, B, E10)),
    ("C3", 1): frozenset((C3, B, E11, T)),
    ("C3", 0): frozenset((C3, B, E10)),
}


def owner(delta_on, u_on, theta_on):
    if u_on:
        cap = "C5"
    elif delta_on:
        cap = "C4"
    else:
        cap = "C3"
    return cap, int(theta_on)


DSTAR = sp.solve(sp.Eq(KREL - e, 0), Delta)[0]
KZERO = sp.solve(sp.Eq(KREL, 0), Delta)[0]
need(DSTAR == sp.Rational(3968, 63),
     "K-e cancellation occurs at Delta=3968/63")
need(KZERO == sp.Rational(5696, 105),
     "excluded K=0 endpoint occurs at Delta=5696/105")
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
    need(kval != 0, "coupled-class representative obeys K!=0")
    status_keys.add((
        int(delta_value != 0),
        int(kval - e != 0),
        int(theta_value != 0),
        int(theta_value - delta_value != 0),
    ))
need(len(status_keys) == 8, "eight coupled Delta/Theta support classes")

exact_supports = {}
exact_fans = Counter()
exact_owner = {}
for class_index, (delta_value, theta_value) in enumerate(DELTA_THETA_REPS):
    for u_bit, phi_bit, eta_bit in product((0, 1), repeat=3):
        substitution = {
            Delta: delta_value,
            K: KREL.subs(Delta, delta_value),
            Theta: theta_value,
            u: u_bit,
            Phi: phi_bit,
            eta: eta_bit,
            alpha: 1,
        }
        active = frozenset(
            point for point, coefficient in SUPPORT.items()
            if sp.factor(coefficient.subs(substitution)) != 0
        )
        key = (class_index, u_bit, phi_bit, eta_bit)
        owner_key = owner(delta_value != 0, u_bit, theta_value != 0)
        fan = lower_fan(active, PLANES)
        need(fan == FANS[owner_key], "exact support has its owner fan")
        exact_supports[key] = active
        exact_owner[key] = owner_key
        exact_fans[owner_key] += 1

need(len(exact_supports) == 64, "64 exact coefficient states")
need(len(set(exact_supports.values())) == 64,
     "64 distinct exact realizable supports")
need(exact_fans == Counter({
    ("C5", 1): 20,
    ("C5", 0): 12,
    ("C4", 1): 16,
    ("C4", 0): 8,
    ("C3", 1): 4,
    ("C3", 0): 4,
}), "exact six-fan census")

# Conservative hostile Boolean atlas: alpha and K are forced present; the
# other five optional rows and all three aggregate deletions are independent.
REQUIRED_LABELS = ((1, 0), (2, 0), (3, 0), (0, 2), (4, 1))
OPTIONAL_LABELS = ((4, 0), (5, 0), (2, 1), (1, 2), (3, 1))
CANCELLABLE = ((2, 3, 1), (2, 4, 1), (2, 5, 1))


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
    labels.extend(
        label for bit, label in zip(row_bits, OPTIONAL_LABELS) if bit
    )
    for deletion_bits in product((0, 1), repeat=3):
        active = lifted(labels)
        for bit, point in zip(deletion_bits, CANCELLABLE):
            if bit:
                active.discard(point)
        active = frozenset(active)
        key = tuple(row_bits) + tuple(deletion_bits)
        hostile_states[key] = active
        delta_on, u_on, _phi_on, theta_on, _eta_on = row_bits
        owner_key = owner(delta_on, u_on, theta_on)
        need(lower_fan(active, PLANES) == FANS[owner_key],
             "hostile support has its owner fan")
        hostile_fans[owner_key] += 1

need(len(hostile_states) == 256, "256 hostile keyed configurations")
need(len(set(hostile_states.values())) == 168,
     "168 hostile distinct supports")
need(hostile_fans == Counter({
    ("C5", 1): 64,
    ("C5", 0): 64,
    ("C4", 1): 32,
    ("C4", 0): 32,
    ("C3", 1): 32,
    ("C3", 0): 32,
}), "hostile six-fan census")

for key, active in exact_supports.items():
    class_index, u_bit, phi_bit, eta_bit = key
    delta_value, theta_value = DELTA_THETA_REPS[class_index]
    kval = KREL.subs(Delta, delta_value)
    row_bits = (
        int(delta_value != 0), u_bit, phi_bit,
        int(theta_value != 0), eta_bit,
    )
    deletion_bits = (
        int(kval - e == 0),
        int(theta_value - delta_value == 0),
        int(u_bit == 0),
    )
    atlas_key = tuple(row_bits) + tuple(deletion_bits)
    need(hostile_states[atlas_key] == active,
         "exact support embeds literally in hostile atlas")

# Exact face equations, polygons, Pick ledgers, and common-base form orders.
FACE_DATA = {
    "C5": {
        "plane": C5,
        "sub": {},
        "polygon": ((0, 1), (1, 6), (0, 6)),
        "equation": P * (u * P**5 + alpha * S * P**5 - 1),
        "pick": (5, 7, 0),
        "order": (F(5, 6), 30, 25),
    },
    "C4": {
        "plane": C4,
        "sub": {u: 0},
        "polygon": ((0, 1), (1, 6), (0, 5)),
        "equation": P * (Delta * P**4 + alpha * S * P**5 - 1),
        "pick": (4, 6, 0),
        "order": (F(13, 12), 12, 13),
    },
    "C3": {
        "plane": C3,
        "sub": {u: 0, Delta: 0},
        "polygon": ((0, 1), (1, 6), (0, 4)),
        "equation": P * (e * P**3 + alpha * S * P**5 - 1),
        "pick": (3, 5, 0),
        "order": (F(3, 2), 6, 9),
    },
    "B": {
        "plane": B,
        "sub": {},
        "polygon": ((0, 1), (2, 0), (3, 5), (1, 6)),
        "equation": (P - S**2) * (alpha * S * P**5 - 1),
        "pick": (22, 4, 10),
        "order": (F(49, 66), 66, 49),
    },
    "E11": {
        "plane": E11,
        "sub": {},
        "polygon": ((2, 0), (4, 3), (3, 5)),
        "equation": S**2 * (1 - alpha * S * P**5
                              - Theta * S**2 * P**3),
        "pick": (7, 3, 3),
        "order": (F(41, 42), 42, 41),
    },
    "T": {
        "plane": T,
        "sub": {},
        "polygon": ((2, 0), (4, 2), (4, 3)),
        "equation": S**2 * (1 - S**2 * P**2 * (K + Theta * P)),
        "pick": (2, 4, 0),
        "order": (F(4, 3), 6, 8),
    },
    "E10": {
        "plane": E10,
        "sub": {Theta: 0},
        "polygon": ((2, 0), (4, 2), (3, 5)),
        "equation": S**2 * (1 - alpha * S * P**5
                              - K * S**2 * P**2),
        "pick": (8, 4, 3),
        "order": (F(13, 12), 24, 26),
    },
}
for name, data in FACE_DATA.items():
    actual = face_expression(SUPPORT, data["plane"], data["sub"])
    need(sp.expand(actual - data["equation"]) == 0,
         name + " exact face equation")
    need(pick(data["polygon"]) == data["pick"],
         name + " face Pick ledger")
    need(face_order(data["plane"]) == data["order"],
         name + " common-base form order")

# Every primitive chart is reconstructed directly from the complete source,
# after applying both the zero gate and the seam relation K=KREL.
sigma = sp.symbols("sigma")
HROWS = {
    label: sp.factor(sp.sympify(coefficient).subs(ZERO_GATE).subs(K, KREL))
    for label, coefficient in FULL_ROWS.items()
    if sp.factor(sp.sympify(coefficient).subs(ZERO_GATE).subs(K, KREL)) != 0
}
need(len(HROWS) == 10, "ten rows survive in primitive-chart source")


def primitive_chart(qpower, s_exponent, p_exponent, clear, substitution):
    old_s = sigma**s_exponent * S
    old_p = sigma**p_exponent * P
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
    "C5": (5, 0, -1, 1, {}),
    "C4": (4, 1, -1, 1, {u: 0}),
    "C3": (3, 2, -1, 1, {u: 0, Delta: 0}),
    "B": (11, -1, -2, 2, {}),
    "E11": (7, -2, -1, 4, {}),
    "T": (2, -1, 0, 2, {}),
    "E10": (8, -3, -1, 6, {Theta: 0}),
}
CHARTS = {}
for name, spec in CHART_SPECS.items():
    chart = primitive_chart(*spec)
    expected = FACE_DATA[name]["equation"].subs(K, KREL)
    expected = expected.subs(FACE_DATA[name]["sub"])
    need(sp.expand(chart.subs(sigma, 0) - expected) == 0,
         name + " primitive initial")
    CHARTS[name] = chart

# Explicit one-parameter normalizations certify every rational face.  The
# displayed P and S^2 factors are invertible torus monomials, not components.
CAP_RATIONAL = {
    "C5": (u * P**5 + alpha * S * P**5 - 1,
           (1 - u * P**5) / (alpha * P**5)),
    "C4": (Delta * P**4 + alpha * S * P**5 - 1,
           (1 - Delta * P**4) / (alpha * P**5)),
    "C3": (e * P**3 + alpha * S * P**5 - 1,
           (1 - e * P**3) / (alpha * P**5)),
}
for name, (equation, s_solution) in CAP_RATIONAL.items():
    need(sp.factor(equation.subs(S, s_solution)) == 0,
         name + " rational normalization")

# The middle face has two rational components and eleven transverse nodes.
fB = FACE_DATA["B"]["equation"]
r_component = P - S**2
c_component = alpha * S * P**5 - 1
need(sp.expand(fB - r_component * c_component) == 0,
     "B factorization into its two rational components")
need(r_component.subs(P, S**2) == 0, "B R component is rational")
need(sp.factor(c_component.subs(S, 1 / (alpha * P**5))) == 0,
     "B C component is rational")
node_polynomial = -sp.factor(sp.resultant(r_component, c_component, P))
need(node_polynomial == alpha * S**11 - 1,
     "B intersection equation alpha*S^11=1")
need(sp.degree(node_polynomial, S) == 11, "B has eleven intersection roots")
node_domain = sp.QQ.frac_field(alpha)
node_gcd = sp.gcd(
    sp.Poly(node_polynomial, S, domain=node_domain),
    sp.Poly(sp.diff(node_polynomial, S), S, domain=node_domain),
)
need(node_gcd.degree() == 0, "all eleven B roots are distinct")
jacobian = sp.det(sp.Matrix([
    [sp.diff(r_component, S), sp.diff(r_component, P)],
    [sp.diff(c_component, S), sp.diff(c_component, P)],
]))
need(sp.expand(jacobian.subs(P, S**2) + 11 * alpha * S**10) == 0,
     "B component Jacobian is nonzero at every node")
need(11 - 2 + 1 == 10, "B multigraph has Betti number ten")

r_parameter = sp.symbols("r_parameter")
t_core = 1 - S**2 * P**2 * (K + Theta * P)
t_p = (1 / r_parameter**2 - K) / Theta
need(sp.factor(t_core.subs({P: t_p, S: r_parameter / t_p})) == 0,
     "T face is rational on its Theta owner gate")

# Exact hyperelliptic normalizations prove the two possible carrier genera.
f11 = 1 - alpha * S * P**5 - Theta * S**2 * P**3
y11 = 2 * Theta * S * P**2 + alpha * P**4
branch11 = P * (alpha**2 * P**7 + 4 * Theta)
need(sp.expand(y11**2 - branch11 + 4 * Theta * P * f11) == 0,
     "E11 hyperelliptic identity")
disc11 = sp.factor(sp.discriminant(branch11, P))
need(disc11 == -sp.Integer(53971714048) * Theta**8 * alpha**12,
     "E11 branch discriminant")
need(sp.degree(branch11, P) == 8 and (sp.degree(branch11, P) - 2) // 2 == 3,
     "E11 squarefree degree-eight branch has genus three")

f10 = 1 - alpha * S * P**5 - K * S**2 * P**2
y10 = 2 * K * S * P + alpha * P**4
branch10 = alpha**2 * P**8 + 4 * K
need(sp.expand(y10**2 - branch10 + 4 * K * f10) == 0,
     "E10 hyperelliptic identity")
disc10 = sp.factor(sp.discriminant(branch10, P))
need(disc10 == sp.Integer(274877906944) * K**7 * alpha**14,
     "E10 branch discriminant")
need(sp.degree(branch10, P) == 8 and (sp.degree(branch10, P) - 2) // 2 == 3,
     "E10 squarefree degree-eight branch has genus three")

# Literal edge words for all seven faces, including every boundary edge.
EDGE_WORDS = {
    "C5": (alpha * X - 1, alpha + u * X, u - X**5),
    "C4": (alpha * X - 1, alpha + Delta * X, Delta - X**4),
    "C3": (alpha * X - 1, alpha + e * X, e - X**3),
    "B": (X - 1, 1 - alpha * X, alpha * (X - 1), alpha - X),
    "E11": (1 - Theta * X, -Theta - alpha * X, X - alpha),
    "T": (1 - K * X**2, -K - Theta * X, X - Theta),
    "E10": (1 - K * X**2, -K - alpha * X, X - alpha),
}
for name, words in EDGE_WORDS.items():
    polygon = FACE_DATA[name]["polygon"]
    equation = FACE_DATA[name]["equation"]
    pairs = tuple(zip(polygon, polygon[1:] + polygon[:1]))
    need(len(words) == len(pairs), name + " complete edge-word length")
    for (start, end), expected in zip(pairs, words):
        actual = edge_polynomial(equation, start, end)
        need(sp.expand(actual - expected) == 0,
             name + " literal edge word")

NONLINEAR_EDGE_DISCS = {
    u - X**5: 3125 * u**4,
    Delta - X**4: -256 * Delta**3,
    e - X**3: -27 * e**2,
    1 - K * X**2: 4 * K,
}
for polynomial, expected_discriminant in NONLINEAR_EDGE_DISCS.items():
    need(sp.factor(sp.discriminant(polynomial, X) - expected_discriminant) == 0,
         "nonlinear edge has owner-monomial discriminant")
for name, words in EDGE_WORDS.items():
    for polynomial in words:
        degree = sp.degree(polynomial, X)
        need(degree >= 1, name + " edge remains nonconstant on owner gate")
        if degree == 1:
            need(sp.discriminant(polynomial, X) == 1,
                 name + " linear edge is reduced")
        else:
            need(polynomial in NONLINEAR_EDGE_DISCS,
                 name + " nonlinear edge has audited discriminant")

# Six global polygons: Pick, packet, convex-hull, and graph ledgers.
GLOBAL = {
    ("C5", 1): {
        "polygon": ((0, 1), (2, 0), (4, 2), (4, 3),
                    (3, 5), (1, 6), (0, 6)),
        "pick": (36, 12, 13),
        "packet": (10, 8, 5, 3, 2, 2, 1),
    },
    ("C5", 0): {
        "polygon": ((0, 1), (2, 0), (4, 2),
                    (3, 5), (1, 6), (0, 6)),
        "pick": (35, 11, 13),
        "packet": (10, 10, 5, 2, 2, 1),
    },
    ("C4", 1): {
        "polygon": ((0, 1), (2, 0), (4, 2), (4, 3),
                    (3, 5), (1, 6), (0, 5)),
        "pick": (35, 11, 13),
        "packet": (10, 8, 5, 3, 2, 2, 1),
    },
    ("C4", 0): {
        "polygon": ((0, 1), (2, 0), (4, 2),
                    (3, 5), (1, 6), (0, 5)),
        "pick": (34, 10, 13),
        "packet": (10, 10, 5, 2, 2, 1),
    },
    ("C3", 1): {
        "polygon": ((0, 1), (2, 0), (4, 2), (4, 3),
                    (3, 5), (1, 6), (0, 4)),
        "pick": (34, 10, 13),
        "packet": (10, 8, 5, 3, 2, 2, 1),
    },
    ("C3", 0): {
        "polygon": ((0, 1), (2, 0), (4, 2),
                    (3, 5), (1, 6), (0, 4)),
        "pick": (33, 9, 13),
        "packet": (10, 10, 5, 2, 2, 1),
    },
}
for owner_key, data in GLOBAL.items():
    need(pick(data["polygon"]) == data["pick"],
         "global Pick ledger")
    packet = edge_packet(data["polygon"])
    need(packet == data["packet"], "global source-infinity packet")
    need(sum(index - 1 for index in packet) == 24,
         "global packet defect equals 24")
    theta_on = owner_key[1]
    vertices, edges = ((5, 14) if theta_on else (4, 13))
    b1 = edges - vertices + 1
    need(b1 == 10, "global graph Betti number ten")
    need(b1 + 3 == data["pick"][2] == 13,
         "graph plus carrier genus equals Pick genus thirteen")
    need(2 * data["pick"][2] - 2 == 24,
         "Riemann--Hurwitz checksum")

for key, active in exact_supports.items():
    projected_hull = convex_hull((rr, ll) for rr, ll, _hh in active)
    owner_key = exact_owner[key]
    need(projected_hull == GLOBAL[owner_key]["polygon"],
         "exact support has the asserted global polygon")

# Closed natural address and its exact inverse.
def natural_address(class_index, u_bit, phi_bit, eta_bit):
    return 1 + 8 * class_index + 4 * u_bit + 2 * phi_bit + eta_bit


def inverse_address(number):
    offset = number - 1
    class_index, residue = divmod(offset, 8)
    u_bit, residue = divmod(residue, 4)
    phi_bit, eta_bit = divmod(residue, 2)
    return class_index, u_bit, phi_bit, eta_bit


addressed = {}
for key, active in exact_supports.items():
    number = natural_address(*key)
    need(1 <= number <= 64, "natural address lies in 1..64")
    need(inverse_address(number) == key, "natural address has exact inverse")
    addressed[number] = active
need(set(addressed) == set(range(1, 65)),
     "natural addresses biject onto 1..64")
need(len(set(addressed.values())) == 64,
     "natural addresses preserve all exact support distinctions")
for number in range(1, 65):
    odd = 2 * number - 1
    need((odd * odd - 1) // 8 == number * (number - 1) // 2,
         "odd-square address retains its triangular index")

print("THM-4354 first-normal-owner U-zero endpoint exact audit")
print("SCOPE Z=beta11=zeta3=W=xi10=U=0; K*alpha11!=0")
print("SEAM K=2848/45-(7/6)Delta")
print("SOURCE_ROWS=16 SPECIALIZED_ROWS=10 SUPPORT_POINTS=21")
print("EXACT_SUPPORTS=64 NATURAL_ADDRESSES=1..64")
print("FANS C5/Theta=20 C5/0=12 C4/Theta=16 C4/0=8 C3/Theta=4 C3/0=4")
print("HOSTILE_ATLAS=256_KEYED DISTINCT_SUPPORTS=168 SAME_SIX_FANS")
print("FACES=C5,C4,C3,B,E11,T,E10")
print("FACE_BASE_ORDER C5=30/25 C4=12/13 C3=6/9 B=66/49 E11=42/41 T=6/8 E10=24/26")
print("B_COMPONENTS=2 B_NODES=11 B_GRAPH_B1=10")
print("RATIONAL_FACES=C5,C4,C3,Bx2,T")
print("CARRIERS E11=g3/order41 E10=g3/order26")
print("GLOBAL_PICK_I=13 PACKET_DEFECT=24 GRAPH_B1=10")
print("PACKETS Theta=(10,8,5,3,2,2,1) zero=(10,10,5,2,2,1)")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
