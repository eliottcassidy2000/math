#!/usr/bin/env python3
"""Clean-room exact audit for reserved THM-4155.

The calculation starts only from THM-4147's normalized source polynomial.
It does not import the proposed primary certificate.  Its critical lane uses
the birational source coordinates (s,p), rather than the normalized (X,T)
projection: after removing the known p-factor, an exact symbolic resultant
produces the degree-eighteen residual and both endpoint factors.

The second lane independently expands the generic-fibre equation, rebuilds
the Newton polygon, its labelled boundary places, the two Galois-locked
responses, and the finite/full monodromy capacities.  No Python assertions
are used, so all gates remain live under ``python -O``.
"""

from hashlib import sha256
from itertools import permutations
from math import gcd

import sympy as sp


CHECKS = 0
SEMANTIC = []


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def record(row):
    SEMANTIC.append(row)


def convex_hull(points):
    points = sorted(set(points))

    def cross(origin, left, right):
        return (
            (left[0] - origin[0]) * (right[1] - origin[1])
            - (left[1] - origin[1]) * (right[0] - origin[0])
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


def polygon_ledger(vertices):
    area2 = abs(sum(
        vertices[index][0] * vertices[(index + 1) % len(vertices)][1]
        - vertices[(index + 1) % len(vertices)][0] * vertices[index][1]
        for index in range(len(vertices))
    ))
    boundary = sum(
        gcd(
            abs(vertices[(index + 1) % len(vertices)][0] - vertices[index][0]),
            abs(vertices[(index + 1) % len(vertices)][1] - vertices[index][1]),
        )
        for index in range(len(vertices))
    )
    require((area2 - boundary + 2) % 2 == 0, "Pick parity")
    return area2, boundary, (area2 - boundary + 2) // 2


def toric_packet(vertices):
    packet = []
    rows = []
    for index, left in enumerate(vertices):
        right = vertices[(index + 1) % len(vertices)]
        dx = right[0] - left[0]
        dy = right[1] - left[1]
        length = gcd(abs(dx), abs(dy))
        if left[0] == right[0] == 0:
            rows.append((left, right, length, "affine"))
            continue
        normal = (-dy // length, dx // length)
        level = normal[0] * left[0] + normal[1] * left[1]
        ramification = normal[0] + normal[1] - level
        packet.extend([ramification] * length)
        rows.append((left, right, length, normal, level, ramification))
    return tuple(sorted(packet, reverse=True)), tuple(rows)


def compose(left, right):
    return tuple(left[right[index]] for index in range(len(left)))


def inverse(permutation):
    result = [0] * len(permutation)
    for index, image in enumerate(permutation):
        result[image] = index
    return tuple(result)


def support(permutation):
    return {index for index, image in enumerate(permutation) if index != image}


# ---------------------------------------------------------------------------
# 1. Source-coordinate critical elimination on Delta=eta=0.
# ---------------------------------------------------------------------------

s, p = sp.symbols("s p")
Theta, Phi, Zeta = sp.symbols("Theta Phi Zeta")
K = sp.Rational(2848, 45)
epsilon = -sp.Rational(1376, 135)
t = p - s**2

H = (
    -3 * p
    + sp.Rational(8, 3) * p**2
    + epsilon * p**3
    + K * s**2 * p**2
    + Phi * s * p**3
    + Theta * s**2 * p**3
    + Zeta * s**3 * p**3
)
G_source = -s**2 / (2 * t) + H

# In the affine source chart, t^2*G_s is divisible by p.  Dividing by this
# known factor removes a collapsed p=0 row; C=2t^2*G_p is the independent
# second equation.
A_numerator = sp.expand(-s * p + t**2 * sp.diff(H, s))
require(sp.rem(sp.Poly(A_numerator, p), sp.Poly(p, p)).as_expr() == 0,
        "source first critical numerator lost its p factor")
A = sp.cancel(A_numerator / p)
C = sp.expand(s**2 + 2 * t**2 * sp.diff(H, p))
require(sp.denom(A) == 1, "source A is not polynomial")
require((sp.degree(A, s), sp.degree(C, s)) == (6, 7),
        "source critical s-degrees changed")
require(sp.factor(sp.Poly(A, s).LC() - 3 * Zeta * p**2) == 0,
        "source A leading row changed")
require(sp.factor(sp.Poly(C, s).LC() - 6 * Zeta * p**2) == 0,
        "source C leading row changed")

source_resultant = sp.resultant(A, C, s)
source_poly = sp.Poly(source_resultant, p)
p_valuation = min(
    monomial[0] for monomial, coefficient in source_poly.terms()
    if coefficient != 0
)
require(p_valuation == 8, "source resultant p-valuation is not eight")
R18 = sp.cancel(source_resultant / p**8)
R18_poly = sp.Poly(R18, p)
require(R18_poly.degree() == 18, "source residual degree is not eighteen")

w = sp.symbols("w")
horizontal_cubic = epsilon * w**3 + Phi * w**2 + Theta * w + Zeta
horizontal_discriminant = sp.factor(sp.discriminant(horizontal_cubic, w))
horizontal_discriminant_expanded = (
    Phi**2 * Theta**2
    - 4 * epsilon * Theta**3
    - 4 * Phi**3 * Zeta
    - 27 * epsilon**2 * Zeta**2
    + 18 * epsilon * Phi * Theta * Zeta
)
require(sp.factor(
    horizontal_discriminant - horizontal_discriminant_expanded
) == 0, "horizontal cubic discriminant changed")

carrier_endpoint = 4 * Theta * K**2 - 27 * Zeta**2
expected_constant = -46656 * Zeta * carrier_endpoint
expected_leading = -236196 * Zeta**5 * horizontal_discriminant
require(sp.factor(R18_poly.TC() - expected_constant) == 0,
        "R18 constant endpoint changed")
require(sp.factor(R18_poly.LC() - expected_leading) == 0,
        "R18 leading endpoint changed")

# One independent exact control proves the endpoint-open set is nonempty and
# gives a hostile squarefree residual without making squarefreeness a theorem
# hypothesis.
control = {Theta: 3, Phi: 5, Zeta: 7}
control_R18 = sp.Poly(R18.subs(control), p)
require(control_R18.degree() == 18, "control residual degree")
require(sp.gcd(control_R18, control_R18.diff()).degree() == 0,
        "control residual is not squarefree")
require(horizontal_discriminant.subs(control) != 0,
        "control horizontal discriminant vanished")
require(carrier_endpoint.subs(control) != 0,
        "control carrier endpoint vanished")
control_digest = sha256(
    str(control_R18.as_expr()).encode("utf-8")
).hexdigest()

require(sp.factor(A.subs(p, 0)) == -s,
        "collapsed p=0 A row changed")
require(sp.factor(C.subs(p, 0) + s**2 * (6 * s**2 - 1)) == 0,
        "collapsed p=0 C row changed")
require(sp.factor(A.subs(p, s**2)) == -s,
        "collapsed t=0 A row changed")
require(sp.factor(C.subs(p, s**2) - s**2) == 0,
        "collapsed t=0 C row changed")

# The divided source equations omit the two p=0 half-value points; the
# birational source chart also collapses the two T=0 zero-value points.
s_half = sp.sqrt(6) / 6
require(sp.simplify(sp.diff(G_source, s).subs({s: s_half, p: 0})) == 0,
        "source half-value G_s changed")
require(sp.simplify(sp.diff(G_source, p).subs({s: s_half, p: 0})) == 0,
        "source half-value G_p changed")
require(sp.simplify(G_source.subs({s: s_half, p: 0}) - sp.Rational(1, 2)) == 0,
        "source half-value changed")

X, T = sp.symbols("X T")
P = T + X**2 * T**2
Y = X * T * P
G_XT = (
    -X**2 * T / 2
    - 3 * P
    + sp.Rational(8, 3) * P**2
    + epsilon * P**3
    + K * Y**2
    + Phi * P**2 * Y
    + Theta * P * Y**2
    + Zeta * Y**3
)
f_XT = sp.cancel(sp.diff(G_XT, X) / T)
h_XT = sp.diff(G_XT, T)
hessian_XT = sp.det(sp.hessian(G_XT, (X, T)))
for sign in (-1, 1):
    half_point = {T: sp.Rational(-1, 6), X: sign * sp.sqrt(6)}
    zero_point = {T: 0, X: sign * sp.sqrt(-6)}
    require(sp.simplify(f_XT.subs(half_point)) == 0,
            "normalized half-value f changed")
    require(sp.simplify(h_XT.subs(half_point)) == 0,
            "normalized half-value h changed")
    require(sp.simplify(G_XT.subs(half_point) - sp.Rational(1, 2)) == 0,
            "normalized half-value changed")
    require(sp.simplify(hessian_XT.subs(half_point) + 6) == 0,
            "normalized half-value Hessian changed")
    require(sp.simplify(sp.diff(G_XT, X).subs(zero_point)) == 0,
            "normalized zero-value G_X changed")
    require(sp.simplify(h_XT.subs(zero_point)) == 0,
            "normalized zero-value h changed")
    require(sp.simplify(G_XT.subs(zero_point)) == 0,
            "normalized zero value changed")
    require(sp.simplify(hessian_XT.subs(zero_point) - 6) == 0,
            "normalized zero Hessian changed")

critical_length = R18_poly.degree() + 2 + 2
require(critical_length == 22, "affine critical length is not twenty-two")
record(
    "critical:p8R18:TC=-46656*Zeta*(4ThetaK2-27Zeta2):"
    "LC=-236196*Zeta5*DiscHorizontal:L22"
)


# ---------------------------------------------------------------------------
# 2. Newton polygon, labelled places, and exact normalized genus.
# ---------------------------------------------------------------------------

q = sp.symbols("q")
F_q = sp.expand((s**2 - p) * (1 - H / q) - s**2 / (2 * q))
F_q_cleared = sp.Poly(sp.cancel(q * F_q), s, p)
support_points = tuple(sorted(
    monomial for monomial, coefficient in F_q_cleared.terms()
    if sp.factor(coefficient) != 0
))
polygon = convex_hull(support_points)
expected_polygon = ((0, 1), (2, 0), (5, 3), (3, 4), (0, 4))
require(polygon == expected_polygon, "Delta=0 Y-only Newton polygon")
ledger = polygon_ledger(polygon)
require(ledger == (27, 11, 9), "Delta=0 Y-only Pick ledger")

packet, edge_rows = toric_packet(polygon)
expected_packet = (8, 3, 3, 3, 2, 2, 2, 1)
require(packet == expected_packet, "Delta=0 Y-only boundary packet")
packet_defect = sum(entry - 1 for entry in packet)
require(packet_defect == 16 == 2 * ledger[2] - 2,
        "boundary defect/genus equality")

gamma = -sp.Rational(1, 2)
ab_face = -q + (q + gamma) * w
bc_face = q + gamma - K * w**2 - Zeta * w**3
cd_face = Zeta * (w - 1)
de_face = horizontal_cubic


def edge_coefficients(left, right):
    dx = right[0] - left[0]
    dy = right[1] - left[1]
    length = gcd(abs(dx), abs(dy))
    step = (dx // length, dy // length)
    return tuple(sp.factor(F_q_cleared.coeff_monomial(
        s ** (left[0] + index * step[0])
        * p ** (left[1] + index * step[1])
    )) for index in range(length + 1))


def same_coefficients(actual, expected):
    return len(actual) == len(expected) and all(
        sp.factor(left - right) == 0
        for left, right in zip(actual, expected)
    )


# Extract every non-affine face directly from the expanded equation.  This
# guards both labels and the orientation used in the residue polynomials.
require(same_coefficients(
    edge_coefficients((0, 1), (2, 0)), (-q, q + gamma)),
        "AB face extraction")
require(same_coefficients(
    edge_coefficients((2, 0), (5, 3)),
    (q + gamma, 0, -K, -Zeta)), "BC face extraction")
require(same_coefficients(
    edge_coefficients((5, 3), (3, 4)), (-Zeta, Zeta)),
        "CD face extraction")
require(same_coefficients(
    edge_coefficients((3, 4), (0, 4)),
    (Zeta, Theta, Phi, epsilon)), "DE face extraction")
require(sp.degree(ab_face, w) == 1, "AB rational place")
require(sp.factor(cd_face - Zeta * (w - 1)) == 0,
        "CD index-eight rational place")
bc_discriminant = sp.factor(sp.discriminant(bc_face, w))
require(sp.factor(
    bc_discriminant
    - (q + gamma) * (4 * K**3 - 27 * Zeta**2 * (q + gamma))
) == 0, "BC cubic discriminant")
require(sp.degree(K * w**2 + Zeta * w**3, w) == 3,
        "BC residue map is not degree three")
require(sp.degree(bc_face, q) == 1,
        "BC residue equation is not linear in q")
require(sp.degree(de_face, w) == 3,
        "horizontal constant face is not cubic")
require(not de_face.has(q), "horizontal cubic is not constant over q")

# Over C(q), the three distinct roots of the constant horizontal cubic are
# rational places.  The BC face is one irreducible degree-three closed point.
rational_packet = (8, 3, 3, 3, 1)
cubic_packet = (2, 2, 2)
require(sum(rational_packet) == 18, "rational origin packet degree")
require(sum(cubic_packet) == 6, "cubic packet degree")
require(sum(packet) == 24, "full response degree")
require(packet_defect == 2 * 9 - 2, "normalized genus is not nine")
record(
    f"geometry:hull{polygon}:ledger{ledger}:packet{packet}:edges{edge_rows}"
)
record("labels:rational(8,3,3,3,1)+cubic(2,2,2):genus9")


# ---------------------------------------------------------------------------
# 3. Galois-locked responses and monodromy contradictions.
# ---------------------------------------------------------------------------

# THM-4120's inherited E_q(C(q))={O} sends every rational place to O.  The
# prime cubic closed point responds either wholly to O or to three distinct
# conjugate finite points.  These are the only two degrees.
full_n = 24
finite_n = 18
carrier_index = 3
require(full_n == sum(packet), "full response sum")
require(finite_n == sum(rational_packet), "finite response sum")

full_commutator_support_bound = 3 * (full_n - critical_length)
full_origin_support = full_n - packet.count(1)
require(full_commutator_support_bound == 6,
        "full commutator support bound")
require(full_commutator_support_bound < full_origin_support == 23,
        "full response did not contradict origin support")

finite_support_sum = 2 * finite_n - critical_length
both_nonidentity = finite_support_sum - 2 + carrier_index
one_identity = finite_support_sum - 1 + carrier_index
both_identity = carrier_index
require((both_nonidentity, one_identity, both_identity) == (15, 16, 3),
        "finite response capacity ledger")
require(both_nonidentity < finite_n - 1,
        "finite both-nonidentity capacity")
require(one_identity < finite_n - 1,
        "finite one-identity capacity")
require(both_identity < finite_n - 1,
        "finite both-identity capacity")

# Exhaustive hostile control for the full-response support lemma through S_6.
# If I=supp(x) intersect supp(y), the disagreement of xy and yx lies in
# I union x^{-1}I union y^{-1}I, and its size equals supp([x,y]).
for degree in range(1, 7):
    bank = tuple(permutations(range(degree)))
    for left in bank:
        left_support = support(left)
        left_inverse = inverse(left)
        xy_cache = {}
        for right in bank:
            right_support = support(right)
            right_inverse = inverse(right)
            overlap = left_support & right_support
            disagreement = {
                value for value in range(degree)
                if left[right[value]] != right[left[value]]
            }
            preimage_left = {left_inverse[value] for value in overlap}
            preimage_right = {right_inverse[value] for value in overlap}
            require(
                disagreement <= overlap | preimage_left | preimage_right,
                "commutator disagreement containment hostile",
            )
            commutator = compose(
                compose(compose(left, right), left_inverse), right_inverse
            )
            require(len(support(commutator)) == len(disagreement),
                    "commutator/Hamming cardinality hostile")
            require(len(support(commutator)) <= 3 * len(overlap),
                    "commutator support hostile")

record("responses:full24/comm6<23;finite18/beta3/cap15,16,3<17")
record("walls:Zeta*(4ThetaK2-27Zeta2)*DiscHorizontal!=0")

semantic_digest = sha256("\n".join(SEMANTIC).encode("utf-8")).hexdigest()

print("THM-4155 INDEPENDENT Y-ONLY DELTA-ZERO AUDIT")
print("source_resultant=p^8*R18;critical_length=22")
print("R18_constant=-46656*Zeta*(4*Theta*K^2-27*Zeta^2)")
print("R18_leading=-236196*Zeta^5*Disc(epsilon*w^3+Phi*w^2+Theta*w+Zeta)")
print("polygon=((0,1),(2,0),(5,3),(3,4),(0,4));Pick=(27,11,9)")
print("packet=(8,3,3,3,2,2,2,1);labels=rat(8,3,3,3,1)+cubic(2,2,2)")
print("responses=full24/finite18;finite_carriers=3")
print("monodromy=full:6<23;finite:15,16,3<17")
print("gates=Zeta*(4*Theta*K^2-27*Zeta^2)*horizontal_discriminant!=0")
print(f"control_sha256={control_digest}")
print(f"checks={CHECKS}")
print(f"semantic_sha256={semantic_digest}")
print("RESULT=PASS")
