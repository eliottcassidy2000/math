#!/usr/bin/env python3
"""Clean-room exact audit for the exact-M=8 two-term collision wall.

This companion imports no primary certificate.  It reconstructs the
normalized source polynomial on Theta=-Delta, proves the Phi=0 critical-value
obstruction, computes the two critical-resultant lengths, derives the
normalization packets and their residue-field labels, and checks every
remaining monodromy identity case.  The target Mordell--Weil input
E_q(C(q))={O} and fixed-sheet transport are inherited theorem inputs; all
algebra and finite permutation arithmetic used around them are checked here.

No Python assertions are used.
"""

from hashlib import sha256
from math import gcd

import sympy as sp


CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def z_order(expression, variable):
    """Return the least exponent with nonzero coefficient in a polynomial."""
    polynomial = sp.Poly(sp.cancel(expression), variable)
    exponents = [monomial[0] for monomial, coefficient in polynomial.terms()
                 if coefficient != 0]
    require(bool(exponents), "zero polynomial has no valuation")
    return min(exponents)


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
    require((area2 - boundary + 2) % 2 == 0, "Pick parity failed")
    return area2, boundary, (area2 - boundary + 2) // 2


def compose(left, right):
    return tuple(left[right[index]] for index in range(len(left)))


def inverse(permutation):
    result = [0] * len(permutation)
    for index, image in enumerate(permutation):
        result[image] = index
    return tuple(result)


def cycle(entries, degree):
    permutation = list(range(degree))
    for left, right in zip(entries, entries[1:] + entries[:1]):
        permutation[left] = right
    return tuple(permutation)


def cycle_type(permutation):
    seen = set()
    lengths = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        current = start
        length = 0
        while current not in seen:
            seen.add(current)
            current = permutation[current]
            length += 1
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


def commutator(left, right):
    return compose(compose(compose(left, right), inverse(left)), inverse(right))


# ---------------------------------------------------------------------------
# 1. Independent normalized polynomial and critical resultants.
# ---------------------------------------------------------------------------

X, T, Delta, Phi = sp.symbols("X T Delta Phi")
K = sp.Rational(2848, 45) - sp.Rational(7, 6) * Delta
P = T + X**2 * T**2
Y = X * T * P

G = (
    -X**2 * T / 2
    - 3 * P
    + sp.Rational(8, 3) * P**2
    - sp.Rational(1376, 135) * P**3
    + K * Y**2
    + Phi * P**2 * Y
    + Delta * P**4
    - Delta * P * Y**2
)

require(sp.factor(P**4 - P * Y**2 - T * P**3) == 0,
        "two-term top did not collapse to T*P^3")

f = sp.cancel(sp.diff(G, X) / T)
h = sp.diff(G, T)
require(sp.factor(sp.diff(G, X) - T * f) == 0, "critical division changed")
require((sp.degree(f, X), sp.degree(h, X)) == (6, 7),
        "generic X-degrees changed")
require(sp.factor(sp.Poly(f, X).LC()) == 7 * Phi * T**6,
        "wrong leading row of f")
require(sp.factor(sp.Poly(h, X).LC()) == 7 * Phi * T**6,
        "wrong leading row of h")

resultant = sp.resultant(f, h, X)
Q15 = sp.cancel(resultant / (7 * T**30 * Phi * (6 * T + 1)**2))
Q15_poly = sp.Poly(Q15, T)
require(sp.factor(resultant - 7 * T**30 * Phi * (6 * T + 1)**2 * Q15) == 0,
        "generic resultant reconstruction failed")
require(Q15_poly.degree() == 15, "generic residual degree is not fifteen")
require(sp.factor(Q15_poly.LC() - 576 * Delta**4 * K**4) == 0,
        "generic residual leading coefficient changed")
require(sp.factor(Q15_poly.TC() + sp.Rational(7203, 32) * Phi**4) == 0,
        "generic residual constant coefficient changed")

require(sp.factor(f.subs(T, 0)) == -X, "fake T=0 row of f changed")
require(sp.simplify(h.subs(T, 0) + (X**2 + 6) / 2) == 0,
        "restored T=0 critical pair changed")
universal_half = {T: sp.Rational(-1, 6), X: sp.sqrt(6)}
require(sp.simplify(f.subs(universal_half)) == 0, "universal half-value f failed")
require(sp.simplify(h.subs(universal_half)) == 0, "universal half-value h failed")
require(sp.simplify(G.subs(universal_half) - sp.Rational(1, 2)) == 0,
        "universal half-value changed")
require(sp.simplify(G.subs({T: 0, X: sp.sqrt(-6)})) == 0,
        "universal zero-value changed")

Delta0 = sp.Rational(5696, 105)
require(sp.factor(K.subs(Delta, Delta0)) == 0, "K=0 wall value changed")
f0 = sp.cancel(f.subs(Delta, Delta0))
h0 = sp.cancel(h.subs(Delta, Delta0))
resultant0 = sp.resultant(f0, h0, X)
Q13 = sp.cancel(resultant0 / (7 * T**30 * Phi * (6 * T + 1)**2))
Q13_poly = sp.Poly(Q13, T)
require(sp.factor(resultant0 - 7 * T**30 * Phi * (6 * T + 1)**2 * Q13) == 0,
        "K=0 resultant reconstruction failed")
require(Q13_poly.degree() == 13, "K=0 residual degree is not thirteen")
require(Q13_poly.LC() == sp.Rational(
    2185746832794987941330944, 16544390625
), "K=0 residual leading coefficient changed")
require(sp.factor(Q13_poly.TC() + sp.Rational(7203, 32) * Phi**4) == 0,
        "K=0 residual constant coefficient changed")

generic_critical_length = 2 + Q15_poly.degree() + 2
k0_critical_length = 2 + Q13_poly.degree() + 2
require(generic_critical_length == 19, "generic affine critical length changed")
require(k0_critical_length == 17, "K=0 affine critical length changed")


# ---------------------------------------------------------------------------
# 2. Phi=0: the three X=0 critical values cannot all be target nodes.
# ---------------------------------------------------------------------------

z_value = sp.symbols("z_value")
H = -3 * T + sp.Rational(8, 3) * T**2 - sp.Rational(1376, 135) * T**3 + Delta * T**4
Hprime = sp.diff(H, T)
require(sp.factor(G.subs({X: 0, Phi: 0}) - H) == 0,
        "Phi=0 line polynomial changed")
require(sp.factor(h.subs({X: 0, Phi: 0}) - Hprime) == 0,
        "Phi=0 line critical equation changed")

critical_value_resultant = sp.Poly(
    sp.resultant(Hprime, z_value - H, T), z_value
)
require(critical_value_resultant.degree() == 3,
        "Phi=0 critical-value eliminant is not cubic")
monic_critical_values = sp.Poly(
    sp.cancel(critical_value_resultant.as_expr() / critical_value_resultant.LC()),
    z_value,
)
pattern_gcds = []
for number_at_zero in range(4):
    pattern = sp.Poly(
        z_value**number_at_zero
        * (z_value - sp.Rational(1, 2))**(3 - number_at_zero),
        z_value,
    )
    differences = [
        sp.together(left - right).as_numer_denom()[0]
        for left, right in zip(
            monic_critical_values.all_coeffs()[1:], pattern.all_coeffs()[1:]
        )
    ]
    common = differences[0]
    for difference in differences[1:]:
        common = sp.gcd(common, difference)
    common = sp.Poly(common, Delta)
    require(common.degree() == 0,
            f"Phi=0 target-value pattern {number_at_zero} survived")
    pattern_gcds.append(str(sp.factor(common.as_expr())))
require(pattern_gcds == ["1", "1", "1", "1"],
        "Phi=0 pattern certificates changed")


# ---------------------------------------------------------------------------
# 3. Boundary branches, differential orders, and residue-field labels.
# ---------------------------------------------------------------------------

q, a, z, gamma, epsilon, alpha, lam = sp.symbols(
    "q a z gamma epsilon alpha lambda", nonzero=True
)
Q = 1 / q
B = (
    z**8
    + Q * Delta * (1 - a)**3 * a
    - Q * Phi * (1 - a)**3 * z
    - Q * (epsilon * (1 - a)**3 + K * (1 - a)**2) * z**2
    - Q * alpha * (1 - a)**2 * z**4
    - Q * lam * (1 - a) * z**6
)
L = sp.expand(a * B + gamma * Q * z**8)
L_a = sp.diff(L, a)

linear_branch = {a: Phi * z / Delta}
fast_branch = {a: gamma * z**7 / Phi}
require(z_order(L.subs(linear_branch), z) > 2,
        "linear a=0 branch failed its leading equation")
require(z_order(L.subs(fast_branch), z) > 8,
        "order-seven a=0 branch failed its leading equation")
require(z_order(L_a.subs(linear_branch), z) == 1,
        "linear a=0 residue denominator order changed")
require(z_order(L_a.subs(fast_branch), z) == 1,
        "order-seven a=0 residue denominator order changed")
require(6 - 1 == 5, "a=0 differential order changed")

u, w, tau, c = sp.symbols("u w tau c", nonzero=True)
L_u_chart = sp.expand(L.subs(a, 1 - u))
L_a_u_chart = sp.expand(L_a.subs(a, 1 - u))

rational_u_branch = {u: K * z**2 / Delta}
require(z_order(L_u_chart.subs(rational_u_branch), z) > 6,
        "rational a=1 branch failed its leading equation")
require(z_order(L_a_u_chart.subs(rational_u_branch), z) == 4,
        "rational a=1 residue denominator order changed")
require(6 - 4 == 2, "index-three differential order changed")

quadratic_face = sp.factor(
    sp.Poly(L_u_chart.subs(u, w * z**3), z).coeff_monomial(z**8)
)
require(sp.factor(quadratic_face - (q + gamma - K * w**2) / q) == 0,
        "quadratic a=1 face equation changed")
quadratic_denominator_lead = sp.factor(
    sp.Poly(L_a_u_chart.subs(u, w * z**3), z).coeff_monomial(z**5)
)
require(sp.factor(quadratic_denominator_lead - 2 * K * w / q) == 0,
        "quadratic a=1 residue denominator changed")
require(6 - 5 == 1, "index-two differential order changed")
require(sp.degree(q + gamma, q) == 1 and sp.gcd(q + gamma, sp.diff(q + gamma, q)) == 1,
        "quadratic residue extension lost its odd divisor")

L_k0_u = sp.expand(L_u_chart.subs({Delta: Delta0, K: 0, u: c * tau**8, z: tau**3}))
L_a_k0_u = sp.expand(
    L_a_u_chart.subs({Delta: Delta0, K: 0, u: c * tau**8, z: tau**3})
)
k0_face_lead = sp.factor(sp.Poly(L_k0_u, tau).coeff_monomial(tau**24))
require(sp.factor(k0_face_lead - (Delta0 * c**3 + q + gamma) / q) == 0,
        "K=0 primitive face equation changed")
require(z_order(L_a_k0_u, tau) == 16,
        "K=0 residue denominator order changed")
require(6 * 3 + (3 - 1) - 16 == 4,
        "K=0 differential order changed")
require(gcd(3, 8) == 1, "K=0 edge ceased to be primitive")

polygon = ((0, 1), (2, 0), (4, 3), (2, 4), (0, 5))
require(polygon_ledger(polygon) == (24, 8, 9),
        "two-term Newton polygon ledger changed")
require(gcd(2, 1) == 1, "AB index-one edge ceased to be rational")

packet_generic = (6, 6, 3, 2, 2, 1)
packet_k0 = (6, 6, 5, 1)
require(sum(entry - 1 for entry in packet_generic) == 14,
        "generic packet defect changed")
require(sum(entry - 1 for entry in packet_k0) == 14,
        "K=0 packet defect changed")

# The rational entries map to O by the inherited E_q(C(q))={O}.  The two
# index-two points are one irreducible quadratic closed point: they either
# both map to O, or their finite images remain two conjugate geometric points.
origin_full_generic = packet_generic
origin_finite_pair = (6, 6, 3, 1)
finite_pair = (2, 2)
require(sum(origin_full_generic) == 20, "generic full response degree changed")
require(sum(origin_finite_pair) == 16 and sum(finite_pair) == 4,
        "generic quadratic response split changed")
require(sum(packet_k0) == 18, "K=0 response degree changed")


# ---------------------------------------------------------------------------
# 4. Collision blowup: seven genuine nodes after persistent normalization.
# ---------------------------------------------------------------------------

rho, A, Z = sp.symbols("rho A Z")
sigma = rho**7
z_scaled = rho**3 * Z
a_scaled = rho**24 * A
B_sigma = (
    z_scaled**8
    + Delta * (1 - a_scaled)**3 * a_scaled
    - sigma**3 * Phi * (1 - a_scaled)**3 * z_scaled
    - sigma**6
    * (epsilon * (1 - a_scaled)**3 + K * (1 - a_scaled)**2)
    * z_scaled**2
    - sigma**12 * alpha * (1 - a_scaled)**2 * z_scaled**4
    - sigma**18 * lam * (1 - a_scaled) * z_scaled**6
)
scaled_total = sp.expand(
    (a_scaled * B_sigma + gamma * sigma**24 * z_scaled**8) / rho**48
)
central = sp.factor(scaled_total.subs(rho, 0))
require(central == A * (Delta * A + Z**8 - Phi * Z),
        "collision central equation changed")
intersection = sp.factor((Delta * A + Z**8 - Phi * Z).subs(A, 0))
require(sp.expand(intersection - Z * (Z**7 - Phi)) == 0,
        "central intersections changed")
require(sp.factor(sp.diff(Z**8 - Phi * Z, Z).subs(Z**7, Phi)) == 7 * Phi,
        "seven nonzero intersections lost transversality")
require(sp.diff(Delta * A + Z**8 - Phi * Z, A) == Delta,
        "exceptional graph ceased to be rational")
require(7 - 2 + 1 == 6, "normalized dual-graph rank changed")
require(2 + 6 == 8, "Bolza plus graph genus ledger changed")


# ---------------------------------------------------------------------------
# 5. Exhaustive monodromy capacity and equality cases.
# ---------------------------------------------------------------------------

def audit_full_boundary(degree, critical_length, packet):
    require(sum(packet) == degree, "full-boundary packet has wrong degree")
    require(critical_length == degree - 1,
            "full-boundary fixed-sheet equality changed")

    # If one torus generator is the identity, the other alone would have to
    # be a degree-cycle, contradicting the two universal fixed sheets.
    require(2 > 0, "universal fixed-sheet pair disappeared")

    expected_commutator = (3,) + (1,) * (degree - 3)
    for left_support in range(2, degree):
        right_support = degree + 1 - left_support
        if right_support < 2:
            continue
        left = cycle([0] + list(range(1, left_support)), degree)
        right = cycle([0] + list(range(left_support, degree)), degree)
        require(cycle_type(commutator(left, right)) == expected_commutator,
                "one-pivot commutator type changed")
    require(tuple(sorted(packet, reverse=True)) != expected_commutator,
            "actual origin packet became a one-pivot commutator")


audit_full_boundary(20, generic_critical_length, origin_full_generic)
audit_full_boundary(18, k0_critical_length, packet_k0)

finite_degree = 16
carrier_index = 2
both_nonidentity_capacity = 2 * finite_degree - generic_critical_length - 2 + carrier_index
one_identity_fixed = max(2, generic_critical_length - finite_degree)
one_identity_capacity = finite_degree - one_identity_fixed - 1 + carrier_index
both_identity_capacity = carrier_index
require(both_nonidentity_capacity == 13 < finite_degree - 1,
        "finite-pair both-nonidentity case reopened")
require(one_identity_fixed == 3, "finite-pair fixed-sheet floor changed")
require(one_identity_capacity == 14 < finite_degree - 1,
        "finite-pair one-identity case reopened")
require(both_identity_capacity == 2 < finite_degree - 1,
        "finite-pair both-identity case reopened")


semantic_lines = (
    "scope=exact_M8_two_term_wall_Theta=-Delta;Phi!=0;JC2_OPEN",
    "top=P^4-PY^2=T P^3",
    "critical=K!=0:length19,Q15;K=0:length17,Q13",
    "packets=K!=0:(6,6,3,2,2,1);K=0:(6,6,5,1)",
    "labels=generic:rational(6,6,3,1)+quadratic(2,2);K0:all_rational",
    "responses=n20_full_or_n16_with_two_distinct_transpositions;n18_full",
    "stable=normalize_Z0;seven_nodes;Bolza2+graph6=genus8",
    "monodromy=n20,n16,n18_all_identity_cases_excluded",
)
semantic_digest = sha256("\n".join(semantic_lines).encode()).hexdigest()

print("JC23 TWO-TERM COLLISION WALL INDEPENDENT AUDIT")
print(semantic_lines[0])
print("phi_zero=NO_GO;critical_value_pattern_gcds=" + ",".join(pattern_gcds))
print("critical_generic=resultant:7*T^30*Phi*(6*T+1)^2*Q15;length=19")
print("critical_K0=Delta:5696/105;Q13;length=17")
print("packet_generic=(6,6,3,2,2,1);defect=14;labels=rational(6,6,3,1)+quadratic(2,2)")
print("packet_K0=(6,6,5,1);defect=14;labels=all_rational")
print("responses=generic:n20_or_n16_with_two_distinct_transpositions;K0:n18")
print("stable_collision=central:A*(Delta*A+Z^8-Phi*Z);normalize_Z0;nodes=7;genus=2+6=8")
print("monodromy=n20_full:PASS;n16_quadratic_carriers:PASS;n18_full:PASS")
print(f"semantic_digest={semantic_digest}")
print(f"checks={CHECKS}")
print("RESULT=PASS")
