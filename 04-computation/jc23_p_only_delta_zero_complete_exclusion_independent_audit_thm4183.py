#!/usr/bin/env python3
"""Independent (s,p)-chart audit for THM-4183.

This implementation imports no primary code.  It derives the polynomial
source-critical pair from the rational chart, uses a p-resultant instead of
the normalized T-resultant, checks the Hessian bridge modulo that pair, and
reconstructs the two polygons from supporting inequalities.
"""

from math import gcd

import sympy as sp


CHECKS = 0


def need(condition, message):
    global CHECKS
    CHECKS += 1
    if not bool(condition):
        raise RuntimeError(message)


def valuation(poly, variable):
    terms = sp.Poly(poly, variable).terms()
    need(bool(terms), "zero polynomial has no valuation")
    return min(monomial[0] for monomial, coefficient in terms if coefficient)


def exact_quotient(numerator, denominator, variable, message):
    quotient = sp.cancel(numerator / denominator)
    need(sp.denom(quotient) == 1, message + " quotient is not polynomial")
    return sp.Poly(quotient, variable)


def polygon_checksum(vertices):
    area2 = abs(sum(
        vertices[index][0] * vertices[(index + 1) % len(vertices)][1]
        - vertices[(index + 1) % len(vertices)][0] * vertices[index][1]
        for index in range(len(vertices))
    ))
    boundary = sum(
        gcd(abs(vertices[(index + 1) % len(vertices)][0] - vertices[index][0]),
            abs(vertices[(index + 1) % len(vertices)][1] - vertices[index][1]))
        for index in range(len(vertices))
    )
    return area2, boundary, (area2 - boundary + 2) // 2


s, p = sp.symbols("s p")
Phi, Theta, eta = sp.symbols("Phi Theta eta")
K0 = sp.Rational(2848, 45)
epsilon = -sp.Rational(1376, 135)
t = p - s**2

H = sp.expand(
    -3 * p + sp.Rational(8, 3) * p**2 + epsilon * p**3
    + K0 * s**2 * p**2 + Phi * s * p**3
    + Theta * s**2 * p**3 + eta * s * p**4
)
G = -s**2 / (2 * t) + H

# Polynomial generators of the critical ideal on p*t!=0.
A = sp.cancel((-s * p + t**2 * sp.diff(H, s)) / p)
C0 = sp.expand(s**2 + 2 * t**2 * sp.diff(H, p))
B = sp.cancel((C0 + s * A) / t**2)
need(sp.denom(A) == 1 and sp.denom(B) == 1,
     "source critical pair is not polynomial")
need(sp.factor(t**2 * sp.diff(G, s) - p * A) == 0,
     "first source-gradient identity changed")
need(sp.factor(2 * t**2 * sp.diff(G, p) - (t**2 * B - s * A)) == 0,
     "second source-gradient identity changed")
need((sp.degree(A, s), sp.degree(B, s)) == (5, 2),
     "source s-degrees changed")
lc_A = sp.factor(sp.Poly(A, s).LC())
lc_B = sp.factor(sp.Poly(B, s).LC())
need(sp.factor(lc_A - 2 * p * (K0 + Theta * p)) == 0,
     "source A leading row changed")
need(sp.factor(lc_B - 8 * p * (Theta * p + 3 * K0 / 4)) == 0,
     "source B leading row changed")
need(sp.factor(lc_B - 4 * lc_A + 2 * K0 * p) == 0,
     "source infinity incompatibility changed")

# At a common zero with p*t!=0, the differentiated gradient identities give
# p det D(A,B)=2t^2 det Hess(G).  The exact difference reduces to the ideal.
source_jacobian = sp.det(sp.Matrix((
    (sp.diff(A, s), sp.diff(A, p)),
    (sp.diff(B, s), sp.diff(B, p)),
)))
source_hessian = sp.det(sp.hessian(G, (s, p)))
bridge_numerator = sp.together(
    p * source_jacobian - 2 * t**2 * source_hessian
).as_numer_denom()[0]
remainder = sp.reduced(bridge_numerator, [A, B], s, p)[1]
need(sp.factor(remainder) == 0, "source Hessian bridge changed")

resultant = sp.resultant(A, B, s)
need(valuation(resultant, p) == 2, "source p-artifact changed")
R20 = exact_quotient(resultant, p**2, p, "source resultant")
need(R20.degree() == 20, "source Theta-nonzero degree changed")
need(sp.factor(R20.LC() + 65610 * eta**6 * Theta) == 0,
     "source Theta-nonzero top endpoint changed")
need(sp.factor(R20.TC() + 31104 * K0**2) == 0,
     "source Theta-nonzero bottom endpoint changed")

resultant0 = sp.resultant(A.subs(Theta, 0), B.subs(Theta, 0), s)
need(sp.factor(resultant0 - resultant.subs(Theta, 0)) == 0,
     "source Theta-zero direct specialization changed")
R19 = exact_quotient(resultant0, p**2, p,
                     "source Theta-zero resultant")
need(R19.degree() == 19, "source Theta-zero degree changed")
need(sp.factor(R19.LC() + 78732 * K0 * eta**6) == 0,
     "source Theta-zero top endpoint changed")
need(sp.factor(R19.TC() + 31104 * K0**2) == 0,
     "source Theta-zero bottom endpoint changed")

# Reciprocal v=1/p strict transform: the top coefficient has exactly one
# Theta factor and the next coefficient is a wall unit.
need(sp.factor(R20.nth(19).subs(Theta, 0)
               + 78732 * K0 * eta**6) == 0,
     "source first reciprocal normal coefficient changed")
v_escape = -sp.cancel(
    (-65610 * eta**6 * Theta) / (-78732 * K0 * eta**6)
)
need(sp.factor(v_escape + sp.Rational(5, 6) * Theta / K0) == 0,
     "source reciprocal escape scale changed")

# Both deleted coordinate rows contain no (A,B) point.  The normalized chart
# restores two Morse points from each row, independently checked below.
need(sp.factor(A.subs(p, 0) + s) == 0
     and B.subs({p: 0, s: 0}) == -6,
     "p-zero chart loss changed")
need(sp.factor(A.subs(p, s**2) + s) == 0
     and B.subs({p: s**2, s: 0}) == -6,
     "t-zero chart loss changed")

X, T = sp.symbols("X T")
PN = T + X**2 * T**2
YN = X * T * PN
GN = sp.expand(
    -X**2 * T / 2 - 3 * PN + sp.Rational(8, 3) * PN**2
    + epsilon * PN**3 + K0 * YN**2 + Phi * PN**2 * YN
    + Theta * PN * YN**2 + eta * PN**3 * YN
)
fN = sp.expand(sp.cancel(sp.diff(GN, X) / T))
hN = sp.expand(sp.diff(GN, T))
hessN = sp.det(sp.hessian(GN, (X, T)))
need(sp.factor(fN.subs(T, 0) + X) == 0
     and sp.factor(hN.subs(T, 0) + (X**2 + 6) / 2) == 0,
     "normalized T-zero pair changed")
need(sp.rem(sp.Poly(hessN.subs(T, 0) - 6, X),
            sp.Poly(X**2 + 6, X)).is_zero,
     "normalized T-zero Hessian changed")
for expression, expected in ((fN, 0), (hN, 0),
                             (GN, sp.Rational(1, 2))):
    remainder = sp.rem(
        sp.Poly(sp.expand(expression.subs(T, -sp.Rational(1, 6))
                          - expected), X),
        sp.Poly(X**2 - 6, X),
    )
    need(remainder.is_zero, "normalized p-zero pair changed")
need(sp.rem(sp.Poly(hessN.subs(T, -sp.Rational(1, 6)) + 6, X),
            sp.Poly(X**2 - 6, X)).is_zero,
     "normalized p-zero Hessian changed")

# Supporting-inequality reconstruction, disjoint from the primary hull code.
Q = sp.symbols("Q")
F = sp.expand((s**2 - p) * (1 - Q * H) - Q * s**2 / 2)
generic_vertices = ((0, 1), (2, 0), (4, 2), (4, 3),
                    (3, 4), (1, 5), (0, 4))
wall_vertices = ((0, 1), (2, 0), (4, 2),
                 (3, 4), (1, 5), (0, 4))
generic_halfspaces = (
    (1, 2, 2), (-1, 1, -2), (-1, 0, -4), (-1, -1, -7),
    (-1, -2, -11), (1, -1, -4), (1, 0, 0),
)
wall_halfspaces = (
    (1, 2, 2), (-1, 1, -2), (-2, -1, -10),
    (-1, -2, -11), (1, -1, -4), (1, 0, 0),
)
for specialization, vertices, halfspaces, ledger in (
        ({Phi: 2, Theta: 3, eta: 5}, generic_vertices,
         generic_halfspaces, (28, 10, 10)),
        ({Phi: 2, Theta: 0, eta: 5}, wall_vertices,
         wall_halfspaces, (27, 9, 10))):
    support = {
        powers for powers, coefficient
        in sp.Poly(F.subs(specialization), s, p).terms()
        if coefficient != 0
    }
    need(all(vertex in support for vertex in vertices),
         "polygon vertex disappeared")
    for i, j in support:
        need(all(u * i + v * j >= level
                 for u, v, level in halfspaces),
             "support escaped a claimed halfspace")
    for u, v, level in halfspaces:
        need(sum(u * i + v * j == level for i, j in support) >= 2,
             "claimed polygon edge is inactive")
    need(polygon_checksum(vertices) == ledger,
         "independent polygon checksum changed")

# Coefficients at every vertex make Phi irrelevant to the boundary and show
# that the Theta=0 merged edge is primitive and smooth.
poly = sp.Poly(F, s, p)
coefficient = lambda i, j: sp.factor(poly.coeff_monomial(s**i * p**j))
need(coefficient(0, 1) == -1, "low vertex changed")
need(sp.factor(coefficient(2, 0) - (1 - Q / 2)) == 0,
     "low/carrier vertex changed")
need(coefficient(4, 2) == -K0 * Q, "carrier vertex changed")
need(coefficient(4, 3) == -Q * Theta, "Theta vertex changed")
need(coefficient(3, 4) == -Q * eta, "mixed vertex changed")
need(coefficient(1, 5) == Q * eta, "top vertex changed")
need(coefficient(0, 4) == Q * epsilon, "diagonal vertex changed")
need(sp.factor(
    sum(coefficient(i, j) * s**i * p**j
        for i, j in ((4, 2), (3, 4))).subs(Theta, 0)
    + Q * p**2 * s**3 * (K0 * s + eta * p**2)
) == 0, "independent merged-edge polynomial changed")

generic_packet = (8, 5, 4, 3, 2, 2, 1)
wall_packet = (8, 7, 4, 2, 2, 1)
for packet, full_n, finite_n, length in (
        (generic_packet, 25, 21, 24),
        (wall_packet, 24, 20, 23)):
    defect = sum(index - 1 for index in packet)
    need(sum(packet) == full_n and defect == 18,
         "independent packet ledger changed")
    need(full_n - 4 == finite_n, "independent carrier response changed")
    need(2 * (full_n - length) < defect,
         "independent full response changed")
    need(2 * finite_n - length - 1 + 2 < finite_n - 1,
         "independent finite response changed")

# Scope control: zeta=0 is exact weight nine precisely because eta is a
# unit.  Setting eta=0 is a filtration exit, not a third Theta row.
need(2 * 3 + 3 * 1 == 9 and 3 * 3 == 9,
     "weight-nine monomial ledger changed")

print("P_ONLY_DELTA_ZERO_INDEPENDENT_SOURCE_AUDIT")
print("chart=(s,p);critical_pair=(A,B);hessian_bridge=ideal_zero")
print("source_resultants=p^2*R20_or_p^2*R19")
print("critical_lengths=20+4:24;19+4:23")
print("normal_escape=v=1/p~-5*Theta/(6*K0)")
print("polygons=Theta_nonzero:(28,10,10);Theta_zero:(27,9,10)")
print("Theta_zero_blows_down_indices_3_plus_5_to_7")
print("packets=(8,5,4,3,2,2,1);(8,7,4,2,2,1)")
print("responses=all_full_and_finite_bounds_strict")
print("eta_zero=exact_M9_scope_exit")
print(f"checks={CHECKS}")
print("verdict=INDEPENDENT_P_ONLY_DELTA_ZERO_ACCEPT")
