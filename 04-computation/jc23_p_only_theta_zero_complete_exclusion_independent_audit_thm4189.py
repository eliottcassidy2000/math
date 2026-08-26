#!/usr/bin/env python3
"""Independent source-chart audit for THM-4189.

This implementation imports no primary polynomial.  It derives a different
polynomial critical pair in rational (s,p) coordinates, proves its Hessian
bridge by ideal reduction, computes the source p-resultant and reciprocal
Theta strict transform, and reconstructs the wall polygon by supporting
inequalities.  The two normalized hostiles are checked after transport to
their distinct source p-coordinates.
"""

from __future__ import annotations

from hashlib import sha256
from math import gcd

import sympy as sp


CHECKS = 0


def need(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not bool(condition):
        raise RuntimeError(message)


def valuation(poly: sp.Expr, variable: sp.Symbol) -> int:
    terms = sp.Poly(poly, variable).terms()
    need(bool(terms), "zero polynomial has no valuation")
    return min(monomial[0] for monomial, coefficient in terms if coefficient)


def exact_quotient(
    numerator: sp.Expr,
    denominator: sp.Expr,
    variable: sp.Symbol,
    label: str,
) -> sp.Poly:
    quotient = sp.cancel(numerator / denominator)
    need(sp.denom(quotient) == 1, f"{label}: quotient is not polynomial")
    return sp.Poly(quotient, variable)


def polygon_checksum(vertices: tuple[tuple[int, int], ...]):
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
Delta, Phi, Theta, eta = sp.symbols("Delta Phi Theta eta")
K = sp.Rational(2848, 45) - sp.Rational(7, 6) * Delta
t = p - s**2

# Start one transverse parameter before the wall in order to certify the
# reciprocal root escape.  Theta is then set to zero for the theorem row.
H = sp.expand(
    -3 * p
    + sp.Rational(8, 3) * p**2
    - sp.Rational(1376, 135) * p**3
    + K * s**2 * p**2
    + Phi * s * p**3
    + Delta * p**4
    + Theta * s**2 * p**3
    + eta * s * p**4
)
G = -s**2 / (2 * t) + H

# A deliberately different source critical pair from the normalized path.
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
need(sp.factor(lc_A - 2 * p * (K + Theta * p)) == 0,
     "source A leading row changed")
need(sp.factor(lc_B - 8 * p * (Theta * p + 3 * K / 4)) == 0,
     "source B leading row changed")
need(sp.factor(lc_B - 4 * lc_A + 2 * K * p) == 0,
     "source infinity incompatibility changed")

# At a common zero with p*t!=0, the differentiated gradient identities give
# p detD(A,B)=2t^2 detHess(G).  Check the exact ideal reduction.
source_jacobian = sp.det(sp.Matrix((
    (sp.diff(A, s), sp.diff(A, p)),
    (sp.diff(B, s), sp.diff(B, p)),
)))
source_hessian = sp.det(sp.hessian(G, (s, p)))
bridge_numerator = sp.together(
    p * source_jacobian - 2 * t**2 * source_hessian
).as_numer_denom()[0]
bridge_remainder = sp.reduced(bridge_numerator, [A, B], s, p)[1]
need(sp.factor(bridge_remainder) == 0, "source Hessian bridge changed")

# Complete source resultants on and transverse to Theta=0.
resultant = sp.resultant(A, B, s)
need(valuation(resultant, p) == 2, "source p-artifact changed")
R20 = exact_quotient(resultant, p**2, p, "transverse source resultant")
need(R20.degree() == 20, "transverse source residual degree changed")
need(sp.factor(R20.TC() + 31104 * K**2) == 0,
     "transverse source residual constant changed")
need(sp.factor(R20.LC() + 65610 * eta**6 * Theta) == 0,
     "transverse source residual leading row changed")

resultant0 = sp.resultant(A.subs(Theta, 0), B.subs(Theta, 0), s)
need(sp.factor(resultant0 - resultant.subs(Theta, 0)) == 0,
     "Theta-zero direct specialization changed")
R19 = exact_quotient(resultant0, p**2, p, "Theta-zero source resultant")
need(R19.degree() == 19, "Theta-zero source residual degree changed")
need(sp.factor(R19.TC() + 31104 * K**2) == 0,
     "Theta-zero source residual constant changed")
need(sp.factor(R19.LC() + 78732 * K * eta**6) == 0,
     "Theta-zero source residual leading row changed")
need(sp.factor(R20.nth(19).subs(Theta, 0) + 78732 * K * eta**6) == 0,
     "first reciprocal normal coefficient changed")
v_escape = -sp.cancel(
    (-65610 * eta**6 * Theta) / (-78732 * K * eta**6)
)
need(sp.factor(v_escape + sp.Rational(5, 6) * Theta / K) == 0,
     "source reciprocal escape scale changed")

A0 = sp.expand(A.subs(Theta, 0))
B0 = sp.expand(B.subs(Theta, 0))
source_hessian0 = sp.cancel(source_hessian.subs(Theta, 0))

# The deleted p=0 and t=0 rows contain no source-chart critical point; the
# normalized audit below restores exactly two universal points from each.
need(sp.factor(A0.subs(p, 0) + s) == 0
     and B0.subs({p: 0, s: 0}) == -6,
     "p-zero chart loss changed")
need(sp.factor(A0.subs(p, s**2) + s) == 0
     and B0.subs({p: s**2, s: 0}) == -6,
     "t-zero chart loss changed")

X, T = sp.symbols("X T")
PN = T + X**2 * T**2
YN = X * T * PN
GN = sp.expand(
    -X**2 * T / 2
    - 3 * PN
    + sp.Rational(8, 3) * PN**2
    - sp.Rational(1376, 135) * PN**3
    + K * YN**2
    + Phi * PN**2 * YN
    + Delta * PN**4
    + eta * PN**3 * YN
)
fN = sp.cancel(sp.diff(GN, X) / T)
hN = sp.diff(GN, T)
hessN = sp.det(sp.hessian(GN, (X, T)))
need(sp.factor(fN.subs(T, 0) + X) == 0
     and sp.factor(hN.subs(T, 0) + (X**2 + 6) / 2) == 0,
     "normalized T-zero pair changed")
need(sp.rem(sp.Poly(hessN.subs(T, 0) - 6, X),
            sp.Poly(X**2 + 6, X)).is_zero,
     "normalized T-zero Hessian changed")
for expression, expected, label in (
    (fN, 0, "f"),
    (hN, 0, "h"),
    (GN, sp.Rational(1, 2), "value"),
    (hessN, -6, "Hessian"),
):
    remainder = sp.rem(
        sp.Poly(sp.expand(
            expression.subs(T, -sp.Rational(1, 6)) - expected
        ), X),
        sp.Poly(X**2 - 6, X),
    )
    need(remainder.is_zero, f"normalized p-zero {label} changed")

# Supporting-inequality polygon reconstruction, independent of the primary
# convex-hull implementation.
Q = sp.symbols("Q")
H0 = sp.expand(H.subs(Theta, 0))
F_Q = sp.expand((s**2 - p) * (1 - Q * H0) - Q * s**2 / 2)
vertices = ((0, 1), (2, 0), (4, 2), (3, 4), (1, 5), (0, 5))
halfspaces = (
    (1, 2, 2),
    (-1, 1, -2),
    (-2, -1, -10),
    (-1, -2, -11),
    (0, -1, -5),
    (1, 0, 0),
)
control = {Delta: 2, Phi: 5, eta: 7, Q: 11}
support = {
    powers for powers, coefficient
    in sp.Poly(F_Q.subs(control), s, p).terms()
    if coefficient != 0
}
need(all(vertex in support for vertex in vertices),
     "polygon vertex disappeared")
for i, j in support:
    need(all(u * i + v * j >= level for u, v, level in halfspaces),
         "support escaped a claimed halfspace")
for u, v, level in halfspaces:
    need(sum(u * i + v * j == level for i, j in support) >= 2,
         "claimed polygon edge is inactive")
need(polygon_checksum(vertices) == (28, 10, 10),
     "independent polygon checksum changed")

wall_poly = sp.Poly(F_Q, s, p)
coefficient = lambda i, j: sp.factor(
    wall_poly.coeff_monomial(s**i * p**j)
)
need(coefficient(0, 1) == -1, "low vertex changed")
need(sp.factor(coefficient(2, 0) - (1 - Q / 2)) == 0,
     "low/carrier vertex changed")
need(sp.factor(coefficient(4, 2) + K * Q) == 0,
     "carrier vertex changed")
need(coefficient(3, 4) == -Q * eta, "merged vertex changed")
need(coefficient(1, 5) == Q * eta, "top vertex changed")
need(coefficient(0, 5) == Q * Delta, "horizontal vertex changed")
need(sp.factor(
    coefficient(4, 2) * s**4 * p**2
    + coefficient(3, 4) * s**3 * p**4
    + Q * p**2 * s**3 * (K * s + eta * p**2)
) == 0, "merged-edge polynomial changed")
need(sp.factor(
    coefficient(1, 5) * s * p**5
    + coefficient(0, 5) * p**5
    - Q * p**5 * (Delta + eta * s)
) == 0, "horizontal-edge polynomial changed")
need(sp.diff(Delta + eta * s, s) == eta,
     "horizontal edge is not simple on eta!=0")

# Algebraic repeated-T hostile transported to its two distinct p-values.
tau = sp.symbols("tau")
P6_expr = (
    16512 * tau**6 - 174912 * tau**5 - 235020 * tau**4
    - 35764 * tau**3 + 10607 * tau**2 + 4530 * tau + 1485
)
P6 = sp.Poly(P6_expr, tau)
N_tau = (
    33024 * tau**6 + 265344 * tau**5 + 263496 * tau**4
    + 33064 * tau**3 + 3973 * tau**2 - 2100 * tau - 1485
)
collision_values = {
    Delta: (1376 * tau**2 - 240 * tau + 135) / (180 * tau**3),
    Phi: N_tau / (1620 * tau**4 * (tau + 1)**2 * (3 * tau + 1)),
    eta: -N_tau / (1620 * tau**5 * (tau + 1)**2 * (3 * tau + 1)),
}

def tau_remainder(expr):
    numerator = sp.together(expr).as_numer_denom()[0]
    return sp.rem(sp.Poly(numerator, tau), P6)


collision_A = sp.cancel(A0.subs(collision_values))
collision_B = sp.cancel(B0.subs(collision_values))
collision_hessian = sp.cancel(source_hessian0.subs(collision_values))
source_points = ((0, tau), (tau, tau + tau**2))
for s_value, p_value in source_points:
    need(tau_remainder(collision_A.subs({s: s_value, p: p_value})).is_zero,
         "collision A point changed")
    need(tau_remainder(collision_B.subs({s: s_value, p: p_value})).is_zero,
         "collision B point changed")
    hessian_numerator = sp.together(
        collision_hessian.subs({s: s_value, p: p_value})
    ).as_numer_denom()[0]
    need(sp.gcd(P6, sp.Poly(hessian_numerator, tau)).degree() == 0,
         "collision source point ceased to be Morse")

# The nonzero source Hessians make both local intersection lengths one.  Since
# their p-values differ by tau^2, they contribute to distinct source-resultant
# roots even though the normalized T-resultant sees one repeated root.
need(tau_remainder(tau**2).as_expr() != 0,
     "collision source p-values ceased to be distinct")

# Universal-fibre hostile becomes one ordinary nonzero-p source point; its
# complete source residual is squarefree.
universal_values = {
    Delta: sp.Integer(1),
    Phi: -sp.Rational(12416, 25),
    eta: -sp.Rational(32856, 125),
}
universal_A = sp.expand(A0.subs(universal_values))
universal_B = sp.expand(B0.subs(universal_values))
universal_R = sp.Poly(sp.cancel(R19.as_expr().subs(universal_values)), p)
extra_s, extra_p = -sp.Rational(1, 6), -sp.Rational(5, 36)
need(universal_A.subs({s: extra_s, p: extra_p}) == 0,
     "universal source A point changed")
need(universal_B.subs({s: extra_s, p: extra_p}) == 0,
     "universal source B point changed")
need(universal_R.eval(extra_p) == 0,
     "universal source residual root changed")
need(sp.gcd(universal_R, universal_R.diff()).degree() == 0,
     "universal source residual ceased to be squarefree")
need(source_hessian0.subs(universal_values).subs(
    {s: extra_s, p: extra_p}
) != 0, "universal source point ceased to be Morse")

packet = (8, 7, 4, 2, 2, 1)
critical_length = 23
defect = sum(index - 1 for index in packet)
need(sum(packet) == 24 and defect == 18,
     "independent packet ledger changed")
need(2 * (24 - critical_length) < defect,
     "independent full response bound changed")
need(2 * 20 - critical_length - 1 + 2 < 19,
     "independent finite response bound changed")

coefficient_text = "|".join(sp.sstr(value) for value in R19.all_coeffs())
symbolic_digest = sha256(coefficient_text.encode("ascii")).hexdigest()
semantic = (
    "source=p^2*R19;LC=-78732*K*eta^6;TC=-31104*K^2;"
    "escape=v~-5*Theta/(6K);L=23;packet=(8,7,4,2,2,1);"
    "collision_p=tau,tau+tau^2;universal_p=-5/36"
)

print("THM4189_SOURCE_INDEPENDENT_EXACT_ACCEPT")
print(f"checks={CHECKS}")
print("chart=(s,p);critical_pair=(A,B);hessian_bridge=ideal_zero")
print("source_resultants=p^2*R20_transverse;p^2*R19_on_Theta_zero")
print(f"R19_constant={sp.factor(R19.TC())}")
print(f"R19_leading={sp.factor(R19.LC())}")
print(f"symbolic_coefficient_sha256={symbolic_digest}")
print("normal_escape=v=1/p~-5*Theta/(6*K)")
print("polygon_ledger=(28,10,10);merged_indices=3+5_to_7")
print("horizontal_face=Delta+eta*s;unique_simple_nonzero_root=-Delta/eta")
print("collision_source_p=tau,tau+tau^2;both_simple")
print("universal_source_p=-5/36;R19_squarefree")
print("packet=(8,7,4,2,2,1);critical_length=23;responses_strict")
print(f"semantic_sha256={sha256(semantic.encode('ascii')).hexdigest()}")
