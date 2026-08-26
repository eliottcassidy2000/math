#!/usr/bin/env python3
"""Independent rational-source exact audit for THM-4192.

It imports no normalized critical resultant from the primary path and uses a
different critical pair.
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
Phi, Theta, eta = sp.symbols("Phi Theta eta")
Delta0 = sp.Rational(5696, 105)
t = p - s**2

H = sp.expand(
    -3 * p
    + sp.Rational(8, 3) * p**2
    - sp.Rational(1376, 135) * p**3
    + Phi * s * p**3
    + Delta0 * p**4
    + Theta * s**2 * p**3
    + eta * s * p**4
)
G = -s**2 / (2 * t) + H

# This source pair differs from the primary normalized pair.  It is designed
# so one polynomial has s-degree at most two on the leading stratum.
A = sp.cancel((-s * p + t**2 * sp.diff(H, s)) / p)
C0 = sp.expand(s**2 + 2 * t**2 * sp.diff(H, p))
B = sp.cancel((C0 + s * A) / t**2)
need(sp.denom(A) == 1 and sp.denom(B) == 1,
     "source critical pair is not polynomial")
need(sp.factor(t**2 * sp.diff(G, s) - p * A) == 0,
     "first source-gradient identity changed")
need(sp.factor(2 * t**2 * sp.diff(G, p) - (t**2 * B - s * A)) == 0,
     "second source-gradient identity changed")

# Exact ideal reduction of p detD(A,B)=2t^2 detHess(G) at common zeros.
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

# The saturated pair has no common point on either deleted coordinate row.
need(sp.factor(A.subs(p, 0) + s) == 0 and B.subs({p: 0, s: 0}) == -6,
     "p=0 chart row changed")
need(sp.factor(A.subs(p, s**2) + s) == 0
     and B.subs({p: s**2, s: 0}) == -6,
     "t=0 chart row changed")

# Stratum 1: Theta!=0, Phi arbitrary.
need((sp.degree(A, s), sp.degree(B, s)) == (5, 2),
     "Theta-nonzero source degrees changed")
need(sp.factor(sp.Poly(A, s).LC() - 2 * Theta * p**2) == 0,
     "Theta-nonzero A leading row changed")
need(sp.factor(sp.Poly(B, s).LC() - 8 * Theta * p**2) == 0,
     "Theta-nonzero B leading row changed")
res18 = sp.resultant(A, B, s)
need(valuation(res18, p) == 4 and sp.degree(res18, p) == 22,
     "Theta-nonzero source artifact changed")
R18 = exact_quotient(res18, p**4, p, "Theta-nonzero source resultant")
need(R18.degree() == 18, "R18 degree changed")
need(sp.factor(R18.TC() + 31104 * Theta**2) == 0,
     "R18 constant changed")
need(sp.factor(R18.LC() + 65610 * Theta * eta**6) == 0,
     "R18 leading row changed")

# Stratum 2: Theta=0, Phi!=0.  The degree drops in s, so compute the direct
# specialized resultant rather than specializing a fixed Sylvester matrix.
A0, B0 = sp.expand(A.subs(Theta, 0)), sp.expand(B.subs(Theta, 0))
need(res18.subs(Theta, 0) == 0,
     "row-A Sylvester determinant no longer exposes the degree-drop trap")
need((sp.degree(A0, s), sp.degree(B0, s)) == (4, 1),
     "Theta-zero source degrees changed")
need(sp.factor(sp.Poly(A0, s).LC() - p**2 * (Phi + eta * p)) == 0,
     "Theta-zero A leading row changed")
need(sp.factor(sp.Poly(B0, s).LC() - p**2 * (7 * Phi + 9 * eta * p)) == 0,
     "Theta-zero B leading row changed")
need(sp.factor(sp.resultant(Phi + eta * p, 7 * Phi + 9 * eta * p, p)
               + 2 * Phi * eta) == 0,
     "Theta-zero infinity incompatibility changed")
res15 = sp.resultant(A0, B0, s)
need(valuation(res15, p) == 2 and sp.degree(res15, p) == 17,
     "Theta-zero source artifact changed")
R15 = exact_quotient(res15, p**2, p, "Theta-zero source resultant")
need(R15.degree() == 15, "R15 degree changed")
need(sp.factor(R15.TC() - 1296 * Phi) == 0,
     "R15 constant changed")
need(sp.factor(R15.LC() - 6561 * eta**5) == 0,
     "R15 leading row changed")

# Stratum 3: Theta=Phi=0.  Recompute once more from the specialized source
# pair; do not infer its p-adic saturation from the row-II endpoint formula.
A00, B00 = sp.expand(A0.subs(Phi, 0)), sp.expand(B0.subs(Phi, 0))
need((sp.degree(A00, s), sp.degree(B00, s)) == (4, 1),
     "terminal source degrees changed")
need(sp.factor(sp.Poly(A00, s).LC() - eta * p**3) == 0,
     "terminal A leading row changed")
need(sp.factor(sp.Poly(B00, s).LC() - 9 * eta * p**3) == 0,
     "terminal B leading row changed")
res14 = sp.resultant(A00, B00, s)
need(sp.factor(res14 - res15.subs(Phi, 0)) == 0,
     "direct Phi-zero recomputation disagrees")
need(valuation(res14, p) == 3 and sp.degree(res14, p) == 17,
     "terminal source artifact changed")
R14 = exact_quotient(res14, p**3, p, "terminal source resultant")
need(R14.degree() == 14, "R14 degree changed")
need(sp.factor(R14.TC() - 1296 * eta) == 0,
     "R14 constant changed")
need(sp.factor(R14.LC() - 6561 * eta**5) == 0,
     "R14 leading row changed")

# Disjoint exact controls attack accidental nonreduced residuals.  The proof
# does not assume their squarefreeness.
source_controls = (
    (sp.Poly(R18.as_expr().subs({Theta: 2, Phi: 3, eta: 5}), p), "R18"),
    (sp.Poly(R15.as_expr().subs({Phi: 2, eta: 3}), p), "R15"),
    (sp.Poly(R14.as_expr().subs(eta, 2), p), "R14"),
)
for residual, label in source_controls:
    need(sp.gcd(residual, residual.diff()).degree() == 0,
         f"{label} source squarefree control failed")

# The normalized chart restores exactly two T=0 and two P=0 universal points.
X, T = sp.symbols("X T")
PN = T + X**2 * T**2
YN = X * T * PN
GN = sp.expand(
    -X**2 * T / 2
    - 3 * PN
    + sp.Rational(8, 3) * PN**2
    - sp.Rational(1376, 135) * PN**3
    + Phi * PN**2 * YN
    + Delta0 * PN**4
    + Theta * PN * YN**2
    + eta * PN**3 * YN
)
fN = sp.cancel(sp.diff(GN, X) / T)
hN = sp.diff(GN, T)
hessN = sp.det(sp.hessian(GN, (X, T)))
need(sp.factor(fN.subs(T, 0) + X) == 0
     and sp.factor(hN.subs(T, 0) + (X**2 + 6) / 2) == 0,
     "normalized T=0 pair changed")
need(sp.rem(sp.Poly(hessN.subs(T, 0) - 6, X),
            sp.Poly(X**2 + 6, X)).is_zero,
     "normalized T=0 Hessian changed")
for expression, expected, label in (
    (fN, 0, "f"),
    (hN, 0, "h"),
    (GN, sp.Rational(1, 2), "value"),
    (hessN, -6, "Hessian"),
):
    need(sp.rem(
        sp.Poly(sp.expand(
            expression.subs(T, -sp.Rational(1, 6)) - expected
        ), X),
        sp.Poly(X**2 - 6, X),
    ).is_zero, f"normalized P=0 {label} changed")

critical_lengths = (R18.degree() + 4, R15.degree() + 4, R14.degree() + 4)
need(critical_lengths == (22, 19, 18), "source critical lengths changed")

# Independent supporting-halfspace reconstruction of the three polygons.
Q = sp.symbols("Q")
F_Q = sp.expand((s**2 - p) * (1 - Q * H) - Q * s**2 / 2)
rows = (
    (
        "Theta_nonzero",
        F_Q,
        {Theta: 7, Phi: 5, eta: 11, Q: 13},
        ((0, 1), (2, 0), (4, 3), (3, 4), (1, 5), (0, 5)),
        ((1, 2, 2), (-3, 2, -6), (-1, -1, -7),
         (-1, -2, -11), (0, -1, -5), (1, 0, 0)),
        (27, 9, 10), (8, 5, 5, 4, 1), 18,
    ),
    (
        "Theta_zero_Phi_nonzero",
        sp.expand(F_Q.subs(Theta, 0)),
        {Phi: 5, eta: 11, Q: 13},
        ((0, 1), (2, 0), (3, 3), (3, 4), (1, 5), (0, 5)),
        ((1, 2, 2), (-3, 1, -6), (-1, 0, -3),
         (-1, -2, -11), (0, -1, -5), (1, 0, 0)),
        (23, 9, 8), (8, 4, 4, 2, 1), 14,
    ),
    (
        "Theta_Phi_zero",
        sp.expand(F_Q.subs({Theta: 0, Phi: 0})),
        {eta: 11, Q: 13},
        ((0, 1), (2, 0), (3, 4), (1, 5), (0, 5)),
        ((1, 2, 2), (-4, 1, -8), (-1, -2, -11),
         (0, -1, -5), (1, 0, 0)),
        (22, 8, 8), (8, 5, 4, 1), 14,
    ),
)

for name, polynomial, control, vertices, halfspaces, checksum, packet, defect in rows:
    support = {
        powers for powers, coefficient
        in sp.Poly(polynomial.subs(control), s, p).terms()
        if coefficient
    }
    need(all(vertex in support for vertex in vertices),
         f"{name} polygon vertex disappeared")
    for i, j in support:
        need(all(u * i + v * j >= level for u, v, level in halfspaces),
             f"{name} support escaped a claimed halfspace")
    for u, v, level in halfspaces:
        need(sum(u * i + v * j == level for i, j in support) >= 2,
             f"{name} claimed edge is inactive")
    need(polygon_checksum(vertices) == checksum,
         f"{name} polygon checksum changed")
    source_indices = tuple(sorted(
        (u + v - level for u, v, level in halfspaces[:-1]), reverse=True
    ))
    need(source_indices == packet, f"{name} packet changed")
    need(sum(index - 1 for index in packet) == defect,
         f"{name} defect changed")

packets = ((8, 5, 5, 4, 1), (8, 4, 4, 2, 1), (8, 5, 4, 1))
degrees = tuple(sum(packet) for packet in packets)
gaps = tuple(n - length for n, length in zip(degrees, critical_lengths))
defects = tuple(sum(index - 1 for index in packet) for packet in packets)
need((degrees, gaps, defects) == ((23, 19, 18), (1, 0, 0), (18, 14, 14)),
     "response consequence ledger changed")
need(all(2 * gap < defect for gap, defect in zip(gaps, defects)),
     "commutator contradictions changed")

# Projection hostile A transported to two distinct source p-fibres.  The
# cubic interval isolates the same real tau used in the normalized audit.
tau = sp.symbols("tau")
C3 = sp.Poly(68352 * tau**3 - 9632 * tau**2 + 1680 * tau - 945, tau)
left, right = sp.Rational(51, 200), sp.Rational(32, 125)
need(C3.eval(left) < 0 and C3.eval(right) > 0
     and C3.count_roots(left, right) == 1,
     "collision interval isolation changed")
B_tau = (
    68352 * tau**6 + 205056 * tau**5 + 195424 * tau**4
    + 49088 * tau**3 - 7952 * tau**2 + 1680 * tau - 315
)
eta_tau = sp.factor(
    -2 * B_tau / (315 * tau**4 * (tau + 1)**2 * (5 * tau + 2))
)
Phi_tau = sp.factor(-eta_tau * tau)
Theta_tau = sp.factor((
    136704 * tau**7 + 410112 * tau**6 + 390848 * tau**5
    + 98176 * tau**4 - 15904 * tau**3 + 3360 * tau**2
    + 945 * tau + 630
) / (630 * tau**4 * (tau + 1)**2 * (5 * tau + 2)))
collision_values = {Phi: Phi_tau, Theta: Theta_tau, eta: eta_tau}
for label, expression in collision_values.items():
    numerator, denominator = sp.together(expression).as_numer_denom()
    need(sp.gcd(C3, sp.Poly(numerator, tau)).degree() == 0,
         f"collision {label} can vanish")
    need(sp.gcd(C3, sp.Poly(denominator, tau)).degree() == 0,
         f"collision {label} denominator can vanish")


def tau_remainder(expr):
    numerator, denominator = sp.together(expr).as_numer_denom()
    need(sp.gcd(C3, sp.Poly(denominator, tau)).degree() == 0,
         "collision denominator meets cubic")
    return sp.rem(sp.Poly(numerator, tau), C3)


collision_A = sp.cancel(A.subs(collision_values))
collision_B = sp.cancel(B.subs(collision_values))
collision_hessian = sp.cancel(source_hessian.subs(collision_values))
source_points = ((0, tau), (tau, tau + tau**2))
for s_value, p_value in source_points:
    need(tau_remainder(collision_A.subs({s: s_value, p: p_value})).is_zero,
         "collision A point changed")
    need(tau_remainder(collision_B.subs({s: s_value, p: p_value})).is_zero,
         "collision B point changed")
    hessian_numerator = sp.together(
        collision_hessian.subs({s: s_value, p: p_value})
    ).as_numer_denom()[0]
    need(sp.gcd(C3, sp.Poly(hessian_numerator, tau)).degree() == 0,
         "collision source point ceased to be Morse")
need(not tau_remainder(tau**2).is_zero,
     "collision source p-values ceased to be distinct")

# Projection hostile B becomes an ordinary nonzero-p source point in row II.
universal_values = {
    Phi: sp.Rational(201808, 1575),
    eta: sp.Rational(2766784, 875),
}
universal_A = sp.expand(A0.subs(universal_values))
universal_B = sp.expand(B0.subs(universal_values))
universal_R15 = sp.Poly(R15.as_expr().subs(universal_values), p)
extra_s, extra_p = -sp.Rational(1, 6), -sp.Rational(5, 36)
need(universal_A.subs({s: extra_s, p: extra_p}) == 0,
     "universal source A point changed")
need(universal_B.subs({s: extra_s, p: extra_p}) == 0,
     "universal source B point changed")
need(universal_R15.eval(extra_p) == 0,
     "universal source residual root changed")
need(sp.gcd(universal_R15, universal_R15.diff()).degree() == 0,
     "universal source residual ceased to be squarefree")
need(source_hessian.subs({Theta: 0, **universal_values}).subs(
    {s: extra_s, p: extra_p}
) == -sp.Rational(27296015083, 26040609),
     "universal source Hessian changed")
need(G.subs({Theta: 0, **universal_values}).subs(
    {s: extra_s, p: extra_p}
) == sp.Rational(42257, 91854),
     "universal source value changed")

coefficient_digest = sha256("|".join(
    sp.sstr(coefficient)
    for polynomial in (R18, R15, R14)
    for coefficient in polynomial.all_coeffs()
).encode("ascii")).hexdigest()
semantic = (
    "source_pair=(A,B);resultants=(p^4R18,p^2R15,p^3R14);"
    "lengths=(22,19,18);packets=((8,5,5,4,1),(8,4,4,2,1),(8,5,4,1));"
    "degrees=(23,19,18);gaps=(1,0,0);defects=(18,14,14);"
    "degree_drop=recompute_Theta0_and_Phi0;hostiles=collision,universal"
)

print("THM4192_K_ZERO_SOURCE_INDEPENDENT_EXACT_ACCEPT")
print(f"checks={CHECKS}")
print("source_pair=(A,B);hessian_bridge=ideal_zero")
print("source_resultants=p^4*R18;p^2*R15;p^3*R14")
print("source_constants=-31104*Theta^2;1296*Phi;1296*eta")
print("source_leading=-65610*Theta*eta^6;6561*eta^5;6561*eta^5")
print("critical_lengths=22,19,18")
print("polygons=(27,9,10);(23,9,8);(22,8,8)")
print("packets=(8,5,5,4,1);(8,4,4,2,1);(8,5,4,1)")
print("degrees=23,19,18;gaps=1,0,0;defects=18,14,14")
print("degree_drop_firewall=res18|Theta=0_is_zero;res15,res14_direct")
print("projection_collision_source_p=tau,tau+tau^2;both_Morse")
print("projection_universal_source_p=-5/36;R15_squarefree")
print(f"symbolic_coefficient_sha256={coefficient_digest}")
print(f"semantic_sha256={sha256(semantic.encode('ascii')).hexdigest()}")
