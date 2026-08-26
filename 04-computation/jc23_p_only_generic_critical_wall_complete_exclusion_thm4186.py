#!/usr/bin/env python3
"""Exact source-chart certificate for THM-4186.

This implementation works in the rational (s,p) chart and does not import
the normalized-resultant code used by THM-4147.  It computes the complete
symbolic P-row source resultant, proves the source Hessian bridge and the
finite-coordinate infinity firewall, and checks two exact hostile points:
one where the normalized T-eliminant has a repeated root because two distinct
Morse points share T, and one where an extra Morse point joins T=-1/6.
"""

from __future__ import annotations

from hashlib import sha256

import sympy as sp


CHECKS = 0


def need(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def valuation(poly: sp.Expr, variable: sp.Symbol) -> int:
    return min(monomial[0] for monomial, _ in sp.Poly(poly, variable).terms())


def exact_quotient(
    numerator: sp.Expr,
    denominator: sp.Expr,
    variable: sp.Symbol,
    label: str,
) -> sp.Poly:
    quotient = sp.cancel(numerator / denominator)
    need(sp.denom(quotient) == 1, f"{label}: quotient is not polynomial")
    return sp.Poly(quotient, variable)


s, p = sp.symbols("s p")
Delta, Phi, Theta, eta = sp.symbols("Delta Phi Theta eta")
t = p - s**2
K = sp.Rational(2848, 45) - sp.Rational(7, 6) * Delta

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

# On p*t != 0 these two polynomials generate the critical ideal.
A = sp.cancel(t**2 * sp.diff(G, s) / p)
C = sp.cancel(2 * t**2 * sp.diff(G, p))
need(sp.denom(A) == sp.denom(C) == 1, "critical pair is not polynomial")
need(sp.factor(t**2 * sp.diff(G, s) - p * A) == 0,
     "first gradient identity changed")
need(sp.factor(2 * t**2 * sp.diff(G, p) - C) == 0,
     "second gradient identity changed")

# The source projection cannot lose a common zero at s=infinity when p*K is
# nonzero.  The displayed linear combination is the useful sidecar.
A_poly = sp.Poly(A, s)
C_poly = sp.Poly(C, s)
need((A_poly.degree(), C_poly.degree()) == (5, 6),
     "source s-degrees changed")
lc_A = sp.factor(A_poly.LC())
lc_C = sp.factor(C_poly.LC())
need(sp.factor(lc_A - 2 * p * (K + Theta * p)) == 0,
     "A leading row changed")
need(sp.factor(lc_C - 2 * p * (2 * K + 3 * Theta * p)) == 0,
     "C leading row changed")
need(sp.factor(3 * lc_A - lc_C - 2 * p * K) == 0,
     "source infinity incompatibility changed")

# Exact ideal-membership form of
#     p det D(A,C) = 2 t^4 det Hess(G)  modulo (A,C).
# This proves reducedness of the full source scheme under Keller--Morse; it
# does not assert squarefreeness of a chosen coordinate eliminant.
A_s, A_p = sp.diff(A, s), sp.diff(A, p)
C_s, C_p = sp.diff(C, s), sp.diff(C, p)
source_jacobian = A_s * C_p - A_p * C_s
source_hessian = sp.det(sp.hessian(G, (s, p)))
bridge_rhs = (
    A * (-4 * C * s + 4 * C_p * p * s + 2 * C_s * p - C_s * t)
    + C * (-4 * A_p * p * s - 2 * A_s * p)
)
need(sp.factor(sp.together(
    t * (2 * t**4 * source_hessian - p * source_jacobian) - bridge_rhs
)) == 0, "source Hessian bridge changed")

# Complete symbolic source resultant.  Its two endpoints are units precisely
# on the theorem's K*Theta*eta != 0 coefficient chamber.
resultant = sp.resultant(A, C, s)
need(valuation(resultant, p) == 4, "source p-artifact changed")
R20 = exact_quotient(resultant, p**4, p, "symbolic source resultant")
need(R20.degree() == 20, "symbolic residual degree changed")
need(sp.factor(R20.TC() + 31104 * K**2) == 0,
     "symbolic bottom endpoint changed")
need(sp.factor(R20.LC() + 65610 * eta**6 * Theta) == 0,
     "symbolic top endpoint changed")

# No residual root lies on p=0, and no p!=0 residual point lies on t=0.
# The p=0 saturation and singular (s,p)=(0,0) collapse are restored in the
# normalized audit as the two universal Morse pairs.
need(sp.factor(A.subs(p, 0) + s) == 0, "p=0 A row changed")
need(sp.factor(C.subs(p, 0) + s**2 * (6 * s**2 - 1)) == 0,
     "p=0 C row changed")
need(sp.factor(A.subs(p, s**2) + s) == 0, "t=0 A row changed")
need(sp.factor(C.subs(p, s**2) - s**2) == 0, "t=0 C row changed")


def specialized_source_residual(values: dict[sp.Symbol, sp.Expr]) -> tuple[sp.Poly, sp.Expr, sp.Expr]:
    specialized_A = sp.expand(A.subs(values))
    specialized_C = sp.expand(C.subs(values))
    specialized_resultant = sp.resultant(specialized_A, specialized_C, s)
    need(valuation(specialized_resultant, p) == 4,
         "specialized p-artifact changed")
    residual = exact_quotient(
        specialized_resultant, p**4, p, "specialized source resultant"
    )
    need(residual.degree() == 20, "specialized residual degree changed")
    need(sp.factor(residual.as_expr() - R20.as_expr().subs(values)) == 0,
         "direct specialization disagrees with symbolic resultant")
    return residual, specialized_A, specialized_C


# Hostile 1: two distinct reduced source points (s,p)=(0,1),(1,2) have the
# same normalized coordinate T=t=1.  The source p-eliminant stays squarefree,
# exposing the normalized repeated root as projection collision rather than
# nonreduced source geometry.
collision_values = {
    Delta: sp.Rational(1271, 180),
    Phi: sp.Rational(1733, 7560),
    Theta: -sp.Rational(206281, 7560),
    eta: -sp.Rational(1733, 7560),
}
collision_K = sp.factor(K.subs(collision_values))
need(collision_K == sp.Rational(11891, 216), "collision K changed")
need(all(value != 0 for value in (
    collision_values[Delta], collision_values[Theta],
    collision_values[eta], collision_K,
)), "collision witness left the coefficient chamber")
collision_R, collision_A, collision_C = specialized_source_residual(
    collision_values
)
need(sp.gcd(collision_R, collision_R.diff()).degree() == 0,
     "collision source residual is not squarefree")
for s_value, p_value in ((sp.Integer(0), sp.Integer(1)),
                         (sp.Integer(1), sp.Integer(2))):
    need(collision_A.subs({s: s_value, p: p_value}) == 0,
         "collision A point changed")
    need(collision_C.subs({s: s_value, p: p_value}) == 0,
         "collision C point changed")
    need(p_value - s_value**2 == 1, "collision points no longer share T=1")
    need(source_hessian.subs(collision_values).subs(
        {s: s_value, p: p_value}
    ) != 0, "collision point ceased to be Morse")
need(collision_R.eval(1) == collision_R.eval(2) == 0,
     "collision p-roots changed")

# Hostile 2: one extra reduced point (s,p)=(-1/6,-5/36) has T=-1/6.
# The source p-coordinate separates it from the two universal p=0 points.
universal_values = {
    Delta: sp.Integer(1),
    Phi: -sp.Rational(1176023, 2700),
    Theta: sp.Rational(32981, 450),
    eta: sp.Integer(1),
}
universal_K = sp.factor(K.subs(universal_values))
need(universal_K == sp.Rational(5591, 90), "universal-wall K changed")
need(all(value != 0 for value in (
    universal_values[Delta], universal_values[Theta],
    universal_values[eta], universal_K,
)), "universal witness left the coefficient chamber")
universal_R, universal_A, universal_C = specialized_source_residual(
    universal_values
)
extra_s = -sp.Rational(1, 6)
extra_p = -sp.Rational(5, 36)
need(extra_p - extra_s**2 == -sp.Rational(1, 6),
     "extra point no longer has T=-1/6")
need(universal_A.subs({s: extra_s, p: extra_p}) == 0,
     "universal-wall A point changed")
need(universal_C.subs({s: extra_s, p: extra_p}) == 0,
     "universal-wall C point changed")
need(universal_R.eval(extra_p) == 0,
     "universal-wall source root changed")
need(sp.gcd(universal_R, universal_R.diff()).degree() == 0,
     "universal-wall source residual is not squarefree")
need(source_hessian.subs(universal_values).subs(
    {s: extra_s, p: extra_p}
) != 0, "universal-wall extra point ceased to be Morse")

# The inherited packet/carrier ledger closes both responses once L=20+4.
packet = (8, 5, 4, 3, 2, 2, 1)
defect = sum(index - 1 for index in packet)
critical_length = 20 + 2 + 2
finite_n, beta, full_n = 21, 2, 25
need((sum(packet), defect, critical_length) == (25, 18, 24),
     "packet/length ledger changed")
need(2 * finite_n - critical_length - 1 + beta < finite_n - 1,
     "finite merger contradiction changed")
need(2 * (full_n - critical_length) < defect,
     "full commutator contradiction changed")

coefficient_text = "|".join(
    sp.sstr(coefficient) for coefficient in R20.all_coeffs()
)
symbolic_digest = sha256(coefficient_text.encode("ascii")).hexdigest()
semantic = (
    f"source=p^4*R20;R20_TC={sp.factor(R20.TC())};"
    f"R20_LC={sp.factor(R20.LC())};L={critical_length};"
    f"packet={packet};finite=({finite_n},{beta});full={full_n};"
    "collision=(T=1;p=1,2);universal_extra=(T=-1/6;p=-5/36)"
)
semantic_digest = sha256(semantic.encode("ascii")).hexdigest()

print("THM4186_SOURCE_CHART_EXACT_ACCEPT")
print(f"checks={CHECKS}")
print("source_resultant=p^4*R20")
print(f"R20_constant={sp.factor(R20.TC())}")
print(f"R20_leading={sp.factor(R20.LC())}")
print(f"symbolic_coefficient_sha256={symbolic_digest}")
print("source_infinity_sidecar=3LC(A)-LC(C)=2pK")
print("source_hessian_bridge=p*detD(A,C)=2*t^4*detHess(G)_mod_(A,C)")
print("collision_hostile=T=1;source_points=(0,1),(1,2);source_R20_squarefree")
print("universal_hostile=T=-1/6;source_point=(-1/6,-5/36);source_R20_squarefree")
print(f"critical_length={critical_length}")
print(f"packet={packet};defect={defect}")
print(f"finite=(n={finite_n},beta={beta});full=n={full_n}")
print(f"semantic_sha256={semantic_digest}")
