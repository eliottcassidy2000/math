#!/usr/bin/env python3
"""Exact companion for THM-3822's SL2/punctured-arm atlas gate."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.cancel(expression) == 0, message)


# The nonlinear Delone--Faddeev completion of THM-3811.
A, C, omega, theta = sp.symbols("A C omega theta")
relations = [
    omega**2 + 7 * A**2 - C * omega + A * theta,
    omega * theta - 3 * A**2 + A * C**2,
    theta**2 - 3 * A * C + C**3 - (C**2 - 3 * A) * omega + 7 * A * theta,
]
S_basis = sp.groebner(
    relations, theta, omega, C, A, order="grevlex", domain=sp.QQ
)


def reduce_S(expression: sp.Expr) -> sp.Expr:
    return sp.factor(S_basis.reduce(sp.expand(expression))[1])


D = C * omega - 3 * A * theta - 14 * A**2
m = 3 * theta + 14 * A
zero(reduce_S(D - (3 * omega**2 - 2 * C * omega + 7 * A**2)),
     "different equals f-prime")
zero(C * omega - m * A - D, "SL2 determinant numerator")

# Clear D^2 from the two intrinsic atlas equations h=A/D, k=omega/D.
eq1_numerator = A * (7 * A**2 + 3 * omega**2) - A * (D + 2 * C * omega)
eq2_numerator = (
    A * (9 * A**2 + 14 * A * omega)
    - (3 * C**2 * A**2 + C * omega**2 - omega * D)
)
zero(reduce_S(eq1_numerator), "first intrinsic atlas equation")
zero(reduce_S(eq2_numerator), "second intrinsic atlas equation")


# Converse presentation.  These three equations recover all three cubic-order
# relations after A=hD, omega=kD, theta=(m-14hD)/3.
h, k, cc, mm, dd = sp.symbols("h k C0 m D0")
E0 = cc * k - mm * h - 1
E1 = dd * (7 * h**2 + 3 * k**2) - 1 - 2 * cc * k
E2 = h * dd * (9 * h + 14 * k) - k * mm - 3 * h * cc**2
atlas_basis = sp.groebner(
    [E0, E1, E2], mm, cc, dd, k, h, order="grevlex", domain=sp.QQ
)


def reduce_atlas(expression: sp.Expr) -> sp.Expr:
    return sp.factor(atlas_basis.reduce(sp.expand(expression))[1])


A_lift = h * dd
omega_lift = k * dd
theta_lift = (mm - 14 * h * dd) / 3
lift_substitution = {A: A_lift, C: cc, omega: omega_lift, theta: theta_lift}
for index, relation in enumerate(relations, start=1):
    zero(reduce_atlas(relation.subs(lift_substitution)),
         f"converse cubic relation {index}")
zero(
    reduce_atlas(D.subs(lift_substitution) - dd),
    "converse different reconstruction",
)
zero(reduce_atlas(2 * cc * k - dd * (3 * k**2 + 7 * h**2) + 1),
     "linear C reconstruction")
zero(
    reduce_atlas(
        2 * h * theta_lift - dd * (k**2 - 7 * h**2) + 1
    ),
    "linear theta reconstruction",
)

# Eliminating D leaves one hypersurface inside SL2.  Every second-row shear
# has the same square-discriminant sidecar H(h,k).
Q = 7 * h**2 + 3 * k**2
Psi = sp.expand(
    h * (14 * k + 9 * h) * (1 + 2 * cc * k)
    - Q * (k * mm + 3 * h * cc**2)
)
zero(reduce_atlas(Psi), "SL2 compatibility follows from the two lift laws")

s = sp.symbols("s")
mm_s = mm + k * s
cc_s = cc + h * s
Psi_s = sp.expand(
    h * (14 * k + 9 * h) * (1 + 2 * cc_s * k)
    - Q * (k * mm_s + 3 * h * cc_s**2)
)
gate(sp.Poly(Psi_s, s).degree() == 2, "second-row shear equation is quadratic")
disc_s = sp.factor(sp.discriminant(Psi_s, s))
H = (
    84 * h**7
    + 36 * h**6 * k**2
    + 196 * h**6 * k
    + 84 * h**5 * k**3
    + 36 * h**5 * k**2
    + 49 * h**4 * k**4
    + 112 * h**4 * k**3
    - 12 * h**3 * k**5
    - 14 * h**2 * k**6
    + 12 * h**2 * k**5
    + k**8
)
det_basis = sp.groebner([E0], mm, cc, k, h, order="grevlex", domain=sp.QQ)
zero(det_basis.reduce(sp.expand(disc_s - 9 * H))[1],
     "universal second-row discriminant")

# Eliminating the second row gives an even smaller quadratic directly in D.
D_quadratic = (
    h**2 * Q**2 * dd**2
    + 2
    * (
        k**5 - 7 * h**2 * k**3 - 3 * h**2 * k**2
        - 6 * h**3 * k**2 - 7 * h**4
    )
    * dd
    + h**2
    - 2 * k**3
)
zero(reduce_atlas(D_quadratic), "intrinsic D quadratic")
zero(sp.discriminant(D_quadratic, dd) - 4 * k**2 * H,
     "intrinsic D quadratic discriminant")

# A useful Pell/conic completion of the same obstruction.
L = k**4 - 7 * h**2 * k**2 + 6 * h**2 * k - 6 * h**3 * k
K = 21 * h**3 + 49 * h**2 * k + 27 * h * k**2 + 49 * k**3 - 9 * k**2
zero(H - L**2 - 4 * h**4 * K, "Pell completion H=L^2+4h^4K")

# The apparently easiest way to make the Pell correction vanish cannot carry
# a nonconstant polynomial row.  For q=k/h=p/r in lowest terms, K=0 gives
# h=9p^2r/F and k=9p^3/F for this squarefree binary cubic F.
p_ratio, r_ratio = sp.symbols("p_ratio r_ratio")
F_ratio = (
    49 * p_ratio**3
    + 27 * p_ratio**2 * r_ratio
    + 49 * p_ratio * r_ratio**2
    + 21 * r_ratio**3
)
gate(
    sp.discriminant(F_ratio.subs(r_ratio, 1), p_ratio) == -27046348,
    "Pell-zero binary cubic has three distinct projective roots",
)
zero(F_ratio.subs(p_ratio, 0) - 21 * r_ratio**3,
     "Pell-zero denominator is coprime to a reduced numerator")
zero(K.subs({h: 0}) - k**2 * (49 * k - 9),
     "Pell-zero h=0 boundary is constant")


# The arm quotient is exactly G_m.  Saturating by k makes the reductions
# explicit: C=k^-1, D=k^-2, m=0 when h=0.
k_inv = sp.symbols("k_inv")
arm_basis = sp.groebner(
    [h, E0, E1, E2, k * k_inv - 1],
    mm,
    cc,
    dd,
    k_inv,
    k,
    h,
    order="grevlex",
    domain=sp.QQ,
)
zero(arm_basis.reduce(cc - k_inv)[1], "arm has C=k inverse")
zero(arm_basis.reduce(dd - k_inv**2)[1], "arm has D=k inverse square")
zero(arm_basis.reduce(mm)[1], "arm has m=0")


# Standard elementary big cell: every determinant-one completion is a lower
# shear of the displayed row, but H specializes at y=0 to an odd-degree
# polynomial and hence cannot be a square.
x, y, t = sp.symbols("x y t")
standard = {h: x, k: 1 + x * y, mm: y, cc: 1}
zero(E0.subs(standard), "standard elementary big-cell determinant")
H_standard = sp.expand(H.subs({h: x, k: 1 + x * y}))
H_standard_axis = sp.Poly(H_standard.subs(y, 0), x)
gate(H_standard_axis.degree() == 7, "standard-axis obstruction has odd degree")
gate(H_standard_axis.LC() == 84, "standard-axis obstruction has nonzero lead")


# Hyperbolic and Cohn rows reduce to one squarefree one-variable polynomial.
u = sp.symbols("u")
P = sp.Poly(H.subs({h: -1, k: u}), u)
expected_P = (
    u**8 - 14 * u**6 + 24 * u**5 + 49 * u**4
    + 28 * u**3 + 196 * u - 84
)
zero(P.as_expr() - expected_P, "hyperbolic boundary polynomial")
gate(sp.gcd(P, P.diff()).degree() == 0, "hyperbolic boundary polynomial squarefree")
gate(P.eval(0) == -84, "hyperbolic boundary polynomial avoids zero")
P_resultant = int(sp.resultant(P.as_expr(), P.diff().as_expr(), u))
gate(P_resultant % 101 == 100, "squarefree resultant mod-101 certificate")

# On the generic hyperbolic fibre xy=tau, the square sidecar is the genus-3
# curve w^2=H(tau-1,Z).  Its t=0 discriminant is already nonzero, so the
# generic polynomial is squarefree over k(tau).
tau, Z = sp.symbols("tau Z")
P_tau = sp.Poly(H.subs({h: tau - 1, k: Z}), Z)
gate(P_tau.degree() == 8 and P_tau.LC() == 1,
     "generic hyperbolic sidecar has monic degree eight")
zero(P_tau.as_expr().subs(tau, 0) - P.as_expr().subs(u, Z),
     "generic hyperbolic sidecar specializes to P")
zero(sp.discriminant(P_tau.as_expr(), Z).subs(tau, 0) - P_resultant,
     "generic hyperbolic discriminant specialization")

# The natural hyperbolic unit row and the first Cohn row both specialize to
# P(c*x^n); the human proof handles every n>=1 by the chain rule.
n_control = 5
h_hyperbolic = x * y - 1
k_hyperbolic = x**n_control
m_hyperbolic = sum((x * y) ** j for j in range(n_control))
C_hyperbolic = y**n_control
zero(k_hyperbolic * C_hyperbolic - h_hyperbolic * m_hyperbolic - 1,
     "hyperbolic unit-row determinant control")
zero(
    H.subs({h: h_hyperbolic, k: k_hyperbolic}).subs(y, 0)
    - P.as_expr().subs(u, x**n_control),
    "hyperbolic row specializes to P(x^n)",
)

h_cohn = 2 * x * t - 1
k_cohn = 4 * t**2
m_cohn = 1 + 2 * x * t
C_cohn = x**2
zero(k_cohn * C_cohn - h_cohn * m_cohn - 1, "Cohn-row determinant")
zero(
    H.subs({h: h_cohn, k: k_cohn}).subs(x, 0)
    - P.as_expr().subs(u, 4 * t**2),
    "Cohn row specializes to P(4t^2)",
)


semantic = {
    "surface": "THM3811 U=Spec S[A/D,omega/D]",
    "matrix": "[[k,h],[m,C]] with h=A/D,k=omega/D,m=3theta+14A and determinant 1",
    "lift": "D(7h^2+3k^2)=1+2Ck; hD(9h+14k)=km+3hC^2",
    "plane": "a plane lift is etale iff Jac(hD,C) is a nonzero scalar",
    "arm": "V(h)=G_m; every etale plane pullback component has a nonconstant unit",
    "discriminant": "Disc_s(Psi)=9H(h,k); H=L^2+4h^4K",
    "D_quadratic": "h^2(3k^2+7h^2)^2D^2+2BD+h^2-2k^3=0; discriminant=4k^2H",
    "closed": "standard elementary first row; every h=xy-1 plane atlas; first Cohn row",
    "open": "nonstandard punctured arms and interacting/non-elementary Cohn words",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()
source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
     "inactive Python assert")

print("theorem=THM-3822-nonlinear-cubic-plane-atlas-sl2-and-punctured-arm-gate")
print("surface=THM3811_U=Spec(S[A/D,omega/D])")
print("sl2=matrix_[k,h;m,C];det=1;m=3theta+14A")
print("lift=D*(7h^2+3k^2)=1+2Ck;hD*(9h+14k)=km+3hC^2")
print("plane_etale=Jac(hD,C)_nonzero_scalar")
print("arm=V(h)=G_m;etale_plane_pullback_requires_nonconstant_units")
print("D_quadratic=intrinsic_degree_2;discriminant=4k^2H")
print("sidecar=Disc_second_row=9H;H=L^2+4h^4K")
print("closed=standard_elementary_big_cell;all_h=xy-1_plane_atlases;first_Cohn_row")
print("open=nonstandard_punctured_arms;interacting_non_elementary_Cohn_words")
print(f"P_resultant_mod101={P_resultant % 101}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
