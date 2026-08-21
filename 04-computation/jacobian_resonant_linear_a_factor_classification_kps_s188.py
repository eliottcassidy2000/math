#!/usr/bin/env python3
"""Exact algebra for THM-3565's resonant linear-a factor classification."""

from __future__ import annotations

import sympy as sp


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(f"failed truth gate: {label}")


a, b, x = sp.symbols("a b x")
A, u = sp.symbols("A u", nonzero=True)
r, v, w = sp.symbols("r v w")


def core(phi: sp.Expr) -> tuple[sp.Expr, sp.Expr]:
    L = sp.expand(27 * a**2 * phi**2 + 18 * a * b * phi + 16 * a - b**3 * phi - b**2)
    E = sp.expand(L * x**3 + (4 + 3 * b * phi) * x + 2 * phi)
    return L, E


# -------------------------------------------------------------------------
# Rational-root branch I: constant numerator, linear denominator.
# -------------------------------------------------------------------------
A_lead = -sp.Rational(2, 27) / u**3
phi_r = A_lead * (a - r)
L_r, _ = core(phi_r)
D1 = a - v
identity_1 = sp.Poly(
    sp.cancel(u**3 * L_r + u * (4 + 3 * b * phi_r) * D1**2 + 2 * phi_r * D1**3),
    a,
)
require(identity_1.degree() == 3, "branch-I leading cancellation")

# The a^3 equation fixes v.  The remaining top equation is a product of two
# candidates; the second candidate is killed by the two lower coefficients.
v_1 = b * u / 2 + r / 3
branch_1_tail = sp.Poly(identity_1.as_expr().subs(v, v_1), a)
coeff_1 = dict(zip(range(branch_1_tail.degree(), -1, -1), branch_1_tail.all_coeffs()))
R = -3 * b * u + 2 * r + 18 * u**2
S = 3 * b * u + 2 * r - 18 * u**2
require(sp.expand(coeff_1[2] * 81 * u**3 + R * S) == 0, "branch-I split equation")

r_S = 9 * u**2 - sp.Rational(3, 2) * b * u
c1_S = sp.factor(coeff_1[1].subs(r, r_S))
c0_S = sp.factor(coeff_1[0].subs(r, r_S))
expected_c1_S = sp.Rational(2, 27) * (-b + 6 * u) ** 2 * (b + 12 * u)
expected_c0_S = b * u * (-b + 6 * u) ** 2 * (b + 6 * u) / 9
require(sp.expand(c1_S - expected_c1_S) == 0, "branch-I false split first lower equation")
require(sp.expand(c0_S - expected_c0_S) == 0, "branch-I false split second lower equation")
gcd_S = sp.Poly(c1_S, b).gcd(sp.Poly(c0_S, b)).monic().as_expr()
require(sp.factor(gcd_S - (b - 6 * u) ** 2) == 0, "branch-I false split gcd")
require(sp.expand(R.subs(r, r_S).subs(b, 6 * u)) == 0, "branch-I false split returns main branch")

r_main = sp.Rational(3, 2) * b * u - 9 * u**2
v_main = b * u - 3 * u**2
phi_main_u = sp.factor(A_lead * (a - r_main))
root_main_u = u / (a - v_main)
L_main_u, E_main_u = core(phi_main_u)
require(sp.factor(E_main_u.subs(x, root_main_u)) == 0, "branch-I rational root")


# -------------------------------------------------------------------------
# Rational-root branch II: numerator phi, quadratic denominator.
# Coprimality with phi is load-bearing.
# -------------------------------------------------------------------------
D2 = a**2 + v * a + w
identity_2 = sp.Poly(
    sp.cancel(
        u**3 * (a - r) ** 3 * L_r
        + u * (a - r) * (4 + 3 * b * phi_r) * D2**2
        + 2 * phi_r * D2**3
    ),
    a,
)
require(identity_2.degree() == 6, "branch-II leading cancellation")

v_2 = -b * u / 2 - 4 * r / 3
w_2 = (9 * b**2 * u**2 + 18 * b * r * u - 108 * b * u**3 + 8 * r**2 + 324 * u**4) / 36
tail_2 = sp.Poly(sp.cancel(identity_2.as_expr().subs({v: v_2, w: w_2})), a)
require(tail_2.degree() == 4, "branch-II sequential top equations")
for coefficient in tail_2.all_coeffs():
    require(sp.rem(sp.together(coefficient).as_numer_denom()[0], R, r) == 0, "branch-II common resonance factor")

# On R=0 the quadratic denominator contains a-r, so this is merely branch I
# written non-primitively and is forbidden by gcd(numerator,denominator)=1.
D2_resonant = sp.factor(D2.subs(v, v_2).subs(w, w_2).subs(r, r_main))
require(
    sp.factor(D2_resonant - (a - r_main) * (a - v_main)) == 0,
    "branch-II R=0 denominator cancellation",
)

# Away from R=0, scale p=b/u and q=r/u^2.  The exact reduced Groebner basis
# of the five residual coefficients is a line and a three-point polynomial.
p, q = sp.symbols("p q")
inner_equations = []
for coefficient in tail_2.all_coeffs():
    quotient = sp.cancel(coefficient / R)
    normalized = sp.cancel(quotient.subs({u: 1, b: p, r: q}))
    inner_equations.append(normalized.as_numer_denom()[0])
G = sp.groebner(inner_equations, q, p, order="lex", domain=sp.QQ)
groebner_monic = [sp.Poly(poly, q, p).monic().as_expr() for poly in G.polys]
expected_G = [q - p**2 / 8 + 3 * p / 4, p**3 - 12 * p**2 - 36 * p + 432]
require(len(groebner_monic) == 2, "branch-II Groebner length")
require(all(sp.expand(g - e) == 0 for g, e in zip(groebner_monic, expected_G)), "branch-II Groebner basis")
require(sp.factor(expected_G[1]) == (p - 12) * (p - 6) * (p + 6), "branch-II three points")

candidate_rows = []
for p_value in (12, 6, -6):
    q_value = sp.Rational(p_value**2 - 6 * p_value, 8)
    R_value = sp.factor(-3 * p_value + 2 * q_value + 18)
    candidate_rows.append((p_value, q_value, R_value))
require(candidate_rows == [(12, 9, 0), (6, 0, 0), (-6, 9, 54)], "branch-II candidate atlas")

# The only coprime branch is p=-6,q=9.  Since b=-6u, its phi has unavoidable
# b-denominators and cannot belong to C[a,b].
u_exceptional = -b / 6
r_exceptional = 9 * u_exceptional**2
phi_exceptional = sp.factor(A_lead.subs(u, u_exceptional) * (a - r_exceptional))
require(sp.cancel(phi_exceptional - (16 * a / b**3 - 4 / b)) == 0, "branch-II exceptional rational graph")


# -------------------------------------------------------------------------
# Polynomial descent and explicit pullback factor.
# -------------------------------------------------------------------------
h = sp.symbols("h", nonzero=True)
phi_h = -2 * h**3 * a + b * h**2 - 2 * h
L_h, E_h = core(phi_h)
D_h = 3 * a * h**2 - b * h + 1
C_h = 12 * a**2 * h**2 - 4 * a * b * h + 16 * a - b**2
quadratic_h = D_h * C_h * x**2 + h * C_h * x + 2 * (2 * a * h**2 - b * h + 2)
require(sp.factor(L_h - D_h**2 * C_h) == 0, "h-family leading coefficient factor")
require(sp.expand(E_h - (D_h * x - h) * quadratic_h) == 0, "h-family core factorization")
quadratic_discriminant = sp.factor(sp.discriminant(quadratic_h, x))
expected_quadratic_discriminant = -(6 * a * h**2 - 3 * b * h + 4) ** 2 * C_h
require(
    sp.expand(quadratic_discriminant - expected_quadratic_discriminant) == 0,
    "h-family residual quadratic discriminant",
)
C_discriminant = sp.factor(sp.discriminant(C_h, a))
require(
    sp.expand(C_discriminant - 64 * (b**2 * h**2 - 2 * b * h + 4)) == 0,
    "h-family residual nonsquare gate",
)

xs, ys, zs, H = sp.symbols("x_s y_s z_s H")
us = 1 + xs * ys
F1 = sp.expand(us**3 * zs + ys**2 * us * (4 + 3 * xs * ys))
F2 = sp.expand(ys + 3 * xs * us**2 * zs + 3 * xs * ys**2 * (4 + 3 * xs * ys))
F3 = sp.expand(2 * xs - 3 * xs**2 * ys - xs**3 * zs)
pull_H = sp.expand(F3 - 2 * H**3 * F1 + H**2 * F2 - 2 * H)
source_linear = xs - us * H
source_cofactor_fraction = sp.cancel(pull_H / source_linear)
source_cofactor, source_denominator = source_cofactor_fraction.as_numer_denom()
require(source_denominator == 1, "universal source cofactor polynomiality")
source_cofactor = sp.factor(source_cofactor)
require(sp.expand(pull_H - source_linear * source_cofactor) == 0, "universal source factorization")

# Collision-compatible total-degree <=2: deg(h)=0, then h=+-2.  These are
# exactly THM-3559's two affine Kummer rows.
h0 = sp.symbols("h0", nonzero=True)
collision_value = sp.factor(phi_h.subs({a: -sp.Rational(1, 4), b: 0, h: h0}))
require(collision_value == h0 * (h0 - 2) * (h0 + 2) / 2, "collision condition")
phi_plus = sp.expand(phi_h.subs(h, 2))
phi_minus = sp.expand(phi_h.subs(h, -2))
require(phi_plus == -16 * a + 4 * b - 4, "positive Kummer row")
require(phi_minus == 16 * a + 4 * b + 4, "negative Kummer row")

print("THM-3565 resonant linear-a target-graph factor audit")
print("branch I: phi=-2*a/(27*u^3)+(b-6*u)/(9*u^2)")
print("branch II normalized Groebner basis:")
for polynomial in expected_G:
    print(sp.factor(polynomial))
print("branch II candidates (p,q,R):", candidate_rows)
print("only coprime candidate: phi=16*a/b^3-4/b, not polynomial in b")
print("polynomial descent: u=1/(3*h), h in C[b] nonzero")
print("factor family: phi=-2*h^3*a+b*h^2-2*h")
print("core linear factor: (3*a*h^2-b*h+1)*x-h")
print("residual discriminant: -(6*a*h^2-3*b*h+4)^2*(12*a^2*h^2-4*a*b*h+16*a-b^2)")
print("for polynomial h, the residual quadratic is irreducible")
print("source pullback factor: x-(1+x*y)*h(F2)")
print("collision-compatible degree<=2 rows: h=2 and h=-2, exactly the affine Kummer pair")
print("all active truth gates passed")
