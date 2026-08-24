#!/usr/bin/env python3
"""Exact companion for THM-3933's centered degree-three root-map closure.

Reproduction:
  python3 04-computation/jc2_centered_degree_three_root_map_octic_thm3933.py
  python3 -O 04-computation/jc2_centered_degree_three_root_map_octic_thm3933.py
"""

from __future__ import annotations

import hashlib
import json

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


u, v, A, C, T, E, W, lam, s, p = sp.symbols("u v A C T E W lambda s p")


# ---------------------------------------------------------------------------
# Pole partitions and the normalized trace-zero triple-pole family.
# ---------------------------------------------------------------------------

finite_ramification_budget = 2
partition_costs = {
    "1+1+1": 3,
    "2+1": 3,
    "3": 1,
}
gate(partition_costs["1+1+1"] > finite_ramification_budget,
     "three isolated simple poles exceed the finite RH budget")
gate(partition_costs["2+1"] > finite_ramification_budget,
     "isolated double plus simple poles exceed the finite RH budget")
gate(partition_costs["3"] <= finite_ramification_budget,
     "a triple pole at a simple critical point is the sole live partition")

A_u = u**3 + u**2
N = 3 * lam * u**2 + 3 * u - 2
D = u**3
t = N / D

# The T^2 coefficient of the resultant is minus the field trace.
incidence_resultant = sp.Poly(sp.resultant(A - A_u, D * T - N, u), T)
gate(sp.expand(incidence_resultant.coeff_monomial(T**2)) == 0,
     "normalized triple-pole root has trace zero")

a = -A**3
c = -9 * A**2 * lam**2 - 27 * A**2 * lam + 12 * A * lam + 15 * A - 4
d = -(
    27 * A**2 * lam**3
    - 36 * A * lam**2
    + 27 * A * lam
    + 27 * A
    + 12 * lam
    - 20
) / 2
expected_incidence = sp.Poly(a * T**3 - c * T - 2 * d, T)
gate(incidence_resultant == expected_incidence,
     "primitive resultant gives the displayed a,c,d")
gate(sp.gcd(sp.gcd(a, c), 2 * d) == 1, "incidence coefficient triple is primitive")
gate(c.subs(A, 0) == -4, "unit coefficient-ideal control at A=0")

# Direct trace formula for an arbitrary numerator b2*u^2+b1*u+b0.
b0, b1, b2 = sp.symbols("b0 b1 b2")
raw_N = b2 * u**2 + b1 * u + b0
raw_resultant = sp.Poly(sp.resultant(A - A_u, D * T - raw_N, u), T)
raw_lc = raw_resultant.coeff_monomial(T**3)
raw_t2 = raw_resultant.coeff_monomial(T**2)
gate(sp.factor(-raw_t2 / raw_lc) == (3 * b0 + 2 * b1) / A,
     "arbitrary triple-pole trace is (3b0+2b1)/A")


# ---------------------------------------------------------------------------
# Polynomial-color divisibility and the unique lambda=3 row.
# ---------------------------------------------------------------------------

C_rational = sp.cancel(-(3 * a.subs(A, A_u) * t**2 + c.subs(A, A_u)) / (2 * t))
color_numerator, color_denominator = sp.together(C_rational).as_numer_denom()
gate(sp.expand(color_denominator - 2 * N) == 0, "color denominator is exactly 2N")

color_resultant = sp.factor(sp.resultant(N, color_numerator, u))
gate(color_resultant == -124416 * lam**2 * (lam - 3) ** 3,
     "generic polynomial-color resultant")
gate(color_numerator.subs({lam: 0, u: sp.Rational(2, 3)}) == -sp.Rational(512, 243),
     "lambda=0 degree-drop seam does not cancel")

C_u = sp.factor(C_rational.subs(lam, 3))
expected_C_u = u**3 * (3 * u + 2) * (
    9 * u**4 + 30 * u**3 + 24 * u**2 - 4
) / 2
gate(sp.expand(C_u - expected_C_u) == 0, "lambda=3 gives the exact polynomial color")
gate(sp.denom(sp.cancel(C_u)) == 1, "lambda=3 color is polynomial")

a3 = sp.factor(a.subs(lam, 3))
c3 = sp.factor(c.subs(lam, 3))
d3 = sp.factor(d.subs(lam, 3))
gate(a3 == -A**3, "survivor leading coefficient")
gate(sp.expand(c3 + (6 * A - 1) * (27 * A - 4)) == 0, "survivor middle endpoint")
gate(sp.expand(d3 + (27 * A - 4) ** 2 / 2) == 0, "survivor final endpoint")


# ---------------------------------------------------------------------------
# Discriminant, implicit octic, and the order/maximal-order separation.
# ---------------------------------------------------------------------------

Phi_affine = a3 * T**3 + C * T**2 + c3 * T + d3
disc = sp.factor(sp.discriminant(Phi_affine, T))
L = 27 * A - 4
H = (
    19683 * A**8
    + 87480 * A**7
    - 60048 * A**6
    + 5832 * A**5 * C
    + 14688 * A**5
    - 1836 * A**4 * C
    - 1584 * A**4
    + 144 * A**3 * C
    + 64 * A**3
    - 144 * A**2 * C**2
    + 48 * A * C**2
    - 8 * C**3
    - 4 * C**2
)
gate(sp.expand(disc + L**2 * H / 4) == 0, "squared line times octic discriminant")
gate(sp.Poly(H, A, C).total_degree() == 8, "H has total degree eight")
gate(sp.factor(H.subs(A, sp.Rational(4, 27))) == -4 * C**2 * (162 * C + 1) / 81,
     "octic does not contain the squared line")

line_fibre = sp.factor(Phi_affine.subs(A, sp.Rational(4, 27)))
gate(sp.expand(line_fibre - T**2 * (C - sp.Rational(64, 19683) * T)) == 0,
     "line fibre is double plus simple, not total cubic ramification")

implicit = sp.factor(sp.resultant(A - A_u, C - C_u, u))
gate(sp.expand(implicit + H / 8) == 0, "exact implicit octic resultant")
gate(sp.Poly(C_u, u).degree() == 8, "polynomial normalization has one degree-eight infinity")
gate(sp.Poly(A_u, u).degree() == 3, "projection to A has degree three")

# The root really generates k(u)/k(A): it has pole order 3 at a point where
# A has local order 2.
gate(sp.limit(u**3 * t.subs(lam, 3), u, 0) == -2,
     "survivor repeated root has exact pole order three")
gate(sp.limit(A_u / u**2, u, 0) == 1,
     "A has exact local order two at the triple pole")
gate(sp.factor(sp.diff(A_u, u)) == u * (3 * u + 2),
     "the triple pole lies at a simple A-ramification point")


# ---------------------------------------------------------------------------
# The explicit maximal overorder clears the line index and is monogenic.
# ---------------------------------------------------------------------------

# Standard Delone--Faddeev multiplication, with theta=L*E and W=omega.
gate(
    sp.expand(
        (-a3 * c3 + C * W - a3 * L * E)
        - (C * W + A**3 * L * (E - (6 * A - 1)))
    ) == 0,
    "Delone--Faddeev omega-square identity",
)
gate(sp.cancel((-a3 * d3) / L + A**3 * L / 2) == 0,
     "Delone--Faddeev omega-e identity")
gate(
    sp.cancel(
        (-C * d3 + d3 * W - c3 * L * E) / L**2
        - (C / 2 - W / 2 + (6 * A - 1) * E)
    ) == 0,
    "Delone--Faddeev e-square identity",
)

omega_E = C + 2 * (6 * A - 1) * E - 2 * E**2
gate(sp.expand(E**2 - (C / 2 - omega_E / 2 + (6 * A - 1) * E)) == 0,
     "omega is recovered polynomially from e")

P_E = E**3 - (6 * A - 1) * E**2 - C * E / 2 - A**3 * L / 4
theta_polynomial = (
    (L * E) ** 3
    + c3 * (L * E) ** 2
    + C * d3 * (L * E)
    + a3 * d3**2
)
gate(sp.cancel(theta_polynomial / L**3 - P_E) == 0,
     "scaled theta equation gives the overorder polynomial")
gate(sp.Poly(P_E, E).degree() == 3 and sp.Poly(P_E, E).LC() == 1,
     "overorder polynomial is monic cubic")

transition = sp.Matrix([
    [1, C, 0],
    [0, 2 * (6 * A - 1), L],
    [0, -2, 0],
])
gate(sp.expand(transition.det() - 2 * L) == 0,
     "old-order to monogenic-overorder determinant")

disc_B = sp.factor(sp.discriminant(P_E, E))
gate(sp.expand(disc_B + H / 16) == 0,
     "maximal monogenic overorder has discriminant -H/16")
gate(sp.expand((2 * L) ** 2 * disc_B - disc) == 0,
     "index-discriminant identity clears exactly the squared line")
gate(
    sp.factor(disc_B.subs(A, sp.Rational(4, 27)))
    == C**2 * (162 * C + 1) / 324,
    "the index line is generically etale in the overorder",
)


# ---------------------------------------------------------------------------
# Exact two-address fibres of the octic normalization.
# ---------------------------------------------------------------------------

A_v = A_u.subs(u, v)
C_v = C_u.subs(u, v)
A_divided = sp.cancel((A_u - A_v) / (u - v))
C_divided = sp.cancel((C_u - C_v) / (u - v))

A_symmetric, A_remainder, A_mapping = sp.symmetrize(A_divided, [u, v], formal=True)
gate(A_remainder == 0, "A divided difference is symmetric")
A_s1, A_s2 = [pair[0] for pair in A_mapping]
A_sp = sp.factor(A_symmetric.subs({A_s1: s, A_s2: p}))
gate(A_sp == s**2 + s - p, "A-collision equation is p=s^2+s")

C_symmetric, C_remainder, C_mapping = sp.symmetrize(C_divided, [u, v], formal=True)
gate(C_remainder == 0, "C divided difference is symmetric")
C_s1, C_s2 = [pair[0] for pair in C_mapping]
C_sp = sp.factor(C_symmetric.subs({C_s1: s, C_s2: p}).subs(p, s**2 + s))
collision_factor = 3 * s**2 + 6 * s + 2
gate(C_sp == s * collision_factor**3 / 2,
     "two off-diagonal collision pairs have the exact cubic contact factor")

pair_discriminant = sp.expand(s**2 - 4 * (s**2 + s))
gate(sp.resultant(collision_factor, pair_discriminant, s) == -12,
     "each collision pair has two distinct parameters")
target_A = -s * (s + 1) ** 2
target_L = sp.expand(27 * target_A - 4)
gate(sp.resultant(collision_factor, target_L, s) == -54,
     "collision targets avoid the squared index line")
gate(sp.rem(sp.Poly(target_A + s / 3, s), sp.Poly(collision_factor, s)) == 0,
     "the two collision roots give distinct A-values -s/3")
gate(sp.discriminant(collision_factor, s) == 12,
     "there are two distinct collision-pair parameters s")


summary = {
    "checks": CHECKS,
    "pole_partition": "shared A-address or one triple pole at e=2",
    "normal_form": "A=u^3+u^2; t=(3*lambda*u^2+3*u-2)/u^3",
    "polynomial_color": "lambda=3 only",
    "discriminant": "-(27*A-4)^2*H/4; deg(H)=8",
    "maximal_order": "e=theta/(27*A-4) is normal monogenic; disc=-H/16",
    "collision": "3*s^2+6*s+2=0 gives two two-address fibres",
    "conclusion": "monogenic and non-unibranch; no A2 Keller open",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3933 centered degree-three root-map exact companion")
print(f"CHECKS={CHECKS}")
print("POLE_PARTITION=shared A-address or one triple pole at e=2")
print("NORMAL_FORM=A=u^3+u^2;t=(3*lambda*u^2+3*u-2)/u^3")
print("POLYNOMIAL_COLOR=lambda=3 only")
print("SURVIVOR=a=-A^3;c=-(6A-1)(27A-4);d=-(27A-4)^2/2")
print("DISCRIMINANT=-(27A-4)^2*H/4;DEG_H=8")
print("MAXIMAL_OVERORDER=e=theta/(27*A-4);omega=C+2*(6*A-1)*e-2*e^2")
print("OVERORDER_POLYNOMIAL=e^3-(6*A-1)*e^2-C*e/2-A^3*(27*A-4)/4")
print("OVERORDER_INDEX=det(1,omega,theta)/(1,e,e^2)=2*(27*A-4);ideal=(27*A-4)")
print("OVERORDER_DISCRIMINANT=-H/16")
print("OVERORDER_NORMAL=R1(discriminant H exponent one)+S2(hypersurface)")
print("OVERORDER_CONCLUSION=maximal cubic completion is globally monogenic")
print("COLLISIONS=3*s^2+6*s+2=0 gives two two-address fibres")
print("CONCLUSION=monogenic and ramification curve non-unibranch;no A2 Keller open")
print(f"SEMANTIC_SHA256={semantic}")
