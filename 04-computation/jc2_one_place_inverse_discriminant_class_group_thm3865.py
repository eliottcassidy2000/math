#!/usr/bin/env python3
"""Assertion-free exact gates for THM-3865.

Reproduction:
  python3 04-computation/jc2_one_place_inverse_discriminant_class_group_thm3865.py
  python3 -O 04-computation/jc2_one_place_inverse_discriminant_class_group_thm3865.py
"""

from __future__ import annotations

import hashlib

import sympy as sp


A, C, W, t, y, ell, lam = sp.symbols("A C W t y ell lambda")
CHECKS = 0


def zero(label: str, value: sp.Expr) -> None:
    global CHECKS
    CHECKS += 1
    value = sp.factor(value)
    if value != 0:
        raise AssertionError(f"{label}: {value}")


def nonzero(label: str, value: sp.Expr) -> None:
    global CHECKS
    CHECKS += 1
    value = sp.factor(value)
    if value == 0:
        raise AssertionError(f"{label}: unexpectedly zero")


def equal(label: str, left: object, right: object) -> None:
    global CHECKS
    CHECKS += 1
    if left != right:
        raise AssertionError(f"{label}: {left!r} != {right!r}")


Delta0 = sp.factor(A * (C + 5 * A) * (4 * C + 19 * A) * (3 * C - 17 * A))
D = sp.factor(t * (5 * t + 1) * (19 * t + 4) * (3 - 17 * t))
delta = sp.expand(Delta0 + lam * C**5)

zero("radial_discriminant", Delta0.subs({A: t * C}) - C**4 * D)
equal("D_degree", sp.degree(D, t), 4)
nonzero("D_squarefree", sp.discriminant(D, t))
zero("special_fibre", delta.subs(C, 0) + 1615 * A**4)

# Four-address affine normalization and its immersive derivative control.
Cpar = -D / lam
Apar = t * Cpar
zero("branch_parametrization", delta.subs({A: Apar, C: Cpar}))
equal("C_parameter_degree", sp.degree(Cpar, t), 4)
equal("A_parameter_degree", sp.degree(Apar, t), 5)
equal("normalization_derivative_gcd", sp.gcd(D, sp.diff(D, t)), 1)

# The quadratic resolvent chart after inverting C.
F = sp.expand(y**2 - D)
Cback = F / lam
Aback = t * Cback
Wback = y * Cback**2
surface = sp.expand(W**2 - delta)
zero("factorial_chart_surface", surface.subs({A: Aback, C: Cback, W: Wback}))
zero("factorial_chart_C_inverse", Cback - F / lam)
zero("factorial_chart_A_inverse", Aback / Cback - t)
zero("factorial_chart_W_inverse", Wback / Cback**2 - y)

# Irreducibility control for F=y^2-D: D has four odd simple valuations.
equal("F_y_degree", sp.degree(F, y), 2)
equal("F_y_linear_coefficient", sp.Poly(F, y).coeff_monomial(y), 0)
equal("F_y_quadratic_coefficient", sp.Poly(F, y).coeff_monomial(y**2), 1)
nonzero("F_constant_not_square_control", sp.discriminant(D, t))

# Delone--Faddeev packet and the two reduced primes over C=0.
a0, b0, c0, d0 = A, C, 7 * A, -3 * A
binary_disc = sp.expand(
    b0**2 * c0**2
    - 4 * a0 * c0**3
    - 4 * b0**3 * d0
    - 27 * a0**2 * d0**2
    + 18 * a0 * b0 * c0 * d0
)
zero("binary_cubic_control", binary_disc - Delta0)
mu = sp.symbols("mu")
special_surface = sp.expand(surface.subs(C, 0))
zero("special_surface_equation", special_surface - (W**2 + 1615 * A**4))
zero(
    "two_prime_special_fibre",
    sp.expand(special_surface - (W - mu * A**2) * (W + mu * A**2)).subs(
        mu**2, -1615
    ),
)
nonzero("two_prime_separation_characteristic_zero", 2 * A**2)

# Independent non-tangent orientation L=A+C, M=A.  This freezes the
# coordinate-free repetition used in the all-line strengthening.
D_sum = sp.factor(
    Delta0.subs({A: t * ell, C: (1 - t) * ell}, simultaneous=True) / ell**4
)
equal(
    "sum_orientation_D",
    D_sum,
    -t * (4 * t + 1) * (15 * t + 4) * (20 * t - 3),
)
equal("sum_orientation_D_degree", sp.degree(D_sum, t), 4)
nonzero("sum_orientation_D_squarefree", sp.discriminant(D_sum, t))

delta_sum = sp.expand(Delta0 + lam * (A + C) ** 5)
Lpar_sum = -D_sum / lam
Apar_sum = t * Lpar_sum
Cpar_sum = (1 - t) * Lpar_sum
zero(
    "sum_orientation_parametrization",
    delta_sum.subs({A: Apar_sum, C: Cpar_sum}, simultaneous=True),
)
equal("sum_orientation_L_degree", sp.degree(Lpar_sum, t), 4)
equal("sum_orientation_M_degree", sp.degree(Apar_sum, t), 5)
equal("sum_orientation_derivative_gcd", sp.gcd(D_sum, sp.diff(D_sum, t)), 1)

F_sum = sp.expand(y**2 - D_sum)
Lback_sum = F_sum / lam
Aback_sum = t * Lback_sum
Cback_sum = (1 - t) * Lback_sum
Wback_sum = y * Lback_sum**2
surface_sum = sp.expand(W**2 - delta_sum)
zero(
    "sum_orientation_factorial_chart",
    surface_sum.subs(
        {A: Aback_sum, C: Cback_sum, W: Wback_sum}, simultaneous=True
    ),
)
zero("sum_orientation_special_fibre", delta_sum.subs(C, -A) + 1200 * A**4)
special_surface_sum = sp.expand(surface_sum.subs(C, -A))
zero(
    "sum_orientation_two_primes",
    sp.expand(
        special_surface_sum - (W - mu * A**2) * (W + mu * A**2)
    ).subs(mu**2, -1200),
)

semantic_packet = (
    "delta_L=Delta0+lambda*L^5 for every non-tangent linear L",
    "four-address A1 normalization in complementary coordinates (M,L)",
    "normal quadratic resolvent S_L=k[A,C,W]/(W^2-delta_L)",
    "(S_L)_L=k[t,y,(y^2-D_L(t))^-1] is factorial",
    "two reduced height-one primes over L with div(L)=P_plus+P_minus",
    "Nagata sequence gives Cl(S_L)=Z and S_L*=k*",
    "Cl(S_L)[3]=0 kills every normal S3 cubic with discriminant delta_L",
    "independent controls L=C and L=A+C",
    "formal THM3855 cubic therefore does not globally algebraize",
)
semantic_sha = hashlib.sha256(repr(semantic_packet).encode()).hexdigest()

print("THM3865_TARGET", "Delta0+lambda*L^5;every non-tangent linear L")
print("THM3865_CONTROLS", "L=C and L=A+C")
print("THM3865_NORMALIZATION", "L=-D_L(t)/lambda,M=t*L;four origins;one infinity")
print("THM3865_UFD_CHART", "(S_L)_L=k[t,y,(y^2-D_L(t))^-1]")
print("THM3865_SPECIAL_FIBRE", "(W-mu*M^2)(W+mu*M^2),mu^2=kappa")
print("THM3865_NAGATA", "Cl(S_L)=Z^2/<(1,1)>=Z;units=k*")
print("THM3865_THREE_TORSION", "NONE")
print("THM3865_GLOBAL_CUBIC", "normal finite-flat S3 algebra with this discriminant=NONE")
print("THM3865_SCOPE", "non-tangent fifth-degree family; other one-place discriminants and JC2 open")
print("SEMANTIC_SHA256", semantic_sha)
print("CHECKS", CHECKS)
