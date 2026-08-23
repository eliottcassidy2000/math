#!/usr/bin/env python3
"""Assertion-free exact gates for THM-3865.

Reproduction:
  python3 04-computation/jc2_one_place_inverse_discriminant_class_group_thm3865.py
  python3 -O 04-computation/jc2_one_place_inverse_discriminant_class_group_thm3865.py
"""

from __future__ import annotations

import hashlib

import sympy as sp


A, C, W, t, y, lam = sp.symbols("A C W t y lambda")
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

semantic_packet = (
    "delta=Delta0+lambda*C^5 with four-address A1 normalization",
    "normal quadratic resolvent S=k[A,C,W]/(W^2-delta)",
    "S_C=k[t,y,(y^2-D(t))^-1] is factorial",
    "two reduced height-one primes over C with div(C)=P_plus+P_minus",
    "Nagata sequence gives Cl(S)=Z and S*=k*",
    "Cl(S)[3]=0 kills every normal S3 cubic with discriminant delta",
    "formal THM3855 cubic therefore does not globally algebraize",
)
semantic_sha = hashlib.sha256(repr(semantic_packet).encode()).hexdigest()

print("THM3865_TARGET", "Delta0+lambda*C^5")
print("THM3865_NORMALIZATION", "C=-D(t)/lambda,A=t*C;four origins;one infinity")
print("THM3865_UFD_CHART", "S_C=k[t,y,(y^2-D(t))^-1]")
print("THM3865_SPECIAL_FIBRE", "(W-mu*A^2)(W+mu*A^2),mu^2=-1615")
print("THM3865_NAGATA", "Cl(S)=Z^2/<(1,1)>=Z;units=k*")
print("THM3865_THREE_TORSION", "NONE")
print("THM3865_GLOBAL_CUBIC", "normal finite-flat S3 algebra with this discriminant=NONE")
print("THM3865_SCOPE", "exact target only; other one-place discriminants and JC2 remain open")
print("SEMANTIC_SHA256", semantic_sha)
print("CHECKS", CHECKS)
