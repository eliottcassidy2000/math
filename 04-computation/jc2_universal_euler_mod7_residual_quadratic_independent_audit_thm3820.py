#!/usr/bin/env python3
"""Cheapest exact test of the RESERVED THM-3820 Euler/mod-seven idea."""

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def check(condition, label):
    global CHECKS
    if not bool(condition):
        raise RuntimeError(f"FAILED: {label}")
    CHECKS += 1


def same(left, right, label):
    check(sp.expand(left - right) == 0, label)


u, t, Y, Z, E, lam = sp.symbols("u t Y Z E lam")

# After G=e^3 Y and D=e^2(Y+Z), the universal source equation P/e^3 is
# minus this quadratic.  Its sign does not affect its discriminant.
p = (Y + 2 * Z) * u**2 + (Y + Z + 4) * u + 2
q = (1 + 2 * u)**3 - 729 * t * Y**3 * u**2 * (1 + u)**2

residual_with_content = sp.factor(sp.resultant(p, q, u))
check(sp.rem(sp.Poly(residual_with_content, Y), sp.Poly(Y**3, Y)) == 0,
      "universal resultant has Y^3 content")
residual = sp.cancel(residual_with_content / Y**3)
check(sp.denom(residual) == 1, "content-free residual is polynomial")
residual = sp.factor(residual)
check(sp.Poly(residual, t).degree() == 2, "residual is quadratic in t")

source_discriminant = sp.factor(sp.discriminant(p, u))
expected_source_discriminant = (Y + Z)**2 - 8 * Z + 16
same(source_discriminant, expected_source_discriminant,
     "source critical quadratic discriminant")

residual_discriminant = sp.factor(sp.discriminant(residual, t))
quotient = sp.factor(sp.cancel(residual_discriminant / source_discriminant))
check(sp.denom(quotient) == 1, "source discriminant divides residual discriminant")

discriminant_scalar, factors = sp.factor_list(quotient)
check(all(exponent % 2 == 0 for _, exponent in factors),
      "discriminant quotient is a square up to scalar")
square_root = sp.Integer(1)
for factor, exponent in factors:
    square_root *= factor ** (exponent // 2)
same(quotient, discriminant_scalar * square_root**2,
     "explicit discriminant square identity")

# The extra square is a genuine seam, not cosmetic.  At this rational jet the
# source quadratic is separable while the residual t-quadratic has a double root.
hostile_jet = {Y: -sp.Rational(9, 5), Z: sp.Integer(1)}
hostile_source_discriminant = sp.factor(source_discriminant.subs(hostile_jet))
hostile_square_factor = sp.factor(square_root.subs(hostile_jet))
hostile_residual = sp.Poly(residual.subs(hostile_jet), t)
check(hostile_source_discriminant == sp.Rational(216, 25),
      "hostile jet source quadratic remains separable")
check(hostile_square_factor == 0, "hostile jet lies on extra square divisor")
check(hostile_residual.degree() == 2 and hostile_residual.LC() != 0,
      "hostile residual remains genuinely quadratic")
check(sp.discriminant(hostile_residual.as_expr(), t) == 0,
      "hostile residual t-quadratic has a double root")
hostile_double_root = sp.factor(
    -hostile_residual.coeff_monomial(t) / (2 * hostile_residual.LC()))

# First genuinely degree-six control: g=1+lam*e^6.  Pull the universal
# identity back and test the repeated-root-safe boundary divisibility directly.
g6 = 1 + lam * E**6
Y6 = g6 / E**3
Z6 = (E * sp.diff(g6, E) - g6) / E**3
H6 = sp.factor(sp.cancel(E**3 * residual.subs({t: E**7, Y: Y6, Z: Z6})))
check(sp.denom(H6) == 1, "sextic-binomial pullback is polynomial")
H6_poly = sp.Poly(H6, E)
check(H6_poly.degree() > 0, "sextic-binomial residual is nonconstant")
check(H6_poly.LC() != 0, "sextic-binomial residual leading coefficient")
boundary_remainder = sp.Poly(
    sp.expand(E * g6 * sp.diff(H6, E)), E,
    domain=sp.QQ.frac_field(lam),
).rem(sp.Poly(H6, E, domain=sp.QQ.frac_field(lam)))
check(not boundary_remainder.is_zero,
      "generic sextic binomial is not boundary-supported")
cleared_boundary_numerator, cleared_boundary_denominator = sp.fraction(
    sp.together(boundary_remainder.as_expr()))
cleared_boundary_poly = sp.Poly(cleared_boundary_numerator, E)
coefficient_gcd_poly = sp.Poly(0, lam)
for boundary_coefficient in cleared_boundary_poly.all_coeffs():
    coefficient_gcd_poly = sp.gcd(
        coefficient_gcd_poly, sp.Poly(boundary_coefficient, lam))
coefficient_gcd = sp.factor(coefficient_gcd_poly.as_expr())
check(sp.Poly(coefficient_gcd, lam).length() == 1,
      "sextic boundary coefficient gcd is a monomial")
check(sp.Poly(coefficient_gcd, lam).terms()[0][0][0] >= 0,
      "sextic boundary gcd has only lam support")

source = Path(__file__).read_text(encoding="utf-8")
check(not any(isinstance(node, ast.Assert)
              for node in ast.walk(ast.parse(source))),
      "assertion-free source")

semantic = {
    "p": sp.sstr(p),
    "q": sp.sstr(q),
    "residual_degree_t": 2,
    "source_discriminant": sp.sstr(source_discriminant),
    "discriminant_scalar": sp.sstr(discriminant_scalar),
    "discriminant_square_root": sp.sstr(sp.factor(square_root)),
    "hostile_extra_square": "Y=-9/5,Z=1: source disc=216/25 but residual disc=0",
    "sextic_control": "g=1+lam*e^6 boundary coefficient gcd is lam-power",
}

print("target=RESERVED-THM-3820-cheapest-exact-test")
print(f"source_discriminant={sp.sstr(source_discriminant)}")
print(f"residual={sp.sstr(residual)}")
print(f"residual_discriminant={sp.sstr(residual_discriminant)}")
print(f"quotient_scalar={sp.sstr(discriminant_scalar)}")
print(f"quotient_square_root={sp.sstr(sp.factor(square_root))}")
print("hostile_jet=Y:-9/5,Z:1;source_disc=216/25;residual_disc=0")
print(f"hostile_double_t_root={sp.sstr(hostile_double_root)}")
print(f"sextic_binomial_H={sp.sstr(H6)}")
print(f"sextic_boundary_remainder_degree={boundary_remainder.degree()}")
print(f"sextic_boundary_denominator={sp.factor(cleared_boundary_denominator)}")
print(f"sextic_boundary_coefficient_gcd={coefficient_gcd}")
print(f"CHECKS={CHECKS}")
print("semantic_sha256=" + hashlib.sha256(json.dumps(
    semantic, sort_keys=True, separators=(",", ":")).encode()).hexdigest())
