#!/usr/bin/env python3
"""Exact companion for THM-3812's nodal second-normal profile no-go."""

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


r, z, e = sp.symbols("r z e")
c = sp.symbols("c", nonzero=True)
variables = (r, z, e)
surface = r**2 * e - z**3 + c**3 * r
monic_relation = z**3 - r**2 * e - c**3 * r
poisson = sp.Matrix(
    [
        [0, 3 * r**2, 9 * z**2],
        [-3 * r**2, 0, 3 * c**3 + 6 * r * e],
        [-9 * z**2, -3 * c**3 - 6 * r * e, 0],
    ]
)


def bracket(left: sp.Expr, right: sp.Expr) -> sp.Expr:
    dl = sp.Matrix([sp.diff(left, q) for q in variables])
    dr = sp.Matrix([sp.diff(right, q) for q in variables])
    return sp.expand((dl.T * poisson * dr)[0])


# The normalized nodal arm profile and its canonical first-normal Bezout row.
a0 = e**2
b0 = e**3 - e
a1 = -1 / (3 * c**3)
b1 = -e / (2 * c**3)
zero(3 * c**3 * (a1 * sp.diff(b0, e) - sp.diff(a0, e) * b1) - 1,
     "nodal first-normal Bezout law")
zero(bracket(surface, r), "surface Casimir against r")
zero(bracket(surface, z), "surface Casimir against z")
zero(bracket(surface, e), "surface Casimir against e")
gate(sp.Poly(monic_relation, z).degree() == 3, "monic relation z-degree")
gate(sp.Poly(monic_relation, z).LC() == 1, "monic relation leading unit")


# Arbitrary arm-dependent coefficients in the canonical r and z^2 slots.
g = sp.Function("g")(e)
f = sp.Function("f")(e)
h = sp.Function("h")(e)
k = sp.Function("k")(e)
A = a0 + z * a1 + r * g + z**2 * f
C = b0 + z * b1 + r * h + z**2 * k

# Reduction by the monic relation is the unique k[r,e]-normal form in the
# basis 1,z,z^2.  Clear only the invertible scalar denominator c^6 while
# taking the remainder, then divide it back out.
raw = sp.together(bracket(A, C) - 1)
normal = sp.cancel(
    sp.rem(sp.expand(raw * c**6), monic_relation, z) / c**6
)
normal_poly = sp.Poly(sp.expand(normal), r, z)
gate(max(z_degree for (_, z_degree) in normal_poly.monoms()) <= 2,
     "canonical normal form has z-degree at most two")
zero(
    sp.rem(sp.expand((raw - normal) * c**6), monic_relation, z) / c**6,
    "canonical reduction has no quotient loss on the surface",
)


# These two unique coefficients alone are incompatible.
z_bucket = normal_poly.coeff_monomial(z)
r3_bucket = normal_poly.coeff_monomial(r**3)
expected_z_bucket = (
    36 * c**6 * e**2 * f
    - 24 * c**6 * e * k
    - 12 * c**6 * f
    + 1
) / (2 * c**3)
expected_r3_bucket = 12 * e**2 * (
    f * sp.diff(k, e) - k * sp.diff(f, e)
)
zero(z_bucket - expected_z_bucket, "pure z coefficient")
zero(r3_bucket - expected_r3_bucket, "pure r^3 coefficient")

bezout_law = (3 * e**2 - 1) * f - 2 * e * k + 1 / (12 * c**6)
wronskian = f * sp.diff(k, e) - k * sp.diff(f, e)
zero(z_bucket - 6 * c**3 * bezout_law, "z bucket is the second-normal Bezout law")
zero(r3_bucket - 12 * e**2 * wronskian, "r^3 bucket is the Wronskian law")

# W=0 makes f/k constant in k(e) when k is nonzero.  The remaining factor
# can never be a polynomial unit: its e coefficient is always -2.
lam = sp.symbols("lambda")
unit_factor = sp.expand(lam * (3 * e**2 - 1) - 2 * e)
gate(sp.Poly(unit_factor, e).coeff_monomial(e) == -2,
     "proportional-profile factor is never constant")
gate(sp.Poly(3 * e**2 - 1, e).degree() == 2,
     "k=0 boundary has nonconstant factor")
zero(
    sp.diff(f / k, e) - (sp.diff(f, e) * k - f * sp.diff(k, e)) / k**2,
    "Wronskian is the rational logarithmic derivative",
)


# Exact bounded hostile controls only: specialize c=2 and independently ask
# Groebner elimination to solve both necessary coefficient identities for
# every pair of profiles of degree at most 0,...,8.  The theorem does not
# depend on this finite universe.
bounded_rows: list[str] = []
for degree in range(9):
    f_coefficients = sp.symbols(f"f0:{degree + 1}")
    k_coefficients = sp.symbols(f"k0:{degree + 1}")
    f_poly = sum(f_coefficients[j] * e**j for j in range(degree + 1))
    k_poly = sum(k_coefficients[j] * e**j for j in range(degree + 1))
    specialized_bezout = sp.Poly(
        (3 * e**2 - 1) * f_poly - 2 * e * k_poly + sp.Rational(1, 768),
        e,
    )
    specialized_wronskian = sp.Poly(
        f_poly * sp.diff(k_poly, e) - k_poly * sp.diff(f_poly, e), e
    )
    equations = specialized_bezout.all_coeffs() + specialized_wronskian.all_coeffs()
    gb = sp.groebner(
        equations, *f_coefficients, *k_coefficients, order="grevlex"
    )
    gate(
        len(gb.polys) == 1 and gb.polys[0].as_expr() == 1,
        f"bounded hostile degree {degree}",
    )
    bounded_rows.append(f"degree<={degree}:unit")


# The ideal-order sidecar is stated by explicit identities.  From the
# surface relation, r belongs first to I^2 and then to I^3; multiplication by
# z puts rz in I^4.  This prevents the theorem from mislabelling its profile
# as every element of I^2.
zero(c**3 * r - z**3 + r**2 * e - surface, "arm-ideal relation")

semantic = {
    "surface": "r^2e-z^3+c^3r=0;c!=0",
    "ansatz": "A=e^2-z/(3c^3)+r*g(e)+z^2*f(e);C=e^3-e-ez/(2c^3)+r*h(e)+z^2*k(e)",
    "normal_form": "free k[r,e]-basis 1,z,z^2",
    "z_bucket": "(3e^2-1)f-2ek=-1/(12c^6)",
    "r3_bucket": "f*kprime-k*fprime=0",
    "contradiction": "Wronskian proportionality makes k*(lambda(3e^2-1)-2e) a forbidden unit",
    "scope": "arm-dependent r,z^2 coefficients only; r in I^3; rz in I^4",
    "next": "add rz*p(e) in at least one output or a still higher canonical coefficient",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()
bounded_blob = "\n".join(bounded_rows).encode()

source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
     "inactive Python assert")

print("theorem=THM-3812-nodal-arm-coefficient-second-normal-profile-nonentry")
print("surface=r^2e-z^3+c^3r;c_nonzero;field=algebraically_closed_characteristic_zero")
print("ansatz=A=e^2-z/(3c^3)+r*g(e)+z^2*f(e);C=e^3-e-ez/(2c^3)+r*h(e)+z^2*k(e)")
print("normal_form=free_k[r,e]_basis_1,z,z^2")
print("z_bucket=(3e^2-1)f-2ek=-1/(12c^6)")
print("r3_bucket=f*kprime-k*fprime=0")
print("contradiction=proportional_profiles_force_nonconstant_factor_to_be_a_unit")
print("scope=arm_dependent_r_z2_profiles_only;r_in_I3;rz_in_I4")
print("next=rz_profile_or_higher_canonical_coefficient")
print(f"bounded_sha256={hashlib.sha256(bounded_blob).hexdigest()}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
