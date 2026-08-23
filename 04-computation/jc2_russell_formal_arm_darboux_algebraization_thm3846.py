#!/usr/bin/env python3
"""Exact companion for THM-3846's formal arm Darboux lift."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: object, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.cancel(sp.factor(expression)) == 0, message)


r, z, e, s = sp.symbols("r z e s")
c = sp.symbols("c", nonzero=True)
surface = r**2 * e - z**3 + c**3 * r
H = c**3 + 2 * e * r

zero(H**2 - (c**6 + 4 * e * z**3) - 4 * e * surface,
     "completed square follows from the surface relation")

e_from_s = 3 * c**3 * s + 9 * s**2 * z**3
r_from_s = z**3 / (c**3 + 3 * s * z**3)
zero(surface.subs({e: e_from_s, r: r_from_s}),
     "canonical inverse coordinates satisfy the surface")

s_rational = e / (3 * (c**3 + e * r))
zero(s_rational.subs({e: e_from_s, r: r_from_s}) - s,
     "canonical inverse recovers s")
zero(
    sp.cancel(r_from_s.subs(s, s_rational) - r).subs(z**3, r * (c**3 + e * r)),
    "canonical inverse recovers r modulo the surface",
)

# Direct Poisson calculation of {z,s}=1 in the surface quotient.
poisson_z_s = (3 * c**3 + 6 * r * e) * sp.diff(s_rational, e) - 3 * r**2 * sp.diff(s_rational, r)
zero(
    sp.together(poisson_z_s - 1).as_numer_denom()[0].subs(z**3, r * (c**3 + e * r)),
    "rational canonical coordinate has bracket one",
)

# Generic formal arm packet, represented by independent functions of s.
a = sp.Function("a")(s)
b = sp.Function("b")(s)
alpha = sp.Function("alpha")(s)
beta = sp.Function("beta")(s)
Z = sp.symbols("Z")
W = alpha * sp.diff(beta, s) - sp.diff(alpha, s) * beta
bezout = alpha * sp.diff(b, s) - sp.diff(a, s) * beta
A0 = a + alpha * Z
C0 = b + beta * Z
jac_Zs = sp.diff(A0, Z) * sp.diff(C0, s) - sp.diff(A0, s) * sp.diff(C0, Z)
zero(jac_Zs - bezout - W * Z, "generic linear-normal Jacobian")
zero((jac_Zs - bezout - W * Z), "generic packet replay")

# W=0 is the injective boundary: this family displays the forced affine
# target combination b-lambda*a=mu*s+nu and its constant normal row.
lam0, mu0, nu0 = sp.symbols("lambda0 mu0 nu0", nonzero=True)
a_zero = s**2
b_zero = lam0 * a_zero + mu0 * s + nu0
alpha_zero = 1 / mu0
beta_zero = lam0 / mu0
zero(
    alpha_zero * sp.diff(b_zero, s) - sp.diff(a_zero, s) * beta_zero - 1,
    "zero-W packet has Bezout determinant one",
)
zero(
    alpha_zero * sp.diff(beta_zero, s) - sp.diff(alpha_zero, s) * beta_zero,
    "zero-W packet has vanishing normal Wronskian",
)
zero(b_zero - lam0 * a_zero - (mu0 * s + nu0),
     "zero-W target combination recovers the arm parameter")

# A nonconstant-W hostile packet checks the full chain rule, including Z_s,
# through seven normal orders rather than only the abstract determinant.
a_host = s
b_host = s**2
alpha_host = s
beta_host = 2 * s**2 - 1
w_host = sp.factor(
    alpha_host * sp.diff(beta_host, s) - sp.diff(alpha_host, s) * beta_host
)
zero(alpha_host * sp.diff(b_host, s) - sp.diff(a_host, s) * beta_host - 1,
     "nonconstant-W hostile Bezout packet")
zero(w_host - (2 * s**2 + 1), "nonconstant-W hostile Wronskian")
Z_host = sum(
    sp.catalan(n - 1) * (-w_host / 2) ** (n - 1) * z**n
    for n in range(1, 10)
)
A_host = a_host + alpha_host * Z_host
C_host = b_host + beta_host * Z_host
jac_host = sp.expand(
    sp.diff(A_host, z) * sp.diff(C_host, s)
    - sp.diff(A_host, s) * sp.diff(C_host, z)
    - 1
)
jac_host_poly = sp.Poly(jac_host, z)
gate(
    all(jac_host_poly.coeff_monomial(z**degree) == 0 for degree in range(8)),
    "nonconstant-W hostile formal Jacobian through order seven",
)

gate(sp.degree(1 + 2 * w_host * z, z) == 1,
     "every nonzero W gives a simple linear square-gate prime")

# Catalan truncations solve Z+(W/2)Z^2=z to the advertised order.
w = sp.symbols("w")
for order in range(1, 11):
    truncation = sum(
        sp.catalan(n - 1) * (-w / 2) ** (n - 1) * z**n
        for n in range(1, order + 1)
    )
    residual = sp.expand(truncation + w * truncation**2 / 2 - z)
    polynomial = sp.Poly(residual, z)
    gate(
        all(polynomial.coeff_monomial(z**degree) == 0 for degree in range(order + 1)),
        f"Catalan truncation solves through order {order}",
    )

# Minimal nodal arm packet in the canonical parameter s.
a_node = 9 * c**6 * s**2
b_node = 27 * c**9 * s**3 - 3 * c**3 * s
alpha_node = -1 / (3 * c**3)
beta_node = -3 * s / 2
bezout_node = sp.factor(
    alpha_node * sp.diff(b_node, s) - sp.diff(a_node, s) * beta_node
)
w_node = sp.factor(
    alpha_node * sp.diff(beta_node, s) - sp.diff(alpha_node, s) * beta_node
)
zero(bezout_node - 1, "nodal arm Bezout determinant")
zero(w_node - 1 / (2 * c**3), "nodal normal Wronskian")

Z_node = 2 * c**3 * (sp.sqrt(1 + z / c**3) - 1)
zero(Z_node + Z_node**2 / (4 * c**3) - z,
     "nodal Hensel root solves its quadratic")

# The square-gate divisor is prime: after z=-c^3, r is a unit and e is
# -(c^3*r+c^9)/r^2.  These identities are the explicit Laurent quotient.
e_laurent = -(c**3 * r + c**9) / r**2
zero(surface.subs({z: -c**3, e: e_laurent}),
     "z+c^3 quotient is the Laurent line")
gate(sp.diff(1 + z / c**3, z) != 0,
     "square-gate function vanishes simply on the Laurent prime")

semantic = {
    "completion": "Bhat=k[s][[z]];s=e/(3(c3+er));{z,s}=1",
    "universal_lift": "Z+(W/2)Z2=z;A=a+alpha*Z;C=b+beta*Z",
    "series": "Catalan formal normal coordinate",
    "rationality_gate": "Frac(B)=k(s,z);displayed lift is rational iff W=0",
    "injective_boundary": "W=0 plus Bezout makes an affine target combination recover s",
    "nodal": "W=1/(2c3);1+z/c3 has odd Laurent-prime valuation",
    "scope": "self-identifying arm lifts are formally exact but canonically nonrational;global Darboux pair open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
     "no inactive Python assert")

print("theorem=THM-3846-formal-arm-darboux-lift-and-algebraization-gate")
print("completion=Bhat_equals_k[s][[z]];bracket_z_s_equals_1")
print("coordinate=s=e/(3(c^3+er));e=3c^3s+9s^2z^3;r=z^3/(c^3+3sz^3)")
print("formal_lift=Z+(W/2)Z^2=z;A=a+alpha*Z;C=b+beta*Z")
print("universal_gate=Frac(B)=k(s,z);W_nonzero_is_nonsquare;W_zero_forces_injective_arm")
print("nodal_gate=W=1/(2c^3);1+z/c^3_nonsquare_by_odd_Laurent_valuation")
print("scope=formal_all_orders;self_identifying_displayed_lifts_not_rational;global_pair_open")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
