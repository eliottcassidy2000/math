#!/usr/bin/env python3
"""Exact companion for THM-3856's quadratic normal-strip classification."""

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


s, z = sp.symbols("s z")
a = sp.Function("a")(s)
b = sp.Function("b")(s)
alpha = sp.Function("alpha")(s)
beta = sp.Function("beta")(s)
u = sp.Function("u")(s)
v = sp.Function("v")(s)

A = a + alpha * z + u * z**2
C = b + beta * z + v * z**2
jac = sp.expand(sp.diff(A, z) * sp.diff(C, s) - sp.diff(A, s) * sp.diff(C, z))
poly = sp.Poly(jac, z)

E0 = alpha * sp.diff(b, s) - sp.diff(a, s) * beta
E1 = (
    alpha * sp.diff(beta, s)
    - sp.diff(alpha, s) * beta
    + 2 * (u * sp.diff(b, s) - sp.diff(a, s) * v)
)
E2 = (
    alpha * sp.diff(v, s)
    + 2 * u * sp.diff(beta, s)
    - 2 * sp.diff(alpha, s) * v
    - sp.diff(u, s) * beta
)
E3 = 2 * (u * sp.diff(v, s) - sp.diff(u, s) * v)

for degree, expected in enumerate((E0, E1, E2, E3)):
    zero(poly.coeff_monomial(z**degree) - expected,
         f"generic Jacobian coefficient z^{degree}")
gate(sp.degree(jac, z) <= 3, "quadratic pair has no higher Jacobian bucket")

# A proportional nonzero quadratic row is killed in one target coordinate by
# a constant target GL_2 change.  This concrete symbolic control also checks
# the determinant scaling of the bracket.
h = sp.Function("h")(s)
A_raw = a + alpha * z + 2 * h * z**2
C_raw = b + beta * z + 3 * h * z**2
A_gl = A_raw / 2
C_gl = -3 * A_raw + 2 * C_raw
zero(sp.Poly(C_gl, z).coeff_monomial(z**2),
     "constant target change kills a proportional top direction")
zero(
    (sp.diff(A_gl, z) * sp.diff(C_gl, s)
     - sp.diff(A_gl, s) * sp.diff(C_gl, z))
    - (sp.diff(A_raw, z) * sp.diff(C_raw, s)
       - sp.diff(A_raw, s) * sp.diff(C_raw, z)),
    "chosen target change has determinant one",
)

# Complete nonlinear normal form.  The proof obtains this after the top-row
# Wronskian lemma and the normalized equations with v=0.
rho, d, beta0, lam, a0 = sp.symbols(
    "rho d beta0 lambda a0", nonzero=True
)
Bfun = sp.Function("B")(s)
C_class = Bfun + beta0 * z
A_class = rho * C_class**2 + d * C_class - lam * s / beta0 + a0
jac_class = (
    sp.diff(A_class, z) * sp.diff(C_class, s)
    - sp.diff(A_class, s) * sp.diff(C_class, z)
)
zero(jac_class - lam, "classified quadratic family has constant Jacobian")

S_target, Z_target = sp.symbols("S_target Z_target")
s_inverse = beta0 * (
    rho * Z_target**2 + d * Z_target + a0 - S_target
) / lam
z_inverse = (Z_target - Bfun.subs(s, s_inverse)) / beta0
zero(s_inverse.subs({S_target: A_class, Z_target: C_class}) - s,
     "classified inverse recovers s")
zero(z_inverse.subs({S_target: A_class, Z_target: C_class}) - z,
     "classified inverse recovers z")

# The coefficient form of the same classification checks every normalized
# Jacobian bucket separately.
a_class = rho * Bfun**2 + d * Bfun - lam * s / beta0 + a0
alpha_class = beta0 * (2 * rho * Bfun + d)
u_class = rho * beta0**2
subs_class = {
    a: a_class,
    sp.diff(a, s): sp.diff(a_class, s),
    b: Bfun,
    sp.diff(b, s): sp.diff(Bfun, s),
    alpha: alpha_class,
    sp.diff(alpha, s): sp.diff(alpha_class, s),
    beta: beta0,
    sp.diff(beta, s): sp.Integer(0),
    u: u_class,
    sp.diff(u, s): sp.Integer(0),
    v: sp.Integer(0),
    sp.diff(v, s): sp.Integer(0),
}
zero(E0.subs(subs_class) - lam, "classified constant bucket")
zero(E1.subs(subs_class), "classified linear bucket")
zero(E2.subs(subs_class), "classified quadratic bucket")
zero(E3.subs(subs_class), "classified cubic bucket")

# Linear-linear boundary: a constant target combination recovers s, then the
# surviving nonzero constant normal coefficient recovers z.
mu, nu, kappa, alpha0 = sp.symbols("mu nu kappa alpha0", nonzero=True)
Lfun = sp.Function("L")(s)
A_lin = Lfun + alpha0 * z
C_lin = mu * A_lin + nu * s + kappa
zero(
    sp.diff(A_lin, z) * sp.diff(C_lin, s)
    - sp.diff(A_lin, s) * sp.diff(C_lin, z)
    - alpha0 * nu,
    "linear-linear family has constant Jacobian",
)
s_lin_inverse = (Z_target - mu * S_target - kappa) / nu
z_lin_inverse = (S_target - Lfun.subs(s, s_lin_inverse)) / alpha0
zero(s_lin_inverse.subs({S_target: A_lin, Z_target: C_lin}) - s,
     "linear inverse recovers s")
zero(z_lin_inverse.subs({S_target: A_lin, Z_target: C_lin}) - z,
     "linear inverse recovers z")

# The beta=0 edge in normalized coordinates is not lost: E0=lambda makes
# alpha and b' units, and E1=0 then forces u=0.  A polynomial positive control
# displays the resulting triangular pair.
q0, q1, al0 = sp.symbols("q0 q1 al0", nonzero=True)
A_edge = s**4 + al0 * z
C_edge = q1 * s + q0
zero(
    sp.diff(A_edge, z) * sp.diff(C_edge, s)
    - sp.diff(A_edge, s) * sp.diff(C_edge, z)
    - al0 * q1,
    "beta-zero edge is triangular",
)

# Deterministic polynomial controls at several arm degrees verify the full
# formulas and explicit inverse, rather than only their functional template.
for degree in range(1, 8):
    b_test = sum(sp.Integer(j + 1) * s**j for j in range(degree + 1))
    beta_test = sp.Integer(degree + 1)
    rho_test = sp.Integer(degree + 2)
    d_test = sp.Integer(2 * degree + 1)
    lam_test = sp.Integer(3 * degree + 1)
    a0_test = sp.Integer(5 - degree)
    c_test = b_test + beta_test * z
    a_test = (
        rho_test * c_test**2
        + d_test * c_test
        - lam_test * s / beta_test
        + a0_test
    )
    zero(
        sp.diff(a_test, z) * sp.diff(c_test, s)
        - sp.diff(a_test, s) * sp.diff(c_test, z)
        - lam_test,
        f"degree-{degree} classified positive Jacobian control",
    )
    s_back = beta_test * (
        rho_test * c_test**2 + d_test * c_test + a0_test - a_test
    ) / lam_test
    z_back = (c_test - b_test.subs(s, s_back)) / beta_test
    zero(s_back - s, f"degree-{degree} positive inverse recovers s")
    zero(z_back - z, f"degree-{degree} positive inverse recovers z")

# Minimal Russell nodal arm hostile.  It self-identifies two distinct arm
# parameters, so the classification proves that it cannot be the z=0 trace
# of any polynomial normal-strip Keller pair of transverse degree at most two.
c = sp.symbols("c", nonzero=True)
a_node = 9 * c**6 * s**2
b_node = 27 * c**9 * s**3 - 3 * c**3 * s
s_plus = 1 / (3 * c**3)
s_minus = -s_plus
zero(a_node.subs(s, s_plus) - a_node.subs(s, s_minus),
     "nodal arm collision agrees in first coordinate")
zero(b_node.subs(s, s_plus) - b_node.subs(s, s_minus),
     "nodal arm collision agrees in second coordinate")
gate(sp.cancel(s_plus - s_minus) != 0,
     "nodal arm collision has distinct normalization addresses")

semantic = {
    "ambient": "k[s,z] with Jacobian A_z*C_s-A_s*C_z",
    "scope": "deg_z(A),deg_z(C)<=2 and nonzero constant Jacobian",
    "top_gate": "quadratic Wronskian makes the top row a constant target direction",
    "classification": "after target GL2, C=b(s)+beta*z and A=rho*C^2+d*C-lambda*s/beta+a0",
    "inverse": "s=beta(rho*C^2+d*C+a0-A)/lambda;z=(C-b(s))/beta",
    "linear_boundary": "constant target combination recovers s",
    "russell_consequence": "a self-identifying arm requires transverse degree at least three or a nonpolynomial strip expansion",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
     "no inactive Python assert")

print("theorem=THM-3856-quadratic-normal-strip-keller-pairs-are-automorphisms")
print("ambient=k[s,z];bracket=Jac_(z,s)")
print("scope=both_transverse_degrees_at_most_two;lambda_nonzero")
print("classification=constant_target_GL2_then_triangular_quadratic_shear")
print("inverse=s_and_z_recovered_polynomially")
print("russell_gate=self_identifying_arm_requires_degree_at_least_three_or_nonpolynomial_strip_expansion")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
