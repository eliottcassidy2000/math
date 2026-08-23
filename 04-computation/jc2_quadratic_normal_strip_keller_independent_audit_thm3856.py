#!/usr/bin/env python3
"""Independent exact audit for THM-3856's quadratic normal-strip theorem.

The proof packet is independent of the Russell surface.  Over a
characteristic-zero field, write a polynomial pair of transverse degree at
most two as

    A=a(s)+alpha(s)z+u(s)z^2,
    C=b(s)+beta(s)z+v(s)z^2.

The exact coefficient equations split according to which of u,v vanish.
Wronskian integration and elementary target shears reduce every branch to a
pair linear in z; the latter has an explicit triangular inverse.  The final
block specializes independently to the THM-3846 minimal nodal arm and checks
that its coefficient ideals are empty in a frozen bounded universe.
"""

from __future__ import annotations

import ast
import hashlib
import json
import sys
from pathlib import Path

import sympy as sp


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


GATES = 0


def gate(condition: object, label: str) -> None:
    global GATES
    GATES += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(label)


def zero(expression: sp.Expr, label: str) -> None:
    gate(sp.cancel(sp.factor(expression)) == 0, label)


def jacobian(left: sp.Expr, right: sp.Expr, z: sp.Symbol, s: sp.Symbol) -> sp.Expr:
    return sp.expand(
        sp.diff(left, z) * sp.diff(right, s)
        - sp.diff(left, s) * sp.diff(right, z)
    )


# ---------------------------------------------------------------------------
# 1. Universal coefficient packet.
# ---------------------------------------------------------------------------
s, z = sp.symbols("s z")
lam = sp.symbols("lambda", nonzero=True)
a = sp.Function("a")(s)
b = sp.Function("b")(s)
alpha = sp.Function("alpha")(s)
beta = sp.Function("beta")(s)
u = sp.Function("u")(s)
v = sp.Function("v")(s)

A = a + alpha * z + u * z**2
C = b + beta * z + v * z**2
remainder = sp.Poly(jacobian(A, C, z, s) - lam, z)

E0 = alpha * sp.diff(b, s) - sp.diff(a, s) * beta - lam
E1 = (
    alpha * sp.diff(beta, s)
    - sp.diff(alpha, s) * beta
    + 2 * u * sp.diff(b, s)
    - 2 * sp.diff(a, s) * v
)
E2 = (
    alpha * sp.diff(v, s)
    - 2 * sp.diff(alpha, s) * v
    + 2 * u * sp.diff(beta, s)
    - sp.diff(u, s) * beta
)
E3 = 2 * (u * sp.diff(v, s) - sp.diff(u, s) * v)

zero(remainder.coeff_monomial(1) - E0, "constant coefficient")
zero(remainder.coeff_monomial(z) - E1, "linear coefficient")
zero(remainder.coeff_monomial(z**2) - E2, "quadratic coefficient")
zero(remainder.coeff_monomial(z**3) - E3, "cubic coefficient")
gate(remainder.degree() <= 3, "complete transverse coefficient packet")


# ---------------------------------------------------------------------------
# 2. Exact Wronskian integrations and target-shear reductions.
# ---------------------------------------------------------------------------
kappa, tau = sp.symbols("kappa tau", nonzero=True)

# If u,v are both nonzero, E3=0 says (v/u)'=0 and hence v=kappa*u.
zero(
    sp.diff(v / u, s)
    - (u * sp.diff(v, s) - sp.diff(u, s) * v) / u**2,
    "top Wronskian is derivative of v/u",
)
D = sp.expand(C - kappa * A)
d = b - kappa * a
gamma = beta - kappa * alpha
zero(D.subs(v, kappa * u) - (d + gamma * z),
     "proportional quadratic direction is removed by C-kappa*A")
zero(jacobian(A, D, z, s) - jacobian(A, C, z, s),
     "first target shear preserves Jacobian")

# The new linear coefficient gamma cannot vanish in this branch.  If it did,
# D=d(s) and the z coefficient of Jac(A,D) would force d'=0 (because u!=0),
# after which the constant coefficient could not be lambda!=0.
gamma_zero_packet = sp.Poly(
    jacobian(a + alpha * z + u * z**2, d, z, s) - lam,
    z,
)
zero(
    gamma_zero_packet.coeff_monomial(1)
    - (alpha * sp.diff(d, s) - lam),
    "gamma-zero constant packet",
)
zero(
    gamma_zero_packet.coeff_monomial(z) - 2 * u * sp.diff(d, s),
    "gamma-zero top packet",
)

# After that shear one has a quadratic first coordinate and a linear second
# coordinate.  Its top equation is 2*u*gamma'-u'*gamma=0, which integrates
# to u=tau*gamma^2.  Subtracting tau*D^2 removes the last quadratic term.
F = a + alpha * z + u * z**2
G = d + gamma * z
FG_remainder = sp.Poly(jacobian(F, G, z, s) - lam, z)
zero(
    FG_remainder.coeff_monomial(z**2)
    - (2 * u * sp.diff(gamma, s) - sp.diff(u, s) * gamma),
    "one-sided top equation",
)
zero(
    sp.diff(u / gamma**2, s)
    - (sp.diff(u, s) * gamma - 2 * u * sp.diff(gamma, s)) / gamma**3,
    "one-sided equation is derivative of u/gamma^2",
)
F_linearized = sp.expand(F.subs(u, tau * gamma**2) - tau * G**2)
gate(sp.degree(F_linearized, z) <= 1,
     "quadratic first coordinate removed by F-tau*G^2")
zero(
    jacobian(F - tau * G**2, G, z, s) - jacobian(F, G, z, s),
    "second target shear preserves Jacobian",
)

# The opposite one-sided branch is symmetric.  If u=0 and v!=0, its top
# equation integrates to v=tau*alpha^2 and C-tau*A^2 is linear in z.
AC_u_zero = sp.Poly(
    (jacobian(A, C, z, s) - lam)
    .subs({u: 0, sp.diff(u, s): 0})
    .doit(),
    z,
)
zero(
    AC_u_zero.coeff_monomial(z**2)
    - (alpha * sp.diff(v, s) - 2 * sp.diff(alpha, s) * v),
    "opposite one-sided top equation",
)
zero(
    sp.diff(v / alpha**2, s)
    - (sp.diff(v, s) * alpha - 2 * v * sp.diff(alpha, s)) / alpha**3,
    "opposite equation is derivative of v/alpha^2",
)
C_linearized = sp.expand(C.subs({u: 0, v: tau * alpha**2}) - tau * (a + alpha * z) ** 2)
gate(sp.degree(C_linearized, z) <= 1,
     "quadratic second coordinate removed by C-tau*A^2")
zero(
    jacobian(A, C - tau * A**2, z, s) - jacobian(A, C, z, s),
    "opposite target shear preserves Jacobian",
)

# The divided linear coefficients in both one-sided integrations are forced
# nonzero.  These packets retain the alpha=0 and beta=0 edge cases instead of
# silently cancelling them.
alpha_zero_packet = sp.Poly(
    jacobian(a, b + beta * z + v * z**2, z, s) - lam,
    z,
)
zero(
    alpha_zero_packet.coeff_monomial(1)
    - (-sp.diff(a, s) * beta - lam),
    "alpha-zero constant packet",
)
zero(
    alpha_zero_packet.coeff_monomial(z)
    - (-2 * sp.diff(a, s) * v),
    "alpha-zero top packet",
)
beta_zero_packet = sp.Poly(
    jacobian(a + alpha * z + u * z**2, b, z, s) - lam,
    z,
)
zero(
    beta_zero_packet.coeff_monomial(1)
    - (alpha * sp.diff(b, s) - lam),
    "beta-zero constant packet",
)
zero(
    beta_zero_packet.coeff_monomial(z) - 2 * u * sp.diff(b, s),
    "beta-zero top packet",
)


# ---------------------------------------------------------------------------
# 3. Linear-normal endpoint and explicit inverses.
# ---------------------------------------------------------------------------
f = sp.Function("f")(s)
g = sp.Function("g")(s)
p = sp.Function("p")(s)
q = sp.Function("q")(s)
L1 = f + p * z
L2 = g + q * z
linear_remainder = sp.Poly(jacobian(L1, L2, z, s) - lam, z)
zero(
    linear_remainder.coeff_monomial(1)
    - (p * sp.diff(g, s) - sp.diff(f, s) * q - lam),
    "linear endpoint constant equation",
)
zero(
    linear_remainder.coeff_monomial(z)
    - (p * sp.diff(q, s) - sp.diff(p, s) * q),
    "linear endpoint Wronskian equation",
)
zero(
    sp.diff(q / p, s)
    - (p * sp.diff(q, s) - sp.diff(p, s) * q) / p**2,
    "linear Wronskian is derivative of q/p",
)

# Three deterministic inverse controls realize all endpoint cases.
X, Y = sp.symbols("X Y")

# p=0: the first coordinate is affine in s and the second is triangular.
F_p0 = 2 * s + 3
G_p0 = s**3 - s + 5 * z
s_p0 = (X - 3) / 2
z_p0 = (Y - (s_p0**3 - s_p0)) / 5
zero(jacobian(F_p0, G_p0, z, s) + 10, "p=0 positive Jacobian")
zero(s_p0.subs({X: F_p0, Y: G_p0}) - s, "p=0 inverse recovers s")
zero(z_p0.subs({X: F_p0, Y: G_p0}) - z, "p=0 inverse recovers z")

# q=0: the second coordinate is affine in s.
F_q0 = s**2 + z
G_q0 = 3 * s - 2
s_q0 = (Y + 2) / 3
z_q0 = X - s_q0**2
zero(jacobian(F_q0, G_q0, z, s) - 3, "q=0 positive Jacobian")
zero(s_q0.subs({X: F_q0, Y: G_q0}) - s, "q=0 inverse recovers s")
zero(z_q0.subs({X: F_q0, Y: G_q0}) - z, "q=0 inverse recovers z")

# p,q nonzero and proportional: L2-kappa*L1 is affine in s.
F_prop = s**3 + s + 2 * z
G_prop = 7 * F_prop + 5 * s - 4
s_prop = (Y - 7 * X + 4) / 5
z_prop = (X - (s_prop**3 + s_prop)) / 2
zero(jacobian(F_prop, G_prop, z, s) - 10, "proportional positive Jacobian")
zero(s_prop.subs({X: F_prop, Y: G_prop}) - s,
     "proportional inverse recovers s")
zero(z_prop.subs({X: F_prop, Y: G_prop}) - z,
     "proportional inverse recovers z")


# ---------------------------------------------------------------------------
# 4. Positive quadratic controls for every top-direction branch.
# ---------------------------------------------------------------------------
# u=0, v!=0: add a square of the first coordinate to the second.
A_second = s + z
C_second = -z + 2 * A_second**2
zero(jacobian(A_second, C_second, z, s) - 1,
     "quadratic-second-coordinate positive control")
gate(sp.degree(A_second, z) == 1 and sp.degree(C_second, z) == 2,
     "quadratic-second support")

# v=0, u!=0: add a square of the second coordinate to the first.
C_first = s + z
A_first = z + 2 * C_first**2
zero(jacobian(A_first, C_first, z, s) - 1,
     "quadratic-first-coordinate positive control")
gate(sp.degree(A_first, z) == 2 and sp.degree(C_first, z) == 1,
     "quadratic-first support")

# u,v both nonzero: compose both target shears and verify an explicit inverse.
D_both = s + z
A_both = z + 2 * D_both**2
C_both = 3 * A_both + D_both
zero(jacobian(A_both, C_both, z, s) - 1,
     "two-sided quadratic positive control")
gate(sp.degree(A_both, z) == 2 and sp.degree(C_both, z) == 2,
     "two-sided quadratic support")
D_inverse = Y - 3 * X
z_inverse = X - 2 * D_inverse**2
s_inverse = D_inverse - z_inverse
zero(s_inverse.subs({X: A_both, Y: C_both}) - s,
     "two-sided inverse recovers s")
zero(z_inverse.subs({X: A_both, Y: C_both}) - z,
     "two-sided inverse recovers z")


# ---------------------------------------------------------------------------
# 5. THM-3846/3849 minimal nodal hostile and independent coefficient scan.
# ---------------------------------------------------------------------------
t, w = sp.symbols("t w")
h = sp.Function("h")(t)
a_node = t**2
b_node = t**3 - t
alpha_node = -1 + 2 * t * h
beta_node = -sp.Rational(3, 2) * t + (3 * t**2 - 1) * h
W_node = sp.expand(
    alpha_node * sp.diff(beta_node, t)
    - sp.diff(alpha_node, t) * beta_node
)
zero(
    alpha_node * sp.diff(b_node, t)
    - sp.diff(a_node, t) * beta_node
    - 1,
    "all nodal Bezout rows",
)
zero(
    W_node
    - (
        sp.diff(h, t)
        - 6 * t * h
        + (6 * t**2 + 2) * h**2
        + sp.Rational(3, 2)
    ),
    "nodal Wronskian formula",
)

# The implicit nodal cubic and its conductor polynomial from THM-3849.
U, V = sp.symbols("U V")
Delta = V**2 - U * (U - 1) ** 2
kappa_node = -(t**2 - 1)
zero(Delta.subs({U: a_node, V: b_node}), "nodal implicit equation")
zero(
    sp.diff(Delta, U).subs({U: a_node, V: b_node})
    - kappa_node * sp.diff(b_node, t),
    "nodal conductor U-gradient",
)
zero(
    sp.diff(Delta, V).subs({U: a_node, V: b_node})
    + kappa_node * sp.diff(a_node, t),
    "nodal conductor V-gradient",
)
gate(
    (a_node.subs(t, 1), b_node.subs(t, 1))
    == (a_node.subs(t, -1), b_node.subs(t, -1)),
    "nodal arm self-identification",
)

# Frozen finite exact scan.  The general theorem already proves emptiness;
# these unit ideals independently test the coefficient equations themselves.
nodal_scan = []
for h_degree, uv_degree in ((0, 0), (0, 1), (1, 1), (1, 2), (2, 2), (2, 3)):
    hs = sp.symbols(f"h0:{h_degree + 1}")
    us = sp.symbols(f"u0:{uv_degree + 1}")
    vs = sp.symbols(f"v0:{uv_degree + 1}")
    hp = sum(hs[index] * t**index for index in range(h_degree + 1))
    up = sum(us[index] * t**index for index in range(uv_degree + 1))
    vp = sum(vs[index] * t**index for index in range(uv_degree + 1))
    alp = alpha_node.subs(h, hp)
    bep = beta_node.subs(h, hp)
    wp = W_node.subs({h: hp, sp.diff(h, t): sp.diff(hp, t)}).doit()
    equations = (
        sp.expand(
            wp
            + 2 * (up * sp.diff(b_node, t) - sp.diff(a_node, t) * vp)
        ),
        sp.expand(
            alp * sp.diff(vp, t)
            - 2 * sp.diff(alp, t) * vp
            + 2 * up * sp.diff(bep, t)
            - sp.diff(up, t) * bep
        ),
        sp.expand(up * sp.diff(vp, t) - sp.diff(up, t) * vp),
    )
    coefficient_equations = []
    for equation in equations:
        coefficient_equations.extend(sp.Poly(equation, t).all_coeffs())
    unknowns = (*hs, *us, *vs)
    basis = sp.groebner(coefficient_equations, *unknowns, order="grevlex")
    unit = len(basis.polys) == 1 and basis.polys[0].as_expr() == 1
    gate(unit, f"nodal unit ideal H={h_degree},D={uv_degree}")
    nodal_scan.append(
        {
            "h_degree_cap": h_degree,
            "uv_degree_cap": uv_degree,
            "unknowns": len(unknowns),
            "equations": len(coefficient_equations),
        }
    )

# A finite Catalan truncation is a hostile, not a solution: it cancels the
# linear normal debt and leaves exact quadratic/cubic residual coefficients.
Z2 = w - sp.Rational(3, 4) * w**2
A_trunc = a_node - Z2
C_trunc = b_node - sp.Rational(3, 2) * t * Z2
truncation_debt = sp.expand(jacobian(A_trunc, C_trunc, w, t) - 1)
zero(
    truncation_debt
    - (-sp.Rational(27, 8) * w**2 + sp.Rational(27, 16) * w**3),
    "quadratic Catalan truncation hostile",
)


# ---------------------------------------------------------------------------
# 6. Frozen semantic packet.
# ---------------------------------------------------------------------------
source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "no inactive Python assert",
)

semantic = {
    "universe": "char(k)=0;A,C in k[s,z];deg_z(A),deg_z(C)<=2;Jac=lambda in k*",
    "coefficients": "E0..E3 exact",
    "classification": "proportional top direction;one-sided square shear;linear-normal triangular inverse",
    "zero_cases": "u=v=0,u=0,v=0,uv!=0 all explicit",
    "division_edges": "gamma=0,alpha=0,beta=0 packets retained and inconsistent with lambda nonzero",
    "positive_controls": "linear p=0/q=0/proportional and quadratic first/second/two-sided",
    "nodal": "arm=(t2,t3-t);all Bezout rows;conductor=-(t2-1);bounded ideals unit",
    "scope": "polynomial canonical strip only;global Russell B algebraization and transverse degree >=3 remain open",
    "nodal_scan": nodal_scan,
}
semantic_sha = hashlib.sha256(
    json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
).hexdigest()

print("THM3856_QUADRATIC_NORMAL_STRIP_INDEPENDENT_AUDIT")
print("status=PASS_PROVED_CLASSIFICATION+VERIFIED_EXACT")
print("universe=char0;A_C_in_k[s,z];deg_z_each_at_most_2;Jac=lambda_nonzero")
print("coefficient_packet=E0_E1_E2_E3_complete")
print("top_cases=both_nonzero_proportional;one_sided_square;both_zero_linear")
print("zero_edges=gamma_zero_alpha_zero_beta_zero_explicitly_inconsistent")
print("reduction=target_shears_to_linear_normal_pair")
print("linear_endpoint=three_explicit_inverse_cases")
print("positive_controls=linear_3;quadratic_first_second_two_sided")
print("nodal_hostile=all_Bezout_rows;conductor=-(t^2-1);self_identifies_at_+-1")
print("nodal_groebner_caps=(H,D)=(0,0),(0,1),(1,1),(1,2),(2,2),(2,3);all_unit")
print("Catalan_quadratic_debt=-27w^2/8+27w^3/16")
print("scope=canonical_polynomial_strip_only;global_Russell_and_deg_z>=3_open")
print(f"semantic_sha256={semantic_sha}")
print(f"GATES={GATES}")
print("RESULT=PASS")
