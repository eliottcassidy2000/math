#!/usr/bin/env python3
"""Exact companion for THM-3837's line--hyperbola contact obstruction."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def check(condition: object, label: str) -> None:
    global CHECKS
    if condition is not True and condition != sp.S.true:
        raise AssertionError(label)
    CHECKS += 1


def same(lhs: object, rhs: object, label: str) -> None:
    check(sp.cancel(sp.expand(lhs - rhs)) == 0, label)  # type: ignore[operator]


alpha = sp.symbols("alpha")
root_relation = 3 * alpha**3 + 7 * alpha**2 + 1


def reduce_alpha(expression: object) -> sp.Expr:
    """Reduce a rational expression modulo the cubic root relation."""

    numerator, denominator = sp.together(expression).as_numer_denom()
    coefficient_symbols = sorted(numerator.free_symbols - {alpha}, key=str)
    domain = (sp.QQ.frac_field(*coefficient_symbols)
              if coefficient_symbols else sp.QQ)
    remainder = sp.rem(
        sp.Poly(numerator, alpha, domain=domain),
        sp.Poly(root_relation, alpha, domain=domain),
    ).as_expr()
    return sp.cancel(remainder / denominator)


def same_mod_root(lhs: object, rhs: object, label: str) -> None:
    check(reduce_alpha(sp.together(lhs - rhs)) == 0, label)  # type: ignore[operator]


z = sp.symbols("z")
r = 3 * z**3 + 7 * z**2 + 1
Qz = 7 * z**2 + 3
az = sp.expand(Qz * r)
bz = sp.expand((z + 1) * (2 * z + 1) * (3 * z - 1))
Uz = 49 * z**2 - 9 * z + 7
collision = 140 * z**2 + 307 * z - 93

check(sp.discriminant(r, z) == -1615, "cubic slopes are simple")
check(sp.resultant(r, Qz, z) == 1615, "Q is a unit at every cubic slope")
check(sp.resultant(r, Uz, z) == 14535, "second-contact multiplier is a unit")
check(sp.resultant(r, collision, z) == -5817545,
      "the two second-contact values are incompatible")

Q = 7 * alpha**2 + 3
B = -Q
rho = 9 * alpha**2 + 14 * alpha
Bprime = 18 * alpha**2 + 14 * alpha
inverse_alpha = -3 * alpha**2 - 7 * alpha
inverse_B = (189 * alpha**2 + 147 * alpha - 767) / 1615
mu = -9 * alpha - 14

same_mod_root(bz.subs(z, alpha), B, "b(alpha)=-Q(alpha)")
same_mod_root(B, alpha * rho, "b(alpha)=alpha*r'(alpha)")
same_mod_root(alpha * inverse_alpha, 1, "explicit inverse of alpha")
same_mod_root(B * inverse_B, 1, "explicit inverse of b(alpha)")
same_mod_root(-B * inverse_alpha**2, mu, "selector boundary scalar mu")

x, y = sp.symbols("x y")
q0, qx, qy, d0, dx, dy = sp.symbols("q0 qx qy d0 dx dy")
q_affine = q0 + qx * x + qy * y
D_affine = d0 + dx * x + dy * y
u = x
v = x * y - 1
f = sp.expand(u * v)
k = sp.expand(x - v + f * q_affine)
h = sp.expand(alpha * k + f)
d = sp.expand(mu * x**2 + f * D_affine)

E0 = 1 + x * y * (y - 1)
m = sp.expand(E0 * q_affine - (y - 1) ** 2)
C = sp.expand(E0 + alpha * m)
same((x * y - 1) - y * x, -1, "line and hyperbola are comaximal")
same(E0 * (x - v) + f * (y - 1) ** 2, 1,
     "base CRT Bezout identity")
same(C * k - m * h, 1, "explicit determinant-one completion")


def A5(hh: sp.Expr, kk: sp.Expr) -> sp.Expr:
    return sp.expand((7 * hh**2 + 3 * kk**2)
                     * (3 * hh**3 + 7 * hh**2 * kk + kk**3))


def B3(hh: sp.Expr, kk: sp.Expr) -> sp.Expr:
    return sp.expand((hh + kk) * (2 * hh + kk) * (3 * hh - kk))


def residual(hh: sp.Expr, kk: sp.Expr, dd: sp.Expr) -> sp.Expr:
    return sp.expand(A5(hh, kk) - dd * (kk * B3(hh, kk) + hh**2 * dd))


same(k.subs(x, 0), 1, "row unit on the selected line")
same(k.subs(y, 1 / x), x, "row unit on the hyperbola")
same(d.subs(x, 0), 0, "selector vanishes on the selected line")
same_mod_root(
    (k * B3(h, k) + h**2 * d).subs(y, 1 / x),
    0,
    "complementary bracket vanishes on the hyperbola",
)

# First and second normal contacts on the selected line x=0.
epsilon, Y, X = sp.symbols("epsilon Y X")
x_u = epsilon
y_u = Y
v_u = x_u * y_u - 1
f_u = x_u * v_u
q_u = q0 + qx * x_u + qy * y_u
D_u = d0 + dx * x_u + dy * y_u
k_u = sp.expand(x_u - v_u + f_u * q_u)
h_u = sp.expand(alpha * k_u + f_u)
d_u = sp.expand(mu * x_u**2 + f_u * D_u)
E_u = residual(h_u, k_u, d_u)
E_u1 = sp.expand(E_u).coeff(epsilon, 1)
E_u2 = sp.expand(E_u).coeff(epsilon, 2)
same_mod_root(E_u1, -Q * (dy * Y + d0 + rho),
              "selected-line first contact")

# Contacts on the hyperbola v=0, in K[X,X^-1][[epsilon]].
v_v = epsilon
x_v = X
y_v = (1 + epsilon) / X
f_v = X * epsilon
q_v = q0 + qx * x_v + qy * y_v
D_v = d0 + dx * x_v + dy * y_v
k_v = sp.expand(x_v - v_v + f_v * q_v)
h_v = sp.expand(alpha * k_v + f_v)
d_v = sp.expand(mu * x_v**2 + f_v * D_v)
E_v = residual(h_v, k_v, d_v)
E_v1 = sp.expand(E_v).coeff(epsilon, 1)
E_v2 = sp.expand(E_v).coeff(epsilon, 2)
q_v0 = q0 + qx * X + qy / X
D_v0 = d0 + dx * X + dy / X
same(sp.diff(k_v, epsilon).subs(epsilon, 0) / X,
     q_v0 - 1 / X,
     "corrected hyperbola derivative k_f=q-1/X")

aprime = Q * rho
expected_v1 = X**5 * (
    aprime
    - mu * (
        2 * B * X * q_v0 - 2 * B
        + X * (Bprime - 2 * B * inverse_alpha)
        + alpha**2 * D_v0
    )
)
same_mod_root(E_v1, expected_v1, "hyperbola first contact")

q0_forced = reduce_alpha(
    (-alpha**2 * dx - Bprime + 2 * B * inverse_alpha)
    * inverse_B / 2
)
first_contact = {
    qx: 0,
    qy: 1 + alpha,
    dy: 0,
    d0: -rho,
    q0: q0_forced,
}
same_mod_root(
    2 * B * q0_forced,
    -alpha**2 * dx - Bprime + 2 * B * inverse_alpha,
    "first-contact q0 relation",
)
same_mod_root(E_u1.subs(first_contact), 0,
              "forced profiles kill the line first contact")
same_mod_root(E_v1.subs(first_contact), 0,
              "forced profiles kill the hyperbola first contact")

L_u = dx + 18 * alpha + 38 + 7 * alpha**2
L_v = dx + 325 * alpha - 55 + 147 * alpha**2
U = 49 * alpha**2 - 9 * alpha + 7
same_mod_root(E_u2.subs(first_contact), -Q * L_u,
              "selected-line second contact is -Q*L_u")
E_v2_forced = sp.expand(E_v2.subs(first_contact))
coefficient_X5 = sp.Poly(E_v2_forced, X).coeff_monomial(X**5)
same_mod_root(coefficient_X5, -U * L_v / 3,
              "hyperbola second contact X^5 coefficient")
same(L_v - L_u, 140 * alpha**2 + 307 * alpha - 93,
     "cross-difference of forced second-contact values")

# Exact characteristic-zero quadratic extension.  Use the universal Taylor
# coefficient of the selector residual to avoid a bounded Groebner inference.
qxx, qxy, qyy, dxx, dxy, dyy = sp.symbols(
    "qxx qxy qyy dxx dxy dyy"
)
tau = sp.symbols("tau")
K0, K1, K2 = sp.symbols("K0 K1 K2")
Delta1, Delta2 = sp.symbols("Delta1 Delta2")
T0, T1, T2 = sp.symbols("T0 T1 T2")
k_taylor = K0 + K1 * tau + K2 * tau**2
delta_taylor = Delta1 * tau + Delta2 * tau**2
h_taylor = alpha * k_taylor + delta_taylor
d_taylor = T0 + T1 * tau + T2 * tau**2
E_taylor = residual(h_taylor, k_taylor, d_taylor)
universal_first = sp.Poly(E_taylor, tau).coeff_monomial(tau)
universal_second = sp.Poly(E_taylor, tau).coeff_monomial(tau**2)

# Jets in the selected-line chart x=tau, y=Y.
quad_q_u0 = q0 + qy * Y + qyy * Y**2
quad_q_u1 = qx + qxy * Y
quad_D_u0 = d0 + dy * Y + dyy * Y**2
quad_D_u1 = dx + dxy * Y
quad_u_jets = {
    K0: 1,
    K1: 1 - Y - quad_q_u0,
    K2: -quad_q_u1 + Y * quad_q_u0,
    Delta1: -1,
    Delta2: Y,
    T0: 0,
    T1: -quad_D_u0,
    T2: mu + Y * quad_D_u0 - quad_D_u1,
}
quadratic_u1 = universal_first.subs(quad_u_jets)
quadratic_u2 = universal_second.subs(quad_u_jets)
same_mod_root(quadratic_u1, -Q * (quad_D_u0 + rho),
              "complete quadratic selected-line first contact")

# Jets in the hyperbola chart v=tau, x=X, y=(1+tau)/X.
quad_q_v0 = (q0 + qx * X + qy / X + qxx * X**2
             + qxy + qyy / X**2)
quad_q_v1 = qy / X + qxy + 2 * qyy / X**2
quad_D_v0 = (d0 + dx * X + dy / X + dxx * X**2
             + dxy + dyy / X**2)
quad_D_v1 = dy / X + dxy + 2 * dyy / X**2
quad_v_jets = {
    K0: X,
    K1: -1 + X * quad_q_v0,
    K2: X * quad_q_v1,
    Delta1: X,
    Delta2: 0,
    T0: mu * X**2,
    T1: X * quad_D_v0,
    T2: X * quad_D_v1,
}
quadratic_v1 = universal_first.subs(quad_v_jets)
quadratic_v2 = universal_second.subs(quad_v_jets)
expected_quadratic_v1 = X**5 * (
    aprime
    - mu * (
        2 * B * X * quad_q_v0 - 2 * B
        + X * (Bprime - 2 * B * inverse_alpha)
        + alpha**2 * quad_D_v0
    )
)
same_mod_root(quadratic_v1, expected_quadratic_v1,
              "complete quadratic hyperbola first contact")

q0_quadratic = q0_forced - qxy
qx_quadratic = reduce_alpha(-alpha**2 * dxx * inverse_B / 2)
qy_quadratic = reduce_alpha(
    1 + alpha - alpha**2 * dxy * inverse_B / 2
)
quadratic_first_contact = {
    qxx: 0,
    qyy: 0,
    d0: -rho,
    dy: 0,
    dyy: 0,
    qx: qx_quadratic,
    qy: qy_quadratic,
    q0: q0_quadratic,
}
same_mod_root(quadratic_u1.subs(quadratic_first_contact), 0,
              "quadratic first system kills the line contact")
same_mod_root(quadratic_v1.subs(quadratic_first_contact), 0,
              "quadratic first system kills the hyperbola contact")

quadratic_u2_forced = quadratic_u2.subs(quadratic_first_contact)
same_mod_root(quadratic_u2_forced, -Q * (dxy * Y + L_u),
              "quadratic selected second contact peels d_xy then L_u")
quadratic_v2_forced = reduce_alpha(
    quadratic_v2.subs(quadratic_first_contact)
)
quadratic_v2_poly = sp.Poly(quadratic_v2_forced, X)
coefficient_X8 = quadratic_v2_poly.coeff_monomial(X**8)
same_mod_root(coefficient_X8, alpha**2 * dxx**2 / 4,
              "quadratic hyperbola X^8 contact peels d_xx")
quadratic_peel = {dxy: 0, dxx: 0}
coefficient_X5_quadratic = sp.Poly(
    quadratic_v2_forced.subs(quadratic_peel), X
).coeff_monomial(X**5)
same_mod_root(coefficient_X5_quadratic, -U * L_v / 3,
              "quadratic hyperbola X^5 contact is the same L_v")
check(sp.diff(coefficient_X5_quadratic, qxy) == 0,
      "the remaining q_xy parameter cannot repair the terminal contact")

source = Path(__file__).read_text(encoding="utf-8")
check(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
      "no inactive Python assert")

semantic = {
    "ansatz": "n=1; u=x; v=xy-1; k=x-v+uv*q; h=alpha*k+uv; d=-b(alpha)x^2/alpha^2+uv*D",
    "profiles": "q,D of total degree at most two over the cubic number field",
    "first_contact": "four parameters dx,dxx,dxy,qxy survive",
    "correction": "on v=0, k_f=q-1/x, not q",
    "second_contact": "dxy=dxx=0, then line forces L_u=0 and hyperbola forces L_v=0; resultant=-5817545",
    "scope": "selector equation only; remaining intrinsic row and planar Jacobian OPEN",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

print("theorem=THM-3837-comaximal-line-hyperbola-affine-bichromatic-contact-nonentry")
print("fibre=u=x;v=xy-1;gcd(u,v)=1")
print("row=k=x-v+uv*q;determinant_completion=exact")
print("selector=d=-b(alpha)*x^2/alpha^2+uv*D")
print("first_contact=four_parameters_survive;corrected_k_f=q-1/x")
print("second_contact=d_xy=d_xx=0_then_L_u_L_v_incompatible;resultant=-5817545")
print("quadratic=characteristic_zero_exact;no_modular_inference")
print("scope=n1_total_degree_at_most2_qD_only;second_row_and_Jacobian_OPEN")
print(f"CHECKS={CHECKS}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
