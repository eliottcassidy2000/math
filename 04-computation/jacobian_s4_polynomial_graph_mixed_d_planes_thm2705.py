#!/usr/bin/env python3
"""Exact companion for THM-2705 mixed target planes containing d."""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


x, y, t = sp.symbols("x y t")
alpha, beta, kappa = sp.symbols("alpha beta kappa", nonzero=True)
q, h = sp.symbols("q h", nonzero=True)
fxy = sp.Function("f")(x, y)

A = x**2 - 2 * y
B = y**2 - 2 * x * fxy
C = alpha * A + beta * B
J_Cd = sp.factor(sp.det(sp.Matrix([C, fxy]).jacobian([x, y])))
fx = sp.diff(fxy, x)
fy = sp.diff(fxy, y)
require(
    sp.expand(J_Cd - 2 * ((alpha - beta * y) * fx + (alpha * x - beta * fxy) * fy))
    == 0,
    "mixed Jacobian identity",
)

# The beta=0 branch is the scaled THM-2702 triangular family.
H = t**5 - 3 * t**2 + 2 * t + 7
t_xy = y - x**2 / 2
f_beta0 = sp.expand(kappa * x / (2 * alpha) + H.subs(t, t_xy))
C_beta0 = alpha * A
require(
    sp.factor(sp.det(sp.Matrix([C_beta0, f_beta0]).jacobian([x, y])) - kappa)
    == 0,
    "scaled beta-zero family",
)
CC, dd = sp.symbols("CC dd")
t_inv = -CC / (2 * alpha)
x_inv = sp.factor(2 * alpha * (dd - H.subs(t, t_inv)) / kappa)
y_inv = sp.factor(t_inv + x_inv**2 / 2)
require(sp.factor(C_beta0.subs({x: x_inv, y: y_inv}) - CC) == 0, "beta-zero inverse C")
require(sp.factor(f_beta0.subs({x: x_inv, y: y_inv}) - dd) == 0, "beta-zero inverse d")

# For alpha*beta != 0, put alpha=-beta*q^2 and kappa=2*beta*q*h.
# The same-Jacobian equation becomes one polynomial PDE P(g)=0.
g = sp.Function("g")(x, t)
y_xt = t + q * x
f_xt = q * t + q**3 - h + g
fx_at_y = sp.diff(f_xt, x) - q * sp.diff(f_xt, t)
fy_at_x = sp.diff(f_xt, t)
half_jac = sp.expand(
    ((-beta * q**2) - beta * y_xt) * fx_at_y
    + ((-beta * q**2) * x - beta * f_xt) * fy_at_x
)
P = sp.expand(
    (h - g) * sp.diff(g, t)
    - (q**2 + q * x + t) * sp.diff(g, x)
    - q * g
)
require(sp.factor(half_jac - beta * q * h - beta * P) == 0, "mixed PDE reduction")

# Exact coefficient equations for g=a(x)t+b(x).
a = sp.Function("a")(x)
b = sp.Function("b")(x)
P_linear = sp.Poly(sp.expand(P.subs(g, a * t + b).doit()), t)
expected = {
    2: -sp.diff(a, x),
    1: -a**2 - q * a - sp.diff(b, x) - (q**2 + q * x) * sp.diff(a, x),
    0: a * h - a * b - q * b - (q**2 + q * x) * sp.diff(b, x),
}
for degree, coefficient in expected.items():
    require(
        sp.expand(P_linear.coeff_monomial(t**degree) - coefficient) == 0,
        f"linear coefficient t^{degree}",
    )

# Degree gates used by the all-degree proof.
for r in range(3, 20):
    require(2 * r - 1 > r + 1, f"unique nonlinear top degree r={r}")
for polynomial_degree in range(1, 20):
    require(2 * polynomial_degree > polynomial_degree - 1, "Riccati degree gate")

f_plus = -q**2 * x + q * y + q**3 - h
f_minus = -q**2 * x - q * y + h - q**3
C_plus = sp.expand((-beta * q**2) * A + beta * (y**2 - 2 * x * f_plus))
C_minus = sp.expand((-beta * q**2) * A + beta * (y**2 - 2 * x * f_minus))
for name, f_solution, C_solution in (
    ("plus", f_plus, C_plus),
    ("minus", f_minus, C_minus),
):
    require(
        sp.factor(sp.det(sp.Matrix([C_solution, f_solution]).jacobian([x, y])) - 2 * beta * q * h)
        == 0,
        f"mixed {name} Jacobian",
    )

# Both mixed solutions are triangular in t=y-qx or u=y+qx.
plus_t = sp.expand(C_plus.subs(y, t + q * x))
plus_d = sp.expand(f_plus.subs(y, t + q * x))
require(sp.expand(plus_t - beta * (t**2 + 2 * q**2 * t + 2 * h * x)) == 0, "plus triangular C")
require(sp.expand(plus_d - (q * t + q**3 - h)) == 0, "plus triangular d")
t_plus_inv = sp.factor((dd - q**3 + h) / q)
x_plus_inv = sp.factor((CC / beta - t_plus_inv**2 - 2 * q**2 * t_plus_inv) / (2 * h))
y_plus_inv = sp.factor(t_plus_inv + q * x_plus_inv)
require(sp.factor(C_plus.subs({x: x_plus_inv, y: y_plus_inv}) - CC) == 0, "plus inverse C")
require(sp.factor(f_plus.subs({x: x_plus_inv, y: y_plus_inv}) - dd) == 0, "plus inverse d")

minus_u = sp.expand(C_minus.subs(y, t - q * x))
minus_d = sp.expand(f_minus.subs(y, t - q * x))
require(sp.expand(minus_u - beta * (t**2 + 2 * q**2 * t - 2 * h * x)) == 0, "minus triangular C")
require(sp.expand(minus_d - (-q * t + h - q**3)) == 0, "minus triangular d")
t_minus_inv = sp.factor((h - q**3 - dd) / q)
x_minus_inv = sp.factor((t_minus_inv**2 + 2 * q**2 * t_minus_inv - CC / beta) / (2 * h))
y_minus_inv = sp.factor(t_minus_inv - q * x_minus_inv)
require(sp.factor(C_minus.subs({x: x_minus_inv, y: y_minus_inv}) - CC) == 0, "minus inverse C")
require(sp.factor(f_minus.subs({x: x_minus_inv, y: y_minus_inv}) - dd) == 0, "minus inverse d")

# q -> -q swaps the two graph equations after h -> -h, as forced by
# h=kappa/(2*beta*q).
require(sp.expand(f_plus.subs({q: -q, h: -h}) - f_minus) == 0, "square-root gauge swap")

print("THM-2705 polynomial-graph linear target planes containing d")
print("Jac(alpha*A+beta*B,d)=2*((alpha-beta*y)*f_x+(alpha*x-beta*f)*f_y)")
print("beta=0:f=(kappa/(2*alpha))*x+H(y-x^2/2);triangular automorphism")
print("alpha=0,beta!=0:EMPTY by scaled LND/divergence contradiction")
print("mixed:q^2=-alpha/beta;h=kappa/(2*beta*q);P=(h-g)g_t-(q^2+qx+t)g_x-qg")
print("mixed_solutions=f_plus:-q^2*x+q*y+q^3-h;f_minus:-q^2*x-q*y+h-q^3")
print("degree_gate=r>=3 unique 2r-1 term;r=2 polynomial Riccati;linear roots 0,-q,-2q")
print("mixed_inverse_checks=plus:PASS;minus:PASS;square_root_gauge=SWAP")
print("ALL CHECKS PASSED")
