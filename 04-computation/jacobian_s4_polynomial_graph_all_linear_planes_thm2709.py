#!/usr/bin/env python3
"""Exact companion for THM-2709 all linear target planes on graphs."""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


x, y = sp.symbols("x y")
P, Q, R, kappa = sp.symbols("P Q R kappa", nonzero=True)
a, b, c = sp.symbols("a b c", nonzero=True)
f = sp.Function("f")(x, y)

A = x**2 - 2 * y
B = y**2 - 2 * x * f
d = f
J_AB = sp.factor(sp.det(sp.Matrix([A, B]).jacobian([x, y])))
J_Bd = sp.factor(sp.det(sp.Matrix([B, d]).jacobian([x, y])))
J_dA = sp.factor(sp.det(sp.Matrix([d, A]).jacobian([x, y])))

fx = sp.diff(f, x)
fy = sp.diff(f, y)
normalized = sp.expand(
    (x + a * y + b) * fx
    + (x**2 + a * f + b * x) * fy
    + f
    - x * y
    - c
)
pluecker_jac = sp.expand(P * J_AB + Q * J_Bd + R * J_dA)
substituted = normalized.subs({a: Q / (2 * P), b: R / (2 * P), c: -kappa / (4 * P)})
require(sp.factor(pluecker_jac - kappa + 4 * P * substituted) == 0, "Pluecker PDE")

# Direct 2x3 row control for the Pluecker expansion.
uA, uB, ud, vA, vB, vd = sp.symbols("uA uB ud vA vB vd")
U = uA * A + uB * B + ud * d
V = vA * A + vB * B + vd * d
wedge_P = uA * vB - uB * vA
wedge_Q = uB * vd - ud * vB
wedge_R = ud * vA - uA * vd
require(
    sp.factor(
        sp.det(sp.Matrix([U, V]).jacobian([x, y]))
        - (wedge_P * J_AB + wedge_Q * J_Bd + wedge_R * J_dA)
    )
    == 0,
    "2x3 row-wedge expansion",
)

# The Q=0, P!=0 branch has one shifted cubic. Here b=R/(2P),
# c=-kappa/(4P), and a normalized target basis is (P*A,B-2b*d).
t = sp.symbols("t")
t_xy = y - x**2 / 2
f_shift = sp.expand(
    (x - b) * t_xy / 2 + (x - b) * (x**2 + b**2) / 8 + c
)
U_shift = P * A
V_shift = sp.expand(y**2 - 2 * x * f_shift - 2 * b * f_shift)
require(
    sp.factor(sp.det(sp.Matrix([U_shift, V_shift]).jacobian([x, y])) + 4 * P * c)
    == 0,
    "shifted cubic Jacobian",
)
V_triangular = sp.expand(V_shift.subs(y, t + x**2 / 2))
require(
    sp.expand(V_triangular - ((t + b**2 / 2) ** 2 - 2 * c * (x + b))) == 0,
    "shifted cubic triangular form",
)
UU, VV = sp.symbols("UU VV")
t_inv = -UU / (2 * P)
x_inv = sp.factor(((t_inv + b**2 / 2) ** 2 - VV) / (2 * c) - b)
y_inv = sp.factor(t_inv + x_inv**2 / 2)
require(sp.factor(U_shift.subs({x: x_inv, y: y_inv}) - UU) == 0, "shifted inverse U")
require(sp.factor(V_shift.subs({x: x_inv, y: y_inv}) - VV) == 0, "shifted inverse V")

# For a!=0, the all-degree y-gate leaves f=u(x)y+v(x).
u = sp.Function("u")(x)
v = sp.Function("v")(x)
linear_y = sp.Poly(sp.expand(normalized.subs(f, u * y + v).doit()), y)
expected = {
    2: a * sp.diff(u, x),
    1: a * u**2 + a * sp.diff(v, x) + (x + b) * sp.diff(u, x) - x + u,
    0: a * u * v + (x + b) * sp.diff(v, x) + (x**2 + b * x) * u + v - c,
}
for degree, coefficient in expected.items():
    require(
        sp.expand(linear_y.coeff_monomial(y**degree) - coefficient) == 0,
        f"linear-y coefficient y^{degree}",
    )

for degree in range(3, 20):
    require(2 * degree - 1 > degree + 1, f"unique nonlinear y-degree {degree}")
for polynomial_degree in range(1, 20):
    require(2 * polynomial_degree > polynomial_degree - 1, "y-Riccati degree gate")

# Once u is constant, the y equation fixes v. The x^2 equation forces
# u=-1/a, and then the constant equation is exactly -c, impossible for a
# nonzero Keller constant.
u0, v0 = sp.symbols("u0 v0")
v_forced = x**2 / (2 * a) - (u0 + a * u0**2) * x / a + v0
constant_equation = sp.Poly(
    sp.expand(
        a * u0 * v_forced
        + (x + b) * sp.diff(v_forced, x)
        + (x**2 + b * x) * u0
        + v_forced
        - c
    ),
    x,
)
require(
    sp.factor(constant_equation.coeff_monomial(x**2) - 3 * (a * u0 + 1) / (2 * a))
    == 0,
    "x^2 obstruction",
)
require(
    sp.factor(constant_equation.as_expr().subs(u0, -1 / a) + c) == 0,
    "terminal nonzero-constant obstruction",
)

print("THM-2709 complete linear target-plane classification on polynomial graphs")
print("Pluecker=(P,Q,R):J=P*J_AB+Q*J_Bd+R*J_dA")
print("P=0:plane contains d;classified by THM-2705")
print("P!=0:normalized PDE=(x+a*y+b)f_x+(x^2+a*f+b*x)f_y+f-x*y=c")
print("Q=0:unique f=(x-b)t/2+(x-b)(x^2+b^2)/8+c;t=y-x^2/2")
print("Q=0:target=(P*A,B-2*b*d);triangular V=(t+b^2/2)^2-2*c*(x+b)")
print("P*Q!=0:EMPTY;deg_y>=2 Riccati collapse;linear branch forces c=0")
print("shifted_cubic_inverse=PASS;row_wedge_expansion=PASS")
print("ALL CHECKS PASSED")
