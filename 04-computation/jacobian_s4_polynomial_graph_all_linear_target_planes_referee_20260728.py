#!/usr/bin/env python3
"""Exact hostile companion for all linear target planes on THM-2696 graphs.

This is an unnumbered proof-candidate companion.  It completes the target
planes excluded by THM-2705: planes whose linear span does not contain d.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


x, y, t = sp.symbols("x y t")
a, b, c = sp.symbols("a b c", nonzero=True)
q0 = sp.symbols("q0")
UU, VV = sp.symbols("UU VV")
fxy = sp.Function("f")(x, y)

A = x**2 - 2 * y
B = y**2 - 2 * x * fxy
d = fxy

# A target two-plane not containing d has projective normal (a,b,1).
# The canonical row basis (U,V)=(A-a*d,B-b*d) has precisely this normal.
rU = sp.Matrix([1, 0, -a])
rV = sp.Matrix([0, 1, -b])
require(rU.cross(rV) == sp.Matrix([a, b, 1]), "canonical target-plane normal")

U = A - a * d
V = B - b * d
J_UV = sp.factor(sp.det(sp.Matrix([U, V]).jacobian([x, y])))
fx = sp.diff(fxy, x)
fy = sp.diff(fxy, y)
bracket = sp.expand(
    (a * y + b + 2 * x) * fx
    + (a * fxy + b * x + 2 * x**2) * fy
    + 2 * fxy
    - 2 * x * y
)
require(sp.expand(J_UV + 2 * bracket) == 0, "normal-form Jacobian equation")

# If a != 0, y-degree r>=3 has unique leading term a*r*p_r^2*y^(2r-1).
for r in range(3, 40):
    require(2 * r - 1 > r + 1, f"unique high y-degree at r={r}")

# At the only tied boundary r=2, the y^3 coefficient is the polynomial
# Riccati equation a*(p2'+2*p2^2)=0.
p2 = sp.Function("p2")(x)
p1 = sp.Function("p1")(x)
p0 = sp.Function("p0")(x)
f_quad = p2 * y**2 + p1 * y + p0
quad_equation = sp.Poly(
    sp.expand(
        (a * y + b + 2 * x) * sp.diff(f_quad, x)
        + (a * f_quad + b * x + 2 * x**2) * sp.diff(f_quad, y)
        + 2 * f_quad
        - 2 * x * y
        - c
    ),
    y,
)
require(
    sp.expand(
        quad_equation.coeff_monomial(y**3)
        - a * (sp.diff(p2, x) + 2 * p2**2)
    )
    == 0,
    "quadratic y-degree Riccati coefficient",
)
for polynomial_degree in range(1, 40):
    require(
        2 * polynomial_degree > polynomial_degree - 1,
        f"polynomial Riccati degree obstruction at degree={polynomial_degree}",
    )

# Therefore f=p(x)y+q(x).  The three coefficients below are exhaustive and
# include the r=0 boundary by allowing p=0.
p = sp.Function("p")(x)
q = sp.Function("q")(x)
f_linear = p * y + q
linear_equation = sp.Poly(
    sp.expand(
        (a * y + b + 2 * x) * sp.diff(f_linear, x)
        + (a * f_linear + b * x + 2 * x**2) * sp.diff(f_linear, y)
        + 2 * f_linear
        - 2 * x * y
        - c
    ),
    y,
)
linear_expected = {
    2: a * sp.diff(p, x),
    1: a * p**2 + a * sp.diff(q, x) + (b + 2 * x) * sp.diff(p, x) - 2 * x + 2 * p,
    0: a * p * q + (b * x + 2 * x**2) * p + (b + 2 * x) * sp.diff(q, x) + 2 * q - c,
}
for degree, expected in linear_expected.items():
    require(
        sp.expand(linear_equation.coeff_monomial(y**degree) - expected) == 0,
        f"linear y-degree coefficient y^{degree}",
    )

# Since a!=0, p is constant and q'=2x/a-2p/a-p^2.  Substitution leaves an
# x^2 coefficient 3*(a*p+2)/a, hence p=-2/a.  At that value the full left
# side is zero, so the equation can hold only on the sharp c=0 boundary.
p_const = sp.symbols("p_const")
q_general = x**2 / a - (2 * p_const / a + p_const**2) * x + q0
constant_row = sp.Poly(
    sp.expand(
        a * p_const * q_general
        + (b * x + 2 * x**2) * p_const
        + (b + 2 * x) * sp.diff(q_general, x)
        + 2 * q_general
        - c
    ),
    x,
)
require(
    sp.factor(constant_row.coeff_monomial(x**2) - 3 * (a * p_const + 2) / a)
    == 0,
    "forced linear slope",
)
require(
    sp.factor(constant_row.as_expr().subs(p_const, -2 / a) + c) == 0,
    "constant collapse after p=-2/a",
)

# The c=0 equality boundary is not a missed Keller family: U is constant and
# the target map has rank at most one.
t_xy = y - x**2 / 2
f_rank_drop = -2 * t_xy / a + q0
U_rank_drop = sp.expand((A - a * fxy).subs(fxy, f_rank_drop))
V_rank_drop = sp.expand((B - b * fxy).subs(fxy, f_rank_drop))
require(sp.expand(U_rank_drop + a * q0) == 0, "zero-Jacobian boundary has constant U")
require(
    sp.factor(sp.det(sp.Matrix([U_rank_drop, V_rank_drop]).jacobian([x, y]))) == 0,
    "zero-Jacobian boundary rank drop",
)

# The remaining non-d planes have a=0, projective normal (0,b,1), and the
# canonical rows (A,B-b*d).  In t=y-x^2/2 the equation is
#   (2x+b) f_x + 2f = 2xt+x^3+c.
# After u=x+b/2 the operator is 2*(u*d_u+1), hence is invertible on every
# polynomial x-degree.  The displayed f is its unique polynomial solution.
f_survivor_xt = sp.expand(
    (x - b / 2) * t / 2
    + (x - b / 2) * (x**2 + b**2 / 4) / 8
    + c / 2
)
ode_residual = sp.expand(
    (2 * x + b) * sp.diff(f_survivor_xt, x)
    + 2 * f_survivor_xt
    - (2 * x * t + x**3 + c)
)
require(ode_residual == 0, "unique a=0 polynomial ODE solution")
for x_degree in range(40):
    require(2 * (x_degree + 1) != 0, f"Euler eigenvalue at x-degree={x_degree}")

f_survivor = sp.expand(f_survivor_xt.subs(t, t_xy))
U_survivor = A
V_survivor = sp.expand((B - b * fxy).subs(fxy, f_survivor))
J_survivor = sp.factor(
    sp.det(sp.Matrix([U_survivor, V_survivor]).jacobian([x, y]))
)
require(sp.expand(J_survivor + 2 * c) == 0, "surviving constant Jacobian")

# The target map is explicitly triangular:
#   U=-2t, V=-c*x+(t+b^2/8)^2-bc/2.
G = sp.expand((t + b**2 / 8) ** 2 - b * c / 2)
require(
    sp.expand(V_survivor.subs(y, t + x**2 / 2) - (-c * x + G)) == 0,
    "surviving triangular target",
)
t_inv = -UU / 2
x_inv = sp.factor((G.subs(t, t_inv) - VV) / c)
y_inv = sp.factor(t_inv + x_inv**2 / 2)
require(
    sp.factor(U_survivor.subs({x: x_inv, y: y_inv}) - UU) == 0,
    "surviving inverse U",
)
require(
    sp.factor(V_survivor.subs({x: x_inv, y: y_inv}) - VV) == 0,
    "surviving inverse V",
)

# Hostile controls: b=0 recovers exactly THM-2702's unique (A,B) graph; a
# concrete shifted member retains the asserted constant and inverse.
kappa = sp.symbols("kappa", nonzero=True)
thm2702_expected = x * t / 2 + x**3 / 8 - kappa / 4
require(
    sp.expand(f_survivor_xt.subs({b: 0, c: -kappa / 2}) - thm2702_expected) == 0,
    "THM-2702 b=0 control",
)
require(
    sp.expand(J_survivor.subs({b: 3, c: 5}) + 10) == 0,
    "shifted nonlinear hostile constant",
)

print("Polynomial graphs: completion of all linear target planes")
print("projective_normal:nu=0 -> THM-2705; nu!=0 -> normalize (a,b,1)")
print("canonical_non-d_pair=(A-a*d,B-b*d)")
print("Jac=-2*((a*y+b+2*x)f_x+(a*f+b*x+2*x^2)f_y+2*f-2*x*y)")
print("a!=0:y-degree r>=3 unique top; r=2 Riccati; r<=1 forces p=-2/a")
print("a!=0:equation collapses to c=0; therefore kappa=-2c!=0 is EMPTY")
print("a!=0,c=0:f=-2*(y-x^2/2)/a+q0 and U=-a*q0 (rank drop)")
print("a=0:unique f=(x-b/2)t/2+(x-b/2)(x^2+b^2/4)/8+c/2")
print("a=0:U=-2t;V=-c*x+(t+b^2/8)^2-b*c/2; explicit inverse PASS")
print("controls=THM2702_b0:PASS;shifted_b3_c5:PASS;zero_boundary:PASS")
print("ALL CHECKS PASSED")
