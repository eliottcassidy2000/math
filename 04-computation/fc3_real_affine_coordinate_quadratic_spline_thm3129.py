#!/usr/bin/env python3
"""Exact controls for THM-3129.

The theorem treats q=Q(r), where r is a real algebraic affine coordinate on
the two-simplex and Q has degree at most two.  When the three vertex values
of r are distinct, its pushforward density is a two-piece linear B-spline.
This script checks that density, the corrected-period decomposition, the
quadratic segment ODE, the signed square source vector, every possible
quadratic endpoint collision, and a sharp transcendental cancellation
hostile.

Reproduce:
  python3 04-computation/fc3_real_affine_coordinate_quadratic_spline_thm3129.py
  python3 -O 04-computation/fc3_real_affine_coordinate_quadratic_spline_thm3129.py
"""

from __future__ import annotations

import sympy as sp


def check(condition: bool, message: str) -> None:
    """Optimization-stable exact check."""
    if not condition:
        raise RuntimeError(message)


def simplex_integral(poly: sp.Expr, u: sp.Symbol, v: sp.Symbol) -> sp.Expr:
    return sp.simplify(sp.integrate(sp.integrate(poly, (v, 0, 1 - u)), (u, 0, 1)))


print("THM-3129 FC(3) REAL AFFINE-COORDINATE QUADRATIC SPLINE AUDIT")

t, u, v = sp.symbols("t u v", real=True)
a = sp.symbols("a", real=True)
p, q = sp.symbols("p q", positive=True)
b = a + p
c = a + p + q
span = p + q

# C1. Exact pushforward density.  The coordinate r has barycentric vertex
# values a<b<c, represented here by r=a*w+b*u+c*v.
rho_left = (t - a) / (p * span)
rho_right = (c - t) / (q * span)
left_mass = sp.integrate(rho_left, (t, a, b))
right_mass = sp.integrate(rho_right, (t, b, c))
check(sp.simplify(left_mass + right_mass - sp.Rational(1, 2)) == 0, "spline mass")

w = 1 - u - v
r_affine = sp.expand(a * w + b * u + c * v)
spline_moment_checks = 0
for m in range(7):
    direct = simplex_integral(r_affine**m, u, v)
    pushed = sp.integrate(rho_left * t**m, (t, a, b)) + sp.integrate(
        rho_right * t**m, (t, b, c)
    )
    check(sp.simplify(direct - pushed) == 0, f"spline moment m={m}")
    spline_moment_checks += 1
print(
    "C1 two-piece spline mass=1/2; "
    f"direct/pushforward moments={spline_moment_checks}"
)

# C2. The quadratic segment system.  For Q(t)=A*t^2+B*t+C and
# H_xy(z)=int_x^y exp(zQ(t))dt, the exact ODE is
#
# 4Az H'+[2A+(B^2-4AC)z]H=Q'(y)e^(zQ(y))-Q'(x)e^(zQ(x)).
# Verify its coefficient recurrence on both spline intervals.
A = sp.symbols("A", nonzero=True)
B = sp.symbols("B")
C = sp.symbols("C")
Q = sp.expand(A * t**2 + B * t + C)


def qvalue(x: sp.Expr) -> sp.Expr:
    return sp.expand(Q.subs(t, x))


def qderiv(x: sp.Expr) -> sp.Expr:
    return sp.expand(2 * A * x + B)


segment_ode_checks = 0
for x0, x1 in ((a, b), (b, c)):
    moments = [sp.integrate(sp.expand(Q**m), (t, x0, x1)) for m in range(7)]
    for m in range(7):
        previous = moments[m - 1] if m else sp.Integer(0)
        lhs = (
            (4 * A * m + 2 * A) * moments[m]
            + (B**2 - 4 * A * C) * m * previous
        )
        rhs = qderiv(x1) * qvalue(x1) ** m - qderiv(x0) * qvalue(x0) ** m
        check(sp.expand(lhs - rhs) == 0, f"segment ODE {(x0, x1, m)}")
        segment_ode_checks += 1
print(f"C2 quadratic segment ODE coefficient recurrences={segment_ode_checks}")

# C3. Subtract the explicit endpoint term from span*K.  The remaining
# corrected period M=alpha*H_ab+beta*H_bc has a signed source vector under
# the common first-order operator.  It simplifies to weighted squares.
h = -B / (2 * A)
alpha = sp.simplify((h - a) / p)
beta = sp.simplify((c - h) / q)
d_a, d_b, d_c = (qderiv(x) for x in (a, b, c))
source_a = sp.factor(-alpha * d_a)
source_b = sp.factor((alpha - beta) * d_b)
source_c = sp.factor(beta * d_c)

expected_a = d_a**2 / (2 * A * p)
expected_b = -span * d_b**2 / (2 * A * p * q)
expected_c = d_c**2 / (2 * A * q)
for label, actual, expected in (
    ("a", source_a, expected_a),
    ("b", source_b, expected_b),
    ("c", source_c, expected_c),
):
    check(sp.factor(actual - expected) == 0, f"square source {label}")

for x, dx in ((a, d_a), (b, d_b), (c, d_c)):
    kappa_numerator = sp.expand(B**2 - 4 * A * C + 4 * A * qvalue(x))
    check(sp.factor(kappa_numerator - dx**2) == 0, f"kappa square at {x}")
print("C3 corrected-period sources=(d_a^2/(2Ap), -span*d_b^2/(2Apq), d_c^2/(2Aq))")

# C4. A nonconstant quadratic can identify at most one pair of knot values.
# In each possible reflection collision, combine the two source coefficients.
# None cancels.
collision_data = (
    ("Q(a)=Q(b)", -A * (a + b), source_a + source_b, -A * p**2 / (2 * q)),
    ("Q(b)=Q(c)", -A * (b + c), source_b + source_c, -A * q**2 / (2 * p)),
    (
        "Q(a)=Q(c)",
        -A * (a + c),
        source_a + source_c,
        A * span**3 / (2 * p * q),
    ),
)
collision_checks = 0
for label, collision_B, combined, expected in collision_data:
    check(
        sp.factor(combined.subs(B, collision_B) - expected) == 0,
        f"collision source {label}",
    )
    collision_checks += 1
print(f"C4 reflection-collision combined sources nonzero: cases={collision_checks}")

# C5. The rational-independence obstruction is the same half-integral pole
# equation as in THM-3116.  Integral pole orders never hit the resonance 1/2.
for pole_order in range(1, 13):
    check(
        -pole_order + sp.Rational(1, 2) != 0,
        f"half-integral pole obstruction {pole_order}",
    )
print("C5 rational coefficient equation: pole multipliers -n+1/2 nonzero for n=1..12")

# C6. Sharp hostile outside the algebraic-coefficient scope.  For knot values
# 0,1,2 the spline transform is (exp(g)-1)^2/(2g^2), which vanishes at
# g=2*pi*i even though the representing measure is positive.
g = sp.symbols("g", nonzero=True)
hostile_transform = sp.integrate(t * sp.exp(g * t) / 2, (t, 0, 1)) + sp.integrate(
    (2 - t) * sp.exp(g * t) / 2, (t, 1, 2)
)
hostile_formula = (sp.exp(g) - 1) ** 2 / (2 * g**2)
check(sp.simplify(hostile_transform - hostile_formula) == 0, "hostile transform formula")
hostile_value = sp.simplify(hostile_formula.subs(g, 2 * sp.pi * sp.I))
check(hostile_value == 0, "transcendental spline cancellation")
print("C6 hostile knots=(0,1,2): transform=(e^g-1)^2/(2g^2), value at g=2*pi*i is 0")

print("ALL THM-3129 CONTROLS PASSED")
