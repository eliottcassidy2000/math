#!/usr/bin/env python3
"""Exact controls for THM-3143.

For a noncollinear algebraic affine coordinate, Cauchy--Green turns the
simplex area period into three oriented complex-edge integrals.  This script
checks the area/boundary orientation, conjugate-line equations, endpoint
correction, closed-cycle relation, common quadratic ODE, the locally derived
turning-square source vector, all endpoint collisions, and the exact
collinear failure boundary.

Reproduce:
  python3 04-computation/fc3_complex_affine_coordinate_cauchy_green_thm3143.py
  python3 -O 04-computation/fc3_complex_affine_coordinate_cauchy_green_thm3143.py
"""

from __future__ import annotations

import sympy as sp


def check(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


print("THM-3143 FC(3) COMPLEX AFFINE-COORDINATE CAUCHY-GREEN AUDIT")

# Normalize an oriented noncollinear algebraic triangle to vertices
# z0=0,z1=1,z2=tau=x+i*y, with real y>0.
u, v, s = sp.symbols("u v s", real=True)
x = sp.symbols("x", real=True)
y = sp.symbols("y", positive=True)
ii = sp.I
tau = x + ii * y
vertices = (sp.Integer(0), sp.Integer(1), tau)
edges = ((0, 1), (1, 2), (2, 0))
W = sp.simplify(tau - sp.conjugate(tau))
check(W == 2 * ii * y, "oriented coordinate Jacobian")

# C1. Cauchy--Green: W int_Delta f(u+tau*v)dudv
# equals int_boundary conjugate(t) f(t)dt for holomorphic monomials.
cauchy_green_checks = 0
for degree in range(7):
    direct = sp.integrate(
        sp.integrate((u + tau * v) ** degree, (v, 0, 1 - u)), (u, 0, 1)
    )
    boundary = sp.Integer(0)
    for j, k in edges:
        edge_point = vertices[j] + s * (vertices[k] - vertices[j])
        edge_derivative = vertices[k] - vertices[j]
        boundary += sp.integrate(
            sp.conjugate(edge_point).expand(complex=True)
            * edge_point**degree
            * edge_derivative,
            (s, 0, 1),
        )
    check(sp.factor(boundary - W * direct) == 0, f"Cauchy-Green degree={degree}")
    cauchy_green_checks += 1
print(f"C1 Cauchy-Green oriented monomial checks={cauchy_green_checks}; W=2*i*y")

# C2. Each edge has conjugate(t)=m*t+n with m=conjugate(e)/e.  Check the
# line identity and the weighted quadratic endpoint correction coefficientwise.
A = sp.symbols("A", nonzero=True)
B = sp.symbols("B")
C = sp.symbols("C")
t = sp.symbols("t")
Q = A * t**2 + B * t + C


def qvalue(z: sp.Expr) -> sp.Expr:
    return sp.expand(Q.subs(t, z))


edge_slopes: list[sp.Expr] = []
edge_checks = 0
for j, k in edges:
    zj, zk = vertices[j], vertices[k]
    e = zk - zj
    m_edge = sp.simplify(sp.conjugate(e) / e)
    n_edge = sp.simplify(sp.conjugate(zj) - m_edge * zj)
    gamma_edge = sp.simplify(n_edge - m_edge * B / (2 * A))
    edge_slopes.append(m_edge)

    edge_point = zj + s * e
    check(
        sp.simplify(
            sp.conjugate(edge_point).expand(complex=True)
            - (m_edge * edge_point + n_edge)
        )
        == 0,
        f"conjugate edge line {j}->{k}",
    )

    # Pointwise identity behind every coefficient of the endpoint correction:
    # conjugate(t)-gamma=m*Q'(t)/(2A) along this edge.
    check(
        sp.factor(
            sp.conjugate(edge_point).expand(complex=True)
            - gamma_edge
            - m_edge * (2 * A * edge_point + B) / (2 * A)
        )
        == 0,
        f"pointwise edge correction {j}->{k}",
    )
    edge_checks += 2
print(f"C2 conjugate-line and endpoint-correction checks={edge_checks}")

# C3. The three edge periods are not independent: every coefficient integral
# is the integral of an entire polynomial around a closed triangle.  Record
# that relation before checking the common segment ODE and its cycle source.
cycle_checks = 0
for degree in range(6):
    primitive = sp.integrate(Q**degree, t)
    cycle_coefficient = sum(
        primitive.subs(t, vertices[k]) - primitive.subs(t, vertices[j])
        for j, k in edges
    )
    check(
        sp.simplify(cycle_coefficient) == 0,
        f"closed-cycle coefficient degree={degree}",
    )
    cycle_checks += 1

# Pin the incoming shifted edge constant at the same vertex before subtracting
# the outgoing one.  Only after that subtraction is the turning-square form a
# valid local corollary.
segment_ode_checks = 0
for degree in range(6):
    previous_term = (
        (B**2 - 4 * A * C) * degree * Q ** (degree - 1)
        if degree
        else sp.Integer(0)
    )
    ode_density = (4 * A * degree + 2 * A) * Q**degree + previous_term
    endpoint_derivative = sp.diff((2 * A * t + B) * Q**degree, t)
    check(
        sp.factor(endpoint_derivative - ode_density) == 0,
        f"complex segment ODE density degree={degree}",
    )
    segment_ode_checks += 1

h = -B / (2 * A)
sources: list[sp.Expr] = []
turns: list[sp.Expr] = []
for j in range(3):
    incoming = (j - 1) % 3
    outgoing = j
    zj = vertices[j]
    d_j = 2 * A * zj + B
    gamma_in_at_j = sp.simplify(
        sp.conjugate(zj) - edge_slopes[incoming] * (zj - h)
    )
    gamma_out_at_j = sp.simplify(
        sp.conjugate(zj) - edge_slopes[outgoing] * (zj - h)
    )
    direct_source = sp.simplify(d_j * (gamma_in_at_j - gamma_out_at_j))
    turn = sp.simplify(edge_slopes[outgoing] - edge_slopes[incoming])
    square_source = sp.simplify(d_j**2 * turn / (2 * A))
    check(
        sp.simplify(direct_source - square_source) == 0,
        f"turning source vertex={j}",
    )
    turns.append(turn)
    sources.append(square_source)
print(
    f"C3 closed-cycle coefficient relations={cycle_checks}; "
    f"complex segment ODE recurrences={segment_ode_checks}; "
    "locally derived turning-square sources=3"
)

# C4. The normalized turning factors are all nonzero multiples of y.  In a
# quadratic reflection collision the shared edge cancels, leaving the two
# outer edge directions; verify all three cases exactly.
expected_turns = (
    2 * ii * y / tau,
    -2 * ii * y / (tau - 1),
    2 * ii * y / (tau * (tau - 1)),
)
for j, expected in enumerate(expected_turns):
    check(sp.simplify(turns[j] - expected) == 0, f"turn factor vertex={j}")

collision_checks = 0
for i, j in ((0, 1), (1, 2), (0, 2)):
    collision_B = -A * (vertices[i] + vertices[j])
    combined = sp.factor((sources[i] + sources[j]).subs(B, collision_B))
    # Under collision the squared derivatives agree.  Divide by that
    # nonzero square factor and compare with the sum of turns, which is the
    # difference of the two outer edge slopes.
    derivative_square = sp.expand(
        (2 * A * vertices[i] + collision_B) ** 2
    )
    check(derivative_square != 0, f"collision derivative pair {(i, j)}")
    check(
        sp.simplify(
            combined - derivative_square * (turns[i] + turns[j]) / (2 * A)
        )
        == 0,
        f"collision outer-edge source {(i, j)}",
    )
    check(
        sp.simplify(turns[i] + turns[j]) != 0,
        f"collision outer-edge slopes remain distinct {(i, j)}",
    )
    check(sp.simplify(combined) != 0, f"collision grouped source {(i, j)}")
    collision_checks += 1
print(f"C4 noncollinear turns and reflection-collision sources: cases={collision_checks}")

# C5. The functional-independence rational equation again has half-integral
# pole coefficient -n+1/2.  No integral pole order resonates.
for pole_order in range(1, 13):
    check(
        -pole_order + sp.Rational(1, 2) != 0,
        f"rational pole order {pole_order}",
    )
print("C5 rational coefficient pole multipliers -n+1/2 nonzero for n=1..12")

# C6. Geometry failure is exact: every turn factor is proportional to y and
# W=2iy.  At y=0 the triangle, Cauchy--Green area normalization, and cycle
# source all collapse together; THM-3142's one-dimensional spline takes over.
check(W.subs(y, 0) == 0, "collinear Jacobian collapse")
for j, turn in enumerate(turns):
    check(sp.simplify(turn.subs(y, 0)) == 0, f"collinear turn collapse {j}")
print("C6 collinear boundary y=0: W and all three turning factors vanish exactly")

print("ALL THM-3143 CONTROLS PASSED")
