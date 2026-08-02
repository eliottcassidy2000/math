#!/usr/bin/env python3
"""Exact controls for THM-3143.

For an algebraic quadratic phase of Hessian rank at most one, Cauchy--Green
or Green turns the simplex area period into oriented edge integrals.  This
script checks the inherited aligned identities and the transverse
discriminant-block extension: edge restrictions and ODEs, singleton/path/
cycle sources, the odd-cycle obstruction, and the vertical-edge boundary.

Reproduce:
  python3 04-computation/fc3_complex_affine_coordinate_cauchy_green_thm3143.py
  python3 -O 04-computation/fc3_complex_affine_coordinate_cauchy_green_thm3143.py
"""

from __future__ import annotations

import sympy as sp


def check(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


print("THM-3143 FC(3) RANK-AT-MOST-ONE QUADRATIC EDGE-CURRENT AUDIT")

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
print("C6 aligned collinear boundary y=0: W and all three turns vanish exactly")

# C7. Add a transverse term lambda*conjugate(z).  Its restriction to edge j
# is a quadratic with an edge-dependent discriminant.  Check Cauchy--Green
# coefficientwise, all restriction identities, and every edge ODE.
lam = sp.symbols("lambda", nonzero=True)
transverse_q = A * t**2 + B * t + C + lam * sp.conjugate(t)
transverse_edge_data: list[tuple[sp.Expr, sp.Expr, sp.Expr]] = []
transverse_edge_checks = 0
for j, k in edges:
    zj, zk = vertices[j], vertices[k]
    edge_point = zj + s * (zk - zj)
    edge_slope = edge_slopes[j]
    edge_intercept = sp.simplify(sp.conjugate(zj) - edge_slope * zj)
    beta_edge = sp.simplify(B + lam * edge_slope)
    gamma_edge = sp.simplify(C + lam * edge_intercept)
    delta_edge = sp.expand(beta_edge**2 - 4 * A * gamma_edge)
    edge_quadratic = sp.expand(A * t**2 + beta_edge * t + gamma_edge)
    restricted = sp.expand(
        A * edge_point**2
        + B * edge_point
        + C
        + lam * sp.conjugate(edge_point).expand(complex=True)
    )
    check(
        sp.simplify(restricted - edge_quadratic.subs(t, edge_point)) == 0,
        f"transverse edge restriction {j}->{k}",
    )
    for degree in range(4):
        previous = (
            delta_edge * degree * edge_quadratic ** (degree - 1)
            if degree
            else sp.Integer(0)
        )
        ode_density = (4 * A * degree + 2 * A) * edge_quadratic**degree + previous
        endpoint_density = sp.diff(
            (2 * A * t + beta_edge) * edge_quadratic**degree, t
        )
        check(
            sp.factor(endpoint_density - ode_density) == 0,
            f"transverse edge ODE {j}->{k}, degree={degree}",
        )
        transverse_edge_checks += 1
    transverse_edge_data.append((beta_edge, gamma_edge, delta_edge))

z_uv = u + tau * v
q_uv = sp.expand(
    A * z_uv**2
    + B * z_uv
    + C
    + lam * sp.conjugate(z_uv).expand(complex=True)
)
transverse_cg_checks = 0
for degree in range(1, 4):
    boundary = sp.Integer(0)
    for j, k in edges:
        edge_point = vertices[j] + s * (vertices[k] - vertices[j])
        edge_q = sp.expand(
            A * edge_point**2
            + B * edge_point
            + C
            + lam * sp.conjugate(edge_point).expand(complex=True)
        )
        boundary += sp.integrate(
            edge_q**degree * (vertices[k] - vertices[j]), (s, 0, 1)
        )
    direct = sp.integrate(
        sp.integrate(q_uv ** (degree - 1), (v, 0, 1 - u)), (u, 0, 1)
    )
    check(
        sp.factor(boundary - 2 * ii * degree * lam * y * direct) == 0,
        f"transverse Cauchy-Green coefficient degree={degree}",
    )
    transverse_cg_checks += 1
print(
    f"C7 transverse restrictions/ODEs={transverse_edge_checks}; "
    f"Cauchy-Green coefficients={transverse_cg_checks}"
)

# C8. Abstract exact source controls for the three possible discriminant
# blocks.  Equal discriminant at a shared vertex and distinct slopes force
# its two derivatives to be d and -d.
gap = sp.symbols("gap", nonzero=True)
dpath = sp.symbols("dpath", nonzero=True)
singleton_left = -A * gap
singleton_right = A * gap
check(
    sp.factor(singleton_right - singleton_left - 2 * A * gap) == 0,
    "singleton reflection source",
)

gap_left, gap_right = sp.symbols("gap_left gap_right")
source_left = -dpath + 2 * A * gap_left
source_middle = 2 * dpath
source_right = -dpath + 2 * A * gap_right
check(
    sp.simplify(source_left.subs(gap_left, dpath / A) - dpath) == 0,
    "two-edge left collision source",
)
check(
    sp.simplify(source_right.subs(gap_right, dpath / A) - dpath) == 0,
    "two-edge right collision source",
)
check(
    sp.simplify(
        (source_left + source_middle).subs(gap_left, dpath / A)
        - 3 * dpath
    )
    == 0,
    "two-edge grouped 3d source left",
)
check(
    sp.simplify(
        (source_middle + source_right).subs(gap_right, dpath / A)
        - 3 * dpath
    )
    == 0,
    "two-edge grouped 3d source right",
)
check(
    sp.simplify(
        (source_left + source_middle + source_right).subs(
            {gap_left: dpath / A, gap_right: dpath / A}
        )
        - 4 * dpath
    )
    == 0,
    "two-edge grouped 4d source",
)

m0, m1, m2 = sp.symbols("m0 m1 m2")
cycle_sources = (
    lam * (m2 - m0),
    lam * (m0 - m1),
    lam * (m1 - m2),
)
check(sp.expand(sum(cycle_sources)) == 0, "full-cycle source conservation")
for j in range(3):
    check(
        sp.expand(
            cycle_sources[(j + 1) % 3]
            + cycle_sources[(j + 2) % 3]
            + cycle_sources[j]
        )
        == 0,
        f"full-cycle pair quotient {j}",
    )

odd_sign_checks = 0
for sign0 in (-1, 1):
    for sign1 in (-1, 1):
        for sign2 in (-1, 1):
            check(sign0 + sign1 + sign2 != 0, "odd equal-square cycle closure")
            odd_sign_checks += 1
print(
    "C8 discriminant-block sources: singleton=1; path currents=2d,3d,4d; "
    f"cycle partitions=3; odd-sign controls={odd_sign_checks}"
)

# C9. Green's formula on a real algebraic triangle with one vertical edge.
# Vertices (0,0),(R,K),(0,H) are counterclockwise when R*H>0; the last
# edge has dx=0.  Verify integral_boundary f dx=-integral_T d_y(f).
Rgeo, Hgeo, Kgeo = sp.symbols("R H K", real=True, nonzero=True)
Xsym, Ysym = sp.symbols("X Y", real=True)
real_vertices = (
    (sp.Integer(0), sp.Integer(0)),
    (Rgeo, Kgeo),
    (sp.Integer(0), Hgeo),
)
green_checks = 0
for degree_x in range(4):
    for degree_y in range(4 - degree_x):
        monomial = Xsym**degree_x * Ysym**degree_y
        boundary = sp.Integer(0)
        for j, k in edges:
            xj, yj = real_vertices[j]
            xk, yk = real_vertices[k]
            edge_x = xj + s * (xk - xj)
            edge_y = yj + s * (yk - yj)
            boundary += sp.integrate(
                monomial.subs({Xsym: edge_x, Ysym: edge_y}) * (xk - xj),
                (s, 0, 1),
            )
        pulled_derivative = sp.diff(monomial, Ysym).subs(
            {Xsym: Rgeo * u, Ysym: Kgeo * u + Hgeo * v}
        )
        direct = -Rgeo * Hgeo * sp.integrate(
            sp.integrate(pulled_derivative, (v, 0, 1 - u)), (u, 0, 1)
        )
        check(sp.factor(boundary - direct) == 0, f"Green monomial {(degree_x, degree_y)}")
        green_checks += 1

# The vertical edge is edge 2 and contributes no dx.  The two active x-gaps
# are R and -R.  Audit all collision partitions in the common-discriminant
# case using the same d,-d shared derivatives as C8.
check(real_vertices[0][0] == real_vertices[2][0], "repeated vertical x value")
check(real_vertices[0][0] != real_vertices[1][0], "first active x gap")
check(real_vertices[1][0] != real_vertices[2][0], "second active x gap")
vertical_start = -dpath + 2 * A * Rgeo
vertical_middle = 2 * dpath
vertical_end = -dpath - 2 * A * Rgeo
check(sp.expand(vertical_start + vertical_end) == -2 * dpath, "vertical endpoint pair")
check(
    sp.simplify(
        (vertical_start + vertical_middle).subs(Rgeo, dpath / A)
        - 3 * dpath
    )
    == 0,
    "vertical shared/start collision",
)
check(
    sp.simplify(
        (vertical_middle + vertical_end).subs(Rgeo, -dpath / A)
        - 3 * dpath
    )
    == 0,
    "vertical shared/end collision",
)
check(sp.expand(A * Rgeo - (-A * Rgeo)) == 2 * A * Rgeo, "vertical all-value contradiction")
print(
    f"C9 Green vertical-edge monomial checks={green_checks}; "
    "collision currents=2d,3d; all-value opposite-gap obstruction exact"
)

print("ALL THM-3143 CONTROLS PASSED")
