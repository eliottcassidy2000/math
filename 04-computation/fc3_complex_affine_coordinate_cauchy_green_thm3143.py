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

import random

import sympy as sp


def check(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


print("THM-3143 FC(3) ALL ALGEBRAIC QUADRATIC EDGE-CURRENT AUDIT")

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

# C10. Full-rank Euler flux and normalized edge-class audit.  These helpers
# use exact algebraic arithmetic only.
def det2(left: tuple[sp.Expr, sp.Expr], right: tuple[sp.Expr, sp.Expr]) -> sp.Expr:
    return sp.expand(left[0] * right[1] - left[1] * right[0])


def vec_sub(
    left: tuple[sp.Expr, sp.Expr], right: tuple[sp.Expr, sp.Expr]
) -> tuple[sp.Expr, sp.Expr]:
    return (sp.expand(left[0] - right[0]), sp.expand(left[1] - right[1]))


def quadratic_form(matrix: sp.Matrix, vector: tuple[sp.Expr, sp.Expr]) -> sp.Expr:
    column = sp.Matrix(vector)
    return sp.expand((column.T * matrix * column)[0])


def polar_form(
    matrix: sp.Matrix,
    left: tuple[sp.Expr, sp.Expr],
    right: tuple[sp.Expr, sp.Expr],
) -> sp.Expr:
    return sp.expand((sp.Matrix(left).T * matrix * sp.Matrix(right))[0])


def same(left: sp.Expr, right: sp.Expr = sp.Integer(0)) -> bool:
    return sp.simplify(left - right) == 0


def add_group(
    groups: list[list[sp.Expr]], key: sp.Expr, coefficient: sp.Expr
) -> None:
    for entry in groups:
        if same(entry[0], key):
            entry[1] = sp.simplify(entry[1] + coefficient)
            return
    groups.append([sp.simplify(key), sp.simplify(coefficient)])


def simplex_polynomial_integral(expression: sp.Expr, jacobian: sp.Expr) -> sp.Expr:
    polynomial = sp.Poly(sp.expand(expression), u, v)
    total = sp.Integer(0)
    for (degree_u, degree_v), coefficient in polynomial.terms():
        total += (
            coefficient
            * sp.factorial(degree_u)
            * sp.factorial(degree_v)
            / sp.factorial(degree_u + degree_v + 2)
        )
    return sp.simplify(jacobian * total)


def analyze_full_rank_case(
    case_name: str,
    raw_vertices: tuple[
        tuple[sp.Expr, sp.Expr],
        tuple[sp.Expr, sp.Expr],
        tuple[sp.Expr, sp.Expr],
    ],
    matrix: sp.Matrix,
    center: tuple[sp.Expr, sp.Expr],
    critical_value: sp.Expr,
) -> dict[str, int]:
    determinant = sp.simplify(matrix.det())
    check(determinant != 0, f"{case_name}: full-rank matrix")

    vertices_real = list(raw_vertices)
    first_edge = vec_sub(vertices_real[1], vertices_real[0])
    second_edge = vec_sub(vertices_real[2], vertices_real[0])
    jacobian = sp.simplify(det2(first_edge, second_edge))
    check(jacobian != 0, f"{case_name}: nondegenerate triangle")
    if bool(jacobian < 0):
        vertices_real[1], vertices_real[2] = vertices_real[2], vertices_real[1]
        first_edge = vec_sub(vertices_real[1], vertices_real[0])
        second_edge = vec_sub(vertices_real[2], vertices_real[0])
        jacobian = sp.simplify(det2(first_edge, second_edge))
    check(jacobian > 0, f"{case_name}: positive orientation")

    pulled_point = (
        sp.expand(vertices_real[0][0] + u * first_edge[0] + v * second_edge[0]),
        sp.expand(vertices_real[0][1] + u * first_edge[1] + v * second_edge[1]),
    )
    pulled_centered = vec_sub(pulled_point, center)
    pulled_q = sp.expand(critical_value + quadratic_form(matrix, pulled_centered))

    edge_records: list[dict[str, sp.Expr | int]] = []
    affine_value_groups: list[list[sp.Expr]] = []
    quadratic_edges = 0
    affine_nonconstant = 0
    affine_constant = 0
    resonant_edges = 0

    for edge_index in range(3):
        start_index = edge_index
        end_index = (edge_index + 1) % 3
        start = vertices_real[start_index]
        end = vertices_real[end_index]
        edge = vec_sub(end, start)
        centered_start = vec_sub(start, center)
        edge_a = sp.simplify(quadratic_form(matrix, edge))
        edge_b = sp.simplify(2 * polar_form(matrix, centered_start, edge))
        edge_c = sp.simplify(
            critical_value + quadratic_form(matrix, centered_start)
        )
        end_value = sp.simplify(
            critical_value + quadratic_form(matrix, vec_sub(end, center))
        )
        edge_weight = sp.simplify(det2(centered_start, edge) / 2)
        edge_polynomial = sp.expand(edge_a * t**2 + edge_b * t + edge_c)
        check(same(edge_polynomial.subs(t, 0), edge_c), f"{case_name}: edge start")
        check(same(edge_polynomial.subs(t, 1), end_value), f"{case_name}: edge end")

        # Binary Gram determinant and Euler flux weight.
        gram_left = sp.expand(
            quadratic_form(matrix, centered_start) * edge_a
            - polar_form(matrix, centered_start, edge) ** 2
        )
        gram_right = sp.expand(determinant * det2(centered_start, edge) ** 2)
        check(same(gram_left, gram_right), f"{case_name}: binary Gram edge={edge_index}")

        if same(edge_a):
            if same(edge_b):
                affine_constant += 1
                check(same(edge_weight), f"{case_name}: constant isotropic flux")
            else:
                affine_nonconstant += 1
                add_group(affine_value_groups, edge_c, -edge_weight / edge_b)
                add_group(affine_value_groups, end_value, edge_weight / edge_b)
            edge_records.append(
                {
                    "kind": 0,
                    "weight": edge_weight,
                    "start_value": edge_c,
                    "end_value": end_value,
                }
            )
            continue

        quadratic_edges += 1
        edge_delta = sp.simplify(edge_b**2 - 4 * edge_a * edge_c)
        edge_kappa = sp.simplify(edge_delta / (4 * edge_a))
        edge_mu = sp.simplify(edge_kappa + critical_value)
        check(
            same(
                edge_weight**2 / edge_a,
                -edge_mu / (4 * determinant),
            ),
            f"{case_name}: flux square edge={edge_index}",
        )
        check(
            same(edge_b**2 / (4 * edge_a), edge_kappa + edge_c),
            f"{case_name}: normalized start derivative edge={edge_index}",
        )
        check(
            same(
                (2 * edge_a + edge_b) ** 2 / (4 * edge_a),
                edge_kappa + end_value,
            ),
            f"{case_name}: normalized end derivative edge={edge_index}",
        )
        if same(edge_mu):
            resonant_edges += 1
            check(same(edge_weight), f"{case_name}: planted resonance flux")
        edge_records.append(
            {
                "kind": 1,
                "a": edge_a,
                "b": edge_b,
                "c": edge_c,
                "end_value": end_value,
                "weight": edge_weight,
                "kappa": edge_kappa,
                "mu": edge_mu,
                "polynomial": edge_polynomial,
            }
        )

    check(
        same(sum(sp.sympify(record["weight"]) for record in edge_records), jacobian / 2),
        f"{case_name}: zero-phase Euler flux",
    )

    # Euler-flux identity coefficientwise through cubic powers.
    previous_integral = sp.Integer(0)
    for degree in range(4):
        direct_integral = simplex_polynomial_integral(pulled_q**degree, jacobian)
        boundary_integral = sp.Integer(0)
        for edge_index, record in enumerate(edge_records):
            start = vertices_real[edge_index]
            end = vertices_real[(edge_index + 1) % 3]
            edge = vec_sub(end, start)
            centered_start = vec_sub(start, center)
            edge_q = sp.expand(
                critical_value
                + quadratic_form(
                    matrix,
                    (
                        centered_start[0] + t * edge[0],
                        centered_start[1] + t * edge[1],
                    ),
                )
            )
            boundary_integral += sp.sympify(record["weight"]) * sp.integrate(
                edge_q**degree, (t, 0, 1)
            )
        expected = (degree + 1) * direct_integral
        if degree:
            expected -= degree * critical_value * previous_integral
        check(
            same(boundary_integral, expected),
            f"{case_name}: Euler coefficient degree={degree}",
        )
        previous_integral = direct_integral

    # Group quadratic edges by kappa, normalize their flux coefficients, and
    # audit source vanishing/noncriticality exactly.
    kappa_classes: list[list[object]] = []
    for record in edge_records:
        if record["kind"] != 1:
            continue
        placed = False
        for entry in kappa_classes:
            if same(sp.sympify(entry[0]), sp.sympify(record["kappa"])):
                entry[1].append(record)
                placed = True
                break
        if not placed:
            kappa_classes.append([record["kappa"], [record]])

    nonzero_flux_blocks = 0
    zero_source_blocks = 0
    nonzero_source_groups = 0
    max_class_size = 0
    for class_kappa, raw_class_records in kappa_classes:
        class_records = list(raw_class_records)
        max_class_size = max(max_class_size, len(class_records))
        class_mu = sp.simplify(sp.sympify(class_kappa) + critical_value)
        if same(class_mu):
            for record in class_records:
                check(same(sp.sympify(record["weight"])), f"{case_name}: resonant class")
            continue

        alpha_square = sp.simplify(-class_mu / (16 * determinant))
        alpha = sp.sqrt(alpha_square)
        source_groups: list[list[sp.Expr]] = []
        normalized_records: list[tuple[sp.Expr, sp.Expr]] = []
        for record in class_records:
            edge_weight = sp.sympify(record["weight"])
            sigma = sp.simplify(edge_weight / (2 * alpha))
            edge_a = sp.sympify(record["a"])
            edge_b = sp.sympify(record["b"])
            check(same(sigma**2, edge_a), f"{case_name}: common-alpha square root")
            normalized_start = sp.simplify(edge_b / (2 * sigma))
            normalized_end = sp.simplify((2 * edge_a + edge_b) / (2 * sigma))
            check(
                same(normalized_start**2, class_kappa + record["c"]),
                f"{case_name}: class start square",
            )
            check(
                same(normalized_end**2, class_kappa + record["end_value"]),
                f"{case_name}: class end square",
            )
            add_group(source_groups, sp.sympify(record["c"]), -normalized_start)
            add_group(
                source_groups,
                sp.sympify(record["end_value"]),
                normalized_end,
            )
            normalized_records.append((sigma, sp.sympify(record["polynomial"])))

        live_sources = [entry for entry in source_groups if not same(entry[1])]
        if not live_sources:
            zero_source_blocks += 1
            for degree in range(5):
                block_coefficient = sum(
                    2 * sigma * sp.integrate(edge_q**degree, (t, 0, 1))
                    for sigma, edge_q in normalized_records
                )
                check(
                    same(block_coefficient),
                    f"{case_name}: zero-source entire block degree={degree}",
                )
            continue

        nonzero_flux_blocks += 1
        for source_value, source_coefficient in live_sources:
            check(source_coefficient != 0, f"{case_name}: live grouped source")
            check(
                not same(class_kappa + source_value),
                f"{case_name}: nonzero source is noncritical",
            )
            nonzero_source_groups += 1

    if nonzero_flux_blocks == 0:
        live_affine_groups = [entry for entry in affine_value_groups if not same(entry[1])]
        check(live_affine_groups, f"{case_name}: affine-only Euler source")

    return {
        "quadratic_edges": quadratic_edges,
        "affine_nonconstant": affine_nonconstant,
        "affine_constant": affine_constant,
        "resonant_edges": resonant_edges,
        "classes": len(kappa_classes),
        "max_class_size": max_class_size,
        "nonzero_flux_blocks": nonzero_flux_blocks,
        "zero_source_blocks": zero_source_blocks,
        "nonzero_source_groups": nonzero_source_groups,
    }


# Planted boundary cases.
standard_triangle = (
    (sp.Integer(0), sp.Integer(0)),
    (sp.Integer(1), sp.Integer(0)),
    (sp.Integer(0), sp.Integer(1)),
)
resonance_case = analyze_full_rank_case(
    "planted-resonance",
    standard_triangle,
    sp.diag(2, 3),
    (sp.Rational(1, 3), sp.Integer(0)),
    sp.Integer(5),
)
check(resonance_case["resonant_edges"] >= 1, "planted resonance detected")

affine_only_case = analyze_full_rank_case(
    "planted-affine-only",
    standard_triangle,
    sp.diag(1, -1),
    (sp.Integer(0), sp.Integer(0)),
    sp.Integer(2),
)
check(affine_only_case["affine_nonconstant"] >= 1, "isotropic affine edge")
check(affine_only_case["nonzero_flux_blocks"] == 0, "affine-only Euler carrier")

two_class_case = analyze_full_rank_case(
    "planted-two-edge-class",
    standard_triangle,
    sp.eye(2),
    (sp.Rational(1, 4), sp.Rational(1, 4)),
    sp.Integer(3),
)
check(two_class_case["max_class_size"] >= 2, "two-edge kappa class")

inradius = (2 - sp.sqrt(2)) / 2
full_class_case = analyze_full_rank_case(
    "planted-full-class",
    standard_triangle,
    sp.eye(2),
    (inradius, inradius),
    sp.Integer(7),
)
check(full_class_case["max_class_size"] == 3, "full tangent kappa class")

complex_case = analyze_full_rank_case(
    "planted-complex-algebraic",
    standard_triangle,
    sp.Matrix([[1 + ii, 2 - ii], [2 - ii, 3 + 2 * ii]]),
    (sp.Rational(1, 3) + ii / 5, -sp.Rational(2, 3) + ii / 7),
    1 + 2 * ii,
)

planted_totals = {
    key: sum(
        case[key]
        for case in (
            resonance_case,
            affine_only_case,
            two_class_case,
            full_class_case,
            complex_case,
        )
    )
    for key in resonance_case
}
print(
    "C10 full-rank planted cases=5 (including complex algebraic); "
    f"quadratic_edges={planted_totals['quadratic_edges']}; "
    f"affine_nonconstant={planted_totals['affine_nonconstant']}; "
    f"resonant_edges={planted_totals['resonant_edges']}; "
    f"max_class_size={max(case['max_class_size'] for case in (resonance_case, affine_only_case, two_class_case, full_class_case, complex_case))}"
)

# C11. Deterministic random exact triangles and full-rank forms.
rng = random.Random(3133)
random_case_count = 18
random_totals = {
    "quadratic_edges": 0,
    "affine_nonconstant": 0,
    "affine_constant": 0,
    "resonant_edges": 0,
    "classes": 0,
    "max_class_size": 0,
    "nonzero_flux_blocks": 0,
    "zero_source_blocks": 0,
    "nonzero_source_groups": 0,
}
for case_index in range(random_case_count):
    while True:
        random_vertices = tuple(
            (sp.Integer(rng.randint(-4, 4)), sp.Integer(rng.randint(-4, 4)))
            for _ in range(3)
        )
        if det2(
            vec_sub(random_vertices[1], random_vertices[0]),
            vec_sub(random_vertices[2], random_vertices[0]),
        ) != 0:
            break
    while True:
        matrix_a = rng.randint(-4, 4)
        matrix_b = rng.randint(-4, 4)
        matrix_c = rng.randint(-4, 4)
        random_matrix = sp.Matrix([[matrix_a, matrix_b], [matrix_b, matrix_c]])
        if random_matrix.det() != 0:
            break
    random_center = (
        sp.Rational(rng.randint(-5, 5), rng.randint(1, 4)),
        sp.Rational(rng.randint(-5, 5), rng.randint(1, 4)),
    )
    random_q0 = sp.Integer(rng.randint(-5, 5))
    case_result = analyze_full_rank_case(
        f"random-{case_index}",
        random_vertices,
        random_matrix,
        random_center,
        random_q0,
    )
    for key in random_totals:
        if key == "max_class_size":
            random_totals[key] = max(random_totals[key], case_result[key])
        else:
            random_totals[key] += case_result[key]
print(
    f"C11 random exact full-rank cases={random_case_count}; "
    f"quadratic_edges={random_totals['quadratic_edges']}; "
    f"affine_edges={random_totals['affine_nonconstant'] + random_totals['affine_constant']}; "
    f"kappa_classes={random_totals['classes']}; "
    f"live_source_groups={random_totals['nonzero_source_groups']}"
)

# C12. The two rational ODE obstructions and zero-source entireness recurrence.
kappa_symbol = sp.symbols("kappa_symbol")
entire_coefficient = sp.Integer(0)
for degree in range(10):
    previous = sp.Integer(0) if degree == 0 else entire_coefficient
    entire_coefficient = sp.simplify(
        -kappa_symbol * previous / (degree + sp.Rational(1, 2))
    )
    check(entire_coefficient == 0, f"zero-source entire recurrence {degree}")

for pole_order in range(2, 13):
    check(1 - pole_order != 0, f"unit-residue pole order {pole_order}")
residue_coefficient, regular0, regular1, nu = sp.symbols(
    "residue_coefficient regular0 regular1 nu"
)
laurent_trial = residue_coefficient / s + regular0 + regular1 * s
unit_residue_left = sp.expand(
    s * sp.diff(laurent_trial, s) + (1 + nu * s) * laurent_trial
)
check(
    sp.simplify(sp.limit(s * unit_residue_left, s, 0)) == 0,
    "simple-pole cancellation cannot make c/s",
)
print(
    "C12 ODE gates: half-residue nonresonance; unit-residue orders=2..12; "
    "simple-pole c/s coefficient=0; zero-source Taylor coefficients=10"
)

print("ALL THM-3143 CONTROLS PASSED")
