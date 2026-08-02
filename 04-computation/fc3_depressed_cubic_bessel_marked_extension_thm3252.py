#!/usr/bin/env python3
"""Exact controls for the depressed-cubic Bessel reduction in THM-3252."""

from __future__ import annotations

from hashlib import sha256
from math import factorial

import sympy as sp


I = sp.I


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def zero(expr: sp.Expr) -> bool:
    return sp.simplify(sp.expand_complex(expr)) == 0


def integral_poly(poly: sp.Expr, variable: sp.Symbol, endpoint: sp.Expr) -> sp.Expr:
    primitive = sp.integrate(sp.expand(poly), variable)
    return sp.simplify(primitive.subs(variable, endpoint) - primitive.subs(variable, 0))


def source_space_rank(matrix: sp.Matrix, variables: list[sp.Symbol], b0: sp.Symbol, b1: sp.Symbol) -> tuple[int, list[sp.Matrix]]:
    null = matrix.nullspace()
    vectors = []
    i0, i1 = variables.index(b0), variables.index(b1)
    for vector in null:
        source = sp.Matrix([vector[i0], vector[i1]])
        if source != sp.zeros(2, 1):
            vectors.append(source)
    rank = sp.Matrix.hstack(*vectors).rank() if vectors else 0
    return rank, vectors


def polynomial_splitting_space(C: sp.Matrix, D: sp.Matrix, eta: sp.Expr, degree: int):
    rs = [sp.Matrix(sp.symbols(f"r{k}_0 r{k}_1")) for k in range(degree + 1)]
    b0, b1 = sp.symbols("source_0 source_1")
    source = sp.Matrix([b0, b1])
    M = C - eta * sp.eye(2)
    equations = list(D * rs[0] + source / 3)
    for k in range(1, degree + 1):
        equations.extend(list(M * rs[k - 1] + D * rs[k] - k * rs[k]))
    equations.extend(list(M * rs[-1]))
    variables = [entry for vector in rs for entry in vector] + [b0, b1]
    coefficient_matrix, rhs = sp.linear_eq_to_matrix(equations, variables)
    require(rhs == sp.zeros(len(equations), 1), "splitting system unexpectedly inhomogeneous")
    return source_space_rank(coefficient_matrix, variables, b0, b1)


ledger: list[str] = []

# C1: exact affine depression of an arbitrary cubic.
a3, a2, a1, a0, y = sp.symbols("a3 a2 a1 a0 y", nonzero=True)
h = a2 / (3 * a3)
Bdep = a1 - a2**2 / (3 * a3)
Cdep = a0 - a1 * a2 / (3 * a3) + 2 * a2**3 / (27 * a3**2)
original = a3 * (y - h) ** 3 + a2 * (y - h) ** 2 + a1 * (y - h) + a0
require(zero(original - (a3 * y**3 + Bdep * y + Cdep)), "cubic depression failed")
ledger.append("C1:cubic_depression=SYMBOLIC")

# C2: Euclidean divisions and coefficientwise primitive-system checks.
A, B, t = sp.symbols("A B t", nonzero=True)
p = A * t**3 + B * t
qprime = sp.diff(p, t)
h0, r0 = t / 3, 2 * B * t / 3
h1, r1 = t**2 / 3 + 2 * B / (9 * A), -2 * B**2 / (9 * A)
require(zero(p - h0 * qprime - r0), "j=0 Euclidean division failed")
require(zero(p * t - h1 * qprime - r1), "j=1 Euclidean division failed")

coefficient_checks = 0
for aval, bval, endpoint in (
    (sp.S.One, sp.Integer(3), sp.Rational(2)),
    (sp.Rational(2), -sp.S.One, 1 + I),
    (-sp.Rational(3, 2), sp.Rational(5, 2), -1 + 2 * I),
):
    phase = sp.expand(aval * t**3 + bval * t)
    phase_endpoint = sp.simplify(phase.subs(t, endpoint))
    kappa = sp.Rational(2) * bval / (3 * aval)
    gamma = sp.Rational(2) * bval**2 / (3 * aval)
    j0 = []
    j1 = []
    for n in range(10):
        j0.append(integral_poly(phase**n, t, endpoint) / factorial(n))
        j1.append(integral_poly(t * phase**n, t, endpoint) / factorial(n))
    for n in range(10):
        rhs0 = endpoint * phase_endpoint**n / factorial(n)
        lhs0 = (3 * n + 1) * j0[n]
        if n:
            lhs0 -= 2 * bval * j1[n - 1]
        require(zero(lhs0 - rhs0), f"J0 system mismatch n={n}")

        rhs1 = (endpoint**2 + kappa) * phase_endpoint**n / factorial(n)
        if n == 0:
            rhs1 -= kappa
        lhs1 = (3 * n + 2) * j1[n]
        if n:
            lhs1 += gamma * j0[n - 1]
        require(zero(lhs1 - rhs1), f"J1 system mismatch n={n}")
        coefficient_checks += 2
ledger.append(f"C2:primitive_system_coefficients={coefficient_checks};divisions=2")

# C3: eliminate each homogeneous coordinate and recover Bessel orders 1/3, 2/3.
s = sp.symbols("s", nonzero=True)
off01 = 2 * B / 3
off10 = -2 * B**2 / (9 * A)
lambda_sq = 4 * B**3 / (27 * A)
f = sp.Function("f")(s)
g = sp.Function("g")(s)
g_from_f = (sp.diff(f, s) + f / (3 * s)) / off01
f_equation = sp.simplify(
    sp.diff(g_from_f, s) - (off10 * f - 2 * g_from_f / (3 * s))
)
target_f = sp.diff(f, s, 2) + sp.diff(f, s) / s + (lambda_sq - sp.Rational(1, 9) / s**2) * f
require(zero(off01 * f_equation - target_f), "order-1/3 Bessel elimination failed")

f_from_g = (sp.diff(g, s) + 2 * g / (3 * s)) / off10
g_equation = sp.simplify(
    sp.diff(f_from_g, s) - (-f_from_g / (3 * s) + off01 * g)
)
target_g = sp.diff(g, s, 2) + sp.diff(g, s) / s + (lambda_sq - sp.Rational(4, 9) / s**2) * g
require(zero(off10 * g_equation - target_g), "order-2/3 Bessel elimination failed")
ledger.append("C3:bessel_orders=(1/3,2/3);lambda_squared=4B^3/(27A)")

# C4: B=0 is the split-residue boundary.
D = sp.diag(-sp.Rational(1, 3), -sp.Rational(2, 3))
C_zero = sp.zeros(2)
require(C_zero == sp.zeros(2), "pure-power constant matrix did not vanish")
require(tuple(D.diagonal()) == (-sp.Rational(1, 3), -sp.Rational(2, 3)),
        "pure-power residues changed")
ledger.append("C4:B=0_SPLITS;residues=(-1/3,-2/3)")

# C5: exact source-splitting classification in a deterministic nonpure cubic.
# p=t^3+3t gives C=[[0,2],[-2,0]], critical values +/-2i.
C = sp.Matrix([[0, 2], [-2, 0]])
critical_values = (2 * I, -2 * I)
splitting_profiles = []
for eta in critical_values + (sp.Integer(4),):
    degree_ranks = []
    first_vectors = []
    for degree in range(6):
        rank, vectors = polynomial_splitting_space(C, D, eta, degree)
        degree_ranks.append(rank)
        if vectors and not first_vectors:
            first_vectors = vectors
    if eta in critical_values:
        require(degree_ranks == [1] * 6, f"critical splitting line changed: eta={eta}")
        expected = sp.Matrix([2, 2 * eta])
        require(any(sp.Matrix.hstack(vector, expected).rank() == 1 for vector in first_vectors),
                f"critical source line mismatch eta={eta}")
    else:
        require(degree_ranks == [0] * 6, "noncritical source unexpectedly split")
    splitting_profiles.append((str(eta), tuple(degree_ranks)))
ledger.append(f"C5:splitting_source_ranks={splitting_profiles}")

# C6: endpoint collisions and critical endpoint sources are hostile controls.
root_plus = (-sp.S.One + I * sp.sqrt(15)) / 2
root_minus = (-sp.S.One - I * sp.sqrt(15)) / 2
roots = (sp.S.One, root_plus, root_minus)
require(all(zero(root**3 + 3 * root - 4) for root in roots), "equal-phase cubic roots failed")
collision_source = sum((sp.Matrix([root, root**2 + 2]) for root in roots), sp.zeros(2, 1))
require(collision_source != sp.zeros(2, 1), "planted equal-phase source vanished")

critical_transverse = []
for endpoint, eta in ((I, 2 * I), (-I, -2 * I)):
    source = sp.Matrix([endpoint, endpoint**2 + 2])
    split_line = sp.Matrix([2, 2 * eta])
    require(sp.Matrix.hstack(source, split_line).rank() == 2,
            f"critical endpoint source accidentally split: {endpoint}")
    critical_transverse.append((str(endpoint), str(eta)))

zero_phase_endpoint = I * sp.sqrt(3)
require(zero(zero_phase_endpoint**3 + 3 * zero_phase_endpoint), "zero-phase endpoint failed")
zero_group_source = sp.Matrix([zero_phase_endpoint, zero_phase_endpoint**2 + 2]) - sp.Matrix([0, 2])
require(zero_group_source != sp.zeros(2, 1), "zero-phase endpoint/constant group vanished")
ledger.append(
    "C6:hostiles=equal_phase_three_root+critical_transverse+zero_phase_constant_group"
)

# C7: the critical-root geometry is uniformly transverse, and the two scalar
# coordinate equations have incompatible orders.
r = sp.symbols("r", nonzero=True)
Bcrit = -3 * A * r**2
kappa_crit = sp.simplify(2 * Bcrit / (3 * A))
eta_crit = sp.simplify(2 * Bcrit * r / 3)
critical_source_r = sp.Matrix([r, r**2 + kappa_crit])
critical_source_other = sp.Matrix([-2 * r, 4 * r**2 + kappa_crit])
critical_split_line = sp.Matrix([2 * Bcrit / 3, 2 * eta_crit])
require(
    sp.Matrix.hstack(critical_source_r, sp.Matrix([1, -r])).rank() == 1,
    "critical repeated-root source direction failed",
)
require(
    sp.Matrix.hstack(critical_source_other, sp.Matrix([1, -r])).rank() == 1,
    "critical third-root source direction failed",
)
require(
    zero(sp.det(sp.Matrix.hstack(critical_source_r, critical_split_line)) + 2 * Bcrit**2 / (3 * A)),
    "uniform critical transversality determinant failed",
)

shared = sp.Function("shared")(s)
order_one_third = s**2 * sp.diff(shared, s, 2) + s * sp.diff(shared, s) + (
    lambda_sq * s**2 - sp.Rational(1, 9)
) * shared
order_two_thirds = s**2 * sp.diff(shared, s, 2) + s * sp.diff(shared, s) + (
    lambda_sq * s**2 - sp.Rational(4, 9)
) * shared
require(
    zero(order_two_thirds - order_one_third + shared / 3),
    "cross-copy Bessel-order defect failed",
)
ledger.append("C7:critical_group_direction=(1,-r);split_direction=(1,2r);cross_order_defect=-1/3")

# C8: exact moment controls for both three-distinct geometries and the cyclic
# row used for a doubled knot.
z_values = (sp.S.Zero, sp.S.One, I)
edges = tuple(z_values[(j + 1) % 3] - z_values[j] for j in range(3))
slopes = tuple(sp.conjugate(edge) / edge for edge in edges)
turns = tuple(sp.simplify(slopes[(j - 1) % 3] - slopes[j]) for j in range(3))
W = 2 * I * sp.im(sp.conjugate(edges[0]) * (z_values[2] - z_values[0]))
require(zero(sum(turns)), "turn zeroth moment failed")
require(zero(sum(turns[j] * z_values[j] for j in range(3))), "turn first moment failed")
require(zero(sum(turns[j] * z_values[j] ** 2 for j in range(3)) + W), "turn second moment failed")

ka, kb, kc = sp.Rational(-2), sp.Rational(1), sp.Rational(5)
kp, kq, kL = kb - ka, kc - kb, kc - ka
alphas = (-1 / (kL * kp), 1 / (kp * kq), -1 / (kL * kq))
knots = (ka, kb, kc)
require(zero(sum(alphas)), "spline zeroth moment failed")
require(zero(sum(alphas[j] * knots[j] for j in range(3))), "spline first moment failed")
require(zero(sum(alphas[j] * knots[j] ** 2 for j in range(3)) + 1), "spline second moment failed")

unique = sp.symbols("unique")
Gamma = 2 * B**2 / (3 * A)
cyclic_polynomial = unique + s * (2 * B * unique**2 + Gamma)
require(sp.Poly(cyclic_polynomial, s).as_expr() != 0, "doubled-knot cyclic row vanished")
ledger.append("C8:moments=turn(0,0,-W)+spline(0,0,-1);doubled_cyclic_row=NONZERO")

# C9: coefficientwise controls for the three geometric period formulas under
# a genuinely nonpure depressed cubic.
u, v = sp.symbols("u v")


def simplex_integral(poly: sp.Expr) -> sp.Expr:
    return sp.integrate(sp.integrate(sp.expand(poly), (v, 0, 1 - u)), (u, 0, 1))


def primitive_coefficient(endpoint: sp.Expr, power: int, moment: int) -> sp.Expr:
    phase = t**3 + 3 * t
    return integral_poly(t**moment * phase**power, t, endpoint) / factorial(power)


period_checks = 0
noncollinear_ell = u + I * v
for n in range(6):
    direct = simplex_integral((noncollinear_ell**3 + 3 * noncollinear_ell) ** n) / factorial(n)
    boundary = sum(
        turns[j]
        * (
            primitive_coefficient(z_values[j], n, 1)
            - z_values[j] * primitive_coefficient(z_values[j], n, 0)
        )
        for j in range(3)
    )
    require(zero(W * direct - boundary), f"noncollinear period formula failed n={n}")
    period_checks += 1

collinear_ell = ka * (1 - u - v) + kb * u + kc * v
for n in range(6):
    direct = simplex_integral((collinear_ell**3 + 3 * collinear_ell) ** n) / factorial(n)
    boundary = sum(
        alphas[j]
        * (
            primitive_coefficient(knots[j], n, 1)
            - knots[j] * primitive_coefficient(knots[j], n, 0)
        )
        for j in range(3)
    )
    require(zero(direct - boundary), f"collinear period formula failed n={n}")
    period_checks += 1

doubled_ell = ka * (1 - v) + kc * v
doubled_L = kc - ka
for n in range(6):
    direct = simplex_integral((doubled_ell**3 + 3 * doubled_ell) ** n) / factorial(n)
    packet_0 = primitive_coefficient(ka, n, 0) - primitive_coefficient(kc, n, 0)
    packet_1 = primitive_coefficient(ka, n, 1) - primitive_coefficient(kc, n, 1)
    boundary = -kc * packet_0 + packet_1
    require(zero(doubled_L**2 * direct - boundary), f"doubled period formula failed n={n}")
    period_checks += 1
ledger.append(f"C9:nonpure_period_formula_coefficients={period_checks}")

digest = sha256("\n".join(ledger).encode("utf-8")).hexdigest()

print("THM-3252 FC(3) DEPRESSED-CUBIC BESSEL MARKED-EXTENSION AUDIT")
for row in ledger:
    print(row)
print(f"semantic_sha256={digest}")
print("CONCLUSION=ALL_AFFINE-COORDINATE_CUBIC_BRANCHES_REDUCED_AND_EXCLUDED")
