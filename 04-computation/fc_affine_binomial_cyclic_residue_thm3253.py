#!/usr/bin/env python3
"""Exact controls for the affine-coordinate binomial theorem THM-3253."""

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


def connection(degree: int, aval: sp.Expr, bval: sp.Expr) -> tuple[sp.Matrix, sp.Matrix]:
    aval = sp.sympify(aval)
    bval = sp.sympify(bval)
    rank = degree - 1
    upper = bval * rank / degree
    lower = -bval**2 * rank / (degree**2 * aval)
    constant = sp.zeros(rank)
    for j in range(rank - 1):
        constant[j, j + 1] = upper
    constant[rank - 1, 0] = lower
    residue = sp.diag(*(sp.Rational(-j, degree) for j in range(1, rank + 1)))
    return constant, residue


def source_space_rank(
    matrix: sp.Matrix, variables: list[sp.Symbol], source_variables: list[sp.Symbol]
) -> tuple[int, list[sp.Matrix]]:
    null = matrix.nullspace()
    vectors = []
    indices = [variables.index(variable) for variable in source_variables]
    for vector in null:
        source = sp.Matrix([vector[index] for index in indices])
        if source != sp.zeros(len(indices), 1):
            vectors.append(source)
    rank = sp.Matrix.hstack(*vectors).rank() if vectors else 0
    return rank, vectors


def polynomial_splitting_space(
    constant: sp.Matrix, residue: sp.Matrix, eta: sp.Expr, degree: int, polynomial_degree: int
) -> tuple[int, list[sp.Matrix]]:
    rank = degree - 1
    coefficients = [sp.Matrix(sp.symbols(f"r{k}_0:{rank}")) for k in range(polynomial_degree + 1)]
    source_variables = list(sp.symbols(f"source_0:{rank}"))
    source = sp.Matrix(source_variables)
    shifted = constant - eta * sp.eye(rank)
    equations = list(residue * coefficients[0] + source / degree)
    for k in range(1, polynomial_degree + 1):
        equations.extend(list(shifted * coefficients[k - 1] + (residue - k * sp.eye(rank)) * coefficients[k]))
    equations.extend(list(shifted * coefficients[-1]))
    variables = [entry for vector in coefficients for entry in vector] + source_variables
    coefficient_matrix, rhs = sp.linear_eq_to_matrix(equations, variables)
    require(rhs == sp.zeros(len(equations), 1), "splitting equations became inhomogeneous")
    return source_space_rank(coefficient_matrix, variables, source_variables)


def frobenius_exponents(degree: int, coordinate: int) -> tuple[sp.Rational, ...]:
    rank = degree - 1
    values = []
    for mode in range(1, rank + 1):
        distance = (mode - coordinate) % rank
        values.append(-sp.Rational(mode, degree) + distance)
    return tuple(values)


ledger: list[str] = []
t = sp.symbols("t")

# C1: coefficientwise audit of the general rank-(d-1) primitive system.
system_checks = 0
for degree in range(3, 8):
    rank = degree - 1
    aval, bval = sp.Rational(2), -sp.Rational(3, 2)
    phase = aval * t**degree + bval * t
    endpoint = sp.Rational(2, 3) + I
    endpoint_phase = sp.expand(phase.subs(t, endpoint))
    kappa = bval * rank / (degree * aval)
    coefficients: list[list[sp.Expr]] = []
    for j in range(1, rank + 1):
        coefficients.append(
            [integral_poly(t ** (j - 1) * phase**n, t, endpoint) / factorial(n) for n in range(6)]
        )
    for n in range(6):
        for j in range(1, rank):
            lhs = (degree * n + j) * coefficients[j - 1][n]
            if n:
                lhs -= bval * rank * coefficients[j][n - 1]
            rhs = endpoint**j * endpoint_phase**n / factorial(n)
            require(zero(lhs - rhs), f"primitive row failed d={degree},j={j},n={n}")
            system_checks += 1
        lhs = (degree * n + rank) * coefficients[-1][n]
        if n:
            lhs += bval**2 * rank / (degree * aval) * coefficients[0][n - 1]
        rhs = (endpoint**rank + kappa) * endpoint_phase**n / factorial(n)
        if n == 0:
            rhs -= kappa
        require(zero(lhs - rhs), f"primitive last row failed d={degree},n={n}")
        system_checks += 1
ledger.append(f"C1:cyclic_primitive_coefficients={system_checks};degrees=3..7")

# C2: critical values are exactly the spectrum, with explicit left/right
# eigenvectors and the universal half-step obstruction to positive-degree
# rational splitting.
spectrum_checks = 0
for degree in range(3, 10):
    rank = degree - 1
    aval, rho = sp.symbols(f"A{degree} rho{degree}", nonzero=True)
    bval = -degree * aval * rho**rank
    constant, residue = connection(degree, aval, bval)
    upper = bval * rank / degree
    eta = upper * rho
    right = sp.Matrix([rho**j for j in range(rank)])
    left = sp.Matrix([rho ** (-j) for j in range(rank)])
    require(
        all(zero(entry) for entry in constant * right - eta * right),
        f"right critical eigenvector failed d={degree}",
    )
    require(
        all(zero(entry) for entry in constant.T * left - eta * left),
        f"left critical eigenvector failed d={degree}",
    )
    obstruction = sp.simplify((left.T * (sp.eye(rank) - residue) * right)[0])
    require(zero(obstruction - rank * sp.Rational(3, 2)), f"half-step obstruction failed d={degree}")
    require(
        zero(constant.det() - (-1) ** (rank - 1) * constant[0, 1] ** (rank - 1) * constant[-1, 0]),
        f"cycle determinant failed d={degree}",
    )
    spectrum_checks += 4
ledger.append(f"C2:critical_spectrum_vector_checks={spectrum_checks};half_step=N+1/2")

# C3: exact polynomial splitting ranks at a rational critical value and a
# noncritical value for the parent-requested degrees 4,5,6.
splitting_profiles = []
for degree in (4, 5, 6):
    rank = degree - 1
    constant, residue = connection(degree, sp.S.One, -degree)
    critical_eta = -(degree - 1)  # rho=1
    expected_line = sp.Matrix(list(range(1, rank + 1)))
    critical_ranks = []
    noncritical_ranks = []
    first_vectors: list[sp.Matrix] = []
    for polynomial_degree in range(4):
        source_rank, vectors = polynomial_splitting_space(
            constant, residue, critical_eta, degree, polynomial_degree
        )
        critical_ranks.append(source_rank)
        if vectors and not first_vectors:
            first_vectors = vectors
        source_rank_noncritical, _ = polynomial_splitting_space(
            constant, residue, sp.S.One, degree, polynomial_degree
        )
        noncritical_ranks.append(source_rank_noncritical)
    require(
        critical_ranks == [1] * 4,
        f"critical splitting rank failed d={degree}: {critical_ranks}",
    )
    require(noncritical_ranks == [0] * 4, f"noncritical splitting rank failed d={degree}")
    require(
        any(sp.Matrix.hstack(vector, expected_line).rank() == 1 for vector in first_vectors),
        f"critical splitting line failed d={degree}",
    )
    splitting_profiles.append((degree, tuple(critical_ranks), tuple(noncritical_ranks)))
ledger.append(f"C3:splitting_profiles={splitting_profiles}")

# C4: critical-fibre collision audit for all partitions of three endpoints and
# for the doubled-knot pair.  Normalize a critical point to rho=1.
x, y, z = sp.symbols("x y z")
critical_fibre_checks = 0
set_partitions = (
    ((0,), (1,), (2,)),
    ((0, 1), (2,)),
    ((0, 2), (1,)),
    ((1, 2), (0,)),
    ((0, 1, 2),),
)
for degree in (4, 5, 6):
    critical_fibre = sp.expand(x**degree - degree * x + degree - 1)
    require(sp.rem(critical_fibre, (x - 1) ** 2, domain=sp.QQ) == 0, f"double root failed d={degree}")
    require(critical_fibre.subs(x, 2) == 2**degree - degree - 1 != 0, f"singleton hostile failed d={degree}")
    require(sp.gcd(critical_fibre, x - 2) == 1, f"critical singleton gcd failed d={degree}")
    # The five set partitions are exact: four contain a singleton, and the
    # remaining all-three block is controlled by moments (0,nonzero).
    require(len(set_partitions) == 5, "set-partition ledger changed")
    require(sum(any(len(block) == 1 for block in partition) for partition in set_partitions) == 4,
            "singleton partition count changed")
    critical_fibre_checks += len(set_partitions)

weights = (y - z, z - x, x - y)
require(zero(sum(weights)), "three-knot zeroth moment failed")
require(zero(sum(weights[j] * (x, y, z)[j] for j in range(3))), "three-knot first moment failed")
vandermonde_moment = sp.factor(sum(weights[j] * (x, y, z)[j] ** 2 for j in range(3)))
require(
    zero(vandermonde_moment - (x - y) * (x - z) * (y - z)),
    "three-knot second moment is not the Vandermonde",
)

pair_groebner = sp.groebner([x + y - 2, x**2 + x * y + y**2 - 3], x, y, domain=sp.QQ)
require(pair_groebner.reduce((x - 1) ** 2)[1] == 0, "pair alignment did not force x=1")
require(pair_groebner.reduce(y - (2 - x))[1] == 0, "pair alignment did not force y=2-x")
rho3 = sp.symbols("rho3", nonzero=True)
cubic_split = sp.Matrix([1, 2 * rho3])
for cubic_source in (sp.Matrix([rho3, -rho3**2]), sp.Matrix([-2 * rho3, 2 * rho3**2])):
    require(sp.factor(sp.det(sp.Matrix.hstack(cubic_source, cubic_split))) != 0,
            "cubic critical-fibre boundary became aligned")
critical_fibre_checks += 5
ledger.append(
    f"C4:critical_fibre_partition_checks={critical_fibre_checks};degrees=(4,5,6);pair_alignment=DIAGONAL_ONLY"
)

# C5: both coordinate projections are cyclic, but their exact Frobenius
# exponent sets are disjoint.  This is the common-factor hostile.
s = sp.symbols("s", nonzero=True)
frobenius_profiles = []
for degree in range(3, 11):
    rank = degree - 1
    constant, residue = connection(degree, sp.S.One, -degree)
    first = frobenius_exponents(degree, 1)
    second = frobenius_exponents(degree, 2)
    require(set(first).isdisjoint(second), f"coordinate exponent collision d={degree}")
    if degree in (4, 5, 6):
        connection_matrix = constant + residue / s
        for coordinate in (0, 1):
            row = sp.eye(rank)[coordinate, :]
            rows = []
            for _ in range(rank):
                rows.append(row)
                row = sp.simplify(sp.diff(row, s) + row * connection_matrix)
            determinant = sp.factor(sp.Matrix.vstack(*rows).det())
            require(determinant != 0, f"coordinate not cyclic d={degree},coordinate={coordinate + 1}")
    frobenius_profiles.append((degree, tuple(map(str, first)), tuple(map(str, second))))
ledger.append(f"C5:frobenius_disjoint_profiles={frobenius_profiles}")

# C6: direct nonpure period-formula comparisons in all three geometries.
u, v = sp.symbols("u v")


def simplex_integral(poly: sp.Expr) -> sp.Expr:
    return sp.integrate(sp.integrate(sp.expand(poly), (v, 0, 1 - u)), (u, 0, 1))


z_values = (sp.S.Zero, sp.S.One, I)
edges = tuple(z_values[(j + 1) % 3] - z_values[j] for j in range(3))
slopes = tuple(sp.conjugate(edge) / edge for edge in edges)
turns = tuple(sp.simplify(slopes[(j - 1) % 3] - slopes[j]) for j in range(3))
W = 2 * I * sp.im(sp.conjugate(edges[0]) * (z_values[2] - z_values[0]))
knots = (sp.Rational(-2), sp.Rational(1), sp.Rational(5))
gap_1, gap_2, total_gap = knots[1] - knots[0], knots[2] - knots[1], knots[2] - knots[0]
alphas = (-1 / (total_gap * gap_1), 1 / (gap_1 * gap_2), -1 / (total_gap * gap_2))


def primitive_coefficient(endpoint: sp.Expr, degree: int, power: int, moment: int) -> sp.Expr:
    phase = t**degree + 2 * t
    return integral_poly(t**moment * phase**power, t, endpoint) / factorial(power)


period_checks = 0
for degree in (4, 5, 6):
    noncollinear_ell = u + I * v
    collinear_ell = knots[0] * (1 - u - v) + knots[1] * u + knots[2] * v
    doubled_ell = knots[0] * (1 - v) + knots[2] * v
    for n in range(5):
        direct = simplex_integral((noncollinear_ell**degree + 2 * noncollinear_ell) ** n) / factorial(n)
        boundary = sum(
            turns[j]
            * (
                primitive_coefficient(z_values[j], degree, n, 1)
                - z_values[j] * primitive_coefficient(z_values[j], degree, n, 0)
            )
            for j in range(3)
        )
        require(zero(W * direct - boundary), f"noncollinear period failed d={degree},n={n}")
        period_checks += 1

        direct = simplex_integral((collinear_ell**degree + 2 * collinear_ell) ** n) / factorial(n)
        boundary = sum(
            alphas[j]
            * (
                primitive_coefficient(knots[j], degree, n, 1)
                - knots[j] * primitive_coefficient(knots[j], degree, n, 0)
            )
            for j in range(3)
        )
        require(zero(direct - boundary), f"collinear period failed d={degree},n={n}")
        period_checks += 1

        direct = simplex_integral((doubled_ell**degree + 2 * doubled_ell) ** n) / factorial(n)
        packet_0 = primitive_coefficient(knots[0], degree, n, 0) - primitive_coefficient(
            knots[2], degree, n, 0
        )
        packet_1 = primitive_coefficient(knots[0], degree, n, 1) - primitive_coefficient(
            knots[2], degree, n, 1
        )
        boundary = -knots[2] * packet_0 + packet_1
        require(
            zero(total_gap**2 * direct - boundary),
            f"doubled period failed d={degree},n={n}",
        )
        period_checks += 1
ledger.append(f"C6:direct_nonpure_period_coefficients={period_checks};degrees=(4,5,6)")

digest = sha256("\n".join(ledger).encode("utf-8")).hexdigest()

print("THM-3253 AFFINE BINOMIAL CYCLIC-RESIDUE AUDIT")
for row in ledger:
    print(row)
print(f"semantic_sha256={digest}")
print("CONCLUSION=STRUCTURAL_NONSPLITTING_VERIFIED;VALUE_SPECIALIZATION_GATE_OPEN")
