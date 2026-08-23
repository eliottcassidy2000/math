#!/usr/bin/env python3
"""Exact affine-bivector and Pluecker gate for THM-3698.

The bracket-to-one equation is linear on the exterior square of a chosen
seven-dimensional collision-ring span.  A single bracket is exactly a
decomposable bivector, so the remaining equations are the 35 Pluecker
quadrics for Gr(2,7).
"""

from __future__ import annotations

import ast
import hashlib
from pathlib import Path

import sympy as sp


CHECKS = 0


def require(condition: bool, detail: object) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(detail)


def bracket(first: sp.Expr, second: sp.Expr, x: sp.Symbol, z: sp.Symbol) -> sp.Expr:
    return sp.expand(
        sp.diff(first, x) * sp.diff(second, z)
        - sp.diff(first, z) * sp.diff(second, x)
    )


def coefficient_equations(polynomial: sp.Expr, x: sp.Symbol, z: sp.Symbol) -> list[sp.Expr]:
    return [coefficient for _monomial, coefficient in sp.Poly(polynomial, x, z).terms()]


def coefficient_matrix(polynomials: list[sp.Expr], x: sp.Symbol, z: sp.Symbol) -> sp.Matrix:
    support = sorted(
        {
            monomial
            for polynomial in polynomials
            for monomial, _coefficient in sp.Poly(polynomial, x, z).terms()
        }
    )
    return sp.Matrix(
        [
            [
                sp.Poly(polynomial, x, z).coeff_monomial(x**i * z**j)
                for polynomial in polynomials
            ]
            for i, j in support
        ]
    )


def main() -> None:
    x, z = sp.symbols("x z")
    t = x**2 * z
    A = sp.expand(3 * z * (2 - t))
    B = sp.expand(2 * x * z * (2 - t))
    C = sp.expand(x * (1 - t))

    U = (C, B * C**2, A * C**3)
    V = (A, B**2, A**2 * C**2, A * B * C)
    functions = U + V
    require(coefficient_matrix(list(functions), x, z).rank() == 7, "function span")
    target_labels = (
        (0, 0, 1), (0, 1, 2), (1, 0, 3),
        (1, 0, 0), (0, 2, 0), (2, 0, 2), (1, 1, 1),
    )
    function_weights = tuple(-2 * a - b_power + c_power for a, b_power, c_power in target_labels)
    require(function_weights == (1, 1, 1, -2, -2, -2, -2), "two-weight profile")
    require(2 * len({*function_weights, sp.Symbol("new_weight")}) == 6 < 7, "adjunction floor")

    left_collision = tuple(value.subs({x: 1, z: 0}) for value in functions)
    right_collision = tuple(value.subs({x: -1, z: 2}) for value in functions)
    require(left_collision == right_collision, "seven-function span lost collision")

    pair_index: dict[tuple[int, int], int] = {}
    wedges: list[sp.Expr] = []
    for first in range(7):
        for second in range(first + 1, 7):
            pair_index[(first, second)] = len(wedges)
            wedges.append(bracket(functions[first], functions[second], x, z))
    require(len(wedges) == 21, "exterior-square dimension")

    wedge_variables = sp.symbols("w0:21")
    generic = sp.expand(sum(a * value for a, value in zip(wedge_variables, wedges)) - 1)
    equations = coefficient_equations(generic, x, z)
    matrix, rhs = sp.linear_eq_to_matrix(equations, wedge_variables)
    augmented = matrix.row_join(rhs)
    require(len(equations) == 29, "coefficient-row count")
    require(matrix.rank() == augmented.rank() == 19, "affine-fibre rank")

    solution_set = sp.linsolve((matrix, rhs), wedge_variables)
    require(solution_set != sp.EmptySet, "affine fibre unexpectedly empty")
    solution = tuple(next(iter(solution_set)))
    parameters = sorted(set().union(*(entry.free_symbols for entry in solution)), key=str)
    require(len(parameters) == 2, "affine-fibre dimension")
    require(parameters == [wedge_variables[11], wedge_variables[14]], "fibre parameters")
    for equation in equations:
        require(sp.expand(equation.subs(dict(zip(wedge_variables, solution)))) == 0, "fibre replay")

    # The whole affine fibre is bipartite between U and V.  In the ordered
    # 3-by-4 cross block its exact two-parameter normal form is the matrix
    # displayed in the theorem.  Two minors already make rank one impossible.
    same_side_pairs = [
        pair_index[pair]
        for pair in (*((i, j) for i in range(3) for j in range(i + 1, 3)),
                     *((i, j) for i in range(3, 7) for j in range(i + 1, 7)))
    ]
    require(all(solution[index] == 0 for index in same_side_pairs), "non-bipartite fibre entry")
    alpha_parameter, beta_parameter = parameters
    cross_block = sp.Matrix(
        [[solution[pair_index[(row, 3 + column)]] for column in range(4)] for row in range(3)]
    )
    expected_cross_block = sp.Matrix(
        [
            [sp.Rational(1, 6), -sp.Rational(1, 4), (1 - 9 * alpha_parameter) / 6, -12 * beta_parameter / 7],
            [8 * beta_parameter / 7, 3 * beta_parameter / 28, -15 * beta_parameter / 14, sp.Rational(7, 4)],
            [alpha_parameter, -sp.Rational(13, 8), 0, beta_parameter],
        ]
    )
    require(cross_block == expected_cross_block, "cross-block normal form")
    first_minor = sp.factor(cross_block.extract((0, 1), (0, 1)).det())
    hostile_minor = sp.factor(cross_block.extract((0, 1), (0, 3)).det())
    require(first_minor == sp.Rational(17, 56) * beta_parameter, "first rank-one minor")
    require(
        sp.expand(hostile_minor - (2304 * beta_parameter**2 + 343) / 1176) == 0,
        "hostile rank-one minor",
    )
    require(hostile_minor.subs(beta_parameter, 0) == sp.Rational(7, 24), "two-minor contradiction")

    pluecker: list[sp.Expr] = []
    for i in range(7):
        for j in range(i + 1, 7):
            for k in range(j + 1, 7):
                for ell in range(k + 1, 7):
                    pluecker.append(
                        sp.factor(
                            solution[pair_index[(i, j)]] * solution[pair_index[(k, ell)]]
                            - solution[pair_index[(i, k)]] * solution[pair_index[(j, ell)]]
                            + solution[pair_index[(i, ell)]] * solution[pair_index[(j, k)]]
                        )
                    )
    require(len(pluecker) == 35, "Pluecker count")
    pluecker_basis = sp.groebner(pluecker, *parameters, order="grevlex")
    require(list(pluecker_basis) == [1], "affine fibre met Gr(2,7)")

    # Independent bipartite chart: P in span(U), Q in span(V).  The affine
    # 3-by-4 coefficient-matrix fibre has dimension two, but all three charts
    # of its rank-one cone are empty.
    table = [[bracket(left, right, x, z) for right in V] for left in U]
    matrix_variables = sp.symbols("m0:12")
    bipartite_generic = sp.expand(
        sum(matrix_variables[4 * row + column] * table[row][column]
            for row in range(3) for column in range(4))
        - 1
    )
    bipartite_equations = coefficient_equations(bipartite_generic, x, z)
    bipartite_matrix, bipartite_rhs = sp.linear_eq_to_matrix(
        bipartite_equations, matrix_variables
    )
    require(len(bipartite_equations) == 11, "bipartite row count")
    require(
        bipartite_matrix.rank()
        == bipartite_matrix.row_join(bipartite_rhs).rank()
        == 10,
        "bipartite fibre rank",
    )

    p = sp.symbols("p0:3")
    q = sp.symbols("q0:4")
    rank_one_expression = sp.expand(
        sum(p[row] * q[column] * table[row][column]
            for row in range(3) for column in range(4))
        - 1
    )
    rank_one_equations = coefficient_equations(rank_one_expression, x, z)
    pivot_bases = []
    for pivot in range(3):
        variables = [entry for index, entry in enumerate(p) if index != pivot] + list(q)
        specialized = [sp.expand(equation.subs(p[pivot], 1)) for equation in rank_one_equations]
        basis = sp.groebner(specialized, *variables, order="grevlex", method="f5b")
        require(list(basis) == [1], ("bipartite rank-one chart", pivot))
        pivot_bases.append(str(list(basis)))

    # Positive controls.  The rational three-bracket identity in THM-3686 is
    # a point of this affine fibre, represented by a rank-six skew matrix.
    rational_wedge = [sp.Integer(0)] * 21
    rational_wedge[pair_index[(0, 3)]] = sp.Rational(1, 6)
    rational_wedge[pair_index[(0, 4)]] = -sp.Rational(1, 4)
    rational_wedge[pair_index[(0, 5)]] = sp.Rational(1, 6)
    rational_wedge[pair_index[(1, 6)]] = sp.Rational(7, 4)
    rational_wedge[pair_index[(2, 4)]] = -sp.Rational(13, 8)
    require(
        sp.expand(sum(a * value for a, value in zip(rational_wedge, wedges))) == 1,
        "additive-bracket control",
    )
    skew = sp.zeros(7)
    for (first, second), index in pair_index.items():
        skew[first, second] = rational_wedge[index]
        skew[second, first] = -rational_wedge[index]
    require(skew.rank() == 6, "additive control rank")
    require(bracket(x, z, x, z) == 1, "ambient bracket orientation")

    source = Path(__file__).read_text(encoding="utf-8")
    require(
        not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
        "inactive Python assert found",
    )

    semantic_rows = [
        "functions=" + ";".join(map(str, functions)),
        "solution=" + ";".join(map(str, solution)),
        f"two_minors={first_minor};{hostile_minor}",
        "pluecker=" + ";".join(map(str, pluecker)),
        "pluecker_basis=" + str(list(pluecker_basis)),
        "bipartite_pivots=" + ";".join(pivot_bases),
        "rational_rank=6",
        "weights=1,-2;one-new-weight-total<=6",
    ]
    semantic_hash = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()

    print("theorem=THM-3698-y0-seven-function-Pluecker-compression-no-go")
    print("span=U:(C,BC^2,AC^3);V:(A,B^2,A^2C^2,ABC);dimension=7")
    print("exterior_square=21;coefficient_rows=29;linear_rank=19;affine_dimension=2")
    print("Pluecker_quadrics=35;substituted_Groebner=[1];one_bracket=EMPTY")
    print("hand_certificate=minor01_01=17*b/56;minor01_03=(2304*b^2+343)/1176;rank1_impossible")
    print("bipartite_3x4=rows:11,rank:10,affine_dimension:2;rank_one_pivot_charts=3/3_unit")
    print("weights=1,-2;one_new_homogeneous_weight_gives_total_support<=6;THM-3695_excludes")
    print("controls=ambient_bracket:PASS;rational_three_bracket_rank6:PASS;collision_retained:PASS")
    print(f"semantic_sha256={semantic_hash}")
    print(f"CHECKS={CHECKS}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
