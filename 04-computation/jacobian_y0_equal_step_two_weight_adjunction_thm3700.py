#!/usr/bin/env python3
"""Exact minimal equal-step two-weight adjunction gate for THM-3700."""

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

    # THM-3698's two-weight span, enlarged by the minimal monomials of weights
    # -5 and +4.  The four weights form the step-three progression.
    functions = (
        C, B * C**2, A * C**3,
        A, B**2, A**2 * C**2, A * B * C,
        A**2 * B, C**4,
    )
    target_labels = (
        (0, 0, 1), (0, 1, 2), (1, 0, 3),
        (1, 0, 0), (0, 2, 0), (2, 0, 2), (1, 1, 1),
        (2, 1, 0), (0, 0, 4),
    )
    weights = tuple(-2 * a - b_power + c_power for a, b_power, c_power in target_labels)
    require(weights == (1, 1, 1, -2, -2, -2, -2, -5, 4), "weight profile")
    require(tuple(sorted(set(weights))) == (-5, -2, 1, 4), "equal-step weights")
    require(coefficient_matrix(list(functions), x, z).rank() == 9, "function span")

    left_collision = tuple(value.subs({x: 1, z: 0}) for value in functions)
    right_collision = tuple(value.subs({x: -1, z: 2}) for value in functions)
    require(left_collision == right_collision, "collision lost")

    pairs = [(i, j) for i in range(9) for j in range(i + 1, 9)]
    pair_index = {pair: index for index, pair in enumerate(pairs)}
    require(len(pairs) == 36, "exterior-square dimension")
    require(pair_index[(2, 6)] == 18 and pair_index[(7, 8)] == 35, "pair indexing")
    wedges = [bracket(functions[i], functions[j], x, z) for i, j in pairs]
    variables = sp.symbols("w0:36")
    expression = sp.expand(sum(a * value for a, value in zip(variables, wedges)) - 1)
    equations = coefficient_equations(expression, x, z)
    matrix, rhs = sp.linear_eq_to_matrix(equations, variables)
    augmented = matrix.row_join(rhs)
    require(len(equations) == 48, "coefficient-row count")
    require(matrix.rank() == augmented.rank() == 29, "affine-fibre rank")

    solution_set = sp.linsolve((matrix, rhs), variables)
    require(solution_set != sp.EmptySet, "affine fibre empty")
    solution = tuple(next(iter(solution_set)))
    parameters = sorted(set().union(*(entry.free_symbols for entry in solution)), key=str)
    expected_parameters = {variables[index] for index in (7, 15, 18, 25, 26, 30, 35)}
    require(set(parameters) == expected_parameters, "affine-fibre parameters")
    for equation in equations:
        require(sp.expand(equation.subs(dict(zip(variables, solution)))) == 0, "fibre replay")

    quadrics: list[sp.Expr] = []
    quadric_by_label: dict[tuple[int, int, int, int], sp.Expr] = {}
    for i in range(9):
        for j in range(i + 1, 9):
            for k in range(j + 1, 9):
                for ell in range(k + 1, 9):
                    quadric = sp.factor(
                        solution[pair_index[(i, j)]] * solution[pair_index[(k, ell)]]
                        - solution[pair_index[(i, k)]] * solution[pair_index[(j, ell)]]
                        + solution[pair_index[(i, ell)]] * solution[pair_index[(j, k)]]
                    )
                    quadrics.append(quadric)
                    quadric_by_label[(i, j, k, ell)] = quadric
    require(len(quadrics) == 126, "Pluecker count")
    unique_nonzero: list[sp.Expr] = []
    for quadric in quadrics:
        if quadric != 0 and quadric not in unique_nonzero and -quadric not in unique_nonzero:
            unique_nonzero.append(quadric)
    require(len(unique_nonzero) == 106, "unique nonzero Pluecker count")
    basis = sp.groebner(unique_nonzero, *parameters, order="grevlex", method="f5b")
    require(list(basis) == [1], "affine fibre met Gr(2,9)")

    # Two explicit Pluecker coordinates already contradict decomposability.
    u = variables[18]
    v = variables[35]
    first_core = sp.factor(quadric_by_label[(0, 1, 3, 4)])
    second_core = sp.factor(quadric_by_label[(0, 1, 3, 6)])
    require(
        sp.expand(first_core + sp.Rational(17, 56) * (u - 2 * v)) == 0,
        "first core minor",
    )
    expected_second = -(
        2304 * u**2 - 9216 * u * v + 9216 * v**2 + 343
    ) / 1176
    require(sp.expand(second_core - expected_second) == 0, "second core minor")
    require(
        sp.simplify(second_core.subs(u, 2 * v)) == -sp.Rational(7, 24),
        "two-minor contradiction",
    )

    # Positive controls: the original rational three-bracket identity remains
    # in the enlarged affine fibre, and the ambient orientation is unchanged.
    Q0 = -B**2 / 4 + A / 6 + A**2 * C**2 / 6
    additive_control = (
        bracket(C, Q0, x, z)
        + sp.Rational(7, 4) * bracket(B * C**2, A * B * C, x, z)
        - sp.Rational(13, 8) * bracket(A * C**3, B**2, x, z)
    )
    require(sp.expand(additive_control) == 1, "additive control")
    require(bracket(x, z, x, z) == 1, "ambient bracket orientation")

    source = Path(__file__).read_text(encoding="utf-8")
    require(
        not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
        "inactive Python assert found",
    )

    semantic_rows = [
        "functions=" + ";".join(map(str, functions)),
        "solution=" + ";".join(map(str, solution)),
        "quadrics=" + ";".join(map(str, quadrics)),
        f"core={first_core};{second_core}",
        "basis=" + str(list(basis)),
    ]
    semantic_hash = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()

    print("theorem=THM-3700-y0-equal-step-two-weight-adjunction-Pluecker-no-go")
    print("span=THM3698_seven_functions_plus_A^2B,C^4;dimension=9;weights=-5,-2,1,4")
    print("exterior_square=36;coefficient_rows=48;linear_rank=29;affine_dimension=7")
    print("Pluecker_quadrics=126;unique_nonzero=106;substituted_Groebner=[1]")
    print("hand_core=Delta0134=-17*(w18-2*w35)/56;Delta0136_at_w18=2w35=-7/24")
    print("conclusion=one_bracket_EMPTY_in_this_nine_function_span")
    print("controls=ambient_bracket:PASS;rational_three_bracket:PASS;collision_retained:PASS")
    print(f"semantic_sha256={semantic_hash}")
    print(f"CHECKS={CHECKS}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
