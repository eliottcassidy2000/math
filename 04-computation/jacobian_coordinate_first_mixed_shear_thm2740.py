#!/usr/bin/env python3
"""Exact symbolic companion for THM-2740.

The theorem is an all-degree coordinate-ring argument.  This script checks
the coordinate chain rule on a generic bounded polynomial, two nonlinear
coordinate controls, the mixed THM-2696 shear pair, both quadratic graph
sections, their inverse, and the characteristic-p boundary.  It uses exact
symbolic arithmetic and no truth-bearing Python assert.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def jacobian(left: sp.Expr, right: sp.Expr, x: sp.Symbol, y: sp.Symbol) -> sp.Expr:
    return sp.expand(sp.diff(left, x) * sp.diff(right, y)
                     - sp.diff(left, y) * sp.diff(right, x))


def main() -> None:
    x, y, t = sp.symbols("x y t")
    kappa = sp.symbols("kappa", nonzero=True)

    # Generic degree-four check of the graph coordinate identity.
    coefficients = {
        (i, j): sp.symbols(f"c_{i}_{j}")
        for i in range(5)
        for j in range(5 - i)
    }
    polynomial_xt = sum(
        coefficient * x**i * t**j
        for (i, j), coefficient in coefficients.items()
    )
    t_xy = y - x**2 / 2
    A = x**2 - 2 * y
    polynomial_xy = sp.expand(polynomial_xt.subs(t, t_xy))
    require(
        sp.expand(jacobian(A, polynomial_xy, x, y)
                  - 2 * sp.diff(polynomial_xt, x).subs(t, t_xy)) == 0,
        "graph coordinate chain rule",
    )
    derivative = sp.Poly(2 * sp.diff(polynomial_xt, x) - kappa, x, t)
    equations = tuple(derivative.coeffs())
    solution = sp.solve(equations, tuple(coefficients.values()), dict=True)
    require(len(solution) == 1, "generic coordinate coefficient solution")
    solved = sp.expand(polynomial_xt.subs(solution[0]) - kappa * x / 2)
    require(sp.diff(solved, x) == 0, "generic survivor depends only on t")
    free_t_coefficients = tuple(
        coefficient for (i, _), coefficient in coefficients.items()
        if i == 0
    )
    require(
        set(solved.free_symbols) <= set(free_t_coefficients) | {t},
        "generic survivor free coefficient census",
    )

    # General nonlinear polynomial-coordinate control U=x+y^2, S=y.
    U = x + y**2
    S = y
    delta = jacobian(U, S, x, y)
    require(delta == 1, "nonlinear coordinate Jacobian")
    G = U**4 - 3 * U**2 + 5
    V = 7 * S + G
    require(jacobian(U, V, x, y) == 7, "nonlinear coordinate target Jacobian")
    u, v = sp.symbols("u v")
    inverse_y = (v - (u**4 - 3 * u**2 + 5)) / 7
    inverse_x = u - inverse_y**2
    require(
        sp.expand(U.subs({x: inverse_x, y: inverse_y}) - u) == 0,
        "nonlinear coordinate inverse first component",
    )
    require(
        sp.expand(V.subs({x: inverse_x, y: inverse_y}) - v) == 0,
        "nonlinear coordinate inverse second component",
    )

    # THM-2696 mixed shear: H(A,d)=-d^2+2d+A^3+2A.
    z = sp.symbols("z")
    B = y**2 - 2 * x * z
    H = -z**2 + 2 * z + A**3 + 2 * A
    sections = (-x - y, y - x + 2)
    target_rows = []
    for section in sections:
        mixed_target = sp.expand((B + H).subs(z, section))
        require(jacobian(A, mixed_target, x, y) == -4,
                "mixed shear constant Jacobian")
        triangular = sp.expand(-2 * x - 2 * t + (-2 * t)**3 + 2 * (-2 * t))
        require(
            sp.expand(mixed_target.subs(y, t + x**2 / 2) - triangular) == 0,
            "mixed shear triangular target",
        )
        target_rows.append(mixed_target)
    require(target_rows[0] == target_rows[1], "two sections share target")

    mixed_u, mixed_v = sp.symbols("mixed_u mixed_v")
    inverse_t = -mixed_u / 2
    inverse_x = sp.expand(
        -(mixed_v - ((-2 * inverse_t)**3 + 2 * (-2 * inverse_t))
          + 2 * inverse_t) / 2
    )
    inverse_y = sp.expand(inverse_t + inverse_x**2 / 2)
    require(
        sp.expand(A.subs({x: inverse_x, y: inverse_y}) - mixed_u) == 0,
        "mixed shear inverse first component",
    )
    require(
        sp.expand(target_rows[0].subs({x: inverse_x, y: inverse_y}) - mixed_v) == 0,
        "mixed shear inverse second component",
    )

    # At kappa=0, V=G(A) has zero Jacobian and is not a coordinate pair.
    zero_response = A**3 + A
    require(jacobian(A, zero_response, x, y) == 0, "zero-response boundary")

    # Characteristic three hostile: d/dy(y+y^3)=1, but its degree prevents a
    # polynomial inverse.  This is only an arithmetic control; the theorem is
    # stated in characteristic zero.
    yy = sp.symbols("yy")
    char_three_derivative = sp.Poly(sp.diff(yy + yy**3, yy), yy, modulus=3)
    require(char_three_derivative.as_expr() == 1, "characteristic-three Jacobian hostile")
    require(sp.degree(yy + yy**3, yy) == 3, "characteristic-three hostile degree")

    print("THM-2740 POLYNOMIAL-COORDINATE FIRST-TARGET AUDIT")
    print(f"generic_degree4_coefficients={len(coefficients)} free_G_coefficients={len(free_t_coefficients)}")
    print("coordinate_identity=Jac(A,V)=2*dV/dx_at_fixed_t")
    print("generic_solution=V=kappa*x/2+G(t)")
    print("nonlinear_coordinate_control=U=x+y^2 kappa=7 inverse=PASS")
    print("mixed_H=-d^2+2d+A^3+2A")
    print("mixed_sections=(-x-y,y-x+2) common_target=YES kappa=-4")
    print("mixed_triangular_inverse=PASS")
    print("zero_response_boundary=PASS")
    print("characteristic3_hostile=(x,y+y^3) Jacobian=1 degree=3")
    print("nongraph_and_noncoordinate_first_targets_not_tested")
    print("PASS")


if __name__ == "__main__":
    main()
