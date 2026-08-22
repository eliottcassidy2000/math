#!/usr/bin/env python3
"""Exact companion for THM-3685.

The calculation keeps the native first two coordinates of the explicit
degree-four weighted Keller map from THM-3438.  It checks

* polynomial cancellation and the ambient determinant ``-6``;
* the top-Jacobian obstruction for general homogeneous graph initials;
* several full graph restrictions, including the graph through the known
  collision pair; and
* the direct affine-target-plane discriminant hostile.

Every truth gate uses ``require`` and therefore remains active under
``python -O``.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def total_initial(poly: sp.Expr, variables: tuple[sp.Symbol, ...]) -> tuple[int, sp.Expr]:
    """Return ordinary total degree and its homogeneous initial form."""
    expanded = sp.Poly(sp.expand(poly), *variables)
    degree = expanded.total_degree()
    initial = sum(
        coefficient * sp.prod(variable**power for variable, power in zip(variables, monomial))
        for monomial, coefficient in expanded.terms()
        if sum(monomial) == degree
    )
    return degree, sp.expand(initial)


def main() -> None:
    x, y, z = sp.symbols("x y z")
    u = 1 + 3 * x * y
    gamma = 1 - 4 * x * y - x**2 * z

    A = sp.cancel((2 * u + u**2 - 3 * u**4 * gamma**2) / x**2)
    B = sp.cancel((1 + u - 2 * u**3 * gamma**2) / x)
    C = x * gamma
    require(sp.denom(A) == 1 and sp.denom(B) == 1, "weighted quotients did not cancel")
    A, B = sp.expand(A), sp.expand(B)

    ambient = sp.factor(
        sp.det(sp.Matrix([[sp.diff(F, variable) for variable in (x, y, z)] for F in (A, B, C)]))
    )
    require(ambient == -6, ("ambient determinant", ambient))

    top_rows: list[tuple[object, ...]] = []

    # Constant graphs, including h=0.  The parameter c is coefficient data,
    # not one of the variables used for ordinary total degree.
    c = sp.symbols("c")
    h0 = c
    A0 = sp.expand(A.subs(z, h0))
    B0 = sp.expand(B.subs(z, h0))
    J0 = sp.expand(sp.diff(A0, x) * sp.diff(B0, y) - sp.diff(A0, y) * sp.diff(B0, x))
    degree0, initial0 = total_initial(J0, (x, y))
    W0 = x**4 * y**3 * (4 * y + c * x) ** 2
    predicted0 = sp.expand(-13122 * W0 * sp.diff(W0, x))
    require((degree0, sp.expand(initial0 - predicted0)) == (17, 0), (degree0, sp.factor(initial0 - predicted0)))
    require(predicted0 != 0, "constant-graph initial vanished")
    top_rows.append(("constant", degree0, sp.factor(predicted0)))

    # General homogeneous initials in degrees 1--5.  Lower graph terms cannot
    # affect the displayed top form, so this checks the structural formula in
    # a coefficient-generic polynomial ring.
    for degree in range(1, 6):
        coefficients = sp.symbols(f"a0:{degree + 1}")
        H = sum(coefficients[i] * x**i * y ** (degree - i) for i in range(degree + 1))
        gamma_top = -x**2 * H
        u_top = 3 * x * y
        A_top = sp.expand(-3 * u_top**4 * gamma_top**2 / x**2)
        B_top = sp.expand(-2 * u_top**3 * gamma_top**2 / x)
        W = x**6 * y**3 * H**2
        require(sp.expand(A_top + 243 * y * W) == 0, ("A top", degree))
        require(sp.expand(B_top + 54 * W) == 0, ("B top", degree))
        J_top = sp.expand(sp.diff(A_top, x) * sp.diff(B_top, y) - sp.diff(A_top, y) * sp.diff(B_top, x))
        predicted = sp.expand(-13122 * W * sp.diff(W, x))
        require(sp.expand(J_top - predicted) == 0, ("Jacobian top", degree))
        require(predicted.subs({coefficient: 1 for coefficient in coefficients}) != 0, ("nonzero top", degree))
        top_rows.append((degree, 4 * degree + 17, len(sp.Poly(predicted, x, y, *coefficients).terms())))

    # Full restrictions exercise lower-order cancellation independently of the
    # homogeneous proof.  The predicted degree is 17 for constants and
    # 4 deg(h)+17 otherwise.
    graph_rows: list[tuple[object, ...]] = []
    graphs = (
        sp.Integer(0),
        sp.Integer(3),
        1 - x,
        x + 2 * y,
        x**2 + x * y + 2 * y**2 + x - 3,
        x**3 - 2 * x * y**2 + y + 1,
    )
    for h in graphs:
        Ah = sp.expand(A.subs(z, h))
        Bh = sp.expand(B.subs(z, h))
        Jh = sp.expand(sp.diff(Ah, x) * sp.diff(Bh, y) - sp.diff(Ah, y) * sp.diff(Bh, x))
        graph_degree = sp.Poly(h, x, y).total_degree()
        expected_degree = 17 if graph_degree == 0 else 4 * graph_degree + 17
        actual_degree, actual_initial = total_initial(Jh, (x, y))
        require(actual_degree == expected_degree and actual_initial != 0, (h, actual_degree, expected_degree))
        graph_rows.append((str(h), actual_degree, len(sp.Poly(Jh, x, y).terms())))

    collision_graph = 1 - x
    collision_points = ((1, 0, 0), (-1, 0, 2))
    collision_images = tuple(tuple(sp.expand(F).subs(dict(zip((x, y, z), point))) for F in (A, B, C)) for point in collision_points)
    require(collision_images == ((0, 0, 1), (0, 0, 1)), collision_images)
    require(tuple(collision_graph.subs({x: point[0], y: point[1]}) for point in collision_points) == (0, 2), "collision graph")

    # Direct affine target-plane control.  If C-aA-bB-c had a graph factor,
    # its quadratic-in-z discriminant would be a square.  In the localization
    # t=xy the discriminant coefficients force a=b=0; the residual C-c still
    # has no polynomial graph factor.
    a, b, c_aff, t = sp.symbols("a b c_aff t")
    plane_pullback = sp.expand(C - a * A - b * B - c_aff)
    discriminant = sp.factor(sp.discriminant(plane_pullback, z))
    discriminant_xt = sp.Poly(sp.expand(discriminant.subs(y, t / x)), x)
    expected_coefficients = (
        36 * a**2 * (t + 1) * (3 * t + 1) ** 5,
        12 * a * b * (3 * t + 1) ** 4 * (5 * t + 4),
        4 * (3 * t + 1) ** 3 * (9 * a * c_aff * t + 3 * a * c_aff + 6 * b**2 * t + 4 * b**2),
        8 * b * c_aff * (3 * t + 1) ** 3,
        sp.Integer(0),
        sp.Integer(0),
        sp.Integer(1),
    )
    actual_coefficients = tuple(sp.factor(discriminant_xt.coeff_monomial(x**power)) for power in range(7))
    require(actual_coefficients == expected_coefficients, actual_coefficients)

    print("theorem=THM-3685-weighted-quartic-native-polynomial-graph-descent-no-go")
    print("ambient_map=THM3438_quartic_G;determinant=-6;generic_degree=4")
    print(f"top_rows={tuple(top_rows)}")
    print(f"full_graph_rows={tuple(graph_rows)}")
    print(f"collision_graph=h=1-x;points={collision_points};common_image={collision_images[0]}")
    print("top_mechanism=A_top=(9/2)yB_top;B_top=-54W;Jac_top=-13122W*d_xW_nonzero")
    print("constant_W=x^4*y^3*(4y+c*x)^2;positive_degree_W=x^6*y^3*h_d^2")
    print(f"affine_plane_discriminant_coefficients={actual_coefficients}")
    print("consequence=no_native_source_graph_z=h_maps_to_any_target_graph_C=g(A,B)")
    print("scope=native_A,B_projection_only;no_arbitrary_target_coordinate_change;no_quartic_C3_exclusion;no_JC2")
    print("commands=python3 -B 04-computation/jacobian_weighted_quartic_native_graph_no_go_thm3685.py;python3 -B -O 04-computation/jacobian_weighted_quartic_native_graph_no_go_thm3685.py")


if __name__ == "__main__":
    main()
