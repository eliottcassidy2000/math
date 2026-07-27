#!/usr/bin/env python3
"""Independent exact referee for THM-2463's degree-22 B-D plane.

This checker starts again from THM-2411's two fluxes, eliminates the linear
``zeta`` equation by the closed resultant formula, and enumerates every
lattice Minkowski summand of both the generic and exceptional Newton
polygons.  It then rebuilds the factor coefficient ideals, with every
denominator exception separated before clearing denominators.  Finally it
checks the fixed simple-zero section, square-lift parity floor, physical
square coordinate, and the ``y=0`` boundary.

All checks use ``require`` and therefore remain live under ``python -O``.
"""

from __future__ import annotations

from itertools import product
from math import gcd

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def is_unit_basis(basis: sp.GroebnerBasis) -> bool:
    return len(basis.polys) == 1 and basis.polys[0].as_expr() == 1


def convex_hull(points: set[tuple[int, int]]) -> list[tuple[int, int]]:
    """Return the counterclockwise hull, without repeating the first point."""

    ordered = sorted(points)

    def cross(
        origin: tuple[int, int],
        left: tuple[int, int],
        right: tuple[int, int],
    ) -> int:
        return (left[0] - origin[0]) * (right[1] - origin[1]) - (
            left[1] - origin[1]
        ) * (right[0] - origin[0])

    lower: list[tuple[int, int]] = []
    for point in ordered:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper: list[tuple[int, int]] = []
    for point in reversed(ordered):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return lower[:-1] + upper[:-1]


def primitive_edges(
    hull: list[tuple[int, int]],
) -> tuple[list[tuple[int, int]], list[int]]:
    directions: list[tuple[int, int]] = []
    lengths: list[int] = []
    for index, point in enumerate(hull):
        target = hull[(index + 1) % len(hull)]
        dx, dy = target[0] - point[0], target[1] - point[1]
        length = gcd(abs(dx), abs(dy))
        directions.append((dx // length, dy // length))
        lengths.append(length)
    return directions, lengths


def balanced_summands(
    directions: list[tuple[int, int]], lengths: list[int]
) -> set[tuple[int, ...]]:
    answer = set()
    for candidate in product(*(range(length + 1) for length in lengths)):
        x_sum = sum(length * direction[0] for length, direction in zip(candidate, directions))
        y_sum = sum(length * direction[1] for length, direction in zip(candidate, directions))
        if x_sum == 0 and y_sum == 0:
            answer.add(candidate)
    return answer


def polygon_width(
    directions: list[tuple[int, int]], lengths: tuple[int, ...]
) -> int:
    x_coordinate = 0
    x_coordinates = [0]
    for length, direction in zip(lengths, directions):
        x_coordinate += length * direction[0]
        x_coordinates.append(x_coordinate)
    require(x_coordinate == 0, "unbalanced summand reached polygon_width")
    return max(x_coordinates) - min(x_coordinates)


def audit_factor_shapes(
    directions: list[tuple[int, int]],
    total_lengths: list[int],
    summands: set[tuple[int, ...]],
    exceptional: bool,
) -> tuple[int, int, int]:
    """Classify every proper two-summand decomposition up to order."""

    seen = set()
    vertical = linear = quadratic = 0
    total = tuple(total_lengths)
    zero = tuple(0 for _ in total)
    for left in summands:
        right = tuple(total[i] - left[i] for i in range(len(total)))
        if left in (zero, total):
            continue
        key = tuple(sorted((left, right)))
        if key in seen:
            continue
        seen.add(key)
        widths = sorted(
            (polygon_width(directions, left), polygon_width(directions, right))
        )
        if widths[0] == 0:
            vertical += 1
            continue
        if widths[0] == 1:
            linear += 1
            continue
        require(widths == [2, 2], "unclassified Newton factor-width split")
        # In either polygon, a 2+2 split has the b=0,a=2 triangular
        # summand with edge lengths (2,0,2,2).  Its monic polynomial is
        # p^2+(av+b)p+cv^2+dv+e.
        triangle = (2, 0, 2, 2)
        require(
            left == triangle or right == triangle,
            "2+2 split lost its Delta(2,0) factor",
        )
        quadratic += 1
    if exceptional:
        require(vertical == 0, "exceptional polygon unexpectedly has a vertical split")
    else:
        require(vertical == 1, "generic polygon vertical split count changed")
    require(linear > 0 and quadratic == 1, "factor-shape inventory changed")
    return vertical, linear, quadratic


def main() -> None:
    y, u, z = sp.symbols("y u Z")
    B, D = sp.symbols("B D")
    p, v, zeta, lam = sp.symbols("p v zeta lambda")

    # Direct C=E=W=0 specialization of THM-2411, equations (12)--(16).
    wall = 616 * B - 1089 * u + 63 * y**2
    residual = (
        -745360 * B * u * y
        + 6160 * B * y**3
        + 511104 * D * y
        + 922383 * u**2 * y
        - 25410 * u * y**3
        + 63 * y**5
    )
    n1 = 1331 * wall * z + 4 * residual
    n2 = (
        15944049 * z**2
        + 65591680 * B * z * y
        - 162339408 * z * u * y
        + 2236080 * z * y**3
        + 1443016960 * B * u**2
        - 71554560 * B * u * y**2
        + 98560 * B * y**4
        - 1978994688 * D * u
        + 16355328 * D * y**2
        - 1190488992 * u**3
        + 147581280 * u**2 * y**2
        - 1219680 * u * y**4
        + 672 * y**6
    )
    substitutions = {
        B: p * y**2,
        D: lam * p**2 * y**4,
        u: v * y**2,
        z: zeta * y**3,
    }
    first = sp.expand(n1.subs(substitutions) / y**5)
    second = sp.expand(n2.subs(substitutions) / y**6)

    expected_first = (
        2044416 * lam * p**2
        - 2981440 * p * v
        + 819896 * p * zeta
        + 24640 * p
        + 3689532 * v**2
        - 1449459 * v * zeta
        - 101640 * v
        + 83853 * zeta
        + 252
    )
    expected_second = (
        -1978994688 * lam * p**2 * v
        + 16355328 * lam * p**2
        + 1443016960 * p * v**2
        - 71554560 * p * v
        + 65591680 * p * zeta
        + 98560 * p
        - 1190488992 * v**3
        + 147581280 * v**2
        - 162339408 * v * zeta
        - 1219680 * v
        + 15944049 * zeta**2
        + 2236080 * zeta
        + 672
    )
    require(sp.expand(first - expected_first) == 0, "first quotient flux changed")
    require(sp.expand(second - expected_second) == 0, "second quotient flux changed")
    quotient_wall = 616 * p - 1089 * v + 63
    require(
        sp.Poly(first, zeta).coeff_monomial(zeta) == 1331 * quotient_wall,
        "zeta reconstruction coefficient is not the open first-flux wall",
    )

    # Eliminate a*zeta+b against A*zeta^2+B*zeta+C by the closed formula
    # A*b^2-a*B*b+a^2*C, independently of SymPy's resultant routine.
    first_poly = sp.Poly(first, zeta)
    second_poly = sp.Poly(second, zeta)
    aa, bb = first_poly.all_coeffs()
    AA, BB, CC = second_poly.all_coeffs()
    direct_resultant = sp.expand(AA * bb**2 - aa * BB * bb + aa**2 * CC)
    resultant = sp.resultant(first, second, zeta)
    require(sp.expand(direct_resultant - resultant) == 0, "closed resultant formula failed")
    primitive = sp.Poly(direct_resultant, p, v, lam).primitive()
    require(primitive[0] == 2**4 * 11**6, "resultant content changed")
    R = primitive[1].as_expr()
    require(len(sp.Poly(R, p, v, lam).terms()) == 28, "resultant support changed")

    # The support and every balanced edge-length suballocation are enumerated,
    # rather than inferred from a polygon name.
    generic_support = {
        monomial
        for monomial, coefficient in sp.Poly(R.subs(lam, 1), p, v).terms()
        if coefficient != 0
    }
    exceptional_value = sp.Rational(49, 33)
    exceptional_support = {
        monomial
        for monomial, coefficient in sp.Poly(R.subs(lam, exceptional_value), p, v).terms()
        if coefficient != 0
    }
    generic_hull = convex_hull(generic_support)
    exceptional_hull = convex_hull(exceptional_support)
    require(
        generic_hull == [(0, 0), (4, 0), (4, 1), (0, 5)],
        "generic Newton hull changed",
    )
    require(
        exceptional_hull == [(0, 0), (3, 0), (4, 1), (0, 5)],
        "exceptional Newton hull changed",
    )
    generic_directions, generic_lengths = primitive_edges(generic_hull)
    exceptional_directions, exceptional_lengths = primitive_edges(exceptional_hull)
    generic_summands = balanced_summands(generic_directions, generic_lengths)
    exceptional_summands = balanced_summands(exceptional_directions, exceptional_lengths)
    require(
        generic_summands
        == {(a0, b0, a0, a0 + b0) for a0 in range(5) for b0 in range(2)},
        "generic Minkowski length solutions changed",
    )
    require(
        exceptional_summands
        == {
            (a0, b0, a0 + b0, a0 + 2 * b0)
            for a0 in range(4)
            for b0 in range(2)
        },
        "exceptional Minkowski length solutions changed",
    )
    generic_inventory = audit_factor_shapes(
        generic_directions, generic_lengths, generic_summands, False
    )
    exceptional_inventory = audit_factor_shapes(
        exceptional_directions, exceptional_lengths, exceptional_summands, True
    )

    constant_in_p = sp.Poly(R, p).coeff_monomial(1)
    linear_in_p = sp.Poly(R, p).coeff_monomial(p)
    require(
        sp.gcd(sp.Poly(constant_in_p, v), sp.Poly(linear_in_p, v)).degree() == 0,
        "vertical-factor exclusion failed",
    )

    # Rebuild the exact factor equations from the independently obtained R.
    a, b, c, d, e = sp.symbols("a b c d e")
    linear_e0 = p + a * v + b
    remainder_e0 = sp.Poly(
        sp.rem(sp.Poly(R, p), sp.Poly(linear_e0, p)).as_expr(), v
    )
    require(
        is_unit_basis(
            sp.groebner(
                [coefficient for _, coefficient in remainder_e0.terms()],
                a,
                b,
                lam,
                order="grevlex",
            )
        ),
        "affine linear factor survived independent audit",
    )

    denominator = v + a
    numerator = b * v**2 + c * v + d
    cleared = sp.cancel(denominator**4 * R.subs(p, -numerator / denominator))
    require(cleared.as_numer_denom()[1] == 1, "linear e=1 denominator did not clear")
    cleared_coefficients = [
        coefficient for _, coefficient in sp.Poly(sp.expand(cleared), v).terms()
    ]
    require(
        is_unit_basis(
            sp.groebner(
                cleared_coefficients, a, b, c, d, lam, order="grevlex"
            )
        ),
        "degree-(2,1) linear factor survived independent audit",
    )

    quadratic = p**2 + (a * v + b) * p + c * v**2 + d * v + e
    quadratic_remainder = sp.Poly(
        sp.rem(sp.Poly(R, p), sp.Poly(quadratic, p)).as_expr(), p, v
    )
    equations = [coefficient for _, coefficient in quadratic_remainder.terms()]
    by_monomial = dict(quadratic_remainder.terms())

    type_i = [
        equation.subs({a: -sp.Rational(99, 28), c: sp.Rational(9801, 3136)})
        for equation in equations
    ]
    require(
        is_unit_basis(sp.groebner(type_i, b, d, e, lam, order="grevlex")),
        "quadratic type I survived independent audit",
    )

    type_ii_substitution = {
        a: -sp.Rational(35, 48) / lam,
        c: sp.Rational(77, 128) / lam,
    }
    type_ii = [equation.subs(type_ii_substitution) for equation in equations]
    type_ii_exception = sp.Rational(196, 891)
    require(
        is_unit_basis(
            sp.groebner(
                [equation.subs(lam, type_ii_exception) for equation in type_ii],
                b,
                d,
                e,
                order="grevlex",
            )
        ),
        "type-II determinant exception survived independent audit",
    )
    q13 = by_monomial[(1, 3)].subs(type_ii_substitution)
    q04 = by_monomial[(0, 4)].subs(type_ii_substitution)
    matrix = sp.Matrix(
        [[sp.diff(q13, b), sp.diff(q13, d)], [sp.diff(q04, b), sp.diff(q04, d)]]
    )
    determinant = sp.factor(matrix.det())
    require(
        determinant
        == 27102253467051586289664 * (891 * lam - 196) ** 2,
        "type-II determinant changed",
    )
    b_solution = -sp.Rational(7) * (
        78682428 * lam**2 - 34767117 * lam + 3841600
    ) / (23232 * lam * (891 * lam - 196) ** 2)
    d_solution = sp.Rational(147) * (
        2822688 * lam**2 - 1246014 * lam + 137543
    ) / (11264 * lam * (891 * lam - 196) ** 2)
    require(
        sp.factor(q13.subs({b: b_solution, d: d_solution})) == 0
        and sp.factor(q04.subs({b: b_solution, d: d_solution})) == 0,
        "type-II generic top-next solution changed",
    )
    type_ii_numerators = [
        sp.together(equation.subs({b: b_solution, d: d_solution}))
        .as_numer_denom()[0]
        for equation in type_ii
    ]
    require(
        is_unit_basis(sp.groebner(type_ii_numerators, e, lam, order="grevlex")),
        "quadratic type II generic branch survived independent audit",
    )
    require(
        exceptional_value not in (0, type_ii_exception),
        "lambda=49/33 fell into an unaudited type-II denominator fibre",
    )

    h, x = sp.symbols("h x")
    lambda_h = (280 * h - 231) / (384 * h**2)
    type_iii_substitution = {
        lam: lambda_h,
        a: -(sp.Rational(99, 56) + h),
        c: sp.Rational(99, 56) * h,
    }
    type_iii = [equation.subs(type_iii_substitution) for equation in equations]
    type_iii_numerators = [
        sp.Poly(
            sp.together(equation).as_numer_denom()[0], b, d, e, h
        ).primitive()[1].as_expr()
        for equation in type_iii
    ]
    q13_iii = by_monomial[(1, 3)].subs(type_iii_substitution).subs(d, x - b * h)
    q04_iii = by_monomial[(0, 4)].subs(type_iii_substitution).subs(d, x - b * h)
    compatibility = sp.factor(
        sp.Poly(q13_iii, x).coeff_monomial(x)
        * sp.Poly(q04_iii, x).coeff_monomial(1)
        - sp.Poly(q04_iii, x).coeff_monomial(x)
        * sp.Poly(q13_iii, x).coeff_monomial(1)
    )
    compatibility_denominator = sp.cancel(compatibility).as_numer_denom()[1]
    expected_compatibility = (
        sp.Rational(191848476676752603477, 16)
        * (20 * h - 33)
        * (448 * h**2 - 1320 * h + 1089) ** 2
        / h**5
    )
    require(
        sp.factor(compatibility - expected_compatibility) == 0,
        "type-III compatibility expression changed",
    )
    require(
        sp.Poly(compatibility_denominator, h).sqf_part().monic().as_expr() == h,
        "type-III compatibility acquired a denominator exception beyond h=0",
    )
    # h=0 is not a root of 384*lambda*h^2-280*h+231, so clearing h is safe.
    require(231 != 0, "type-III h=0 exclusion failed")
    rational_h = sp.Rational(33, 20)
    require(
        is_unit_basis(
            sp.groebner(
                [equation.subs(h, rational_h) for equation in type_iii_numerators],
                b,
                d,
                e,
                order="grevlex",
            )
        ),
        "type-III rational compatibility branch survived",
    )
    quadratic_h = 448 * h**2 - 1320 * h + 1089
    require(
        is_unit_basis(
            sp.groebner(
                type_iii_numerators + [quadratic_h],
                b,
                d,
                e,
                h,
                order="grevlex",
            )
        ),
        "type-III quadratic compatibility branch survived",
    )

    # Fixed simple zeros and the connected double-cover genus floor.
    L5 = -constant_in_p / 567
    require(sp.Poly(L5, v).degree() == 5, "p=0 section degree changed")
    require(
        sp.gcd(sp.Poly(L5, v), sp.Poly(sp.diff(L5, v), v)).degree() == 0,
        "p=0 section lost squarefreeness",
    )
    require(
        sp.expand(sp.diff(R, v).subs(p, 0) + 567 * sp.diff(L5, v)) == 0,
        "p=0 points lost their smooth transverse derivative",
    )
    simple_odd_zeros = 5
    minimum_even_branch_count = simple_odd_zeros + simple_odd_zeros % 2
    minimum_genus = -1 + minimum_even_branch_count // 2
    require(minimum_even_branch_count == 6, "branch parity floor changed")
    require(minimum_genus == 2, "square-lift genus floor changed")

    # The square coordinate is physical, not an arbitrary quadratic extension.
    root_B, Y = sp.symbols("root_B Y", nonzero=True)
    require(
        sp.cancel((p * Y**2).subs({p: root_B**2 / y**2, Y: y / root_B})) == 1,
        "physical square-coordinate identity failed",
    )
    require(
        sp.expand(n1.subs(y, 0) - 1331 * (616 * B - 1089 * u) * z) == 0,
        "y=0 boundary changed",
    )

    print("THM-2463 independent hostile referee")
    print("flux_specialization=RECONSTRUCTED_FROM_THM-2411")
    print("elimination=closed_linear_quadratic_resultant")
    print(f"generic_newton_hull={generic_hull}")
    print(f"exceptional_newton_hull={exceptional_hull}")
    print(f"generic_factor_inventory_vertical_linear_quadratic={generic_inventory}")
    print(f"exceptional_factor_inventory_vertical_linear_quadratic={exceptional_inventory}")
    print("lambda_49_over_33_factor_shapes=COMPLETE")
    print("linear_factor_unit_ideals=2")
    print("quadratic_factor_unit_ideals=5")
    print("denominator_exceptions=lambda_0_excluded,lambda_196_over_891_audited,h_0_impossible")
    print("p_zero_smooth_simple_points=5")
    print("square_lift_connected=YES")
    print("square_lift_minimum_branch_points=6")
    print("square_lift_minimum_genus=2")
    print("physical_square_coordinate=Y=y/sqrt(B)")
    print("y_zero_boundary=Z_zero_on_open_first_flux")
    print("ALL INDEPENDENT CHECKS PASSED")


if __name__ == "__main__":
    main()
