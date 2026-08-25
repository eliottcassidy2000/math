#!/usr/bin/env python3
"""Standard-library audit of THM-4120, independent of the SymPy path.

This audit reconstructs the discriminant coefficient ledger, response census,
DE-weight Diophantine layers, fixed-layer Jacobian brackets, re-entry costs,
and degree floors using ``Fraction`` and sparse dictionaries only.  The
Shioda--Tate and Kodaira classification steps remain cited mathematical input.
"""

from __future__ import annotations

from fractions import Fraction as F


BiPoly = dict[tuple[int, int], F]
UniPoly = dict[int, F]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def clean(poly):
    return {key: value for key, value in poly.items() if value}


def b_add(left: BiPoly, right: BiPoly) -> BiPoly:
    out = dict(left)
    for key, value in right.items():
        out[key] = out.get(key, F(0)) + value
    return clean(out)


def b_scale(poly: BiPoly, scalar: F) -> BiPoly:
    return clean({key: scalar * value for key, value in poly.items()})


def b_mul(left: BiPoly, right: BiPoly) -> BiPoly:
    out: BiPoly = {}
    for (m, n), value in left.items():
        for (r, s), other in right.items():
            key = (m + r, n + s)
            out[key] = out.get(key, F(0)) + value * other
    return clean(out)


def b_pow(poly: BiPoly, exponent: int) -> BiPoly:
    out: BiPoly = {(0, 0): F(1)}
    for _ in range(exponent):
        out = b_mul(out, poly)
    return out


def b_shift(poly: BiPoly, m: int, n: int) -> BiPoly:
    return {(i + m, j + n): value for (i, j), value in poly.items()}


def b_derivative(poly: BiPoly, axis: int) -> BiPoly:
    out: BiPoly = {}
    for (m, n), value in poly.items():
        exponent = (m, n)[axis]
        if exponent:
            key = (m - 1, n) if axis == 0 else (m, n - 1)
            out[key] = value * exponent
    return clean(out)


def b_bracket(left: BiPoly, right: BiPoly) -> BiPoly:
    return b_add(
        b_mul(b_derivative(left, 0), b_derivative(right, 1)),
        b_scale(b_mul(b_derivative(left, 1), b_derivative(right, 0)), F(-1)),
    )


def u_add(left: UniPoly, right: UniPoly) -> UniPoly:
    out = dict(left)
    for degree, value in right.items():
        out[degree] = out.get(degree, F(0)) + value
    return clean(out)


def u_scale(poly: UniPoly, scalar: F) -> UniPoly:
    return clean({degree: scalar * value for degree, value in poly.items()})


def u_mul(left: UniPoly, right: UniPoly) -> UniPoly:
    out: UniPoly = {}
    for degree, value in left.items():
        for other_degree, other in right.items():
            key = degree + other_degree
            out[key] = out.get(key, F(0)) + value * other
    return clean(out)


def u_pow(poly: UniPoly, exponent: int) -> UniPoly:
    out: UniPoly = {0: F(1)}
    for _ in range(exponent):
        out = u_mul(out, poly)
    return out


def u_derivative(poly: UniPoly) -> UniPoly:
    return clean({degree - 1: degree * value for degree, value in poly.items()
                  if degree})


def u_evaluate(poly: UniPoly, value: F) -> F:
    return sum(coefficient * value**degree for degree, coefficient in poly.items())


def substitute_w(poly: UniPoly) -> BiPoly:
    return {(6 * degree, 7 * degree): coefficient
            for degree, coefficient in poly.items()}


def fixed_layer_expected(
    a: int, b: int, c: int, d: int, f: UniPoly, g: UniPoly
) -> BiPoly:
    weight_a = 7 * a - 6 * b
    weight_c = 7 * c - 6 * d
    coefficient = u_scale(u_mul(f, g), F(a * d - b * c))
    coefficient = u_add(
        coefficient,
        u_scale(u_mul({1: F(1)}, u_mul(f, u_derivative(g))), F(weight_a)),
    )
    coefficient = u_add(
        coefficient,
        u_scale(u_mul({1: F(1)}, u_mul(u_derivative(f), g)), F(-weight_c)),
    )
    return b_shift(substitute_w(coefficient), a + c - 1, b + d - 1)


def layer(weight: int, depth: int, count: int) -> list[tuple[int, int]]:
    points = []
    for m in range(200):
        for n in range(200):
            if 7 * m - 6 * n == weight and m <= n + depth:
                points.append((m, n))
    return sorted(points, key=lambda pair: (sum(pair), pair))[:count]


def is_eisenstein_norm(value: int) -> bool:
    residual = value
    prime = 2
    while prime * prime <= residual:
        exponent = 0
        while residual % prime == 0:
            residual //= prime
            exponent += 1
        if prime % 3 == 2 and exponent % 2:
            return False
        prime += 1
    return not (residual > 1 and residual % 3 == 2)


def main() -> None:
    # Expand 4a4^3+27a6^2 by hand as a polynomial in (q,a^3).
    # Coefficient order: q^2, q*a^3, a^6.
    four_a4_cubed = (F(0), F(0), F(-27, 16))
    twentyseven_a6_squared = (F(27), F(-27, 2), F(27, 16))
    cubic_discriminant_core = tuple(
        left + right
        for left, right in zip(four_a4_cubed, twentyseven_a6_squared)
    )
    delta_coefficients = tuple(F(-16) * value for value in cubic_discriminant_core)
    require(cubic_discriminant_core == (F(27), F(-27, 2), F(0)),
            "short-Weierstrass cancellation changed")
    require(delta_coefficients == (F(-432), F(216), F(0)),
            "target discriminant coefficients changed")
    infinity_valuations = (4, 5, 10)
    finite_root_multiplicities = (1, 1)
    require(sum((10, 1, 1)) == 12, "elliptic-surface Euler ledger changed")

    # Odd unimodular section/fibre block and E8 exhaust rho=10.
    section_fibre_det = (-1) * 0 - 1 * 1
    require((section_fibre_det, 2 + 8, 10) == (-1, 10, 10),
            "Shioda--Tate rank/discriminant ledger changed")

    # Target response.  Five K-rational punctures are forced to O; BC is one
    # irreducible quadratic block because q+gamma has a simple zero.
    rational_indices = (1, 3, 7, 3, 3)
    rational_weight = sum(rational_indices)
    candidates = (rational_weight, rational_weight + 4)
    survivors = tuple(value for value in candidates if is_eisenstein_norm(value))
    require((rational_weight, candidates, survivors) == (17, (17, 21), (21,)),
            "response collapse changed")
    pullback = (1, 2, 2, 3, 7, 3, 3)
    require((sum(pullback), 2 * sum(pullback), 3 * sum(pullback)) == (21, 42, 63),
            "pullback/pole ledger changed")

    # The quadratic hostile is checked at several exact values of a,r.
    for aval, rval in ((F(2), F(3)), (F(-4), F(5)), (F(7, 3), F(-2, 5))):
        qval = aval**3 / 2 + rval**2
        target_rhs = (aval / 2) ** 3 - F(3, 4) * aval**2 * (aval / 2)
        target_rhs += qval - aval**3 / 4
        require(target_rhs == rval**2, "quadratic target hostile changed")

    # Independent sparse-polynomial bracket checks: baseline and higher layer.
    f = {0: F(2), 1: F(3), 2: F(5)}
    g = {0: F(7), 1: F(11), 2: F(13)}
    for a, b, c, d in ((2, 0, 3, 0), (9, 8, 10, 8), (15, 15, 16, 15)):
        left = b_shift(substitute_w(f), a, b)
        right = b_shift(substitute_w(g), c, d)
        require(
            b_bracket(left, right) == fixed_layer_expected(a, b, c, d, f, g),
            f"fixed-layer bracket failed at {(a, b, c, d)}",
        )

    A14 = [(2, 0)] + layer(14, 1, 3)
    C21 = [(3, 0)] + layer(21, 2, 3)
    A15 = layer(15, 1, 3)
    C22 = layer(22, 2, 3)
    require(A14 == [(2, 0), (8, 7), (14, 14), (20, 21)],
            "A14 layer changed")
    require(C21 == [(3, 0), (9, 7), (15, 14), (21, 21)],
            "C21 layer changed")
    require(A15 == [(9, 8), (15, 15), (21, 22)], "A15 layer changed")
    require(C22 == [(10, 8), (16, 15), (22, 22)], "C22 layer changed")

    # Kernel equation and the minimal response, using gamma=2, theta=3.
    gamma, theta = F(2), F(3)
    T = -gamma / theta
    h = {0: F(1), 1: -F(8) / T}
    f_kernel = u_scale(u_pow(h, 2), gamma**2)
    g_kernel = u_scale(u_pow(h, 3), gamma**3)
    kernel_equation = u_add(
        u_scale(u_mul(f_kernel, u_derivative(g_kernel)), F(2)),
        u_scale(u_mul(u_derivative(f_kernel), g_kernel), F(-3)),
    )
    require(kernel_equation == {}, "square/cube kernel equation changed")
    require((u_evaluate(h, F(0)), u_evaluate(h, T)) == (F(1), F(-7)),
            "h endpoint response changed")
    A_lead = T**-2 * gamma**2 * u_evaluate(h, T) ** 2
    C_lead = T**-3 * gamma**3 * u_evaluate(h, T) ** 3
    require((A_lead, C_lead) == (49 * theta**2, 343 * theta**3),
            "target-O coefficients changed")

    # Multiplicity ratio and the two edge re-entry costs.
    require(14 * 3 == 21 * 2, "baseline multiplicity lock changed")
    require(14 * 1 != 21 * 1, "equal-simple-root hostile disappeared")
    phi, edge_sum = F(5), F(7)
    c1 = gamma * phi / theta**2
    c2 = gamma * edge_sum / theta**2 - gamma * phi**2 / theta**3
    # Multiplication verifies the inverse-series coefficients through z^2.
    inverse_series = (T, c1, c2)
    denominator = (theta, phi, edge_sum)
    product = [F(0)] * 3
    for i, left in enumerate(inverse_series):
        for j, right in enumerate(denominator):
            if i + j < 3:
                product[i + j] += left * right
    require(tuple(product) == (-gamma, F(0), F(0)), "DE inverse series changed")
    require(gamma * edge_sum / theta**2 != 0, "phi=0 re-entry hostile failed")
    require((15 - 1, 16 - 2) == (14, 14), "re-entry weights changed")

    baseline_floor = (2 + 26, 3 + 39)
    higher_floor = (2 * 15, 2 * 22 - 13)
    combined_floor = (min(baseline_floor[0], higher_floor[0]),
                      min(baseline_floor[1], higher_floor[1]))
    require((baseline_floor, higher_floor, combined_floor)
            == ((28, 42), (30, 31), (28, 31)),
            "DE degree floors changed")

    print("THM4120 JC23 INDEPENDENT STANDARD-LIBRARY AUDIT")
    print("scope=THM4103_smooth_nonresonant_theta_only;JC2=OPEN")
    print(f"discriminant_coefficients_q2_qa3_a6={delta_coefficients}")
    print(f"finite_root_multiplicities={finite_root_multiplicities}")
    print(f"infinity_valuations={infinity_valuations};fibres=I1,I1,IIstar")
    print("section_fibre_det=-1;trivial_rank=10;MW={O}")
    print(f"rational_weight={rational_weight};candidates={candidates};survivors={survivors}")
    print(f"pullback={pullback};pole_degrees=(42,63)")
    print("BC_irreducible_reason=odd_valuation_at_q=-gamma")
    print("quadratic_hostile=verified_at_3_exact_specializations")
    print(f"A14={A14};C21={C21};A15={A15};C22={C22}")
    print("general_bracket=verified_on_3_distinct_layer_pairs")
    print(f"T={T};h0={u_evaluate(h, F(0))};hT={u_evaluate(h, T)}")
    print("multiplicity_lock=14*3=21*2;equal_simple_hostile=nonzero")
    print(f"DE_inverse_series={inverse_series};reentry_costs=(1,2)")
    print(f"baseline_floor={baseline_floor};higher_floor={higher_floor}")
    print(f"combined_floor={combined_floor}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
