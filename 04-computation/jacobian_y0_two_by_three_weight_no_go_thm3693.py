#!/usr/bin/env python3
"""Exact factor companion for THM-3693.

THM-3693 is an all-degree support/centralizer/UFD proof.  This companion
checks the arithmetic support grammar and all four decisive ODE/factor
calculations (two possible constant buckets, with either endpoint parameter
equal to two).  It also checks the cube-root exceptional homogeneous terms.

Every truth gate uses ``require`` and therefore survives ``python -O``.
"""

from __future__ import annotations

from itertools import combinations
from math import gcd

import sympy as sp


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def weight_bracket(
    first_weight: int,
    first_coefficient: sp.Expr,
    second_weight: int,
    second_coefficient: sp.Expr,
    u: sp.Symbol,
) -> sp.Expr:
    """Coefficient in {x^r p(u),x^s q(u)} after removing x^(r+s+1)."""
    return sp.expand(
        second_weight * sp.diff(first_coefficient, u) * second_coefficient
        - first_weight * first_coefficient * sp.diff(second_coefficient, u)
    )


def require_polynomial_quotient(numerator: sp.Expr, factor: sp.Expr, detail: object) -> sp.Expr:
    quotient = sp.cancel(numerator / factor)
    require(sp.denom(quotient) == 1, (detail, sp.factor(numerator), sp.factor(factor), quotient))
    return sp.factor(quotient)


def main() -> None:
    u = sp.symbols("u")
    kappa = sp.symbols("kappa")
    H = sp.Function("H")(u)
    K = sp.Function("K")(u)
    T = sp.Function("T")(u)

    # Two-by-three sumset grammar on a hostile window.  A duplicated bucket
    # requires the P-gap to be one Q interval; non-arithmetic supports have
    # exactly one duplicate, while equal gaps have the two middle duplicates.
    support_rows = {"nonarithmetic": 0, "arithmetic": 0}
    for first_support in combinations(range(-6, 7), 2):
        first_gap = first_support[1] - first_support[0]
        for second_support in combinations(range(-6, 7), 3):
            first_q_gap = second_support[1] - second_support[0]
            second_q_gap = second_support[2] - second_support[1]
            buckets: dict[int, list[tuple[int, int]]] = {}
            for first_weight in first_support:
                for second_weight in second_support:
                    buckets.setdefault(first_weight + second_weight + 1, []).append((first_weight, second_weight))
            duplicated = [weight for weight, pairs in buckets.items() if len(pairs) > 1]
            if duplicated:
                require(first_gap in (first_q_gap, second_q_gap, first_q_gap + second_q_gap), (first_support, second_support, duplicated))
                if first_gap == first_q_gap == second_q_gap:
                    require(len(duplicated) == 2, ("arithmetic duplicate count", first_support, second_support, duplicated))
                    support_rows["arithmetic"] += 1
                else:
                    require(len(duplicated) == 1, ("nonarithmetic duplicate count", first_support, second_support, duplicated))
                    support_rows["nonarithmetic"] += 1

    first_a2_rows = []
    # First middle constant bucket, a=2 and even b.  The high extreme common
    # power has e=gcd(b-1,b+2) in {1,3}.
    for b in range(2, 22, 2):
        m = b // 2
        e = gcd(b - 1, b + 2)
        small_power = (b - 1) // e
        large_power = (b + 2) // e
        particular = sp.Rational(large_power, small_power) * H * K ** (3 // e)
        cancellation = weight_bracket(-2, H, b + 2, K**large_power, u) + weight_bracket(
            b - 1, K**small_power, 1, particular, u
        )
        require(sp.simplify(cancellation) == 0, ("first bucket a=2 cancellation", b, e, sp.factor(cancellation)))
        constant_core = weight_bracket(-2, H, 1, particular, u) + weight_bracket(
            b - 1, K**small_power, -b, H**m, u
        )
        if e == 1:
            divisor = sp.diff(H, u) * K + 2 * H * sp.diff(K, u)
            require_polynomial_quotient(constant_core, divisor, ("first bucket a=2,e=1", b))
            with_homogeneous = kappa * K + particular
            require(
                sp.simplify(
                    weight_bracket(-2, H, b + 2, K**large_power, u)
                    + weight_bracket(b - 1, K**small_power, 1, with_homogeneous, u)
                )
                == 0,
                ("first bucket a=2,e=1 homogeneous solution", b),
            )
            require_polynomial_quotient(
                weight_bracket(-2, H, 1, with_homogeneous, u)
                + weight_bracket(b - 1, K**small_power, -b, H**m, u),
                divisor,
                ("first bucket a=2,e=1 full factor", b),
            )
            factor_kind = "HprimeK+2HKprime"
        else:
            require(e == 3, ("unexpected gcd", b, e))
            divisor = 3 * sp.diff(H, u) * K + 2 * H * sp.diff(K, u)
            require_polynomial_quotient(constant_core, divisor, ("first bucket a=2,e=3", b))

            # A polynomial homogeneous solution exists only on the cube
            # branch K=T^3.  There every term factors through H'T+2HT'.
            cube_K = T**3
            cube_particular = sp.Rational(large_power, small_power) * H * cube_K
            cube_middle = kappa * T + cube_particular
            cube_cancellation = weight_bracket(-2, H, b + 2, cube_K**large_power, u) + weight_bracket(
                b - 1, cube_K**small_power, 1, cube_middle, u
            )
            require(sp.simplify(cube_cancellation) == 0, ("first bucket cube cancellation", b))
            cube_constant = weight_bracket(-2, H, 1, cube_middle, u) + weight_bracket(
                b - 1, cube_K**small_power, -b, H**m, u
            )
            cube_divisor = sp.diff(H, u) * T + 2 * H * sp.diff(T, u)
            require_polynomial_quotient(cube_constant, cube_divisor, ("first bucket cube factor", b))
            factor_kind = "3HprimeK+2HKprime_or_cube_root"
        first_a2_rows.append((b, e, small_power, large_power, factor_kind))

    first_b2_rows = []
    # First middle bucket, b=2 and even a.  Odd a is killed by endpoint
    # valuation in the theorem.  The cancellation ODE and unit factor are:
    for a in range(2, 22, 2):
        m = a // 2
        middle = kappa * K ** (a - 1) + 2 * a * H**m * K ** (2 * a - 1)
        cancellation = weight_bracket(-a, H**m, 2 * a, K ** (2 * a), u) + weight_bracket(
            1, K, a - 1, middle, u
        )
        require(sp.simplify(cancellation) == 0, ("first bucket b=2 cancellation", a, sp.factor(cancellation)))
        constant_core = weight_bracket(-a, H**m, a - 1, middle, u) + weight_bracket(1, K, -2, H, u)
        divisor = sp.diff(H, u) * K + 2 * H * sp.diff(K, u)
        require_polynomial_quotient(constant_core, divisor, ("first bucket b=2 factor", a))
        first_b2_rows.append((a, m, "HprimeK+2HKprime"))

    second_a2_rows = []
    # Second middle constant bucket, a=2.  This calculation is deliberately
    # separate: the ring has no symmetry exchanging h and u boundaries.
    for b in range(2, 13):
        far_weight = 2 * b + 1
        middle = kappa * H**b + sp.Rational(far_weight, 2) * H ** (2 * b - 1) * K ** (b - 1)
        cancellation = weight_bracket(-2, H**2, -b, middle, u) + weight_bracket(
            b - 1, K ** (b - 1), -far_weight, H**far_weight, u
        )
        require(sp.simplify(cancellation) == 0, ("second bucket a=2 cancellation", b, sp.factor(cancellation)))
        constant_core = weight_bracket(-2, H**2, 1, K, u) + weight_bracket(
            b - 1, K ** (b - 1), -b, middle, u
        )
        divisor = H * (sp.diff(H, u) * K + H * sp.diff(K, u))
        require_polynomial_quotient(constant_core, divisor, ("second bucket a=2 factor", b))
        second_a2_rows.append((b, far_weight, "H*(HprimeK+HKprime)"))

    second_b2_rows = []
    # Second middle bucket, b=2.  Here e=gcd(a,a+3) is 1 or 3.
    for a in range(2, 19):
        e = gcd(a, a + 3)
        small_power = a // e
        large_power = (a + 3) // e
        particular = sp.Rational(large_power, small_power) * H ** (3 // e) * K
        cancellation = weight_bracket(-a, H**small_power, -2, particular, u) + weight_bracket(
            1, K, -(a + 3), H**large_power, u
        )
        require(sp.simplify(cancellation) == 0, ("second bucket b=2 cancellation", a, e, sp.factor(cancellation)))
        constant_core = weight_bracket(-a, H**small_power, a - 1, K ** (a - 1), u) + weight_bracket(
            1, K, -2, particular, u
        )
        if e == 1:
            divisor = sp.diff(H, u) * K + H * sp.diff(K, u)
            require_polynomial_quotient(constant_core, divisor, ("second bucket b=2,e=1", a))
            with_homogeneous = kappa * H**2 + particular
            require(
                sp.simplify(
                    weight_bracket(-a, H**small_power, -2, with_homogeneous, u)
                    + weight_bracket(1, K, -(a + 3), H**large_power, u)
                )
                == 0,
                ("second bucket b=2,e=1 homogeneous solution", a),
            )
            require_polynomial_quotient(
                weight_bracket(-a, H**small_power, a - 1, K ** (a - 1), u)
                + weight_bracket(1, K, -2, with_homogeneous, u),
                divisor,
                ("second bucket b=2,e=1 full factor", a),
            )
            factor_kind = "HprimeK+HKprime"
        else:
            require(e == 3, ("unexpected second gcd", a, e))
            divisor = sp.diff(H, u) * K + 3 * H * sp.diff(K, u)
            require_polynomial_quotient(constant_core, divisor, ("second bucket b=2,e=3", a))

            # If the homogeneous solution exists, H=T^3 and it is kappa*T^2.
            cube_H = T**3
            cube_particular = sp.Rational(large_power, small_power) * T**3 * K
            cube_middle = kappa * T**2 + cube_particular
            cube_cancellation = weight_bracket(-a, cube_H**small_power, -2, cube_middle, u) + weight_bracket(
                1, K, -(a + 3), cube_H**large_power, u
            )
            require(sp.simplify(cube_cancellation) == 0, ("second bucket cube cancellation", a))
            cube_constant = weight_bracket(-a, cube_H**small_power, a - 1, K ** (a - 1), u) + weight_bracket(
                1, K, -2, cube_middle, u
            )
            cube_divisor = sp.diff(T, u) * K + T * sp.diff(K, u)
            require_polynomial_quotient(cube_constant, cube_divisor, ("second bucket cube factor", a))
            factor_kind = "HprimeK+3HKprime_or_cube_root"
        second_b2_rows.append((a, e, small_power, large_power, factor_kind))

    print("theorem=THM-3693-y0-collision-ring-two-by-three-weight-no-go")
    print(f"support_sum_window=-6..6;rows={support_rows}")
    print("nonarithmetic_mechanism=unique_duplicate_bucket_plus_Poisson_centralizer_transitivity")
    print(f"first_middle_a2_rows={tuple(first_a2_rows)}")
    print(f"first_middle_b2_even_rows={tuple(first_b2_rows)}")
    print(f"second_middle_a2_rows={tuple(second_a2_rows)}")
    print(f"second_middle_b2_rows={tuple(second_b2_rows)}")
    print("edge_mechanism=a_or_b_nonpositive_or_one_reduces_to_opposite_weight,weight_zero,or_endpoint_no_go")
    print("consequence=no_Darboux_pair_with_support_sizes_at_most_2_and_3_in_either_order")
    print("scope=all_degree_in_R;every_survivor_needs_at_least_three_active_weights_in_each_output")
    print("commands=python3 -B 04-computation/jacobian_y0_two_by_three_weight_no_go_thm3693.py;python3 -B -O 04-computation/jacobian_y0_two_by_three_weight_no_go_thm3693.py")


if __name__ == "__main__":
    main()
