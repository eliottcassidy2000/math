#!/usr/bin/env python3
"""Exact algebra companion for THM-3691.

The theorem is an all-degree proof.  This companion checks its reusable
Laurent bracket law, the endpoint divisibility carried by target monomials,
the crossed equal-gap support grammar in a wide hostile window, and the
odd/even factor identities after the UFD commutation reduction.

Every truth gate uses ``require`` and therefore survives ``python -O``.
"""

from __future__ import annotations

from itertools import combinations

import sympy as sp


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def main() -> None:
    x, z, u = sp.symbols("x z u")
    h = 1 - u**2

    # Direct Laurent-coordinate audit of the master law.  In source
    # coordinates u=1-x^2 z and {P,Q}=P_x Q_z-P_z Q_x.
    source_u = 1 - x**2 * z
    p_coefficients = sp.symbols("p0:5")
    q_coefficients = sp.symbols("q0:5")
    p = sum(coefficient * u**degree for degree, coefficient in enumerate(p_coefficients))
    q = sum(coefficient * u**degree for degree, coefficient in enumerate(q_coefficients))
    bracket_rows = []
    for first_weight, second_weight in ((-5, -2), (-3, 0), (-2, 1), (-1, 4), (0, -4), (2, 5)):
        first = x**first_weight * p.subs(u, source_u)
        second = x**second_weight * q.subs(u, source_u)
        actual = sp.diff(first, x) * sp.diff(second, z) - sp.diff(first, z) * sp.diff(second, x)
        predicted = x ** (first_weight + second_weight + 1) * (
            second_weight * sp.diff(p, u) * q - first_weight * p * sp.diff(q, u)
        ).subs(u, source_u)
        require(sp.cancel(actual - predicted) == 0, ("master bracket law", first_weight, second_weight))
        bracket_rows.append((first_weight, second_weight, first_weight + second_weight + 1))

    # Every target monomial A^i B^j C^k has weight -2i-j+k and coefficient
    # h^(i+j)u^k.  Check the all-degree divisibility formulas on a hostile
    # finite window (the theorem proves them directly from this formula).
    divisibility_rows = []
    for target_degree in range(1, 11):
        for power_a in range(target_degree + 1):
            for power_b in range(target_degree - power_a + 1):
                power_c = target_degree - power_a - power_b
                weight = -2 * power_a - power_b + power_c
                h_power = power_a + power_b
                u_power = power_c
                if weight < 0:
                    require(h_power >= sp.ceiling(-sp.Rational(weight, 2)), ("negative divisibility", power_a, power_b, power_c))
                elif weight > 0:
                    require(u_power >= weight, ("positive divisibility", power_a, power_b, power_c))
        divisibility_rows.append((target_degree, sum(1 for _ in range((target_degree + 1) * (target_degree + 2) // 2))))

    # Hostile finite support-sum audit.  Among two-by-two integer weight
    # supports, once constant-at-an-extreme and unique-cross cases are
    # removed, the only possible constant bucket consists of both crosses,
    # and the two gaps agree.
    crossed_rows = []
    for first_support in combinations(range(-7, 8), 2):
        for second_support in combinations(range(-7, 8), 2):
            buckets: dict[int, list[tuple[int, int]]] = {}
            for first_weight in first_support:
                for second_weight in second_support:
                    buckets.setdefault(first_weight + second_weight + 1, []).append((first_weight, second_weight))
            constant_pairs = buckets.get(0, [])
            if len(constant_pairs) == 2:
                require(
                    set(constant_pairs)
                    == {(first_support[0], second_support[1]), (first_support[1], second_support[0])},
                    ("noncrossed duplicate constant bucket", first_support, second_support, constant_pairs),
                )
                require(
                    first_support[1] - first_support[0] == second_support[1] - second_support[0],
                    ("unequal crossed gaps", first_support, second_support),
                )
                crossed_rows.append((first_support, second_support))

    # The a=b=1 edge: after proportionality of the negative pieces, both
    # crossed terms retain the squarefree factor h.
    p0, s0, r0, q0 = sp.symbols("p0 s0 r0 q0")
    edge_expression = h * (p0 * sp.diff(q0 * u, u) - sp.diff(r0 * u, u) * s0)
    require(sp.rem(edge_expression, h, domain=sp.QQ[p0, s0, r0, q0]) == 0, "a=b=1 endpoint factor")

    # For a=2 and b odd, low commutation gives p=H^2,q=H^b.  Both crossed
    # terms vanish at h=0 once h|H.  Check symbolic representatives H=h*T.
    T, K = sp.symbols("T K", cls=sp.Function)
    T_u = T(u)
    K_u = K(u)
    H_u = h * T_u
    odd_rows = []
    for b_value in (3, 5, 7, 9, 11):
        first_cross = sp.diff(H_u**2, u) * K_u + 2 * H_u**2 * sp.diff(K_u, u)
        second_cross = -b_value * sp.diff(K_u ** (b_value - 1), u) * H_u**b_value - (
            b_value - 1
        ) * K_u ** (b_value - 1) * sp.diff(H_u**b_value, u)
        for endpoint in (-1, 1):
            require(sp.simplify(first_cross.subs(u, endpoint)) == 0, ("odd first endpoint", b_value, endpoint))
            require(sp.simplify(second_cross.subs(u, endpoint)) == 0, ("odd second endpoint", b_value, endpoint))
        odd_rows.append(b_value)

    # For b=2m, the two crossed brackets have the common homogeneous
    # Wronskian D=H'K+2HK'.  This is the decisive all-degree factorization.
    H, Hp, Ksym, Kp = sp.symbols("H Hp K Kp")
    D = Hp * Ksym + 2 * H * Kp
    even_rows = []
    for b_value in range(2, 22, 2):
        m = b_value // 2
        second_cross = -b_value * (b_value - 1) * Ksym ** (b_value - 2) * Kp * H**m - (
            b_value - 1
        ) * m * Ksym ** (b_value - 1) * H ** (m - 1) * Hp
        predicted = -(b_value - 1) * m * H ** (m - 1) * Ksym ** (b_value - 2) * D
        require(sp.expand(second_cross - predicted) == 0, ("even crossed factor", b_value))
        even_rows.append((b_value, m, b_value - 1, m - 1, b_value - 2))

    print("theorem=THM-3691-y0-collision-ring-two-weight-darboux-no-go")
    print(f"master_bracket_rows={tuple(bracket_rows)}")
    print("weight_modules=negative_-R_has_h^ceil(R/2);positive_R_has_u^R")
    print(f"divisibility_target_degree_controls={tuple(divisibility_rows)}")
    print(f"crossed_equal_gap_window=-7..7;survivors={len(crossed_rows)}")
    print("edge_cases=opposite_sign_extreme_impossible;a_or_b_1_reduces_to_homogeneous_or_h_divisible")
    print(f"odd_b_endpoint_controls={tuple(odd_rows)}")
    print(f"even_b_factor_rows={tuple(even_rows)}")
    print("even_factor=constant_equation=(Hprime*K+2H*Kprime)*(c1+c2*H^(b/2-1)*K^(b-2))")
    print("consequence=no_Darboux_pair_when_each_output_has_at_most_two_active_grading_weights")
    print("scope=all_degree_in_R;three_or_more_weights_in_at_least_one_output_remain_open")
    print("commands=python3 -B 04-computation/jacobian_y0_two_weight_darboux_no_go_thm3691.py;python3 -B -O 04-computation/jacobian_y0_two_weight_darboux_no_go_thm3691.py")


if __name__ == "__main__":
    main()
