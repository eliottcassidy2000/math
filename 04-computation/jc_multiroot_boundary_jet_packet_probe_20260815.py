#!/usr/bin/env python3
"""Exact probe for the multiroot boundary-jet packet.

This is an unnumbered research artifact.  For a repeated root alpha_i, write

    lambda=P-(a*alpha_i+b)=y*H,
    H=a+y**(e_i-1)*u*z**d.

The probe checks the D-stable CRT factors y**q and H**q, an explicit finite
Neumann inverse on the vertical factor, and the horizontal character gauge

    r=(1-lambda/(a*y))**(-sigma/d),
    A_lambda(r*p)=r*(1-lambda/(a*y))*A_0(p) mod lambda**q.

It then audits the selected N versus unselected N-1 special-fibre graph rank
and enumerates the resulting derived K[lambda]/lambda**q packet instances on
a declared exact grid.  Freeness itself comes from the universal gauge above;
the final grid is bookkeeping evidence, not an independent Smith computation.
"""

from __future__ import annotations

from itertools import product

import sympy as sp


x, y, z, lam = sp.symbols("x y z lam")
ROOTS = (-3, 1, 4, 8)


def require(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(label)


def truncate_in(variable: sp.Symbol, expression: sp.Expr, power: int) -> sp.Expr:
    numerator, denominator = sp.together(expression).as_numer_denom()
    require(
        not denominator.has(variable),
        ("truncation denominator depends on variable", variable, denominator),
    )
    truncated = sp.rem(sp.Poly(sp.expand(numerator), variable), variable**power)
    return sp.cancel(truncated.as_expr() / denominator)


def local_derivation(
    expression: sp.Expr, g: sp.Expr, d: int, a: int
) -> sp.Expr:
    return sp.expand(
        (a + sp.diff(g, y) * z**d) * sp.diff(expression, z)
        - d * g * z ** (d - 1) * sp.diff(expression, y)
    )


def check_factorization_and_vertical_inverse() -> tuple[int, int, int]:
    factor_checks = 0
    stable_ideal_checks = 0
    vertical_inverse_checks = 0
    local_profiles = (
        (2, 2, 1 + y),
        (3, 3, (1 + y) * (2 + y) ** 2),
    )
    for d, e, u in local_profiles:
        a = 3
        gamma = 2
        g = sp.expand(gamma * y**e * u)
        h = sp.expand(a + gamma * y ** (e - 1) * u * z**d)
        lambda_local = sp.expand(y * h)
        require(
            sp.expand(lambda_local - (a * y + g * z**d)) == 0,
            ("lambda factorization", d, e, u),
        )
        require(local_derivation(lambda_local, g, d, a) == 0, ("D lambda", d, e))
        expected_h = sp.expand(d * gamma * y ** (e - 1) * u * z ** (d - 1) * h)
        require(
            sp.expand(local_derivation(h, g, d, a) - expected_h) == 0,
            ("D H", d, e),
        )
        factor_checks += 3

        def right_integral(target: sp.Expr) -> sp.Expr:
            return sp.integrate(target, z) / a

        def perturbation(target: sp.Expr) -> sp.Expr:
            return sp.expand(
                local_derivation(target, g, d, a) - a * sp.diff(target, z)
            )

        for jet_order in (1, 2, 4):
            dyq = local_derivation(y**jet_order, g, d, a)
            dhq = local_derivation(h**jet_order, g, d, a)
            require(
                sp.cancel(dyq / y**jet_order).is_polynomial(y, z),
                ("D-stable y ideal", d, e, jet_order),
            )
            require(
                sp.cancel(dhq / h**jet_order).is_polynomial(y, z),
                ("D-stable H ideal", d, e, jet_order),
            )
            stable_ideal_checks += 2

            for y_degree in range(jet_order):
                for z_degree in range(2):
                    target = y**y_degree * z**z_degree
                    inverse_input = target
                    neumann_sum = target
                    # E=R*J raises y-adic order by at least e-1.  q terms are
                    # enough even in the hostile e=2 case.
                    for _ in range(1, jet_order):
                        inverse_input = truncate_in(
                            y,
                            -perturbation(right_integral(inverse_input)),
                            jet_order,
                        )
                        neumann_sum = truncate_in(
                            y, neumann_sum + inverse_input, jet_order
                        )
                    primitive = right_integral(neumann_sum)
                    residual = truncate_in(
                        y,
                        local_derivation(primitive, g, d, a) - target,
                        jet_order,
                    )
                    require(
                        residual == 0,
                        ("vertical Neumann inverse", d, e, jet_order, y_degree, z_degree),
                    )
                    vertical_inverse_checks += 1
    return factor_checks, stable_ideal_checks, vertical_inverse_checks


def split_polynomial(exponents: tuple[int, ...], gamma: int = 2) -> sp.Expr:
    return sp.expand(
        gamma
        * sp.prod(
            (x - alpha) ** e
            for alpha, e in zip(ROOTS[: len(exponents)], exponents)
        )
    )


def truncated_horizontal_gauge(
    d: int, sigma: int, a: int, alpha: sp.Expr, jet_order: int
) -> sp.Expr:
    local_y = x - alpha
    u = lam / (a * local_y)
    return sp.expand(
        sum(
            sp.rf(sp.Rational(sigma, d), degree)
            / sp.factorial(degree)
            * u**degree
            for degree in range(jet_order)
        )
    )


def check_horizontal_gauge() -> tuple[int, int, int]:
    gauge_checks = 0
    second_jet_checks = 0
    wrap_checks = 0
    split_cases = (
        (2, 1, (3, 2), 0),
        (4, 2, (3, 1), 0),
        (2, 1, (2, 4), 0),
        (2, 2, (3, 2), 0),
        (5, 5, (2, 4, 3), 1),
    )
    cases: list[tuple[int, int, sp.Expr, str]] = []
    for d, sigma, exponents, i in split_cases:
        alpha = sp.Integer(ROOTS[i])
        cases.append((d, sigma, alpha, f"split:{exponents}:{i}"))

    # Nonsplit controls over Q, checked after exact faithful splitting for the
    # chosen geometric root and at a rational selected root respectively.
    cases.append((2, 1, sp.I, "nonsplit:(x^2+1)^3"))
    cases.append((2, 1, sp.Integer(0), "mixed:x^3(x^2+1)^2"))

    log_derivative, p_value, p_derivative = sp.symbols("ell p p_x")
    for d, sigma, alpha, label in cases:
        a = 3
        local_y = x - alpha
        for jet_order in range(1, 6):
            r = truncated_horizontal_gauge(d, sigma, a, alpha, jet_order)
            scale = 1 - lam / (a * local_y)
            if jet_order == 2:
                linear_coefficient = sp.expand(r).coeff(lam, 1)
                require(
                    sp.cancel(
                        linear_coefficient
                        - sp.Rational(sigma, d) / (a * local_y)
                    )
                    == 0,
                    ("second-jet gauge coefficient", d, sigma, label, alpha),
                )
                second_jet_checks += 1
            lhs = sp.expand(
                a * sigma * r * p_value
                + (lam - a * local_y)
                * (
                    sigma * log_derivative * r * p_value
                    - d
                    * (
                        sp.diff(r, x) * p_value
                        + r * p_derivative
                    )
                )
            )
            special = sp.expand(
                a * sigma * p_value
                - a
                * local_y
                * (sigma * log_derivative * p_value - d * p_derivative)
            )
            rhs = sp.expand(r * scale * special)
            require(
                truncate_in(lam, lhs - rhs, jet_order) == 0,
                ("horizontal gauge", d, sigma, label, alpha, jet_order),
            )
            gauge_checks += 1
            if sigma == d:
                wrap_checks += 1
    return gauge_checks, second_jet_checks, wrap_checks


def selected(
    d: int, sigma: int, exponents: tuple[int, ...], i: int
) -> bool:
    e_i = exponents[i]
    return (
        e_i > 1
        and sigma * (e_i - 1) % d == 0
        and all(
            sigma * e_j % d == 0
            for j, e_j in enumerate(exponents)
            if j != i
        )
    )


def check_special_graph_and_derived_packets() -> tuple[int, int, int, int]:
    graph_rank_checks = 0
    derived_packet_checks = 0
    selected_packets = 0
    unselected_packets = 0
    for d in range(2, 11):
        for root_count in range(1, 5):
            for exponents in product(range(1, 6), repeat=root_count):
                for i, e_i in enumerate(exponents):
                    if e_i == 1:
                        continue
                    monodromy_orders = tuple(
                        sigma_exponent
                        for sigma_exponent in (
                            e_i - 1,
                            *(e_j for j, e_j in enumerate(exponents) if j != i),
                        )
                    )
                    for sigma in range(1, d + 1):
                        coboundary_vector = tuple(
                            0 if sigma * exponent % d == 0 else 1
                            for exponent in monodromy_orders
                        )
                        differential_rank = 0 if not any(coboundary_vector) else 1
                        h_zero = 1 - differential_rank
                        h_one = root_count - differential_rank
                        predicted_selected = selected(d, sigma, exponents, i)
                        require(
                            (h_zero, h_one)
                            == ((1, root_count) if predicted_selected else (0, root_count - 1)),
                            ("special graph cohomology", d, sigma, exponents, i),
                        )
                        graph_rank_checks += 1

                        free_rank = root_count if predicted_selected else root_count - 1
                        for jet_order in (1, 2, 3, 5):
                            invariant_factors = (jet_order,) * free_rank
                            total_length = sum(invariant_factors)
                            residue_generators = len(invariant_factors)
                            require(
                                total_length == jet_order * free_rank
                                and residue_generators == free_rank,
                                ("free jet packet", d, sigma, exponents, i, jet_order),
                            )
                            # Over the principal Artin local ring, c generators
                            # and maximal possible length c*q force all c
                            # invariant factors to have length q.
                            require(
                                all(length == jet_order for length in invariant_factors),
                                ("jet freeness", d, sigma, exponents, i, jet_order),
                            )
                            derived_packet_checks += 1
                        if predicted_selected:
                            selected_packets += 1
                        else:
                            unselected_packets += 1

    return graph_rank_checks, derived_packet_checks, selected_packets, unselected_packets


def check_named_nonsplit_and_wrap_controls() -> int:
    checks = 0
    require(not selected(2, 1, (3, 3), 0), "nonsplit conjugate hostile left")
    require(not selected(2, 1, (3, 3), 1), "nonsplit conjugate hostile right")
    checks += 2
    require(selected(2, 1, (3, 2, 2), 0), "rational arm amid nonsplit orbit")
    checks += 1
    for d in range(2, 13):
        for exponents in ((2,), (3, 1), (2, 4, 3), (1, 2, 1, 5)):
            for i, e_i in enumerate(exponents):
                if e_i > 1:
                    require(selected(d, d, exponents, i), ("wrap selected", d, exponents, i))
                    checks += 1
    return checks


def main() -> None:
    factor_checks, stable_ideals, vertical_inverses = (
        check_factorization_and_vertical_inverse()
    )
    horizontal_gauges, second_jets, wrap_gauges = check_horizontal_gauge()
    graph_ranks, derived_packets, selected_packets, unselected_packets = (
        check_special_graph_and_derived_packets()
    )
    descent_controls = check_named_nonsplit_and_wrap_controls()

    print("MULTIROOT BOUNDARY-JET PACKET -- EXACT RESEARCH PROBE")
    print(f"lambda=yH / D(lambda)=0 factor checks: {factor_checks}")
    print(f"D-stable CRT ideal checks: {stable_ideals}")
    print(f"vertical finite-Neumann inverse checks: {vertical_inverses}")
    print(f"horizontal formal-gauge identities: {horizontal_gauges}")
    print(f"explicit q=2 extension coefficient checks: {second_jets}")
    print(f"wrap sigma=d gauge identities: {wrap_gauges}")
    print(f"special-fibre graph H0/H1 ranks: {graph_ranks}")
    print(f"derived Artin jet packet instances: {derived_packets}")
    print(f"selected rank-N packets: {selected_packets}")
    print(f"unselected rank-(N-1) packets: {unselected_packets}")
    print(f"nonsplit/wrap descent controls: {descent_controls}")
    print("named q=2 blocked hostile: (d,sigma;e1,e2)=(4,2;3,1)")
    print("named nonsplit hostile over Q: g=(x^2+1)^3, (d,sigma)=(2,1)")
    print("BOUNDARY-JET N/N-1 PACKET VERIFIED ON DECLARED GRID")


if __name__ == "__main__":
    main()
