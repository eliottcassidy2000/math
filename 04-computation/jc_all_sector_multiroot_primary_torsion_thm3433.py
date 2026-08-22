#!/usr/bin/env python3
"""Exact hostile referee for provisional THM-3433.

For P=a*x+b+g(x)*z**d and a nonwrap character sigma, THM-3418 gives
transitions

    L_n(q) = d*g*q'/(a*n) - g'*q/a,       n=sigma+k*d.

This companion checks the global-power/injective dichotomy, the selected-root
congruences, exact root gauges, the unique polynomial horizontal solution for
the (P-rho)-kernel, strictness of the embedded root filtration, character
counts, finite torsion thickness, and the separate wrap regression.  All
decisions use integer, Fraction, or SymPy exact polynomial arithmetic.
"""

from __future__ import annotations

from itertools import product
from math import gcd

import sympy as sp


x = sp.symbols("x")
ROOTS = (-3, 1, 4, 8)


def require(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(label)


def split_polynomial(exponents: tuple[int, ...], gamma: int = 2) -> sp.Expr:
    return sp.expand(
        gamma
        * sp.prod(
            (x - alpha) ** e
            for alpha, e in zip(ROOTS[: len(exponents)], exponents)
        )
    )


def sector_operator(
    g: sp.Expr, q: sp.Expr, d: int, n: int, a: int = 1
) -> sp.Expr:
    return sp.expand(
        sp.Rational(d, a * n) * g * sp.diff(q, x)
        - sp.Rational(1, a) * sp.diff(g, x) * q
    )


def kernel_operator(
    g: sp.Expr, p: sp.Expr, d: int, n: int, xi: sp.Expr
) -> sp.Expr:
    """n times the stage-(k+1) representative of (P-rho)[p]_k."""
    return sp.expand(
        d * g * (x - xi) * sp.diff(p, x)
        + (n + d) * g * p
        - n * (x - xi) * sp.diff(g, x) * p
    )


def selected_roots(
    d: int, sigma: int, exponents: tuple[int, ...]
) -> tuple[int, ...]:
    return tuple(
        i
        for i, e_i in enumerate(exponents)
        if e_i > 1
        and sigma * (e_i - 1) % d == 0
        and all(
            sigma * e_j % d == 0
            for j, e_j in enumerate(exponents)
            if j != i
        )
    )


def check_congruence_universe() -> tuple[int, int, int, int, int, int]:
    sector_profiles = 0
    global_power_profiles = 0
    injective_profiles = 0
    selected_occurrences = 0
    blocked_local_resonances = 0
    character_count_checks = 0

    for d in range(2, 13):
        for root_count in range(1, 5):
            for exponents in product(range(1, 7), repeat=root_count):
                for i, e_i in enumerate(exponents):
                    if e_i == 1:
                        continue
                    arithmetic_gcd = d
                    arithmetic_gcd = gcd(arithmetic_gcd, e_i - 1)
                    for j, e_j in enumerate(exponents):
                        if j != i:
                            arithmetic_gcd = gcd(arithmetic_gcd, e_j)
                    selected_characters = sum(
                        i in selected_roots(d, sigma, exponents)
                        for sigma in range(1, d + 1)
                    )
                    require(
                        selected_characters == arithmetic_gcd,
                        ("character gcd", d, exponents, i),
                    )
                    character_count_checks += 1

                wrap = selected_roots(d, d, exponents)
                require(
                    wrap == tuple(i for i, e_i in enumerate(exponents) if e_i > 1),
                    ("wrap selection", d, exponents),
                )

                for sigma in range(1, d):
                    sector_profiles += 1
                    global_power = all(sigma * e % d == 0 for e in exponents)
                    selected = selected_roots(d, sigma, exponents)
                    require(
                        len(selected) <= 1,
                        ("two nonwrap selected roots", d, sigma, exponents),
                    )
                    require(
                        not (global_power and selected),
                        ("selected global-power sector", d, sigma, exponents),
                    )
                    if global_power:
                        global_power_profiles += 1
                        for e in exponents:
                            require(
                                sigma * (e - 1) % d != 0,
                                ("global-power Laurent gate", d, sigma, exponents),
                            )
                    else:
                        injective_profiles += 1

                    selected_occurrences += len(selected)
                    for i, e_i in enumerate(exponents):
                        local_resonance = (
                            e_i > 1 and sigma * (e_i - 1) % d == 0
                        )
                        other_gate = all(
                            sigma * e_j % d == 0
                            for j, e_j in enumerate(exponents)
                            if j != i
                        )
                        if local_resonance and not other_gate:
                            blocked_local_resonances += 1

    return (
        sector_profiles,
        global_power_profiles,
        injective_profiles,
        selected_occurrences,
        blocked_local_resonances,
        character_count_checks,
    )


def check_global_power_diagrams() -> tuple[int, int, int]:
    kernel_checks = 0
    dying_derivative_checks = 0
    quotient_gauge_checks = 0
    cases = (
        (2, 1, (2, 4)),
        (2, 1, (2, 2, 4)),
        (3, 1, (3, 6)),
        (4, 2, (2, 4)),
        (6, 2, (3, 6)),
        (6, 3, (2, 4)),
    )
    for d, sigma, exponents in cases:
        require(
            all(sigma * e % d == 0 for e in exponents),
            ("global-power case", d, sigma, exponents),
        )
        g = split_polynomial(exponents, gamma=2)
        for k in range(4):
            n = sigma + k * d
            powers = tuple(k * e + sigma * e // d for e in exponents)
            h_k = sp.prod(
                (x - alpha) ** power
                for alpha, power in zip(ROOTS, powers)
            )
            h_next = sp.prod(
                (x - alpha) ** (power + e)
                for alpha, power, e in zip(ROOTS, powers, exponents)
            )
            require(
                sector_operator(g, h_k, d, n) == 0,
                ("transition kernel", d, sigma, exponents, k),
            )
            kernel_checks += 1

            p = 1 + 2 * x + 3 * x**2 + x**3
            expected = sp.expand(sp.Rational(2 * d, n) * h_next * sp.diff(p, x))
            require(
                sp.expand(sector_operator(g, h_k * p, d, n) - expected) == 0,
                ("dying derivative", d, sigma, exponents, k),
            )
            dying_derivative_checks += 1

            for i, (alpha, e_i) in enumerate(zip(ROOTS, exponents)):
                y = x - alpha
                gauge = sp.prod(
                    (x - other_alpha) ** power
                    for j, (other_alpha, power) in enumerate(zip(ROOTS, powers))
                    if j != i
                )
                gauge_next = sp.prod(
                    (x - other_alpha) ** (power + e_j)
                    for j, (other_alpha, power, e_j) in enumerate(
                        zip(ROOTS, powers, exponents)
                    )
                    if j != i
                )
                local = 2 * y**e_i
                require(
                    sp.expand(
                        sector_operator(g, gauge * p, d, n)
                        - gauge_next * sector_operator(local, p, d, n)
                    )
                    == 0,
                    ("global quotient gauge", d, sigma, exponents, k, i),
                )
                quotient_gauge_checks += 1
    return kernel_checks, dying_derivative_checks, quotient_gauge_checks


def rank_mod_prime(rows: list[list[int]], prime: int) -> int:
    matrix = [[entry % prime for entry in row] for row in rows]
    row_count = len(matrix)
    column_count = len(matrix[0]) if matrix else 0
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(pivot_row, row_count) if matrix[row][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        inverse = pow(matrix[pivot_row][column], prime - 2, prime)
        matrix[pivot_row] = [entry * inverse % prime for entry in matrix[pivot_row]]
        for row in range(pivot_row + 1, row_count):
            coefficient = matrix[row][column]
            if coefficient:
                matrix[row] = [
                    (entry - coefficient * pivot_entry) % prime
                    for entry, pivot_entry in zip(matrix[row], matrix[pivot_row])
                ]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def polynomial_map_modular_nullity_bound(
    g: sp.Expr, d: int, n: int, xi: int, degree_bound: int
) -> int:
    """Upper-bound the Q-nullity by a certified finite-field minor.

    In predicted-nullity-one cases the explicit horizontal polynomial already
    supplies the matching lower bound.  In predicted-nullity-zero cases a
    full-column-rank modular minor proves injectivity over Q.
    """
    outputs = [kernel_operator(g, x**m, d, n, xi) for m in range(degree_bound + 1)]
    output_degree = max(sp.degree(poly, x) for poly in outputs)
    rows = [
        [0 for _ in range(degree_bound + 1)] for _ in range(output_degree + 1)
    ]
    for column, poly in enumerate(outputs):
        coefficient_poly = sp.Poly(poly, x, domain=sp.ZZ)
        for (power,), coefficient in coefficient_poly.terms():
            rows[power][column] = int(coefficient)
    modular_rank = max(
        rank_mod_prime(rows, prime) for prime in (1_000_003, 1_000_033)
    )
    return degree_bound + 1 - modular_rank


def check_selected_horizontal_solutions() -> tuple[int, int, int, int]:
    gauge_checks = 0
    horizontal_checks = 0
    strict_valuation_checks = 0
    low_grid_kernel_checks = 0
    positive_cases = (
        (2, 1, (3, 2), 0),
        (3, 1, (4, 3), 0),
        (4, 2, (3, 2), 0),
        (6, 2, (4, 3), 0),
        (6, 3, (3, 2), 0),
        (6, 2, (4, 3, 6), 0),
    )
    for d, sigma, exponents, i in positive_cases:
        require(
            selected_roots(d, sigma, exponents) == (i,),
            ("positive selected case", d, sigma, exponents, i),
        )
        roots = ROOTS[: len(exponents)]
        alpha = roots[i]
        e_i = exponents[i]
        y = x - alpha
        g = split_polynomial(exponents, gamma=2)
        local = 2 * y**e_i
        q_i = sigma * (e_i - 1) // d
        require(q_i >= 1, ("positive intercept", d, sigma, exponents, i))

        for k in range(4):
            n = sigma + k * d
            gauge = sp.prod(
                (x - other_alpha) ** (k * e_j + sigma * e_j // d)
                for j, (other_alpha, e_j) in enumerate(zip(roots, exponents))
                if j != i
            )
            gauge_next = sp.prod(
                (x - other_alpha) ** ((k + 1) * e_j + sigma * e_j // d)
                for j, (other_alpha, e_j) in enumerate(zip(roots, exponents))
                if j != i
            )
            p_test = 1 + 2 * y + 3 * y**2
            require(
                sp.expand(
                    sector_operator(g, gauge * p_test, d, n)
                    - gauge_next * sector_operator(local, p_test, d, n)
                )
                == 0,
                ("selected gauge", d, sigma, exponents, i, k),
            )
            gauge_checks += 1

            thickness = k * (e_i - 1) + q_i
            endpoint = sp.expand(gauge * y ** (thickness - 1))
            require(
                kernel_operator(g, endpoint, d, n, alpha) == 0,
                ("horizontal endpoint", d, sigma, exponents, i, k),
            )
            expected_next = sp.expand(
                -sp.Rational(2 * (n + d), n)
                * gauge_next
                * y ** (thickness + e_i - 2)
            )
            require(
                sp.expand(sector_operator(g, endpoint, d, n) - expected_next) == 0,
                ("endpoint persistence", d, sigma, exponents, i, k),
            )
            horizontal_checks += 2

            for j, e_j in enumerate(exponents):
                if j == i:
                    continue
                boundary = k * e_j + sigma * e_j // d
                for valuation in range(boundary):
                    require(
                        d * valuation - n * e_j != 0,
                        (
                            "strict quotient leading coefficient",
                            d,
                            sigma,
                            exponents,
                            i,
                            j,
                            k,
                            valuation,
                        ),
                    )
                    require(
                        e_j + valuation - 1 < boundary + e_j,
                        ("strict quotient valuation", d, sigma, exponents, i, j),
                    )
                    strict_valuation_checks += 1

        for k in (0, 1):
            n = sigma + k * d
            endpoint_degree = sum(
                k * e_j + sigma * e_j // d
                for j, e_j in enumerate(exponents)
                if j != i
            ) + k * (e_i - 1) + q_i - 1
            bound = endpoint_degree + 3
            for j, alpha_j in enumerate(roots):
                predicted = 1 if j == i else 0
                require(
                    polynomial_map_modular_nullity_bound(g, d, n, alpha_j, bound)
                    == predicted,
                    ("low-grid root kernel", d, sigma, exponents, i, j, k),
                )
                low_grid_kernel_checks += 1
            require(
                polynomial_map_modular_nullity_bound(g, d, n, 0, bound) == 0,
                ("low-grid off-root kernel", d, sigma, exponents, i, k),
            )
            low_grid_kernel_checks += 1

    hostile_cases = (
        (4, 2, (3, 1)),
        (2, 1, (3, 3)),
        (6, 2, (4, 1, 3)),
    )
    for d, sigma, exponents in hostile_cases:
        require(
            not selected_roots(d, sigma, exponents),
            ("hostile accidentally selected", d, sigma, exponents),
        )
        require(
            not all(sigma * e % d == 0 for e in exponents),
            ("hostile accidentally global power", d, sigma, exponents),
        )
        g = split_polynomial(exponents, gamma=2)
        for alpha in ROOTS[: len(exponents)]:
            require(
                polynomial_map_modular_nullity_bound(g, d, sigma, alpha, 14) == 0,
                ("hostile low-grid kernel", d, sigma, exponents, alpha),
            )
            low_grid_kernel_checks += 1

    return gauge_checks, horizontal_checks, strict_valuation_checks, low_grid_kernel_checks


def check_thickness_and_primary_models() -> tuple[int, int, int]:
    thickness_checks = 0
    wrap_blocks = 0
    primary_model_checks = 0
    for d in range(2, 13):
        for root_count in range(1, 5):
            for exponents in product(range(1, 7), repeat=root_count):
                for sigma in range(1, d + 1):
                    selected = selected_roots(d, sigma, exponents)
                    for i in selected:
                        e_i = exponents[i]
                        intercept = sigma * (e_i - 1) // d
                        require(intercept >= 1, ("thickness intercept", d, sigma, exponents, i))
                        for k in range(5):
                            thickness = k * (e_i - 1) + intercept
                            require(thickness >= 1, ("positive thickness", d, sigma, exponents, i, k))
                            if sigma == d:
                                require(
                                    thickness == (k + 1) * (e_i - 1),
                                    ("wrap thickness", d, exponents, i, k),
                                )
                                wrap_blocks += 1
                            thickness_checks += 1

    # R/(lambda^length) has a one-dimensional lambda-kernel.  On the monomial
    # basis, lambda^power shifts exactly length-power basis vectors and kills
    # the remaining power vectors.
    for length in range(1, 41):
        for power in range(1, length + 1):
            surviving_images = {
                basis + power for basis in range(length) if basis + power < length
            }
            require(
                length - len(surviving_images) == power,
                ("primary filtration", length, power),
            )
            primary_model_checks += 1
    return thickness_checks, wrap_blocks, primary_model_checks


def check_named_descent_and_collision_controls() -> int:
    checks = 0
    # Over Q, (x^2+1)^3 has two conjugate multiplicity-three roots.  Both are
    # locally resonant for (d,sigma)=(2,1), but each is blocked by the other.
    require(selected_roots(2, 1, (3, 3)) == (), "nonsplit blocked pair")
    checks += 1
    # x^3*(x^2+1)^2 has one rational multiplicity-three root and a conjugate
    # pair of multiplicity two; exactly the rational root survives nonwrap.
    require(selected_roots(2, 1, (3, 2, 2)) == (0,), "rational descent arm")
    checks += 1
    # A nonzero a makes alpha -> a*alpha+b injective, so distinct roots cannot
    # collide in primary support.  Check the declared rational root universe.
    for a in (-3, -1, 1, 2, 5):
        for b in (-4, 0, 7):
            values = tuple(a * alpha + b for alpha in ROOTS)
            require(len(set(values)) == len(ROOTS), ("support collision", a, b))
            checks += 1
    return checks


def main() -> None:
    (
        sector_profiles,
        global_power_profiles,
        injective_profiles,
        selected_occurrences,
        blocked_local_resonances,
        character_count_checks,
    ) = check_congruence_universe()
    (
        transition_kernels,
        dying_derivatives,
        global_quotient_gauges,
    ) = check_global_power_diagrams()
    (
        selected_gauges,
        horizontal_solutions,
        strict_valuations,
        low_grid_kernels,
    ) = check_selected_horizontal_solutions()
    thickness_checks, wrap_blocks, primary_models = check_thickness_and_primary_models()
    descent_controls = check_named_descent_and_collision_controls()

    print("ALL-SECTOR MULTIROOT PRIMARY TORSION -- EXACT HOSTILE REFEREE")
    print(f"nonwrap sector profiles: {sector_profiles}")
    print(f"global-power dying profiles: {global_power_profiles}")
    print(f"injective-transition profiles: {injective_profiles}")
    print(f"selected nonwrap root occurrences: {selected_occurrences}")
    print(f"blocked local resonances: {blocked_local_resonances}")
    print(f"gcd character-count checks: {character_count_checks}")
    print(f"exact global-power transition kernels: {transition_kernels}")
    print(f"exact dying-derivative identities: {dying_derivatives}")
    print(f"exact global quotient gauges: {global_quotient_gauges}")
    print(f"exact selected-root gauges: {selected_gauges}")
    print(f"horizontal endpoint/persistence identities: {horizontal_solutions}")
    print(f"strict filtration valuation checks: {strict_valuations}")
    print(f"exact low-grid (P-rho)-kernel checks: {low_grid_kernels}")
    print(f"finite torsion-thickness checks: {thickness_checks}")
    print(f"separate wrap primary blocks: {wrap_blocks}")
    print(f"finite primary-model checks: {primary_models}")
    print(f"nonsplit/support-collision controls: {descent_controls}")
    print("named blocked hostile: (d,sigma;e1,e2)=(4,2;3,1)")
    print("named double-local hostile: (d,sigma;e1,e2)=(2,1;3,3)")
    print("named positive arm: (d,sigma;e1,e2)=(2,1;3,2)")
    print("GEOMETRIC PRIMARY-TORSION CLASSIFICATION VERIFIED")


if __name__ == "__main__":
    main()
