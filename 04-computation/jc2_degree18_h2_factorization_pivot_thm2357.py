#!/usr/bin/env python3
"""Exact factorization-pivot companion for THM-2357.

This uses the same moving total-ramification gauge as
``jc2_degree18_h2_moving_root_scratch.py`` but avoids the degree-20--37
terminal subresultants.  Write

    R_10 = H_2 S_4^2,
    S_4 = y^4+s3*y^3+s2*y^2+s1*y+s0,
    H_2 = L*y^2+h1*y+h0,

where L=lc(R_10).  Coefficients y^9,y^8,y^7,y^6 solve successively for
h1,h0,s1,s0.  Only six residual equations in B,C,s3,s2 remain.

"""

from __future__ import annotations

import argparse
import time

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def primitive(expression: sp.Expr, variables: tuple[sp.Symbol, ...]) -> sp.Expr:
    numerator, denominator = sp.together(expression).as_numer_denom()
    require(
        denominator.is_Integer and denominator != 0,
        "residual acquired a nonconstant denominator",
    )
    content, result = sp.Poly(numerator, *variables).primitive()
    require(content != 0, "primitive residual met zero")
    return result.as_expr()


def solve_linear(expression: sp.Expr, variable: sp.Symbol) -> sp.Expr:
    expanded = sp.expand(expression)
    coefficient = expanded.coeff(variable)
    remainder = sp.expand(expanded.subs(variable, 0))
    require(
        coefficient != 0
        and not coefficient.has(variable)
        and sp.expand(expanded - coefficient * variable - remainder) == 0,
        f"equation stopped being linear in {variable}",
    )
    return sp.cancel(-remainder / coefficient)


def coefficient(expression: sp.Expr, y: sp.Symbol, degree: int) -> sp.Expr:
    return sp.Poly(expression, y).coeff_monomial(y**degree)


def build() -> tuple[
    sp.Symbol,
    sp.Symbol,
    sp.Symbol,
    sp.Symbol,
    sp.Expr,
    sp.Expr,
    list[sp.Expr],
]:
    y = sp.symbols("y")
    bvar, cvar, dvar, wvar = sp.symbols("B C D W")
    s3, s2, s1, s0 = sp.symbols("s3 s2 s1 s0")
    h1, h0 = sp.symbols("h1 h0")

    quartic_p = (
        245 * y**4
        + 1890 * bvar * y**2
        - 24300 * bvar**2
        + 122472 * dvar
    )
    sextic_q = (
        539 * y**6
        + 11340 * bvar * y**4
        + 183708 * cvar * y**3
        + (72900 * bvar**2 - 367416 * dvar) * y**2
        + (2361960 * bvar * cvar + 2480058 * wvar) * y
    )
    d_incidence = (
        24300 * bvar**2 - 1890 * bvar - 245
    ) / sp.Integer(122472)
    w_incidence = (
        -sp.Rational(20, 21) * bvar * cvar
        - sp.Rational(2, 27) * cvar
        - (91 + 1215 * bvar) / sp.Integer(177147)
    )
    p_incidence = sp.expand(quartic_p.subs(dvar, d_incidence))
    q_incidence = sp.expand(
        sextic_q.subs({dvar: d_incidence, wvar: w_incidence})
    )
    smooth = 245 + 2835 * bvar + 26244 * cvar
    require(
        p_incidence.subs(y, 1) == 0
        and q_incidence.subs(y, 1) == 0
        and sp.expand(sp.diff(q_incidence, y).subs(y, 1) - 14 * smooth)
        == 0,
        "moving-root incidence changed",
    )

    mordell = sp.expand(4 * p_incidence**3 + 49 * q_incidence**2)
    residual, remainder = sp.div(mordell, (y - 1) ** 2, y)
    require(remainder == 0 and sp.degree(residual, y) == 10, "bad R10")
    residual = sp.expand(residual)
    leading = sp.Poly(residual, y).LC()
    require(leading == 73060029, "R10 leading coefficient changed")

    square_factor = y**4 + s3 * y**3 + s2 * y**2 + s1 * y + s0
    branch_factor = leading * y**2 + h1 * y + h0
    difference = sp.expand(branch_factor * square_factor**2 - residual)

    h1_solution = solve_linear(coefficient(difference, y, 9), h1)
    difference = sp.expand(difference.subs(h1, h1_solution))
    h0_solution = solve_linear(coefficient(difference, y, 8), h0)
    difference = sp.expand(difference.subs(h0, h0_solution))
    s1_solution = solve_linear(coefficient(difference, y, 7), s1)
    difference = sp.expand(difference.subs(s1, s1_solution))
    s0_solution = solve_linear(coefficient(difference, y, 6), s0)
    difference = sp.expand(difference.subs(s0, s0_solution))

    require(
        all(coefficient(difference, y, degree) == 0 for degree in range(6, 11)),
        "top-down coefficient recursion failed",
    )
    variables = (bvar, cvar, s3, s2)
    residual_equations = [
        primitive(coefficient(difference, y, degree), variables)
        for degree in range(5, -1, -1)
    ]
    require(len(residual_equations) == 6, "residual equation count changed")

    branch_discriminant = sp.factor(
        sp.together(h1_solution**2 - 4 * leading * h0_solution)
    )
    branch_discriminant = primitive(branch_discriminant, variables)

    common_wall = 54 * bvar + 7
    double_zero_wall = 13122 * cvar + 1215 * bvar + 91
    localization = sp.expand(
        cvar
        * smooth
        * common_wall
        * double_zero_wall
        * branch_discriminant
    )
    return (
        bvar,
        cvar,
        s3,
        s2,
        smooth,
        localization,
        residual_equations,
    )


def signatures(
    expressions: list[sp.Expr],
    variables: tuple[sp.Symbol, ...],
) -> list[tuple[int, int, int, int, int, int]]:
    bvar, cvar, s3, s2 = variables
    return [
        (
            sp.total_degree(item),
            sp.degree(item, bvar),
            sp.degree(item, cvar),
            sp.degree(item, s3),
            sp.degree(item, s2),
            len(sp.Poly(item, *variables).terms()),
        )
        for item in expressions
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode",
        choices=(
            "profile",
            "reduced-profile",
            "reduced-modular",
            "bruteforce",
            "modular-unsaturated",
            "modular",
            "full",
        ),
        default="reduced-profile",
    )
    parser.add_argument("--primes", default="29,31,37,41,43")
    args = parser.parse_args()

    started = time.perf_counter()
    (
        bvar,
        cvar,
        s3,
        s2,
        smooth,
        localization,
        residual_equations,
    ) = build()
    variables = (bvar, cvar, s3, s2)
    print("THM-2357 degree-18 H2 factorization-pivot exact companion")
    print("gauge: smooth e=3 branch at y=1")
    print(f"smoothness factor: {smooth}")
    print("top-down solved variables: h1,h0,s1,s0")
    print(f"residual equation signatures: {signatures(residual_equations, variables)}")
    print(
        "localization factors: "
        "C*K*(54B+7)*(13122C+1215B+91)*Disc(H2)"
    )
    print("factorization build completed exactly")
    if args.mode == "profile":
        return

    if args.mode in ("reduced-profile", "reduced-modular"):
        first_equation = sp.Poly(
            residual_equations[0],
            cvar,
        )
        require(first_equation.degree() == 1, "first residual is no longer C-linear")
        c_leading = first_equation.LC()
        c_constant = first_equation.TC()
        c_solution = sp.cancel(-c_constant / c_leading)
        reduced_variables = (bvar, s3, s2)
        reduced_equations: list[sp.Expr] = []
        reduced_denominators: list[sp.Expr] = []
        for equation in residual_equations[1:]:
            numerator, denominator = sp.together(
                equation.subs(cvar, c_solution)
            ).as_numer_denom()
            content, primitive_numerator = sp.Poly(
                numerator,
                *reduced_variables,
            ).primitive()
            require(content != 0, "C-eliminated residual met zero")
            reduced_equations.append(primitive_numerator.as_expr())
            reduced_denominators.append(sp.factor(denominator))
        c_factorization = sp.factor_list(c_leading, bvar, s3, s2)
        boundary_b = sp.Rational(1771, 10260) * (s2 - 1)
        boundary_numerator = sp.together(
            c_constant.subs(bvar, boundary_b)
        ).as_numer_denom()[0]
        boundary_content, boundary_polynomial = sp.Poly(
            boundary_numerator,
            s3,
            s2,
        ).primitive()
        boundary_factorization = sp.factor_list(
            boundary_polynomial.as_expr(),
            s3,
            s2,
        )
        reduced_signatures = [
            (
                sp.total_degree(item),
                sp.degree(item, bvar),
                sp.degree(item, s3),
                sp.degree(item, s2),
                len(sp.Poly(item, *reduced_variables).terms()),
            )
            for item in reduced_equations
        ]
        print(
            "C-leading factor signatures/exponents: "
            + str(
                [
                    (
                        (
                            sp.total_degree(factor),
                            sp.degree(factor, bvar),
                            sp.degree(factor, s3),
                            sp.degree(factor, s2),
                            len(sp.Poly(factor, *reduced_variables).terms()),
                        ),
                        exponent,
                    )
                    for factor, exponent in c_factorization[1]
                ]
            )
        )
        print(f"C-leading content: {c_factorization[0]}")
        print(
            "C-leading wall: -10260*B+1771*s2-1771=0; "
            "on it the first residual constant has "
            "factor signatures/exponents="
            + str(
                [
                    (
                        (
                            sp.total_degree(factor),
                            sp.degree(factor, s3),
                            sp.degree(factor, s2),
                            len(sp.Poly(factor, s3, s2).terms()),
                        ),
                        exponent,
                    )
                    for factor, exponent in boundary_factorization[1]
                ]
            )
            + f", content={boundary_content * boundary_factorization[0]}"
        )
        print(
            "C-leading-wall residual polynomial: "
            f"{boundary_polynomial.as_expr()}"
        )
        print(f"C-eliminated signatures: {reduced_signatures}")
        for index, equation in enumerate(reduced_equations[:2]):
            factor_started = time.perf_counter()
            factorization = sp.factor_list(
                equation,
                *reduced_variables,
            )
            print(
                f"reduced factor {index}: content={factorization[0]}, "
                "factor signatures/exponents="
                + str(
                    [
                        (
                            (
                                sp.total_degree(factor),
                                sp.degree(factor, bvar),
                                sp.degree(factor, s3),
                                sp.degree(factor, s2),
                                len(
                                    sp.Poly(
                                        factor,
                                        *reduced_variables,
                                    ).terms()
                                ),
                            ),
                            exponent,
                        )
                        for factor, exponent in factorization[1]
                    ]
                )
                + ", factorization exact"
            )
        print(
            "distinct reduced denominators: "
            + str(list(dict.fromkeys(map(str, reduced_denominators))))
        )
        if args.mode == "reduced-modular":
            for prime_text in args.primes.split(","):
                prime = int(prime_text)
                require(prime > 7 and sp.isprime(prime), f"bad prime {prime}")
                basis_started = time.perf_counter()
                basis = sp.groebner(
                    reduced_equations,
                    bvar,
                    s3,
                    s2,
                    order="grevlex",
                    method="f5b",
                    modulus=prime,
                )
                expressions = [item.as_expr() for item in basis.polys]
                unit = len(expressions) == 1 and expressions[0] == 1
                print(
                    f"C-eliminated mod {prime}: "
                    f"basis count={len(expressions)}, unit={unit}, "
                    f"seconds={time.perf_counter() - basis_started:.3f}"
                )
        return

    inverse = sp.symbols("U")
    if args.mode == "bruteforce":
        import numpy as np

        def evaluate_mod(
            expression: sp.Expr,
            prime: int,
            grids: tuple[np.ndarray, ...],
        ) -> np.ndarray:
            polynomial = sp.Poly(expression, *variables)
            maximum_powers = [
                max(monomial[index] for monomial, _ in polynomial.terms())
                for index in range(4)
            ]
            power_tables = [
                [
                    np.ones_like(grid)
                    if exponent == 0
                    else np.mod(np.power(grid, exponent), prime)
                    for exponent in range(maximum + 1)
                ]
                for grid, maximum in zip(grids, maximum_powers)
            ]
            value = np.zeros_like(grids[0])
            for monomial, rational_coefficient in polynomial.terms():
                require(
                    rational_coefficient.q % prime != 0,
                    f"denominator vanished modulo {prime}",
                )
                term = np.full_like(
                    value,
                    (
                        int(rational_coefficient.p)
                        * pow(int(rational_coefficient.q), -1, prime)
                    )
                    % prime,
                )
                for index, exponent in enumerate(monomial):
                    term = np.mod(
                        term * power_tables[index][exponent],
                        prime,
                    )
                value = np.mod(value + term, prime)
            return value

        for prime_text in args.primes.split(","):
            prime = int(prime_text)
            require(prime > 7 and sp.isprime(prime), f"bad prime {prime}")
            scan_started = time.perf_counter()
            axes = np.arange(prime, dtype=np.int64)
            grids = tuple(
                item.reshape(-1)
                for item in np.meshgrid(
                    axes,
                    axes,
                    axes,
                    axes,
                    indexing="ij",
                )
            )
            common_mask = np.ones(grids[0].shape, dtype=bool)
            for equation in residual_equations:
                common_mask &= evaluate_mod(equation, prime, grids) == 0
            raw_count = int(np.count_nonzero(common_mask))
            localization_values = evaluate_mod(localization, prime, grids)
            localized_mask = common_mask & (localization_values != 0)
            localized_count = int(np.count_nonzero(localized_mask))
            sample_indices = np.flatnonzero(localized_mask)[:8]
            samples = [
                tuple(int(grid[index]) for grid in grids)
                for index in sample_indices
            ]
            print(
                f"bruteforce mod {prime}: raw={raw_count}, "
                f"localized={localized_count}, samples={samples}"
            )
        return

    if args.mode == "modular-unsaturated":
        for prime_text in args.primes.split(","):
            prime = int(prime_text)
            require(prime > 7 and sp.isprime(prime), f"bad prime {prime}")
            basis_started = time.perf_counter()
            basis = sp.groebner(
                residual_equations,
                bvar,
                cvar,
                s3,
                s2,
                order="grevlex",
                method="f5b",
                modulus=prime,
            )
            expressions = [item.as_expr() for item in basis.polys]
            unit = len(expressions) == 1 and expressions[0] == 1
            print(
                f"unsaturated mod {prime}: basis count={len(expressions)}, "
                f"unit={unit}, seconds={time.perf_counter() - basis_started:.3f}"
            )
        return

    if args.mode == "modular":
        for prime_text in args.primes.split(","):
            prime = int(prime_text)
            require(prime > 7 and sp.isprime(prime), f"bad prime {prime}")
            basis_started = time.perf_counter()
            basis = sp.groebner(
                residual_equations + [inverse * localization - 1],
                inverse,
                bvar,
                cvar,
                s3,
                s2,
                order="grevlex",
                method="f5b",
                modulus=prime,
            )
            expressions = [item.as_expr() for item in basis.polys]
            unit = len(expressions) == 1 and expressions[0] == 1
            print(
                f"mod {prime}: basis count={len(expressions)}, "
                f"unit={unit}, seconds={time.perf_counter() - basis_started:.3f}"
            )
        return

    basis_started = time.perf_counter()
    basis = sp.groebner(
        residual_equations + [inverse * localization - 1],
        inverse,
        bvar,
        cvar,
        s3,
        s2,
        order="grevlex",
        method="f5b",
    )
    expressions = [item.as_expr() for item in basis.polys]
    unit = len(expressions) == 1 and expressions[0] == 1
    print(f"exact basis count: {len(expressions)}")
    print(f"exact localized ideal is unit: {unit}")
    print(f"exact basis seconds: {time.perf_counter() - basis_started:.3f}")


if __name__ == "__main__":
    main()
