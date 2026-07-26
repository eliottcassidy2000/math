#!/usr/bin/env python3
"""Exact moving-root and subresultant companion for THM-2357.

THM-2332 shows that an H_2 S_5^2 survivor has branch signature
(3,2,2).  Infinity is unramified, so its unique three-cycle branch value
is finite.  At that value the depressed cubic has a smooth triple fibre,
and hence P=Q=0 and Q'!=0.

Use weighted scaling to move that branch value to y=1.  Solving P(1)=Q(1)=0
eliminates D and W, leaving a degree-ten polynomial

    R_10 = (4 P^3 + 49 Q^2)/(y-1)^2

in the two parameters B,C.  The H_2 condition forces

    R_10 = H_2 S_4^2,

so gcd(R_10,R_10') has degree at least four.  This script forms all ten
coefficients of Sres_0,...,Sres_3 and probes their ideal after localizing
away from the smoothness factor and the already-closed parameter walls.

No executable check uses Python assert.
"""

from __future__ import annotations

import argparse
import time

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def primitive(expression: sp.Expr, variables: tuple[sp.Symbol, ...]) -> sp.Expr:
    """Clear rational content and return a primitive integral polynomial."""

    numerator, denominator = sp.together(expression).as_numer_denom()
    require(
        denominator.is_Integer and denominator != 0,
        "coefficient acquired a nonconstant denominator",
    )
    content, result = sp.Poly(numerator, *variables).primitive()
    require(content != 0, "primitive part met the zero polynomial")
    return result.as_expr()


def signatures(
    expressions: list[sp.Expr],
    variables: tuple[sp.Symbol, ...],
) -> list[tuple[int, int, int, int]]:
    """Return (total degree, B degree, C degree, term count)."""

    bvar, cvar = variables
    return [
        (
            sp.total_degree(item),
            sp.degree(item, bvar),
            sp.degree(item, cvar),
            len(sp.Poly(item, bvar, cvar).terms()),
        )
        for item in expressions
    ]


def strip_local_factors(
    expression: sp.Expr,
    factors: tuple[sp.Expr, ...],
    variables: tuple[sp.Symbol, ...],
) -> tuple[sp.Expr, tuple[int, ...]]:
    """Strip maximal powers of factors inverted in the localized chart."""

    polynomial = sp.Poly(expression, *variables)
    orders: list[int] = []
    for factor in factors:
        factor_polynomial = sp.Poly(factor, *variables)
        order = 0
        while True:
            quotient, remainder = sp.div(polynomial, factor_polynomial)
            if not remainder.is_zero:
                break
            polynomial = quotient
            order += 1
        orders.append(order)
    content, result = polynomial.primitive()
    require(content != 0, "localized factor stripping met zero")
    return result.as_expr(), tuple(orders)


def build() -> tuple[
    sp.Symbol,
    sp.Symbol,
    sp.Symbol,
    sp.Expr,
    sp.Expr,
    sp.Expr,
    list[sp.Expr],
]:
    """Build R_10, localization factors, and the ten terminal coefficients."""

    y = sp.symbols("y")
    bvar, cvar, dvar, wvar = sp.symbols("B C D W")

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

    # The unique smooth three-cycle branch is moved to y=1.
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
        sp.expand(p_incidence.subs(y, 1)) == 0
        and sp.expand(q_incidence.subs(y, 1)) == 0,
        "moving-root incidence stopped imposing P(1)=Q(1)=0",
    )
    require(
        sp.expand(sp.diff(q_incidence, y).subs(y, 1) - 14 * smooth) == 0,
        "smooth triple-fibre factor changed",
    )

    mordell = sp.expand(4 * p_incidence**3 + 49 * q_incidence**2)
    residual, remainder = sp.div(mordell, (y - 1) ** 2, y)
    require(remainder == 0, "Mordell polynomial lost its exact double zero")
    residual = sp.expand(residual)
    require(
        sp.degree(residual, y) == 10
        and sp.expand(residual.subs(y, 1) - 49 * (14 * smooth) ** 2) == 0,
        "residual degree or smooth-fibre value changed",
    )

    subresultants = sp.subresultants(residual, sp.diff(residual, y), y)
    degree_profile = [sp.degree(item, y) for item in subresultants]
    require(
        degree_profile == list(range(10, -1, -1)),
        f"unexpected generic subresultant degree profile: {degree_profile}",
    )
    by_degree = {
        sp.degree(item, y): sp.Poly(item, y)
        for item in subresultants
    }

    terminal_coefficients: list[sp.Expr] = []
    for degree in range(4):
        coefficient_list = by_degree[degree].all_coeffs()
        require(
            len(coefficient_list) == degree + 1,
            f"Sres_{degree} changed shape",
        )
        terminal_coefficients.extend(
            primitive(item, (bvar, cvar)) for item in coefficient_list
        )
    require(
        len(terminal_coefficients) == 10,
        "terminal coefficient count changed",
    )

    # In the y=1 gauge the previously closed walls become linear:
    #
    # 126D-25B^2 = -(35/972)(54B+7),
    # 20BC+21W   = -(7/59049)(13122C+1215B+91).
    common_wall = 54 * bvar + 7
    double_zero_wall = 13122 * cvar + 1215 * bvar + 91
    require(
        sp.expand(
            126 * d_incidence
            - 25 * bvar**2
            + sp.Rational(35, 972) * common_wall
        )
        == 0,
        "common-root-wall pullback changed",
    )
    require(
        sp.expand(
            20 * bvar * cvar
            + 21 * w_incidence
            + sp.Rational(7, 59049) * double_zero_wall
        )
        == 0,
        "double-zero-wall pullback changed",
    )

    # C=0 is kept out of the main run as a requested boundary chart.
    localization = sp.expand(
        cvar * smooth * common_wall * double_zero_wall
    )
    return (
        y,
        bvar,
        cvar,
        residual,
        smooth,
        localization,
        terminal_coefficients,
    )


def gcd_degree(
    residual: sp.Expr,
    y: sp.Symbol,
    bvar: sp.Symbol,
    cvar: sp.Symbol,
    bvalue: sp.Expr,
    cvalue: sp.Expr,
) -> int:
    specialization = sp.Poly(
        residual.subs({bvar: bvalue, cvar: cvalue}),
        y,
        domain=sp.QQ,
    )
    return sp.gcd(specialization, specialization.diff()).degree()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode",
        choices=(
            "profile",
            "factors",
            "modular",
            "unsaturated",
            "selected",
            "full",
        ),
        default="profile",
    )
    parser.add_argument(
        "--selected-count",
        type=int,
        default=6,
        help="number of smallest terminal coefficients in selected mode",
    )
    parser.add_argument(
        "--primes",
        default="29,31,37,41,43",
        help="comma-separated primes for modular mode",
    )
    args = parser.parse_args()

    started = time.perf_counter()
    (
        y,
        bvar,
        cvar,
        residual,
        smooth,
        localization,
        terminal_coefficients,
    ) = build()
    built_seconds = time.perf_counter() - started

    coefficient_signatures = signatures(
        terminal_coefficients,
        (bvar, cvar),
    )
    common_wall = 54 * bvar + 7
    double_zero_wall = 13122 * cvar + 1215 * bvar + 91
    local_factors = (cvar, smooth, common_wall, double_zero_wall)
    stripped_with_orders = [
        strip_local_factors(
            item,
            local_factors,
            (bvar, cvar),
        )
        for item in terminal_coefficients
    ]
    stripped_coefficients = [item[0] for item in stripped_with_orders]
    stripped_orders = [item[1] for item in stripped_with_orders]
    stripped_signatures = signatures(
        stripped_coefficients,
        (bvar, cvar),
    )
    ranked = sorted(
        zip(
            stripped_coefficients,
            stripped_signatures,
            range(len(terminal_coefficients)),
        ),
        key=lambda item: (item[1][3], item[1][0], item[2]),
    )

    generic_controls = (
        gcd_degree(residual, y, bvar, cvar, 0, 1),
        gcd_degree(residual, y, bvar, cvar, 1, 1),
        gcd_degree(residual, y, bvar, cvar, 2, -1),
    )
    require(
        generic_controls == (0, 0, 0),
        f"generic hostile gcd controls changed: {generic_controls}",
    )

    # K=0 is outside the smooth e=3 chart.  After removing (y-1)^2, the
    # generic order-three discriminant zero leaves a *simple* y=1 factor,
    # so this hostile control should not falsely pass the gcd-four test.
    k_zero_c = -sp.Rational(245, 26244)
    k_zero_degree = gcd_degree(
        residual,
        y,
        bvar,
        cvar,
        0,
        k_zero_c,
    )
    require(
        sp.expand(
            residual.subs(
                {
                    bvar: 0,
                    cvar: k_zero_c,
                    y: 1,
                }
            )
        )
        == 0
        and k_zero_degree == 0,
        "K=0 hostile control changed its simple residual collision",
    )

    print("THM-2357 degree-18 H2 moving-root exact companion")
    print("gauge: unique smooth e=3 branch moved to y=1")
    print("incidence: P(1)=Q(1)=0 and Q'(1)=14*K")
    print(f"smoothness factor K: {smooth}")
    print("residual: R10=(4P^3+49Q^2)/(y-1)^2")
    print("subresultant target: all coefficients of Sres_0,...,Sres_3")
    print(f"terminal coefficient count: {len(terminal_coefficients)}")
    print(f"terminal signatures: {coefficient_signatures}")
    print(
        "stripped factor order columns: "
        "(C,K,54B+7,13122C+1215B+91)"
    )
    print(f"stripped factor orders: {stripped_orders}")
    print(f"stripped signatures: {stripped_signatures}")
    print(
        "ranked coefficient indices/signatures: "
        + str([(item[2], item[1]) for item in ranked])
    )
    print(
        "localization: C*K*(54B+7)*(13122C+1215B+91)"
    )
    print(f"generic hostile gcd degrees: {generic_controls}")
    print(f"K=0 hostile gcd degree: {k_zero_degree}")
    print("profile build completed exactly")

    if args.mode == "profile":
        return

    if args.mode == "factors":
        for expression, signature, index in ranked[:4]:
            factor_started = time.perf_counter()
            content, factor_list = sp.factor_list(expression, bvar, cvar)
            print(
                f"factor index {index}: content={content}, "
                f"factor signatures/exponents="
                f"{[(signatures([factor], (bvar, cvar))[0], exponent) for factor, exponent in factor_list]}, "
                f"seconds={time.perf_counter() - factor_started:.3f}"
            )
        return

    inverse = sp.symbols("U")
    if args.mode == "modular":
        prime_list = [int(item) for item in args.primes.split(",") if item]
        for prime in prime_list:
            require(prime > 7 and sp.isprime(prime), f"bad probe prime {prime}")
            prime_started = time.perf_counter()
            smallest_unit: int | None = None
            for count in range(len(ranked), 1, -1):
                modular_basis = sp.groebner(
                    [item[0] for item in ranked[:count]]
                    + [inverse * localization - 1],
                    inverse,
                    bvar,
                    cvar,
                    order="grevlex",
                    method="f5b",
                    modulus=prime,
                )
                modular_expressions = [
                    item.as_expr() for item in modular_basis.polys
                ]
                if (
                    len(modular_expressions) == 1
                    and modular_expressions[0] == 1
                ):
                    smallest_unit = count
                    continue
                break
            print(
                f"mod {prime}: smallest tested localized unit prefix="
                f"{smallest_unit}, "
                f"seconds={time.perf_counter() - prime_started:.3f}"
            )
        return

    if args.mode == "unsaturated":
        groebner_started = time.perf_counter()
        unsaturated_basis = sp.groebner(
            stripped_coefficients,
            bvar,
            cvar,
            order="grevlex",
            method="f5b",
        )
        unsaturated_expressions = [
            item.as_expr() for item in unsaturated_basis.polys
        ]
        unit = (
            len(unsaturated_expressions) == 1
            and unsaturated_expressions[0] == 1
        )
        print(f"unsaturated basis count: {len(unsaturated_expressions)}")
        print(
            "unsaturated basis signatures: "
            + str(
                [
                    (
                        sp.total_degree(item),
                        sp.degree(item, bvar),
                        sp.degree(item, cvar),
                        len(sp.Poly(item, bvar, cvar).terms()),
                    )
                    for item in unsaturated_expressions
                ]
            )
        )
        print(f"unsaturated ideal is unit: {unit}")
        print(
            "unsaturated Groebner seconds: "
            f"{time.perf_counter() - groebner_started:.3f}"
        )
        return

    if args.mode == "selected":
        count = max(1, min(args.selected_count, len(ranked)))
        generators = [item[0] for item in ranked[:count]]
        print(
            "selected generator indices: "
            + str([item[2] for item in ranked[:count]])
        )
    else:
        generators = stripped_coefficients

    groebner_started = time.perf_counter()
    basis = sp.groebner(
        generators + [inverse * localization - 1],
        inverse,
        bvar,
        cvar,
        order="grevlex",
        method="f5b",
    )
    groebner_seconds = time.perf_counter() - groebner_started
    basis_expressions = [item.as_expr() for item in basis.polys]
    unit = len(basis_expressions) == 1 and basis_expressions[0] == 1
    basis_signatures = [
        (
            sp.total_degree(item),
            sp.degree(item, inverse),
            sp.degree(item, bvar),
            sp.degree(item, cvar),
            len(sp.Poly(item, inverse, bvar, cvar).terms()),
        )
        for item in basis_expressions
    ]
    print(f"Groebner generator count: {len(generators) + 1}")
    print(f"Groebner basis count: {len(basis_expressions)}")
    print(f"Groebner basis signatures: {basis_signatures}")
    print(f"localized ideal is unit: {unit}")
    print(f"Groebner seconds: {groebner_seconds:.3f}")


if __name__ == "__main__":
    main()
