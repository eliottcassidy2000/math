#!/usr/bin/env python3
"""Exact saturated accessory and Nielsen controls for THM-2817."""

import ast
from collections import defaultdict
from itertools import permutations
from pathlib import Path

import sympy as sp


def require(condition, message):
    if not bool(condition):
        raise RuntimeError(message)


def has_asserts(path):
    return any(
        isinstance(node, ast.Assert)
        for node in ast.walk(ast.parse(Path(path).read_text(encoding="utf-8")))
    )


def reduce_mod(expr, variable, modulus):
    numerator, denominator = sp.cancel(expr).as_numer_denom()
    modulus_poly = sp.Poly(modulus, variable, domain=sp.QQ)
    denominator_poly = sp.Poly(denominator, variable, domain=sp.QQ)
    require(sp.gcd(denominator_poly, modulus_poly).degree() == 0,
            "nonunit denominator")
    inverse = sp.invert(denominator_poly, modulus_poly)
    return sp.rem(
        sp.Poly(numerator, variable, domain=sp.QQ) * inverse,
        modulus_poly,
    ).as_expr()


def reduce_x(expr, x, parameter, modulus):
    return all(
        reduce_mod(coefficient, parameter, modulus) == 0
        for coefficient in sp.Poly(sp.expand(expr), x).all_coeffs()
    )


def compose(left, right):
    return tuple(left[right[index]] for index in range(len(left)))


def cycles(permutation):
    seen = set()
    answer = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        cycle = []
        point = start
        while point not in seen:
            seen.add(point)
            cycle.append(point)
            point = permutation[point]
        answer.append(tuple(cycle))
    return tuple(answer)


def noncrossing_matchings(points):
    points = tuple(points)
    if not points:
        yield ()
        return
    first = points[0]
    for partner_index in range(1, len(points), 2):
        partner = points[partner_index]
        inside = points[1:partner_index]
        outside = points[partner_index + 1 :]
        for left in noncrossing_matchings(inside):
            for right in noncrossing_matchings(outside):
                yield ((first, partner),) + left + right


def rotate_matching(matching, shift, degree):
    return tuple(
        sorted(
            tuple(sorted(((left + shift) % degree, (right + shift) % degree)))
            for left, right in matching
        )
    )


def labelled_key(matching, labelled_cycles, degree):
    return min(
        (
            rotate_matching(matching, shift, degree),
            tuple(
                tuple(sorted((point + shift) % degree for point in cycle))
                for cycle in labelled_cycles
            ),
        )
        for shift in range(degree)
    )


def nielsen_sextic_atlas():
    degree = 6
    full_cycle = tuple((index + 1) % degree for index in range(degree))
    atlas = defaultdict(set)
    for matching in noncrossing_matchings(range(degree)):
        involution = list(range(degree))
        for left, right in matching:
            involution[left], involution[right] = (
                involution[right],
                involution[left],
            )
        pole_cycles = cycles(compose(tuple(involution), full_cycle))
        require(len(pole_cycles) == 4, "four pole cycles")
        for labelled_cycles in permutations(pole_cycles):
            passport = tuple(len(cycle) for cycle in labelled_cycles)
            atlas[passport].add(
                labelled_key(matching, labelled_cycles, degree)
            )
    require(len(atlas) == 10, "ten ordered sextic passports")
    require(all(len(classes) == 2 for classes in atlas.values()),
            "two Nielsen classes per passport")
    return atlas


def accessory_data(parts, x, lam, mu, inverse):
    a, b, c, d = parts
    degree = 6
    D = sp.expand(
        x**a * (x - 1) ** b * (x - lam) ** c * (x - mu) ** d
    )
    T = x * (x - 1) * (x - lam) * (x - mu)
    K = sp.expand(
        a * (x - 1) * (x - lam) * (x - mu)
        + b * x * (x - lam) * (x - mu)
        + c * x * (x - 1) * (x - mu)
        + d * x * (x - 1) * (x - lam)
    )
    E = sp.expand(K / degree)
    remainder = sp.rem(D, K, x)
    remainder_poly = sp.Poly(remainder, x)
    R2 = remainder_poly.coeff_monomial(x**2)
    R1 = remainder_poly.coeff_monomial(x)
    r0 = remainder_poly.coeff_monomial(1)
    f2 = sp.Poly(R2, lam, mu).clear_denoms()[1].as_expr()
    f1 = sp.Poly(R1, lam, mu).clear_denoms()[1].as_expr()
    A = sp.expand(-r0)
    collision = sp.discriminant(K, x)
    gate = sp.together(
        lam
        * mu
        * (lam - 1)
        * (mu - 1)
        * (lam - mu)
        * collision
        * A
    ).as_numer_denom()[0]

    basis = sp.groebner(
        (f2, f1, inverse * gate - 1),
        inverse,
        mu,
        lam,
        order="lex",
    )
    elimination = [
        polynomial.as_expr()
        for polynomial in basis.polys
        if not polynomial.as_expr().has(inverse)
    ]
    require(len(elimination) == 2, "two-generator saturated eliminant")
    linear = next(poly for poly in elimination if sp.degree(poly, mu) == 1)
    quadratic = next(poly for poly in elimination if sp.degree(poly, mu) == 0)
    linear_poly = sp.Poly(linear, mu)
    mu_expression = sp.cancel(
        -linear_poly.coeff_monomial(1)
        / linear_poly.coeff_monomial(mu)
    )
    Q = sp.Poly(quadratic, lam, domain=sp.QQ).monic().as_expr()
    linear_normal = sp.expand(mu - mu_expression)

    require(sp.degree(Q, lam) == 2, "quadratic accessory field")
    require(sp.discriminant(Q, lam) != 0, "two simple accessory points")
    reduced_gate = sp.cancel(gate.subs(mu, mu_expression))
    gate_numerator = reduced_gate.as_numer_denom()[0]
    require(
        sp.gcd(
            sp.Poly(gate_numerator, lam, domain=sp.QQ),
            sp.Poly(Q, lam, domain=sp.QQ),
        ).degree()
        == 0,
        "saturated gate remains a unit",
    )

    B = sp.expand(D + A)
    S, square_remainder = sp.div(B, E**2, x)
    require(
        reduce_x(square_remainder.subs(mu, mu_expression), x, lam, Q),
        "critical cubic square division",
    )
    require(
        reduce_x((S - 1).subs(mu, mu_expression), x, lam, Q),
        "minimal sextic simple factor",
    )
    pole_core = (
        x ** (a - 1)
        * (x - 1) ** (b - 1)
        * (x - lam) ** (c - 1)
        * (x - mu) ** (d - 1)
    )
    C = sp.expand(-degree * A)
    require(
        sp.expand(
            sp.diff(B, x) * D
            - B * sp.diff(D, x)
            - C * E * pole_core
        )
        == 0,
        "response derivative",
    )
    return {
        "linear": linear_normal,
        "Q": Q,
        "mu_expression": mu_expression,
        "D": D,
        "E": E,
        "A": A,
    }


def main():
    require(not has_asserts(Path(__file__)), "truth-bearing Python assert found")
    x, lam, mu, inverse = sp.symbols("x lambda mu inverse")
    print("THM-2817 SEXTIC E=3 MAXIMAL-POLE ACCESSORY AUDIT")
    print("(a,b,c,d) | saturated linear relation | saturated quadratic")

    passports = sorted(
        set(permutations((3, 1, 1, 1)))
        | set(permutations((2, 2, 1, 1)))
    )
    data = {}
    for parts in passports:
        item = accessory_data(parts, x, lam, mu, inverse)
        data[parts] = item
        print(
            f"{parts} | {sp.factor(item['linear'])}"
            f" | {sp.factor(item['Q'])}"
        )

    atlas = nielsen_sextic_atlas()
    require(set(atlas) == set(passports), "algebra/Nielsen passport agreement")

    power = data[(3, 1, 1, 1)]
    require(
        reduce_x(
            (
                (x - 1) * (x - lam) * (x - power["mu_expression"])
                - (x**3 - 1)
            ),
            x,
            lam,
            power["Q"],
        ),
        "cubic-power pole identity",
    )
    require(
        reduce_x(power["E"].subs(mu, power["mu_expression"]) - (x**3 - sp.Rational(1, 2)),
                 x, lam, power["Q"]),
        "cubic-power critical polynomial",
    )
    require(reduce_mod(power["A"].subs(mu, power["mu_expression"]),
                       lam, power["Q"]) == sp.Rational(1, 4),
            "cubic-power constant")

    y = sp.symbols("y")
    cheb_parts = (2, 2, 1, 1)
    a, b, c, d = cheb_parts
    lam_value = sp.Rational(3, 2)
    mu_value = sp.Rational(-1, 2)
    D_cheb = sp.expand(
        x**a
        * (x - 1) ** b
        * (x - lam_value) ** c
        * (x - mu_value) ** d
    )
    K_cheb = sp.expand(
        a * (x - 1) * (x - lam_value) * (x - mu_value)
        + b * x * (x - lam_value) * (x - mu_value)
        + c * x * (x - 1) * (x - mu_value)
        + d * x * (x - 1) * (x - lam_value)
    )
    E_cheb = sp.expand(K_cheb / 6)
    T3 = 4 * y**3 - 3 * y
    require(
        sp.expand(
            D_cheb.subs(x, y + sp.Rational(1, 2))
            - (y**2 - sp.Rational(1, 4)) ** 2 * (y**2 - 1)
        )
        == 0,
        "Chebyshev pole identity",
    )
    require(
        sp.expand(E_cheb.subs(x, y + sp.Rational(1, 2)) - T3 / 4) == 0,
        "Chebyshev critical polynomial",
    )
    require(sp.expand(E_cheb**2 - D_cheb) == sp.Rational(1, 16),
            "Chebyshev constant")

    print("ordered passports=10, accessory points/passport=2, Nielsen classes=2")
    print("power carrier: P3=2*x^3-1, F=P3^2/(P3^2-1)")
    print("Chebyshev carrier: P3=T3(x-1/2), F=P3^2/(P3^2-1)")
    print("saturated ideals, response identities, and ternary atlas: PASS")
    print("assert_nodes=0")


if __name__ == "__main__":
    main()
