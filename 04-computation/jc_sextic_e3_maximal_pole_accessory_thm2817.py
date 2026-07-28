#!/usr/bin/env python3
"""Exact controls for THM-2817's sextic e=3 maximal-pole atlas.

The universe is the ten ordered positive four-part compositions of six.
For each passport we form the two critical-value remainder equations,
saturate every collision/zero-value boundary, and compare the exact
elimination ideal with the claimed two-point table.  A separate permutation
enumeration checks the marked Nielsen count.
"""

import ast
from itertools import permutations
from pathlib import Path

import sympy as sp


def require(condition, message):
    if not bool(condition):
        raise RuntimeError(message)


def has_asserts(path):
    tree = ast.parse(Path(path).read_text(encoding="utf-8"))
    return any(isinstance(node, ast.Assert) for node in ast.walk(tree))


X, LAM, MU, INV = sp.symbols("x lambda mu inv")


EXPECTED = {
    (1, 1, 1, 3): (3 * MU - LAM - 1, LAM**2 - LAM + 1),
    (1, 1, 2, 2): (MU + LAM - 1, 16 * LAM**2 - 16 * LAM + 3),
    (1, 1, 3, 1): (MU - 3 * LAM + 1, 3 * LAM**2 - 3 * LAM + 1),
    (1, 2, 1, 2): (MU - LAM + 1, 3 * LAM**2 - 16 * LAM + 16),
    (1, 2, 2, 1): (MU - LAM - 1, 3 * LAM**2 - 10 * LAM + 3),
    (1, 3, 1, 1): (MU + LAM - 3, LAM**2 - 3 * LAM + 3),
    (2, 1, 1, 2): (MU - LAM - 1, 3 * LAM**2 + 10 * LAM + 3),
    (2, 1, 2, 1): (MU - LAM + 1, 3 * LAM**2 + 4 * LAM - 4),
    (2, 2, 1, 1): (MU + LAM - 1, 4 * LAM**2 - 4 * LAM - 3),
    (3, 1, 1, 1): (MU + LAM + 1, LAM**2 + LAM + 1),
}


def exact_data(parts):
    a, b, c, d = parts
    denominator = sp.expand(
        X**a * (X - 1) ** b * (X - LAM) ** c * (X - MU) ** d
    )
    forced_pole_factor = sp.expand(
        X ** (a - 1)
        * (X - 1) ** (b - 1)
        * (X - LAM) ** (c - 1)
        * (X - MU) ** (d - 1)
    )
    critical = sp.cancel(sp.diff(denominator, X) / (6 * forced_pole_factor))
    require(sp.Poly(critical, X).degree() == 3, "critical cubic degree")
    require(sp.Poly(critical, X).LC() == 1, "critical cubic monic")

    field = sp.QQ.frac_field(LAM, MU)
    remainder = sp.rem(
        sp.Poly(denominator, X, domain=field),
        sp.Poly(critical, X, domain=field),
    )
    r2 = sp.cancel(remainder.coeff_monomial(X**2))
    r1 = sp.cancel(remainder.coeff_monomial(X))
    value = sp.cancel(remainder.coeff_monomial(1))
    r2_num = sp.factor(r2.as_numer_denom()[0])
    r1_num = sp.factor(r1.as_numer_denom()[0])

    boundary = sp.factor(
        LAM
        * MU
        * (LAM - 1)
        * (MU - 1)
        * (LAM - MU)
        * sp.discriminant(critical, X)
        * value
    )
    return {
        "D": denominator,
        "P": forced_pole_factor,
        "E": sp.expand(critical),
        "r2": r2_num,
        "r1": r1_num,
        "v": value,
        "boundary": boundary,
    }


def zero_in_expected_quotient(expr, linear, quadratic):
    mu_value = sp.solve(linear, MU)[0]
    substituted = sp.cancel(expr.subs(MU, mu_value))
    numerator, denominator = substituted.as_numer_denom()
    qpoly = sp.Poly(quadratic, LAM, domain=sp.QQ)
    require(
        sp.gcd(sp.Poly(denominator, LAM, domain=sp.QQ), qpoly).degree() == 0,
        "quotient denominator is not a unit",
    )
    return sp.rem(
        sp.Poly(numerator, LAM, domain=sp.QQ)
        * sp.invert(sp.Poly(denominator, LAM, domain=sp.QQ), qpoly),
        qpoly,
    ).is_zero


def polynomial_zero_in_expected_quotient(expr, linear, quadratic):
    return all(
        zero_in_expected_quotient(coefficient, linear, quadratic)
        for coefficient in sp.Poly(sp.expand(expr), X).all_coeffs()
    )


def check_accessory_table():
    rows = []
    checks = 0
    for parts, (linear, quadratic) in sorted(EXPECTED.items()):
        data = exact_data(parts)

        saturated = sp.groebner(
            [data["r2"], data["r1"], 1 - INV * data["boundary"]],
            INV,
            MU,
            LAM,
            order="lex",
            domain=sp.QQ,
        )
        elimination = [
            poly.as_expr()
            for poly in saturated.polys
            if not poly.as_expr().has(INV)
        ]
        claimed = sp.groebner(
            [linear, quadratic], MU, LAM, order="lex", domain=sp.QQ
        )
        require(elimination, f"empty elimination basis for {parts}")
        for generator in elimination:
            require(
                claimed.reduce(generator)[1] == 0,
                f"saturated ideal not contained in table for {parts}",
            )
            checks += 1
        for generator in (linear, quadratic):
            require(
                saturated.reduce(generator)[1] == 0,
                f"table not contained in saturated ideal for {parts}",
            )
            checks += 1

        qpoly = sp.Poly(quadratic, LAM, domain=sp.QQ)
        require(qpoly.degree() == 2, f"table length for {parts}")
        require(sp.discriminant(qpoly.as_expr(), LAM) != 0, f"radical for {parts}")
        checks += 2

        mu_value = sp.solve(linear, MU)[0]
        boundary_sub = sp.cancel(data["boundary"].subs(MU, mu_value))
        boundary_num = boundary_sub.as_numer_denom()[0]
        require(
            sp.gcd(sp.Poly(boundary_num, LAM, domain=sp.QQ), qpoly).degree()
            == 0,
            f"boundary survives at all roots for {parts}",
        )
        checks += 1

        require(
            sp.expand(sp.diff(data["D"], X) - 6 * data["P"] * data["E"])
            == 0,
            f"logarithmic derivative for {parts}",
        )
        require(
            polynomial_zero_in_expected_quotient(
                data["E"] ** 2 - data["D"] + data["v"],
                linear,
                quadratic,
            ),
            f"response reconstruction for {parts}",
        )
        checks += 2
        rows.append((parts, sp.factor(linear), sp.factor(quadratic)))
    return rows, checks


def all_matchings(vertices):
    if not vertices:
        yield ()
        return
    first = vertices[0]
    for index in range(1, len(vertices)):
        second = vertices[index]
        rest = vertices[1:index] + vertices[index + 1 :]
        for tail in all_matchings(rest):
            yield tuple(sorted(((first, second),) + tail))


def chords_cross(edge_one, edge_two):
    a, b = sorted(edge_one)
    c, d = sorted(edge_two)
    return (a < c < b < d) or (c < a < d < b)


def noncrossing(matching):
    return all(
        not chords_cross(matching[i], matching[j])
        for i in range(len(matching))
        for j in range(i)
    )


def permutation_cycles(permutation):
    seen = set()
    cycles = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        cycle = []
        current = start
        while current not in seen:
            seen.add(current)
            cycle.append(current)
            current = permutation[current]
        cycles.append(tuple(cycle))
    return cycles


def rotated_state(matching, labels, shift):
    degree = 6
    shifted_matching = tuple(
        sorted(
            tuple(sorted(((a + shift) % degree, (b + shift) % degree)))
            for a, b in matching
        )
    )
    shifted_labels = tuple(labels[(i - shift) % degree] for i in range(degree))
    return shifted_matching, shifted_labels


def check_marked_nielsen():
    degree = 6
    matchings = [
        matching
        for matching in all_matchings(tuple(range(degree)))
        if noncrossing(matching)
    ]
    require(len(matchings) == 5, "Catalan perfect matching count")

    states = {}
    for matching in matchings:
        tau = {}
        for a, b in matching:
            tau[a] = b
            tau[b] = a
        sigma = tuple(tau[(i + 1) % degree] for i in range(degree))
        cycles = permutation_cycles(sigma)
        require(len(cycles) == 4, "noncrossing matching has four pole cycles")
        for order in permutations(range(4)):
            vertex_labels = [None] * degree
            profile = []
            for label, cycle_index in enumerate(order):
                cycle = cycles[cycle_index]
                profile.append(len(cycle))
                for vertex in cycle:
                    vertex_labels[vertex] = label
            canonical = min(
                rotated_state(matching, tuple(vertex_labels), shift)
                for shift in range(degree)
            )
            states[canonical] = tuple(profile)

    counts = {profile: 0 for profile in EXPECTED}
    for profile in states.values():
        require(profile in counts, "unexpected ordered pole profile")
        counts[profile] += 1
    require(len(states) == 20, "total marked Nielsen state count")
    require(all(count == 2 for count in counts.values()), "two charts per passport")
    return counts


def check_symmetric_carriers():
    y = sp.symbols("y")
    d_power = sp.expand(X**3 * (X**3 - 1))
    e_power = X**3 - sp.Rational(1, 2)
    v_power = -sp.Rational(1, 4)
    require(sp.expand(e_power**2 - d_power + v_power) == 0, "cubic power carrier")

    d_cheb = sp.expand((y**2 - sp.Rational(1, 4)) ** 2 * (y**2 - 1))
    e_cheb = y**3 - sp.Rational(3, 4) * y
    v_cheb = -sp.Rational(1, 16)
    require(sp.expand(e_cheb**2 - d_cheb + v_cheb) == 0, "Chebyshev carrier")
    require(
        sp.expand(4 * e_cheb - (4 * y**3 - 3 * y)) == 0,
        "T3 normalization",
    )


def main():
    path = Path(__file__)
    require(not has_asserts(path), "truth-bearing assert node")
    rows, algebra_checks = check_accessory_table()
    nielsen_counts = check_marked_nielsen()
    check_symmetric_carriers()

    print("THM-2817 sextic e=3 maximal-pole exact controls")
    for parts, linear, quadratic in rows:
        print(f"{parts}: {linear}; {quadratic}; charts={nielsen_counts[parts]}")
    print("symmetric carriers: cubic power and T3/4 PASS")
    print(f"PASS passports={len(rows)} algebra_gates={algebra_checks} marked_charts=20")


if __name__ == "__main__":
    main()
