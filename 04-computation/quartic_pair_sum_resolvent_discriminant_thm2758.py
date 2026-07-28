#!/usr/bin/env python3
"""Exact quartic pair-sum/resolvent discriminant audit for THM-2758.

The six pair sums split into three opposite pairs around e1/2.  Their centered
squares are a translate of the three standard cubic-resolvent roots.  The
pair-sum Vandermonde has twelve adjacent factors and three opposite factors,
giving disc(G)=disc(f)^2*(e1^3-4*e1*e2+8*e3)^2.  All checks are exact, use a
small custom polynomial ring plus integer/Fraction arithmetic, and use no
truth-bearing Python assertions.
"""

from fractions import Fraction
from itertools import combinations, combinations_with_replacement


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


N_VARIABLES = 4
ZERO_EXPONENT = (0,) * N_VARIABLES


def poly_constant(value):
    return {} if value == 0 else {ZERO_EXPONENT: value}


def poly_variable(index):
    exponent = [0] * N_VARIABLES
    exponent[index] = 1
    return {tuple(exponent): 1}


def poly_add(left, right):
    out = dict(left)
    for monomial, coefficient in right.items():
        out[monomial] = out.get(monomial, 0) + coefficient
        if out[monomial] == 0:
            del out[monomial]
    return out


def poly_scale(poly, scalar):
    return {
        monomial: scalar * coefficient
        for monomial, coefficient in poly.items()
        if scalar * coefficient != 0
    }


def poly_subtract(left, right):
    return poly_add(left, poly_scale(right, -1))


def poly_multiply(left, right):
    out = {}
    for monomial_left, coefficient_left in left.items():
        for monomial_right, coefficient_right in right.items():
            monomial = tuple(
                a + b for a, b in zip(monomial_left, monomial_right)
            )
            out[monomial] = (
                out.get(monomial, 0)
                + coefficient_left * coefficient_right
            )
            if out[monomial] == 0:
                del out[monomial]
    return out


def poly_power(poly, exponent):
    out = poly_constant(1)
    base = dict(poly)
    while exponent:
        if exponent % 2:
            out = poly_multiply(out, base)
        base = poly_multiply(base, base)
        exponent //= 2
    return out


def poly_sum(polys):
    out = {}
    for poly in polys:
        out = poly_add(out, poly)
    return out


ROOTS = tuple(poly_variable(index) for index in range(4))
EDGES = tuple(combinations(range(4), 2))
MATCHINGS = (
    ((0, 1), (2, 3)),
    ((0, 2), (1, 3)),
    ((0, 3), (1, 2)),
)


def root_linear(coefficients):
    return poly_sum(
        poly_scale(ROOTS[index], coefficient)
        for index, coefficient in enumerate(coefficients)
    )


def elementary_polynomials():
    e1 = poly_sum(ROOTS)
    e2 = poly_sum(
        poly_multiply(ROOTS[i], ROOTS[j])
        for i, j in combinations(range(4), 2)
    )
    e3 = poly_sum(
        poly_multiply(poly_multiply(ROOTS[i], ROOTS[j]), ROOTS[k])
        for i, j, k in combinations(range(4), 3)
    )
    return e1, e2, e3


def integer_elementaries(values):
    e1 = sum(values)
    e2 = sum(values[i] * values[j]
             for i, j in combinations(range(4), 2))
    e3 = sum(values[i] * values[j] * values[k]
             for i, j, k in combinations(range(4), 3))
    e4 = values[0] * values[1] * values[2] * values[3]
    return e1, e2, e3, e4


def discriminant_from_roots(values):
    out = 1
    for i, j in combinations(range(len(values)), 2):
        out *= (values[i] - values[j]) ** 2
    return out


def vandermonde(values):
    out = 1
    for i, j in combinations(range(len(values)), 2):
        out *= values[i] - values[j]
    return out


def pair_sums(values):
    return tuple(values[i] + values[j] for i, j in EDGES)


def resolvent_roots(values):
    return tuple(
        values[a] * values[b] + values[c] * values[d]
        for (a, b), (c, d) in MATCHINGS
    )


def opposite_differences(values):
    return tuple(
        values[a] + values[b] - values[c] - values[d]
        for (a, b), (c, d) in MATCHINGS
    )


def main():
    # The fifteen edge pairs split as twelve adjacent pairs and three opposite
    # pairs.  Each of the six sheet differences occurs twice among the twelve.
    adjacent_counts = {pair: 0 for pair in EDGES}
    opposite_pairs = []
    for first_index, second_index in combinations(range(6), 2):
        first = set(EDGES[first_index])
        second = set(EDGES[second_index])
        intersection = first & second
        if intersection:
            require(len(intersection) == 1,
                    "two distinct K4 edges acquired a larger intersection")
            left = next(iter(first - intersection))
            right = next(iter(second - intersection))
            adjacent_counts[tuple(sorted((left, right)))] += 1
        else:
            opposite_pairs.append((EDGES[first_index], EDGES[second_index]))
    require(set(adjacent_counts.values()) == {2},
            "the twelve adjacent factors stopped doubling the Vandermonde")
    require(len(opposite_pairs) == 3
            and {frozenset(pair) for pair in opposite_pairs}
            == {frozenset(pair) for pair in MATCHINGS},
            "the three opposite factors stopped being the perfect matchings")

    # Exact polynomial identities in Z[r1,r2,r3,r4].
    e1, e2, e3 = elementary_polynomials()
    differences = (
        root_linear((1, 1, -1, -1)),
        root_linear((1, -1, 1, -1)),
        root_linear((1, -1, -1, 1)),
    )
    opposite_product = poly_multiply(
        poly_multiply(differences[0], differences[1]), differences[2]
    )
    symmetric_t = poly_add(
        poly_subtract(poly_power(e1, 3),
                      poly_scale(poly_multiply(e1, e2), 4)),
        poly_scale(e3, 8),
    )
    require(opposite_product == symmetric_t,
            "T stopped being e1^3-4e1e2+8e3")

    z_polys = tuple(
        poly_add(poly_multiply(ROOTS[a], ROOTS[b]),
                 poly_multiply(ROOTS[c], ROOTS[d]))
        for (a, b), (c, d) in MATCHINGS
    )
    for difference, z_poly in zip(differences, z_polys):
        require(
            poly_power(difference, 2)
            == poly_add(
                poly_subtract(poly_power(e1, 2), poly_scale(e2, 4)),
                poly_scale(z_poly, 4),
            ),
            "a centered pair-square stopped being a resolvent-root translate",
        )

    resolvent_difference_factors = (
        ((0, 3), (1, 2)),
        ((0, 2), (1, 3)),
        ((0, 1), (2, 3)),
    )
    for (i, j), factor_pairs in zip(
        combinations(range(3), 2), resolvent_difference_factors
    ):
        first_pair, second_pair = factor_pairs
        expected = poly_multiply(
            poly_subtract(ROOTS[first_pair[0]], ROOTS[first_pair[1]]),
            poly_subtract(ROOTS[second_pair[0]], ROOTS[second_pair[1]]),
        )
        require(poly_subtract(z_polys[i], z_polys[j]) == expected,
                "a resolvent difference stopped factoring into root differences")

    # Exact integer census, including repeated-root tuples.
    distinct_total = 0
    distinct_opposite_walls = 0
    for values in combinations(range(-4, 6), 4):
        distinct_total += 1
        e1_value, e2_value, e3_value, _ = integer_elementaries(values)
        t_value = e1_value ** 3 - 4 * e1_value * e2_value + 8 * e3_value
        quartic_discriminant = discriminant_from_roots(values)
        sextic_vandermonde = vandermonde(pair_sums(values))
        sextic_discriminant = discriminant_from_roots(pair_sums(values))
        require(sextic_vandermonde == quartic_discriminant * t_value,
                "the lexicographic Vandermonde sign convention changed")
        require(
            sextic_discriminant
            == quartic_discriminant ** 2 * t_value ** 2,
            "the pair-sum discriminant square failed on a distinct-root row",
        )
        z_values = resolvent_roots(values)
        require(discriminant_from_roots(z_values) == quartic_discriminant,
                "the resolvent discriminant stopped equalling the quartic one")
        shift = Fraction(e1_value ** 2, 4) - e2_value
        d_values = opposite_differences(values)
        require(
            tuple(Fraction(d * d, 4) for d in d_values)
            == tuple(Fraction(z, 1) + shift for z in z_values),
            "the centered quadratic-pullback roots changed",
        )
        if t_value == 0:
            distinct_opposite_walls += 1

    repeated_total = 0
    for values in combinations_with_replacement(range(-2, 3), 4):
        if len(set(values)) == 4:
            continue
        repeated_total += 1
        e1_value, e2_value, e3_value, _ = integer_elementaries(values)
        t_value = e1_value ** 3 - 4 * e1_value * e2_value + 8 * e3_value
        quartic_discriminant = discriminant_from_roots(values)
        sextic_discriminant = discriminant_from_roots(pair_sums(values))
        require(quartic_discriminant == 0 and sextic_discriminant == 0,
                "a repeated quartic root stopped forcing pair-sum collision")
        require(sextic_discriminant
                == quartic_discriminant ** 2 * t_value ** 2,
                "the polynomial identity failed on a repeated-root row")

    hostile = (-3, -1, 1, 3)
    hostile_e1, hostile_e2, hostile_e3, hostile_e4 = integer_elementaries(
        hostile
    )
    hostile_t = (
        hostile_e1 ** 3
        - 4 * hostile_e1 * hostile_e2
        + 8 * hostile_e3
    )
    hostile_quartic_discriminant = discriminant_from_roots(hostile)
    hostile_pair_sums = pair_sums(hostile)
    require(hostile_quartic_discriminant != 0
            and hostile_t == 0
            and discriminant_from_roots(hostile_pair_sums) == 0
            and hostile_pair_sums.count(0) == 2,
            "the separable-quartic/multiple-sextic hostile changed")

    generic = (0, 1, 3, 7)
    generic_e1, generic_e2, generic_e3, _ = integer_elementaries(generic)
    generic_t = generic_e1 ** 3 - 4 * generic_e1 * generic_e2 + 8 * generic_e3
    generic_quartic_discriminant = discriminant_from_roots(generic)
    generic_sextic_discriminant = discriminant_from_roots(pair_sums(generic))

    print("QUARTIC PAIR-SUM SEXTIC / RESOLVENT PULLBACK AUDIT")
    print("edge_pair_partition=12_adjacent+3_opposite")
    print("adjacent_factor_counts=each_of_6_root_differences_twice")
    print("T=product_opposite_differences=e1^3-4*e1*e2+8*e3")
    print("centered_roots=+-d_m/2 t_m=d_m^2/4")
    print("t_m=z_m+e1^2/4-e2 standard_resolvent_translate=exact")
    print("disc(resolvent)=disc(quartic)")
    print("lex_edge_vandermonde=+disc(quartic)*T")
    print("disc(pair_sum_sextic)=disc(quartic)^2*T^2")
    print(
        f"distinct_census={distinct_total} opposite_sum_walls="
        f"{distinct_opposite_walls} repeated_controls={repeated_total}"
    )
    print(
        "generic_roots=0,1,3,7 "
        f"disc4={generic_quartic_discriminant} T={generic_t} "
        f"disc6={generic_sextic_discriminant}"
    )
    print(
        "hostile_roots=-3,-1,1,3 quartic="
        f"x^4{hostile_e2:+d}x^2{(-hostile_e3):+d}x{hostile_e4:+d} "
        f"disc4={hostile_quartic_discriminant} T={hostile_t}"
    )
    print(f"hostile_pair_sums={','.join(map(str, hostile_pair_sums))} disc6=0")
    print("coefficient_wall_for_x4+a*x3+b*x2+c*x+d=a^3-4*a*b+8*c=0")
    print("SCOPE: exact quartic/resolvent polynomial identity; no Keller exclusion")
    print("FAILED CHECKS: NONE")


if __name__ == "__main__":
    main()
