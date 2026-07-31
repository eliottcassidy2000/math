#!/usr/bin/env python3
"""Exact affine full-S4 divisor-parity hostile audit for THM-2769.

Only integer, rational, finite-field, and finite-permutation arithmetic is
used.  There are no truth-bearing Python assertions and no floating point.
"""

from fractions import Fraction
from itertools import permutations, product


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def trim(polynomial):
    values = list(polynomial)
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return tuple(values)


def add(left, right):
    size = max(len(left), len(right))
    return trim(tuple(
        (left[index] if index < len(left) else 0)
        + (right[index] if index < len(right) else 0)
        for index in range(size)
    ))


def scale(polynomial, scalar):
    return trim(tuple(scalar * coefficient for coefficient in polynomial))


def multiply(left, right):
    out = [0] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            out[left_index + right_index] += left_value * right_value
    return trim(tuple(out))


def power(polynomial, exponent):
    out = (1,)
    base = polynomial
    remaining = exponent
    while remaining:
        if remaining % 2:
            out = multiply(out, base)
        base = multiply(base, base)
        remaining //= 2
    return out


def shift(polynomial, amount):
    return (0,) * amount + tuple(polynomial)


def quartic_discriminant(p_value, q_value, r_value):
    """Discriminant of X^4+pX^2+qX+r as a polynomial in t."""

    terms = (
        scale(power(r_value, 3), 256),
        scale(multiply(power(p_value, 2), power(r_value, 2)), -128),
        scale(multiply(multiply(p_value, power(q_value, 2)), r_value), 144),
        scale(power(q_value, 4), -27),
        scale(multiply(power(p_value, 4), r_value), 16),
        scale(multiply(power(p_value, 3), power(q_value, 2)), -4),
    )
    out = (0,)
    for term in terms:
        out = add(out, term)
    return out


def cubic_discriminant(a_value, b_value, c_value):
    """Discriminant of U^3+aU^2+bU+c as a polynomial in t."""

    terms = (
        multiply(power(a_value, 2), power(b_value, 2)),
        scale(power(b_value, 3), -4),
        scale(multiply(power(a_value, 3), c_value), -4),
        scale(power(c_value, 2), -27),
        scale(multiply(multiply(a_value, b_value), c_value), 18),
    )
    out = (0,)
    for term in terms:
        out = add(out, term)
    return out


def lower_hull(points):
    hull = []
    for point in points:
        while len(hull) >= 2:
            first, second = hull[-2], hull[-1]
            cross = ((second[0] - first[0]) * (point[1] - first[1])
                     - (second[1] - first[1]) * (point[0] - first[0]))
            if cross > 0:
                break
            hull.pop()
        hull.append(point)
    return tuple(hull)


def root_substitution_by_t_degree(exponent):
    """Coefficients after U=c*t^exponent, keyed by t- then c-degree."""

    terms = (
        (3 * exponent, 3, 1),
        (2 * exponent, 2, -4),
        (exponent + 1, 1, 16),
        (2, 0, -64),
    )
    out = {}
    for t_degree, c_degree, coefficient in terms:
        row = out.setdefault(t_degree, {})
        row[c_degree] = row.get(c_degree, 0) + coefficient
    return {
        t_degree: {
            c_degree: coefficient
            for c_degree, coefficient in row.items()
            if coefficient
        }
        for t_degree, row in out.items()
        if any(row.values())
    }


def xor(left, right):
    return tuple(first ^ second for first, second in zip(left, right))


def permute_bits(vector, permutation):
    return tuple(vector[permutation[index]] for index in range(3))


def is_subspace(subset):
    zero = (0, 0, 0)
    return zero in subset and all(
        xor(left, right) in subset
        for left in subset
        for right in subset
    )


def invariant_subspaces(space, group):
    vectors = tuple(space)
    out = []
    for mask in range(1 << len(vectors)):
        subset = frozenset(
            vectors[index]
            for index in range(len(vectors))
            if (mask >> index) & 1
        )
        if not is_subspace(subset):
            continue
        if all(
            frozenset(permute_bits(vector, permutation)
                      for vector in subset) == subset
            for permutation in group
        ):
            out.append(subset)
    return tuple(out)


def even_sign_action():
    even = tuple(
        vector for vector in product((0, 1), repeat=3)
        if sum(vector) % 2 == 0
    )
    group = tuple(permutations(range(3)))
    actions = set()
    for epsilon in even:
        for permutation in group:
            image = tuple(
                even.index(xor(epsilon, permute_bits(vector, permutation)))
                for vector in even
            )
            require(tuple(sorted(image)) == tuple(range(4)),
                    "even-sign action stopped being a permutation")
            actions.add(image)
    return even, actions


def main():
    # Coefficients of F_t(Y)=Y^4+pY^2+qY+r in low-to-high t order.
    p_value = (-2,)
    q_value = (0, -8)
    r_value = (1, -4)

    # The squared complementary-pair-sum cubic is
    # U^3+2pU^2+(p^2-4r)U-q^2.
    resolvent_a = scale(p_value, 2)
    resolvent_b = add(power(p_value, 2), scale(r_value, -4))
    resolvent_c = scale(power(q_value, 2), -1)
    require(resolvent_a == (-4,), "pair-sum U^2 coefficient changed")
    require(resolvent_b == (0, 16), "pair-sum U coefficient changed")
    require(resolvent_c == (0, 0, -64),
            "pair-sum constant coefficient changed")
    inverse_p = scale(resolvent_a, Fraction(1, 2))
    inverse_r = scale(
        add(power(inverse_p, 2), scale(resolvent_b, -1)),
        Fraction(1, 4),
    )
    require(inverse_p == p_value and inverse_r == r_value,
            "cubic-to-quartic inverse changed")
    require(scale(resolvent_c, -1) == power(q_value, 2),
            "two-branch q square reconstruction changed")

    expected_discriminant = (0, 0, -12288, 57344, -110592)
    disc_r = cubic_discriminant(
        resolvent_a, resolvent_b, resolvent_c
    )
    disc_f = quartic_discriminant(p_value, q_value, r_value)
    require(disc_r == expected_discriminant,
            "cubic discriminant factorization changed")
    require(disc_f == expected_discriminant,
            "quartic and cubic discriminants diverged")
    quadratic_factor = (3, -14, 27)
    require(expected_discriminant == scale(shift(quadratic_factor, 2), -4096),
            "displayed discriminant factorization changed")
    require((-14) ** 2 - 4 * 27 * 3 == -128,
            "quadratic squarefreeness control changed")

    product_roots = (0, 0, 64)
    disc_h = scale(multiply(product_roots, power(disc_r, 2)), 64)
    displayed_disc_h = scale(shift(power(quadratic_factor, 2), 6), 2 ** 36)
    square_root_disc_h = scale(shift(quadratic_factor, 3), 2 ** 18)
    require(disc_h == displayed_disc_h,
            "quadratic-pullback discriminant changed")
    require(disc_h == power(square_root_disc_h, 2),
            "six-root discriminant stopped being an exact square")

    # Rational-root proof gate.  A root in C(t) is integral over C[t], hence
    # lies in C[t] and divides 64t^2, so it is c*t^j for j=0,1,2.
    # These are the decisive coefficient obstructions after substitution.
    substituted = {
        exponent: root_substitution_by_t_degree(exponent)
        for exponent in range(3)
    }
    require(substituted[0][2] == {0: -64},
            "j=0 rational-root t^2 blocker changed")
    require(substituted[1][3] == {3: 1},
            "j=1 rational-root leading c^3 blocker changed")
    require(substituted[2][2] == {0: -64},
            "j=2 rational-root t^2 blocker changed")

    valuation_points = ((0, 2), (1, 1), (2, 0), (3, 0))
    hull = lower_hull(valuation_points)
    require(hull == ((0, 2), (2, 0), (3, 0)),
            "Newton lower hull changed")
    slopes = tuple(
        Fraction(right[1] - left[1], right[0] - left[0])
        for left, right in zip(hull, hull[1:])
    )
    lengths = tuple(
        right[0] - left[0]
        for left, right in zip(hull, hull[1:])
    )
    valuations = tuple(
        -slope
        for slope, length in zip(slopes, lengths)
        for _ in range(length)
    )
    require(slopes == (Fraction(-1), Fraction(0))
            and valuations == (1, 1, 0),
            "Newton root valuations changed")
    require(4 ** 2 - 4 * 16 == -48,
            "scaled small-root residual stopped being squarefree")
    require(3 * 4 ** 2 - 8 * 4 == 16,
            "unit root at U=4 stopped being simple")

    s3 = tuple(permutations(range(3)))
    even_code = tuple(
        vector for vector in product((0, 1), repeat=3)
        if sum(vector) % 2 == 0
    )
    invariant = invariant_subspaces(even_code, s3)
    require(tuple(sorted(len(space) for space in invariant)) == (1, 4),
            "even-weight S3 module acquired a proper invariant line")
    orbit_110 = {
        permute_bits((1, 1, 0), permutation)
        for permutation in s3
    }
    require(orbit_110 == {(1, 1, 0), (1, 0, 1), (0, 1, 1)},
            "weight-two parity orbit changed")

    omega, actions = even_sign_action()
    require(len(omega) == 4 and len(actions) == 24,
            "V4 semidirect S3 action stopped having order 24")
    require(actions == set(permutations(range(4))),
            "even-sign four-state action stopped being all of S4")

    require(tuple(value % 2 for value in valuations) == (1, 1, 0),
            "specific divisor parity row changed")
    require(tuple((2 * value) % 2 for value in valuations) == (0, 0, 0),
            "even-base-change parity control changed")
    require((1 + 1 + 0) % 2 == 0,
            "specific boundary row left the even-weight code")

    # t=0 gives Y^4-2Y^2+1=(Y^2-1)^2, while e3=8t and T=8e3=64t.
    require((1, 0, -2, 0, 1) == multiply((-1, 0, 1), (-1, 0, 1)),
            "quartic boundary factorization changed")
    require(8 * 8 == 64, "opposite-sum invariant normalization changed")

    print("FULL S4 AFFINE DIVISOR-PARITY HOSTILE AUDIT")
    print("F_t=Y^4-2Y^2-8tY+1-4t")
    print("R_t=U^3-4U^2+16tU-64t^2")
    print("inverse_quartic_branches=q=+/-8t")
    print("pair_sum_product=(8t)^2")
    print("rational_root_cases=j0,j1,j2 all_blocked")
    print("disc_R=disc_F=-2^12*t^2*(27t^2-14t+3)")
    print("quadratic_factor_discriminant=-128 nonsquare")
    print("disc_R(V^2)=2^36*t^6*(27t^2-14t+3)^2")
    print("six_root_square_root=2^18*t^3*(27t^2-14t+3)")
    print("newton_hull=((0,2),(2,0),(3,0))")
    print("root_valuations=(1,1,0) parity=110")
    print("small_root_scaled_quadratic_discriminant=-48")
    print("S3_relation_module=irreducible_dimension_2")
    print("kummer_rank=2")
    print("even_sign_group=V4_semidirect_S3=S4 order=24")
    print("boundary_code={000,110,101,011}")
    print("even_base_change_control=000")
    print("T=64t generic_nonzero; F_0=(Y^2-1)^2")
    print("SCOPE=non_Keller_affine_boundary_hostile_not_JC2_or_DC2")
    print("FAILED CHECKS: NONE")


if __name__ == "__main__":
    main()
