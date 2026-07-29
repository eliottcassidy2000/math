#!/usr/bin/env python3
"""Exact quartic edge/orientation sextic-resolvent audit for THM-2864.

The companion is dependency free.  It implements its own sparse integer
polynomial arithmetic, with every verification routed through explicit
exceptions that remain active in optimized mode.
"""

from collections import Counter
from itertools import combinations, permutations


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


class SparsePolynomial:
    """Sparse polynomial over Z in a fixed number of variables."""

    __slots__ = ("variable_count", "terms")

    def __init__(self, variable_count, terms=None):
        self.variable_count = variable_count
        cleaned = {}
        for monomial, coefficient in (terms or {}).items():
            require(
                len(monomial) == variable_count,
                "polynomial monomial has wrong arity",
            )
            if coefficient:
                cleaned[tuple(monomial)] = (
                    cleaned.get(tuple(monomial), 0) + int(coefficient)
                )
        self.terms = {
            monomial: coefficient
            for monomial, coefficient in cleaned.items()
            if coefficient
        }

    @classmethod
    def constant(cls, variable_count, coefficient):
        if coefficient == 0:
            return cls(variable_count)
        return cls(
            variable_count,
            {(0,) * variable_count: int(coefficient)},
        )

    @classmethod
    def variable(cls, variable_count, position):
        exponent = [0] * variable_count
        exponent[position] = 1
        return cls(variable_count, {tuple(exponent): 1})

    def coerce(self, other):
        if isinstance(other, SparsePolynomial):
            require(
                other.variable_count == self.variable_count,
                "mixed polynomial arities",
            )
            return other
        return SparsePolynomial.constant(self.variable_count, other)

    def __add__(self, other):
        other = self.coerce(other)
        terms = dict(self.terms)
        for monomial, coefficient in other.terms.items():
            terms[monomial] = terms.get(monomial, 0) + coefficient
        return SparsePolynomial(self.variable_count, terms)

    def __radd__(self, other):
        return self + other

    def __neg__(self):
        return SparsePolynomial(
            self.variable_count,
            {
                monomial: -coefficient
                for monomial, coefficient in self.terms.items()
            },
        )

    def __sub__(self, other):
        return self + (-self.coerce(other))

    def __rsub__(self, other):
        return self.coerce(other) - self

    def __mul__(self, other):
        other = self.coerce(other)
        terms = {}
        for left_monomial, left_coefficient in self.terms.items():
            for right_monomial, right_coefficient in other.terms.items():
                monomial = tuple(
                    left + right
                    for left, right
                    in zip(left_monomial, right_monomial)
                )
                terms[monomial] = (
                    terms.get(monomial, 0)
                    + left_coefficient * right_coefficient
                )
        return SparsePolynomial(self.variable_count, terms)

    def __rmul__(self, other):
        return self * other

    def __pow__(self, exponent):
        require(
            isinstance(exponent, int) and exponent >= 0,
            "polynomial exponent must be a nonnegative integer",
        )
        result = SparsePolynomial.constant(self.variable_count, 1)
        base = self
        remaining = exponent
        while remaining:
            if remaining & 1:
                result = result * base
            base = base * base
            remaining //= 2
        return result

    def __eq__(self, other):
        other = self.coerce(other)
        return self.terms == other.terms

    def is_zero(self):
        return not self.terms

    def term_count(self):
        return len(self.terms)


def quartic_discriminant(p_value, q_value, r_value):
    return (
        256 * r_value**3
        - 128 * p_value**2 * r_value**2
        + 144 * p_value * q_value**2 * r_value
        - 27 * q_value**4
        + 16 * p_value**4 * r_value
        - 4 * p_value**3 * q_value**2
    )


def orientation_wall(p_value, q_value, r_value):
    return (
        1024 * r_value**3
        + 768 * p_value**2 * r_value**2
        - 288 * p_value * q_value**2 * r_value
        + 27 * q_value**4
    )


def matching_resolvent(u_value, p_value, q_value, r_value):
    return (
        u_value**3
        + 2 * p_value * u_value**2
        + (p_value**2 - 4 * r_value) * u_value
        - q_value**2
    )


def orientation_radicand(u_value, p_value, r_value):
    return (
        16
        * (
            16 * r_value
            - 4 * p_value * u_value
            - 3 * u_value**2
        )
        * (
            u_value**2
            + 2 * p_value * u_value
            + p_value**2
            - 4 * r_value
        )
    )


def matching_derivative(u_value, p_value, r_value):
    return (
        3 * u_value**2
        + 4 * p_value * u_value
        + p_value**2
        - 4 * r_value
    )


def monic_cubic_discriminant(
    coefficient_two,
    coefficient_one,
    coefficient_zero,
):
    return (
        coefficient_two**2 * coefficient_one**2
        - 4 * coefficient_one**3
        - 4 * coefficient_two**3 * coefficient_zero
        - 27 * coefficient_zero**2
        + 18
        * coefficient_two
        * coefficient_one
        * coefficient_zero
    )


def orientation_coefficients(p_value, q_value, r_value):
    discriminant = quartic_discriminant(
        p_value,
        q_value,
        r_value,
    )
    coefficient_four = -32 * (
        8 * p_value**2 * r_value
        - 3 * p_value * q_value**2
        - 32 * r_value**2
    )
    coefficient_two = (
        -256
        * q_value**2
        * (p_value**2 + 12 * r_value)
        * (32 * p_value * r_value - 9 * q_value**2)
    )
    coefficient_zero = (
        -4096 * q_value**4 * discriminant
    )
    return (
        coefficient_zero,
        0,
        coefficient_two,
        0,
        coefficient_four,
        0,
        1,
    )


def edge_coefficients(p_value, q_value, r_value):
    return (
        -q_value**2,
        0,
        p_value**2 - 4 * r_value,
        0,
        2 * p_value,
        0,
        1,
    )


def polynomial_value(coefficients, value):
    result = 0
    for coefficient in reversed(coefficients):
        result = result * value + coefficient
    return result


def polynomial_from_roots(roots):
    coefficients = [1]
    for root in roots:
        product_coefficients = [0] * (len(coefficients) + 1)
        for position, coefficient in enumerate(coefficients):
            product_coefficients[position] -= coefficient * root
            product_coefficients[position + 1] += coefficient
        coefficients = product_coefficients
    return tuple(coefficients)


def coefficient_product(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for left_position, left_coefficient in enumerate(left):
        for right_position, right_coefficient in enumerate(right):
            result[left_position + right_position] += (
                left_coefficient * right_coefficient
            )
    return tuple(result)


def coefficient_derivative(coefficients):
    if len(coefficients) <= 1:
        return (0,)
    return tuple(
        degree * coefficients[degree]
        for degree in range(1, len(coefficients))
    )


def monic_division_remainder(dividend, divisor):
    require(
        divisor[-1] == 1,
        "exact division helper requires a monic divisor",
    )
    remainder = list(dividend)
    while len(remainder) >= len(divisor):
        leading = remainder[-1]
        shift = len(remainder) - len(divisor)
        for position, coefficient in enumerate(divisor):
            remainder[shift + position] -= leading * coefficient
        while remainder and remainder[-1] == 0:
            remainder.pop()
    return tuple(remainder or [0])


def quotient_zero(variable_count):
    zero = SparsePolynomial.constant(variable_count, 0)
    return (zero, zero, zero)


def quotient_constant(polynomial):
    zero = SparsePolynomial.constant(polynomial.variable_count, 0)
    return (polynomial, zero, zero)


def quotient_add(left, right):
    return tuple(
        left_coordinate + right_coordinate
        for left_coordinate, right_coordinate in zip(left, right)
    )


def quotient_negate(value):
    return tuple(-coordinate for coordinate in value)


def quotient_subtract(left, right):
    return quotient_add(left, quotient_negate(right))


def quotient_scale(value, scalar):
    return tuple(coordinate * scalar for coordinate in value)


def quotient_multiply(left, right, p_value, q_value, r_value):
    variable_count = p_value.variable_count
    zero = SparsePolynomial.constant(variable_count, 0)
    coefficients = [zero for _ in range(5)]
    for left_degree, left_coefficient in enumerate(left):
        for right_degree, right_coefficient in enumerate(right):
            coefficients[left_degree + right_degree] = (
                coefficients[left_degree + right_degree]
                + left_coefficient * right_coefficient
            )
    relation = (
        q_value**2,
        -(p_value**2 - 4 * r_value),
        -2 * p_value,
    )
    for degree in range(4, 2, -1):
        leading = coefficients[degree]
        coefficients[degree] = zero
        shift = degree - 3
        for relation_degree, relation_coefficient in enumerate(
            relation
        ):
            coefficients[shift + relation_degree] = (
                coefficients[shift + relation_degree]
                + leading * relation_coefficient
            )
    return tuple(coefficients[:3])


def determinant_three(matrix):
    return (
        matrix[0][0]
        * (
            matrix[1][1] * matrix[2][2]
            - matrix[1][2] * matrix[2][1]
        )
        - matrix[0][1]
        * (
            matrix[1][0] * matrix[2][2]
            - matrix[1][2] * matrix[2][0]
        )
        + matrix[0][2]
        * (
            matrix[1][0] * matrix[2][1]
            - matrix[1][1] * matrix[2][0]
        )
    )


def permutation_compose(left, right):
    return tuple(left[right[position]] for position in range(4))


def permutation_inverse(permutation):
    inverse = [0] * 4
    for position, image in enumerate(permutation):
        inverse[image] = position
    return tuple(inverse)


def permutation_order(permutation):
    identity = tuple(range(4))
    power = identity
    for exponent in range(1, 13):
        power = permutation_compose(power, permutation)
        if power == identity:
            return exponent
    raise RuntimeError("permutation order exceeded exact cap")


def permutation_cycle_type(permutation):
    visited = set()
    lengths = []
    for start in range(4):
        if start in visited:
            continue
        current = start
        length = 0
        while current not in visited:
            visited.add(current)
            length += 1
            current = permutation[current]
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


def generated_permutation_group(generators):
    identity = tuple(range(4))
    group = {identity}
    pending = [identity]
    while pending:
        element = pending.pop()
        for generator in generators:
            candidate = permutation_compose(element, generator)
            if candidate not in group:
                group.add(candidate)
                pending.append(candidate)
    return frozenset(group)


def require_permutation_subgroup(subgroup, label):
    identity = tuple(range(4))
    require(identity in subgroup, f"{label}: identity missing")
    for left in subgroup:
        require(
            permutation_inverse(left) in subgroup,
            f"{label}: inverse missing",
        )
        for right in subgroup:
            require(
                permutation_compose(left, right) in subgroup,
                f"{label}: closure failed",
            )


def quartic_coefficients_from_roots(roots):
    p_value = sum(
        roots[left] * roots[right]
        for left, right in combinations(range(4), 2)
    )
    q_value = -sum(
        roots[first] * roots[second] * roots[third]
        for first, second, third in combinations(range(4), 3)
    )
    r_value = 1
    for root in roots:
        r_value *= root
    return p_value, q_value, r_value


def vandermonde_square(roots):
    result = 1
    for left, right in combinations(range(4), 2):
        result *= (roots[left] - roots[right]) ** 2
    return result


def edge_value(roots, permutation):
    return (
        roots[permutation[0]]
        + roots[permutation[2]]
    )


def matching_value(roots, permutation):
    return edge_value(roots, permutation) ** 2


def omega_value(roots, permutation):
    first, second, third, fourth = (
        roots[position]
        for position in permutation
    )
    left_difference = first - third
    right_difference = second - fourth
    return (
        left_difference
        * right_difference
        * (
            left_difference**2
            - right_difference**2
        )
    )


def orbit_representatives(group, value_function):
    values = {}
    for permutation in group:
        values.setdefault(
            value_function(permutation),
            permutation,
        )
    return tuple(
        values[value]
        for value in sorted(values)
    )


def character_with_kernel(group, kernel):
    return {
        element: 1 if element in kernel else -1
        for element in group
    }


def check_character(group, character, label):
    for left in group:
        for right in group:
            require(
                character[permutation_compose(left, right)]
                == character[left] * character[right],
                f"{label}: character law failed",
            )


def check_specialization(
    roots,
    edge_representatives,
    orientation_representatives,
):
    p_value, q_value, r_value = quartic_coefficients_from_roots(roots)
    discriminant = quartic_discriminant(
        p_value,
        q_value,
        r_value,
    )
    require(
        discriminant == vandermonde_square(roots),
        "specialization discriminant failed",
    )
    edge_roots = tuple(
        edge_value(roots, representative)
        for representative in edge_representatives
    )
    orientation_roots = tuple(
        omega_value(roots, representative)
        for representative in orientation_representatives
    )
    require(
        polynomial_from_roots(edge_roots)
        == edge_coefficients(p_value, q_value, r_value),
        "specialized edge polynomial failed",
    )
    require(
        polynomial_from_roots(orientation_roots)
        == orientation_coefficients(
            p_value,
            q_value,
            r_value,
        ),
        "specialized orientation polynomial failed",
    )
    return {
        "p": p_value,
        "q": q_value,
        "r": r_value,
        "discriminant": discriminant,
        "edge_distinct": len(set(edge_roots)),
        "orientation_distinct": len(set(orientation_roots)),
    }


def main():
    # Root-coordinate symbolic proof in Z[a,b,c], with d=-a-b-c.
    root_variable_count = 3
    root_a = SparsePolynomial.variable(root_variable_count, 0)
    root_b = SparsePolynomial.variable(root_variable_count, 1)
    root_c = SparsePolynomial.variable(root_variable_count, 2)
    root_d = -(root_a + root_b + root_c)
    roots_symbolic = (root_a, root_b, root_c, root_d)
    p_root = sum(
        roots_symbolic[left] * roots_symbolic[right]
        for left, right in combinations(range(4), 2)
    )
    q_root = -sum(
        roots_symbolic[first]
        * roots_symbolic[second]
        * roots_symbolic[third]
        for first, second, third in combinations(range(4), 3)
    )
    r_root = (
        root_a * root_b * root_c * root_d
    )
    u_root = (root_a + root_c) ** 2
    omega_root = (
        (root_a - root_c)
        * (root_b - root_d)
        * (
            (root_a - root_c) ** 2
            - (root_b - root_d) ** 2
        )
    )
    require(
        matching_resolvent(
            u_root,
            p_root,
            q_root,
            r_root,
        ).is_zero(),
        "symbolic matching-root identity failed",
    )
    f_root = orientation_radicand(
        u_root,
        p_root,
        r_root,
    )
    require(
        omega_root**2 == f_root,
        "symbolic Omega-square identity failed",
    )
    delta_root = quartic_discriminant(
        p_root,
        q_root,
        r_root,
    )
    vandermonde_root = SparsePolynomial.constant(
        root_variable_count,
        1,
    )
    for left, right in combinations(range(4), 2):
        vandermonde_root = (
            vandermonde_root
            * (
                roots_symbolic[left]
                - roots_symbolic[right]
            ) ** 2
        )
    require(
        delta_root == vandermonde_root,
        "symbolic quartic discriminant failed",
    )

    # Coefficient-coordinate proof in Z[p,q,r,Z].
    coefficient_variable_count = 4
    p_value = SparsePolynomial.variable(
        coefficient_variable_count,
        0,
    )
    q_value = SparsePolynomial.variable(
        coefficient_variable_count,
        1,
    )
    r_value = SparsePolynomial.variable(
        coefficient_variable_count,
        2,
    )
    z_value = SparsePolynomial.variable(
        coefficient_variable_count,
        3,
    )
    zero = SparsePolynomial.constant(
        coefficient_variable_count,
        0,
    )
    one = SparsePolynomial.constant(
        coefficient_variable_count,
        1,
    )
    quotient_one = quotient_constant(one)
    quotient_u = (zero, one, zero)
    quotient_u_squared = quotient_multiply(
        quotient_u,
        quotient_u,
        p_value,
        q_value,
        r_value,
    )
    first_factor = quotient_add(
        quotient_constant(16 * r_value),
        quotient_add(
            quotient_scale(
                quotient_u,
                -4 * p_value,
            ),
            quotient_scale(
                quotient_u_squared,
                -3,
            ),
        ),
    )
    second_factor = quotient_add(
        quotient_u_squared,
        quotient_add(
            quotient_scale(
                quotient_u,
                2 * p_value,
            ),
            quotient_constant(
                p_value**2 - 4 * r_value,
            ),
        ),
    )
    quotient_f = quotient_scale(
        quotient_multiply(
            first_factor,
            second_factor,
            p_value,
            q_value,
            r_value,
        ),
        16,
    )
    quotient_c = quotient_add(
        quotient_scale(
            quotient_u_squared,
            3,
        ),
        quotient_add(
            quotient_scale(
                quotient_u,
                4 * p_value,
            ),
            quotient_constant(
                p_value**2 - 4 * r_value,
            ),
        ),
    )
    basis = (
        quotient_one,
        quotient_u,
        quotient_u_squared,
    )
    multiplication_matrix = [[zero for _ in range(3)] for _ in range(3)]
    for column, basis_vector in enumerate(basis):
        image = quotient_multiply(
            quotient_f,
            basis_vector,
            p_value,
            q_value,
            r_value,
        )
        for row in range(3):
            multiplication_matrix[row][column] = image[row]
    norm_matrix = [
        [
            (
                z_value
                if row == column
                else zero
            )
            - multiplication_matrix[row][column]
            for column in range(3)
        ]
        for row in range(3)
    ]
    norm_polynomial = determinant_three(norm_matrix)
    delta_coefficient = quartic_discriminant(
        p_value,
        q_value,
        r_value,
    )
    matching_discriminant = monic_cubic_discriminant(
        2 * p_value,
        p_value**2 - 4 * r_value,
        -q_value**2,
    )
    require(
        matching_discriminant == delta_coefficient,
        "matching cubic discriminant is not quartic Delta",
    )
    orientation_coefficient_four = -32 * (
        8 * p_value**2 * r_value
        - 3 * p_value * q_value**2
        - 32 * r_value**2
    )
    orientation_coefficient_two = (
        -256
        * q_value**2
        * (p_value**2 + 12 * r_value)
        * (32 * p_value * r_value - 9 * q_value**2)
    )
    orientation_coefficient_zero = (
        -4096 * q_value**4 * delta_coefficient
    )
    expected_norm = (
        z_value**3
        + orientation_coefficient_four * z_value**2
        + orientation_coefficient_two * z_value
        + orientation_coefficient_zero
    )
    require(
        norm_polynomial == expected_norm,
        "coefficient orientation resultant failed",
    )
    product_remainder = quotient_subtract(
        quotient_multiply(
            quotient_multiply(
                quotient_c,
                quotient_c,
                p_value,
                q_value,
                r_value,
            ),
            quotient_multiply(
                quotient_u,
                quotient_f,
                p_value,
                q_value,
                r_value,
            ),
            p_value,
            q_value,
            r_value,
        ),
        quotient_constant(
            16 * q_value**2 * delta_coefficient,
        ),
    )
    require(
        all(coordinate.is_zero() for coordinate in product_remainder),
        "D8 radicand product modulo S failed",
    )
    orientation_collision_invariant = monic_cubic_discriminant(
        orientation_coefficient_four,
        orientation_coefficient_two,
        orientation_coefficient_zero,
    )
    require(
        not orientation_collision_invariant.is_zero(),
        "universal orientation collision invariant vanished",
    )
    orientation_wall_polynomial = orientation_wall(
        p_value,
        q_value,
        r_value,
    )
    require(
        orientation_collision_invariant
        == (
            2**24
            * q_value**4
            * orientation_wall_polynomial**2
            * delta_coefficient
        ),
        "Theta factorization failed",
    )

    # General root-product proof for Disc(P(Y^2)).
    even_root_variable_count = 3
    even_roots = tuple(
        SparsePolynomial.variable(
            even_root_variable_count,
            position,
        )
        for position in range(3)
    )
    even_root_product = (
        even_roots[0] * even_roots[1] * even_roots[2]
    )
    cubic_root_discriminant = (
        (even_roots[0] - even_roots[1]) ** 2
        * (even_roots[0] - even_roots[2]) ** 2
        * (even_roots[1] - even_roots[2]) ** 2
    )
    even_sextic_root_discriminant = (
        4**3
        * even_root_product
        * cubic_root_discriminant**2
    )
    require(
        even_sextic_root_discriminant
        == (
            2**6
            * even_root_product
            * cubic_root_discriminant**2
        ),
        "general even-sextic root discriminant failed",
    )
    sextic_discriminant = (
        -2**6
        * orientation_coefficient_zero
        * orientation_collision_invariant**2
    )
    require(
        sextic_discriminant
        == (
            2**18
            * q_value**4
            * delta_coefficient
            * orientation_collision_invariant**2
        )
        == (
            2**66
            * q_value**12
            * orientation_wall_polynomial**4
            * delta_coefficient**3
        ),
        "orientation sextic discriminant identity failed",
    )

    # A separated exact control realizes the generic S4 orbit/stabilizers.
    generic_roots = (-6, -2, 3, 5)
    require(
        len(set(generic_roots)) == 4
        and sum(generic_roots) == 0,
        "generic root control failed",
    )
    group = tuple(permutations(range(4)))
    p_generic, q_generic, r_generic = (
        quartic_coefficients_from_roots(generic_roots)
    )
    delta_generic = quartic_discriminant(
        p_generic,
        q_generic,
        r_generic,
    )
    require(
        (p_generic, q_generic, r_generic)
        == (-37, 24, 180)
        and delta_generic == 768398400
        and delta_generic == vandermonde_square(generic_roots),
        "generic coefficient control failed",
    )
    generic_orientation_coefficients = orientation_coefficients(
        p_generic,
        q_generic,
        r_generic,
    )
    generic_orientation_collision = monic_cubic_discriminant(
        generic_orientation_coefficients[4],
        generic_orientation_coefficients[2],
        generic_orientation_coefficients[0],
    )
    require(
        generic_orientation_collision
        == (
            (4064256 - 9216) ** 2
            * (27878400 - 9216) ** 2
            * (27878400 - 4064256) ** 2
        )
        and generic_orientation_collision != 0,
        "generic orientation collision invariant failed",
    )
    generic_orientation_wall = orientation_wall(
        p_generic,
        q_generic,
        r_generic,
    )
    require(
        generic_orientation_collision
        == (
            2**24
            * q_generic**4
            * generic_orientation_wall**2
            * delta_generic
        ),
        "generic Theta factorization failed",
    )
    edge_values = {
        permutation: edge_value(generic_roots, permutation)
        for permutation in group
    }
    matching_values = {
        permutation: matching_value(generic_roots, permutation)
        for permutation in group
    }
    omega_values = {
        permutation: omega_value(generic_roots, permutation)
        for permutation in group
    }
    require(
        Counter(edge_values.values())
        == Counter(
            {
                -8: 4,
                -3: 4,
                -1: 4,
                1: 4,
                3: 4,
                8: 4,
            }
        ),
        "generic edge orbit failed",
    )
    require(
        Counter(matching_values.values())
        == Counter({1: 8, 9: 8, 64: 8}),
        "generic matching orbit failed",
    )
    require(
        Counter(omega_values.values())
        == Counter(
            {
                -5280: 4,
                -2016: 4,
                -96: 4,
                96: 4,
                2016: 4,
                5280: 4,
            }
        ),
        "generic orientation orbit failed",
    )
    for permutation in group:
        edge = edge_values[permutation]
        matching = matching_values[permutation]
        omega = omega_values[permutation]
        require(
            edge**2 == matching,
            "edge square does not equal matching root",
        )
        require(
            matching_resolvent(
                matching,
                p_generic,
                q_generic,
                r_generic,
            ) == 0,
            "generic matching root fails S",
        )
        require(
            omega**2
            == orientation_radicand(
                matching,
                p_generic,
                r_generic,
            ),
            "generic orientation square fails F",
        )
        derivative = matching_derivative(
            matching,
            p_generic,
            r_generic,
        )
        require(
            derivative**2
            * matching
            * orientation_radicand(
                matching,
                p_generic,
                r_generic,
            )
            == 16 * q_generic**2 * delta_generic,
            "generic D8 radicand product failed",
        )
        require(
            polynomial_value(
                edge_coefficients(
                    p_generic,
                    q_generic,
                    r_generic,
                ),
                edge,
            ) == 0,
            "generic edge sextic root failed",
        )
        require(
            polynomial_value(
                orientation_coefficients(
                    p_generic,
                    q_generic,
                    r_generic,
                ),
                omega,
            ) == 0,
            "generic orientation sextic root failed",
        )

    edge_representatives = orbit_representatives(
        group,
        lambda permutation: edge_values[permutation],
    )
    matching_representatives = orbit_representatives(
        group,
        lambda permutation: matching_values[permutation],
    )
    orientation_representatives = orbit_representatives(
        group,
        lambda permutation: omega_values[permutation],
    )
    require(
        (
            len(edge_representatives),
            len(matching_representatives),
            len(orientation_representatives),
        ) == (6, 3, 6),
        "generic representative counts failed",
    )
    require(
        polynomial_from_roots(
            tuple(
                edge_values[representative]
                for representative in edge_representatives
            )
        )
        == edge_coefficients(
            p_generic,
            q_generic,
            r_generic,
        ),
        "generic edge orbit polynomial failed",
    )
    require(
        polynomial_from_roots(
            tuple(
                omega_values[representative]
                for representative in orientation_representatives
            )
        )
        == orientation_coefficients(
            p_generic,
            q_generic,
            r_generic,
        ),
        "generic orientation orbit polynomial failed",
    )

    identity = tuple(range(4))
    edge_base = edge_values[identity]
    matching_base = matching_values[identity]
    omega_base = omega_values[identity]
    edge_stabilizer = frozenset(
        permutation
        for permutation in group
        if edge_values[permutation] == edge_base
    )
    matching_stabilizer = frozenset(
        permutation
        for permutation in group
        if matching_values[permutation] == matching_base
    )
    orientation_stabilizer = frozenset(
        permutation
        for permutation in group
        if omega_values[permutation] == omega_base
    )
    for label, subgroup in (
        ("edge stabilizer", edge_stabilizer),
        ("matching stabilizer", matching_stabilizer),
        ("orientation stabilizer", orientation_stabilizer),
    ):
        require_permutation_subgroup(subgroup, label)
    require(
        len(edge_stabilizer) == 4
        and Counter(
            permutation_order(element)
            for element in edge_stabilizer
        ) == Counter({2: 3, 1: 1}),
        "generic edge stabilizer is not V4",
    )
    require(
        len(orientation_stabilizer) == 4
        and Counter(
            permutation_order(element)
            for element in orientation_stabilizer
        ) == Counter({4: 2, 2: 1, 1: 1}),
        "generic orientation stabilizer is not C4",
    )
    require(
        len(matching_stabilizer) == 8
        and Counter(
            permutation_order(element)
            for element in matching_stabilizer
        ) == Counter({2: 5, 4: 2, 1: 1}),
        "generic matching stabilizer is not D8",
    )
    normal_klein = frozenset(
        permutation
        for permutation in group
        if permutation_cycle_type(permutation)
        in ((1, 1, 1, 1), (2, 2))
    )
    require(
        normal_klein.issubset(matching_stabilizer)
        and len(normal_klein) == 4,
        "normal quartic V4 failed",
    )
    index_two_bank = {
        edge_stabilizer,
        orientation_stabilizer,
        normal_klein,
    }
    require(
        len(index_two_bank) == 3
        and all(
            len(subgroup) == 4
            and subgroup.issubset(matching_stabilizer)
            for subgroup in index_two_bank
        ),
        "generic D8 index-two bank failed",
    )
    require(
        generated_permutation_group(
            tuple(edge_stabilizer)
            + tuple(orientation_stabilizer)
        ) == matching_stabilizer,
        "edge and orientation stabilizers do not generate D8",
    )
    character_edge = character_with_kernel(
        matching_stabilizer,
        edge_stabilizer,
    )
    character_orientation = character_with_kernel(
        matching_stabilizer,
        orientation_stabilizer,
    )
    character_discriminant = character_with_kernel(
        matching_stabilizer,
        normal_klein,
    )
    check_character(
        matching_stabilizer,
        character_edge,
        "edge character",
    )
    check_character(
        matching_stabilizer,
        character_orientation,
        "orientation character",
    )
    check_character(
        matching_stabilizer,
        character_discriminant,
        "discriminant character",
    )
    for element in matching_stabilizer:
        require(
            character_edge[element]
            * character_orientation[element]
            == character_discriminant[element],
            "D8 edge/orientation/discriminant character product failed",
        )

    # The compact separated control used by the theorem.
    theorem_control_roots = (-3, 0, 1, 2)
    theorem_control = check_specialization(
        theorem_control_roots,
        edge_representatives,
        orientation_representatives,
    )
    theorem_edge_values = Counter(
        edge_value(theorem_control_roots, permutation)
        for permutation in group
    )
    theorem_orientation_values = Counter(
        omega_value(theorem_control_roots, permutation)
        for permutation in group
    )
    require(
        (
            theorem_control["p"],
            theorem_control["q"],
            theorem_control["r"],
        ) == (-7, 6, 0)
        and theorem_edge_values
        == Counter(
            {
                -3: 4,
                -2: 4,
                -1: 4,
                1: 4,
                2: 4,
                3: 4,
            }
        )
        and theorem_orientation_values
        == Counter(
            {
                -120: 4,
                -96: 4,
                -24: 4,
                24: 4,
                96: 4,
                120: 4,
            }
        ),
        "theorem separated-root control failed",
    )
    theorem_control_wall = orientation_wall(-7, 6, 0)
    require(
        theorem_control_wall == 34992,
        "theorem separated-root J_or failed",
    )

    # Sharp specialization controls.
    even_control = check_specialization(
        (-5, -2, 2, 5),
        edge_representatives,
        orientation_representatives,
    )
    require(
        even_control["q"] == 0
        and even_control["discriminant"] != 0
        and even_control["edge_distinct"] == 5
        and even_control["orientation_distinct"] == 3,
        "q=0 specialization boundary failed",
    )
    repeated_control = check_specialization(
        (-4, -1, -1, 6),
        edge_representatives,
        orientation_representatives,
    )
    require(
        repeated_control["q"] != 0
        and repeated_control["discriminant"] == 0
        and repeated_control["edge_distinct"] == 4
        and repeated_control["orientation_distinct"] == 3,
        "Delta=0 specialization boundary failed",
    )
    orientation_collision = check_specialization(
        (-30, 2, 10, 18),
        edge_representatives,
        orientation_representatives,
    )
    require(
        orientation_collision["q"] != 0
        and orientation_collision["discriminant"] != 0
        and orientation_collision["edge_distinct"] == 6
        and orientation_collision["orientation_distinct"] == 4,
        "orientation-only collision boundary failed",
    )
    collision_orientation_coefficients = orientation_coefficients(
        orientation_collision["p"],
        orientation_collision["q"],
        orientation_collision["r"],
    )
    require(
        monic_cubic_discriminant(
            collision_orientation_coefficients[4],
            collision_orientation_coefficients[2],
            collision_orientation_coefficients[0],
        ) == 0,
        "orientation-only collision invariant is nonzero",
    )
    require(
        orientation_wall(
            orientation_collision["p"],
            orientation_collision["q"],
            orientation_collision["r"],
        ) == 0,
        "orientation-only collision does not lie on J_or=0",
    )

    # Coefficient-only wall control with the exact repeated quadratic factor.
    wall_p, wall_q, wall_r = (1, 4, -3)
    wall_delta = quartic_discriminant(
        wall_p,
        wall_q,
        wall_r,
    )
    wall_jor = orientation_wall(
        wall_p,
        wall_q,
        wall_r,
    )
    wall_edge_coefficients = edge_coefficients(
        wall_p,
        wall_q,
        wall_r,
    )
    wall_orientation_coefficients = orientation_coefficients(
        wall_p,
        wall_q,
        wall_r,
    )
    wall_edge_factorization = coefficient_product(
        (-1, 0, 1),
        (16, 0, 3, 0, 1),
    )
    repeated_orientation_factor = (-1280, 0, 1)
    wall_orientation_factorization = coefficient_product(
        coefficient_product(
            repeated_orientation_factor,
            repeated_orientation_factor,
        ),
        (14080, 0, 1),
    )
    require(
        wall_delta == -22000
        and wall_jor == 0
        and wall_edge_coefficients == wall_edge_factorization
        and wall_orientation_coefficients
        == wall_orientation_factorization,
        "explicit coefficient-wall factorization failed",
    )
    require(
        monic_division_remainder(
            wall_orientation_coefficients,
            repeated_orientation_factor,
        ) == (0,)
        and monic_division_remainder(
            coefficient_derivative(
                wall_orientation_coefficients,
            ),
            repeated_orientation_factor,
        ) == (0,),
        "orientation gcd collision factor failed",
    )

    print("QUARTIC SEXTIC-RESOLVENT TRIANGLE AUDIT - exact")
    print(
        "symbolic roots: "
        "S((a+c)^2)=0 Omega^2=F(u) Delta=Vandermonde^2"
    )
    print(
        "symbolic coefficient ring: "
        f"Norm(Z-F) terms={norm_polynomial.term_count()} "
        "disc(S)=Delta "
        "C(u)^2*u*F(u)-16*q^2*Delta mod S=(0,0,0)"
    )
    print(
        "orientation collision invariant: "
        f"Theta=2^24*q^4*J_or^2*Delta terms={orientation_collision_invariant.term_count()}"
    )
    print(
        "orientation sextic discriminant: "
        "Disc_Y(O)=2^18*q^4*Delta*Theta^2="
        "2^66*q^12*J_or^4*Delta^3"
    )
    print(
        "generic control roots=(-6,-2,3,5): "
        f"(p,q,r)=({p_generic},{q_generic},{r_generic}) "
        f"Delta={delta_generic}"
    )
    print(
        "generic edge orbit: "
        f"values={sorted(set(edge_values.values()))} "
        "multiplicity=4 stabilizer=V4 orders=[1,2,2,2]"
    )
    print(
        "generic matching orbit: "
        f"values={sorted(set(matching_values.values()))} "
        "multiplicity=8 stabilizer=D8 orders=[1,2,2,2,2,2,4,4]"
    )
    print(
        "generic orientation orbit: "
        f"values={sorted(set(omega_values.values()))} "
        "multiplicity=4 stabilizer=C4 orders=[1,2,4,4]"
    )
    print(
        "orbit polynomials: "
        "E(Y)=S(Y^2); O(Y)=Norm(Y^2-F)"
    )
    print(
        "generic edge coefficients: "
        f"{list(edge_coefficients(p_generic,q_generic,r_generic))}"
    )
    print(
        "generic orientation coefficients: "
        f"{list(generic_orientation_coefficients)}"
    )
    print(
        "generic radicand product: "
        f"16*q^2*Delta={16*q_generic**2*delta_generic} "
        f"J_or={generic_orientation_wall} "
        f"Theta={generic_orientation_collision}"
    )
    print(
        "D8 character triangle: "
        "chi_edge*chi_orientation=chi_discriminant on 8/8"
    )
    print(
        "theorem root control=(-3,0,1,2): "
        "(p,q,r)=(-7,6,0) "
        "edge=[-3,-2,-1,1,2,3] "
        "Omega=[-120,-96,-24,24,96,120] "
        f"J_or={theorem_control_wall}"
    )
    print(
        "q=0 boundary roots=(-5,-2,2,5): "
        f"Delta={even_control['discriminant']} "
        f"edge_distinct={even_control['edge_distinct']} "
        f"orientation_distinct={even_control['orientation_distinct']}"
    )
    print(
        "Delta=0 boundary roots=(-4,-1,-1,6): "
        f"q={repeated_control['q']} "
        f"edge_distinct={repeated_control['edge_distinct']} "
        f"orientation_distinct={repeated_control['orientation_distinct']}"
    )
    print(
        "orientation-only collision roots=(-30,2,10,18): "
        f"q={orientation_collision['q']} "
        f"Delta={orientation_collision['discriminant']} "
        "J_or=0 "
        f"edge_distinct={orientation_collision['edge_distinct']} "
        f"orientation_distinct={orientation_collision['orientation_distinct']}"
    )
    print(
        "coefficient wall (p,q,r)=(1,4,-3): "
        "Delta=-22000 J_or=0 "
        "E=(Y^2-1)(Y^4+3Y^2+16) "
        "O=(Y^2-1280)^2(Y^2+14080) gcd=Y^2-1280"
    )
    print(
        "specialization scope: "
        "edge separable iff q*Delta!=0; "
        "orientation separable iff q*Delta*J_or!=0 (characteristic zero)"
    )
    print("PASS")


if __name__ == "__main__":
    main()
