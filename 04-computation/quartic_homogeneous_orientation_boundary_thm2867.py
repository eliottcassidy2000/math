#!/usr/bin/env python3
"""Exact homogeneous quartic boundary audit for THM-2867.

The companion is dependency-free.  It implements sparse multivariate
integer-polynomial arithmetic directly, and every truth-bearing check uses
an explicit exception so optimized execution retains the full audit.
"""

from collections import Counter
from itertools import combinations, permutations


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


class Poly:
    __slots__ = ("arity", "terms")

    def __init__(self, arity, terms=None):
        self.arity = arity
        cleaned = {}
        for monomial, coefficient in (terms or {}).items():
            require(len(monomial) == arity, "wrong monomial arity")
            monomial = tuple(monomial)
            coefficient = int(coefficient)
            if coefficient:
                cleaned[monomial] = (
                    cleaned.get(monomial, 0) + coefficient
                )
        self.terms = {
            monomial: coefficient
            for monomial, coefficient in cleaned.items()
            if coefficient
        }

    @classmethod
    def constant(cls, arity, coefficient):
        if not coefficient:
            return cls(arity)
        return cls(arity, {(0,) * arity: coefficient})

    @classmethod
    def variable(cls, arity, position):
        monomial = [0] * arity
        monomial[position] = 1
        return cls(arity, {tuple(monomial): 1})

    def coerce(self, value):
        if isinstance(value, Poly):
            require(value.arity == self.arity, "mixed polynomial arities")
            return value
        return Poly.constant(self.arity, value)

    def __add__(self, other):
        other = self.coerce(other)
        terms = dict(self.terms)
        for monomial, coefficient in other.terms.items():
            terms[monomial] = terms.get(monomial, 0) + coefficient
        return Poly(self.arity, terms)

    def __radd__(self, other):
        return self + other

    def __neg__(self):
        return Poly(
            self.arity,
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
        return Poly(self.arity, terms)

    def __rmul__(self, other):
        return self * other

    def __pow__(self, exponent):
        require(
            isinstance(exponent, int) and exponent >= 0,
            "invalid polynomial exponent",
        )
        result = Poly.constant(self.arity, 1)
        base = self
        remaining = exponent
        while remaining:
            if remaining & 1:
                result = result * base
            base = base * base
            remaining //= 2
        return result

    def __eq__(self, other):
        return self.terms == self.coerce(other).terms

    def is_zero(self):
        return not self.terms


def divide_by_variable_power(polynomial, position, power):
    terms = {}
    for monomial, coefficient in polynomial.terms.items():
        require(
            monomial[position] >= power,
            "polynomial lacks requested variable power",
        )
        reduced = list(monomial)
        reduced[position] -= power
        terms[tuple(reduced)] = coefficient
    return Poly(polynomial.arity, terms)


def divide_coefficients(polynomial, divisor):
    require(divisor != 0, "zero coefficient divisor")
    terms = {}
    for monomial, coefficient in polynomial.terms.items():
        quotient, remainder = divmod(coefficient, divisor)
        require(remainder == 0, "inexact coefficient division")
        if quotient:
            terms[monomial] = quotient
    return Poly(polynomial.arity, terms)


def variable_order_and_lead(polynomial, position):
    require(not polynomial.is_zero(), "zero has no variable order")
    order = min(
        monomial[position] for monomial in polynomial.terms
    )
    terms = {}
    for monomial, coefficient in polynomial.terms.items():
        if monomial[position] == order:
            reduced = list(monomial)
            reduced[position] -= order
            terms[tuple(reduced)] = coefficient
    return order, Poly(polynomial.arity, terms)


def factorized_order_and_lead(factors, position, constant=1):
    require(factors, "factorized order needs at least one factor")
    arity = factors[0][0].arity
    order = 0
    lead = Poly.constant(arity, constant)
    for polynomial, exponent in factors:
        require(
            polynomial.arity == arity and exponent >= 0,
            "invalid factorized order input",
        )
        factor_order, factor_lead = variable_order_and_lead(
            polynomial, position
        )
        order += exponent * factor_order
        lead = lead * factor_lead**exponent
    return order, lead


def specialize_variable_zero(polynomial, position):
    return Poly(
        polynomial.arity,
        {
            monomial: coefficient
            for monomial, coefficient in polynomial.terms.items()
            if monomial[position] == 0
        },
    )


def monic_cubic_discriminant(coefficient_two, coefficient_one, constant):
    return (
        coefficient_two**2 * coefficient_one**2
        - 4 * coefficient_one**3
        - 4 * coefficient_two**3 * constant
        - 27 * constant**2
        + 18 * coefficient_two * coefficient_one * constant
    )


def quartic_discriminant(a, b, c, d, e):
    return (
        256 * a**3 * e**3
        - 192 * a**2 * b * d * e**2
        - 128 * a**2 * c**2 * e**2
        + 144 * a**2 * c * d**2 * e
        - 27 * a**2 * d**4
        + 144 * a * b**2 * c * e**2
        - 6 * a * b**2 * d**2 * e
        - 80 * a * b * c**2 * d * e
        + 18 * a * b * c * d**3
        + 16 * a * c**4 * e
        - 4 * a * c**3 * d**2
        - 27 * b**4 * e**2
        + 18 * b**3 * c * d * e
        - 4 * b**3 * d**3
        - 4 * b**2 * c**3 * e
        + b**2 * c**2 * d**2
    )


def cubic_discriminant(b, c, d, e):
    return (
        c**2 * d**2
        - 4 * b * d**3
        - 4 * c**3 * e
        - 27 * b**2 * e**2
        + 18 * b * c * d * e
    )


def edge_invariants(a, b, c, d, e):
    p_value = 8 * a * c - 3 * b**2
    q_value = b**3 - 4 * a * b * c + 8 * a**2 * d
    r_value = (
        -3 * b**4
        + 16 * a * b**2 * c
        - 64 * a**2 * b * d
        + 256 * a**3 * e
    )
    k_value = (
        3 * b**4
        - 16 * a * b**2 * c
        + 16 * a**2 * (c**2 + b * d)
        - 64 * a**3 * e
    )
    return p_value, q_value, r_value, k_value


def raw_j(p_value, q_value, r_value):
    return (
        r_value**3
        + 3 * p_value**2 * r_value**2
        - 36 * p_value * q_value**2 * r_value
        + 108 * q_value**4
    )


def polynomial_product(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for left_index, left_coefficient in enumerate(left):
        for right_index, right_coefficient in enumerate(right):
            result[left_index + right_index] += (
                left_coefficient * right_coefficient
            )
    return tuple(result)


def polynomial_from_roots(roots):
    coefficients = (1,)
    for root in roots:
        coefficients = polynomial_product(coefficients, (-root, 1))
    return coefficients


def polynomial_value(coefficients, value):
    result = 0
    for coefficient in reversed(coefficients):
        result = result * value + coefficient
    return result


def permutation_sign(permutation):
    inversions = sum(
        permutation[left] > permutation[right]
        for left in range(len(permutation))
        for right in range(left + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def cycle_lengths(mapping):
    require(
        sorted(mapping) == list(range(len(mapping))),
        "cycle input is not a permutation",
    )
    unseen = set(range(len(mapping)))
    lengths = []
    while unseen:
        start = min(unseen)
        current = start
        length = 0
        while current in unseen:
            unseen.remove(current)
            current = mapping[current]
            length += 1
        require(current == start, "permutation cycle did not close")
        lengths.append(length)
    return tuple(sorted(lengths))


def quartic_coefficients_from_roots(leading, roots):
    return (
        leading,
        -leading * sum(roots),
        leading
        * sum(
            roots[left] * roots[right]
            for left, right in combinations(range(4), 2)
        ),
        -leading
        * sum(
            roots[first] * roots[second] * roots[third]
            for first, second, third in combinations(range(4), 3)
        ),
        leading * roots[0] * roots[1] * roots[2] * roots[3],
    )


def root_v(leading, b_value, roots, permutation):
    return (
        2 * leading
        * (roots[permutation[0]] + roots[permutation[2]])
        + b_value
    )


def root_omega(roots, permutation):
    first, second, third, fourth = (
        roots[position] for position in permutation
    )
    left = first - third
    right = second - fourth
    return left * right * (left**2 - right**2)


def root_natural_w(leading, roots, permutation):
    return leading**3 * root_omega(roots, permutation)


def natural_numeric_invariants(a, b, c, d, e, p_value, q_value, r_value):
    i_value = c**2 - 3 * b * d + 12 * a * e
    l_value = (
        128 * a**2 * c * e
        - 36 * a**2 * d**2
        - 48 * a * b**2 * e
        + 4 * a * b * c * d
        + 3 * b**3 * d
        - b**2 * c**2
    )
    h_numerator = (
        r_value**2
        - p_value**2 * r_value
        + 12 * p_value * q_value**2
    )
    j_numerator = raw_j(p_value, q_value, r_value)
    h_value, h_remainder = divmod(h_numerator, 128 * a**2)
    j_value, j_remainder = divmod(j_numerator, 256 * a**3)
    require(
        h_remainder == 0 and j_remainder == 0,
        "numeric orientation quotient is not integral",
    )
    return i_value, l_value, h_value, j_value


def main():
    arity = 6
    a, b, c, d, e, t = tuple(
        Poly.variable(arity, position) for position in range(arity)
    )
    zero = Poly.constant(arity, 0)
    one = Poly.constant(arity, 1)

    p_value, q_value, r_value, k_value = edge_invariants(
        a, b, c, d, e
    )
    require(
        4 * k_value == p_value**2 - r_value,
        "K=(P^2-R)/4 failed",
    )
    discriminant = quartic_discriminant(a, b, c, d, e)
    delta_three = cubic_discriminant(b, c, d, e)
    i_zero = c**2 - 3 * b * d
    tau_three = 2 * c**3 - 9 * b * c * d + 27 * b**2 * e
    require(
        specialize_variable_zero(discriminant, 0)
        == b**2 * delta_three,
        "Disc(f)|A=0=B^2*Delta3 failed",
    )
    require(
        4 * i_zero**3 - tau_three**2
        == 27 * b**2 * delta_three,
        "4*I0^3-tau3^2=27*B^2*Delta3 failed",
    )

    # Matching cubic T(v), edge sextic T(V^2), and the A=0 blow-up.
    t_discriminant = monic_cubic_discriminant(
        p_value, k_value, -q_value**2
    )
    require(
        t_discriminant == 2**12 * a**6 * discriminant,
        "matching cubic discriminant failed",
    )
    edge_discriminant = 2**6 * q_value**2 * t_discriminant**2
    expected_edge_discriminant = (
        2**30 * a**12 * q_value**2 * discriminant**2
    )
    require(
        edge_discriminant == expected_edge_discriminant,
        "edge sextic discriminant failed",
    )

    u_two = 8 * c
    u_one = 16 * (c**2 + b * d - 4 * a * e)
    u_zero = 64 * (b * c * d - b**2 * e - a * d**2)
    substituted_t_coefficients = (
        b**6 + p_value * b**4 + k_value * b**2 - q_value**2,
        a * (3 * b**4 + 2 * p_value * b**2 + k_value),
        a**2 * (3 * b**2 + p_value),
        a**3,
    )
    expected_substitution = (
        a**3 * u_zero,
        a**3 * u_one,
        a**3 * u_two,
        a**3,
    )
    require(
        substituted_t_coefficients == expected_substitution,
        "T(B^2+A*t)=A^3*U(t) failed",
    )
    require(
        monic_cubic_discriminant(u_two, u_one, u_zero)
        == 2**12 * discriminant,
        "Disc(U)=2^12*Disc(f) failed",
    )
    edge_at_zero = (
        -specialize_variable_zero(q_value, 0) ** 2,
        zero,
        specialize_variable_zero(k_value, 0),
        zero,
        specialize_variable_zero(p_value, 0),
        zero,
        one,
    )
    require(
        edge_at_zero
        == (-b**6, zero, 3 * b**4, zero, -3 * b**2, zero, one),
        "edge at A=0 is not (V^2-B^2)^3",
    )

    # Integral natural coordinate W=A^3*Omega.
    i_value = c**2 - 3 * b * d + 12 * a * e
    l_value = (
        128 * a**2 * c * e
        - 36 * a**2 * d**2
        - 48 * a * b**2 * e
        + 4 * a * b * c * d
        + 3 * b**3 * d
        - b**2 * c**2
    )
    h_numerator = (
        r_value**2
        - p_value**2 * r_value
        + 12 * p_value * q_value**2
    )
    h_value = divide_coefficients(
        divide_by_variable_power(h_numerator, 0, 2),
        128,
    )
    j_numerator = raw_j(p_value, q_value, r_value)
    j_value = divide_coefficients(
        divide_by_variable_power(j_numerator, 0, 3),
        256,
    )
    require(
        p_value**2 + 3 * r_value == 64 * a**2 * i_value,
        "P^2+3R=64A^2*I failed",
    )
    require(
        p_value * r_value - 9 * q_value**2
        == 16 * a**2 * l_value,
        "PR-9Q^2=16A^2*L failed",
    )
    require(
        h_numerator == 128 * a**2 * h_value,
        "integral H quotient failed",
    )
    require(
        j_numerator == 256 * a**3 * j_value,
        "integral J quotient failed",
    )

    natural_four = 2 * h_value
    natural_two = -q_value**2 * i_value * l_value
    natural_zero = -q_value**4 * discriminant
    natural_cubic_discriminant = monic_cubic_discriminant(
        natural_four, natural_two, natural_zero
    )
    require(
        natural_cubic_discriminant
        == q_value**4 * j_value**2 * discriminant,
        "natural W^2-cubic discriminant failed",
    )

    # Verify once, universally, the cubic/even-sextic scaling exponents.
    scale_arity = 4
    cubic_a, cubic_b, cubic_c, scale = tuple(
        Poly.variable(scale_arity, position)
        for position in range(scale_arity)
    )
    abstract_cubic_discriminant = monic_cubic_discriminant(
        cubic_a, cubic_b, cubic_c
    )
    scaled_cubic_discriminant = monic_cubic_discriminant(
        scale**2 * cubic_a,
        scale**4 * cubic_b,
        scale**6 * cubic_c,
    )
    require(
        scaled_cubic_discriminant
        == scale**12 * abstract_cubic_discriminant,
        "monic cubic coordinate-scaling exponent failed",
    )
    abstract_even_discriminant = (
        -2**6 * cubic_c * abstract_cubic_discriminant**2
    )
    scaled_even_discriminant = (
        -2**6
        * scale**6
        * cubic_c
        * scaled_cubic_discriminant**2
    )
    require(
        scaled_even_discriminant
        == scale**30 * abstract_even_discriminant,
        "even sextic coordinate-scaling exponent failed",
    )

    # Raw THM-2864 g-orientation G=Omega_Y=256*A*W.
    raw_scale = 256 * a
    raw_four = raw_scale**2 * natural_four
    raw_two = raw_scale**4 * natural_two
    raw_zero = raw_scale**6 * natural_zero
    require(
        raw_four == 1024 * h_numerator
        and raw_two
        == (
            -2**22
            * q_value**2
            * (p_value**2 + 3 * r_value)
            * (p_value * r_value - 9 * q_value**2)
        )
        and raw_zero
        == -2**48 * a**6 * q_value**4 * discriminant,
        "G=256*A*W coefficient formulas failed",
    )
    raw_j_g = 1024 * j_numerator
    require(
        raw_j_g == 2**18 * a**3 * j_value,
        "raw g-orientation collision invariant failed",
    )

    # Universal A-orders and leading cubic residual.
    leading_expectations = (
        (p_value, 0, -3 * b**2, "P"),
        (q_value, 0, b**3, "Q"),
        (r_value, 0, -3 * b**4, "R"),
        (k_value, 0, 3 * b**4, "K"),
        (discriminant, 0, b**2 * delta_three, "Disc(f)"),
        (i_value, 0, i_zero, "I"),
        (l_value, 0, -b**2 * i_zero, "L"),
        (h_value, 0, -b**4 * i_zero, "H"),
        (j_value, 0, -b**6 * tau_three, "J"),
        (
            natural_four,
            0,
            -2 * b**4 * i_zero,
            "natural W4 coefficient",
        ),
        (
            natural_two,
            0,
            b**8 * i_zero**2,
            "natural W2 coefficient",
        ),
        (
            natural_zero,
            0,
            -b**14 * delta_three,
            "natural constant",
        ),
        (
            raw_four,
            2,
            -2**17 * b**4 * i_zero,
            "raw G4 coefficient",
        ),
        (
            raw_two,
            4,
            2**32 * b**8 * i_zero**2,
            "raw G2 coefficient",
        ),
        (
            raw_zero,
            6,
            -2**48 * b**14 * delta_three,
            "raw constant",
        ),
        (
            raw_j_g,
            3,
            -2**18 * b**6 * tau_three,
            "raw g J",
        ),
    )
    orders = {}
    for polynomial, expected_order, expected_lead, label in (
        leading_expectations
    ):
        order, lead = variable_order_and_lead(polynomial, 0)
        require(
            order == expected_order and lead == expected_lead,
            f"{label} A-order or leading term failed",
        )
        orders[label] = order

    factored_order_checks = (
        (
            (
                (q_value, 12),
                (j_value, 4),
                (discriminant, 3),
            ),
            2**6,
            0,
            2**6 * b**66 * tau_three**4 * delta_three**3,
            "natural orientation discriminant",
        ),
        (
            ((a, 12), (q_value, 2), (discriminant, 2)),
            2**30,
            12,
            2**30 * b**10 * delta_three**2,
            "edge discriminant",
        ),
        (
            (
                (a, 12),
                (q_value, 4),
                (j_value, 2),
                (discriminant, 1),
            ),
            2**96,
            12,
            2**96 * b**26 * tau_three**2 * delta_three,
            "raw orientation cubic discriminant",
        ),
        (
            (
                (a, 30),
                (q_value, 12),
                (j_value, 4),
                (discriminant, 3),
            ),
            2**246,
            30,
            2**246 * b**66 * tau_three**4 * delta_three**3,
            "raw orientation sextic discriminant",
        ),
    )
    for factors, constant, expected_order, expected_lead, label in (
        factored_order_checks
    ):
        order, lead = factorized_order_and_lead(
            factors, 0, constant
        )
        require(
            order == expected_order and lead == expected_lead,
            f"{label} factorized A-order or lead failed",
        )
        orders[label] = order

    natural_at_zero = (
        specialize_variable_zero(natural_zero, 0),
        zero,
        specialize_variable_zero(natural_two, 0),
        zero,
        specialize_variable_zero(natural_four, 0),
        zero,
        one,
    )
    expected_natural_at_zero = (
        -b**14 * delta_three,
        zero,
        b**8 * i_zero**2,
        zero,
        -2 * b**4 * i_zero,
        zero,
        one,
    )
    require(
        natural_at_zero == expected_natural_at_zero,
        "leading cubic orientation residual failed",
    )

    u_at_zero_coefficients = (
        specialize_variable_zero(u_zero, 0),
        specialize_variable_zero(u_one, 0),
        specialize_variable_zero(u_two, 0),
        one,
    )
    conjugate_u_polynomial = (
        (t + 4 * c) ** 3
        - 4 * c * (t + 4 * c) ** 2
        + 16 * b * d * (t + 4 * c)
        - 64 * b**2 * e
    )
    explicit_u_polynomial = (
        t**3
        + 8 * c * t**2
        + 16 * (c**2 + b * d) * t
        + 64 * (b * c * d - b**2 * e)
    )
    require(
        conjugate_u_polynomial == explicit_u_polynomial,
        "U0=-64B^2*h(-(t+4C)/(4B)) failed",
    )

    # Generic quartic control: all 24 labels and both sextic orbits.
    generic_roots = (-4, -1, 2, 6)
    generic_coefficients = quartic_coefficients_from_roots(
        3, generic_roots
    )
    ga, gb, gc, gd, ge = generic_coefficients
    gp, gq, gr, gk = edge_invariants(ga, gb, gc, gd, ge)
    gdisc = quartic_discriminant(ga, gb, gc, gd, ge)
    gi, gl, gh, gj = natural_numeric_invariants(
        ga, gb, gc, gd, ge, gp, gq, gr
    )
    require(
        generic_coefficients == (3, -9, -72, 84, 144)
        and (gp, gq, gr, gk)
        == (-1971, -2457, 1131165, 688419)
        and (gi, gl, gh, gj)
        == (
            12636,
            -15860124,
            -2827817244,
            2437528497221664,
        )
        and gdisc == 166659897600,
        "generic quartic invariant control failed",
    )
    s_four = tuple(permutations(range(4)))
    v_values = Counter(
        root_v(ga, gb, generic_roots, permutation)
        for permutation in s_four
    )
    w_values = Counter(
        root_natural_w(ga, generic_roots, permutation)
        for permutation in s_four
    )
    raw_values = Counter(
        16 * ga * root_natural_w(ga, generic_roots, permutation)
        for permutation in s_four
    )
    require(
        v_values
        == Counter(
            {-39: 4, -21: 4, -3: 4, 3: 4, 21: 4, 39: 4}
        ),
        "generic edge orbit failed",
    )
    require(
        w_values
        == Counter(
            {
                -73710: 4,
                -14742: 4,
                -2268: 4,
                2268: 4,
                14742: 4,
                73710: 4,
            }
        ),
        "generic natural orientation orbit failed",
    )
    require(
        raw_values
        == Counter(
            {
                -3538080: 4,
                -707616: 4,
                -108864: 4,
                108864: 4,
                707616: 4,
                3538080: 4,
            }
        ),
        "generic raw orientation orbit failed",
    )
    edge_coefficients = (-gq**2, 0, gk, 0, gp, 0, 1)
    natural_coefficients = (
        -gq**4 * gdisc,
        0,
        -gq**2 * gi * gl,
        0,
        2 * gh,
        0,
        1,
    )
    raw_coefficients = (
        (16 * ga) ** 6 * natural_coefficients[0],
        0,
        (16 * ga) ** 4 * natural_coefficients[2],
        0,
        (16 * ga) ** 2 * natural_coefficients[4],
        0,
        1,
    )
    require(
        polynomial_from_roots(sorted(v_values)) == edge_coefficients
        and polynomial_from_roots(sorted(w_values))
        == natural_coefficients
        and polynomial_from_roots(sorted(raw_values))
        == raw_coefficients,
        "generic sextic orbit polynomial failed",
    )
    for permutation in s_four:
        require(
            polynomial_value(
                edge_coefficients,
                root_v(ga, gb, generic_roots, permutation),
            )
            == 0
            and polynomial_value(
                natural_coefficients,
                root_natural_w(ga, generic_roots, permutation),
            )
            == 0,
            "generic labelled root substitution failed",
        )

    # A=0 squarefree cubic hostile: edge collapses, orientation does not.
    boundary_b, boundary_c, boundary_d, boundary_e = (1, 0, -1, 1)
    boundary_i = (
        boundary_c**2 - 3 * boundary_b * boundary_d
    )
    boundary_tau = (
        2 * boundary_c**3
        - 9 * boundary_b * boundary_c * boundary_d
        + 27 * boundary_b**2 * boundary_e
    )
    boundary_delta = cubic_discriminant(
        boundary_b, boundary_c, boundary_d, boundary_e
    )
    boundary_orientation = (
        -boundary_b**14 * boundary_delta,
        0,
        boundary_b**8 * boundary_i**2,
        0,
        -2 * boundary_b**4 * boundary_i,
        0,
        1,
    )
    boundary_discriminant = (
        2**6
        * boundary_b**66
        * boundary_tau**4
        * boundary_delta**3
    )
    require(
        (boundary_i, boundary_tau, boundary_delta) == (3, 27, -23)
        and boundary_orientation == (23, 0, 9, 0, -6, 0, 1)
        and boundary_discriminant != 0,
        "A=0 squarefree cubic hostile failed",
    )

    # Split cubic control: the ordered-pair S3 action is regular.
    cubic_roots = (-3, 1, 6)
    split_b, split_c, split_d, split_e = (1, -4, -15, 18)
    split_delta = cubic_discriminant(
        split_b, split_c, split_d, split_e
    )
    split_i = split_c**2 - 3 * split_b * split_d
    split_tau = (
        2 * split_c**3
        - 9 * split_b * split_c * split_d
        + 27 * split_b**2 * split_e
    )
    ordered_pairs = tuple(
        (left, right)
        for left in range(3)
        for right in range(3)
        if left != right
    )
    ordered_pair_index = {
        pair: index for index, pair in enumerate(ordered_pairs)
    }
    ordered_pair_roots = tuple(
        split_b**3
        * (cubic_roots[left] - cubic_roots[right])
        for left, right in ordered_pairs
    )
    split_orientation = (
        -split_b**14 * split_delta,
        0,
        split_b**8 * split_i**2,
        0,
        -2 * split_b**4 * split_i,
        0,
        1,
    )
    require(
        (split_i, split_tau, split_delta) == (61, -182, 32400)
        and sorted(ordered_pair_roots) == [-9, -5, -4, 4, 5, 9]
        and polynomial_from_roots(ordered_pair_roots)
        == split_orientation,
        "leading cubic ordered-pair root formula failed",
    )
    u0_split = (
        64
        * (
            split_b * split_c * split_d
            - split_b**2 * split_e
        ),
        16 * (split_c**2 + split_b * split_d),
        8 * split_c,
        1,
    )
    conjugate_roots = tuple(
        -4 * (split_c + split_b * root) for root in cubic_roots
    )
    require(
        sorted(conjugate_roots) == [-8, 12, 28]
        and polynomial_from_roots(conjugate_roots) == u0_split,
        "leading cubic U0 root conjugacy failed",
    )

    s_three = tuple(permutations(range(3)))
    base_pair = (0, 1)
    base_stabilizer = 0
    image_counts = Counter()
    cycle_census = Counter()
    for permutation in s_three:
        image_pair = (
            permutation[base_pair[0]],
            permutation[base_pair[1]],
        )
        image_counts[image_pair] += 1
        if image_pair == base_pair:
            base_stabilizer += 1
        induced = tuple(
            ordered_pair_index[
                (permutation[left], permutation[right])
            ]
            for left, right in ordered_pairs
        )
        lengths = cycle_lengths(induced)
        cycle_census[lengths] += 1
        ambient_sign = -1 if (6 - len(lengths)) % 2 else 1
        require(
            ambient_sign == permutation_sign(permutation),
            "ordered-pair ambient sign differs from cubic sign",
        )
        if permutation_sign(permutation) == -1:
            require(
                lengths == (2, 2, 2),
                "cubic transposition is not cycle type 2^3",
            )
        elif permutation != (0, 1, 2):
            require(
                lengths == (3, 3),
                "cubic 3-cycle is not cycle type 3^2",
            )
    require(
        base_stabilizer == 1
        and set(image_counts) == set(ordered_pairs)
        and set(image_counts.values()) == {1}
        and cycle_census
        == Counter({(2, 2, 2): 3, (3, 3): 2, (1, 1, 1, 1, 1, 1): 1}),
        "ordered-pair S3 action is not regular",
    )

    # Interior primitive-value collision: Disc(f) != 0 but J=0.
    collision_coefficients = (1, 0, 1, 4, -3)
    ca, cb, cc, cd, ce = collision_coefficients
    cp, cq, cr, ck = edge_invariants(ca, cb, cc, cd, ce)
    cdisc = quartic_discriminant(ca, cb, cc, cd, ce)
    ci, cl, ch, cj = natural_numeric_invariants(
        ca, cb, cc, cd, ce, cp, cq, cr
    )
    collision_orientation = (
        -cq**4 * cdisc,
        0,
        -cq**2 * ci * cl,
        0,
        2 * ch,
        0,
        1,
    )
    collision_factorization = polynomial_product(
        polynomial_product(
            (-1280, 0, 1),
            (-1280, 0, 1),
        ),
        (14080, 0, 1),
    )
    require(
        (cp, cq, cr, ck) == (8, 32, -768, 208)
        and cdisc == -22000
        and cj == 0
        and collision_orientation == collision_factorization,
        "interior J=0 collision hostile failed",
    )

    print("QUARTIC HOMOGENEOUS ORIENTATION BOUNDARY AUDIT - exact")
    print(
        "edge invariants: P=8AC-3B^2; "
        "Q=B^3-4ABC+8A^2D; "
        "R=-3B^4+16AB^2C-64A^2BD+256A^3E; "
        "4K=P^2-R"
    )
    print(
        "edge sextic: V^6+P V^4+K V^2-Q^2; "
        "Disc=2^30*A^12*Q^2*Disc(f)^2"
    )
    print(
        "matching blow-up: T(B^2+A*t)=A^3*U(t); "
        "U=t^3+8Ct^2+16(C^2+BD-4AE)t"
        "+64(BCD-B^2E-AD^2)"
    )
    print(
        "leading conjugacy: U0(t)=-64B^2*h(-(t+4C)/(4B)); "
        "roots t_i=-4(C+B*alpha_i)"
    )
    print(
        "natural coordinate: W=A^3*Omega; "
        "O=W^6+2H W^4-Q^2 I L W^2-Q^4 Disc(f)"
    )
    print(
        "natural quotients: P^2+3R=64A^2 I; "
        "PR-9Q^2=16A^2 L; "
        "R^2-P^2R+12PQ^2=128A^2 H; rawJ=256A^3 J"
    )
    print(
        "natural discriminants: "
        "Disc_(W^2)=Q^4*J^2*Disc(f); "
        "Disc_W(O)=2^6*Q^12*J^4*Disc(f)^3"
    )
    print(
        "raw g-coordinate: G=Omega_Y=256A*W; "
        "Disc_G=(256A)^30*Disc_W; generic A-order=30"
    )
    print(
        "raw A-orders: coefficients=(2,4,6); "
        "rawJ=3; cubicDisc=12; sexticDisc=30; edgeDisc=12"
    )
    print(
        "A=0 edge=(V^2-B^2)^3; "
        "O0=W^6-2B^4I0 W^4+B^8I0^2 W^2-B^14 Delta3"
    )
    print(
        "A=0 invariants: J|0=-B^6*tau3; "
        "4I0^3-tau3^2=27B^2*Delta3"
    )
    print(
        "generic quartic: roots=(-4,-1,2,6); "
        "(A,B,C,D,E)=(3,-9,-72,84,144); "
        f"(P,Q,R,K)=({gp},{gq},{gr},{gk}); Disc={gdisc}"
    )
    print(
        "generic edge roots: "
        f"{sorted(v_values)} multiplicity=4"
    )
    print(
        "generic natural orientation roots: "
        f"{sorted(w_values)} multiplicity=4"
    )
    print(
        "cubic ordered-pair control: roots=(-3,1,6); "
        "orientation roots=(-9,-5,-4,4,5,9); "
        "S3 regular; transposition=2^3; 3-cycle=3^2; signs agree"
    )
    print(
        "boundary hostile x^3-x+1: "
        "edge triple-collapses; "
        "O0=W^6-6W^4+9W^2+23 is squarefree"
    )
    print(
        "collision hostile x^4+x^2+4x-3: "
        "Disc=-22000; J=0; "
        "O=(W^2-1280)^2(W^2+14080)"
    )
    print(
        "scope: A-orders are universal generic orders and may rise on "
        "B*I0*tau3*Delta3=0; the blow-up and residual are coordinate "
        "identities, not a Keller divisor gate"
    )
    print("PASS")


if __name__ == "__main__":
    main()
