#!/usr/bin/env python3
"""Exact companion for THM-2871.

The script uses only the Python standard library.  It checks the polynomial
identities in a small exact multivariate ring, computes the universal
quartic resultant, verifies the two Newton-polygon controls, and certifies
the S4 specialization by exact finite-field factor degrees.

Run:
    python 04-computation/quartic_leading_face_integral_depression_thm2871.py
    python -O 04-computation/quartic_leading_face_integral_depression_thm2871.py
"""

from fractions import Fraction
from itertools import combinations, permutations


VARS = ("a", "p", "q", "b", "c", "d", "e", "x", "y", "u", "w")
NVARS = len(VARS)


class Poly:
    """Sparse exact polynomial in the fixed variables VARS."""

    def __init__(self, terms=None):
        cleaned = {}
        if terms:
            for monomial, coefficient in terms.items():
                value = Fraction(coefficient)
                if value:
                    cleaned[tuple(monomial)] = (
                        cleaned.get(tuple(monomial), Fraction(0)) + value
                    )
        self.terms = {m: c for m, c in cleaned.items() if c}

    @staticmethod
    def constant(value):
        value = Fraction(value)
        if not value:
            return Poly()
        return Poly({(0,) * NVARS: value})

    @staticmethod
    def variable(name):
        if name not in VARS:
            raise RuntimeError(f"unknown polynomial variable {name}")
        exponents = [0] * NVARS
        exponents[VARS.index(name)] = 1
        return Poly({tuple(exponents): Fraction(1)})

    def __add__(self, other):
        other = as_poly(other)
        result = dict(self.terms)
        for monomial, coefficient in other.terms.items():
            result[monomial] = result.get(monomial, Fraction(0)) + coefficient
            if not result[monomial]:
                del result[monomial]
        return Poly(result)

    def __radd__(self, other):
        return self + other

    def __neg__(self):
        return Poly({m: -c for m, c in self.terms.items()})

    def __sub__(self, other):
        return self + (-as_poly(other))

    def __rsub__(self, other):
        return as_poly(other) - self

    def __mul__(self, other):
        other = as_poly(other)
        result = {}
        for left_monomial, left_coefficient in self.terms.items():
            for right_monomial, right_coefficient in other.terms.items():
                monomial = tuple(
                    left_monomial[index] + right_monomial[index]
                    for index in range(NVARS)
                )
                result[monomial] = (
                    result.get(monomial, Fraction(0))
                    + left_coefficient * right_coefficient
                )
        return Poly(result)

    def __rmul__(self, other):
        return self * other

    def __pow__(self, exponent):
        if not isinstance(exponent, int) or exponent < 0:
            raise RuntimeError("polynomial exponent must be a nonnegative integer")
        result = Poly.constant(1)
        base = self
        power = exponent
        while power:
            if power & 1:
                result = result * base
            base = base * base
            power >>= 1
        return result

    def __truediv__(self, scalar):
        scalar = Fraction(scalar)
        if not scalar:
            raise RuntimeError("division by zero")
        return Poly({m: c / scalar for m, c in self.terms.items()})

    def substitute(self, substitutions):
        result = Poly()
        for monomial, coefficient in self.terms.items():
            term = Poly.constant(coefficient)
            for index, exponent in enumerate(monomial):
                if not exponent:
                    continue
                name = VARS[index]
                replacement = substitutions.get(name, Poly.variable(name))
                term = term * (as_poly(replacement) ** exponent)
            result = result + term
        return result

    def monomial_divide(self, name, exponent=1):
        index = VARS.index(name)
        result = {}
        for monomial, coefficient in self.terms.items():
            if monomial[index] < exponent:
                raise RuntimeError(
                    f"term {monomial} is not divisible by {name}^{exponent}"
                )
            divided = list(monomial)
            divided[index] -= exponent
            result[tuple(divided)] = coefficient
        return Poly(result)

    def __eq__(self, other):
        return self.terms == as_poly(other).terms

    def is_zero(self):
        return not self.terms


def as_poly(value):
    if isinstance(value, Poly):
        return value
    return Poly.constant(value)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sign_of_permutation(permutation):
    inversions = 0
    for left in range(len(permutation)):
        for right in range(left + 1, len(permutation)):
            if permutation[left] > permutation[right]:
                inversions += 1
    return -1 if inversions % 2 else 1


def determinant(matrix):
    size = len(matrix)
    if any(len(row) != size for row in matrix):
        raise RuntimeError("determinant requires a square matrix")
    total = Poly()
    for permutation in permutations(range(size)):
        term = Poly.constant(sign_of_permutation(permutation))
        for row in range(size):
            entry = matrix[row][permutation[row]]
            if entry.is_zero():
                term = Poly()
                break
            term = term * entry
        total = total + term
    return total


def sylvester_resultant(first_high, second_high):
    m = len(first_high) - 1
    n = len(second_high) - 1
    size = m + n
    zero = Poly()
    matrix = []
    for shift in range(n):
        row = [zero for _ in range(size)]
        for index, coefficient in enumerate(first_high):
            row[shift + index] = coefficient
        matrix.append(row)
    for shift in range(m):
        row = [zero for _ in range(size)]
        for index, coefficient in enumerate(second_high):
            row[shift + index] = coefficient
        matrix.append(row)
    return determinant(matrix)


def poly_trim_mod(coefficients, prime):
    result = [coefficient % prime for coefficient in coefficients]
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def poly_add_mod(left, right, prime):
    length = max(len(left), len(right))
    result = [0] * length
    for index in range(length):
        result[index] = (
            (left[index] if index < len(left) else 0)
            + (right[index] if index < len(right) else 0)
        ) % prime
    return poly_trim_mod(result, prime)


def poly_sub_mod(left, right, prime):
    length = max(len(left), len(right))
    result = [0] * length
    for index in range(length):
        result[index] = (
            (left[index] if index < len(left) else 0)
            - (right[index] if index < len(right) else 0)
        ) % prime
    return poly_trim_mod(result, prime)


def poly_mul_mod(left, right, prime):
    result = [0] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            result[left_index + right_index] = (
                result[left_index + right_index] + left_value * right_value
            ) % prime
    return poly_trim_mod(result, prime)


def poly_divmod_mod(numerator, denominator, prime):
    numerator = poly_trim_mod(numerator, prime)
    denominator = poly_trim_mod(denominator, prime)
    if denominator == [0]:
        raise RuntimeError("polynomial division by zero")
    if len(numerator) < len(denominator):
        return [0], numerator
    quotient = [0] * (len(numerator) - len(denominator) + 1)
    inverse = pow(denominator[-1], -1, prime)
    remainder = list(numerator)
    while remainder != [0] and len(remainder) >= len(denominator):
        shift = len(remainder) - len(denominator)
        coefficient = remainder[-1] * inverse % prime
        quotient[shift] = coefficient
        subtraction = [0] * shift + [
            coefficient * value % prime for value in denominator
        ]
        remainder = poly_sub_mod(remainder, subtraction, prime)
    return poly_trim_mod(quotient, prime), remainder


def poly_mod_mod(numerator, denominator, prime):
    return poly_divmod_mod(numerator, denominator, prime)[1]


def poly_gcd_mod(left, right, prime):
    left = poly_trim_mod(left, prime)
    right = poly_trim_mod(right, prime)
    while right != [0]:
        left, right = right, poly_mod_mod(left, right, prime)
    inverse = pow(left[-1], -1, prime)
    return [(value * inverse) % prime for value in left]


def poly_powmod_mod(base, exponent, modulus, prime):
    result = [1]
    base = poly_mod_mod(base, modulus, prime)
    power = exponent
    while power:
        if power & 1:
            result = poly_mod_mod(poly_mul_mod(result, base, prime), modulus, prime)
        base = poly_mod_mod(poly_mul_mod(base, base, prime), modulus, prime)
        power >>= 1
    return result


def factor_type_quartic_mod(coefficients, prime):
    polynomial = poly_trim_mod(coefficients, prime)
    require(len(polynomial) == 5, f"degree dropped modulo {prime}")
    derivative = [
        index * polynomial[index] % prime
        for index in range(1, len(polynomial))
    ]
    require(
        len(poly_gcd_mod(polynomial, derivative, prime)) == 1,
        f"specialization is not squarefree modulo {prime}",
    )
    x_polynomial = [0, 1]
    x_to_p = poly_powmod_mod(x_polynomial, prime, polynomial, prime)
    degree_one = len(
        poly_gcd_mod(
            polynomial,
            poly_sub_mod(x_to_p, x_polynomial, prime),
            prime,
        )
    ) - 1
    x_to_p2 = poly_powmod_mod(x_polynomial, prime * prime, polynomial, prime)
    degree_one_two = len(
        poly_gcd_mod(
            polynomial,
            poly_sub_mod(x_to_p2, x_polynomial, prime),
            prime,
        )
    ) - 1
    if degree_one == 0 and degree_one_two == 0:
        return (4,)
    if degree_one == 1 and degree_one_two == 1:
        return (1, 3)
    if degree_one == 2 and degree_one_two == 4:
        return (1, 1, 2)
    raise RuntimeError(
        "unexpected quartic factor census "
        f"mod {prime}: deg1={degree_one}, deg1or2={degree_one_two}"
    )


def compose_permutations(left, right):
    return tuple(left[right[index]] for index in range(len(left)))


def permutation_inverse(permutation):
    result = [0] * len(permutation)
    for index, image in enumerate(permutation):
        result[image] = index
    return tuple(result)


def generated_subgroup(generators):
    identity = tuple(range(4))
    group = {identity}
    frontier = [identity]
    expanded_generators = list(generators) + [
        permutation_inverse(generator) for generator in generators
    ]
    while frontier:
        current = frontier.pop()
        for generator in expanded_generators:
            product = compose_permutations(current, generator)
            if product not in group:
                group.add(product)
                frontier.append(product)
    return frozenset(group)


def cycle_type(permutation):
    seen = set()
    lengths = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        length = 0
        value = start
        while value not in seen:
            seen.add(value)
            length += 1
            value = permutation[value]
        lengths.append(length)
    return tuple(sorted(lengths))


def transitive(group):
    return {permutation[0] for permutation in group} == {0, 1, 2, 3}


def enumerate_subgroups_s4():
    all_permutations = list(permutations(range(4)))
    subgroups = {generated_subgroup([])}
    changed = True
    while changed:
        changed = False
        for group in list(subgroups):
            for permutation in all_permutations:
                enlarged = generated_subgroup(list(group) + [permutation])
                if enlarged not in subgroups:
                    subgroups.add(enlarged)
                    changed = True
    return subgroups


def lower_newton_hull(points):
    """Return lower hull with increasing abscissa."""
    hull = []
    for point in sorted(points):
        while len(hull) >= 2:
            x0, y0 = hull[-2]
            x1, y1 = hull[-1]
            x2, y2 = point
            left_slope = Fraction(y1 - y0, x1 - x0)
            right_slope = Fraction(y2 - y1, x2 - x1)
            if left_slope >= right_slope:
                hull.pop()
            else:
                break
        hull.append(point)
    return hull


def root_valuations_from_hull(points):
    hull = lower_newton_hull(points)
    valuations = []
    for left, right in zip(hull, hull[1:]):
        slope = Fraction(right[1] - left[1], right[0] - left[0])
        valuations.extend([-slope] * (right[0] - left[0]))
    return hull, tuple(valuations)


def main():
    a = Poly.variable("a")
    p = Poly.variable("p")
    q = Poly.variable("q")
    b = Poly.variable("b")
    c = Poly.variable("c")
    d = Poly.variable("d")
    e = Poly.variable("e")
    x = Poly.variable("x")
    y = Poly.variable("y")
    u = Poly.variable("u")
    w = Poly.variable("w")

    # The universal connected S4 hostile family.
    quartic_high = [a, Poly.constant(1), Poly.constant(0), p, q]
    quartic_derivative_high = [4 * a, Poly.constant(3), Poly.constant(0), p]
    resultant = sylvester_resultant(quartic_high, quartic_derivative_high)
    discriminant = resultant.monomial_divide("a")
    expected_discriminant = (
        256 * a**3 * q**3
        - 27 * a**2 * p**4
        - 192 * a**2 * p * q**2
        - 6 * a * p**2 * q
        - 4 * p**3
        - 27 * q**2
    )
    require(
        discriminant == expected_discriminant,
        "universal quartic discriminant/resultant mismatch",
    )
    leading_cubic_discriminant = -4 * p**3 - 27 * q**2
    require(
        discriminant.substitute({"a": 0}) == leading_cubic_discriminant,
        "leading-cubic discriminant specialization failed",
    )

    # THM-2867's leading-face matching blow-up, cross-multiplied so that
    # no localization or symbolic rational-function package is used.
    h = b * x**3 + c * x**2 + d * x + e
    u0 = (
        u**3
        + 8 * c * u**2
        + 16 * (c**2 + b * d) * u
        + 64 * (b * c * d - b**2 * e)
    )
    affine_u = -4 * (c + b * x)
    require(
        u0.substitute({"u": affine_u}) + 64 * b**2 * h == 0,
        "U0 is not the claimed affine copy of the leading cubic",
    )

    # Integral cubic depression at the second leading flag.
    h_flag = b * x**3 + b * c * x**2 + d * x + e
    shifted_h = h_flag.substitute({"x": y - c / 3})
    p_flag = d - b * c**2 / 3
    q_flag = e - c * d / 3 + 2 * b * c**3 / 27
    require(
        shifted_h == b * y**3 + p_flag * y + q_flag,
        "integral depression identity failed",
    )
    raw_delta_flag = (
        (b * c) ** 2 * d**2
        - 4 * b * d**3
        - 4 * (b * c) ** 3 * e
        - 27 * b**2 * e**2
        + 18 * b * (b * c) * d * e
    )
    factored_delta_flag = -4 * b * p_flag**3 - 27 * b**2 * q_flag**2
    require(
        raw_delta_flag == factored_delta_flag,
        "cubic discriminant flag factorization failed",
    )

    # THM-2473's extra square is an additional identity, not a consequence
    # of integral depression.  Here (a,p,q) stand for its target (a,b,g).
    sporadic_l = (
        27 * a**2 * q**2
        - 18 * a * p * q
        + 16 * a
        + p**3 * q
        - p**2
    )
    sporadic_p = 4 - 3 * p * q
    sporadic_square = 27 * a * q**2 - 9 * p * q + 8
    require(
        sporadic_p**3 + 27 * sporadic_l * q**2 == sporadic_square**2,
        "sporadic trisection square identity failed",
    )

    # The two sharp one-parameter controls.
    delta_c_unit = (-4 - 27 * b**2)
    delta_paired = (-4 * b - 27 * b**2)
    cubic_formula = (
        c**2 * d**2
        - 4 * b * d**3
        - 4 * c**3 * e
        - 27 * b**2 * e**2
        + 18 * b * c * d * e
    )
    require(
        cubic_formula.substitute({"c": 1, "d": 0, "e": 1})
        == delta_c_unit,
        "C-unit hostile discriminant failed",
    )
    require(
        cubic_formula.substitute({"c": 0, "d": 1, "e": 1})
        == delta_paired,
        "paired-escape hostile discriminant failed",
    )

    # Boundary specialization (p,q)=(-1,1).
    hostile_u0 = (u**3 + 16 * p * u - 64 * q).substitute({"p": -1, "q": 1})
    hostile_orientation = (
        w**6 + 6 * p * w**4 + 9 * p**2 * w**2 + 4 * p**3 + 27 * q**2
    ).substitute({"p": -1, "q": 1})
    require(
        hostile_u0 == u**3 - 16 * u - 64,
        "hostile U0 specialization failed",
    )
    require(
        hostile_orientation == w**6 - 6 * w**4 + 9 * w**2 + 23,
        "hostile orientation specialization failed",
    )

    # Newton polygons: coefficient points are (degree, valuation).
    stage_one_hull, stage_one_values = root_valuations_from_hull(
        [(0, 0), (1, 0), (3, 0), (4, 1)]
    )
    c_unit_hull, c_unit_values = root_valuations_from_hull(
        [(0, 0), (2, 0), (3, 1)]
    )
    paired_hull, paired_values = root_valuations_from_hull(
        [(0, 0), (1, 0), (3, 1)]
    )
    require(
        stage_one_values == (Fraction(0), Fraction(0), Fraction(0), Fraction(-1)),
        "stage-one Newton slopes failed",
    )
    require(
        c_unit_values == (Fraction(0), Fraction(0), Fraction(-1)),
        "C-unit Newton slopes failed",
    )
    require(
        paired_values == (Fraction(0), Fraction(-1, 2), Fraction(-1, 2)),
        "paired-escape Newton slopes failed",
    )

    # Exact Frobenius factor types for the S4 specialization
    # 2*x^4+x^3-x+1.
    specialized_quartic_low = [1, -1, 0, 1, 2]
    factor_types = {
        prime: factor_type_quartic_mod(specialized_quartic_low, prime)
        for prime in (5, 13, 17)
    }
    require(factor_types[5] == (4,), "mod-5 irreducibility certificate failed")
    require(
        factor_types[13] == (1, 1, 2),
        "mod-13 transposition certificate failed",
    )
    require(
        factor_types[17] == (1, 3),
        "mod-17 three-cycle certificate failed",
    )
    specialized_discriminant = 2673
    for prime in factor_types:
        require(
            specialized_discriminant % prime != 0,
            f"bad-prime factor type used at {prime}",
        )

    # Independently enumerate the subgroup boundary in S4.
    subgroups = enumerate_subgroups_s4()
    qualifying = []
    for group in subgroups:
        types = {cycle_type(permutation) for permutation in group}
        if (
            transitive(group)
            and (4,) in types
            and (1, 1, 2) in types
            and (1, 3) in types
        ):
            qualifying.append(group)
    require(
        len(qualifying) == 1 and len(qualifying[0]) == 24,
        "S4 cycle-type uniqueness failed",
    )

    print("THM-2871 exact companion")
    print(
        "universal_D4="
        "256*a^3*q^3-27*a^2*p^4-192*a^2*p*q^2"
        "-6*a*p^2*q-4*p^3-27*q^2"
    )
    print("D4_at_a0=-4*p^3-27*q^2")
    print("U0_affine_cubic=PASS")
    print("integral_depression=PASS")
    print("delta_flag=-b*(4*P^3+27*b*Q^2)")
    print("sporadic_square_law=P^3+27*L*g^2=S^2")
    print("hostile_C_unit_delta=-4-27*b^2")
    print("hostile_paired_delta=-b*(4+27*b)")
    print(f"newton_stage1_hull={stage_one_hull} roots={stage_one_values}")
    print(f"newton_C_unit_hull={c_unit_hull} roots={c_unit_values}")
    print(f"newton_paired_hull={paired_hull} roots={paired_values}")
    print(f"S4_specialization_factor_types={factor_types}")
    print(
        "S4_subgroup_gate="
        f"{len(subgroups)} subgroups; unique qualifying order {len(qualifying[0])}"
    )
    print("hostile_U0=u^3-16*u-64")
    print("hostile_orientation=w^6-6*w^4+9*w^2+23")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
