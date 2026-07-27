#!/usr/bin/env python3
"""Finite-place witnesses for the algebraic B--C ratio factors.

For each non-linear parameter polynomial in THM-2311's B--C bank, find a
prime where:

* the parameter polynomial remains irreducible, giving a finite residue
  field for the whole conjugacy packet;
* the u-discriminant of the rationalized trigonal curve has degree 12
  and gcd degree exactly 1 with its x-derivative;
* the leading cubic at infinity is separable; and
* one small integral x-specialization is an irreducible cubic in u.

The gcd degree can only increase on reduction.  Since THM-2311 already
puts each parameter packet on the repeated-discriminant resultant, these
witnesses force characteristic-zero gcd degree exactly one, hence ten
simple finite branch values.  The last two bullets certify a degree-three
geometrically irreducible cover after the rational smooth-point sidecar.
"""

from __future__ import annotations

from itertools import product


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def trim_prime(poly: list[int], prime: int) -> list[int]:
    result = [value % prime for value in poly]
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def div_prime(
    first: list[int],
    second: list[int],
    prime: int,
) -> tuple[list[int], list[int]]:
    numerator = trim_prime(first, prime)
    denominator = trim_prime(second, prime)
    quotient = [0] * max(1, len(numerator) - len(denominator) + 1)
    inverse = pow(denominator[-1], -1, prime)
    while numerator != [0] and len(numerator) >= len(denominator):
        shift = len(numerator) - len(denominator)
        coefficient = numerator[-1] * inverse % prime
        quotient[shift] = coefficient
        for index, value in enumerate(denominator):
            numerator[shift + index] -= coefficient * value
        numerator = trim_prime(numerator, prime)
    return trim_prime(quotient, prime), numerator


def irreducible_mod_prime(poly: list[int], prime: int) -> bool:
    reduced = trim_prime(poly, prime)
    degree = len(reduced) - 1
    if degree != len(poly) - 1:
        return False
    inverse = pow(reduced[-1], -1, prime)
    monic = [value * inverse % prime for value in reduced]
    for divisor_degree in range(1, degree // 2 + 1):
        for coefficients in product(range(prime), repeat=divisor_degree):
            divisor = list(coefficients) + [1]
            _, remainder = div_prime(monic, divisor, prime)
            if remainder == [0]:
                return False
    return True


class FiniteField:
    def __init__(self, prime: int, modulus: list[int]) -> None:
        self.prime = prime
        reduced = trim_prime(modulus, prime)
        inverse = pow(reduced[-1], -1, prime)
        self.modulus = tuple(value * inverse % prime for value in reduced)
        self.degree = len(self.modulus) - 1
        self.order = prime**self.degree
        self.zero = (0,) * self.degree
        self.one = (1,) + (0,) * (self.degree - 1)

    def element(self, coefficients: list[int] | tuple[int, ...]) -> tuple[int, ...]:
        values = [value % self.prime for value in coefficients]
        while len(values) > self.degree:
            leading = values[-1]
            shift = len(values) - self.degree - 1
            if leading:
                for index in range(self.degree):
                    values[shift + index] -= leading * self.modulus[index]
            values.pop()
            values = [value % self.prime for value in values]
        values += [0] * (self.degree - len(values))
        return tuple(value % self.prime for value in values)

    def scalar(self, value: int) -> tuple[int, ...]:
        return self.element([value])

    def add(
        self,
        first: tuple[int, ...],
        second: tuple[int, ...],
    ) -> tuple[int, ...]:
        return tuple(
            (left + right) % self.prime
            for left, right in zip(first, second, strict=True)
        )

    def neg(self, value: tuple[int, ...]) -> tuple[int, ...]:
        return tuple((-coefficient) % self.prime for coefficient in value)

    def sub(
        self,
        first: tuple[int, ...],
        second: tuple[int, ...],
    ) -> tuple[int, ...]:
        return self.add(first, self.neg(second))

    def mul(
        self,
        first: tuple[int, ...],
        second: tuple[int, ...],
    ) -> tuple[int, ...]:
        raw = [0] * (2 * self.degree - 1)
        for i, left in enumerate(first):
            for j, right in enumerate(second):
                raw[i + j] += left * right
        return self.element(raw)

    def pow(self, value: tuple[int, ...], exponent: int) -> tuple[int, ...]:
        result = self.one
        base = value
        remaining = exponent
        while remaining:
            if remaining & 1:
                result = self.mul(result, base)
            base = self.mul(base, base)
            remaining >>= 1
        return result

    def inv(self, value: tuple[int, ...]) -> tuple[int, ...]:
        require(value != self.zero, "attempted finite-field division by zero")
        return self.pow(value, self.order - 2)

    def elements(self):
        for coefficients in product(range(self.prime), repeat=self.degree):
            yield tuple(coefficients)


FieldElement = tuple[int, ...]
FieldPoly = list[FieldElement]


def trim_field(poly: FieldPoly, field: FiniteField) -> FieldPoly:
    result = list(poly)
    while len(result) > 1 and result[-1] == field.zero:
        result.pop()
    return result


def add_field(first: FieldPoly, second: FieldPoly, field: FiniteField) -> FieldPoly:
    size = max(len(first), len(second))
    return trim_field(
        [
            field.add(
                first[index] if index < len(first) else field.zero,
                second[index] if index < len(second) else field.zero,
            )
            for index in range(size)
        ],
        field,
    )


def neg_field(poly: FieldPoly, field: FiniteField) -> FieldPoly:
    return [field.neg(coefficient) for coefficient in poly]


def sub_field(first: FieldPoly, second: FieldPoly, field: FiniteField) -> FieldPoly:
    return add_field(first, neg_field(second, field), field)


def scale_field(
    poly: FieldPoly,
    scalar: FieldElement,
    field: FiniteField,
) -> FieldPoly:
    return trim_field([field.mul(scalar, coefficient) for coefficient in poly], field)


def mul_field(first: FieldPoly, second: FieldPoly, field: FiniteField) -> FieldPoly:
    result = [field.zero] * (len(first) + len(second) - 1)
    for i, left in enumerate(first):
        for j, right in enumerate(second):
            result[i + j] = field.add(
                result[i + j],
                field.mul(left, right),
            )
    return trim_field(result, field)


def power_field(poly: FieldPoly, exponent: int, field: FiniteField) -> FieldPoly:
    result = [field.one]
    base = list(poly)
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = mul_field(result, base, field)
        base = mul_field(base, base, field)
        remaining >>= 1
    return result


def derivative_field(poly: FieldPoly, field: FiniteField) -> FieldPoly:
    if len(poly) == 1:
        return [field.zero]
    return trim_field(
        [
            field.mul(field.scalar(index), coefficient)
            for index, coefficient in enumerate(poly)
        ][1:],
        field,
    )


def divide_field(
    first: FieldPoly,
    second: FieldPoly,
    field: FiniteField,
) -> tuple[FieldPoly, FieldPoly]:
    numerator = trim_field(first, field)
    denominator = trim_field(second, field)
    quotient = [field.zero] * max(1, len(numerator) - len(denominator) + 1)
    inverse = field.inv(denominator[-1])
    while numerator != [field.zero] and len(numerator) >= len(denominator):
        shift = len(numerator) - len(denominator)
        coefficient = field.mul(numerator[-1], inverse)
        quotient[shift] = coefficient
        numerator = sub_field(
            numerator,
            [field.zero] * shift + scale_field(denominator, coefficient, field),
            field,
        )
    return trim_field(quotient, field), trim_field(numerator, field)


def monic_field(poly: FieldPoly, field: FiniteField) -> FieldPoly:
    return scale_field(poly, field.inv(poly[-1]), field)


def gcd_field(first: FieldPoly, second: FieldPoly, field: FiniteField) -> FieldPoly:
    left, right = trim_field(first, field), trim_field(second, field)
    while right != [field.zero]:
        _, remainder = divide_field(left, right, field)
        left, right = right, remainder
    return monic_field(left, field)


def power_mod_field(
    base: FieldPoly,
    exponent: int,
    modulus: FieldPoly,
    field: FiniteField,
) -> FieldPoly:
    result = [field.one]
    current = list(base)
    remaining = exponent
    while remaining:
        if remaining & 1:
            _, result = divide_field(
                mul_field(result, current, field),
                modulus,
                field,
            )
        _, current = divide_field(
            mul_field(current, current, field),
            modulus,
            field,
        )
        remaining >>= 1
    return result


def discriminant_cubic(
    a: FieldPoly,
    b: FieldPoly,
    c: FieldPoly,
    d: FieldPoly,
    field: FiniteField,
) -> FieldPoly:
    return add_field(
        sub_field(
            sub_field(
                mul_field(power_field(b, 2, field), power_field(c, 2, field), field),
                scale_field(
                    mul_field(a, power_field(c, 3, field), field),
                    field.scalar(4),
                    field,
                ),
                field,
            ),
            scale_field(
                mul_field(power_field(b, 3, field), d, field),
                field.scalar(4),
                field,
            ),
            field,
        ),
        add_field(
            scale_field(
                mul_field(power_field(a, 2, field), power_field(d, 2, field), field),
                field.scalar(-27),
                field,
            ),
            scale_field(
                mul_field(mul_field(mul_field(a, b, field), c, field), d, field),
                field.scalar(18),
                field,
            ),
            field,
        ),
        field,
    )


def evaluate_field(
    poly: FieldPoly,
    value: FieldElement,
    field: FiniteField,
) -> FieldElement:
    result = field.zero
    for coefficient in reversed(poly):
        result = field.add(field.mul(result, value), coefficient)
    return result


def scalar_cubic_discriminant(
    a: FieldElement,
    b: FieldElement,
    c: FieldElement,
    d: FieldElement,
    field: FiniteField,
) -> FieldElement:
    def product_of(*values: FieldElement) -> FieldElement:
        result = field.one
        for value in values:
            result = field.mul(result, value)
        return result

    result = product_of(b, b, c, c)
    result = field.sub(
        result,
        field.mul(field.scalar(4), product_of(a, c, c, c)),
    )
    result = field.sub(
        result,
        field.mul(field.scalar(4), product_of(b, b, b, d)),
    )
    result = field.sub(
        result,
        field.mul(field.scalar(27), product_of(a, a, d, d)),
    )
    return field.add(
        result,
        field.mul(field.scalar(18), product_of(a, b, c, d)),
    )


def curve_coefficients(
    field: FiniteField,
    t: FieldElement,
    plane: str,
) -> tuple[FieldPoly, FieldPoly, FieldPoly, FieldPoly]:
    """Return a(x),b(x),c(x),d(x) for H=a u^3+b u^2+c u+d."""

    s = field.scalar
    t2, t3 = field.pow(t, 2), field.pow(t, 3)
    if plane == "BC":
        return (
            [field.mul(s(-26040609), t3)],
            [
                field.mul(s(49601160), t3),
                field.zero,
                field.mul(s(1607445), t2),
            ],
            [
                field.mul(s(-20995200), t3),
                field.zero,
                field.mul(s(-2857680), t2),
                field.zero,
                field.mul(s(-138915), t),
            ],
            [
                field.zero,
                field.mul(s(-5598720), t3),
                field.mul(s(777600), t2),
                field.mul(s(-435456), t2),
                field.mul(s(78120), t),
                field.zero,
                s(1127),
            ],
        )
    if plane == "BW":
        return (
            [field.mul(s(-26040609), t3)],
            [
                field.mul(s(49601160), t3),
                field.zero,
                field.mul(s(1607445), t2),
            ],
            [
                field.mul(s(-20995200), t3),
                field.zero,
                field.mul(s(-2857680), t2),
                field.zero,
                field.mul(s(-138915), t),
            ],
            [
                field.zero,
                field.mul(s(-5878656), t3),
                field.mul(s(777600), t2),
                field.zero,
                field.mul(s(78120), t),
                field.zero,
                s(1127),
            ],
        )
    if plane == "CD":
        return (
            [field.mul(s(-26040609), t)],
            [field.zero, field.zero, field.mul(s(1607445), t)],
            [
                s(-52907904),
                field.zero,
                field.zero,
                field.zero,
                field.mul(s(-138915), t),
            ],
            [
                field.zero,
                field.zero,
                s(1959552),
                s(-435456),
                field.zero,
                field.zero,
                field.mul(s(1127), t),
            ],
        )
    if plane == "DW":
        return (
            [field.mul(s(-26040609), t)],
            [field.zero, field.zero, field.mul(s(1607445), t)],
            [
                s(-52907904),
                field.zero,
                field.zero,
                field.zero,
                field.mul(s(-138915), t),
            ],
            [
                field.zero,
                s(-5878656),
                s(1959552),
                field.zero,
                field.zero,
                field.zero,
                field.mul(s(1127), t),
            ],
        )
    raise ValueError(f"unknown plane {plane}")


def witness_for(
    parameter_polynomial: list[int],
    plane: str = "BC",
) -> tuple[int, FiniteField, int, int, bool, bool, int]:
    primes = (5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43)
    diagnostics: list[tuple[int, int, int, bool, bool, bool, bool]] = []
    for prime in primes:
        if not irreducible_mod_prime(parameter_polynomial, prime):
            continue
        field = FiniteField(prime, parameter_polynomial)
        t = field.element([0, 1])
        parameter_value = field.zero
        for coefficient in reversed(parameter_polynomial):
            parameter_value = field.add(
                field.mul(parameter_value, t),
                field.scalar(coefficient),
            )
        require(parameter_value == field.zero, "parameter root reduction failed")
        if t == field.zero:
            continue
        require(
            field.pow(t, field.order) == t,
            "finite-field Frobenius control failed",
        )
        a, b, c, d = curve_coefficients(field, t, plane)
        branch = discriminant_cubic(a, b, c, d, field)
        branch_degree = len(branch) - 1
        branch_derivative = derivative_field(branch, field)
        if branch == [field.zero] or branch_derivative == [field.zero]:
            continue
        branch_gcd = gcd_field(branch, branch_derivative, field)
        gcd_degree = len(branch_gcd) - 1

        infinity_discriminant = scalar_cubic_discriminant(
            a[0],
            b[2],
            c[4],
            d[6],
            field,
        )
        infinity_separable = infinity_discriminant != field.zero

        no_root = False
        specialization = 0
        for candidate in range(1, 13):
            base_value = field.scalar(candidate)
            specialized = tuple(
                evaluate_field(coefficient, base_value, field)
                for coefficient in (a, b, c, d)
            )
            specialized_poly = [
                specialized[3],
                specialized[2],
                specialized[1],
                specialized[0],
            ]
            frobenius_minus_identity = sub_field(
                power_mod_field(
                    [field.zero, field.one],
                    field.order,
                    specialized_poly,
                    field,
                ),
                [field.zero, field.one],
                field,
            )
            root_gcd = gcd_field(
                specialized_poly,
                frobenius_minus_identity,
                field,
            )
            if len(root_gcd) == 1:
                no_root = True
                specialization = candidate
                break
        smooth_origin = c[0] != field.zero
        cubic_leading_unit = a[0] != field.zero
        diagnostics.append(
            (
                prime,
                branch_degree,
                gcd_degree,
                infinity_separable,
                no_root,
                smooth_origin,
                cubic_leading_unit,
            )
        )

        if (
            cubic_leading_unit
            and branch_degree == 12
            and gcd_degree == 1
            and infinity_separable
            and no_root
            and smooth_origin
        ):
            return (
                prime,
                field,
                branch_degree,
                gcd_degree,
                infinity_separable,
                no_root,
                specialization,
            )
    raise RuntimeError(f"no finite-place witness found: {diagnostics}")


def main() -> None:
    require(
        irreducible_mod_prime([1, 0, 1], 3),
        "positive irreducibility control failed",
    )
    require(
        not irreducible_mod_prime([1, 0, 1], 5),
        "reducible-polynomial hostile control failed",
    )
    packets = {
        "BC_linear_1": [
            2000,
            15309,
        ],
        "BC_linear_2": [
            125,
            1134,
        ],
        "BC_quadratic": [
            3321125,
            -161754894,
            -2812385772,
        ],
        "BC_quintic": [
            410644531250000,
            -18114791748046875,
            -545436228093750000,
            -4951165276923468750,
            -18946967714644599000,
            -26529827304546537363,
        ],
    }
    print("finite_field_controls=pass")
    total_roots = 0
    for name, polynomial in packets.items():
        (
            prime,
            field,
            branch_degree,
            gcd_degree,
            infinity_separable,
            no_root,
            specialization,
        ) = witness_for(polynomial)
        print(f"packet={name}")
        degree = len(polynomial) - 1
        total_roots += degree
        print(f"parameter_degree={degree}")
        print(f"witness_prime={prime}")
        print(
            "monic_parameter_mod_prime="
            + ",".join(str(value) for value in field.modulus)
        )
        print(f"residue_field_order={field.order}")
        print("cubic_leading_unit=1")
        print(f"branch_degree={branch_degree}")
        print(f"branch_derivative_gcd_degree={gcd_degree}")
        print(f"infinity_separable={int(infinity_separable)}")
        print(f"irreducible_specialization_x={specialization}")
        print(f"specialized_cubic_no_residue_root={int(no_root)}")
    print(f"candidate_roots_covered={total_roots}")
    print("status=BC_ALGEBRAIC_FINITE_PLACE_GENUS_FLOOR_EXACT")


if __name__ == "__main__":
    main()
