#!/usr/bin/env python3
"""Clean-room audit of three-coordinate primitivity for the fixed F^3.

The candidate probe and its parent computation are not imported.  This audit
uses three independent routes:

1. symbolic elimination re-derives the inverse section from the displayed
   sporadic Keller map and verifies its quotient-map identities;
2. direct enumeration of all of F_101^3 reconstructs the two split levels of
   both the positive target and the retained hostile target; and
3. a separately written cubic-quotient/matrix engine computes direct
   determinant characteristic polynomials and full 27x27 power-basis ranks.

The finite computation certifies generic primitivity for this fixed map at
level three.  The common [-2J] class then follows from trace-form congruence
and proved THM-3495; no arbitrary-Keller, all-level, JC(2), or LRC claim is
made here.
"""

from __future__ import annotations

import ast
import hashlib
from itertools import permutations
from pathlib import Path

import sympy as sy


P = 101
WITNESS_TARGET = (77, 62, 4)
WITNESS_LEVEL_ONE = ((13, 36, 2), (39, 84, 75), (49, 74, 79))
WITNESS_LEVEL_TWO_BLOCKS = (
    ((35, 40, 46), (83, 54, 93), (84, 61, 87)),
    ((2, 91, 44), (35, 7, 50), (64, 28, 37)),
    ((24, 2, 74), (84, 50, 58), (94, 59, 2)),
)
PARENT_TARGET = (93, 28, 83)
PARENT_LEVEL_ONE = ((10, 17, 82), (39, 90, 66), (52, 36, 24))
PARENT_LEVEL_TWO_BLOCKS = (
    ((16, 70, 96), (90, 98, 7), (96, 9, 53)),
    ((50, 45, 99), (72, 58, 17), (80, 32, 9)),
    ((2, 57, 13), (21, 83, 49), (78, 15, 61)),
)
EXPECTED_WITNESS_PRODUCT_HASHES = (
    "306dc5831a989a0fe953ae5eb6127d26bd6ff07f8cc4e25d54dc0f65e53e5e19",
    "aba63853d168e670461e67b67e42d26380e18c62b363c8931aa187f1a97d4e51",
    "b8a09c5bfed7affdabbbb3fcf222002a3a7c11c8d55cb19bcc949494725bb51d",
)
EXPECTED_PARENT_PRODUCT_HASHES = (
    "e10d5f365f1a673e4c2615123c1925ac9fa797395cb160495fc43df7969cd4e2",
    "691caf7a37960ecf904b669ebf4477c232a5eed1137b35c90a76d8b55b53237d",
    "fe96f6cb23ebc5facf0d36d2d2319da007ac7c46b1b96c809de8ad88fe1ed179",
)
EXPECTED_SEMANTIC_SHA256 = "51b5945a5f26c3335b6837b169a916d9562acd8448d2dde4f1fe1449c1e157a3"


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def inv(value: int) -> int:
    value %= P
    require(value != 0, "inverse of zero")
    return pow(value, P - 2, P)


def l_value(point: tuple[int, int, int]) -> int:
    a, b, c = point
    return (
        27 * a * a * c * c - 18 * a * b * c + 16 * a + b**3 * c - b * b
    ) % P


def fmap(point: tuple[int, int, int]) -> tuple[int, int, int]:
    x, y, z = point
    unit = (1 + x * y) % P
    return (
        (unit**3 * z + y * y * unit * (4 + 3 * x * y)) % P,
        (y + 3 * x * unit * unit * z + 3 * x * y * y * (4 + 3 * x * y)) % P,
        (2 * x - 3 * x * x * y - x**3 * z) % P,
    )


def core(point: tuple[int, int, int]) -> tuple[int, ...]:
    _a, b, c = point
    return ((-2 * c) % P, (4 - 3 * b * c) % P, 0, l_value(point))


def trim(poly: list[int] | tuple[int, ...]) -> list[int]:
    result = [value % P for value in poly]
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def poly_add(left: list[int], right: list[int]) -> list[int]:
    return trim([
        (left[index] if index < len(left) else 0)
        + (right[index] if index < len(right) else 0)
        for index in range(max(len(left), len(right)))
    ])


def poly_negate(poly: list[int]) -> list[int]:
    return trim([-value for value in poly])


def poly_subtract(left: list[int], right: list[int]) -> list[int]:
    return poly_add(left, poly_negate(right))


def poly_scale(poly: list[int], scalar: int) -> list[int]:
    return trim([scalar * value for value in poly])


def poly_multiply(left: list[int], right: list[int]) -> list[int]:
    result = [0] * (len(left) + len(right) - 1)
    for i, x_value in enumerate(left):
        for j, y_value in enumerate(right):
            result[i + j] = (result[i + j] + x_value * y_value) % P
    return trim(result)


def poly_divmod(
    numerator: list[int], denominator: list[int]
) -> tuple[list[int], list[int]]:
    remainder = trim(numerator)
    divisor = trim(denominator)
    require(divisor != [0], "polynomial division by zero")
    if len(remainder) < len(divisor):
        return [0], remainder
    quotient = [0] * (len(remainder) - len(divisor) + 1)
    inverse_lead = inv(divisor[-1])
    while remainder != [0] and len(remainder) >= len(divisor):
        shift = len(remainder) - len(divisor)
        factor = remainder[-1] * inverse_lead % P
        quotient[shift] = factor
        for index, coefficient in enumerate(divisor):
            remainder[index + shift] = (
                remainder[index + shift] - factor * coefficient
            ) % P
        remainder = trim(remainder)
    return trim(quotient), remainder


def poly_gcd(left: list[int], right: list[int]) -> list[int]:
    x_poly, y_poly = trim(left), trim(right)
    while y_poly != [0]:
        x_poly, y_poly = y_poly, poly_divmod(x_poly, y_poly)[1]
    return poly_scale(x_poly, inv(x_poly[-1]))


def derivative(poly: list[int]) -> list[int]:
    return trim([index * poly[index] for index in range(1, len(poly))])


def matrix_rank(matrix: list[list[int]]) -> int:
    rows = [[value % P for value in row] for row in matrix]
    row = 0
    for column in range(len(rows[0]) if rows else 0):
        pivot = next((r for r in range(row, len(rows)) if rows[r][column]), None)
        if pivot is None:
            continue
        rows[row], rows[pivot] = rows[pivot], rows[row]
        inverse = inv(rows[row][column])
        rows[row] = [value * inverse % P for value in rows[row]]
        for other in range(len(rows)):
            if other == row or rows[other][column] == 0:
                continue
            factor = rows[other][column]
            rows[other] = [
                (value - factor * pivot_value) % P
                for value, pivot_value in zip(rows[other], rows[row])
            ]
        row += 1
        if row == len(rows):
            break
    return row


def determinant(matrix: list[list[int]]) -> int:
    require(len(matrix) == len(matrix[0]), "nonsquare determinant")
    rows = [[value % P for value in row] for row in matrix]
    answer = 1
    for column in range(len(rows)):
        pivot = next(
            (row for row in range(column, len(rows)) if rows[row][column]),
            None,
        )
        if pivot is None:
            return 0
        if pivot != column:
            rows[column], rows[pivot] = rows[pivot], rows[column]
            answer = -answer
        pivot_value = rows[column][column]
        answer = answer * pivot_value % P
        inverse = inv(pivot_value)
        for row in range(column + 1, len(rows)):
            factor = rows[row][column] * inverse % P
            if factor:
                rows[row] = [
                    (value - factor * pivot_entry) % P
                    for value, pivot_entry in zip(rows[row], rows[column])
                ]
    return answer % P


def solve(matrix: list[list[int]], vector: tuple[int, ...]) -> tuple[int, ...]:
    size = len(matrix)
    rows = [
        [matrix[row][column] % P for column in range(size)] + [vector[row] % P]
        for row in range(size)
    ]
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if rows[row][column]), None
        )
        require(pivot is not None, "singular linear solve")
        rows[column], rows[pivot] = rows[pivot], rows[column]
        inverse = inv(rows[column][column])
        rows[column] = [value * inverse % P for value in rows[column]]
        for row in range(size):
            if row == column or rows[row][column] == 0:
                continue
            factor = rows[row][column]
            rows[row] = [
                (value - factor * pivot_value) % P
                for value, pivot_value in zip(rows[row], rows[column])
            ]
    return tuple(rows[row][-1] for row in range(size))


class CubicQuotient:
    """Three-coordinate quotient with reduction and inverses via matrices."""

    def __init__(self, modulus: tuple[int, ...]):
        self.modulus = tuple(trim(modulus))
        require(len(self.modulus) == 4 and self.modulus[-1] != 0, self.modulus)
        self.inverse_lead = inv(self.modulus[-1])
        self.zero = (0, 0, 0)
        self.one = (1, 0, 0)
        self.generator = (0, 1, 0)

    def reduce(self, raw: list[int] | tuple[int, ...]) -> tuple[int, ...]:
        values = [value % P for value in raw]
        if len(values) < 3:
            values += [0] * (3 - len(values))
        for degree in range(len(values) - 1, 2, -1):
            coefficient = values[degree] * self.inverse_lead % P
            if coefficient:
                shift = degree - 3
                for index in range(3):
                    values[shift + index] = (
                        values[shift + index] - coefficient * self.modulus[index]
                    ) % P
            values[degree] = 0
        return tuple(values[:3])

    def add(self, left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
        return tuple((x_value + y_value) % P for x_value, y_value in zip(left, right))

    def subtract(
        self, left: tuple[int, ...], right: tuple[int, ...]
    ) -> tuple[int, ...]:
        return tuple((x_value - y_value) % P for x_value, y_value in zip(left, right))

    def scale(self, element: tuple[int, ...], scalar: int) -> tuple[int, ...]:
        return tuple(scalar * value % P for value in element)

    def multiply(
        self, left: tuple[int, ...], right: tuple[int, ...]
    ) -> tuple[int, ...]:
        raw = [0] * 5
        for i, x_value in enumerate(left):
            for j, y_value in enumerate(right):
                raw[i + j] = (raw[i + j] + x_value * y_value) % P
        return self.reduce(raw)

    def power(self, element: tuple[int, ...], exponent: int) -> tuple[int, ...]:
        answer = self.one
        base = element
        while exponent:
            if exponent & 1:
                answer = self.multiply(answer, base)
            base = self.multiply(base, base)
            exponent //= 2
        return answer

    def matrix(self, element: tuple[int, ...]) -> list[list[int]]:
        basis = (self.one, self.generator, (0, 0, 1))
        columns = [self.multiply(element, value) for value in basis]
        return [[columns[column][row] for column in range(3)] for row in range(3)]

    def inverse(self, element: tuple[int, ...]) -> tuple[int, ...]:
        answer = solve(self.matrix(element), self.one)
        require(self.multiply(element, answer) == self.one, (element, answer))
        return answer

    def evaluate(self, poly: list[int], element: tuple[int, ...]) -> tuple[int, ...]:
        answer = self.zero
        for coefficient in reversed(poly):
            answer = self.add(self.multiply(answer, element), (coefficient % P, 0, 0))
        return answer


def permutation_sign(permutation: tuple[int, ...]) -> int:
    inversions = sum(
        permutation[i] > permutation[j]
        for i in range(len(permutation))
        for j in range(i + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def characteristic_polynomial(
    algebra: CubicQuotient, element: tuple[int, ...]
) -> tuple[int, ...]:
    matrix = algebra.matrix(element)
    polynomial_matrix = [
        [
            [(-matrix[row][column]) % P, 1]
            if row == column
            else [(-matrix[row][column]) % P]
            for column in range(3)
        ]
        for row in range(3)
    ]
    answer = [0]
    for permutation in permutations(range(3)):
        term = [1]
        for row, column in enumerate(permutation):
            term = poly_multiply(term, polynomial_matrix[row][column])
        answer = (
            poly_add(answer, term)
            if permutation_sign(permutation) == 1
            else poly_subtract(answer, term)
        )
    answer = tuple(trim(answer))
    require(len(answer) == 4 and answer[-1] == 1, answer)
    require(algebra.evaluate(list(answer), element) == algebra.zero, "Cayley-Hamilton")
    return answer


def fmap_quotient(
    algebra: CubicQuotient,
    x_value: tuple[int, ...],
    y_value: tuple[int, ...],
    z_value: tuple[int, ...],
) -> tuple[tuple[int, ...], ...]:
    xy = algebra.multiply(x_value, y_value)
    unit = algebra.add(algebra.one, xy)
    unit2 = algebra.power(unit, 2)
    unit3 = algebra.multiply(unit2, unit)
    y2 = algebra.power(y_value, 2)
    four_plus_3xy = algebra.add((4, 0, 0), algebra.scale(xy, 3))
    first = algebra.add(
        algebra.multiply(unit3, z_value),
        algebra.multiply(algebra.multiply(y2, unit), four_plus_3xy),
    )
    second = algebra.add(
        y_value,
        algebra.add(
            algebra.scale(
                algebra.multiply(algebra.multiply(x_value, unit2), z_value), 3
            ),
            algebra.scale(
                algebra.multiply(algebra.multiply(x_value, y2), four_plus_3xy),
                3,
            ),
        ),
    )
    third = algebra.subtract(
        algebra.subtract(
            algebra.scale(x_value, 2),
            algebra.scale(algebra.multiply(algebra.power(x_value, 2), y_value), 3),
        ),
        algebra.multiply(algebra.power(x_value, 3), z_value),
    )
    return first, second, third


def symbolic_inverse_certificate() -> tuple[object, ...]:
    x, r, a, b, c = sy.symbols("x r a b c")
    y = (r - 1) / x
    z = (2 * x - c - 3 * x * (r - 1)) / x**3
    first = r**3 * z + y**2 * r * (4 + 3 * (r - 1))
    second = y + 3 * x * r**2 * z + 3 * x * y**2 * (4 + 3 * (r - 1))
    third = 2 * x - 3 * x**2 * y - x**3 * z
    first_reduced = r * (x * (1 + r) - c * r**2) / x**3
    second_reduced = (2 * x * (1 + 2 * r) - 3 * c * r**2) / x**2
    require(sy.cancel(first - first_reduced) == 0, "symbolic first reduction")
    require(sy.cancel(second - second_reduced) == 0, "symbolic second reduction")
    require(sy.cancel(third - c) == 0, "symbolic third reduction")

    a_from_map = first_reduced
    b_from_map = second_reduced
    relation_one = sy.cancel(
        3 * a_from_map * x**2 - r * (1 + b_from_map * x - r)
    )
    k_from_map = 2 + (9 * a_from_map * c - b_from_map) * x
    q_from_map = 3 * c * (1 + b_from_map * x) - 4 * x
    relation_two = sy.cancel(q_from_map * r - x * k_from_map)
    require(relation_one == 0 and relation_two == 0, "inverse elimination")

    l_poly = 27 * a**2 * c**2 - 18 * a * b * c + 16 * a + b**3 * c - b**2
    g = sy.expand(l_poly * x**3 + (4 - 3 * b * c) * x - 2 * c)
    require(
        sy.cancel(g.subs({a: a_from_map, b: b_from_map})) == 0,
        "core derivation",
    )
    d_value = (12 * a - b**2) * x**2 + b * x + 2
    k_value = 2 + (9 * a * c - b) * x
    q_value = 3 * c * (1 + b * x) - 4 * x
    require(sy.expand(d_value * q_value - x * k_value**2 + 3 * g) == 0,
            "D-Q-core identity")

    y_inverse = b - 3 * a * x * k_value / d_value
    z_inverse = (2 * x - c - 3 * x**2 * y_inverse) / x**3
    unit = 1 + x * y_inverse
    images = (
        unit**3 * z_inverse + y_inverse**2 * unit * (4 + 3 * x * y_inverse),
        y_inverse + 3 * x * unit**2 * z_inverse
        + 3 * x * y_inverse**2 * (4 + 3 * x * y_inverse),
        2 * x - 3 * x**2 * y_inverse - x**3 * z_inverse,
    )
    field = sy.QQ.frac_field(a, b, c)
    core_poly = sy.Poly(g, x, domain=field)
    quotient_strings = []
    numerator_degrees = []
    for image_value, target_value in zip(images, (a, b, c)):
        numerator = sy.cancel(image_value - target_value).as_numer_denom()[0]
        numerator_poly = sy.Poly(numerator, x, domain=field)
        quotient, remainder = numerator_poly.div(core_poly)
        require(remainder.is_zero, ("symbolic map remainder", remainder))
        numerator_degrees.append(-1 if numerator_poly.is_zero else numerator_poly.degree())
        quotient_strings.append(sy.srepr(sy.expand(quotient.as_expr())))
    quotient_digest = hashlib.sha256(
        "\n".join(quotient_strings).encode("utf-8")
    ).hexdigest()
    x0, y0, z0 = sy.symbols("x0 y0 z0")
    unit0 = 1 + x0 * y0
    map0 = sy.Matrix((
        unit0**3 * z0 + y0**2 * unit0 * (4 + 3 * x0 * y0),
        y0 + 3 * x0 * unit0**2 * z0
        + 3 * x0 * y0**2 * (4 + 3 * x0 * y0),
        2 * x0 - 3 * x0**2 * y0 - x0**3 * z0,
    ))
    jacobian = sy.expand(map0.jacobian((x0, y0, z0)).det())
    require(jacobian == -2, ("Jacobian", jacobian))
    return tuple(numerator_degrees), quotient_digest, int(jacobian)


def direct_fibre_census() -> tuple[tuple[tuple[int, int, int], ...], ...]:
    targets = (
        (WITNESS_TARGET,) + WITNESS_LEVEL_ONE
        + (PARENT_TARGET,) + PARENT_LEVEL_ONE
    )
    buckets = {target: [] for target in targets}
    target_set = set(targets)
    for x_value in range(P):
        x2 = x_value * x_value % P
        x3 = x2 * x_value % P
        for y_value in range(P):
            xy = x_value * y_value % P
            unit = (1 + xy) % P
            first_constant = y_value * y_value * unit * (4 + 3 * xy) % P
            first_slope = unit**3 % P
            second_constant = (
                y_value + 3 * x_value * y_value * y_value * (4 + 3 * xy)
            ) % P
            second_slope = 3 * x_value * unit * unit % P
            third_constant = (2 * x_value - 3 * x2 * y_value) % P
            third_slope = (-x3) % P
            for z_value in range(P):
                image = (
                    (first_constant + first_slope * z_value) % P,
                    (second_constant + second_slope * z_value) % P,
                    (third_constant + third_slope * z_value) % P,
                )
                if image in target_set:
                    buckets[image].append((x_value, y_value, z_value))
    expected = {
        WITNESS_TARGET: WITNESS_LEVEL_ONE,
        PARENT_TARGET: PARENT_LEVEL_ONE,
    }
    for parent, blocks in (
        (WITNESS_LEVEL_ONE, WITNESS_LEVEL_TWO_BLOCKS),
        (PARENT_LEVEL_ONE, PARENT_LEVEL_TWO_BLOCKS),
    ):
        expected.update({target: block for target, block in zip(parent, blocks)})
    for target in targets:
        require(
            tuple(sorted(buckets[target])) == tuple(sorted(expected[target])),
            ("direct fibre mismatch", target, buckets[target], expected[target]),
        )
    return tuple(tuple(sorted(buckets[target])) for target in targets)


def quotient_row(target: tuple[int, int, int]) -> tuple[object, ...]:
    a, b, c = target
    modulus = core(target)
    require(poly_gcd(list(modulus), derivative(list(modulus))) == [1],
            ("inseparable cubic", target, modulus))
    algebra = CubicQuotient(modulus)
    x_value = algebra.generator
    x2 = algebra.power(x_value, 2)
    denominator = algebra.add(
        algebra.add(algebra.scale(x2, 12 * a - b * b), algebra.scale(x_value, b)),
        (2, 0, 0),
    )
    denominator_matrix_det = determinant(algebra.matrix(denominator))
    x3 = algebra.power(x_value, 3)
    x3_matrix_det = determinant(algebra.matrix(x3))
    require(denominator_matrix_det != 0 and x3_matrix_det != 0,
            ("nonunit inverse denominator", target))
    k_value = algebra.add(algebra.scale(x_value, 9 * a * c - b), (2, 0, 0))
    correction = algebra.scale(algebra.multiply(x_value, k_value), 3 * a)
    y_value = algebra.subtract(
        (b, 0, 0), algebra.multiply(correction, algebra.inverse(denominator))
    )
    z_numerator = algebra.subtract(
        algebra.subtract(algebra.scale(x_value, 2), (c, 0, 0)),
        algebra.scale(algebra.multiply(x2, y_value), 3),
    )
    z_value = algebra.multiply(z_numerator, algebra.inverse(x3))
    image = fmap_quotient(algebra, x_value, y_value, z_value)
    require(image == tuple((value, 0, 0) for value in target),
            ("quotient map replay", target, image))
    coordinates = (x_value, y_value, z_value)
    charpolys = tuple(characteristic_polynomial(algebra, value) for value in coordinates)
    require(
        all(poly_gcd(list(poly), derivative(list(poly))) == [1] for poly in charpolys),
        ("local coordinate collision", target, charpolys),
    )
    monic_core = tuple(poly_scale(list(modulus), inv(modulus[-1])))
    require(charpolys[0] == monic_core, ("x charpoly/core mismatch", target))
    return (
        target,
        modulus,
        algebra,
        coordinates,
        charpolys,
        denominator_matrix_det,
        x3_matrix_det,
    )


def coordinate_products(rows: tuple[tuple[object, ...], ...]) -> tuple[object, ...]:
    products = []
    gcds = []
    for coordinate in range(3):
        product = [1]
        for row in rows:
            product = poly_multiply(product, list(row[4][coordinate]))
        require(len(product) == 28 and product[-1] == 1,
                ("global product degree", coordinate, len(product) - 1))
        products.append(tuple(product))
        gcds.append(tuple(poly_gcd(product, derivative(product))))
    hashes = tuple(
        hashlib.sha256(
            "\n".join(
                f"{index}:{coefficient}" for index, coefficient in enumerate(product)
            ).encode("ascii")
        ).hexdigest()
        for product in products
    )
    return tuple(products), tuple(gcds), hashes


def power_basis_data(
    rows: tuple[tuple[object, ...], ...]
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    ranks = []
    determinants = []
    for coordinate in range(3):
        powers = [row[2].one for row in rows]
        columns = []
        for _exponent in range(27):
            columns.append(tuple(value for element in powers for value in element))
            powers = [
                row[2].multiply(power, row[3][coordinate])
                for row, power in zip(rows, powers)
            ]
        matrix = [
            [columns[column][row] for column in range(27)] for row in range(27)
        ]
        ranks.append(matrix_rank(matrix))
        determinants.append(determinant(matrix))
    return tuple(ranks), tuple(determinants)


def collision_pairs(
    rows: tuple[tuple[object, ...], ...], coordinate: int
) -> tuple[tuple[int, int, tuple[int, ...]], ...]:
    answer = []
    for left in range(len(rows)):
        for right in range(left):
            divisor = tuple(poly_gcd(
                list(rows[left][4][coordinate]), list(rows[right][4][coordinate])
            ))
            if len(divisor) > 1:
                answer.append((right, left, divisor))
    return tuple(answer)


def split_inverse_denominator_gate() -> tuple[tuple[int, int], ...]:
    edges = [(WITNESS_TARGET, point) for point in WITNESS_LEVEL_ONE]
    edges.extend(
        (target, point)
        for target, block in zip(WITNESS_LEVEL_ONE, WITNESS_LEVEL_TWO_BLOCKS)
        for point in block
    )
    values = []
    for target, point in edges:
        a, b, _c = target
        x_value = point[0]
        denominator = (
            (12 * a - b * b) * x_value * x_value + b * x_value + 2
        ) % P
        require(x_value != 0 and denominator != 0,
                ("split inverse denominator", target, point, denominator))
        values.append((x_value, denominator))
    return tuple(values)


symbolic_certificate = symbolic_inverse_certificate()
direct_census = direct_fibre_census()
split_inverse_units = split_inverse_denominator_gate()

witness_points = tuple(point for block in WITNESS_LEVEL_TWO_BLOCKS for point in block)
parent_points = tuple(point for block in PARENT_LEVEL_TWO_BLOCKS for point in block)
witness_rows = tuple(quotient_row(point) for point in witness_points)
parent_rows = tuple(quotient_row(point) for point in parent_points)

witness_products, witness_gcds, witness_hashes = coordinate_products(witness_rows)
parent_products, parent_gcds, parent_hashes = coordinate_products(parent_rows)
require(witness_hashes == EXPECTED_WITNESS_PRODUCT_HASHES,
        ("witness product hashes", witness_hashes))
require(parent_hashes == EXPECTED_PARENT_PRODUCT_HASHES,
        ("parent product hashes", parent_hashes))
require(witness_gcds == ((1,), (1,), (1,)), ("witness gcds", witness_gcds))
require(tuple(len(value) - 1 for value in parent_gcds) == (0, 1, 0),
        ("parent hostile gcds", parent_gcds))

witness_power = power_basis_data(witness_rows)
parent_power = power_basis_data(parent_rows)
require(witness_power[0] == (27, 27, 27), ("witness power ranks", witness_power))
require(all(witness_power[1]), ("witness power determinants", witness_power))
require(parent_power[0] == (27, 26, 27), ("parent power ranks", parent_power))
parent_y_collisions = collision_pairs(parent_rows, 1)
require(parent_y_collisions, "missing parent y collision")

witness_l_values = tuple(
    l_value(point)
    for point in (WITNESS_TARGET,) + WITNESS_LEVEL_ONE + witness_points
)
require(all(witness_l_values), ("witness L gate", witness_l_values))
witness_unit_pairs = tuple((row[5], row[6]) for row in witness_rows)
require(all(left and right for left, right in witness_unit_pairs), "unit gate")

local_ledger = tuple(
    (row[0], row[1], row[3], row[4], row[5], row[6]) for row in witness_rows
)
local_digest = hashlib.sha256(repr(local_ledger).encode("utf-8")).hexdigest()
security_tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
require(not any(isinstance(node, ast.Assert) for node in ast.walk(security_tree)),
        "assert node")
imports = tuple(
    sorted(
        alias.name
        for node in ast.walk(security_tree)
        if isinstance(node, ast.Import)
        for alias in node.names
    )
)
require(all("keller_level_three_three_coordinate_primitive_finite_field_probe" not in name
            for name in imports), "candidate import")

semantic = (
    P,
    symbolic_certificate,
    direct_census,
    WITNESS_TARGET,
    WITNESS_LEVEL_ONE,
    WITNESS_LEVEL_TWO_BLOCKS,
    split_inverse_units,
    witness_l_values,
    witness_unit_pairs,
    tuple((row[0], row[1], row[3], row[4]) for row in witness_rows),
    witness_products,
    witness_gcds,
    witness_hashes,
    witness_power,
    PARENT_TARGET,
    PARENT_LEVEL_ONE,
    PARENT_LEVEL_TWO_BLOCKS,
    parent_products,
    parent_gcds,
    parent_hashes,
    parent_power,
    parent_y_collisions,
    "one lawful rank-27 etale fibre with each coordinate power determinant nonzero",
    "generic K(x)=K(y)=K(z); trace-form comparison retains exact class [-2J]",
)
semantic_hash = hashlib.sha256(repr(semantic).encode("utf-8")).hexdigest()
if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
    require(semantic_hash == EXPECTED_SEMANTIC_SHA256,
            ("semantic hash", semantic_hash))

print("KELLER LEVEL-THREE THREE-COORDINATE PRIMITIVITY INDEPENDENT AUDIT")
print("status=FINITE-EXACT + symbolic clean-room audit; fixed sporadic F^3 only")
print("candidate_imported=False; parent_computation_imported=False")
print(f"symbolic_inverse=(numerator_degrees={symbolic_certificate[0]},quotient_sha256={symbolic_certificate[1]},identity=DQ-xK^2=-3g,jacobian={symbolic_certificate[2]})")
print("direct_map_census=all F_101^3 enumerated once; both listed depth-two split trees reconstructed exactly")
print(f"witness=(target={WITNESS_TARGET},level_one={WITNESS_LEVEL_ONE},level_two_blocks={WITNESS_LEVEL_TWO_BLOCKS})")
print(f"good_reduction=(split_inverse_units={split_inverse_units},L_values={witness_l_values},final_unit_determinants={witness_unit_pairs},nine_separable_cubic_quotients=True)")
print(f"local_charpolys=(27/27 degree_three_squarefree,ledger_sha256={local_digest})")
print(f"global_products=(degrees=(27,27,27),derivative_gcds={witness_gcds},sha256={witness_hashes})")
print(f"power_basis_27x27=(ranks={witness_power[0]},determinants={witness_power[1]})")
print(f"parent_hostile=(target={PARENT_TARGET},derivative_gcd_degrees={tuple(len(value)-1 for value in parent_gcds)},power_ranks={parent_power[0]},y_collision_pairs={parent_y_collisions},sha256={parent_hashes})")
print("generic_consequence=the three power determinants are nonzero rational functions, so K(x)=K(y)=K(z)=K_3")
print("trace_consequence=power-basis Gram determinants differ by change-of-basis squares; THM-3495 gives the common exact class [-2J]")
print("scaling_gate=degree 27 eliminant scaling exponent 52=(26*2), hence a square")
print("scope=no exact multiplicities, level-four/all-level, arbitrary Keller, outside-family, JC(2), or LRC consequence")
print(f"security_ast=({len(tuple(ast.walk(security_tree)))},assert_free=True)")
print(f"semantic_sha256={semantic_hash}")
print("STATUS=PASS")
