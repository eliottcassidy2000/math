#!/usr/bin/env python3
"""Exact companion for THM-3253's universal positive owner-mass cyclicity."""

import ast
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from math import comb, gcd
from pathlib import Path


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = {
    ROOT / "01-canon/theorems/THM-3234-singer-owner-compactification-and-pointed-heisenberg-carrier-gate.md":
        "ef77a1f8fce16eb851eb38d5110a61ab73aa693f2d0ee9e11a912aa4fc302c87",
    ROOT / "01-canon/theorems/THM-3246-all-dilation-second-owner-seam-stabilization-and-sign-word.md":
        "6badc0c9aba09b56d3d055a96cb8ef8b619d8492508bf21476eba5f624b13055",
    ROOT / "04-computation/lrc_second_owner_all_dilation_seam_thm3246.py":
        "e23b098b38aa2199a348f48f8ab4ac0ce5913c870ead972bd31296494fc25a4b",
    ROOT / "05-knowledge/results/lrc_second_owner_all_dilation_seam_thm3246.out":
        "d7f7dd96b01c597113e78f903cad36246cb47b10e9a1758cb831aa0e83e8cebc",
    ROOT / "01-canon/theorems/THM-3250-charged-heisenberg-blowup-address-intertwiner-and-pointed-multiplicity-gate.md":
        "7e91a07e38d6869e6621b9057594cdf4745827a98de9618681d09427102b27ea",
}


def lf_bytes(path):
    return path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")


for dependency, expected in DEPENDENCIES.items():
    require(sha256(lf_bytes(dependency)).hexdigest() == expected,
            ("dependency hash drift", dependency.name))

tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
float_literals = sum(
    isinstance(node, ast.Constant) and isinstance(node.value, float)
    for node in ast.walk(tree)
)
require(assert_nodes == 0, "optimization-sensitive assert")
require(float_literals == 0, "floating literal")


def numerator_polynomial(cell):
    """Return (g^2,g,1) coefficients from THM-3246's exact table."""
    if cell <= 5 or cell >= 162:
        return 12096, -1032, 2
    if 6 <= cell <= 23 or 144 <= cell <= 161:
        return 12096, -24, 0
    if 24 <= cell <= 71:
        return 16044 - 168 * cell, 48, 0
    if 72 <= cell <= 95:
        return 4032, 96, 0
    if 96 <= cell <= 143:
        return 168 * cell - 12012, 48, 0
    raise RuntimeError(("cell outside table", cell))


polynomials = tuple(numerator_polynomial(cell) for cell in range(168))
require(all(polynomials[167 - cell] == polynomials[cell]
            for cell in range(168)), "owner reflection")
require(all(a > 0 and a + b + c > 0 and 2 * a + b > 0
            for a, b, c in polynomials), "positive increasing numerators")


def numerator(cell, dilation):
    a, b, c = polynomials[cell]
    return a * dilation * dilation + b * dilation + c


def field_mul(x, y):
    a, b = x
    c, d = y
    return ((a * c + 2 * b * d) % 13, (a * d + b * c) % 13)


alpha = (1, 2)
points = []
point = (1, 0)
for _ in range(168):
    points.append(point)
    point = field_mul(point, alpha)
require(point == (1, 0) and len(set(points)) == 168, "Singer orbit")
require(points[7] == (0, 9) and points[14] == (6, 0),
        "anti-diagonal and scalar phase steps")
require(all(points[(index + 7) % 168]
            == (5 * points[index][1] % 13, 9 * points[index][0] % 13)
            for index in range(168)), "anti-diagonal phase action")
require(all(points[(index + 14) % 168]
            == (6 * points[index][0] % 13, 6 * points[index][1] % 13)
            for index in range(168)), "scalar phase action")
require(all(points[13 * index % 168]
            == (points[index][0], -points[index][1] % 13)
            for index in range(168)), "Frobenius exponent action")


def permutation_sign(permutation):
    inversions = sum(permutation[left] > permutation[right]
                     for left in range(len(permutation))
                     for right in range(left + 1, len(permutation)))
    return -1 if inversions % 2 else 1


require(permutation_sign(tuple(5 * value % 13 for value in range(13))) == -1,
        "row anti-diagonal sign")
require(permutation_sign(tuple(9 * value % 13 for value in range(13))) == 1,
        "column anti-diagonal sign")
require(permutation_sign(tuple(6 * value % 13 for value in range(13))) ** 2
        == 1, "scalar row-column sign")
require(permutation_sign(tuple(-value % 13 for value in range(13))) == 1,
        "Frobenius column sign")


def cleared_matrix(multiplier, shift, dilation):
    matrix = [[0] * 13 for _ in range(13)]
    for owner in range(168):
        x, y = points[(shift + multiplier * owner) % 168]
        matrix[x][y] = numerator(owner, dilation)
    return matrix


def determinant_bareiss(matrix):
    a = [row[:] for row in matrix]
    previous = 1
    sign = 1
    for column in range(len(a) - 1):
        pivot_row = next(
            (row for row in range(column, len(a)) if a[row][column]), None)
        if pivot_row is None:
            return 0
        if pivot_row != column:
            a[column], a[pivot_row] = a[pivot_row], a[column]
            sign = -sign
        pivot = a[column][column]
        for row in range(column + 1, len(a)):
            for index in range(column + 1, len(a)):
                value = (a[row][index] * pivot
                         - a[row][column] * a[column][index])
                require(value % previous == 0,
                        ("Bareiss divisibility", column, row, index))
                a[row][index] = value // previous
        previous = pivot
        for row in range(column + 1, len(a)):
            a[row][column] = 0
    return sign * a[-1][-1]


@lru_cache(maxsize=None)
def determinant(multiplier, shift, dilation):
    return determinant_bareiss(cleared_matrix(multiplier, shift, dilation))


def forward_coefficients(values):
    values = list(values)
    coefficients = []
    while values:
        coefficients.append(values[0])
        values = [values[index + 1] - values[index]
                  for index in range(len(values) - 1)]
    return tuple(coefficients)


def newton_value(coefficients, offset):
    return sum(value * comb(offset, degree)
               for degree, value in enumerate(coefficients))


# Multiplication of a primitive exponent by -1 is owner reflection, and
# multiplication by 13 is Frobenius.  The 48 primitive multipliers therefore
# split into twelve four-element orbits.  The phase symmetries reduce 168
# shifts to seven determinant classes up to sign.
units = tuple(multiplier for multiplier in range(168)
              if gcd(multiplier, 168) == 1)
seen = set()
multiplier_representatives = []
for multiplier in units:
    if multiplier in seen:
        continue
    orbit = {(factor * multiplier) % 168
             for factor in (1, -1, 13, -13)}
    require(len(orbit) == 4, ("multiplier orbit size", multiplier))
    seen.update(orbit)
    multiplier_representatives.append(min(orbit))
multiplier_representatives = tuple(multiplier_representatives)
require(seen == set(units), "primitive multiplier orbit cover")
require(multiplier_representatives
        == (1, 5, 11, 17, 19, 23, 29, 31, 43, 47, 59, 71),
        "primitive multiplier representatives")

for multiplier in multiplier_representatives:
    for shift in range(7):
        for dilation in (1, 2, 57, 58):
            require(determinant(multiplier, shift + 7, dilation)
                    == -determinant(multiplier, shift, dilation),
                    ("phase-seven determinant symmetry",
                     multiplier, shift, dilation))


# Every representative determinant has degree at most 26.  A uniform base 58
# works for all 84 classes.  The finite head is essential: sixteen classes
# change sign before that base, including one between g=57 and g=58.
head_rows = []
newton_rows = []
tail_signs = Counter()
transition_counts = Counter()
base_one_same_sign = 0
minimum_signed_coefficient = None
maximum_signed_coefficient = None
representative_class_count = 0
for multiplier in multiplier_representatives:
    for shift in range(7):
        values = tuple(determinant(multiplier, shift, dilation)
                       for dilation in range(1, 101))
        head = values[:57]
        require(all(head), ("zero finite-head determinant", multiplier, shift))

        base_one = forward_coefficients(values[:27])
        base_one_sign = 1 if base_one[0] > 0 else -1
        if all(base_one_sign * value > 0 for value in base_one):
            base_one_same_sign += 1

        coefficients = forward_coefficients(values[57:84])
        require(len(coefficients) == 27 and coefficients[-1] != 0,
                ("determinant degree", multiplier, shift))
        sign = 1 if coefficients[0] > 0 else -1
        require(all(sign * value > 0 for value in coefficients),
                ("mixed Newton tail", multiplier, shift))
        require(newton_value(coefficients, 27) == values[84],
                ("first Newton extrapolation", multiplier, shift))
        require(newton_value(coefficients, 42) == values[99],
                ("distant Newton extrapolation", multiplier, shift))

        for dilation in range(1, 58):
            if values[dilation - 1] * values[dilation] < 0:
                transition_counts[dilation] += 1

        signed = tuple(sign * value for value in coefficients)
        local_minimum = min(signed)
        local_maximum = max(signed)
        if minimum_signed_coefficient is None:
            minimum_signed_coefficient = local_minimum
            maximum_signed_coefficient = local_maximum
        else:
            minimum_signed_coefficient = min(
                minimum_signed_coefficient, local_minimum)
            maximum_signed_coefficient = max(
                maximum_signed_coefficient, local_maximum)

        tail_signs[sign] += 1
        head_rows.append("%d,%d|%s" %
                         (multiplier, shift, ",".join(map(str, head))))
        newton_rows.append("%d,%d|%s" %
                           (multiplier, shift,
                            ",".join(map(str, coefficients))))
        representative_class_count += 1

require(representative_class_count == 12 * 7 == 84,
        "representative class count")
require(base_one_same_sign == 68, "base-one Newton census")
require(tail_signs == Counter({1: 44, -1: 40}), "tail sign census")
require(transition_counts == Counter(
    {1: 4, 2: 4, 3: 1, 6: 1, 7: 1, 9: 2, 17: 2, 57: 1}),
    "finite-head sign-transition census")
require(minimum_signed_coefficient ==
        31583339331193330095097135851612043896766825662862153422608271409152000000,
        "minimum signed Newton coefficient")

head_digest = sha256("\n".join(head_rows).encode("ascii")).hexdigest()
newton_digest = sha256("\n".join(newton_rows).encode("ascii")).hexdigest()
require(head_digest ==
        "97260646b6d140268a649f6f5d4e3a8c0d6af6849945e81b172a2170e310a63f",
        "finite-head digest")
require(newton_digest ==
        "31fe23e30c06e0ddde4bb19f19c5d44e5795f00d0456e6a5e71c846d88a9c008",
        "Newton-tail digest")


# Independent modular elimination over every one of the 8,064 gauges at two
# primes and four hostile dilations.  The 168 values, each repeated 48 times,
# are the expected signed phase/reflection/Frobenius/scalar orbit census.
def determinant_mod(matrix, prime):
    a = [[entry % prime for entry in row] for row in matrix]
    result = 1
    for column in range(len(a)):
        pivot_row = next(
            (row for row in range(column, len(a)) if a[row][column]), None)
        if pivot_row is None:
            return 0
        if pivot_row != column:
            a[column], a[pivot_row] = a[pivot_row], a[column]
            result = -result
        pivot = a[column][column]
        result = result * pivot % prime
        inverse = pow(pivot, prime - 2, prime)
        for row in range(column + 1, len(a)):
            multiplier = a[row][column] * inverse % prime
            for index in range(column + 1, len(a)):
                a[row][index] = (
                    a[row][index] - multiplier * a[column][index]) % prime
    return result % prime


modular_rows = []
modular_controls = []
for prime in (1000000007, 1000000009):
    for dilation in (1, 57, 58, 1000):
        values = tuple(
            determinant_mod(cleared_matrix(multiplier, shift, dilation), prime)
            for multiplier in units
            for shift in range(168)
        )
        frequencies = Counter(values)
        require(len(values) == 8064 and 0 not in frequencies,
                ("singular full-gauge modular control", prime, dilation))
        require(len(frequencies) == 168
                and set(frequencies.values()) == {48},
                ("full-gauge determinant orbit census", prime, dilation))
        require(sum(values) % prime == 0,
                ("signed phase cancellation", prime, dilation))
        if dilation in (1, 57, 58):
            for multiplier in multiplier_representatives:
                for shift in range(7):
                    exact = determinant(multiplier, shift, dilation) % prime
                    modular = determinant_mod(
                        cleared_matrix(multiplier, shift, dilation), prime)
                    require(exact == modular,
                            ("exact/modular determinant mismatch",
                             prime, dilation, multiplier, shift))
        modular_rows.extend(map(str, values))
        modular_controls.append(
            (prime, dilation, len(frequencies), min(values), max(values)))

modular_digest = sha256("\n".join(modular_rows).encode("ascii")).hexdigest()
require(modular_digest ==
        "4c7cb86935def79998900ee220ca758ef9788d2d5b4d16b420995802f60ff12d",
        "full-gauge modular digest")


charged_dimension = 12 * 169
neutral_common_dimension = 13
common_dimension = charged_dimension + neutral_common_dimension
require((charged_dimension, neutral_common_dimension, common_dimension)
        == (2028, 13, 2041), "common-module dimensions")
require(common_dimension == 13 ** 3 - 13 ** 2 + 13,
        "THM-3250 maximal common rank")

print("THM-3253 UNIVERSAL POSITIVE OWNER-MASS NEWTON CYCLICITY EXACT AUDIT")
print("dependency_hash_checks=%d" % len(DEPENDENCIES))
print("assert_nodes=%d,float_literals=%d" % (assert_nodes, float_literals))
print("owner_numerator_polynomials=168,degree=2,all_positive_for_g>=1")
print("singer_plane=F13[u]/(u^2-2),alpha=(1,2),orbit=168")
print("primitive_multiplier_representatives=%s" %
      (multiplier_representatives,))
print("phase_classes=7,multiplier_classes=12,determinant_classes=84")
print("finite_head=g1..57,nonzero_checks=84*57=%d,base1_one_sign=68" %
      (84 * 57))
print("head_sign_transitions=%s,base58_is_sharp=PASS" %
      (tuple(sorted(transition_counts.items())),))
print("finite_head_digest=%s" % head_digest)
print("newton_base=58,degree26_classes=84,strict_one_sign_coefficients=84*27")
print("newton_tail_signs=%s,minimum_signed_coefficient=%d" %
      ((tail_signs[1], tail_signs[-1]), minimum_signed_coefficient))
print("newton_tail_digest=%s" % newton_digest)
print("full_gauge_count=8064,modular_controls=2*4*8064=%d" %
      (2 * 4 * 8064))
print("modular_control_summaries=%s" % (tuple(modular_controls),))
print("full_gauge_modular_digest=%s" % modular_digest)
print("all_integer_dilations_all_primitive_gauges_nonsingular=PASS")
print("nonnegative_delta0_packet_orbit_span=12*169+13=2041")
print("scope=abstract-positive-owner-mass-relocation-not-canonical-endpoint-current")
print("all_exact_checks=PASS")
