#!/usr/bin/env python3
"""Exact companion for THM-3253's positive owner-mass cyclicity theorem."""

import ast
from fractions import Fraction
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
    """Return (g^2,g,1) coefficients from THM-3246's table."""
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


def cleared_matrix(multiplier, shift, dilation):
    matrix = [[0] * 13 for _ in range(13)]
    for owner in range(168):
        x, y = points[(shift + multiplier * owner) % 168]
        matrix[x][y] = numerator(owner, dilation)
    return matrix


@lru_cache(maxsize=None)
def determinant(multiplier, shift, dilation):
    return determinant_bareiss(cleared_matrix(multiplier, shift, dilation))


def forward_coefficients(multiplier, shift, base):
    values = [determinant(multiplier, shift, dilation)
              for dilation in range(base, base + 27)]
    coefficients = []
    while values:
        coefficients.append(values[0])
        values = [values[index + 1] - values[index]
                  for index in range(len(values) - 1)]
    return tuple(coefficients)


def newton_value(coefficients, offset):
    return sum(value * comb(offset, degree)
               for degree, value in enumerate(coefficients))


# Determinants have degree at most 2*13=26.  Five phase representatives have
# one-sign Newton expansions from g=1.  The two exceptional representatives
# have finite negative heads and one-sign positive tails.
regular_coefficients = tuple(
    forward_coefficients(1, shift, 1) for shift in range(5)
)
require(all(all(value < 0 for value in coefficients)
            for coefficients in regular_coefficients),
        "regular Newton signs")

phase_five_head = tuple(determinant(1, 5, dilation)
                        for dilation in range(1, 18))
phase_five_tail = forward_coefficients(1, 5, 18)
require(all(value < 0 for value in phase_five_head), "phase-five head")
require(all(value > 0 for value in phase_five_tail), "phase-five tail")

phase_six_head = determinant(1, 6, 1)
phase_six_tail = forward_coefficients(1, 6, 2)
require(phase_six_head < 0, "phase-six head")
require(all(value > 0 for value in phase_six_tail), "phase-six tail")

# Independent extrapolation controls beyond every interpolation window.
for shift, coefficients in enumerate(regular_coefficients):
    require(newton_value(coefficients, 79)
            == determinant(1, shift, 80), ("regular extrapolation", shift))
require(newton_value(phase_five_tail, 62) == determinant(1, 5, 80),
        "phase-five extrapolation")
require(newton_value(phase_six_tail, 78) == determinant(1, 6, 80),
        "phase-six extrapolation")

# alpha^7=9u swaps the two affine axes, with row multiplier 5 and column
# multiplier 9.  Their permutation signs are -1 and +1, respectively.
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
for shift in range(7):
    for dilation in range(1, 28):
        require(determinant(1, shift + 7, dilation)
                == -determinant(1, shift, dilation),
                ("phase-seven determinant symmetry", shift, dilation))

# alpha^14=6 is scalar, so all 12 radial translates give simultaneous row
# and column permutations.  Frobenius is (x,y)->(x,-y), and owner reflection
# turns multiplier -a into +a with a phase correction.  These identities
# first expand the fourteen Newton-certified phases to 672 gauges.
restricted_multipliers = (1, 13, 155, 167)  # +/-1,+/-13 modulo 168
restricted_gauge_count = len(restricted_multipliers) * 168
require(restricted_gauge_count == 672, "restricted gauge census")

# Direct finite controls at two head and two tail dilations guard the index
# conventions used in the algebraic scalar/Frobenius/reflection reductions.
for multiplier in restricted_multipliers:
    for shift in range(168):
        for dilation in (1, 2, 18, 80):
            require(determinant(multiplier, shift, dilation) != 0,
                    ("restricted gauge direct control",
                     multiplier, shift, dilation))


# The remaining primitive Singer generators need not share the one-sign
# Newton tails above.  Reduce all 48 unit multipliers by the exact
# {+/-1,+/-13} symmetry and all phases by the scalar/anti-diagonal symmetries.
# This leaves 12*7=84 determinant polynomials.  We recover each polynomial in
# the ordinary monomial basis over Z, factor its exact zero at g=0, and find a
# prime for which the remaining factor has no root.  An integer zero of the
# remaining factor would reduce to a root in every finite field, so one such
# prime is an all-integer certificate.
unit_multipliers = tuple(value for value in range(168)
                         if gcd(value, 168) == 1)
symmetry_multipliers = (1, 13, 155, 167)
multiplier_orbits = []
unseen = set(unit_multipliers)
while unseen:
    representative = min(unseen)
    orbit = {representative * symmetry % 168
             for symmetry in symmetry_multipliers}
    require(len(orbit) == 4 and orbit <= unseen,
            ("multiplier orbit", representative, orbit, unseen & orbit))
    multiplier_orbits.append(tuple(sorted(orbit)))
    unseen -= orbit
multiplier_representatives = tuple(orbit[0]
                                   for orbit in multiplier_orbits)
require(multiplier_representatives ==
        (1, 5, 11, 17, 19, 23, 29, 31, 43, 47, 59, 71),
        "multiplier quotient representatives")
require(len(unit_multipliers) == 48 and len(multiplier_orbits) == 12,
        "unit multiplier quotient census")


def monomial_coefficients(values):
    """Interpolate a degree<=26 integer polynomial exactly from g=0..26."""
    rows = list(values)
    forward = []
    while rows:
        forward.append(rows[0])
        rows = [rows[index + 1] - rows[index]
                for index in range(len(rows) - 1)]

    polynomial = [Fraction(0) for _ in range(len(values))]
    binomial_basis = [Fraction(1)]
    for degree, coefficient in enumerate(forward):
        for exponent, value in enumerate(binomial_basis):
            polynomial[exponent] += coefficient * value
        if degree + 1 < len(values):
            next_basis = [Fraction(0)
                          for _ in range(len(binomial_basis) + 1)]
            for exponent, value in enumerate(binomial_basis):
                next_basis[exponent] -= value * degree / (degree + 1)
                next_basis[exponent + 1] += value / (degree + 1)
            binomial_basis = next_basis
    require(all(value.denominator == 1 for value in polynomial),
            "integer monomial coefficients")
    return tuple(value.numerator for value in polynomial)


def evaluate_mod(coefficients, value, prime):
    result = 0
    for coefficient in reversed(coefficients):
        result = (result * value + coefficient) % prime
    return result


def monomial_coefficients_mod(values, prime):
    """Independent Newton interpolation directly in F_prime."""
    rows = [value % prime for value in values]
    forward = []
    while rows:
        forward.append(rows[0])
        rows = [(rows[index + 1] - rows[index]) % prime
                for index in range(len(rows) - 1)]

    polynomial = [0 for _ in range(len(values))]
    binomial_basis = [1]
    for degree, coefficient in enumerate(forward):
        for exponent, value in enumerate(binomial_basis):
            polynomial[exponent] = (
                polynomial[exponent] + coefficient * value) % prime
        if degree + 1 < len(values):
            inverse = pow(degree + 1, -1, prime)
            next_basis = [0 for _ in range(len(binomial_basis) + 1)]
            for exponent, value in enumerate(binomial_basis):
                next_basis[exponent] = (
                    next_basis[exponent] - degree * value * inverse) % prime
                next_basis[exponent + 1] = (
                    next_basis[exponent + 1] + value * inverse) % prime
            binomial_basis = next_basis
    return tuple(polynomial)


certificate_primes = (29, 31, 37, 41, 43, 47, 53, 59, 61)
certificate_table = []
exact_polynomials = {}
order_counts = {}
prime_counts = {}
for multiplier in multiplier_representatives:
    for shift in range(7):
        values = tuple(determinant(multiplier, shift, dilation)
                       for dilation in range(27))
        coefficients = monomial_coefficients(values)
        exact_polynomials[(multiplier, shift)] = coefficients
        order = next((exponent for exponent, coefficient
                      in enumerate(coefficients) if coefficient), None)
        require(order is not None and order > 0,
                ("exact positive vanishing order", multiplier, shift, order))
        require(sum(coefficient * (80 ** exponent)
                    for exponent, coefficient in enumerate(coefficients))
                == determinant(multiplier, shift, 80),
                ("all-gauge polynomial extrapolation", multiplier, shift))
        quotient = coefficients[order:]
        certificate = next(
            (prime for prime in certificate_primes
             if all(evaluate_mod(quotient, value, prime) != 0
                    for value in range(prime))), None)
        require(certificate is not None,
                ("missing root-free prime", multiplier, shift, order))
        direct_modular = monomial_coefficients_mod(values, certificate)
        require(direct_modular == tuple(coefficient % certificate
                                        for coefficient in coefficients),
                ("exact/modular interpolation comparison",
                 multiplier, shift, certificate))
        certificate_table.append((multiplier, shift, certificate, order))
        order_counts[order] = order_counts.get(order, 0) + 1
        prime_counts[certificate] = prime_counts.get(certificate, 0) + 1

require(len(certificate_table) == 84, "modular certificate census")
require(order_counts == {4: 2, 5: 18, 6: 47, 7: 15, 8: 2},
        ("exact g-adic order census", order_counts))
require(prime_counts == {29: 34, 31: 21, 37: 9, 41: 5,
                         43: 7, 47: 4, 53: 2, 59: 1, 61: 1},
        ("root-free prime cover census", prime_counts))
certificate_digest = sha256("\n".join(
    ",".join(map(str, row)) for row in certificate_table
).encode("ascii")).hexdigest()
require(certificate_digest ==
        "34771ce3677e1fa9cb1932324b54aa743e6e9b97304805e65459071d3b9fc132",
        "modular certificate table digest")
exceptional_lead = exact_polynomials[(17, 3)][7]
require(exceptional_lead == -266920406286336,
        "exceptional exact seventh coefficient")
require(exceptional_lead % 31 == 0 and exceptional_lead % 59 == 44,
        "exceptional modular-order hostile")


def digest_nested(rows):
    return sha256("\n".join(
        ",".join(map(str, row)) for row in rows
    ).encode("ascii")).hexdigest()


regular_digest = digest_nested(regular_coefficients)
phase_five_head_digest = sha256(
    "\n".join(map(str, phase_five_head)).encode("ascii")).hexdigest()
phase_five_tail_digest = sha256(
    "\n".join(map(str, phase_five_tail)).encode("ascii")).hexdigest()
phase_six_tail_digest = sha256(
    "\n".join(map(str, phase_six_tail)).encode("ascii")).hexdigest()
require(regular_digest ==
        "1aeb6d4070908447584ff0fee52c2ab7e7f0d2287f766bf69b962ef6f2815e16",
        "regular Newton digest")
require(phase_five_head_digest ==
        "b228dcdcca8a65cc11b61894e5cf07921238596ff03680873d88fe427d444274",
        "phase-five head digest")
require(phase_five_tail_digest ==
        "f7a36db8a325e9bc576343143455035dab7e8bc1b16341159406c7427b0d52af",
        "phase-five tail digest")
require(phase_six_tail_digest ==
        "5fb6da1e11014bec9e2e97daa588f982ff318d6a437ebfc08da15141c2da40b4",
        "phase-six tail digest")

print("THM-3253 POSITIVE OWNER-MASS ALL-GAUGE CYCLICITY EXACT AUDIT")
print("dependency_hash_checks=%d" % len(DEPENDENCIES))
print("assert_nodes=%d,float_literals=%d" % (assert_nodes, float_literals))
print("owner_numerator_polynomials=168,degree=2,all_positive_for_g>=1")
print("singer_plane=F13[u]/(u^2-2),alpha=(1,2),orbit=168")
print("regular_phase_classes=0..4,base=1,newton_coefficients=5*27,all_negative")
print("regular_newton_digest=%s" % regular_digest)
print("phase5=head_1_to_17_negative,tail_base18_27_coefficients_positive")
print("phase5_head_digest=%s" % phase_five_head_digest)
print("phase5_tail_digest=%s" % phase_five_tail_digest)
print("phase6=head_g1_negative,tail_base2_27_coefficients_positive")
print("phase6_tail_digest=%s" % phase_six_tail_digest)
print("phase7_symmetry=D_(b+7)=-D_b,scalar_phase_period=14")
print("restricted_multipliers=(1,13,155,167),gauge_count=672")
print("direct_restricted_gauge_controls=672*4=%d" %
      (restricted_gauge_count * 4))
print("unit_multiplier_orbits=48/4=12,phase_classes=7,modular_polynomials=84")
print("exact_g_orders=%s" % sorted(order_counts.items()))
print("root_free_prime_cover=%s" % sorted(prime_counts.items()))
print("modular_certificate_digest=%s" % certificate_digest)
print("exceptional_(a,b)=(17,3),exact_r=7,c7=-266920406286336,c7_mod31=0,c7_mod59=44")
print("all_primitive_singer_gauges=48*168=8064")
print("all_integer_dilations_nonsingular=PASS")
print("nonnegative_delta0_packet_orbit_span=12*169+13=2041")
print("scope=abstract_positive-owner-mass-relocation-not-canonical-endpoint-current")
print("all_exact_checks=PASS")
