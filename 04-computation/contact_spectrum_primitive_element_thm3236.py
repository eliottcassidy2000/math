"""Exact companion for THM-3236's contact-spectrum primitive-element gate."""

from fractions import Fraction
from itertools import permutations
from math import comb
import ast
import hashlib
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCY = ROOT / (
    "01-canon/theorems/"
    "THM-3232-root-free-contact-stratum-norm-and-discriminant-power.md"
)
DEPENDENCY_SHA256 = "dcf42b70cad1bc940aaeea13d7d0a4b981678b8d6ca87e828b06624c242d77a3"


def lf_sha256(path):
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


require(lf_sha256(DEPENDENCY) == DEPENDENCY_SHA256, "THM-3232 dependency hash")
syntax_tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax_tree))
float_literals = sum(
    isinstance(node, ast.Constant) and isinstance(node.value, float)
    for node in ast.walk(syntax_tree)
)
require(assert_nodes == 0, "assert statements are optimization-sensitive")
require(float_literals == 0, "floating literals are forbidden")


def trim(poly):
    out = [Fraction(value) for value in poly]
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def add(left, right):
    degree = max(len(left), len(right))
    out = [Fraction(0)] * degree
    for index in range(degree):
        out[index] = (left[index] if index < len(left) else 0) + (
            right[index] if index < len(right) else 0
        )
    return trim(out)


def sub(left, right):
    return add(left, [-Fraction(value) for value in right])


def scale(poly, scalar):
    return trim([Fraction(scalar) * value for value in poly])


def mul(left, right):
    out = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, left_i in enumerate(left):
        for j, right_j in enumerate(right):
            out[i + j] += left_i * right_j
    return trim(out)


def power(poly, exponent):
    out = [Fraction(1)]
    for _ in range(exponent):
        out = mul(out, poly)
    return out


def monic(poly):
    poly = trim(poly)
    return scale(poly, Fraction(1, 1) / poly[-1])


def derivative(poly):
    if len(poly) <= 1:
        return [Fraction(0)]
    return trim([index * poly[index] for index in range(1, len(poly))])


def divmod_poly(dividend, divisor):
    dividend = trim(dividend)
    divisor = trim(divisor)
    require(divisor != [0], "division by zero polynomial")
    if len(dividend) < len(divisor):
        return [Fraction(0)], dividend
    quotient = [Fraction(0)] * (len(dividend) - len(divisor) + 1)
    remainder = list(dividend)
    while remainder != [0] and len(remainder) >= len(divisor):
        shift = len(remainder) - len(divisor)
        coefficient = remainder[-1] / divisor[-1]
        quotient[shift] += coefficient
        remainder = sub(
            remainder,
            [Fraction(0)] * shift + [coefficient * value for value in divisor],
        )
    return trim(quotient), trim(remainder)


def remainder(poly, modulus):
    return divmod_poly(poly, modulus)[1]


def quotient_exact(dividend, divisor):
    quotient, rem = divmod_poly(dividend, divisor)
    require(rem == [0], "expected exact division")
    return quotient


def gcd_poly(left, right):
    left = trim(left)
    right = trim(right)
    while right != [0]:
        left, right = right, remainder(left, right)
    return monic(left)


def inverse_mod(poly, modulus):
    old_r, new_r = trim(poly), trim(modulus)
    old_s, new_s = [Fraction(1)], [Fraction(0)]
    while new_r != [0]:
        quotient, next_r = divmod_poly(old_r, new_r)
        old_r, new_r = new_r, next_r
        old_s, new_s = new_s, sub(old_s, mul(quotient, new_s))
    require(len(old_r) == 1 and old_r[0] != 0, "class is not invertible")
    return remainder(scale(old_s, Fraction(1, 1) / old_r[0]), modulus)


def hasse_poly(poly, order):
    if order >= len(poly):
        return [Fraction(0)]
    return trim(
        [
            Fraction(poly[index]) * comb(index, order)
            for index in range(order, len(poly))
        ]
    )


def omega(poly_f, poly_g, order):
    return sub(
        mul(hasse_poly(poly_f, order), derivative(poly_g)),
        mul(hasse_poly(poly_g, order), derivative(poly_f)),
    )


def determinant(matrix):
    size = len(matrix)
    work = [[Fraction(value) for value in row] for row in matrix]
    require(all(len(row) == size for row in work), "square matrix required")
    value = Fraction(1)
    sign = 1
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if work[row][column] != 0),
            None,
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            sign = -sign
        pivot_value = work[column][column]
        value *= pivot_value
        for row in range(column + 1, size):
            factor = work[row][column] / pivot_value
            for entry in range(column, size):
                work[row][entry] -= factor * work[column][entry]
    return sign * value


def solve_linear(matrix, target):
    size = len(matrix)
    work = [
        [Fraction(value) for value in row] + [Fraction(target[index])]
        for index, row in enumerate(matrix)
    ]
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if work[row][column] != 0),
            None,
        )
        require(pivot is not None, "singular reconstruction matrix")
        work[column], work[pivot] = work[pivot], work[column]
        pivot_value = work[column][column]
        for entry in range(column, size + 1):
            work[column][entry] /= pivot_value
        for row in range(size):
            if row == column:
                continue
            factor = work[row][column]
            for entry in range(column, size + 1):
                work[row][entry] -= factor * work[column][entry]
    return [work[index][-1] for index in range(size)]


def resultant(poly_f, poly_g):
    poly_f = trim(poly_f)
    poly_g = trim(poly_g)
    degree_f = len(poly_f) - 1
    degree_g = len(poly_g) - 1
    require(degree_f >= 1 and degree_g >= 0, "resultant degree boundary")
    if degree_g == 0:
        return poly_g[0] ** degree_f
    width = degree_f + degree_g
    rows = []
    descending_f = list(reversed(poly_f))
    descending_g = list(reversed(poly_g))
    for shift in range(degree_g):
        row = [Fraction(0)] * width
        row[shift : shift + degree_f + 1] = descending_f
        rows.append(row)
    for shift in range(degree_f):
        row = [Fraction(0)] * width
        row[shift : shift + degree_g + 1] = descending_g
        rows.append(row)
    return determinant(rows)


def discriminant(monic_poly):
    monic_poly = monic(monic_poly)
    degree = len(monic_poly) - 1
    sign = -1 if degree * (degree - 1) // 2 % 2 else 1
    return sign * resultant(monic_poly, derivative(monic_poly))


def class_multiply(left, right, modulus):
    return remainder(mul(left, right), modulus)


def class_power(poly, exponent, modulus):
    out = [Fraction(1)]
    for _ in range(exponent):
        out = class_multiply(out, poly, modulus)
    return out


def contact_class(poly_f, poly_g, stratum, order):
    numerator = scale(omega(poly_f, poly_g, order), -1)
    denominator = mul(derivative(poly_f), derivative(poly_g))
    return remainder(mul(numerator, inverse_mod(denominator, stratum)), stratum)


def multiplication_matrix(poly, modulus):
    degree = len(modulus) - 1
    columns = []
    for shift in range(degree):
        product = class_multiply(poly, [Fraction(0)] * shift + [Fraction(1)], modulus)
        columns.append(product + [Fraction(0)] * (degree - len(product)))
    return [[columns[column][row] for column in range(degree)] for row in range(degree)]


def permutation_sign(permutation):
    inversions = sum(
        permutation[left] > permutation[right]
        for left in range(len(permutation))
        for right in range(left + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def determinant_polynomial(matrix):
    size = len(matrix)
    out = [Fraction(0)]
    for permutation in permutations(range(size)):
        term = [Fraction(permutation_sign(permutation))]
        for row in range(size):
            term = mul(term, matrix[row][permutation[row]])
        out = add(out, term)
    return trim(out)


def characteristic_polynomial(poly, modulus):
    matrix = multiplication_matrix(poly, modulus)
    size = len(matrix)
    polynomial_matrix = []
    for row in range(size):
        entries = []
        for column in range(size):
            entry = [Fraction(-matrix[row][column])]
            if row == column:
                entry = add(entry, [Fraction(0), Fraction(1)])
            entries.append(entry)
        polynomial_matrix.append(entries)
    return monic(determinant_polynomial(polynomial_matrix))


def interpolate(nodes, values):
    out = [Fraction(0)]
    for index, node in enumerate(nodes):
        basis = [Fraction(1)]
        denominator = Fraction(1)
        for other_index, other in enumerate(nodes):
            if index == other_index:
                continue
            basis = mul(basis, [Fraction(-other), Fraction(1)])
            denominator *= Fraction(node - other)
        out = add(out, scale(basis, Fraction(values[index], 1) / denominator))
    return trim(out)


def spectrum_from_resultants(poly_f, poly_g, stratum, order):
    degree = len(stratum) - 1
    denominator = mul(derivative(poly_f), derivative(poly_g))
    omega_order = omega(poly_f, poly_g, order)
    denominator_resultant = resultant(stratum, denominator)
    require(denominator_resultant != 0, "simple-root denominator")
    nodes = list(range(degree + 1))
    values = [
        resultant(stratum, add(scale(denominator, node), omega_order))
        / denominator_resultant
        for node in nodes
    ]
    return monic(interpolate(nodes, values))


def krylov_matrix(poly, modulus):
    degree = len(modulus) - 1
    columns = []
    for exponent in range(degree):
        value = class_power(poly, exponent, modulus)
        columns.append(value + [Fraction(0)] * (degree - len(value)))
    return [[columns[column][row] for column in range(degree)] for row in range(degree)]


def reconstruct_x(poly, modulus):
    matrix = krylov_matrix(poly, modulus)
    degree = len(modulus) - 1
    x_class = remainder([Fraction(0), Fraction(1)], modulus)
    target = x_class + [Fraction(0)] * (degree - len(x_class))
    coefficients = solve_linear(matrix, target)
    reconstructed = [Fraction(0)]
    for exponent, coefficient in enumerate(coefficients):
        reconstructed = add(
            reconstructed,
            scale(class_power(poly, exponent, modulus), coefficient),
        )
    require(remainder(reconstructed, modulus) == x_class, "root reconstruction")
    return trim(coefficients)


def spectrum_gate(poly_f, poly_g, stratum, order):
    contact = contact_class(poly_f, poly_g, stratum, order)
    spectrum_matrix = characteristic_polynomial(contact, stratum)
    spectrum_resultant = spectrum_from_resultants(poly_f, poly_g, stratum, order)
    require(spectrum_matrix == spectrum_resultant, "spectrum resultant formula")
    repeated = gcd_poly(spectrum_matrix, derivative(spectrum_matrix))
    minimal = monic(quotient_exact(spectrum_matrix, repeated))
    krylov = krylov_matrix(contact, stratum)
    kappa = determinant(krylov)
    require(
        discriminant(spectrum_matrix) == kappa * kappa * discriminant(stratum),
        "spectral discriminant/Krylov identity",
    )
    generator = kappa != 0
    require(generator == (len(minimal) == len(stratum)), "primitive-element iff")
    reconstruction = reconstruct_x(contact, stratum) if generator else None
    return contact, spectrum_matrix, minimal, kappa, reconstruction


# Quadratic power family: parity alternates primitive and collapsed spectra.
quadratic_p = [Fraction(-2), Fraction(0), Fraction(1)]
quadratic_spectra = []
quadratic_generator_pattern = []
quadratic_checks = 0
for order in range(2, 7):
    poly_f = quadratic_p
    poly_g = add(quadratic_p, power(quadratic_p, order))
    contact, spectrum, minimal, kappa, reconstruction = spectrum_gate(
        poly_f, poly_g, quadratic_p, order
    )
    expected_contact = remainder(
        power(derivative(quadratic_p), order - 1), quadratic_p
    )
    require(contact == expected_contact, "quadratic power contact")
    generator = kappa != 0
    require(generator == (order % 2 == 0), "quadratic parity gate")
    if generator:
        require(reconstruction is not None, "primitive reconstruction exists")
    else:
        require(len(minimal) == 2, "collapsed spectrum has one value")
    quadratic_spectra.append(spectrum)
    quadratic_generator_pattern.append("P" if generator else "C")
    quadratic_checks += 1

require(quadratic_spectra[0] == [Fraction(-8), 0, 1], "m2 spectrum")
require(quadratic_spectra[1] == [64, -16, 1], "m3 doubled spectrum")


# Delayed resplitting of the collapsed m=3 quadratic spectrum by the next raw jet.
quadratic_g3 = add(quadratic_p, power(quadratic_p, 3))
quadratic_c3 = contact_class(quadratic_p, quadratic_g3, quadratic_p, 3)
quadratic_d4 = contact_class(quadratic_p, quadratic_g3, quadratic_p, 4)
require(quadratic_c3 == [8] and quadratic_d4 == [0, 6], "quadratic delayed jets")
quadratic_combined = add(quadratic_c3, quadratic_d4)
combined_spectrum = characteristic_polynomial(quadratic_combined, quadratic_p)
combined_kappa = determinant(krylov_matrix(quadratic_combined, quadratic_p))
combined_reconstruction = reconstruct_x(quadratic_combined, quadratic_p)
require(combined_spectrum == [-8, -16, 1], "quadratic delayed spectrum")
require(discriminant(combined_spectrum) == 288, "quadratic delayed discriminant")
require(combined_kappa == 6, "quadratic delayed Krylov determinant")
require(combined_reconstruction == [Fraction(-4, 3), Fraction(1, 6)], "x=(h-8)/6")


# Cubic collision: c2 glues two roots, while (c2,d3) jointly separates them.
cubic_p = [Fraction(0), Fraction(-3), Fraction(0), Fraction(1)]
cubic_f = cubic_p
cubic_g = add(cubic_p, power(cubic_p, 2))
cubic_c2, cubic_spectrum, cubic_minimal, cubic_kappa, _ = spectrum_gate(
    cubic_f, cubic_g, cubic_p, 2
)
expected_cubic_spectrum = mul([3, 1], power([-6, 1], 2))
expected_cubic_minimal = mul([3, 1], [-6, 1])
require(cubic_c2 == [-3, 0, 3], "cubic c2=P'")
require(cubic_spectrum == expected_cubic_spectrum, "cubic doubled eigenvalue")
require(cubic_minimal == expected_cubic_minimal, "cubic minimal polynomial")
require(cubic_kappa == 0, "cubic c2 is not primitive")

cubic_d3 = contact_class(cubic_f, cubic_g, cubic_p, 3)
require(cubic_d3 == [0, 6], "cubic next raw jet")
cubic_combined = add(cubic_c2, cubic_d3)
cubic_combined_spectrum = characteristic_polynomial(cubic_combined, cubic_p)
cubic_combined_kappa = determinant(krylov_matrix(cubic_combined, cubic_p))
cubic_reconstruction = reconstruct_x(cubic_combined, cubic_p)
require(cubic_combined_kappa != 0, "cubic delayed primitive element")
require(
    discriminant(cubic_combined_spectrum)
    == cubic_combined_kappa * cubic_combined_kappa * discriminant(cubic_p),
    "cubic delayed discriminant identity",
)


# Exact exceptional-slope polynomial for c2+lambda*d3.
lambda_nodes = list(range(4))
lambda_kappas = []
for value in lambda_nodes:
    combined = add(cubic_c2, scale(cubic_d3, value))
    lambda_kappas.append(determinant(krylov_matrix(combined, cubic_p)))
exceptional_polynomial = interpolate(lambda_nodes, lambda_kappas)
expected_shape = mul([0, 1], [-3, 0, 4])
exceptional_scale = exceptional_polynomial[-1] / expected_shape[-1]
require(
    exceptional_polynomial == scale(expected_shape, exceptional_scale),
    "three pair-collision slopes",
)
require(exceptional_scale != 0, "exceptional polynomial is nonzero")


# Abstract small-field boundary: a jointly separating pair need not contain a
# primitive base-field linear combination when the field is too small.
f2_pairs = ((0, 0), (0, 1), (1, 0))
require(len(set(f2_pairs)) == 3, "F2 pair jointly separates")
f2_primitive_counts = []
for lam in (0, 1):
    values = {(left + lam * right) % 2 for left, right in f2_pairs}
    f2_primitive_counts.append(len(values))
require(f2_primitive_counts == [2, 2], "F2 has no primitive pencil member")


def format_poly(poly):
    return "[" + ",".join(str(value) for value in poly) + "]"


print("dependency_thm3232_sha256=%s" % DEPENDENCY_SHA256)
print("assert_nodes=%d,float_literals=%d" % (assert_nodes, float_literals))
print("quadratic_power_orders=2..6,checks=%d" % quadratic_checks)
print("quadratic_generator_pattern=%s" % "".join(quadratic_generator_pattern))
print("quadratic_m2_spectrum=%s" % format_poly(quadratic_spectra[0]))
print("quadratic_m3_spectrum=%s" % format_poly(quadratic_spectra[1]))
print("quadratic_delayed=spectrum%s,kappa%s" % (format_poly(combined_spectrum), combined_kappa))
print("cubic_collapsed_spectrum=%s" % format_poly(cubic_spectrum))
print("cubic_delayed_spectrum=%s" % format_poly(cubic_combined_spectrum))
print("cubic_exceptional_kappa=%s" % format_poly(exceptional_polynomial))
print("discriminant_identity=Disc(chi)=Disc(S)*kappa^2")
print("f2_joint_pair_boundary=primitive_value_counts_%s" % ",".join(str(v) for v in f2_primitive_counts))
print("scope=coefficient-label-reconstruction-not-root-order-or-physical-carrier")
print("all_exact_checks=PASS")
