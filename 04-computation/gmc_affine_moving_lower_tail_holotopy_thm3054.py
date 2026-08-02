#!/usr/bin/env python3
"""Exact referee for THM-3054 using dependency-free rational arithmetic."""

from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations, product
from math import comb, factorial


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def rising(x, length):
    out = Fraction(1)
    for step in range(length):
        out *= x + step
    return out


def determinant(matrix):
    a = [[Fraction(entry) for entry in row] for row in matrix]
    size = len(a)
    det = Fraction(1)
    for column in range(size):
        pivot = next((row for row in range(column, size) if a[row][column]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            a[column], a[pivot] = a[pivot], a[column]
            det = -det
        value = a[column][column]
        det *= value
        for j in range(column, size):
            a[column][j] /= value
        for row in range(column + 1, size):
            scale = a[row][column]
            if scale:
                for j in range(column, size):
                    a[row][j] -= scale * a[column][j]
    return det


def binary_resultant(f, g):
    """Sylvester resultant for descending coefficient lists."""
    degree_f = len(f) - 1
    degree_g = len(g) - 1
    matrix = []
    for shift in range(degree_g):
        matrix.append([0] * shift + list(f) + [0] * (degree_g - 1 - shift))
    for shift in range(degree_f):
        matrix.append([0] * shift + list(g) + [0] * (degree_f - 1 - shift))
    return determinant(matrix)


def character(slot_count):
    harmonic = sum(Fraction(1, j) for j in range(1, slot_count + 1))
    a_count = factorial(slot_count) * (harmonic - 1)
    b_count = factorial(slot_count) * (slot_count + 1 - 2 * harmonic)
    interior = a_count + b_count
    require(a_count.denominator == b_count.denominator == interior.denominator == 1,
            "character lost integrality")
    return int(a_count), int(b_count), int(interior)


def width_flag(slot_count, depth, width):
    a_count, b_count, interior = character(slot_count)
    t = Fraction(1, depth)
    out = Fraction(1)
    for s in range(1, width):
        out *= (1 + s * t) ** interior
    if width:
        out *= (1 + width * t) ** b_count
    return out


def top_coefficient(slot_count, depth, moving):
    return rising(slot_count * depth + 1, slot_count * moving) / rising(depth + 1, moving) ** slot_count


def dominant_carrier(slot_count, depth, gap, moving):
    mu = factorial(slot_count - 1)
    return width_flag(slot_count, depth, moving + gap) * top_coefficient(slot_count, depth, moving) ** mu


def k3_parameters(depth, lower):
    denominator = rising(depth + 1, lower)
    a = rising(2 * depth + 1, lower) / denominator
    b = rising(2 * depth + 1, 2 * lower) / denominator**2
    p = rising(3 * depth + 1, lower) / denominator
    q = rising(3 * depth + 1, 2 * lower) / denominator**2
    r = rising(3 * depth + 1, 3 * lower) / denominator**3
    return a, b, p, q, r


def k3_formula(a, b, p, q, r):
    return (r**2 - 6 * a * q * r - 6 * b * p * r + 6 * a * b * r
            - 8 * a**3 * r + 12 * a**2 * p * r + 12 * a**2 * b * q
            - 6 * a * b**2 * p - 18 * a * b * p * q + b**3
            + 9 * b**2 * p**2 - 6 * b**2 * q + 9 * b * q**2)


def k3_resultant(depth, lower):
    a, b, p, q, r = k3_parameters(depth, lower)
    return k3_formula(a, b, p, q, r)


def k3_corner(depth, gap, lower):
    return width_flag(3, depth, lower + gap) * k3_resultant(depth, lower)


print("THM-3054 AFFINE MOVING-LOWER TROPICAL BETA-GAMMA TAIL HOLOTOPY")

# The generic 13-term quadratic/cubic resultant identity.
generic_cells = 0
test_values = (Fraction(1, 3), Fraction(2, 5), Fraction(3, 2))
for a, b, p, q, r in product(test_values, repeat=5):
    direct = binary_resultant((b, 2 * a, 1), (r, 3 * q, 3 * p, 1))
    require(direct == k3_formula(a, b, p, q, r), "generic 13-term resultant failed")
    generic_cells += 1
print(f"generic_resultant_cells={generic_cells} terms=13")

# Exact all-depth/all-lower positivity, including the division remainder and
# the factorized termwise inequality which makes C0 positive.
positivity_cells = 0
termwise_cells = 0
for depth in range(1, 17):
    for lower in range(1, 21):
        a, b, p, q, r = k3_parameters(depth, lower)
        direct = binary_resultant((b, 2 * a, 1), (r, 3 * q, 3 * p, 1))
        formula = k3_formula(a, b, p, q, r)
        c0 = b**2 + 2 * a * r - 3 * b * q
        c1 = 4 * a**2 * r - 6 * a * b * q + 3 * b**2 * p - b * r
        remainder_norm = ((c1 - a * c0) ** 2 + (b - a**2) * c0**2) / b**2
        require(direct == formula == remainder_norm, "k3 remainder identity failed")
        require(b > a**2 and b * q**3 < a**2 * r**2 and c0 > 0 and direct > 0,
                "k3 positivity certificate failed")
        positivity_cells += 1
        for i in range(1, lower + 1):
            left = (2 * depth + i) * (3 * depth + 2 * lower + i) ** 2
            left -= (2 * depth + lower + i) * (3 * depth + i) * (3 * depth + lower + i)
            expected = lower * (3 * lower * i + 5 * lower * depth + 2 * i**2
                                + 9 * i * depth + 9 * depth**2)
            require(left == expected and expected > 0, "termwise positivity factor failed")
            termwise_cells += 1
print(f"k3_positive_cells={positivity_cells} termwise_factor_cells={termwise_cells}")

# Translation of arbitrary distinct lower pairs to (0,delta) at shifted depth.
translation_cells = 0
for depth in range(1, 7):
    for low in range(5):
        for high in range(low + 1, 8):
            delta = high - low
            shifted_depth = depth + low
            for order in (2, 3):
                scalar = rising(order * depth + 1, order * low) / rising(depth + 1, low) ** order
                for j in range(order + 1):
                    raw = (Fraction(comb(order, j))
                           * rising(order * depth + 1, (order - j) * low + j * high)
                           / (rising(depth + 1, low) ** (order - j)
                              * rising(depth + 1, high) ** j))
                    reduced = (scalar * comb(order, j)
                               * rising(order * shifted_depth + 1, j * delta)
                               / rising(shifted_depth + 1, delta) ** j)
                    require(raw == reduced, "lower-pair translation failed")
                    translation_cells += 1
print(f"lower_pair_translation_cells={translation_cells}")

# Exact Beta-Gamma factorization of the top tropical carrier, first on the
# unit clock and then on affine subsequences C=a_clock*c+b_clock.
carrier_cells = 0
affine_cells = 0
for slot_count in range(2, 6):
    a_count, b_count, interior = character(slot_count)
    mu = factorial(slot_count - 1)
    for depth in (1, 2, 3):
        for gap in (1, 2):
            base = dominant_carrier(slot_count, depth, gap, 0)
            scale = Fraction(slot_count ** factorial(slot_count), depth**interior)
            for moving in range(5):
                factored = base * scale**moving
                factored *= rising(depth + gap, moving) ** a_count
                factored *= rising(depth + gap + 1, moving) ** b_count
                for q_index in range(1, slot_count):
                    factored *= (rising(Fraction(depth * slot_count + q_index, slot_count), moving)
                                 / rising(depth + 1, moving)) ** mu
                require(factored == dominant_carrier(slot_count, depth, gap, moving),
                        "unit-clock carrier factorization failed")
                carrier_cells += 1

            for clock_slope in (1, 2, 3):
                for clock_offset in (0, 1, 2):
                    anchor = dominant_carrier(slot_count, depth, gap, clock_offset)
                    affine_scale = (Fraction(slot_count ** (clock_slope * factorial(slot_count)),
                                             depth ** (clock_slope * interior))
                                    * clock_slope ** (clock_slope * interior))
                    for index in range(4):
                        factored = anchor * affine_scale**index
                        for residue in range(clock_slope):
                            alpha = Fraction(depth + gap + clock_offset + residue, clock_slope)
                            factored *= rising(alpha, index) ** a_count
                            alpha = Fraction(depth + gap + 1 + clock_offset + residue, clock_slope)
                            factored *= rising(alpha, index) ** b_count
                            for q_index in range(1, slot_count):
                                beta = Fraction(depth + clock_offset + residue, clock_slope)
                                beta += Fraction(q_index, slot_count * clock_slope)
                                gamma = Fraction(depth + 1 + clock_offset + residue, clock_slope)
                                factored *= (rising(beta, index) / rising(gamma, index)) ** mu
                        moving = clock_slope * index + clock_offset
                        require(factored == dominant_carrier(slot_count, depth, gap, moving),
                                "affine-clock carrier factorization failed")
                        affine_cells += 1
print(f"carrier_factor_cells={carrier_cells} affine_subsequence_cells={affine_cells}")

# Sharp universal tropical redistribution gap.  Dynamic programming maximizes
# the exponential base after L missing pure top coefficients.
tropical_cells = 0
for slot_count in range(2, 11):
    rho = Fraction((slot_count - 1) ** (slot_count - 1), slot_count**slot_count)
    max_weight = slot_count * (slot_count - 1)
    best = [0] * (max_weight + 1)
    best[0] = 1
    for weight in range(1, max_weight + 1):
        best[weight] = max((best[weight - part] * part**part
                            for part in range(1, min(slot_count - 1, weight) + 1)), default=0)
    for lost in range(1, slot_count):
        ratio = Fraction(best[slot_count * lost], slot_count ** (slot_count * lost))
        require(ratio <= rho, "tropical gap bound failed")
        if lost == 1:
            require(ratio == rho, "sharp tropical gap equality failed")
        tropical_cells += 1
print(f"tropical_gap_cells={tropical_cells} sharp_rho_k=(k-1)^(k-1)/k^k")


# Polynomial alternant leading coefficient for product-Gamma carriers.
def poly_add(left, right, sign=1):
    out = [Fraction(0)] * max(len(left), len(right))
    for i, value in enumerate(left):
        out[i] += value
    for i, value in enumerate(right):
        out[i] += sign * value
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def poly_mul(left, right):
    out = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            out[i + j] += x * y
    return out


def poly_shift(poly, shift):
    out = [Fraction(0)] * len(poly)
    for degree, coefficient in enumerate(poly):
        for power in range(degree + 1):
            out[power] += coefficient * comb(degree, power) * shift ** (degree - power)
    return out


def permutation_sign(permutation):
    inversions = sum(permutation[i] > permutation[j]
                     for i in range(len(permutation)) for j in range(i + 1, len(permutation)))
    return -1 if inversions % 2 else 1


def polynomial_determinant(matrix):
    size = len(matrix)
    out = [Fraction(0)]
    for perm in permutations(range(size)):
        term = [Fraction(permutation_sign(perm))]
        for row in range(size):
            term = poly_mul(term, matrix[row][perm[row]])
        out = poly_add(out, term)
    return out


def vandermonde(values):
    out = 1
    for i in range(len(values)):
        for j in range(i + 1, len(values)):
            out *= values[j] - values[i]
    return out


alternant_cells = 0
for factor_count in (2, 3, 7):
    shapes = [Fraction(2 * j + 1, 3) for j in range(factor_count)]
    for rows, columns in (((0, 1), (0, 2)), ((0, 2, 3), (0, 1, 3))):
        size = len(rows)
        matrix = []
        for row in rows:
            matrix_row = []
            for column in columns:
                polynomial = [Fraction(1)]
                for shape in shapes:
                    for step in range(column):
                        polynomial = poly_mul(polynomial, [shape + step, 1])
                matrix_row.append(poly_shift(polynomial, row))
            matrix.append(matrix_row)
        det_poly = polynomial_determinant(matrix)
        degree = factor_count * sum(columns) - size * (size - 1) // 2
        factorial_denominator = 1
        for derivative_order in range(size):
            factorial_denominator *= factorial(derivative_order)
        expected_leader = (factor_count ** (size * (size - 1) // 2)
                           * vandermonde(rows) * vandermonde(columns)
                           / factorial_denominator)
        require(len(det_poly) - 1 == degree and det_poly[-1] == expected_leader,
                "alternant leading coefficient failed")
        alternant_cells += 1
print(f"alternant_leader_cells={alternant_cells}")

# Exact finite-window recovery controls for the actual k=3 moving corner and
# its straight carrier interpolation.  These illustrate, but do not prove,
# the all-tail theorem.
recovery_cells = 0
recovery_thresholds = []
for depth in (1, 2, 3):
    for gap in (1, 2):
        maximum_base = 18
        maximum_offset = 5
        actual = {c: k3_corner(depth, gap, c) for c in range(1, maximum_base + maximum_offset + 1)}
        carrier = {c: dominant_carrier(3, depth, gap, c) for c in actual}
        for size in (2, 3):
            passing = []
            for base_index in range(1, maximum_base + 1):
                okay = True
                for theta in (Fraction(0), Fraction(1, 2), Fraction(1)):
                    sequence = {c: (1 - theta) * carrier[c] + theta * actual[c] for c in actual}
                    block = [[sequence[base_index + i + j] for j in range(size)] for i in range(size)]
                    if determinant(block) <= 0:
                        okay = False
                    recovery_cells += 1
                passing.append(okay)
            threshold = next((base for base in range(1, maximum_base + 1)
                              if all(passing[base - 1:])), None)
            require(threshold is not None, "finite recovery window did not stabilize")
            recovery_thresholds.append((depth, gap, size, threshold))
print(f"finite_recovery_cells={recovery_cells} thresholds={recovery_thresholds}")

# Correction-only and scope hostiles.
depth = 1
correction = {}
for lower in range(1, 30):
    _, _, _, _, r = k3_parameters(depth, lower)
    correction[lower] = k3_resultant(depth, lower) / r**2
require((correction[1], correction[2], correction[3])
        == (Fraction(1, 25), Fraction(53, 140), Fraction(28621, 38115)),
        "correction values changed")
correction_h2 = correction[1] * correction[3] - correction[2] ** 2
require(correction_h2 == Fraction(-12089453, 106722000), "correction H2 hostile changed")

inverse_difference = sum((-1) ** (19 - j) * comb(19, j) / correction[5 + j]
                         for j in range(20))
require(inverse_difference > 0, "inverse correction order-19 hostile changed")
inverse_hash = sha256((str(inverse_difference) + "\n").encode()).hexdigest()

nonaffine = dominant_carrier(3, 1, 1, 1) * dominant_carrier(3, 1, 1, 4)
nonaffine -= dominant_carrier(3, 1, 1, 3) ** 2
varying_gap = dominant_carrier(3, 1, 1, 1) * dominant_carrier(3, 1, 1, 3)
varying_gap -= dominant_carrier(3, 1, 2, 2) ** 2
require(nonaffine == -14623245661460508379818491904000000000,
        "nonaffine-clock hostile changed")
require(varying_gap == -80810821433972705093222400000000,
        "varying-gap hostile changed")
print(f"correction_h2={correction_h2} inverse_delta19_sha256={inverse_hash}")
print(f"scope_hostiles=nonaffine:{nonaffine},varying_gap:{varying_gap}")
print("all_exact_checks=PASS")
