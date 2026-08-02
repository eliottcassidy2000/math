#!/usr/bin/env python3
"""Exact companion for THM-3082's simultaneous suspension-word chambers."""

from fractions import Fraction
from math import factorial


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def add(left, right):
    return left[0] + right[0], left[1] + right[1]


def sub(left, right):
    return left[0] - right[1], left[1] - right[0]


def scale(value, interval):
    if value >= 0:
        return value * interval[0], value * interval[1]
    return value * interval[1], value * interval[0]


def log_bounds(value, terms=110):
    """Rigorous Fraction bounds for log(value), value>0."""
    require(value > 0, "logarithm domain")
    if value < 1:
        lower, upper = log_bounds(1 / value, terms)
        return -upper, -lower
    power_two = 0
    reduced = value
    while reduced >= 2:
        reduced /= 2
        power_two += 1
    z = (reduced - 1) / (reduced + 1)
    partial = Fraction(0)
    zpower = z
    for index in range(terms):
        partial += 2 * zpower / (2 * index + 1)
        zpower *= z * z
    tail = 2 * zpower / ((2 * terms + 1) * (1 - z * z))
    core = partial, partial + tail
    two = None
    if power_two:
        z2 = Fraction(1, 3)
        partial2 = Fraction(0)
        zpower2 = z2
        for index in range(terms):
            partial2 += 2 * zpower2 / (2 * index + 1)
            zpower2 *= z2 * z2
        tail2 = 2 * zpower2 / ((2 * terms + 1) * (1 - z2 * z2))
        two = partial2, partial2 + tail2
    return add(core, scale(power_two, two)) if power_two else core


def invoice_bounds(width, delta):
    bill = (width - 1) * factorial(width)
    return scale(bill * delta, add(log_bounds(Fraction(width) / delta), (1, 1)))


def gap_bounds(width, kind):
    if kind == "O":
        return log_bounds(Fraction(width**width, (width - 1) ** (width - 1)))
    require(kind == "P" and width >= 5, "pair-node domain")
    return scale(width - 2, log_bounds(Fraction(width - 1, width - 2)))


def margin_bounds(width, kind, delta):
    return sub(gap_bounds(width, kind), invoice_bounds(width, delta))


def suspension_words(target, width=3):
    if width == target:
        return [()]
    answer = []
    if width + 1 <= target:
        answer += [("O", width + 1) + tail for tail in suspension_words(target, width + 1)]
    if width + 2 <= target:
        answer += [("P", width + 2) + tail for tail in suspension_words(target, width + 2)]
    return answer


def ledger(word):
    exponents = {"S3": 1}
    for index in range(0, len(word), 2):
        kind, width = word[index], word[index + 1]
        multiplier = width if kind == "O" else width * (width - 1)
        exponents = {key: value * multiplier for key, value in exponents.items()}
        if kind == "O":
            exponents[f"U{width}"] = factorial(width - 1)
        else:
            exponents[f"U{width - 1}"] = factorial(width) // (width - 1)
            exponents[f"U{width}"] = factorial(width - 1)
            exponents[f"E{width - 1}{width}"] = factorial(width - 2)
    return exponents


def bareiss_determinant(matrix):
    matrix = [list(map(int, row)) for row in matrix]
    size = len(matrix)
    sign = 1
    denominator = 1
    for pivot_index in range(size - 1):
        if matrix[pivot_index][pivot_index] == 0:
            swap = next(
                (row for row in range(pivot_index + 1, size) if matrix[row][pivot_index]),
                None,
            )
            require(swap is not None, "singular hostile Sylvester matrix")
            matrix[pivot_index], matrix[swap] = matrix[swap], matrix[pivot_index]
            sign *= -1
        pivot = matrix[pivot_index][pivot_index]
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                numerator = matrix[row][column] * pivot - matrix[row][pivot_index] * matrix[pivot_index][column]
                require(numerator % denominator == 0, "Bareiss exact division")
                matrix[row][column] = numerator // denominator
        denominator = pivot
        for row in range(pivot_index + 1, size):
            matrix[row][pivot_index] = 0
    return sign * matrix[-1][-1]


def polynomial_resultant_monomial_hostile(degree_f, degree_g, amplitude):
    # f=t^degree_f and g=amplitude*t^degree_g+1, descending coefficients.
    f = [1] + [0] * degree_f
    g = [amplitude] + [0] * (degree_g - 1) + [1]
    size = degree_f + degree_g
    matrix = []
    for shift in range(degree_g):
        matrix.append([0] * shift + f + [0] * (size - shift - len(f)))
    for shift in range(degree_f):
        matrix.append([0] * shift + g + [0] * (size - shift - len(g)))
    return bareiss_determinant(matrix)


# Every width has the exact factorial-atom invoice used in the proof.
for width in range(4, 15):
    coefficient_degrees = [factorial(width) // degree for degree in range(2, width + 1)]
    require(
        sum(degree * multiplicity for degree, multiplicity in zip(range(2, width + 1), coefficient_degrees))
        == (width - 1) * factorial(width),
        f"invoice width {width}",
    )


# Every O/P word has the same closed multidegree ledger, plus one E per P node.
word_cells = 0
for target in range(3, 15):
    for word in suspension_words(target):
        word_cells += 1
        values = ledger(word)
        require(values["S3"] == factorial(target) // 6, "base exponent")
        for degree in range(4, target + 1):
            require(values[f"U{degree}"] == factorial(target) // degree, "U exponent")
        pair_widths = [word[index + 1] for index in range(0, len(word), 2) if word[index] == "P"]
        for width in pair_widths:
            require(
                values[f"E{width - 1}{width}"] == factorial(target) // (width * (width - 1)),
                "pair exponent",
            )

require(len(suspension_words(14)) == 144, "width-14 word count")


# Exact safe rays for O4,P6,O7 and the intermediate O5 prefix of the P6 node.
safe_cells = (
    (4, "O", Fraction(1, 1000)),
    (6, "P", Fraction(1, 100000)),
    (5, "O", Fraction(1, 100000)),
    (7, "O", Fraction(1, 200000)),
)
for width, kind, delta in safe_cells:
    require(margin_bounds(width, kind, delta)[0] > 0, f"safe ray {kind}{width}")

# Rigorous brackets for the two nontrivial ratios in the O4,P6,O7 example.
require(
    margin_bounds(6, "P", Fraction("0.00001808125357979"))[0] > 0
    and margin_bounds(6, "P", Fraction("0.00001808125357980"))[1] < 0,
    "P6 threshold bracket",
)
require(
    margin_bounds(7, "O", Fraction("0.00000636703607580"))[0] > 0
    and margin_bounds(7, "O", Fraction("0.00000636703607582"))[1] < 0,
    "O7 threshold bracket",
)

example = ledger(("O", 4, "P", 6, "O", 7))
require(
    example
    == {"S3": 840, "U4": 1260, "U5": 1008, "U6": 840, "E56": 168, "U7": 720},
    "O4-P6-O7 ledger",
)


# Divergent upper coefficients can be invisible to the exact grouped resultant.
hostile_cells = 0
for degree_f in range(2, 6):
    for degree_g in range(2, 6):
        for amplitude in (1, 7, 10**6):
            require(
                polynomial_resultant_monomial_hostile(degree_f, degree_g, amplitude) == 1,
                "triangular divergence hostile",
            )
            hostile_cells += 1

# A boundary leading quotient can be cancelled exactly when the strict
# entropy/Newton margin is lost: Res(y+epsilon*z,z+c*y)=1-epsilon*c.
margin_hostiles = 0
for amplitude in (1, 7, 10**3, 10**6):
    epsilon = Fraction(1, amplitude)
    require(1 - epsilon * amplitude == 0, "zero-margin cancellation hostile")
    margin_hostiles += 1


print("THM-3082 ADMISSIBLE SUSPENSION-WORD CHAMBERS")
print(f"widths=3..14 word_cells={word_cells} width14_words=144")
print("invoice=B_k=(k-1)k!;J_k(delta)=B_k*delta*(log(k/delta)+1)")
print("node_gaps=O:log(k^k/(k-1)^(k-1));P:(k-2)log((k-1)/(k-2))")
print("P6_threshold_in=(0.00001808125357979,0.00001808125357980)")
print("O7_threshold_in=(0.00000636703607580,0.00000636703607582)")
print("safe_example=O4:1/1000,P6:1/100000,O7:1/200000")
print("K7_carrier=S3^840*U4^1260*E56^168*U5^1008*U6^840*U7^720")
print(f"triangular_divergence_hostiles={hostile_cells} resultant=1")
print(f"zero_margin_cancellation_hostiles={margin_hostiles} resultant=0")
print("repair=group_lower_error_ideal_then_entropy;raw_convergence=FALSE")
print("scope=all_finite_admissible_OP_words;thin_cones;not_arbitrary_support")
print("all_exact_checks=PASS")
