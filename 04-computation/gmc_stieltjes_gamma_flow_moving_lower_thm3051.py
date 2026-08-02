#!/usr/bin/env python3
"""Exact referee for THM-3051 using dependency-free rational arithmetic."""

from fractions import Fraction
from hashlib import sha256
from math import comb, factorial, gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def rising(x, length):
    out = Fraction(1)
    for step in range(length):
        out *= x + step
    return out


def moments(atoms, count):
    return [sum(weight * point**m for point, weight in atoms) for m in range(count)]


def determinant_fraction(matrix):
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


def lcm(left, right):
    return abs(left * right) // gcd(left, right) if left and right else 0


def determinant_bareiss_fraction(matrix):
    """Clear denominators rowwise, then use fraction-free Bareiss elimination."""
    integer_matrix = []
    denominator = 1
    for row in matrix:
        row_lcm = 1
        fractions = [Fraction(entry) for entry in row]
        for entry in fractions:
            row_lcm = lcm(row_lcm, entry.denominator)
        integer_matrix.append([int(entry * row_lcm) for entry in fractions])
        denominator *= row_lcm

    a = integer_matrix
    size = len(a)
    sign = 1
    previous = 1
    for column in range(size - 1):
        pivot_row = next((row for row in range(column, size) if a[row][column]), None)
        if pivot_row is None:
            return Fraction(0)
        if pivot_row != column:
            a[column], a[pivot_row] = a[pivot_row], a[column]
            sign = -sign
        pivot = a[column][column]
        for row in range(column + 1, size):
            for j in range(column + 1, size):
                numerator = a[row][j] * pivot - a[row][column] * a[column][j]
                require(numerator % previous == 0, "Bareiss division lost exactness")
                a[row][j] = numerator // previous
            a[row][column] = 0
        previous = pivot
    return Fraction(sign * a[-1][-1], denominator)


def width_flag(slot_count, t, m):
    harmonic = sum(Fraction(1, j) for j in range(1, slot_count + 1))
    a_count = factorial(slot_count) * (harmonic - 1)
    b_count = factorial(slot_count) * (slot_count + 1 - 2 * harmonic)
    interior = a_count + b_count
    require(a_count.denominator == b_count.denominator == 1, "nonintegral character")
    out = Fraction(1)
    for s in range(1, m):
        out *= (1 + s * t) ** int(interior)
    if m:
        out *= (1 + m * t) ** int(b_count)
    return out


print("THM-3051 STIELTJES MULTIPLIERS, GAMMA FLOW, AND MOVING-LOWER BOUNDARY")

# Finite-atomic exact controls for the universal Hadamard multiplier mechanism.
multiplier_cells = 0
base_measures = (
    ((Fraction(1, 2), Fraction(2, 3)), (Fraction(3), Fraction(1, 3))),
    ((Fraction(2), Fraction(1, 4)), (Fraction(5), Fraction(3, 4))),
    ((Fraction(1, 3), Fraction(1, 5)), (Fraction(4, 3), Fraction(2, 5)), (Fraction(7), Fraction(2, 5))),
)
factor_measures = (
    ((Fraction(0), Fraction(1, 5)), (Fraction(2), Fraction(4, 5))),
    ((Fraction(1), Fraction(1, 2)), (Fraction(3), Fraction(1, 2))),
    ((Fraction(5, 2), Fraction(1)),),
)
for mu in base_measures:
    f = moments(mu, 13)
    for lam in factor_measures:
        ell = moments(lam, 13)
        product_atoms = tuple((x * y, wx * wy) for x, wx in mu for y, wy in lam)
        product_moments = moments(product_atoms, 13)
        require([f[m] * ell[m] for m in range(13)] == product_moments, "product moment identity failed")
        for shift in range(4):
            for size in (1, 2, 3):
                block = [[product_moments[shift + i + j] for j in range(size)] for i in range(size)]
                # Finite atomic controls can have fewer than `size` distinct
                # product nodes, so semidefiniteness (not strictness) is the
                # correct executable assertion here.
                require(determinant_fraction(block) >= 0, "positive atomic product lost Hankel positivity")
                multiplier_cells += 1
print(f"multiplier_exact_cells={multiplier_cells}")

# Adjacent Gamma-shape flow and a non-Stieltjes multiplier preserved by F.
flow_cells = 0
beta_escape_cells = 0
for slot_count in (2, 3, 4):
    for t in (Fraction(1, 2), Fraction(2, 3), Fraction(1)):
        a = 1 / t
        harmonic = sum(Fraction(1, j) for j in range(1, slot_count + 1))
        a_count = int(factorial(slot_count) * (harmonic - 1))
        b_count = int(factorial(slot_count) * (slot_count + 1 - 2 * harmonic))
        for m in range(9):
            base = width_flag(slot_count, t, m)
            transfer = (a + m) / a
            shifted = t ** ((a_count + b_count) * m) * rising(a, m) ** (a_count - 1) * rising(a + 1, m) ** (b_count + 1)
            require(base * transfer == shifted, "one-step Gamma flow failed")
            flow_cells += 1
            if slot_count == 2:
                inverse_transfer = a / (a + m)
                beta_product = t**m * rising(a, m) * a / (a + m)
                require(base * inverse_transfer == beta_product,
                        "negative-inventory Beta escape failed")
                beta_escape_cells += 1
non_stieltjes_t = Fraction(2, 3)
non_stieltjes_h2 = (1 + 2 * non_stieltjes_t) - (1 + non_stieltjes_t) ** 2
require(non_stieltjes_h2 == -(non_stieltjes_t**2), "non-Stieltjes hostile changed")
print(f"gamma_flow_cells={flow_cells} beta_escape_cells={beta_escape_cells} nonstieltjes_multiplier_h2={non_stieltjes_h2}")

# Universal order-three curvature identity and a strict-log-convex hostile.
curvature_cells = 0
for a in (Fraction(9, 8), Fraction(5, 4), Fraction(3, 2)):
    for b in (Fraction(10, 9), Fraction(4, 3)):
        for c in (Fraction(17, 16), Fraction(7, 5)):
            q = Fraction(5, 3)
            g = Fraction(7, 4)
            values = [g]
            ratios = [q, q * a, q * a * b, q * a * b * c]
            for ratio in ratios:
                values.append(values[-1] * ratio)
            direct = determinant_fraction([[values[i + j] for j in range(3)] for i in range(3)])
            psi = a * b * b * c - a * b * b - b * b * c + 2 * b - 1
            predicted = g**3 * q**6 * a**3 * psi
            require(direct == predicted, "order-three curvature identity failed")
            curvature_cells += 1

strict_log_convex = [factorial(m) * (m + 1) ** 2 for m in range(5)]
for m in range(1, 4):
    require(strict_log_convex[m - 1] * strict_log_convex[m + 1] > strict_log_convex[m] ** 2,
            "strict log-convex hostile lost strictness")
strict_log_convex_h3 = determinant_fraction([[strict_log_convex[i + j] for j in range(3)] for i in range(3)])
require(strict_log_convex_h3 == -24, "strict-log-convex H3 hostile changed")
print(f"curvature_identity_cells={curvature_cells} strict_logconvex_h3={strict_log_convex_h3}")


def binary_form(order, lower, depth):
    """Descending X,Y coefficients of the normalized low factorial form."""
    coefficients = []
    for j in range(order + 1):
        value = Fraction(comb(order, j))
        value *= rising(order * depth + 1, j * lower)
        value /= rising(depth + 1, lower) ** j
        coefficients.append(value)
    return coefficients


def binary_resultant(f, g):
    m = len(f) - 1
    n = len(g) - 1
    matrix = []
    for shift in range(n):
        matrix.append([Fraction(0)] * shift + list(f) + [Fraction(0)] * (n - 1 - shift))
    for shift in range(m):
        matrix.append([Fraction(0)] * shift + list(g) + [Fraction(0)] * (m - 1 - shift))
    return determinant_fraction(matrix)


def eval_desc(coefficients, x):
    out = 0
    for coefficient in coefficients:
        out = out * x + coefficient
    return out


A_DESC = [1, 22, 69, 44, 8]
B_DESC = [15625, 427500, 2850950, 8097488, 11893269, 9683740, 4378980, 1022112, 95616]
C_DESC = [47045881, 1797908516, 23919687803, 167453540986, 722042412723,
          2070797182104, 4110657369269, 5751069261826, 5683612606852,
          3918712637960, 1829514912288, 546467611968, 93389808384, 6890365440]
P_SHIFT_ASC = [
    3268881760112640000, 37643475128286720000, 146291966980162617344,
    311943061093012269312, 433931551733794554880, 429003646190355919616,
    316601230236031112592, 179745211720414744752, 79988104378771836568,
    28208693320266071200, 7923650163754635817, 1772672109839816584,
    314264085331109200, 43683300848187864, 4678174175000898,
    375655828575392, 21687531118480, 840873878616, 19348674701, 197094744,
]


def polynomial_add(left, right, sign=1):
    out = [0] * max(len(left), len(right))
    for i, coefficient in enumerate(left):
        out[i] += coefficient
    for i, coefficient in enumerate(right):
        out[i] += sign * coefficient
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def polynomial_multiply(left, right):
    out = [0] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            out[i + j] += x * y
    return out


def shift_polynomial(coefficients, shift):
    """Ascending coefficients of P(m+shift) from ascending coefficients of P(n)."""
    out = [0] * len(coefficients)
    for degree, coefficient in enumerate(coefficients):
        for power in range(degree + 1):
            out[power] += coefficient * comb(degree, power) * shift ** (degree - power)
    return out


# Three-slot literal moving-lower family: exact resultants and all-depth sign.
k3_cells = 0
for depth in range(1, 33):
    direct = [binary_resultant(binary_form(2, lower, depth), binary_form(3, lower, depth))
              for lower in (1, 2, 3)]
    expected = [
        Fraction(eval_desc(A_DESC, depth), (depth + 1) ** 4),
        Fraction(eval_desc(B_DESC, depth), (depth + 1) ** 4 * (depth + 2) ** 4),
        Fraction(eval_desc(C_DESC, depth), (depth + 1) ** 4 * (depth + 2) ** 5 * (depth + 3) ** 4),
    ]
    require(direct == expected, "three-slot moving-lower resultant formula failed")
    k3_cells += 3

a_asc = list(reversed(A_DESC))
b_asc = list(reversed(B_DESC))
c_asc = list(reversed(C_DESC))
denominator_polynomial = polynomial_multiply(
    polynomial_multiply([4, 4, 1], [3, 1]), polynomial_multiply(b_asc, b_asc))
numerator_polynomial = polynomial_multiply(polynomial_multiply([16, 8, 1], a_asc), c_asc)
p_polynomial = polynomial_add(denominator_polynomial, numerator_polynomial, sign=-1)
require(shift_polynomial(p_polynomial, 2) == P_SHIFT_ASC, "all-depth shifted polynomial certificate failed")
require(all(coefficient > 0 for coefficient in P_SHIFT_ASC), "shifted positivity certificate lost a coefficient")
require(sum(coefficient for coefficient in p_polynomial) == -36470140909977600, "depth-one boundary changed")

k_depth_one = Fraction((1 + 4) ** 2 * eval_desc(A_DESC, 1) * eval_desc(C_DESC, 1),
                       (1 + 2) ** 2 * (1 + 3) * eval_desc(B_DESC, 1) ** 2)
require(k_depth_one == Fraction(3577625, 2123604) and k_depth_one > 1,
        "depth-one positive boundary changed")
for lower in (1, 2, 3):
    require((3**lower - 2**lower) ** 6 > 0, "three-slot t=0 resultant control failed")

depth = 2
r1, r2, r3 = [binary_resultant(binary_form(2, lower, depth), binary_form(3, lower, depth))
              for lower in (1, 2, 3)]
h1 = width_flag(3, Fraction(1, depth), 2) * r1
h2 = width_flag(3, Fraction(1, depth), 3) * r2
h3 = width_flag(3, Fraction(1, depth), 4) * r3
require((r1, r2, r3) == (Fraction(188, 27), Fraction(1025675, 27), Fraction(44055606)),
        "three-slot depth-two resultants changed")
require(h1 * h3 - h2 * h2 == Fraction(-9090480053413828125, 512),
        "three-slot depth-two Hankel hostile changed")
p_checksum = sha256((",".join(str(value) for value in P_SHIFT_ASC) + "\n").encode()).hexdigest()
require(p_checksum == "9104ea7c56f38d91905d0baeaf39bb0a7afc17a705232df076171523e8408c96",
        "shifted polynomial checksum changed")
print(f"k3_resultant_cells={k3_cells} depth1_curvature=positive depth_ge_2_curvature=negative")
print(f"k3_shift_certificate_sha256={p_checksum} depth2_hankel={h1 * h3 - h2 * h2}")


def compositions(total, variables=3):
    if variables == 1:
        return [(total,)]
    out = []
    for first in range(total + 1):
        for suffix in compositions(total - first, variables - 1):
            out.append((first,) + suffix)
    return out


MON = {degree: compositions(degree) for degree in range(8)}


def multinomial(alpha):
    out = factorial(sum(alpha))
    for entry in alpha:
        out //= factorial(entry)
    return out


def ternary_form(order, lower, depth=2):
    out = {}
    for alpha in MON[order]:
        offset_sum = alpha[1] + lower * alpha[2]
        coefficient = Fraction(multinomial(alpha) * rising(order * depth + 1, offset_sum),
                               rising(depth + 1, 1) ** alpha[1] * rising(depth + 1, lower) ** alpha[2])
        out[alpha] = coefficient
    return out


SELECTED_ROWS = tuple(range(20)) + tuple(range(21, 30)) + (35,) + tuple(range(36, 42))


def ternary_resultant_selected(lower):
    columns = {monomial: index for index, monomial in enumerate(MON[7])}
    all_rows = []
    forms = {order: ternary_form(order, lower) for order in (2, 3, 4)}
    for order in (2, 3, 4):
        for beta in MON[7 - order]:
            row = [Fraction(0)] * len(MON[7])
            for alpha, coefficient in forms[order].items():
                target = tuple(alpha[j] + beta[j] for j in range(3))
                row[columns[target]] += coefficient
            all_rows.append(row)
    matrix = [all_rows[index] for index in SELECTED_ROWS]
    delta = determinant_bareiss_fraction(matrix)
    q = forms[2]
    cubic = forms[3]
    q200 = q[(2, 0, 0)]
    q110 = q[(1, 1, 0)]
    q020 = q[(0, 2, 0)]
    c300 = cubic[(3, 0, 0)]
    c210 = cubic[(2, 1, 0)]
    c120 = cubic[(1, 2, 0)]
    flag = c120 * q200**2 - c210 * q110 * q200 - c300 * q020 * q200 + c300 * q110**2
    require(q200 and c300 and flag, "selected Macaulay chart hit its extraneous wall")
    return delta / (q200**6 * c300 * flag), delta, flag


EXPECTED_K4 = (
    Fraction(671696427641384054000000000, 1162261467),
    Fraction(4720147255045226732121521309597418358046720, 19683),
    Fraction(23877237441067576157124399040421454340333007132221361662211000369152,
             2017815046875),
)
k4_values = []
for lower, expected in zip((2, 3, 4), EXPECTED_K4):
    value, delta, flag = ternary_resultant_selected(lower)
    require(delta and flag and value == expected, "four-slot selected resultant changed")
    k4_values.append(value)

response_det = determinant_fraction([[Fraction(order**offset) for offset in (0, 1, 2)]
                                     for order in (2, 3, 4)])
require(response_det == 2 and response_det**24 == 2**24, "t=0 resultant control changed")
universal_curvature = Fraction(6, 5) ** 26 * Fraction(7, 6) ** 20
total_curvature = universal_curvature * k4_values[0] * k4_values[2] / (k4_values[1] ** 2)
expected_curvature = Fraction(
    44087099546187941338084870318252428808920867535835042464970681136337702301141312674129,
    148395225288083619662436115015878402556351144600145942926978139965896606445312500000000,
)
require(total_curvature == expected_curvature and total_curvature < 1,
        "four-slot moving-lower curvature hostile changed")
deficit = expected_curvature.denominator - expected_curvature.numerator
require(deficit == 104308125741895678324351244697625973747430277064310900462007458829558904144171187325871,
        "four-slot curvature deficit changed")
print("k4_macaulay_resultants=3 selected_matrix=36x36 t0_resultant=2^24")
print(f"k4_n2_total_curvature={total_curvature} curvature_lt_1={total_curvature < 1}")
print("all_exact_checks=PASS")
