#!/usr/bin/env python3
"""Exact controls for THM-3186's exterior path convolution."""

from itertools import combinations

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def compound(matrix):
    """Second compound in the ordered pairs (01,02,12)."""

    pairs = ((0, 1), (0, 2), (1, 2))
    return sp.Matrix([
        [sp.det(matrix.extract(rows, columns)) for columns in pairs]
        for rows in pairs
    ])


def path_continuant(alpha, beta, scalar_d, length):
    """Matching expansion on vertices 1,...,length."""

    if length == 0:
        return sp.Integer(1), 1
    total = sp.Integer(0)
    count = 0
    # A dimer with right endpoint r covers vertices r-1,r.
    endpoints = tuple(range(2, length + 1))
    for size in range(len(endpoints) + 1):
        for chosen in combinations(endpoints, size):
            if any(right + 1 in chosen for right in chosen):
                continue
            covered = {vertex for right in chosen
                       for vertex in (right - 1, right)}
            weight = sp.prod(scalar_d * beta[right] for right in chosen)
            weight *= sp.prod(alpha[vertex]
                              for vertex in range(1, length + 1)
                              if vertex not in covered)
            total += weight
            count += 1
    return sp.expand(total), count


# Abstract direct matrix products versus the recurrence and exit convolution.
MAX_STEPS = 7
scalar_d = sp.symbols("D")
u = sp.symbols("u0:" + str(MAX_STEPS))
c = sp.symbols("c0:" + str(MAX_STEPS))
alpha = sp.symbols("alpha0:" + str(MAX_STEPS))
beta = sp.symbols("beta0:" + str(MAX_STEPS))
matrices = tuple(
    sp.Matrix([[u[index], -c[index], 0],
               [0, alpha[index], beta[index]],
               [0, scalar_d, 0]])
    for index in range(MAX_STEPS)
)
source = sp.Matrix([0, 0, 1])

t = [sp.Integer(0), beta[0]]
for index in range(1, MAX_STEPS):
    t.append(sp.expand(alpha[index] * t[index]
                       + beta[index] * scalar_d * t[index - 1]))

product = sp.eye(3)
ABSTRACT_PRODUCT_CHECKS = 0
visible_by_length = [sp.Integer(0)]
for length in range(1, MAX_STEPS + 1):
    product = matrices[length - 1] * product
    direct = product * source
    visible = -sum(
        c[exit_time] * t[exit_time]
        * sp.prod(u[tail] for tail in range(exit_time + 1, length))
        for exit_time in range(1, length)
    )
    require(sp.expand(direct[0] - visible) == 0,
            ("abstract visible convolution", length))
    require(sp.expand(direct[1] - t[length]) == 0,
            ("abstract transverse continuant", length))
    require(sp.expand(direct[2] - scalar_d * t[length - 1]) == 0,
            ("abstract returning coordinate", length))
    visible_by_length.append(sp.expand(direct[0]))
    ABSTRACT_PRODUCT_CHECKS += 1

# Eliminate the transverse state without dividing by an exit coefficient.
# This leaves a scalar polynomial-coefficient recurrence that remains valid
# at the isolated indices where one c-factor vanishes.
VISIBLE_RECURRENCE_CHECKS = 0
for index in range(1, MAX_STEPS - 1):
    left = c[index] * c[index - 1] * (
        u[index + 1] * visible_by_length[index + 1]
        - visible_by_length[index + 2])
    right = (
        c[index + 1] * c[index - 1] * alpha[index]
        * (u[index] * visible_by_length[index]
           - visible_by_length[index + 1])
        + c[index + 1] * c[index] * beta[index] * scalar_d
        * (u[index - 1] * visible_by_length[index - 1]
           - visible_by_length[index])
    )
    require(sp.expand(left - right) == 0,
            ("cleared visible recurrence", index))
    VISIBLE_RECURRENCE_CHECKS += 1

MATCHING_COUNTS = []
for length in range(MAX_STEPS):
    matching, count = path_continuant(alpha, beta, scalar_d, length)
    require(sp.expand(t[length + 1] - beta[0] * matching) == 0,
            ("matching-continuant expansion", length))
    MATCHING_COUNTS.append(count)
require(tuple(MATCHING_COUNTS) == (1, 1, 2, 3, 5, 8, 13),
        "path matching census")

require((matrices[0] * source)[0] == 0, "length-one boundary")
require(sp.expand((matrices[1] * matrices[0] * source)[0]
                  + c[1] * beta[0]) == 0,
        "length-two boundary")
require(sp.expand((matrices[2] * matrices[1] * matrices[0] * source)[0]
                  + beta[0] * (c[1] * u[2] + c[2] * alpha[1])) == 0,
        "length-three two-path formula")
transverse_step = sp.Matrix([[alpha[1], beta[1]], [scalar_d, 0]])
require(sp.factor(sp.det(transverse_step) + scalar_d * beta[1]) == 0,
        "transverse-step determinant")


# Factorial specialization.
n, d, v = sp.symbols("n d v")
Delta = 1 - 4 * d * v


def scalar_transfer(index):
    return sp.Matrix([
        [2 * (index + 1) * (2 * index + 1) * v,
         index * (index + 1) * Delta,
         d - index - 1],
        [1, 0, 0],
        [0, 0, d],
    ])


def factorial_data(index):
    a_i = 2 * (index + 1) * (2 * index + 1) * v
    b_i = index * (index + 1) * Delta
    c_i = d - index - 1
    return -b_i, c_i, a_i * d, b_i * d


# The factorial transverse continuant has a smaller normalized recurrence.
# Define G_r without division, so the identity remains meaningful on d=0.
NORMALIZED_CONTINUANT_CHECKS = 0
NORMALIZED_OGF_COEFFICIENT_CHECKS = 0
normalized_g = [sp.Integer(1)]
factorial_continuant = [sp.Integer(1)]
for length in range(1, MAX_STEPS + 1):
    previous_two = (normalized_g[length - 2]
                    if length >= 2 else sp.Integer(0))
    normalized_g.append(sp.expand(
        2 * v * (2 * n + 2 * length + 1) * normalized_g[length - 1]
        + Delta * previous_two))

    if length == 1:
        next_continuant = factorial_data(n + 1)[2]
    else:
        next_continuant = (
            factorial_data(n + length)[2]
            * factorial_continuant[length - 1]
            + d * factorial_data(n + length)[3]
            * factorial_continuant[length - 2]
        )
    factorial_continuant.append(sp.expand(next_continuant))
    expected = d**length * sp.rf(n + 2, length) * normalized_g[length]
    require(sp.expand(factorial_continuant[length] - expected) == 0,
            ("normalized factorial continuant", length))
    NORMALIZED_CONTINUANT_CHECKS += 1

# Coefficientwise audit of
# 4v z^2 G' + [2v(2n+3)z + Delta z^2 - 1]G + 1 = 0.
for degree in range(MAX_STEPS + 1):
    coefficient = -normalized_g[degree]
    if degree == 0:
        coefficient += 1
    if degree >= 1:
        coefficient += (4 * v * (degree - 1)
                        + 2 * v * (2 * n + 3)) \
            * normalized_g[degree - 1]
    if degree >= 2:
        coefficient += Delta * normalized_g[degree - 2]
    require(sp.expand(coefficient) == 0,
            ("normalized OGF differential equation", degree))
    NORMALIZED_OGF_COEFFICIENT_CHECKS += 1


index = sp.symbols("i", integer=True)
u_i, c_i, alpha_i, beta_i = factorial_data(index)
expected_exterior = sp.Matrix([
    [u_i, -c_i, 0],
    [0, alpha_i, beta_i],
    [0, d, 0],
])
require((compound(scalar_transfer(index)) - expected_exterior)
        .applyfunc(sp.expand) == sp.zeros(3),
        "factorial exterior reconstruction")

u_n, c_n, alpha_n, beta_n = factorial_data(n)
u_n1, c_n1, alpha_n1, beta_n1 = factorial_data(n + 1)
u_n2, c_n2, alpha_n2, beta_n2 = factorial_data(n + 2)
factorial_v2 = sp.factor(-c_n1 * beta_n)
require(sp.factor(factorial_v2
                  - (n + 2 - d) * n * (n + 1) * Delta * d) == 0,
        "THM-3183 two-step specialization")
factorial_v3 = sp.factor(
    -beta_n * (c_n1 * u_n2 + c_n2 * alpha_n1))
direct_factorial_v3 = (
    compound(scalar_transfer(n + 2))
    * compound(scalar_transfer(n + 1))
    * compound(scalar_transfer(n))
    * source
)[0]
require(sp.factor(direct_factorial_v3 - factorial_v3) == 0,
        "factorial length-three direct product")


# Sharp actual-family cancellation hostile.
HOSTILE = {n: 1, d: 5, v: sp.Rational(4, 105)}
hostile_delta = sp.factor(Delta.subs(HOSTILE))
hostile_beta = sp.factor(beta_n.subs(HOSTILE))
hostile_c1 = sp.factor(c_n1.subs(HOSTILE))
hostile_u2 = sp.factor(u_n2.subs(HOSTILE))
hostile_c2 = sp.factor(c_n2.subs(HOSTILE))
hostile_alpha1 = sp.factor(alpha_n1.subs(HOSTILE))
path_one = sp.factor((-beta_n * c_n1 * u_n2).subs(HOSTILE))
path_two = sp.factor((-beta_n * c_n2 * alpha_n1).subs(HOSTILE))
hostile_v2 = sp.factor(factorial_v2.subs(HOSTILE))
hostile_v3 = sp.factor(factorial_v3.subs(HOSTILE))
require(hostile_delta == sp.Rational(5, 21), "hostile discriminant")
require((hostile_beta, hostile_c1, hostile_u2,
         hostile_c2, hostile_alpha1)
        == (sp.Rational(50, 21), 2, sp.Rational(-20, 7),
            1, sp.Rational(40, 7)),
        "hostile local factors")
require((path_one, path_two)
        == (sp.Rational(2000, 147), sp.Rational(-2000, 147)),
        "hostile cancelling path values")
require(hostile_v2 == sp.Rational(-100, 21) and hostile_v3 == 0,
        "hostile visibility profile")
hostile_full_v3 = (
    compound(scalar_transfer(3))
    * compound(scalar_transfer(2))
    * compound(scalar_transfer(1))
    * source
).subs(HOSTILE)
require(hostile_full_v3[0] == 0
        and (hostile_full_v3[1] != 0 or hostile_full_v3[2] != 0),
        "hostile selected-chart rather than exterior-death boundary")


def valuation_rational(value, prime):
    value = sp.Rational(value)
    require(value != 0, "valuation of zero")
    numerator = abs(int(value.p))
    denominator = abs(int(value.q))
    answer = 0
    while numerator % prime == 0:
        numerator //= prime
        answer += 1
    while denominator % prime == 0:
        denominator //= prime
        answer -= 1
    return answer


local_weights = (hostile_beta, hostile_c1, hostile_u2,
                 hostile_c2, hostile_alpha1, path_one, path_two)
require(all(valuation_rational(value, 11) == 0 for value in local_weights),
        "hostile path weight lost 11-adic unit status")
hostile_scalar_matrices = tuple(
    scalar_transfer(step).subs(HOSTILE) for step in (1, 2, 3)
)
hostile_exterior_matrices = tuple(
    compound(matrix) for matrix in hostile_scalar_matrices
)
require(
    all(valuation_rational(entry, 11) >= 0
        for matrix in hostile_scalar_matrices + hostile_exterior_matrices
        for entry in matrix if entry != 0),
    "hostile transfer lost 11-adic integrality",
)
hostile_determinants = tuple(
    sp.factor(sp.det(matrix)) for matrix in hostile_scalar_matrices
)
require(all(valuation_rational(value, 11) == 0
            for value in hostile_determinants),
        "hostile local transfer lost 11-adic unimodularity")
require(all(valuation_rational(value**2, 11) == 0
            for value in hostile_determinants),
        "hostile exterior transfer lost 11-adic unimodularity")

# A same-local-Smith positive control is needed for the non-determination
# conclusion; unit local profiles alone are not a comparison.
POSITIVE = {n: 1, d: 5, v: 1}
positive_scalar_matrices = tuple(
    scalar_transfer(step).subs(POSITIVE) for step in (1, 2, 3)
)
positive_exterior_matrices = tuple(
    compound(matrix) for matrix in positive_scalar_matrices
)
require(
    all(valuation_rational(entry, 11) >= 0
        for matrix in positive_scalar_matrices + positive_exterior_matrices
        for entry in matrix if entry != 0),
    "positive control lost 11-adic integrality",
)
positive_determinants = tuple(
    sp.factor(sp.det(matrix)) for matrix in positive_scalar_matrices
)
require(all(valuation_rational(value, 11) == 0
            for value in positive_determinants),
        "positive control scalar Smith profile")
require(all(valuation_rational(value**2, 11) == 0
            for value in positive_determinants),
        "positive control exterior Smith profile")
positive_v3 = sp.factor(factorial_v3.subs(POSITIVE))
require(positive_v3 == 115140
        and valuation_rational(positive_v3, 11) == 0,
        "same-Smith positive visible control")


def support_pattern(matrix):
    return tuple(entry != 0 for entry in matrix)


hostile_scalar_support = tuple(
    support_pattern(matrix) for matrix in hostile_scalar_matrices
)
positive_scalar_support = tuple(
    support_pattern(matrix) for matrix in positive_scalar_matrices
)
hostile_exterior_support = tuple(
    support_pattern(matrix) for matrix in hostile_exterior_matrices
)
positive_exterior_support = tuple(
    support_pattern(matrix) for matrix in positive_exterior_matrices
)
require(hostile_scalar_support == positive_scalar_support,
        "scalar support-pattern mismatch")
require(hostile_exterior_support == positive_exterior_support,
        "exterior support-pattern mismatch")


print("THM-3186 FULL EXTERIOR CONTINUANT PATH CONVOLUTION EXACT CONTROL")
print("abstract_direct_product_checks=" + str(ABSTRACT_PRODUCT_CHECKS))
print("cleared_order3_visibility_recurrence_checks="
      + str(VISIBLE_RECURRENCE_CHECKS))
print("normalized_factorial_continuant_checks="
      + str(NORMALIZED_CONTINUANT_CHECKS))
print("normalized_ogf_differential_coefficient_checks="
      + str(NORMALIZED_OGF_COEFFICIENT_CHECKS))
print("continuant_matching_counts=" + repr(tuple(MATCHING_COUNTS)))
print("boundary_visible_coefficients=(0,-c_(n+1)beta_n)")
print("length3_exit_polynomial=c_(n+1)u_(n+2)+c_(n+2)alpha_(n+1)")
print("factorial_hostile=(n=1,d=5,v=4/105,Delta=5/21)")
print("hostile_path_contributions=(2000/147,-2000/147)")
print("hostile_visibility=(V2=-100/21,V3=0)")
print("hostile_V3_zero_but_transverse_chart_nonzero=PASS")
print("hostile_local_11adic_profiles=all_unimodular")
print("same_smith_positive_control=(n=1,d=5,v=1,V3=115140_unit)")
print("hostile_positive_scalar_exterior_support_patterns=identical")
print("scope=bare_scalar_exterior_tail_not_PRS_or_GMC_closure")
print("ALL EXACT CHECKS PASSED")
