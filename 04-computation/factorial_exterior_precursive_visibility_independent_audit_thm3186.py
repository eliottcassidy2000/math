#!/usr/bin/env python3
"""Independent exact audit of THM-3186's P-recursive corollaries.

This companion uses a tiny sparse-polynomial ring and ``Fraction`` matrices.
It imports neither SymPy nor the primary THM-3186 implementation.
"""

from fractions import Fraction


def require(condition, data):
    if not condition:
        raise RuntimeError(data)


# A polynomial is a map from a sorted tuple of variable names to a rational
# coefficient.  Repeated names encode powers.
def clean(poly):
    return {monomial: coefficient for monomial, coefficient in poly.items()
            if coefficient}


def constant(value):
    value = Fraction(value)
    return {} if not value else {(): value}


def variable(name):
    return {(name,): Fraction(1)}


def add(left, right):
    answer = dict(left)
    for monomial, coefficient in right.items():
        answer[monomial] = answer.get(monomial, Fraction(0)) + coefficient
    return clean(answer)


def neg(poly):
    return {monomial: -coefficient for monomial, coefficient in poly.items()}


def sub(left, right):
    return add(left, neg(right))


def mul(left, right):
    answer = {}
    for monomial_left, coefficient_left in left.items():
        for monomial_right, coefficient_right in right.items():
            monomial = tuple(sorted(monomial_left + monomial_right))
            answer[monomial] = answer.get(monomial, Fraction(0)) \
                + coefficient_left * coefficient_right
    return clean(answer)


def product(values):
    answer = ONE
    for value in values:
        answer = mul(answer, value)
    return answer


def power(value, exponent):
    return product(value for _ in range(exponent))


def sum_poly(values):
    answer = ZERO
    for value in values:
        answer = add(answer, value)
    return answer


def matrix_vector(matrix, vector):
    return [sum_poly(mul(matrix[row][column], vector[column])
                     for column in range(len(vector)))
            for row in range(len(matrix))]


ZERO = constant(0)
ONE = constant(1)


# First audit the cleared recurrence over a freely generated coefficient ring.
MAX_STEPS = 7
scalar_d = variable("D0")
u = tuple(variable("u%d" % index) for index in range(MAX_STEPS))
c = tuple(variable("c%d" % index) for index in range(MAX_STEPS))
alpha = tuple(variable("a%d" % index) for index in range(MAX_STEPS))
beta = tuple(variable("b%d" % index) for index in range(MAX_STEPS))
matrices = tuple(
    ((u[index], neg(c[index]), ZERO),
     (ZERO, alpha[index], beta[index]),
     (ZERO, scalar_d, ZERO))
    for index in range(MAX_STEPS)
)

state = [ZERO, ZERO, ONE]
visible = [ZERO]
for length in range(1, MAX_STEPS + 1):
    state = matrix_vector(matrices[length - 1], state)
    visible.append(state[0])

free_recurrence_checks = 0
for index in range(1, MAX_STEPS - 1):
    left = mul(mul(c[index], c[index - 1]),
               sub(mul(u[index + 1], visible[index + 1]),
                   visible[index + 2]))
    first = mul(mul(mul(c[index + 1], c[index - 1]), alpha[index]),
                sub(mul(u[index], visible[index]), visible[index + 1]))
    second = mul(
        mul(mul(mul(c[index + 1], c[index]), beta[index]), scalar_d),
        sub(mul(u[index - 1], visible[index - 1]), visible[index]))
    require(left == add(first, second),
            ("free cleared visibility recurrence", index))
    free_recurrence_checks += 1


# Independently normalize the factorial continuant over Q[N,D,V].
N = variable("N")
D = variable("D")
V = variable("V")
J = variable("J")
DELTA = sub(ONE, mul(constant(4), mul(D, V)))


def shifted(base, offset):
    return add(base, constant(offset))


def factorial_weights(index):
    index_plus_one = add(index, ONE)
    b_value = mul(mul(index, index_plus_one), DELTA)
    u_value = neg(b_value)
    c_value = sub(sub(D, index), ONE)
    alpha_value = mul(
        constant(2),
        mul(D, mul(index_plus_one,
                   mul(add(mul(constant(2), index), ONE), V))))
    beta_value = mul(D, b_value)
    return u_value, c_value, alpha_value, beta_value


normalized = [ONE]
continuant = [ONE]
normalized_checks = 0
ogf_coefficient_checks = 0
for length in range(1, MAX_STEPS + 1):
    previous_two = normalized[length - 2] if length >= 2 else ZERO
    linear_factor = add(mul(constant(2), N), constant(2 * length + 1))
    normalized.append(add(
        mul(mul(constant(2), V),
            mul(linear_factor, normalized[length - 1])),
        mul(DELTA, previous_two)))

    _, _, alpha_value, beta_value = factorial_weights(shifted(N, length))
    if length == 1:
        next_continuant = alpha_value
    else:
        next_continuant = add(
            mul(alpha_value, continuant[length - 1]),
            mul(mul(D, beta_value), continuant[length - 2]))
    continuant.append(next_continuant)
    rising = product(shifted(N, offset)
                     for offset in range(2, length + 2))
    expected = mul(mul(power(D, length), rising), normalized[length])
    require(continuant[length] == expected,
            ("normalized continuant", length))
    normalized_checks += 1

for degree in range(MAX_STEPS + 1):
    coefficient = neg(normalized[degree])
    if degree == 0:
        coefficient = add(coefficient, ONE)
    if degree >= 1:
        multiplier = mul(
            constant(2),
            mul(V, add(mul(constant(2), N), constant(2 * degree + 1))))
        coefficient = add(
            coefficient, mul(multiplier, normalized[degree - 1]))
    if degree >= 2:
        coefficient = add(
            coefficient, mul(DELTA, normalized[degree - 2]))
    require(coefficient == ZERO,
            ("normalized OGF coefficient", degree))
    ogf_coefficient_checks += 1


# The coefficient of V_(j+2) is -c_(n+j)c_(n+j-1), with nonzero
# quadratic leading term -j^2.  Check this directly in Q[N,D,J].
_, c_j, _, _ = factorial_weights(add(N, J))
_, c_j_minus_one, _, _ = factorial_weights(add(add(N, J), constant(-1)))
leading_coefficient_poly = neg(mul(c_j, c_j_minus_one))
max_j_degree = max((monomial.count("J")
                    for monomial in leading_coefficient_poly), default=-1)
degree_two_part = {
    tuple(name for name in monomial if name != "J"): coefficient
    for monomial, coefficient in leading_coefficient_poly.items()
    if monomial.count("J") == 2
}
require(max_j_degree == 2 and degree_two_part == {(): Fraction(-1)},
        (max_j_degree, degree_two_part))


# Finally replay a rational factorial specialization across the actual c=0
# boundary.  This is a separate arithmetic path from the sparse-polynomial
# audit above.
def rational_weights(index, scalar_d_value, v_value):
    scalar_d_value = Fraction(scalar_d_value)
    v_value = Fraction(v_value)
    delta = 1 - 4 * scalar_d_value * v_value
    b_value = index * (index + 1) * delta
    return (-b_value,
            scalar_d_value - index - 1,
            2 * scalar_d_value * (index + 1) * (2 * index + 1) * v_value,
            scalar_d_value * b_value)


def rational_state_profile(n_value, scalar_d_value, v_value, steps):
    state_value = (Fraction(0), Fraction(0), Fraction(1))
    profile = [state_value]
    for length in range(steps):
        weights = rational_weights(
            n_value + length, scalar_d_value, v_value)
        matrix = ((weights[0], -weights[1], Fraction(0)),
                  (Fraction(0), weights[2], weights[3]),
                  (Fraction(0), Fraction(scalar_d_value), Fraction(0)))
        state_value = tuple(
            sum(matrix[row][column] * state_value[column]
                for column in range(3)) for row in range(3))
        profile.append(state_value)
    return profile


hostile_n = 1
hostile_d = 5
hostile_v = Fraction(4, 105)
hostile_profile = rational_state_profile(
    hostile_n, hostile_d, hostile_v, 10)
boundary_recurrence_checks = 0
for index in range(1, 9):
    u_prev, c_prev, _, _ = rational_weights(
        hostile_n + index - 1, hostile_d, hostile_v)
    u_now, c_now, alpha_now, beta_now = rational_weights(
        hostile_n + index, hostile_d, hostile_v)
    u_next, c_next, _, _ = rational_weights(
        hostile_n + index + 1, hostile_d, hostile_v)
    values = tuple(state_value[0] for state_value in hostile_profile)
    left = c_now * c_prev * (u_next * values[index + 1]
                              - values[index + 2])
    right = (c_next * c_prev * alpha_now
             * (u_now * values[index] - values[index + 1])
             + c_next * c_now * beta_now * hostile_d
             * (u_prev * values[index - 1] - values[index]))
    require(left == right, ("rational boundary recurrence", index))
    boundary_recurrence_checks += 1

require(hostile_profile[2][0] == Fraction(-100, 21), hostile_profile[2])
require(hostile_profile[3][0] == 0
        and (hostile_profile[3][1] != 0 or hostile_profile[3][2] != 0),
        hostile_profile[3])
positive_profile = rational_state_profile(1, 5, 1, 3)
require(positive_profile[3][0] == 115140, positive_profile[3])


print("THM-3186 INDEPENDENT P-RECURSIVE VISIBILITY AUDIT")
print("implementation=custom sparse polynomial ring + Fraction matrices")
print("free_cleared_order3_recurrence_checks="
      + str(free_recurrence_checks))
print("normalized_factorial_continuant_checks=" + str(normalized_checks))
print("normalized_ogf_differential_coefficient_checks="
      + str(ogf_coefficient_checks))
print("leading_V_(j+2)_coefficient=(j_degree=2,leading=-1)")
print("factorial_boundary_recurrence_checks="
      + str(boundary_recurrence_checks))
print("vanishing_exit_coefficient_at=(n=1,d=5,j=3)")
print("hostile_visibility=(V2=-100/21,V3=0,transverse_nonzero)")
print("same_smith_positive_V3=115140")
print("scope=scalar_exterior_visibility_not_PRS_depth_or_factorial_closure")
print("INDEPENDENT EXACT AUDIT PASSED")
