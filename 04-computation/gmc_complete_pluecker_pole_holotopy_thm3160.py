#!/usr/bin/env python3
"""Exact controls for THM-3160's complete Pluecker pole holotopy.

The finite degree-two basis is ``(1,h1,m2,m11)``.  It already contains the
minimal same-degree selector projection failure, while the theorem proves the
general exterior-square identity formally.
"""

from itertools import product


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def matvec(matrix, vector):
    return tuple(sum(a * b for a, b in zip(row, vector)) for row in matrix)


def matmul(left, right):
    columns = tuple(zip(*right))
    return tuple(tuple(sum(a * b for a, b in zip(row, column))
                       for column in columns) for row in left)


def transpose(matrix):
    return tuple(zip(*matrix))


def wedge(left, right):
    return tuple(tuple(left[i] * right[j] - right[i] * left[j]
                       for j in range(len(left)))
                 for i in range(len(left)))


def pole_matrix(m):
    """Matrix of f[X] -> f[X-m] on (1,h1,m2,m11)."""

    return (
        (1, 0, 0, 0),
        (-m, 1, 0, 0),
        (-m * m, 0, 1, 0),
        (m * m, -m, 0, 1),
    )


def exterior_transport(matrix, omega):
    return matmul(matmul(matrix, omega), transpose(matrix))


ONE = (1, 0, 0, 0)
H1 = (0, 1, 0, 0)
M2 = (0, 0, 1, 0)
H2 = (0, 0, 1, 1)


def evaluate(vector, function):
    return sum(a * b for a, b in zip(vector, function))


def minor(left, right, first, second):
    return (evaluate(left, first) * evaluate(right, second)
            - evaluate(left, second) * evaluate(right, first))


def top_current(left, right):
    return minor(left, right, H2, M2)


# Pole squares commute, so prefix transport depends only on the multiset.
commuting_checks = 0
for m, n in product(range(-3, 4), repeat=2):
    require(matmul(pole_matrix(m), pole_matrix(n))
            == matmul(pole_matrix(n), pole_matrix(m)),
            f"pole square lost commutativity at {(m, n)}")
    commuting_checks += 1


# Directly verify Omega' = K Omega K^T on a deterministic endpoint bank.
endpoint_bank = tuple(
    (1, a, b, c)
    for a in range(-2, 3)
    for b in range(-1, 2)
    for c in range(-1, 2)
)
exterior_checks = 0
recurrence_checks = 0
for index in range(0, len(endpoint_bank), 7):
    left = endpoint_bank[index]
    right = endpoint_bank[-index - 1]
    for m in range(-3, 4):
        matrix = pole_matrix(m)
        transported_left = matvec(matrix, left)
        transported_right = matvec(matrix, right)
        require(wedge(transported_left, transported_right)
                == exterior_transport(matrix, wedge(left, right)),
                f"exterior-square transport failed at {(index, m)}")
        exterior_checks += 1

        # N=2 instance of the four-term recurrence (9).
        expected = (
            top_current(left, right)
            - m * m * minor(left, right, H2, ONE)
            - m * minor(left, right, H1, M2)
            + m * m * m * minor(left, right, H1, ONE)
        )
        require(top_current(transported_left, transported_right) == expected,
                f"top-coordinate recurrence failed at {(index, m)}")
        recurrence_checks += 1


def partitions(total, maximum=None):
    if total == 0:
        return ((),)
    if maximum is None or maximum > total:
        maximum = total
    answer = []
    for first in range(maximum, 0, -1):
        for tail in partitions(total - first, first):
            answer.append((first,) + tail)
    return tuple(answer)


def one_letter_monomial(shape, letter):
    return letter ** shape[0] if len(shape) == 1 else 0


def one_letter_h(degree, letter):
    return letter ** degree


# Every same-degree current of two one-letter endpoints vanishes.
parent_zero_checks = 0
for x, y in product(range(-3, 5), repeat=2):
    for degree in range(1, 9):
        for shape in partitions(degree):
            value = (
                one_letter_h(degree, x) * one_letter_monomial(shape, y)
                - one_letter_monomial(shape, x) * one_letter_h(degree, y)
            )
            require(value == 0,
                    f"one-letter parent current survived at {(x, y, shape)}")
            parent_zero_checks += 1


# The child degree-two top coordinate has the exact factorization (11).
factorization_checks = 0
for m, x, y in product(range(-3, 4), range(-2, 6), range(-2, 6)):
    left = (1, x, x * x, 0)
    right = (1, y, y * y, 0)
    child = top_current(matvec(pole_matrix(m), left),
                        matvec(pole_matrix(m), right))
    expected = m * (x - m) * (y - m) * (x - y)
    require(child == expected,
            f"one-letter child factorization failed at {(m, x, y)}")
    factorization_checks += 1


hostile_parent = top_current((1, 2, 4, 0), (1, 3, 9, 0))
control_parent = top_current((1, 2, 4, 0), (1, 2, 4, 0))
hostile_child = top_current(
    matvec(pole_matrix(1), (1, 2, 4, 0)),
    matvec(pole_matrix(1), (1, 3, 9, 0)),
)
control_child = top_current(
    matvec(pole_matrix(1), (1, 2, 4, 0)),
    matvec(pole_matrix(1), (1, 2, 4, 0)),
)
hostile_cross_minor = minor((1, 2, 4, 0), (1, 3, 9, 0), H1, ONE)
control_cross_minor = minor((1, 2, 4, 0), (1, 2, 4, 0), H1, ONE)
require((hostile_parent, control_parent) == (0, 0),
        "hostile parents left the common zero selector profile")
require((hostile_child, control_child) == (-2, 0),
        "hostile children lost selector separation")
require((hostile_cross_minor, control_cross_minor) == (-1, 0),
        "cross-degree missing coordinate drift")


def one_letter_child_top(degree, m, x, y):
    left_h = x ** (degree - 1) * (x - m)
    right_h = y ** (degree - 1) * (y - m)
    left_m = x ** degree - m ** degree
    right_m = y ** degree - m ** degree
    return left_h * right_m - left_m * right_h


hostile_child_sequence = tuple(
    one_letter_child_top(degree, 1, 2, 3) for degree in range(2, 9)
)
control_child_sequence = tuple(
    one_letter_child_top(degree, 1, 2, 2) for degree in range(2, 9)
)
require(hostile_child_sequence
        == (-2, -22, -170, -1150, -7322, -45262, -275690),
        "all-degree hostile child sequence drift")
require(control_child_sequence == (0,) * 7,
        "equal-endpoint child control became nonzero")


# No fixed degree-lag cap sees the degree-N/degree-zero minor uniformly.
bounded_lag_children = []
for lag in range(9):
    degree = lag + 2
    m = 2
    # A extracts m_(degree), B extracts 1.  The only nonzero minor has
    # degree gap `degree>lag`; the four-term recurrence gives this child.
    parent_top = 0
    omega_hn_one = 1
    omega_hprev_mn = 0
    omega_hprev_one = 0
    child_top = (parent_top - m ** degree * omega_hn_one
                 - m * omega_hprev_mn
                 + m ** (degree + 1) * omega_hprev_one)
    require(degree > lag and child_top == -(2 ** degree),
            f"bounded-lag hostile failed at {lag}")
    bounded_lag_children.append(child_top)


print("THM-3160 complete Pluecker pole holotopy")
print("degree_two_basis=(1,h1,m2,m11)")
print(f"commuting_prefix_square_checks={commuting_checks}")
print(f"exterior_square_transport_checks={exterior_checks}")
print(f"four_term_top_recurrence_checks={recurrence_checks}")
print(f"one_letter_parent_zero_checks_through_degree8={parent_zero_checks}")
print(f"one_letter_child_factorization_checks={factorization_checks}")
print("hostile_parent_top_and_control=(0,0)")
print("hostile_child_top_and_control=(-2,0)")
print("hostile_child_top_N2_N8=" + repr(hostile_child_sequence))
print("hostile_cross_degree_minor_and_control=(-1,0)")
print("bounded_lag_L0_L8_children=" + repr(tuple(bounded_lag_children)))
print("bifiltration_maps=(horizon_inclusion,depth_zero_extension)")
print("survivor=full_cross_degree_endpoint_Pluecker_tensor")
print("scope=universal_endpoint_transport_not_fixed_bank_positivity")
print("all_exact_checks=PASS")
