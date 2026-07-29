#!/usr/bin/env python3
"""Exact group-theoretic audit for the THM-2859 V4/carry hinge.

THM-2851 proves that the natural residue lifts live on the nonsplit cyclic
extension C_169, not on a split F_13^2 address plane.  This companion checks
the sharp realization distinction:

* the natural lift L_1 has order 169 and L_1^13 is the nontrivial carry T;
* every lift of the residue generator to C_169 has order 169, so the
  extension has no section;
* every p-element of AGL_2(F_13) has exponent 13, so an affine action on the
  split 169-point plane cannot realize L_1;
* the standard q3/q11/q7 triangle has exactly one carry; and
* the formal V4 degrees supplied by the horn diagonal and even collar add to
  the missing degree, but degree addition alone is not arrow composition.

No Python ``assert`` statement is used.
"""

from collections import Counter


P = 13
P2 = P * P


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def pair_add(left, right):
    return ((left[0] + right[0]) % P, (left[1] + right[1]) % P)


def matrix_multiply(left, right):
    return (
        (
            (left[0][0] * right[0][0]
             + left[0][1] * right[1][0]) % P,
            (left[0][0] * right[0][1]
             + left[0][1] * right[1][1]) % P,
        ),
        (
            (left[1][0] * right[0][0]
             + left[1][1] * right[1][0]) % P,
            (left[1][0] * right[0][1]
             + left[1][1] * right[1][1]) % P,
        ),
    )


IDENTITY = ((1, 0), (0, 1))


def matrix_power(matrix, exponent):
    result = IDENTITY
    base = matrix
    power = exponent
    while power:
        if power & 1:
            result = matrix_multiply(result, base)
        base = matrix_multiply(base, base)
        power //= 2
    return result


def matrix_vector(matrix, vector):
    return (
        (matrix[0][0] * vector[0] + matrix[0][1] * vector[1]) % P,
        (matrix[1][0] * vector[0] + matrix[1][1] * vector[1]) % P,
    )


def matrix_add(left, right):
    return (
        (
            (left[0][0] + right[0][0]) % P,
            (left[0][1] + right[0][1]) % P,
        ),
        (
            (left[1][0] + right[1][0]) % P,
            (left[1][1] + right[1][1]) % P,
        ),
    )


ZERO_MATRIX = ((0, 0), (0, 0))


def matrix_geometric_sum(matrix, terms):
    total = ZERO_MATRIX
    power = IDENTITY
    for _index in range(terms):
        total = matrix_add(total, power)
        power = matrix_multiply(power, matrix)
    return total


def determinant(matrix):
    return (
        matrix[0][0] * matrix[1][1]
        - matrix[0][1] * matrix[1][0]
    ) % P


def affine_power(linear, shift, exponent):
    result_linear = IDENTITY
    result_shift = (0, 0)
    base_linear = linear
    base_shift = shift
    power = exponent
    while power:
        if power & 1:
            result_shift = pair_add(
                matrix_vector(result_linear, base_shift),
                result_shift,
            )
            result_linear = matrix_multiply(
                result_linear, base_linear
            )
        doubled_shift = pair_add(
            matrix_vector(base_linear, base_shift),
            base_shift,
        )
        base_linear = matrix_multiply(base_linear, base_linear)
        base_shift = doubled_shift
        power //= 2
    return result_linear, result_shift


def digit_pair(number):
    return number // P, number % P


def pair_number(pair):
    ancestry, residue = pair
    return (P * ancestry + residue) % P2


def natural_lift(pair, step):
    ancestry, residue = pair
    return (
        (ancestry + (residue + step) // P) % P,
        (residue + step) % P,
    )


def carry(pair):
    return ((pair[0] + 1) % P, pair[1])


def iterate(function, value, times):
    answer = value
    for _index in range(times):
        answer = function(answer)
    return answer


def permutation_order(permutation):
    seen = set()
    order = 1

    def gcd(left, right):
        while right:
            left, right = right, left % right
        return left

    for start in range(len(permutation)):
        if start in seen:
            continue
        cursor = start
        length = 0
        while cursor not in seen:
            seen.add(cursor)
            cursor = permutation[cursor]
            length += 1
        order = order * length // gcd(order, length)
    return order


def main():
    states = tuple((ancestry, residue)
                   for ancestry in range(P) for residue in range(P))

    carry_cocycle_checks = 0
    for pair in states:
        for first in range(P):
            for second in range(P):
                left = natural_lift(
                    natural_lift(pair, first), second
                )
                epsilon = (first + second) // P
                right = natural_lift(
                    pair, (first + second) % P
                )
                if epsilon:
                    right = carry(right)
                require(left == right, "carry cocycle identity failed")
                carry_cocycle_checks += 1

    l1_permutation = tuple(
        pair_number(natural_lift(digit_pair(number), 1))
        for number in range(P2)
    )
    t_permutation = tuple(
        pair_number(carry(digit_pair(number)))
        for number in range(P2)
    )
    require(
        permutation_order(l1_permutation) == P2,
        "natural lifted generator lost order 169",
    )
    require(
        permutation_order(t_permutation) == P,
        "carry lost order 13",
    )
    for pair in states:
        require(
            iterate(lambda value: natural_lift(value, 1), pair, P)
            == carry(pair),
            "L1^13 stopped equalling T",
        )

    lift_orders = Counter()
    for high_digit in range(P):
        lift = 1 + P * high_digit
        cursor = 0
        order = 0
        while True:
            cursor = (cursor + lift) % P2
            order += 1
            if cursor == 0:
                break
        lift_orders[order] += 1
    require(
        lift_orders == Counter({P2: P}),
        "a residue-generator lift unexpectedly split C169",
    )

    matrices = tuple(
        ((a, b), (c, d))
        for a in range(P)
        for b in range(P)
        for c in range(P)
        for d in range(P)
        if (a * d - b * c) % P
    )
    require(
        len(matrices) == (P**2 - 1) * (P**2 - P),
        "GL2(F13) census changed",
    )
    p_power_linear = tuple(
        matrix for matrix in matrices
        if matrix_power(matrix, P2) == IDENTITY
    )
    require(
        len(p_power_linear) == P**2,
        "GL2(F13) p-power candidate census changed",
    )
    linear_orders = Counter()
    affine_candidate_count = 0
    affine_killed_by_p = 0
    for matrix in p_power_linear:
        linear_order = (
            1 if matrix == IDENTITY
            else P if matrix_power(matrix, P) == IDENTITY
            else P2
        )
        linear_orders[linear_order] += 1
        geometric = matrix_geometric_sum(matrix, P)
        require(
            geometric == ZERO_MATRIX,
            "a p-power GL2 candidate acquired nonzero p-step trace",
        )
        for shift in ((x, y) for x in range(P) for y in range(P)):
            affine_candidate_count += 1
            power = affine_power(matrix, shift, P)
            if power == (IDENTITY, (0, 0)):
                affine_killed_by_p += 1
    require(
        affine_candidate_count == P**4
        and affine_killed_by_p == affine_candidate_count,
        "AGL2(F13) acquired an affine element of order 169",
    )

    q3 = (0, 3)
    q11 = natural_lift(q3, 8)
    q7_via_q11 = natural_lift(q11, 9)
    q7_direct = natural_lift(q3, 4)
    require(
        (q11, q7_via_q11, q7_direct)
        == ((0, 11), (1, 7), (0, 7)),
        "q3/q11/q7 carry triangle changed",
    )
    require(
        q7_via_q11 == carry(q7_direct),
        "q3/q11/q7 triangle lost its central carry",
    )

    # Exact physical half-step hinge singled out in the common horn cells.
    half_step = 401080680
    root = (142004190428100, 142004216872980)
    m1 = (142004591508780, 142004617953660)
    m2 = (142004992589460, 142005019034340)
    require(
        tuple(right - left for left, right in zip(root, m1))
        == (half_step, half_step)
        and tuple(right - left for left, right in zip(root, m2))
        == (2 * half_step, 2 * half_step),
        "distinguished collar hinge changed scale",
    )

    horn_degree = (1, 1)
    even_collar_degree = (0, 1)
    formal_missing_degree = (
        horn_degree[0] ^ even_collar_degree[0],
        horn_degree[1] ^ even_collar_degree[1],
    )
    require(
        formal_missing_degree == (1, 0),
        "horn/collar degree sum stopped selecting the missing V4 degree",
    )

    print("THM-2859 V4 / WITT-CARRY HINGE AUDIT")
    print(
        f"carry_cocycle_checks={carry_cocycle_checks};"
        f"L1_order={permutation_order(l1_permutation)};"
        f"T_order={permutation_order(t_permutation)};L1^13=T"
    )
    print(
        f"C169_residue_generator_lift_orders={tuple(sorted(lift_orders.items()))};"
        "section=no"
    )
    print(
        "split_plane="
        f"|GL2(F13)|:{len(matrices)},"
        f"p_power_linear_candidates:{len(p_power_linear)},"
        f"linear_orders:{tuple(sorted(linear_orders.items()))},"
        f"affine_candidates:{affine_candidate_count},"
        f"killed_by_13:{affine_killed_by_p};"
        "AGL2(F13)_order169=no"
    )
    print(
        "mechanism=an affine p-element embeds as a unipotent 3x3 matrix;"
        "for p=13 its nilpotence index is <13, hence its exponent is 13"
    )
    print(
        "sharp_realization=C169=(W2(F13),+) has the same 169 underlying "
        "states as F13^2 but exponent169 rather than exponent13"
    )
    print(
        f"q3_q11_q7={q3}->{q11}->{q7_via_q11};"
        f"direct={q7_direct};two_step=T*direct"
    )
    print(
        f"distinguished_hinge=R:{root}->M1:{m1}->M2:{m2};"
        f"h={half_step};"
        "the 20-label collapse is an inherited concurrent-scout fact, "
        "not reconstructed by this script"
    )
    print(
        f"formal_V4_degree=horn:{horn_degree}+"
        f"even_collar:{even_collar_degree}={formal_missing_degree}"
    )
    print(
        "composition_boundary=the horn q3/q11/q7 atoms and the collar q0 "
        "triple are not a common object; adding degree labels is not "
        "composition.  A faithful comparison must first lift them to the "
        "same nonsplit C169 carry torsor."
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
