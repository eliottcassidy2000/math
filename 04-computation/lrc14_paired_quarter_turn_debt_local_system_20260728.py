#!/usr/bin/env python3
"""Exact paired quarter-turn audit for the recurrent LRC14 debt carriers.

The BABA and relative-Segal four-cycles both have unsigned incidence
u+u^{-1}, which becomes multiplication by two after reflection quotient.
This companion verifies the minimal integral quarter-turn local system and
the matching mod-17 phase / mod-13 coefficient shadows.  It does not build a
nonnegative semantic current or a physical coefficient transport.
"""

from __future__ import annotations

from math import lcm


P = 13


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def matrix_mul(left, right):
    return tuple(
        tuple(sum(left[i][k] * right[k][j]
                  for k in range(len(right)))
              for j in range(len(right[0])))
        for i in range(len(left))
    )


def matrix_add(left, right):
    return tuple(
        tuple(x + y for x, y in zip(left_row, right_row))
        for left_row, right_row in zip(left, right)
    )


def poly_mul(left, right, modulus=None):
    result = [0] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            result[i + j] += x * y
    if modulus is not None:
        result = [value % modulus for value in result]
    return tuple(result)


def phi7_mul(left, right):
    result = list(poly_mul(left, right, P))
    result.extend([0] * (11 - len(result)))
    for degree in range(10, 5, -1):
        coefficient = result[degree] % P
        for index in range(6):
            result[degree - 6 + index] = (
                result[degree - 6 + index] - coefficient
            ) % P
    return tuple(result[:6])


def phi7_power(value, exponent):
    result = (1, 0, 0, 0, 0, 0)
    while exponent:
        if exponent & 1:
            result = phi7_mul(result, value)
        value = phi7_mul(value, value)
        exponent //= 2
    return result


def quadratic_mul(left, right, linear_coefficient):
    """Multiply modulo z^2+linear_coefficient*z+1 over F_13."""
    constant = (left[0] * right[0]) % P
    linear = (left[0] * right[1] + left[1] * right[0]) % P
    quadratic = (left[1] * right[1]) % P
    return (
        (constant - quadratic) % P,
        (linear - linear_coefficient * quadratic) % P,
    )


def quadratic_reduce(value, linear_coefficient):
    result = (0, 0)
    for coefficient in reversed(value):
        result = quadratic_mul(
            result, (0, 1), linear_coefficient
        )
        result = ((result[0] + coefficient) % P, result[1])
    return result


def quadratic_power(value, exponent, linear_coefficient):
    result = (1, 0)
    while exponent:
        if exponent & 1:
            result = quadratic_mul(
                result, value, linear_coefficient
            )
        value = quadratic_mul(value, value, linear_coefficient)
        exponent //= 2
    return result


def multiplicative_order(value, linear_coefficient, bound=168):
    for exponent in range(1, bound + 1):
        if quadratic_power(value, exponent, linear_coefficient) == (1, 0):
            return exponent
    raise RuntimeError("quadratic component order exceeded declared bound")


def matrix_rank_mod_p(rows, modulus):
    """Row rank over the prime field F_modulus."""
    work = [list(value % modulus for value in row) for row in rows]
    if not work:
        return 0
    row = 0
    for column in range(len(work[0])):
        pivot = next(
            (index for index in range(row, len(work))
             if work[index][column]),
            None,
        )
        if pivot is None:
            continue
        work[row], work[pivot] = work[pivot], work[row]
        inverse = pow(work[row][column], -1, modulus)
        work[row] = [(inverse * value) % modulus for value in work[row]]
        for index in range(len(work)):
            if index == row:
                continue
            coefficient = work[index][column]
            if coefficient:
                work[index] = [
                    (left - coefficient * right) % modulus
                    for left, right in zip(work[index], work[row])
                ]
        row += 1
        if row == len(work):
            break
    return row


def main():
    zero2 = ((0, 0), (0, 0))
    identity2 = ((1, 0), (0, 1))
    minus_identity2 = ((-1, 0), (0, -1))
    quarter_turn = ((0, -1), (1, 0))
    inverse_quarter_turn = ((0, 1), (-1, 0))
    require(matrix_mul(quarter_turn, quarter_turn) == minus_identity2,
            "integral quarter-turn stopped squaring to -I")
    require(matrix_mul(quarter_turn, inverse_quarter_turn) == identity2,
            "displayed inverse quarter-turn changed")
    require(matrix_add(quarter_turn, inverse_quarter_turn) == zero2,
            "symmetric debt did not cancel in the rank-two local system")
    rank_one_sums = tuple(sign + sign for sign in (-1, 1))
    require(rank_one_sums == (-2, 2),
            "rank-one integral hostile controls changed")

    phase_orbit = [4]
    for _ in range(4):
        phase_orbit.append(13 * phase_orbit[-1] % 17)
    require(tuple(phase_orbit) == (4, 1, 13, 16, 4),
            "C17 quarter-turn orbit changed")
    require(pow(13, 2, 17) == 16 and pow(13, 4, 17) == 1
            and (13 + pow(13, -1, 17)) % 17 == 0,
            "mod-17 quarter-turn identities changed")
    endpoints = (4, 13)
    ghosts = (1, 16)
    incidence = tuple(
        tuple(
            sum(
                ghost == pow(13, exponent, 17) * endpoint % 17
                for exponent in (1, 3)
            )
            for endpoint in endpoints
        )
        for ghost in ghosts
    )
    require(incidence == ((1, 1), (1, 1)),
            "endpoint-to-ghost incidence changed")

    factors = ((1, 3, 1), (1, 5, 1), (1, 6, 1))
    factor_product = poly_mul(poly_mul(factors[0], factors[1], P),
                              factors[2], P)
    require(factor_product == (1, 1, 1, 1, 1, 1, 1),
            "Phi7 factorization over F13 changed")
    quadratic_discriminants = tuple(
        (factor[1] * factor[1] - 4) % P for factor in factors
    )
    quadratic_residues = {value * value % P for value in range(P)}
    require(quadratic_discriminants == (5, 8, 6)
            and all(value not in quadratic_residues
                    for value in quadratic_discriminants),
            "a displayed Phi7 quadratic stopped being irreducible")
    unit_rows = (
        (11, 0, 0, 0, 0, 6),
        (0, 0, 0, 5, 0, 0),
        (11, 0, 0, 0, 1, 2),
        (0, 0, 0, 8, 0, 0),
    )
    unit_product = (1, 0, 0, 0, 0, 0)
    for row in unit_rows:
        unit_product = phi7_mul(unit_product, row)
    require(unit_product == (9, 2, 8, 7, 6, 9),
            "BABA conditional unit product changed")

    component_rows = []
    for linear_coefficient in (3, 5, 6):
        component = quadratic_reduce(unit_product, linear_coefficient)
        order = multiplicative_order(component, linear_coefficient)
        power42 = quadratic_power(component, 42, linear_coefficient)
        inverse42 = quadratic_power(
            component, order - (42 % order), linear_coefficient
        )
        symmetric42 = tuple(
            (x + y) % P for x, y in zip(power42, inverse42)
        )
        component_rows.append(
            (linear_coefficient, component, order, power42, symmetric42)
        )
    component_rows = tuple(component_rows)
    require(component_rows == (
        (3, (7, 0), 12, (12, 0), (11, 0)),
        (5, (4, 5), 168, (8, 0), (0, 0)),
        (6, (5, 4), 28, (12, 0), (11, 0)),
    ), "conditional unit component decomposition changed")
    require(lcm(*(row[2] for row in component_rows)) == 168
            and phi7_power(unit_product, 168)
            == (1, 0, 0, 0, 0, 0),
            "global conditional unit order changed")
    global_power42 = phi7_power(unit_product, 42)
    global_inverse42 = phi7_power(unit_product, 168 - 42)
    global_symmetric42 = tuple(
        (left + right) % P
        for left, right in zip(global_power42, global_inverse42)
    )
    symmetric_components = tuple(
        quadratic_reduce(global_symmetric42, coefficient)
        for coefficient in (3, 5, 6)
    )
    require(global_power42 == (8, 0, 10, 8, 8, 10)
            and global_inverse42 == (5, 0, 11, 1, 1, 11)
            and global_symmetric42 == (0, 0, 8, 9, 9, 8),
            "global symmetric U^42 debt changed")
    require(symmetric_components == ((11, 0), (0, 0), (11, 0)),
            "CRT image of the global symmetric debt changed")
    middle_basis_images = tuple(
        quadratic_reduce(
            tuple(1 if index == basis else 0 for index in range(6)), 5
        )
        for basis in range(6)
    )
    middle_projection_matrix = tuple(
        tuple(image[coordinate] for image in middle_basis_images)
        for coordinate in range(2)
    )
    middle_projection_rank = matrix_rank_mod_p(
        middle_projection_matrix, P
    )
    middle_kernel_dimension = 6 - middle_projection_rank
    require(middle_projection_rank == 2 and middle_kernel_dimension == 4,
            "middle CRT projection dimension ledger changed")
    require(global_symmetric42 != (0, 0, 0, 0, 0, 0)
            and quadratic_reduce(global_symmetric42, 5) == (0, 0),
            "middle projection no longer kills exactly the displayed debt")
    require(pow(8, 2, 13) == 12
            and (8 + pow(8, -1, 13)) % 13 == 0,
            "mod-13 coefficient quarter-turn changed")
    roots13 = tuple(value for value in range(13)
                    if value * value % 13 == 12)
    roots17 = tuple(value for value in range(17)
                    if value * value % 17 == 16)
    roots221 = tuple(value for value in range(221)
                     if value * value % 221 == 220)
    joint_root = tuple(
        value for value in roots221
        if value % 13 == 8 and value % 17 == 13
    )
    require(roots13 == (5, 8) and roots17 == (4, 13)
            and roots221 == (21, 47, 174, 200)
            and joint_root == (47,),
            "CRT quarter-turn alignment changed")

    print("LRC14 PAIRED QUARTER-TURN DEBT LOCAL-SYSTEM AUDIT")
    print("status=PROVED-ELEMENTARY ALGEBRA + VERIFIED-EXACT CONNECTION")
    print(f"integral_quarter_turn={quarter_turn} square={minus_identity2} "
          f"symmetric_sum={zero2}")
    print(f"rank_one_GL1Z_symmetric_sums={rank_one_sums}")
    print(f"C17_phase_orbit={tuple(phase_orbit)} "
          "i17=13 i17^2=-1 inverse_sum=0")
    print(f"endpoint_ghost_unsigned_incidence={incidence} "
          "reflection_quotient_multiplier=2")
    print(f"Phi7_F13_factors={factors}")
    print(f"conditional_BABA_unit_product={unit_product}")
    print(f"unit_component_data={component_rows}")
    print("Phi7_CRT_product=A=K3xK5xK6 component_dimensions=(2,2,2) "
          f"middle_projection_rank={middle_projection_rank} "
          f"kernel_dimension={middle_kernel_dimension}")
    print(f"global_U42={global_power42} "
          f"global_Uminus42={global_inverse42}")
    print(f"global_symmetric_debt={global_symmetric42} "
          f"CRT_components={symmetric_components}")
    print("middle_component_quarter_turn=i13=8 i13^2=-1 "
          "inverse_sum=0")
    print(f"integral_J_reductions=roots13={roots13} roots17={roots17} "
          f"roots221={roots221} selected_i221={joint_root[0]}")
    print("loss_ledger=quarter cancellation uses signed/local coefficients; "
          "middle projection has a 4-dimensional K3xK6 kernel and kills "
          "the nonzero global debt; U is a conditional four-row product "
          "and U^42 is cycle holonomy, not an edgewise gain; the four unit "
          "rows have no proved transport")
    print("SCOPE: no nonnegative current, semantic cospan, row exclusion, "
          "or LRC14 conclusion")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
