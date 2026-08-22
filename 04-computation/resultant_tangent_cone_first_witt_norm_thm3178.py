#!/usr/bin/env python3
"""Exact controls for THM-3178's resultant tangent-cone formula."""

from fractions import Fraction


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def trim(polynomial):
    answer = list(polynomial)
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


def add(left, right):
    answer = [0] * max(len(left), len(right))
    for index, coefficient in enumerate(left):
        answer[index] += coefficient
    for index, coefficient in enumerate(right):
        answer[index] += coefficient
    return trim(answer)


def scale(polynomial, scalar):
    return trim([scalar * coefficient for coefficient in polynomial])


def multiply(left, right):
    answer = [0] * (len(left) + len(right) - 1)
    for left_index, left_coefficient in enumerate(left):
        for right_index, right_coefficient in enumerate(right):
            answer[left_index + right_index] += left_coefficient * right_coefficient
    return trim(answer)


def derivative(polynomial):
    if len(polynomial) == 1:
        return [0]
    return [index * coefficient for index, coefficient in enumerate(polynomial)][1:]


def evaluate(polynomial, value):
    answer = 0
    for coefficient in reversed(polynomial):
        answer = answer * value + coefficient
    return answer


def remainder_monic(dividend, divisor):
    """Remainder over Z for a monic divisor."""
    divisor = trim(divisor)
    require(divisor[-1] == 1, "divisor is not monic")
    answer = trim(dividend)
    while len(answer) >= len(divisor):
        coefficient = answer[-1]
        shift = len(answer) - len(divisor)
        for index, entry in enumerate(divisor):
            answer[index + shift] -= coefficient * entry
        answer = trim(answer)
    return answer


def bareiss_determinant(matrix):
    if not matrix:
        return 1
    work = [row[:] for row in matrix]
    sign = 1
    previous = 1
    for column in range(len(work) - 1):
        pivot_row = next(
            (row for row in range(column, len(work)) if work[row][column]),
            None,
        )
        if pivot_row is None:
            return 0
        if pivot_row != column:
            work[column], work[pivot_row] = work[pivot_row], work[column]
            sign = -sign
        pivot = work[column][column]
        for row in range(column + 1, len(work)):
            for entry in range(column + 1, len(work)):
                numerator = (
                    work[row][entry] * pivot
                    - work[row][column] * work[column][entry]
                )
                require(numerator % previous == 0, "Bareiss division failed")
                work[row][entry] = numerator // previous
            work[row][column] = 0
        previous = pivot
    return sign * work[-1][-1]


def resultant(left, right):
    left = trim(left)
    right = trim(right)
    m = len(left) - 1
    n = len(right) - 1
    require(not (m == 0 and left[0] == 0), "zero left polynomial")
    require(not (n == 0 and right[0] == 0), "zero right polynomial")
    left_descending = list(reversed(left))
    right_descending = list(reversed(right))
    matrix = []
    for shift in range(n):
        matrix.append([0] * shift + left_descending + [0] * (n - 1 - shift))
    for shift in range(m):
        matrix.append([0] * shift + right_descending + [0] * (m - 1 - shift))
    return bareiss_determinant(matrix)


def binomial_basis(maximum):
    """Return coefficient arrays of binom(x,k), 0<=k<=maximum."""
    answer = [[Fraction(1)]]
    for k in range(1, maximum + 1):
        factor = [Fraction(-(k - 1), k), Fraction(1, k)]
        answer.append(multiply(answer[-1], factor))
    return answer


def interpolate(values):
    """Interpolate P(0),...,P(D) in the Newton forward basis."""
    differences = [Fraction(value) for value in values]
    leading_differences = []
    while differences:
        leading_differences.append(differences[0])
        differences = [
            differences[index + 1] - differences[index]
            for index in range(len(differences) - 1)
        ]
    basis = binomial_basis(len(values) - 1)
    polynomial = [Fraction(0)] * len(values)
    for difference, term in zip(leading_differences, basis):
        for index, coefficient in enumerate(term):
            polynomial[index] += difference * coefficient
    require(all(coefficient.denominator == 1 for coefficient in polynomial),
            "nonintegral interpolation")
    return trim([int(coefficient) for coefficient in polynomial])


def deformation_resultant(h, f, g, phi_f, phi_g):
    """Return Res(Hf+e phi_F,Hg+e phi_G) as an exact Z[e] polynomial."""
    base_f = multiply(h, f)
    base_g = multiply(h, g)
    degree_bound = (len(base_f) - 1) + (len(base_g) - 1)
    values = [
        resultant(add(base_f, scale(phi_f, epsilon)),
                  add(base_g, scale(phi_g, epsilon)))
        for epsilon in range(degree_bound + 1)
    ]
    return interpolate(values)


def tangent_data(h, f, g, phi_f, phi_g):
    r = len(h) - 1
    m = len(multiply(h, f)) - 1
    theta = remainder_monic(
        add(multiply(phi_f, g), scale(multiply(phi_g, f), -1)), h
    )
    norm_theta = resultant(h, theta)
    residual_resultant = resultant(f, g)
    predicted = (-1) ** (r * (m - r + 1)) * residual_resultant * norm_theta
    deformation = deformation_resultant(h, f, g, phi_f, phi_g)
    require(all(coefficient == 0 for coefficient in deformation[:r]),
            "resultant order below common-factor degree")
    actual = deformation[r] if r < len(deformation) else 0
    require(actual == predicted, "tangent-cone coefficient mismatch")
    return theta, norm_theta, residual_resultant, actual, deformation


def squarefree_data(h, theta):
    r = len(h) - 1
    h_prime = derivative(h)
    discriminant = (-1) ** (r * (r - 1) // 2) * resultant(h, h_prime)
    omega = remainder_monic(multiply(h_prime, theta), h)
    norm_omega = resultant(h, omega)
    require(
        norm_omega
        == (-1) ** (r * (r - 1) // 2) * discriminant * resultant(h, theta),
        "first-Witt/discriminant norm mismatch",
    )
    return discriminant, omega, norm_omega


def subtract(left, right):
    return add(left, scale(right, -1))


def quotient_multiply(left, right, h):
    return remainder_monic(multiply(left, right), h)


def main():
    cases = [
        ("linear", [-2, 1], [1, 1], [3, 2], [1, 1], [-1, 3], True),
        (
            "split2",
            [-2, 1, 1],
            [3, 1],
            [5, 2],
            [1, 0, 1],
            [-2, 0, 3],
            True,
        ),
        (
            "irreducible2",
            [1, 0, 1],
            [2, 1],
            [2, 1, 1],
            [3, -1, 1],
            [1, 2, -1, 1],
            True,
        ),
        (
            "split3",
            [8, -6, -3, 1],
            [5, 1],
            [7, 2],
            [1, 0, 1],
            [-1, 3],
            True,
        ),
        (
            "repeated3",
            [-1, 3, -3, 1],
            [2, 1],
            [3, 2],
            [1, 0, 1],
            [-1, 3],
            False,
        ),
    ]
    rows = []
    saved = {}
    for name, h, f, g, phi_f, phi_g, squarefree in cases:
        theta, norm_theta, residual_resultant, actual, deformation = tangent_data(
            h, f, g, phi_f, phi_g
        )
        first_order = next(
            index for index, coefficient in enumerate(deformation) if coefficient
        )
        require(first_order == len(h) - 1, "unexpected tangent order")

        # Changing the monic lift H~ to H~+pi*u subtracts u*f and u*g
        # from the two pi-columns and must not change theta.
        u = [2, -1, 1]
        shifted_phi_f = remainder_monic(subtract(phi_f, multiply(u, f)), h)
        shifted_phi_g = remainder_monic(subtract(phi_g, multiply(u, g)), h)
        shifted_theta = remainder_monic(
            subtract(multiply(shifted_phi_f, g), multiply(shifted_phi_g, f)), h
        )
        require(shifted_theta == theta, "conormal lift dependence")

        if squarefree:
            discriminant, omega, norm_omega = squarefree_data(h, theta)
            require(discriminant != 0 and norm_omega != 0,
                    "squarefree first-Witt unit control")
        else:
            discriminant = (-1) ** ((len(h) - 1) * (len(h) - 2) // 2) * resultant(
                h, derivative(h)
            )
            require(discriminant == 0, "repeated-root hostile is squarefree")
            omega = remainder_monic(multiply(derivative(h), theta), h)
            norm_omega = resultant(h, omega) if omega != [0] else 0
            require(norm_omega == 0 and norm_theta != 0,
                    "repeated-root first-Witt boundary missing")
        rows.append(
            (name, len(h) - 1, len(multiply(h, f)) - 1,
             len(multiply(h, g)) - 1, theta, norm_theta,
             residual_resultant, actual, discriminant, norm_omega)
        )
        saved[name] = (h, f, g, phi_f, phi_g, theta, norm_theta)

    # Direct branchwise first-Witt checks on the split degree-three case.
    h, f, g, phi_f, phi_g, theta, _ = saved["split3"]
    roots = [-2, 1, 4]
    branch_omegas = []
    h_prime = derivative(h)
    bar_f = multiply(h, f)
    bar_g = multiply(h, g)
    for root in roots:
        omega_direct = (
            evaluate(phi_f, root) * evaluate(derivative(bar_g), root)
            - evaluate(phi_g, root) * evaluate(derivative(bar_f), root)
        )
        omega_conormal = evaluate(h_prime, root) * evaluate(theta, root)
        require(omega_direct == omega_conormal != 0,
                "branch first-Witt mismatch")
        branch_omegas.append(omega_direct)

    # Polynomial-frame covariance in A=Q[v]/(v^2+1).
    h, f, g, phi_f, phi_g, theta, norm_theta = saved["irreducible2"]
    m11, m12, m21, m22 = [1, 1], [1], [0, 1], [2]
    phi_new_1 = remainder_monic(
        add(multiply(m11, phi_f), multiply(m12, phi_g)), h
    )
    phi_new_2 = remainder_monic(
        add(multiply(m21, phi_f), multiply(m22, phi_g)), h
    )
    f_new = remainder_monic(add(multiply(m11, f), multiply(m12, g)), h)
    g_new = remainder_monic(add(multiply(m21, f), multiply(m22, g)), h)
    theta_new = remainder_monic(
        subtract(multiply(phi_new_1, g_new), multiply(phi_new_2, f_new)), h
    )
    determinant = remainder_monic(
        subtract(multiply(m11, m22), multiply(m12, m21)), h
    )
    expected_theta_new = quotient_multiply(determinant, theta, h)
    require(theta_new == expected_theta_new, "polynomial-frame covariance")
    norm_determinant = resultant(h, determinant)
    require(
        resultant(h, theta_new) == norm_determinant * norm_theta,
        "frame norm scaling",
    )

    # A determinant that vanishes on one split component kills the norm.
    h, _, _, _, _, theta, norm_theta = saved["split2"]
    wall_determinant = [-1, 1]
    wall_theta = quotient_multiply(wall_determinant, theta, h)
    require(norm_theta != 0, "wall source theta is not a unit")
    require(resultant(h, wall_determinant) == 0, "determinant wall absent")
    require(resultant(h, wall_theta) == 0, "determinant wall did not kill norm")

    # If H is not the full gcd, Res(f,g)=0 and the first possible order rises.
    h = [-2, 1]
    shared = [1, 1]
    deformation = deformation_resultant(h, shared, shared, [1], [2])
    first_nonzero = next(
        index for index, coefficient in enumerate(deformation) if coefficient
    )
    require(first_nonzero > 1, "non-full-gcd hostile did not raise order")

    # THM-2598's cubic 1+2 discriminant gate is the r=1 specialization
    # for (S_A,S_A').  Use p=3,q=5,r=7 as an exact nondegenerate control.
    cubic_p, cubic_q, cubic_r = 3, 5, 7
    cubic_h = [cubic_p, 1]
    cubic_f = multiply([0, 1], cubic_h)
    cubic_g = [cubic_p, 3]
    cubic_phi = [-cubic_q**2, -4 * cubic_r]
    cubic_phi_derivative = [-4 * cubic_r]
    (
        cubic_theta,
        _,
        cubic_residual_resultant,
        cubic_resultant_tangent,
        _,
    ) = tangent_data(
        cubic_h,
        cubic_f,
        cubic_g,
        cubic_phi,
        cubic_phi_derivative,
    )
    cubic_omega = -2 * cubic_p * (4 * cubic_p * cubic_r - cubic_q**2)
    cubic_discriminant_tangent = -cubic_resultant_tangent
    require(cubic_theta == [cubic_omega], "cubic first-Witt coordinate")
    require(cubic_residual_resultant == -2 * cubic_p**2,
            "cubic residual cofactor")
    require(
        cubic_discriminant_tangent
        == -2 * cubic_p**2 * cubic_omega
        == 4 * cubic_p**3 * (4 * cubic_p * cubic_r - cubic_q**2),
        "cubic discriminant tangent",
    )

    print("RESULTANT TANGENT CONE / CONORMAL NORM EXACT CONTROL")
    for row in rows:
        print(
            "case=" + row[0]
            + f" rmn={row[1]},{row[2]},{row[3]}"
            + f" theta={row[4]} norm={row[5]} resfg={row[6]}"
            + f" tangent={row[7]} disc={row[8]} normomega={row[9]}"
        )
    print(f"split3_branch_omegas={branch_omegas}")
    print(
        "frame_covariance="
        f"det:{determinant},normdet:{norm_determinant},"
        f"theta_new:{theta_new},norm_new:{resultant(saved['irreducible2'][0], theta_new)}"
    )
    print(
        "determinant_wall="
        f"det:{wall_determinant},theta:{wall_theta},norm:0"
    )
    print(f"nonfull_gcd_first_nonzero_order={first_nonzero}")
    print(
        "cubic_1plus2="
        f"p:{cubic_p},q:{cubic_q},r:{cubic_r},omega:{cubic_omega},"
        f"res_over_A:{cubic_resultant_tangent},"
        f"disc_over_A:{cubic_discriminant_tangent}"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
