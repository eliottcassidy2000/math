#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2524.

The referee works with the integer circulants A_tau, exact rational
autocorrelations, the cyclotomic relation Phi_13, and F_547 only as an
independent check of the Fourier intertwining.  Every check remains active
under ``python -O``.
"""

from fractions import Fraction


P = 13
Q = 7
MOD = 547
ZETA = 475
CHI = (0, 1, 1, -1, 1, -1, -1)


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


def zero_matrix():
    return [[0 for _ in range(P)] for _ in range(P)]


def identity_matrix():
    return [[int(row == column) for column in range(P)] for row in range(P)]


def ones_matrix():
    return [[1 for _ in range(P)] for _ in range(P)]


def matrix_add(*terms):
    return [
        [sum(term[row][column] for term in terms) for column in range(P)]
        for row in range(P)
    ]


def matrix_scale(scale, matrix):
    return [[scale * entry for entry in row] for row in matrix]


def matrix_product(first, second):
    return [
        [
            sum(first[row][middle] * second[middle][column] for middle in range(P))
            for column in range(P)
        ]
        for row in range(P)
    ]


def matrix_vector(matrix, vector):
    return tuple(sum(entry * value for entry, value in zip(row, vector)) for row in matrix)


def matrix_power(matrix, exponent):
    result = identity_matrix()
    base = matrix
    while exponent:
        if exponent & 1:
            result = matrix_product(result, base)
        base = matrix_product(base, base)
        exponent //= 2
    return result


def a_matrix(tau):
    """Matrix of (A_tau x)_t=-sum_s chi(s)(x_(t+2tau s)+x_(t-2tau s))."""

    matrix = zero_matrix()
    for target in range(P):
        for s in range(1, Q):
            difference = 2 * tau * s
            matrix[target][(target + difference) % P] -= CHI[s]
            matrix[target][(target - difference) % P] -= CHI[s]
    return matrix


def centered_basis():
    return tuple(
        tuple(-1 if position == 0 else int(position == index) for position in range(P))
        for index in range(1, P)
    )


def normalized_transform_mod(values, frequency):
    inverse_p = pow(P, MOD - 2, MOD)
    return (
        inverse_p
        * sum(
            (value % MOD) * pow(ZETA, (-frequency * position) % P, MOD)
            for position, value in enumerate(values)
        )
    ) % MOD


def lambda_mod(tau, frequency):
    return (
        -sum(
            CHI[s]
            * (
                pow(ZETA, (2 * frequency * tau * s) % P, MOD)
                + pow(ZETA, (-2 * frequency * tau * s) % P, MOD)
            )
            for s in range(1, Q)
        )
    ) % MOD


def exact_primitive_mode_nonzero(values, frequency):
    """Test sum_t values[t] zeta^(-frequency*t) != 0 in Q(zeta_13)."""

    coefficients = [0] * P
    for position, value in enumerate(values):
        coefficients[(-frequency * position) % P] += value
    # Modulo Phi_13, X^12=-(1+...+X^11).  The remainder is zero exactly
    # when all thirteen coefficients are equal.
    return any(coefficients[index] != coefficients[P - 1] for index in range(P - 1))


def autocorrelation(profile):
    return tuple(
        sum(profile[root] * profile[(root + shift) % P] for root in range(P))
        for shift in range(P)
    )


def cross_correlation(first, second):
    return tuple(
        sum(first[root] * second[(root + shift) % P] for root in range(P))
        for shift in range(P)
    )


def translated_bank_from_correlation(correlation, tau):
    # If B_t=correlation_t/13, this is
    # -13 sum_s chi(s)(B_(t+2tau s)+B_(t-2tau s)).
    return tuple(
        -sum(
            CHI[s]
            * (
                correlation[(translation + 2 * tau * s) % P]
                + correlation[(translation - 2 * tau * s) % P]
            )
            for s in range(1, Q)
        )
        for translation in range(P)
    )


def translated_gradient_polarization(first, second, tau):
    return tuple(
        sum(
            CHI[s]
            * sum(
                (first[(root - tau * s) % P] - first[(root + tau * s) % P])
                * (
                    second[(root + translation - tau * s) % P]
                    - second[(root + translation + tau * s) % P]
                )
                for root in range(P)
            )
            for s in range(1, Q)
        )
        for translation in range(P)
    )


def audit_matrices_and_inverse():
    identity = identity_matrix()
    ones = ones_matrix()
    basis = centered_basis()
    inverse_entries = 0
    fourier_entries = 0
    nonzero_lambdas = 0

    for tau in range(1, P):
        matrix = a_matrix(tau)
        require(matrix == [list(row) for row in zip(*matrix)], ("symmetric", tau))
        require(all(sum(row) == 0 for row in matrix), ("constant kernel", tau))
        require(all(matrix[index][index] == 0 for index in range(P)), ("zero diagonal", tau))
        require(
            {matrix[row][column] for row in range(P) for column in range(P) if row != column}
            == {-1, 1},
            ("signed complete graph", tau),
        )

        square = matrix_power(matrix, 2)
        cube = matrix_power(matrix, 3)
        fourth = matrix_power(matrix, 4)
        fifth = matrix_power(matrix, 5)
        sixth = matrix_power(matrix, 6)
        polynomial = matrix_add(
            sixth,
            matrix_scale(-39, fourth),
            matrix_scale(299, square),
            matrix_scale(-325, identity),
            matrix_scale(25, ones),
        )
        require(polynomial == zero_matrix(), ("degree-six polynomial", tau))

        inverse_numerator = matrix_add(fifth, matrix_scale(-39, cube), matrix_scale(299, matrix))
        expected_inverse_product = matrix_add(matrix_scale(325, identity), matrix_scale(-25, ones))
        require(
            matrix_product(matrix, inverse_numerator) == expected_inverse_product,
            ("left inverse", tau),
        )
        require(
            matrix_product(inverse_numerator, matrix) == expected_inverse_product,
            ("right inverse", tau),
        )

        for vector in basis:
            response = matrix_vector(matrix, vector)
            recovered = matrix_vector(inverse_numerator, response)
            require(recovered == tuple(325 * value for value in vector), ("basis inverse", tau))
            inverse_entries += P

            for frequency in range(1, P):
                multiplier = lambda_mod(tau, frequency)
                require(multiplier != 0, ("lambda mod prime", tau, frequency))
                actual = normalized_transform_mod(response, frequency)
                expected = multiplier * normalized_transform_mod(vector, frequency) % MOD
                require(actual == expected, ("Fourier factorization", tau, frequency))
                fourier_entries += 1

        nonzero_lambdas += P - 1

    return inverse_entries, fourier_entries, nonzero_lambdas


def audit_binary_self_correlations():
    profiles = 0
    slope_profiles = 0
    primitive_modes = 0
    zero_diagonal = 0

    for mask in range(1, (1 << P) - 1):
        profile = tuple((mask >> root) & 1 for root in range(P))
        correlation = autocorrelation(profile)
        total = sum(profile)
        drift = Fraction(P * sum(value * value for value in profile) - total * total, P * P)
        require(drift > 0, ("positive binary drift", mask))

        for tau in range(1, Q):
            matrix = a_matrix(tau)
            response = translated_bank_from_correlation(correlation, tau)
            require(response == matrix_vector(matrix, correlation), ("translated formula", mask, tau))
            require(sum(response) == 0 and any(response), ("centred nonzero bank", mask, tau))
            if response[0] == 0:
                zero_diagonal += 1
            for frequency in range(1, P):
                require(
                    exact_primitive_mode_nonzero(response, frequency),
                    ("rational primitive mode", mask, tau, frequency),
                )
                primitive_modes += 1
            slope_profiles += 1
        profiles += 1

    return profiles, slope_profiles, primitive_modes, zero_diagonal


def audit_delta_and_cross_controls():
    delta = (1,) + (0,) * (P - 1)
    correlation = autocorrelation(delta)
    collision_profile = tuple(Fraction(value, P) for value in correlation)
    collision_mean = sum(collision_profile) / P
    centered_collision = tuple(value - collision_mean for value in collision_profile)
    require(centered_collision[0] == Fraction(12, 169), "delta drift")
    delta_responses = []
    for tau in range(1, Q):
        response = translated_bank_from_correlation(correlation, tau)
        require(
            response
            == tuple(13 * value for value in matrix_vector(a_matrix(tau), centered_collision)),
            ("delta R=13Ab", tau),
        )
        require(response[0] == 0, ("delta isotropic diagonal", tau))
        require(sorted(response[1:]) == [-1] * 6 + [1] * 6, ("delta translated bank", tau))
        delta_responses.append(response)

    first = (P - 1,) + (-1,) * (P - 1)
    second = tuple(first[(root + 1) % P] for root in range(P))
    cross = cross_correlation(first, second)
    reverse_cross = cross_correlation(second, first)
    response = translated_bank_from_correlation(cross, 1)
    reverse_response = translated_bank_from_correlation(reverse_cross, 1)
    gradient_response = translated_gradient_polarization(first, second, 1)
    require(response == gradient_response, "cross gradient polarization")
    require(
        reverse_response == tuple(response[-translation % P] for translation in range(P)),
        "cross reversal",
    )
    require(any(response[t] != response[-t % P] for t in range(1, Q)), "oriented odd component")

    self_response = translated_bank_from_correlation(autocorrelation(first), 1)
    require(
        self_response == tuple(self_response[-translation % P] for translation in range(P)),
        "self correlation remains even",
    )
    return len(delta_responses), len(response), len(reverse_response)


def main():
    require(pow(ZETA, P, MOD) == 1 and ZETA != 1, "primitive thirteenth root")
    inverse_entries, fourier_entries, lambdas = audit_matrices_and_inverse()
    profiles, slope_profiles, primitive_modes, zero_diagonal = audit_binary_self_correlations()
    delta_slopes, cross_entries, reverse_entries = audit_delta_and_cross_controls()

    print("THM-2524 exact translated chi7 Hamilton-polarization referee")
    print(
        "matrix_audit:",
        "slopes=12",
        "centered_rank=12",
        "inverse_polynomial=A5-39A3+299A",
        f"inverse_entries={inverse_entries}",
    )
    print(
        "fourier_audit:",
        f"nonzero_lambdas={lambdas}",
        f"factor_entries={fourier_entries}",
        "inverse_denominator=325",
    )
    print(
        "binary_self_audit:",
        f"profiles={profiles}",
        f"slope_profiles={slope_profiles}",
        f"primitive_modes={primitive_modes}",
        f"zero_diagonal_controls={zero_diagonal}",
    )
    print(
        "delta_hostile:",
        "D13=12/169",
        "R_tau(0)=0",
        f"translated_slopes={delta_slopes}",
        "off_diagonal=6*(-1)+6*(+1)",
    )
    print(
        "cross_orientation:",
        f"direct_entries={cross_entries}",
        f"reverse_entries={reverse_entries}",
        "odd_component=preserved",
        "self_component=even",
    )
    print(
        "VERIFIED: the full translated chi7 bank recovers the centred collision profile; "
        "a diagonal contrast may vanish, and no oriented Boolean owner current is created."
    )


if __name__ == "__main__":
    main()
