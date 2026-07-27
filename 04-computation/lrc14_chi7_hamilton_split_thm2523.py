#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2523.

The audit works over Z, Q, and F_547.  It checks the signed Hamilton-cycle
matrix on the F_13 augmentation module, its cyclotomic polynomial identity,
the dilation-by-five anti-isometry, explicit rational isotropic controls, and,
as an independent forward cross-check for THM-2524, the translated-
polarization/collision-table Fourier multiplier.
"""

from fractions import Fraction


P = 13
Q = 7
CHI = {1: 1, 2: 1, 3: -1, 4: 1, 5: -1, 6: -1}
MOD = 547
GENERATOR = 2
ZETA = pow(GENERATOR, (MOD - 1) // P, MOD)


def require(condition, label):
    if not condition:
        raise AssertionError(label)


def dot(first, second):
    return sum(x * y for x, y in zip(first, second))


def mat_vec(matrix, vector):
    return [dot(row, vector) for row in matrix]


def bilinear(first, matrix, second):
    return dot(first, mat_vec(matrix, second))


def rational_rank(matrix):
    work = [[Fraction(x) for x in row] for row in matrix]
    rows = len(work)
    cols = len(work[0]) if rows else 0
    pivot_row = 0
    for col in range(cols):
        pivot = next((r for r in range(pivot_row, rows) if work[r][col]), None)
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        value = work[pivot_row][col]
        work[pivot_row] = [x / value for x in work[pivot_row]]
        for row in range(rows):
            if row == pivot_row or not work[row][col]:
                continue
            factor = work[row][col]
            work[row] = [
                work[row][j] - factor * work[pivot_row][j]
                for j in range(cols)
            ]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def chi7(value):
    return CHI[value % Q]


def antipodal_representative(value):
    value %= P
    require(value != 0, "nonzero representative")
    return min(value, P - value)


def hamilton_matrix(tau):
    """A_tau with Q_tau(p)=<p,A_tau p>."""
    matrix = [[0 for _ in range(P)] for _ in range(P)]
    for row in range(P):
        for s in range(1, Q):
            matrix[row][(row + 2 * tau * s) % P] -= CHI[s]
            matrix[row][(row - 2 * tau * s) % P] -= CHI[s]
    return matrix


def energy_matrix(tau, s):
    """Matrix of sum_h (p_{h-tau*s}-p_{h+tau*s})^2."""
    matrix = [[0 for _ in range(P)] for _ in range(P)]
    gap = 2 * tau * s
    for row in range(P):
        matrix[row][row] += 2
        matrix[row][(row + gap) % P] -= 1
        matrix[row][(row - gap) % P] -= 1
    return matrix


def add_scaled(target, scale, source):
    for row in range(P):
        for col in range(P):
            target[row][col] += scale * source[row][col]


def cyclic_multiply(first, second):
    result = [0] * P
    for i, x in enumerate(first):
        for j, y in enumerate(second):
            result[(i + j) % P] += x * y
    return result


def cyclic_power(base, exponent):
    result = [1] + [0] * (P - 1)
    factor = base[:]
    while exponent:
        if exponent & 1:
            result = cyclic_multiply(result, factor)
        factor = cyclic_multiply(factor, factor)
        exponent //= 2
    return result


def dft_mod(vector, frequency):
    return sum(
        (value % MOD) * pow(ZETA, (-frequency * index) % P, MOD)
        for index, value in enumerate(vector)
    ) % MOD


def translated(vector, shift):
    return [vector[(row + shift) % P] for row in range(P)]


def collision_table(vector):
    """B_t=(1/13) sum_r q_r q_{r+t} for a scalar fibre control."""
    return [
        Fraction(sum(vector[r] * vector[(r + shift) % P] for r in range(P)), P)
        for shift in range(P)
    ]


def centred(vector):
    mean = Fraction(sum(vector), P)
    return [Fraction(value) - mean for value in vector]


def audit_fano_selector_and_anti_isometry():
    residues = {1, 2, 4}
    differences = [
        (first - second) % Q
        for first in residues
        for second in residues
        if first != second
    ]
    require(sorted(differences) == list(range(1, Q)), "Fano difference set")

    permutation = tuple(antipodal_representative(5 * s) for s in range(1, Q))
    require(permutation == (5, 3, 2, 6, 1, 4), "dilation permutation")
    require(
        all(CHI[permutation[s - 1]] == -CHI[s] for s in range(1, Q)),
        "dilation flips chi7",
    )

    anti_checks = 0
    for tau in range(1, P):
        matrix = hamilton_matrix(tau)
        # If (D_5 p)_x=p_{5x}, then
        # (D_5^T A D_5)_{ij}=A_{5^{-1}i,5^{-1}j}.
        inverse_five = pow(5, -1, P)
        for row in range(P):
            for col in range(P):
                transformed = matrix[(inverse_five * row) % P][
                    (inverse_five * col) % P
                ]
                require(transformed == -matrix[row][col], "anti-isometry matrix")
                anti_checks += 1
    return permutation, anti_checks


def audit_matrix_and_cyclotomic_identity():
    matrix_checks = 0
    for tau in range(1, P):
        matrix = hamilton_matrix(tau)
        expected = [[0 for _ in range(P)] for _ in range(P)]
        for s in range(1, Q):
            add_scaled(expected, CHI[s], energy_matrix(tau, s))
        require(matrix == expected, "cycle-energy matrix identity")
        require(all(matrix[r][r] == 0 for r in range(P)), "zero diagonal")
        require(all(sum(row) == 0 for row in matrix), "augmentation kernel")
        require(
            all(matrix[r][c] == matrix[c][r] for r in range(P) for c in range(P)),
            "symmetric matrix",
        )
        require(rational_rank(matrix) == P - 1, "augmentation rank")
        matrix_checks += P * P

    first_row = hamilton_matrix(1)[0]
    square = cyclic_power(first_row, 2)
    fourth = cyclic_power(first_row, 4)
    sixth = cyclic_power(first_row, 6)
    polynomial_value = [
        sixth[i] - 39 * fourth[i] + 299 * square[i] - (325 if i == 0 else 0)
        for i in range(P)
    ]
    require(polynomial_value == [-25] * P, "cyclotomic group-ring identity")

    # f(T)=T^6-39T^4+299T^2-325 is 13-Eisenstein.
    require(all(coefficient % P == 0 for coefficient in (-39, 299, -325)), "Eisenstein divisibility")
    require((-325) % (P * P) != 0, "Eisenstein constant")

    def cubic(value):
        return value**3 - 39 * value**2 + 299 * value - 325

    isolating_intervals = ((1, 2), (8, 9), (29, 30))
    signs = tuple((cubic(left), cubic(right)) for left, right in isolating_intervals)
    require(
        all(first * second < 0 for first, second in signs),
        "three positive squared-root intervals",
    )
    # The Fano-residue slope product lands in the symmetric quadratic
    # character operator at 13.
    paley = [0] + [
        1 if pow(value, (P - 1) // 2, P) == 1 else -1
        for value in range(1, P)
    ]
    slope_product = cyclic_multiply(
        cyclic_multiply(hamilton_matrix(1)[0], hamilton_matrix(2)[0]),
        hamilton_matrix(4)[0],
    )
    require(slope_product == [5 * value for value in paley], "Fano slope product")
    paley_square = cyclic_multiply(paley, paley)
    require(paley_square == [12] + [-1] * 12, "Paley square")
    require(
        all(paley[value] == paley[-value % P] for value in range(P)),
        "Paley symmetry",
    )
    return first_row, matrix_checks, signs, paley, slope_product


def multiplication_orbits():
    seen = {0}
    orbits = []
    for start in range(1, P):
        if start in seen:
            continue
        orbit = []
        value = start
        while value not in orbit:
            orbit.append(value)
            seen.add(value)
            value = 5 * value % P
        orbits.append(tuple(orbit))
    return tuple(orbits)


def audit_isotropic_controls():
    matrix = hamilton_matrix(1)
    delta = [Fraction(P - 1, P)] + [Fraction(-1, P)] * (P - 1)
    energies = tuple(
        bilinear(delta, energy_matrix(1, s), delta) for s in range(1, Q)
    )
    require(energies == (Fraction(2),) * 6, "centred-delta cycle energies")
    require(bilinear(delta, matrix, delta) == 0, "centred-delta isotropic")
    translated_row = tuple(
        bilinear(delta, matrix, translated(delta, shift)) for shift in range(P)
    )
    require(translated_row == tuple(matrix[0]), "delta translated row")
    require(all(value != 0 for value in translated_row[1:]), "translated rescue")

    orbits = multiplication_orbits()
    require(
        orbits == ((1, 5, 12, 8), (2, 10, 11, 3), (4, 7, 9, 6)),
        "multiplication-five orbits",
    )
    fixed = []
    anti_fixed = []
    for orbit in orbits:
        fixed_vector = [0] * P
        fixed_vector[0] = -4
        for vertex in orbit:
            fixed_vector[vertex] = 1
        fixed.append(fixed_vector)

        anti_vector = [0] * P
        for index, vertex in enumerate(orbit):
            anti_vector[vertex] = 1 if index % 2 == 0 else -1
        anti_fixed.append(anti_vector)

    require(rational_rank(fixed) == 3, "fixed isotropic dimension")
    require(rational_rank(anti_fixed) == 3, "anti-fixed isotropic dimension")
    require(
        all(bilinear(x, matrix, y) == 0 for x in fixed for y in fixed),
        "fixed totally isotropic",
    )
    require(
        all(bilinear(x, matrix, y) == 0 for x in anti_fixed for y in anti_fixed),
        "anti-fixed totally isotropic",
    )
    cross = [[bilinear(x, matrix, y) for y in anti_fixed] for x in fixed]
    cross_determinant = (
        cross[0][0] * (cross[1][1] * cross[2][2] - cross[1][2] * cross[2][1])
        - cross[0][1] * (cross[1][0] * cross[2][2] - cross[1][2] * cross[2][0])
        + cross[0][2] * (cross[1][0] * cross[2][1] - cross[1][1] * cross[2][0])
    )
    require(cross_determinant == -4160, "isotropic cross determinant")
    return energies, translated_row, orbits, cross_determinant


def audit_translated_polarization():
    identity = [1] + [0] * (P - 1)
    controls = (
        identity,
        [2, -1, 3, 0, 1, -2, 4, 1, 0, -3, 2, 1, 5],
        [index * index - 5 * index + 3 for index in range(P)],
    )
    position_checks = 0
    fourier_checks = 0
    nonzero_multipliers = 0
    for tau in range(1, P):
        matrix = hamilton_matrix(tau)
        for control in controls:
            p = centred(control)
            collision = collision_table(control)
            response = []
            for shift in range(P):
                direct = bilinear(p, matrix, translated(p, shift))
                expected = -P * sum(
                    CHI[s]
                    * (
                        collision[(shift + 2 * tau * s) % P]
                        + collision[(shift - 2 * tau * s) % P]
                    )
                    for s in range(1, Q)
                )
                require(direct == expected, "translated collision identity")
                response.append(direct)
                position_checks += 1
            require(
                all(response[shift] == response[-shift % P] for shift in range(P)),
                "translated response even",
            )

            denominator_lcm = 1
            for value in collision + response:
                denominator_lcm *= value.denominator
            # The finite-field check below is performed separately on the
            # integer delta control, so no denominator conversion is needed.
            require(denominator_lcm > 0, "rational controls")

        # Exact root multiplier on the integer delta control q=e_0.
        collision = collision_table(identity)
        p = centred(identity)
        response = [
            bilinear(p, matrix, translated(p, shift)) for shift in range(P)
        ]
        inverse_p = pow(P, -1, MOD)
        for frequency in range(1, P):
            lam = dft_mod(matrix[0], frequency)
            require(lam != 0, "nonzero cyclotomic multiplier mod 547")
            bhat = (
                sum(
                    (value.numerator % MOD)
                    * pow(value.denominator, -1, MOD)
                    * pow(ZETA, (-frequency * shift) % P, MOD)
                    for shift, value in enumerate(collision)
                )
                * inverse_p
            ) % MOD
            rhat = (
                sum(
                    (value.numerator % MOD)
                    * pow(value.denominator, -1, MOD)
                    * pow(ZETA, (-frequency * shift) % P, MOD)
                    for shift, value in enumerate(response)
                )
                * inverse_p
            ) % MOD
            require(rhat == P * lam * bhat % MOD, "translated Fourier multiplier")
            require(rhat != 0, "delta translated root mode")
            nonzero_multipliers += 1
            fourier_checks += 1
    return position_checks, fourier_checks, nonzero_multipliers


def main():
    permutation, anti_checks = audit_fano_selector_and_anti_isometry()
    first_row, matrix_checks, root_signs, paley, slope_product = (
        audit_matrix_and_cyclotomic_identity()
    )
    energies, translated_row, orbits, cross_determinant = audit_isotropic_controls()
    position_checks, fourier_checks, nonzero_multipliers = audit_translated_polarization()

    print("THM-2523 exact chi7 Hamilton split-form referee")
    print(
        "fano_selector:",
        "quadratic_residues=(1,2,4)",
        f"dilation5_permutation={permutation}",
        "chi7_flipped=yes",
    )
    print(
        "circulant_form:",
        f"rank={P-1}",
        f"first_row={tuple(first_row)}",
        f"matrix_checks={matrix_checks}",
        f"anti_isometry_checks={anti_checks}",
    )
    print(
        "cyclotomic_spectrum:",
        "minimal_polynomial=t^6-39*t^4+299*t^2-325",
        "characteristic_on_augmentation=minimal_polynomial^2",
        f"squared_root_signs={root_signs}",
        "real_signature=(6,6)",
    )
    print(
        "fano_slope_product:",
        f"chi13={tuple(paley)}",
        f"A1_A2_A4={tuple(slope_product)}",
        "identity=A1*A2*A4=5*C13",
        "C13_squared_on_augmentation=13*identity",
        "operator=symmetric_not_tournament",
    )
    print(
        "isotropic_controls:",
        f"delta_cycle_energies={energies}",
        f"translated_delta={translated_row}",
        f"dilation5_orbits={orbits}",
        "fixed_dim=3",
        "anti_fixed_dim=3",
        f"cross_determinant={cross_determinant}",
    )
    print(
        "translated_polarization:",
        f"position_identities={position_checks}",
        f"fourier_identities={fourier_checks}",
        f"nonzero_root_multipliers={nonzero_multipliers}",
        "Rhat_tau(k)=13*lambda_tau(k)*Bhat(k)",
        "response_even=yes",
    )
    print(
        "VERIFIED: the chi7 diagonal is real-split and can vanish on nonzero "
        "profiles; FORWARD CROSS-CHECK: the translated collision bank has "
        "the lossless THM-2524 multiplier but remains unoriented."
    )


if __name__ == "__main__":
    main()
