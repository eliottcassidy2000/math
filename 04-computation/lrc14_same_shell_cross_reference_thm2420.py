#!/usr/bin/env python3
"""Dependency-free exact referee for candidate THM-2420."""

from itertools import product


P = 13
ZERO = (0, 0)
ONE = (1, 0)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def gadd(left, right):
    return left[0] + right[0], left[1] + right[1]


def gneg(value):
    return -value[0], -value[1]


def gsub(left, right):
    return gadd(left, gneg(right))


def gmul(left, right):
    return (
        left[0] * right[0] - left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def gconj(value):
    return value[0], -value[1]


def gscale(scale, value):
    return scale * value[0], scale * value[1]


def gsum(values):
    total = ZERO
    for value in values:
        total = gadd(total, value)
    return total


def gnorm2(value):
    return value[0] ** 2 + value[1] ** 2


def dot(left, right):
    return sum(a * b for a, b in zip(left, right))


def add_mod(left, right, modulus):
    return tuple((a + b) % modulus for a, b in zip(left, right))


def kernel(row, modulus):
    return tuple(
        point
        for point in product(range(modulus), repeat=len(row))
        if dot(row, point) % modulus == 0
    )


def dictionary_sum(function):
    return gsum(function.values())


def correlate(group, charged, reference, modulus):
    result = {}
    for target in group:
        total = ZERO
        for base, value in reference.items():
            shifted = add_mod(target, base, modulus)
            if shifted in charged:
                total = gadd(
                    total,
                    gmul(charged[shifted], gconj(value)),
                )
        if total != ZERO:
            result[target] = total
    return result


def transform_polynomial(function, character, modulus):
    """Unnormalized sum f(x) T^(-character.x) in Z[i][T]/(T^M-1)."""

    polynomial = [ZERO] * modulus
    for point, value in function.items():
        exponent = (-dot(character, point)) % modulus
        polynomial[exponent] = gadd(polynomial[exponent], value)
    return tuple(polynomial)


def conjugated_transform_polynomial(function, character, modulus):
    """Conjugate of the preceding transform, exactly."""

    polynomial = [ZERO] * modulus
    for point, value in function.items():
        exponent = dot(character, point) % modulus
        polynomial[exponent] = gadd(
            polynomial[exponent],
            gconj(value),
        )
    return tuple(polynomial)


def polynomial_multiply(left, right):
    modulus = len(left)
    output = [ZERO] * modulus
    for first, left_value in enumerate(left):
        if left_value == ZERO:
            continue
        for second, right_value in enumerate(right):
            if right_value == ZERO:
                continue
            exponent = (first + second) % modulus
            output[exponent] = gadd(
                output[exponent],
                gmul(left_value, right_value),
            )
    return tuple(output)


def integer_polynomial_multiply(left, right, modulus):
    output = [0] * modulus
    for first, left_value in enumerate(left):
        for second, right_value in enumerate(right):
            output[(first + second) % modulus] += left_value * right_value
    return tuple(output)


def cyclotomic_reduce_13(polynomial):
    require(len(polynomial) == P, "wrong Phi_13 polynomial length")
    top = polynomial[-1]
    return tuple(value - top for value in polynomial[:-1])


def integer_monomial(exponent, coefficient=1):
    polynomial = [0] * P
    polynomial[exponent % P] = coefficient
    return tuple(polynomial)


def deterministic_value(point, salt):
    first = 1 + (sum((index + salt) * value for index, value in enumerate(point)) % 4)
    second = (
        sum((2 * index + salt + 1) * value for index, value in enumerate(point))
        % 5
    ) - 2
    return first, second


def finite_shell_checks():
    cases = (
        ((1, 0, 0), 13),
        ((1, 2, 3), 26),
        ((2, 3, 5), 39),
        ((3, 4, 6), 26),
    )
    fibre_checks = 0
    correlation_checks = 0
    fourier_checks = 0
    sharp_checks = 0

    for row, modulus in cases:
        group = kernel(row, modulus)
        prime_group = kernel(tuple(value % P for value in row), P)
        fibres = {}
        for point in group:
            residue = tuple(value % P for value in point)
            fibres.setdefault(residue, []).append(point)
        require(set(fibres) == set(prime_group), "kernel reduction not surjective")
        expected_size = (modulus // P) ** (len(row) - 1)
        require(
            {len(fibre) for fibre in fibres.values()} == {expected_size},
            "kernel reduction fibre size failed",
        )
        fibre_checks += len(fibres)

        nonzero_residues = [
            residue
            for residue in prime_group
            if any(residue)
        ]
        charged_residue = nonzero_residues[0]
        charged = {
            point: deterministic_value(point, 2)
            for point in fibres[charged_residue]
        }
        reference = {
            point: deterministic_value(point, 3)
            for point in fibres[(0, 0, 0)]
        }
        charged_sum = dictionary_sum(charged)
        reference_sum = dictionary_sum(reference)
        require(charged_sum != ZERO, "charged total vanished in control")
        require(reference_sum != ZERO, "reference total vanished in control")

        correlation = correlate(group, charged, reference, modulus)
        require(
            all(
                tuple(value % P for value in point) == charged_residue
                for point in correlation
            ),
            "correlation escaped charged residue fibre",
        )
        product_total = gmul(charged_sum, gconj(reference_sum))
        require(
            dictionary_sum(correlation) == product_total,
            "correlation total identity failed",
        )
        maximum = max(gnorm2(value) for value in correlation.values())
        energy = sum(gnorm2(value) for value in correlation.values())
        require(
            maximum * expected_size**2 >= gnorm2(product_total),
            "sharp maximum lower bound failed",
        )
        require(
            energy * expected_size >= gnorm2(product_total),
            "sharp energy lower bound failed",
        )
        correlation_checks += 1

        characters = tuple(product(range(min(modulus, 4)), repeat=len(row)))
        for character in characters:
            direct = transform_polynomial(correlation, character, modulus)
            factored = polynomial_multiply(
                transform_polynomial(charged, character, modulus),
                conjugated_transform_polynomial(
                    reference,
                    character,
                    modulus,
                ),
            )
            require(direct == factored, "Fourier product identity failed")
            fourier_checks += 1

        zero_point = fibres[(0, 0, 0)][0]
        sharp_charged = {
            point: ONE
            for point in fibres[charged_residue]
        }
        sharp_reference = {zero_point: ONE}
        sharp_correlation = correlate(
            group,
            sharp_charged,
            sharp_reference,
            modulus,
        )
        require(
            set(sharp_correlation.values()) == {ONE}
            and len(sharp_correlation) == expected_size,
            "sharp uniform correlation control failed",
        )
        sharp_total = dictionary_sum(sharp_correlation)
        sharp_maximum = max(
            gnorm2(value)
            for value in sharp_correlation.values()
        )
        sharp_energy = sum(
            gnorm2(value)
            for value in sharp_correlation.values()
        )
        require(
            sharp_maximum * expected_size**2 == gnorm2(sharp_total),
            "maximum constant was not sharp",
        )
        require(
            sharp_energy * expected_size == gnorm2(sharp_total),
            "energy constant was not sharp",
        )
        sharp_checks += 1

    return fibre_checks, correlation_checks, fourier_checks, sharp_checks


def address_pushforward_checks():
    parameters = tuple(range(-4, 5))
    charged = {
        parameter: (parameter * parameter + 2, parameter - 1)
        for parameter in parameters
    }
    reference = {
        parameter: (2 * parameter + 3, 1 - parameter * parameter)
        for parameter in parameters
    }
    exact = {}
    for difference in range(-8, 9):
        value = ZERO
        for base in parameters:
            shifted = base + difference
            if shifted in charged:
                value = gadd(
                    value,
                    gmul(charged[shifted], gconj(reference[base])),
                )
        if value != ZERO:
            exact[difference] = value

    l1_exact = sum(
        abs(value[0]) + abs(value[1])
        for value in exact.values()
    )
    l1_charged = sum(
        abs(value[0]) + abs(value[1])
        for value in charged.values()
    )
    l1_reference = sum(
        abs(value[0]) + abs(value[1])
        for value in reference.values()
    )
    require(
        l1_exact <= l1_charged * l1_reference,
        "finite Gaussian coefficient-l1 control failed",
    )

    charged_push = {}
    reference_push = {}
    for parameter, value in charged.items():
        residue = parameter % P
        charged_push[residue] = gadd(
            charged_push.get(residue, ZERO),
            value,
        )
    for parameter, value in reference.items():
        residue = parameter % P
        reference_push[residue] = gadd(
            reference_push.get(residue, ZERO),
            value,
        )

    exact_push = {}
    for difference, value in exact.items():
        residue = difference % P
        exact_push[residue] = gadd(exact_push.get(residue, ZERO), value)

    finite_push = {}
    for target in range(P):
        value = ZERO
        for base, reference_value in reference_push.items():
            shifted = (target + base) % P
            if shifted in charged_push:
                value = gadd(
                    value,
                    gmul(
                        charged_push[shifted],
                        gconj(reference_value),
                    ),
                )
        if value != ZERO:
            finite_push[target] = value
    require(exact_push == finite_push, "exact-address pushforward failed")

    exact_total = dictionary_sum(exact)
    require(
        exact_total
        == gmul(dictionary_sum(charged), gconj(dictionary_sum(reference))),
        "exact-address total failed",
    )
    return len(exact), len(exact_push)


def hostile_checks():
    speeds = (1, 2, 26)
    safe_character = (1, 6, 2)
    deepest_character = 12
    address = (
        safe_character[0],
        safe_character[1],
        (safe_character[2] + deepest_character) % P,
    )
    require(dot(speeds, safe_character) % P == 0, "safe character not closed")
    require(address == (1, 6, 1), "all-unit address mismatch")
    require(all(address), "address lost a coordinate unit")
    require(
        deepest_character != (-safe_character[2]) % P,
        "deepest character hit the ineligible diagonal",
    )
    require(dot(speeds, address) % P == 0, "affine address not periodic")

    phase_word = []
    sigma_bank = []
    for cell in range(2 * P):
        sigma_1 = 0
        sigma_2 = cell // P
        sigma_26 = cell % P
        phase = (
            safe_character[0] * sigma_1
            + safe_character[1] * sigma_2
            + (safe_character[2] + deepest_character) * sigma_26
        ) % P
        phase_word.append(phase)
        sigma_bank.append((sigma_1, sigma_2, sigma_26))
    phase_word = tuple(phase_word)
    require(
        phase_word
        == tuple(range(P))
        + tuple((residue + 6) % P for residue in range(P)),
        "hostile phase word derivation failed",
    )
    phase_counts = tuple(phase_word.count(residue) for residue in range(P))
    require(phase_counts == (2,) * P, "hostile mean phase census failed")

    multiplier_counts = [0] * P
    for cell, phase in enumerate(phase_word):
        multiplier_counts[(phase - cell) % P] += 1
    expected_multiplier = [0] * P
    expected_multiplier[0] = P
    expected_multiplier[6] = P
    require(
        tuple(multiplier_counts) == tuple(expected_multiplier),
        "sideband multiplier census failed",
    )

    first_factor = [0] * P
    first_factor[0] = 1
    first_factor[P - 1] = -1
    second_factor = [0] * P
    second_factor[0] = 1
    second_factor[6] = 1
    numerator = integer_polynomial_multiply(
        first_factor,
        second_factor,
        P,
    )
    require(
        any(cyclotomic_reduce_13(numerator)),
        "sideband cyclotomic numerator vanished",
    )

    zero_reference_constants = []
    for residue in range(P):
        cell_values = []
        for _, _, sigma_26 in sigma_bank:
            deepest = integer_monomial(residue * sigma_26)
            safe_1 = integer_monomial(0, 12)
            safe_2 = integer_monomial(0, 12)
            safe_26_character = (-residue) % P
            if safe_26_character == 0:
                safe_26 = integer_monomial(0, 12)
            else:
                safe_26 = integer_monomial(
                    safe_26_character * sigma_26,
                    -1,
                )
            value = integer_polynomial_multiply(
                integer_polynomial_multiply(deepest, safe_1, P),
                integer_polynomial_multiply(safe_2, safe_26, P),
                P,
            )
            cell_values.append(cyclotomic_reduce_13(value))
        require(
            len(set(cell_values)) == 1,
            "zero-reference sector was not constant across cells",
        )
        reduced = cell_values[0]
        require(
            not any(reduced[1:]),
            "zero-reference constant retained a cyclotomic phase",
        )
        zero_reference_constants.append(reduced[0])
    require(
        zero_reference_constants[0] == 12**3
        and set(zero_reference_constants[1:]) == {-(12**2)},
        "complete zero-reference bank constants failed",
    )
    require(
        all(value != 0 for value in zero_reference_constants),
        "zero-reference sector unexpectedly vanished pointwise",
    )
    self_gram_cells = []
    for phase in phase_word:
        packet = integer_monomial(phase, -1)
        conjugate_packet = integer_monomial(-phase, -1)
        gram = integer_polynomial_multiply(packet, conjugate_packet, P)
        self_gram_cells.append(cyclotomic_reduce_13(gram))
    require(
        len(set(self_gram_cells)) == 1,
        "self-Gram varied across quotient cells",
    )
    require(
        self_gram_cells[0] == (1,) + (0,) * (P - 2),
        "self-Gram numerator was not one",
    )

    return (
        phase_counts,
        tuple(multiplier_counts),
        tuple(zero_reference_constants),
    )


def main():
    (
        fibre_checks,
        correlation_checks,
        fourier_checks,
        sharp_checks,
    ) = finite_shell_checks()
    exact_addresses, pushed_fibres = address_pushforward_checks()
    phase_counts, multiplier_counts, zero_constants = hostile_checks()

    print("THM-2420 SAME-SHELL CROSS-REFERENCE -- exact audit")
    print(f"primitive kernel reduction fibres checked={fibre_checks}")
    print(f"finite correlation packets checked={correlation_checks}")
    print(f"exact group-ring Fourier products checked={fourier_checks}")
    print(f"sharp max/energy equality controls={sharp_checks}")
    print("support / total / max / energy identities: PASS")
    print(
        "finite exact-address differences="
        f"{exact_addresses}, pushed fibres={pushed_fibres}"
    )
    print("coefficient-l1 / mod-M pushforward / aggregate amplitude: PASS")
    print(f"26-cell phase counts={phase_counts}")
    print(f"mode-2 multiplier counts={multiplier_counts}")
    print("Fhat(0)=0 and cyclotomic numerator at X=26 is nonzero: PASS")
    print(f"complete zero-reference numerators={zero_constants}")
    print("all 13 zero sectors and self-Gram are sideband-constant: PASS")
    print("canonical same-shell physical reference remains OPEN")
    print("THM-2420 exact companion PASS")


if __name__ == "__main__":
    main()
