#!/usr/bin/env python3
"""Dependency-free exact referee for candidate THM-2418."""

from fractions import Fraction
from itertools import product
from math import floor


P = 7


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def fractional_part(value):
    return value - floor(value)


def circle_distance(value):
    residue = fractional_part(value)
    return min(residue, 1 - residue)


def word(speed, width, base):
    return {
        root
        for root in range(P)
        if circle_distance(Fraction(speed) * (base + root) / P)
        < Fraction(width, 14)
    }


def affine_image(values, sign, carry):
    return {(sign * (value - carry)) % P for value in values}


def matrix_multiply(left, right):
    size = len(left)
    return [
        [
            sum(
                (left[row][middle] * right[middle][column]
                 for middle in range(size)),
                Fraction(0),
            )
            for column in range(size)
        ]
        for row in range(size)
    ]


def transpose(matrix):
    return [list(column) for column in zip(*matrix)]


def matrix_vector(matrix, vector):
    return [
        sum(
            (entry * value for entry, value in zip(row, vector)),
            Fraction(0),
        )
        for row in matrix
    ]


def matrix_rank(matrix):
    work = [list(row) for row in matrix]
    rows = len(work)
    columns = len(work[0]) if rows else 0
    pivot_row = 0
    for column in range(columns):
        pivot = next(
            (row for row in range(pivot_row, rows) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        scale = work[pivot_row][column]
        work[pivot_row] = [value / scale for value in work[pivot_row]]
        for row in range(rows):
            if row == pivot_row:
                continue
            factor = work[row][column]
            if factor:
                work[row] = [
                    value - factor * base
                    for value, base in zip(work[row], work[pivot_row])
                ]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def residue_count(length, residue, modulus=P):
    quotient, remainder = divmod(length, modulus)
    return quotient + int(residue < remainder)


def cyclotomic_reduce(coefficients):
    require(len(coefficients) == P, "wrong cyclotomic vector length")
    top = coefficients[-1]
    return tuple(value - top for value in coefficients[:-1])


def digit_fourier_vector(length, character):
    counts = [residue_count(length, residue) for residue in range(P)]
    vector = [0] * P
    for residue, count in enumerate(counts):
        vector[(character * residue) % P] += count
    return vector


def value_fourier_vector(values, character):
    vector = [Fraction(0)] * P
    for residue, value in enumerate(values):
        vector[(character * residue) % P] += value
    return vector


def jump_count(values):
    return sum(
        values[index] != values[index - 1]
        for index in range(len(values))
    )


def multiplicative_order(base, modulus):
    require(modulus > 1, "order modulus must exceed one")
    value = 1
    for exponent in range(1, modulus + 1):
        value = value * base % modulus
        if value == 1:
            return exponent
    raise RuntimeError("multiplicative order not found")


def quotient_phase_prefix(limit, denominator_tail):
    """Counts floor(j/D0) modulo seven for 0<=j<limit."""

    period = P * denominator_tail
    cycles, remainder = divmod(limit, period)
    histogram = [cycles * denominator_tail] * P
    for index in range(remainder):
        histogram[(index // denominator_tail) % P] += 1
    return histogram


def quotient_phase_interval(start, stop, denominator_tail):
    left = quotient_phase_prefix(start, denominator_tail)
    right = quotient_phase_prefix(stop, denominator_tail)
    return [right[index] - left[index] for index in range(P)]


def weighted_carry_masses(source_values, terminal_values, dilation):
    """Exact masses for G(y)Q({Ry}) on the seven carry classes."""

    source_denominator = len(source_values)
    terminal_denominator = len(terminal_values)
    breakpoints = {Fraction(0), Fraction(1)}
    breakpoints.update(
        Fraction(index, terminal_denominator)
        for index in range(1, terminal_denominator)
    )
    for index in range(1, source_denominator):
        scaled = Fraction(dilation * index, source_denominator)
        breakpoint = scaled - floor(scaled)
        if breakpoint:
            breakpoints.add(breakpoint)
    ordered = sorted(breakpoints)

    masses = [Fraction(0)] * P
    for left, right in zip(ordered, ordered[1:]):
        midpoint = (left + right) / 2
        terminal_index = floor(midpoint * terminal_denominator)
        terminal_value = terminal_values[terminal_index]
        length = right - left
        for prefix in range(dilation):
            source_index = floor(
                source_denominator * (prefix + midpoint) / dilation
            )
            masses[prefix % P] += (
                length
                * source_values[source_index]
                * terminal_value
                / dilation
            )
    return tuple(masses)


def centered_scaled_histogram(masses, dilation):
    mean = sum(masses, Fraction(0)) / P
    return tuple(dilation * (mass - mean) for mass in masses)


def main():
    cocycle_checks = 0
    for denominator in (17, 19, 31, 101):
        for numerator in range(1, denominator):
            base = Fraction(numerator, denominator)
            for depth in range(1, 6):
                dilation = 13**depth
                reduced = fractional_part(dilation * base)
                carry = floor(dilation * base)
                sign = (-1) ** depth
                for speed in (1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12):
                    for width in (1, 2):
                        require(
                            word(dilation * speed, width, base)
                            == affine_image(
                                word(speed, width, reduced),
                                sign,
                                carry % P,
                            ),
                            "affine word cocycle failed",
                        )
                        cocycle_checks += 1
    require(cocycle_checks == 18_040, "cocycle census mismatch")

    digit_words = 0
    for depth in range(1, 6):
        for digits in product(range(13), repeat=depth):
            prefix = sum(
                digit * 13 ** (depth - 1 - index)
                for index, digit in enumerate(digits)
            )
            alternating = sum(
                digit * (-1) ** (depth - 1 - index)
                for index, digit in enumerate(digits)
            )
            require(
                prefix % P == alternating % P,
                "alternating digit carry failed",
            )
            digit_words += 1
    require(digit_words == 402_233, "digit-word census mismatch")

    pi = [[Fraction(1, P) for _ in range(P)] for _ in range(P)]
    identity = [
        [Fraction(row == column) for column in range(P)]
        for row in range(P)
    ]
    reflection = [
        [Fraction((row + column) % P == 6) for column in range(P)]
        for row in range(P)
    ]
    transition = [
        [
            Fraction(2 - int((row + column) % P == 6), 13)
            for column in range(P)
        ]
        for row in range(P)
    ]

    power = identity
    matrix_depths = 0
    terminal_support_checks = 0
    terminal_charged_checks = 0
    for depth in range(1, 9):
        dilation = 13**depth
        sign = (-1) ** depth
        power = matrix_multiply(power, transition)
        direct = [
            [
                Fraction(
                    residue_count(
                        dilation,
                        (column - sign * row) % P,
                    ),
                    dilation,
                )
                for column in range(P)
            ]
            for row in range(P)
        ]
        closed = [
            [
                pi[row][column]
                + (
                    Fraction(1, dilation)
                    * (identity[row][column] - pi[row][column])
                    if depth % 2 == 0
                    else -Fraction(1, dilation)
                    * (reflection[row][column] - pi[row][column])
                )
                for column in range(P)
            ]
            for row in range(P)
        ]
        require(power == direct == closed, "carry matrix formula failed")
        gram = matrix_multiply(closed, transpose(closed))
        expected_gram = [
            [
                pi[row][column]
                + Fraction(1, dilation**2)
                * (identity[row][column] - pi[row][column])
                for column in range(P)
            ]
            for row in range(P)
        ]
        require(gram == expected_gram, "carry Gram formula failed")
        require(matrix_rank(closed) == P, "raw carry matrix lost rank")

        for character in range(1, P):
            actual = cyclotomic_reduce(
                digit_fourier_vector(dilation, character)
            )
            expected_vector = [0] * P
            if depth % 2 == 0:
                expected_vector[0] = 1
            else:
                expected_vector[(-character) % P] = -1
            expected = cyclotomic_reduce(expected_vector)
            require(actual == expected, "raw digit Fourier law failed")

        for mask in range(1, 1 << P):
            terminal = tuple(
                Fraction((mask >> column) & 1)
                for column in range(P)
            )
            filtered = [
                [
                    closed[row][column] * terminal[column]
                    for column in range(P)
                ]
                for row in range(P)
            ]
            require(
                matrix_rank(filtered) == sum(bool(value) for value in terminal),
                "terminal support/rank identity failed",
            )
            marginal = matrix_vector(closed, terminal)
            require(
                (len(set(marginal)) > 1)
                == (len(set(terminal)) > 1),
                "terminal nonflatness equivalence failed",
            )
            if len(set(terminal)) > 1:
                for character in range(1, P):
                    require(
                        any(
                            cyclotomic_reduce(
                                value_fourier_vector(terminal, character)
                            )
                        ),
                        "nonconstant rational terminal lost a colour",
                    )
                    require(
                        any(
                            cyclotomic_reduce(
                                value_fourier_vector(marginal, character)
                            )
                        ),
                        "filtered marginal lost a charged colour",
                    )
                    terminal_charged_checks += 1
            terminal_support_checks += 1
        require(
            matrix_vector(closed, (Fraction(1, 7),) * P)
            == [Fraction(1, 7)] * P,
            "flat D_7 terminal hostile failed",
        )
        matrix_depths += 1

    flat_terminal = (Fraction(1, 7),) * P
    require(matrix_vector(power, flat_terminal) == list(flat_terminal),
            "final flat terminal check failed")

    rank_one_terminal = [Fraction(0)] * P
    rank_one_terminal[3] = Fraction(1, 2)
    first_step_masses = matrix_vector(transition, rank_one_terminal)
    require(
        first_step_masses
        == [Fraction(1, 13)] * 3
        + [Fraction(1, 26)]
        + [Fraction(1, 13)] * 3,
        "one-cylinder hostile masses failed",
    )
    require(
        matrix_rank(
            [
                [
                    transition[row][column] * rank_one_terminal[column]
                    for column in range(P)
                ]
                for row in range(P)
            ]
        )
        == 1,
        "one-cylinder terminal matrix did not have rank one",
    )

    fixed_source_checks = 0
    for depth in range(1, 9):
        dilation = 13**depth
        unit = 13 ** (depth - 1)
        lower = 3 * unit
        upper = 10 * unit
        require(
            lower + upper == dilation,
            "fixed source block is not reflection invariant",
        )
        counts = [
            residue_count(upper, residue) - residue_count(lower, residue)
            for residue in range(P)
        ]
        require(
            counts == [unit] * P,
            "fixed source did not equalize carries",
        )
        require(
            Fraction(upper - lower, dilation) == Fraction(7, 13),
            "fixed source mass failed",
        )
        source_kernel = [
            [Fraction(unit, dilation) for _ in range(P)]
            for _ in range(P)
        ]
        require(matrix_rank(source_kernel) == 1, "source kernel rank failed")
        require(
            Fraction(unit, dilation) == Fraction(1, 13),
            "fixed source kernel normalization failed",
        )
        require(
            jump_count(tuple(int(3 <= cell < 10) for cell in range(13)))
            == 2,
            "fixed source BV-two control failed",
        )
        fixed_source_checks += 1

    danger_source_checks = 0
    require(
        Fraction(81, 169) - Fraction(13, 28)
        == Fraction(15, 28) - Fraction(88, 169)
        == Fraction(71, 4732),
        "D2 source clearance failed",
    )
    require(
        floor(Fraction(13 * 81, 169))
        == floor(Fraction(13 * 88, 169))
        == 6,
        "D2 source left its one predecessor sheet",
    )
    require(
        jump_count(tuple(int(81 <= cell < 88) for cell in range(169)))
        == 2,
        "D2 source variation failed",
    )
    for depth in range(2, 8):
        unit = 13 ** (depth - 2)
        lower = 81 * unit
        upper = 88 * unit
        counts = [
            residue_count(upper, residue) - residue_count(lower, residue)
            for residue in range(P)
        ]
        require(counts == [unit] * P, "D2 source carry balance failed")
        require(
            Fraction(upper - lower, 13**depth) == Fraction(7, 169),
            "D2 source mass failed",
        )
        danger_source_checks += 1

    universal_block_checks = 0
    for base, prime in ((5, 3), (7, 5), (11, 3), (13, 7), (17, 5)):
        anchor = (base - prime) // 2
        for depth in range(1, 5):
            dilation = base**depth
            unit = base ** (depth - 1)
            lower = anchor * unit
            upper = (anchor + prime) * unit
            counts = [
                residue_count(upper, residue, prime)
                - residue_count(lower, residue, prime)
                for residue in range(prime)
            ]
            require(
                counts == [unit] * prime,
                f"universal carry block failed at B={base}, p={prime}, k={depth}",
            )
            require(
                lower + upper == dilation,
                f"universal centred block failed at B={base}, p={prime}, k={depth}",
            )
            require(
                Fraction(upper - lower, dilation) == Fraction(prime, base),
                "universal block mass failed",
            )
            universal_block_checks += 1

    cylinder_profiles = 0
    uniform_cylinder_profiles = 0
    for mask in range(1 << 13):
        values = [(mask >> cell) & 1 for cell in range(13)]
        histogram = [
            sum(values[cell] for cell in range(13) if cell % P == residue)
            for residue in range(P)
        ]
        uniform = len(set(histogram)) == 1
        uniform_cylinder_profiles += int(uniform)
        for character in range(1, P):
            reduced = cyclotomic_reduce(
                value_fourier_vector(histogram, character)
            )
            require(
                (not any(reduced)) == uniform,
                "rational cylinder all-or-flat classification failed",
            )
        cylinder_profiles += 1

    tail_formula_checks = 0
    tail_controls = (
        [int(3 <= cell < 10) for cell in range(13)],
        [int(cell in (0, 1, 4, 9)) for cell in range(13)],
        [Fraction((cell * cell + 3) % 11, 11) for cell in range(13)],
    )
    for values in tail_controls:
        for extra_depth in range(6):
            tail_length = 13**extra_depth
            for character in range(1, P):
                direct = [Fraction(0)] * P
                for cell, value in enumerate(values):
                    for residue in range(P):
                        count = residue_count(tail_length, residue)
                        exponent = (
                            character * (tail_length * cell + residue)
                        ) % P
                        direct[exponent] += value * count

                predicted = [Fraction(0)] * P
                if extra_depth % 2 == 0:
                    for cell, value in enumerate(values):
                        predicted[(character * cell) % P] += value
                else:
                    for cell, value in enumerate(values):
                        predicted[(-character * (cell + 1)) % P] -= value
                require(
                    cyclotomic_reduce(direct)
                    == cyclotomic_reduce(predicted),
                    "fixed-cylinder tail Fourier formula failed",
                )
                tail_formula_checks += 1

    rational_period_checks = 0
    rational_controls = (
        (0, 2, (Fraction(1), Fraction(0))),
        (0, 5, tuple(Fraction((3 * cell + 1) % 7, 7) for cell in range(5))),
        (
            1,
            3,
            tuple(
                Fraction(int(cell in (0, 5, 8, 9, 13, 16, 18)))
                for cell in range(39)
            ),
        ),
        (
            1,
            2,
            tuple(Fraction((cell * cell + 2) % 5, 5) for cell in range(26)),
        ),
        (
            1,
            5,
            tuple(Fraction((2 * cell + 3) % 11, 11) for cell in range(65)),
        ),
    )
    for base_depth, denominator_tail, values in rational_controls:
        require(
            len(values) == 13**base_depth * denominator_tail,
            "rational control denominator mismatch",
        )
        period = multiplicative_order(13, P * denominator_tail)
        period_bank = {}
        for extra_depth in range(2 * period):
            descendants = 13**extra_depth
            masses = [Fraction(0)] * P
            for cell, value in enumerate(values):
                histogram = quotient_phase_interval(
                    cell * descendants,
                    (cell + 1) * descendants,
                    denominator_tail,
                )
                for residue, count in enumerate(histogram):
                    masses[residue] += value * count / denominator_tail

            uniform = len(set(masses)) == 1
            if base_depth == 1 and denominator_tail == 3:
                alternating_expected = (
                    (Fraction(1, 3),) * P,
                    (
                        Fraction(5),
                        Fraction(19, 3),
                        Fraction(6),
                        Fraction(17, 3),
                        Fraction(4),
                        Fraction(4, 3),
                        Fraction(2),
                    ),
                    (Fraction(169, 3),) * P,
                    (
                        Fraction(733),
                        Fraction(2203, 3),
                        Fraction(734),
                        Fraction(2201, 3),
                        Fraction(732),
                        Fraction(2188, 3),
                        Fraction(730),
                    ),
                )
                require(
                    tuple(masses) == alternating_expected[extra_depth],
                    "Boolean period-two carry histogram failed",
                )
                require(
                    uniform == (extra_depth % 2 == 0),
                    "Boolean flat/all-six alternation failed",
                )
            charged = []
            for character in range(1, P):
                vector = value_fourier_vector(masses, character)
                reduced = cyclotomic_reduce(vector)
                charged.append(reduced)
                require(
                    (not any(reduced)) == uniform,
                    "rational-step one-colour/all-flat gate failed",
                )
                rational_period_checks += 1

            key = extra_depth % period
            signature = tuple(charged)
            if key in period_bank:
                require(
                    signature == period_bank[key],
                    "rational-step charged scale periodicity failed",
                )
            else:
                period_bank[key] = signature

    weighted_tail_checks = 0
    fixed_source = tuple(Fraction(int(3 <= cell < 10)) for cell in range(13))
    terminal_profile = (
        Fraction(1),
        Fraction(0),
        Fraction(2, 3),
        Fraction(1, 5),
        Fraction(4, 7),
    )
    for depth in range(1, 5):
        masses = weighted_carry_masses(
            fixed_source,
            terminal_profile,
            13**depth,
        )
        require(
            len(set(masses)) == 1,
            "fixed source failed against a nonconstant terminal weight",
        )
        weighted_tail_checks += 1

    general_source = tuple(
        Fraction((cell * cell + 3 * cell + 2) % 7, 7)
        for cell in range(26)
    )
    general_terminal = tuple(
        Fraction((2 * cell + 1) % 5, 5)
        for cell in range(5)
    )
    parity_signatures = {}
    for depth in range(1, 5):
        dilation = 13**depth
        masses = weighted_carry_masses(
            general_source,
            general_terminal,
            dilation,
        )
        signature = centered_scaled_histogram(masses, dilation)
        parity = depth % 2
        if parity in parity_signatures:
            require(
                signature == parity_signatures[parity],
                "weighted source-terminal period-two law failed",
            )
        else:
            parity_signatures[parity] = signature
        weighted_tail_checks += 1

    raw_scaled = []
    centered_scaled = []
    for depth in (1, 3):
        dilation = 13**depth
        masses = weighted_carry_masses(
            (Fraction(1),),
            (Fraction(1),),
            dilation,
        )
        raw_scaled.append(tuple(dilation * mass for mass in masses))
        centered_scaled.append(centered_scaled_histogram(masses, dilation))
        weighted_tail_checks += 1
    require(raw_scaled[0] != raw_scaled[1], "raw scaled masses became periodic")
    require(
        centered_scaled[0] == centered_scaled[1],
        "centered scaled masses lost periodicity",
    )

    print("THM-2418 SEPTIMAL CARRY MATRIX -- exact audit")
    print(f"rational affine word cocycle checks={cocycle_checks}")
    print(f"base-thirteen digit words checked={digit_words}")
    print(f"carry matrix depths checked={matrix_depths} (k=1..8)")
    print("K_k=P^k and alternating closed forms: PASS")
    print("K_k K_k^T=Pi+13^(-2k)(I-Pi): PASS")
    print("raw rank=7 / charged singular values=13^(-k): PASS")
    print("raw digit Fourier phases even=1 / odd=-zeta^(-e): PASS")
    print(f"terminal support/rank checks={terminal_support_checks}")
    print(f"terminal charged-colour checks={terminal_charged_checks}")
    print("terminal nonconstant iff every rational charged colour survives: PASS")
    print("flat D_7 terminal profile=(1/7)^7: PASS")
    print("one-cylinder terminal profile=(1/2)e_3 / rank=1: PASS")
    print(f"fixed [3/13,10/13) source depths checked={fixed_source_checks}")
    print("fixed source variation=2 / mass=7/13 / kernel=(1/13)J")
    print(
        "fixed one-sheet D2 source depths checked="
        f"{danger_source_checks} (k=2..7)"
    )
    print("D2 source variation=2 / mass=7/169 / kernel=(1/169)J")
    print(f"universal odd (B,p) equal-carry blocks checked={universal_block_checks}")
    print(
        "depth-one Boolean cylinder profiles="
        f"{cylinder_profiles} (flat tails={uniform_cylinder_profiles})"
    )
    print(f"fixed-cylinder tail Fourier formulas checked={tail_formula_checks}")
    print("one histogram classifies all deeper clocks: PASS")
    print(f"rational finite-step scale-period checks={rational_period_checks}")
    print("one ord_(7D0)(13) period classifies every rational tail: PASS")
    print("D0=3 Boolean period-2 flat/all-six alternation: PASS")
    print(f"fixed source-terminal weighted checks={weighted_tail_checks}")
    print("Q denominator-free period / centred-defect typing: PASS")
    print("canonical source-terminal correlation remains OPEN")
    print("THM-2418 exact companion PASS")


if __name__ == "__main__":
    main()
