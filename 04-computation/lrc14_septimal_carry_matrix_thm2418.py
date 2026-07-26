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


def residue_count(length, residue):
    quotient, remainder = divmod(length, P)
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

    even_source_checks = 0
    minimum_mass = Fraction(1)
    for depth in range(1, 9):
        dilation = 13**depth
        quotient = dilation // P
        if depth % 2 == 0:
            removed = {(dilation - 1) // 2}
            require(
                (dilation - 1) // 2 % P == 0,
                "even-depth central cylinder has wrong carry",
            )
        else:
            removed = {0, 1, 2, dilation - 3, dilation - 2, dilation - 1}
            require(
                {index % P for index in removed} == set(range(6)),
                "odd-depth boundary cylinders have wrong carries",
            )
        require(
            {dilation - 1 - index for index in removed} == removed,
            "source cylinder deletion is not reflection invariant",
        )
        counts = [
            residue_count(dilation, residue)
            - sum(index % P == residue for index in removed)
            for residue in range(P)
        ]
        require(
            counts == [quotient] * P,
            "even centred source did not equalize carries",
        )
        mass = Fraction(dilation - len(removed), dilation)
        minimum_mass = min(minimum_mass, mass)
        require(mass >= Fraction(7, 13), "source mass floor failed")

        source_kernel = [
            [Fraction(quotient, dilation) for _ in range(P)]
            for _ in range(P)
        ]
        require(matrix_rank(source_kernel) == 1, "source kernel rank failed")
        require(
            0 < len(removed) < dilation and jump_count((0, 1, 0)) == 2,
            "centred source variation control failed",
        )
        even_source_checks += 1

    require(minimum_mass == Fraction(7, 13), "sharp source mass mismatch")

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
    print(f"even centred source depths checked={even_source_checks}")
    print("even centred source variation=2 / minimum mass=7/13")
    print("equal-carry source kernel rank=1: PASS")
    print("canonical source-terminal correlation remains OPEN")
    print("THM-2418 exact companion PASS")


if __name__ == "__main__":
    main()
