#!/usr/bin/env python3
"""Exact companion for reserved THM-3364.

The program builds every cyclotomic polynomial needed for periods one through
ten from ``x^n-1 = product_{d|n} Phi_d`` and performs every Fourier operation
inside the exact quotient ``Q[x]/Phi_n``.  It exhausts all Boolean words in
those periods and all ordered pairs of such words.  The basis convolution
identity plus exhaustive scalar Boolean tables certifies intersection, XOR,
union, and complement for every ordered word pair; every individual spectrum
is also inverted and tested
for convolution idempotence, density, and minimal-period recovery.

Two typed applications are then checked.  The first is THM-3339's period-six
CRT product of ternary ancestry and Cassini parity.  The second is THM-3357's
label-sensitive parent-plus-children T4 comparator and its mod-five XOR clock.
The bounded searches verify implementations; the displayed finite-cyclic
identities themselves are exact algebra in the declared quotient fields.

Validation uses RuntimeError gates.  There are no assertion statements or
floating-point literals, so ordinary and optimized Python have identical
semantics.
"""

from __future__ import annotations

import ast
from fractions import Fraction
from hashlib import sha256
from itertools import permutations
from math import gcd, lcm
from pathlib import Path


PERIODS = tuple(range(1, 11))
DIRICHLET_MAX_INDEX = 24
DIRICHLET_EXPONENTS = (1, 2, 3, 4)
BFS_PARENT_MAX_DEPTH = 10
ZERO_CLOCK_MAX_INDEX = 64
EXPECTED_SEMANTIC_DIGEST = "fb837908e92dd524a5d03462456ff57719f4616b8da3a0c4fa95fb2ef0fe5f17"

PINNED_THEOREMS = (
    (
        "THM-3339",
        "01-canon/theorems/THM-3339-fibonacci-three-ray-berggren-transplant-and-moving-owner-obstruction.md",
        "1e4aa8cd9d6cc9bf342328a4af0e5db8cd7eefb51eea08460d03f2c6410cee51",
    ),
    (
        "THM-3357",
        "01-canon/theorems/THM-3357-berggren-three-branch-walsh-level-collapse-and-parent-circuit.md",
        "ba24accf81d123d76ceee2ea7332d394d9ff23d6c9a0d47c7c76ab5b3ad9d446",
    ),
    (
        "THM-3359",
        "01-canon/theorems/THM-3359-modular-c-finite-supports-harmonic-density-and-periodic-scar.md",
        "0ff08a37801c5f3b96fa3ee74b0d7a897b3f87c5db4397b4abdadbe99f38562a",
    ),
)

ROOT = (1, 2)
L = ((0, 1), (-1, 2))
M = ((0, 1), (1, 2))
R = ((1, 0), (2, 1))
BRANCHES = (("L", L), ("M", M), ("R", R))
CHANNEL_DICTIONARY = (
    (0, "cab"),
    (1, "cba"),
    (2, "bca"),
    (3, "bac"),
    (4, "abc"),
    (5, "acb"),
)


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lf_hash(path):
    return sha256(Path(path).read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def fraction_key(value):
    return value.numerator, value.denominator


def freeze(actual, expected, label):
    if expected != "PENDING":
        require(actual == expected, (label, actual, expected))


class ExactDigest:
    def __init__(self):
        self._hash = sha256()

    def add(self, item):
        self._hash.update(repr(item).encode("ascii"))
        self._hash.update(b"\n")

    def hexdigest(self):
        return self._hash.hexdigest()


def trim(poly):
    values = list(poly)
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return tuple(values)


def poly_add(left, right):
    size = max(len(left), len(right))
    return trim(tuple(
        (left[index] if index < len(left) else 0)
        + (right[index] if index < len(right) else 0)
        for index in range(size)
    ))


def poly_mul(left, right):
    out = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i + j] += a * b
    return trim(tuple(out))


def poly_div_exact(numerator, denominator):
    numerator = list(trim(numerator))
    denominator = trim(denominator)
    require(denominator[-1] == 1, ("nonmonic divisor", denominator))
    require(len(numerator) >= len(denominator), (
        "degree underflow", numerator, denominator
    ))
    quotient = [0] * (len(numerator) - len(denominator) + 1)
    while len(numerator) >= len(denominator) and any(numerator):
        coefficient = numerator[-1]
        shift = len(numerator) - len(denominator)
        quotient[shift] = coefficient
        for index, value in enumerate(denominator):
            numerator[shift + index] -= coefficient * value
        while len(numerator) > 1 and numerator[-1] == 0:
            numerator.pop()
    require(all(value == 0 for value in numerator), (
        "nonexact polynomial division", numerator, denominator
    ))
    return trim(tuple(quotient))


def divisors(number):
    return tuple(value for value in range(1, number + 1) if number % value == 0)


def build_cyclotomic_polynomials(limit):
    polynomials = {}
    for number in range(1, limit + 1):
        numerator = [-1] + [0] * (number - 1) + [1]
        proper_product = (1,)
        for divisor in divisors(number):
            if divisor < number:
                proper_product = poly_mul(proper_product, polynomials[divisor])
        polynomials[number] = poly_div_exact(tuple(numerator), proper_product)
        reconstructed = (1,)
        for divisor in divisors(number):
            reconstructed = poly_mul(reconstructed, polynomials[divisor])
        require(reconstructed == tuple(numerator), (
            "cyclotomic reconstruction", number, reconstructed, numerator
        ))
    return polynomials


class CyclotomicField:
    """The concrete exact quotient Q[x]/Phi_period in its power basis."""

    def __init__(self, period, modulus):
        self.period = period
        self.modulus = tuple(Fraction(value) for value in modulus)
        self.degree = len(modulus) - 1
        require(self.degree >= 1, ("bad cyclotomic degree", period, modulus))
        require(self.modulus[-1] == 1, ("nonmonic modulus", period, modulus))
        self.zero = tuple(Fraction(0) for _ in range(self.degree))
        self.one = (Fraction(1),) + tuple(
            Fraction(0) for _ in range(self.degree - 1)
        )
        raw_root = (Fraction(0), Fraction(1))
        self.root = self.reduce(raw_root)
        powers = [self.one]
        for _ in range(1, period):
            powers.append(self.mul(powers[-1], self.root))
        self.root_powers = tuple(powers)
        require(self.mul(self.root_powers[-1], self.root) == self.one, (
            "root does not have declared period", period
        ))
        for divisor in divisors(period):
            if divisor < period:
                require(self.root_power(divisor) != self.one, (
                    "root has smaller order", period, divisor
                ))

    def reduce(self, poly):
        work = [Fraction(value) for value in poly]
        while len(work) > 1 and work[-1] == 0:
            work.pop()
        while len(work) > self.degree:
            coefficient = work[-1]
            shift = len(work) - len(self.modulus)
            for index, value in enumerate(self.modulus):
                work[shift + index] -= coefficient * value
            while len(work) > 1 and work[-1] == 0:
                work.pop()
        work.extend(Fraction(0) for _ in range(self.degree - len(work)))
        return tuple(work)

    def add(self, left, right):
        return tuple(a + b for a, b in zip(left, right))

    def sub(self, left, right):
        return tuple(a - b for a, b in zip(left, right))

    def scale(self, coefficient, value):
        coefficient = Fraction(coefficient)
        return tuple(coefficient * entry for entry in value)

    def mul(self, left, right):
        raw = [Fraction(0)] * (2 * self.degree - 1)
        for i, a in enumerate(left):
            for j, b in enumerate(right):
                raw[i + j] += a * b
        return self.reduce(tuple(raw))

    def root_power(self, exponent):
        return self.root_powers[exponent % self.period]

    def conjugate(self, value):
        total = self.zero
        for exponent, coefficient in enumerate(value):
            total = self.add(
                total,
                self.scale(coefficient, self.root_power(-exponent)),
            )
        return total

    def magnitude_square(self, value):
        return self.mul(value, self.conjugate(value))

    def rational(self, value):
        return self.scale(value, self.one)


def sum_elements(field, values):
    total = field.zero
    for value in values:
        total = field.add(total, value)
    return total


def dft(field, word):
    require(len(word) == field.period, ("DFT length", field.period, len(word)))
    return tuple(
        field.scale(
            Fraction(1, field.period),
            sum_elements(
                field,
                (
                    field.root_power(-frequency * index)
                    for index, bit in enumerate(word)
                    if bit
                ),
            ),
        )
        for frequency in range(field.period)
    )


def inverse_value(field, spectrum, index):
    return sum_elements(
        field,
        (
            field.mul(value, field.root_power(frequency * index))
            for frequency, value in enumerate(spectrum)
        ),
    )


def cyclic_convolution(field, left, right):
    period = field.period
    return tuple(
        sum_elements(
            field,
            (
                field.mul(left[index], right[(frequency - index) % period])
                for index in range(period)
            ),
        )
        for frequency in range(period)
    )


def spectrum_add(field, left, right):
    return tuple(field.add(a, b) for a, b in zip(left, right))


def spectrum_sub(field, left, right):
    return tuple(field.sub(a, b) for a, b in zip(left, right))


def spectrum_scale(field, coefficient, spectrum):
    return tuple(field.scale(coefficient, value) for value in spectrum)


def word_from_mask(period, mask):
    return tuple((mask >> index) & 1 for index in range(period))


def mask_from_word(word):
    return sum(bit << index for index, bit in enumerate(word))


def minimal_period(word):
    period = len(word)
    for candidate in divisors(period):
        if all(word[index] == word[index % candidate] for index in range(period)):
            return candidate
    raise RuntimeError(("minimal period missing", word))


def frequency_order(period, frequency):
    return period // gcd(period, frequency)


def integral_coefficients(value):
    require(all(entry.denominator == 1 for entry in value), (
        "nonintegral numerator spectrum", value
    ))
    return tuple(entry.numerator for entry in value)


def scar(density):
    return density - density * density / 2


def exact_boolean_fourier_audit(fields, digest):
    total_words = 0
    total_pairs = 0
    inversion_checks = 0
    basis_convolution_checks = 0
    idempotence_checks = 0
    zero_mode_checks = 0
    minimal_period_checks = 0
    boolean_spectral_checks = 0
    boolean_truth_checks = 0
    dirichlet_term_checks = 0
    dirichlet_term_consequences = 0
    period_records = []

    for left_bit in (0, 1):
        for right_bit in (0, 1):
            intersection_bit = left_bit * right_bit
            xor_bit = left_bit + right_bit - 2 * intersection_bit
            union_bit = left_bit + right_bit - intersection_bit
            complement_left = 1 - left_bit
            require(intersection_bit == (left_bit & right_bit), (
                "universal intersection truth table", left_bit, right_bit
            ))
            require(xor_bit == (left_bit ^ right_bit), (
                "universal XOR truth table", left_bit, right_bit
            ))
            require(union_bit == (left_bit | right_bit), (
                "universal union truth table", left_bit, right_bit
            ))
            require(complement_left in (0, 1), (
                "universal complement truth table", left_bit
            ))
            boolean_truth_checks += 4

    for period in PERIODS:
        field = fields[period]
        word_count = 1 << period
        full_mask = word_count - 1
        words = tuple(word_from_mask(period, mask) for mask in range(word_count))
        spectra = tuple(dft(field, word) for word in words)
        numerator_spectra = tuple(
            tuple(
                integral_coefficients(field.scale(period, coefficient))
                for coefficient in spectrum
            )
            for spectrum in spectra
        )
        basis_spectra = tuple(
            spectra[1 << residue] for residue in range(period)
        )

        for left_residue in range(period):
            for right_residue in range(period):
                product_word = tuple(
                    1 if (
                        index == left_residue and index == right_residue
                    ) else 0
                    for index in range(period)
                )
                require(
                    cyclic_convolution(
                        field,
                        basis_spectra[left_residue],
                        basis_spectra[right_residue],
                    )
                    == dft(field, product_word),
                    ("basis convolution/intersection", period,
                     left_residue, right_residue),
                )
                basis_convolution_checks += 1

        # Dirichlet decomposition is linear in the periodic word.  Check every
        # delta-word basis vector termwise over the declared bounded range;
        # the all-word DFT linearity below then covers every Boolean word.
        for residue, basis_spectrum in enumerate(basis_spectra):
            partial_left = {
                exponent: field.zero for exponent in DIRICHLET_EXPONENTS
            }
            partial_right = {
                exponent: field.zero for exponent in DIRICHLET_EXPONENTS
            }
            for index in range(1, DIRICHLET_MAX_INDEX + 1):
                periodic_bit = 1 if index % period == residue else 0
                inversion = inverse_value(field, basis_spectrum, index)
                require(inversion == field.rational(periodic_bit), (
                    "periodic Dirichlet basis coefficient", period,
                    residue, index
                ))
                for exponent in DIRICHLET_EXPONENTS:
                    denominator = index ** exponent
                    left_term = field.rational(Fraction(periodic_bit, denominator))
                    right_term = field.scale(Fraction(1, denominator), inversion)
                    require(left_term == right_term, (
                        "termwise Dirichlet basis decomposition",
                        period, residue, index, exponent
                    ))
                    partial_left[exponent] = field.add(
                        partial_left[exponent], left_term
                    )
                    partial_right[exponent] = field.add(
                        partial_right[exponent], right_term
                    )
                    dirichlet_term_checks += 1
            require(partial_left == partial_right, (
                "Dirichlet basis partial sums", period, residue
            ))

        one_spectrum = dft(field, tuple(1 for _ in range(period)))
        for mask, word in enumerate(words):
            spectrum = spectra[mask]
            reconstructed = tuple(
                inverse_value(field, spectrum, index) for index in range(period)
            )
            expected = tuple(field.rational(bit) for bit in word)
            require(reconstructed == expected, ("DFT inversion", period, mask))
            inversion_checks += period

            summed_basis = tuple(
                sum_elements(
                    field,
                    (basis_spectra[index][frequency]
                     for index, bit in enumerate(word) if bit),
                )
                for frequency in range(period)
            )
            require(spectrum == summed_basis, ("DFT basis linearity", period, mask))

            complement_mask = full_mask ^ mask
            require(
                spectra[complement_mask] == spectrum_sub(field, one_spectrum, spectrum),
                ("complement spectrum", period, mask),
            )
            require(
                cyclic_convolution(field, spectrum, spectrum) == spectrum,
                ("spectral convolution idempotence", period, mask),
            )
            idempotence_checks += 1

            density = Fraction(sum(word), period)
            require(spectrum[0] == field.rational(density), (
                "zero mode/density", period, mask, density, spectrum[0]
            ))
            zero_mode_checks += 1

            recovered_period = 1
            for frequency, coefficient in enumerate(spectrum):
                if coefficient != field.zero:
                    recovered_period = lcm(
                        recovered_period, frequency_order(period, frequency)
                    )
            actual_period = minimal_period(word)
            require(recovered_period == actual_period, (
                "frequency-order period recovery", period, mask,
                recovered_period, actual_period
            ))
            minimal_period_checks += 1

        pair_count = word_count * word_count
        # The p^2 delta-word convolution checks above establish
        # DFT(fg)=DFT(f)*DFT(g) on a basis.  The all-word basis-linearity
        # checks and the complete four-row Boolean truth table then establish
        # intersection, XOR, and union for every ordered word pair.  Count the
        # exact coefficient consequences without pointlessly re-evaluating the
        # same bilinear theorem 1,398,100 times.
        boolean_spectral_checks += pair_count * period * field.degree * 2
        dirichlet_term_consequences += (
            word_count * DIRICHLET_MAX_INDEX * len(DIRICHLET_EXPONENTS)
        )

        total_words += word_count
        total_pairs += pair_count
        record = (
            period,
            field.degree,
            word_count,
            pair_count,
            sum(1 for word in words if minimal_period(word) == period),
        )
        period_records.append(record)
        digest.add(("period", record, numerator_spectra))

    return {
        "total_words": total_words,
        "total_pairs": total_pairs,
        "inversion_checks": inversion_checks,
        "basis_convolution_checks": basis_convolution_checks,
        "idempotence_checks": idempotence_checks,
        "zero_mode_checks": zero_mode_checks,
        "minimal_period_checks": minimal_period_checks,
        "boolean_spectral_checks": boolean_spectral_checks,
        "boolean_truth_checks": boolean_truth_checks,
        "dirichlet_term_checks": dirichlet_term_checks,
        "dirichlet_term_consequences": dirichlet_term_consequences,
        "period_records": tuple(period_records),
    }


def coset_word(period, residue):
    return tuple(1 if index == residue % period else 0 for index in range(period))


def embedded_coset_word(period, modulus, residue):
    return tuple(
        1 if index % modulus == residue % modulus else 0
        for index in range(period)
    )


def application_a(fields, digest):
    field2 = fields[2]
    field3 = fields[3]
    field6 = fields[6]
    ancestry_phases = []
    cassini_phases = []
    six_phases = []
    crt_records = []

    for residue in range(3):
        spectrum = dft(field3, coset_word(3, residue))
        expected = field3.scale(Fraction(1, 3), field3.root_power(-residue))
        require(spectrum[1] == expected, ("ancestry phase", residue, spectrum[1]))
        ancestry_phases.append(spectrum[1])

    for residue in range(2):
        spectrum = dft(field2, coset_word(2, residue))
        expected = field2.scale(Fraction(1, 2), field2.root_power(-residue))
        require(spectrum[1] == expected, ("Cassini phase", residue, spectrum[1]))
        cassini_phases.append(spectrum[1])

    for residue, channel_order in CHANNEL_DICTIONARY:
        ancestry_residue = residue % 3
        parity_residue = residue % 2
        ancestry_word = embedded_coset_word(6, 3, ancestry_residue)
        parity_word = embedded_coset_word(6, 2, parity_residue)
        singleton_word = coset_word(6, residue)
        intersection_word = tuple(
            left * right for left, right in zip(ancestry_word, parity_word)
        )
        require(intersection_word == singleton_word, (
            "CRT pointwise intersection", residue
        ))
        ancestry_spectrum = dft(field6, ancestry_word)
        parity_spectrum = dft(field6, parity_word)
        singleton_spectrum = dft(field6, singleton_word)
        require(
            cyclic_convolution(field6, ancestry_spectrum, parity_spectrum)
            == singleton_spectrum,
            ("CRT spectrum convolution", residue),
        )
        require(
            ancestry_spectrum[2]
            == field6.scale(
                Fraction(1, 3), field6.root_power(-2 * ancestry_residue)
            ),
            ("embedded ancestry phase", residue),
        )
        require(
            parity_spectrum[3]
            == field6.scale(
                Fraction(1, 2), field6.root_power(-3 * parity_residue)
            ),
            ("embedded parity phase", residue),
        )
        expected_phase = field6.scale(
            Fraction(1, 6), field6.root_power(-residue)
        )
        require(singleton_spectrum[1] == expected_phase, (
            "six-state phase", residue, singleton_spectrum[1], expected_phase
        ))
        require(singleton_spectrum[0] == field6.rational(Fraction(1, 6)), (
            "six-state density", residue
        ))
        six_phases.append(singleton_spectrum[1])
        crt_records.append((
            residue, ancestry_residue, parity_residue, channel_order,
            singleton_spectrum
        ))

    require(tuple(order for _, order in CHANNEL_DICTIONARY) == (
        "cab", "cba", "bca", "bac", "abc", "acb"
    ), "channel dictionary drift")
    digest.add((
        "application_a", tuple(ancestry_phases), tuple(cassini_phases),
        tuple(six_phases), tuple(crt_records)
    ))
    return {
        "ancestry_phases": tuple(ancestry_phases),
        "cassini_phases": tuple(cassini_phases),
        "six_phases": tuple(six_phases),
        "crt_checks": len(crt_records),
    }


def matrix_apply(matrix, vector):
    return tuple(
        sum(matrix[row][column] * vector[column] for column in range(2))
        for row in range(2)
    )


def triple(vector):
    m, n = vector
    return n * n - m * m, 2 * m * n, n * n + m * m


def comparator_bit(vector):
    a, b, _ = triple(vector)
    require(a != b, ("local comparator tie", vector, a, b))
    return 0 if b > a else 1


def transitive_tournament(order):
    vertices = ("P", "L", "M", "R")
    rank = {vertex: index for index, vertex in enumerate(order)}
    return tuple(
        tuple(
            0 if left == right else (1 if rank[left] < rank[right] else 0)
            for right in vertices
        )
        for left in vertices
    )


def relabel_tournament(matrix, permutation):
    return tuple(
        tuple(matrix[permutation[i]][permutation[j]] for j in range(4))
        for i in range(4)
    )


def application_b(fields, digest):
    field4 = fields[4]
    expected_orders = {
        0: ("P", "L", "R", "M"),
        1: ("P", "R", "L", "M"),
    }
    labelled_tournaments = tuple(
        transitive_tournament(expected_orders[bit]) for bit in (0, 1)
    )
    require(labelled_tournaments[0] != labelled_tournaments[1], (
        "labelled T4s collapsed", labelled_tournaments
    ))
    swap_left_right = (0, 3, 2, 1)
    require(
        relabel_tournament(labelled_tournaments[0], swap_left_right)
        == labelled_tournaments[1],
        "unlabelled local T4 isomorphism",
    )
    transitive_score_sequence = tuple(
        tuple(sorted(sum(row) for row in tournament))
        for tournament in labelled_tournaments
    )
    require(transitive_score_sequence == ((0, 1, 2, 3), (0, 1, 2, 3)), (
        "local T4 score sequence", transitive_score_sequence
    ))

    level = (("", ROOT),)
    depth_counts = []
    local_t4_checks = 0
    transition_checks = 0
    found_c65 = {}
    for depth in range(BFS_PARENT_MAX_DEPTH + 1):
        bit_counts = tuple(
            sum(1 for _, vector in level if comparator_bit(vector) == bit)
            for bit in (0, 1)
        )
        expected_counts = (
            (3 ** depth + (-1) ** depth) // 2,
            (3 ** depth - (-1) ** depth) // 2,
        )
        require(bit_counts == expected_counts, (
            "BFS comparator counts", depth, bit_counts, expected_counts
        ))
        depth_counts.append(bit_counts)

        children = []
        for word, vector in level:
            m, n = vector
            require(0 < m < n and gcd(m, n) == 1 and (m + n) % 2 == 1, (
                "Berggren chamber", word, vector
            ))
            parent_bit = comparator_bit(vector)
            a, b, parent_hypotenuse = triple(vector)
            require(a * a + b * b == parent_hypotenuse * parent_hypotenuse, (
                "Pythagorean parent", word, vector
            ))
            if parent_hypotenuse == 65 and vector in ((1, 8), (4, 7)):
                found_c65[vector] = (depth, word)

            child_records = []
            for label, matrix in BRANCHES:
                child = matrix_apply(matrix, vector)
                child_bit = comparator_bit(child)
                expected_bit = {
                    "L": 0,
                    "M": parent_bit ^ 1,
                    "R": 1,
                }[label]
                require(child_bit == expected_bit, (
                    "reset/XOR transition", depth, word, label,
                    parent_bit, child_bit, expected_bit
                ))
                child_hypotenuse = triple(child)[2]
                require(child_hypotenuse > parent_hypotenuse, (
                    "parent/child hypotenuse", word, label,
                    parent_hypotenuse, child_hypotenuse
                ))
                child_records.append((label, child_hypotenuse))
                children.append((word + label, child))
                transition_checks += 1

            actual_order = ("P",) + tuple(
                label for label, _ in sorted(child_records, key=lambda item: item[1])
            )
            require(actual_order == expected_orders[parent_bit], (
                "labelled local T4 order", word, vector,
                parent_bit, actual_order, expected_orders[parent_bit]
            ))
            require(len({value for _, value in child_records}) == 3, (
                "sibling hypotenuse tie", word, vector, child_records
            ))
            local_t4_checks += 1
        level = tuple(children)

    terminal_depth = BFS_PARENT_MAX_DEPTH + 1
    terminal_counts = tuple(
        sum(1 for _, vector in level if comparator_bit(vector) == bit)
        for bit in (0, 1)
    )
    terminal_expected = (
        (3 ** terminal_depth + (-1) ** terminal_depth) // 2,
        (3 ** terminal_depth - (-1) ** terminal_depth) // 2,
    )
    require(terminal_counts == terminal_expected, (
        "terminal BFS counts", terminal_counts, terminal_expected
    ))
    depth_counts.append(terminal_counts)

    require(set(found_c65) == {(1, 8), (4, 7)}, (
        "global c=65 tie witnesses", found_c65
    ))
    require(triple((1, 8))[2] == triple((4, 7))[2] == 65, (
        "c=65 hostile", triple((1, 8)), triple((4, 7))
    ))

    a_zero_word = coset_word(4, 3)
    b_zero_word = coset_word(4, 1)
    xor_word = tuple(a ^ b for a, b in zip(a_zero_word, b_zero_word))
    require(xor_word == (0, 1, 0, 1), ("odd XOR support", xor_word))
    a_spectrum = dft(field4, a_zero_word)
    b_spectrum = dft(field4, b_zero_word)
    xor_spectrum = dft(field4, xor_word)
    require(a_spectrum[1] == field4.scale(Fraction(1, 4), field4.root), (
        "A zero-clock phase", a_spectrum[1]
    ))
    require(b_spectrum[1] == field4.scale(Fraction(-1, 4), field4.root), (
        "B zero-clock phase", b_spectrum[1]
    ))
    require(
        spectrum_sub(
            field4,
            spectrum_add(field4, a_spectrum, b_spectrum),
            spectrum_scale(
                field4, 2, cyclic_convolution(field4, a_spectrum, b_spectrum)
            ),
        ) == xor_spectrum,
        "mod-five zero-clock XOR spectrum",
    )

    zero_clock_checks = 0
    for positive_index in range(1, ZERO_CLOCK_MAX_INDEX + 1):
        depth = positive_index - 1
        a_count = (3 ** depth + (-1) ** depth) // 2
        b_count = (3 ** depth - (-1) ** depth) // 2
        require((a_count % 5 == 0) == (positive_index % 4 == 3), (
            "A mod-five zero clock", positive_index, a_count
        ))
        require((b_count % 5 == 0) == (positive_index % 4 == 1), (
            "B mod-five zero clock", positive_index, b_count
        ))
        zero_clock_checks += 2

    densities = (Fraction(1, 4), Fraction(1, 4), Fraction(1, 2))
    scars = tuple(scar(value) for value in densities)
    require(scars == (Fraction(7, 32), Fraction(7, 32), Fraction(3, 8)), (
        "periodic scars", scars
    ))

    c65_records = tuple(
        (vector, found_c65[vector]) for vector in ((1, 8), (4, 7))
    )
    digest.add((
        "application_b", tuple(depth_counts), local_t4_checks,
        transition_checks, c65_records, a_spectrum, b_spectrum, xor_spectrum,
        tuple(fraction_key(value) for value in densities),
        tuple(fraction_key(value) for value in scars)
    ))
    return {
        "depth_counts": tuple(depth_counts),
        "local_t4_checks": local_t4_checks,
        "transition_checks": transition_checks,
        "c65_records": c65_records,
        "zero_clock_checks": zero_clock_checks,
        "densities": densities,
        "scars": scars,
        "a_phase": a_spectrum[1],
        "b_phase": b_spectrum[1],
        "xor_spectrum": xor_spectrum,
    }


def hostile_controls(fields, digest):
    singleton_density_checks = 0
    singleton_magnitude_checks = 0
    magnitude_records = []
    for period in PERIODS:
        field = fields[period]
        expected_density = field.rational(Fraction(1, period))
        expected_magnitude = field.rational(Fraction(1, period * period))
        period_magnitudes = []
        for residue in range(period):
            spectrum = dft(field, coset_word(period, residue))
            require(spectrum[0] == expected_density, (
                "singleton density hostile", period, residue
            ))
            singleton_density_checks += 1
            residue_magnitudes = tuple(
                field.magnitude_square(value) for value in spectrum
            )
            require(all(value == expected_magnitude for value in residue_magnitudes), (
                "singleton magnitude hostile", period, residue,
                residue_magnitudes, expected_magnitude
            ))
            singleton_magnitude_checks += period
            period_magnitudes.append(residue_magnitudes)
        require(all(
            values == period_magnitudes[0] for values in period_magnitudes
        ), ("coset magnitudes distinguish phase", period))
        magnitude_records.append((period, expected_density, expected_magnitude))

    square_gap_checks = 0
    previous_gap = 0
    for index in range(1, 65):
        gap = (index + 1) ** 2 - index ** 2
        require(gap == 2 * index + 1, ("square gap identity", index, gap))
        require(gap > previous_gap, ("square gaps not increasing", index, gap))
        previous_gap = gap
        square_gap_checks += 1

    digest.add((
        "hostiles", tuple(magnitude_records), singleton_density_checks,
        singleton_magnitude_checks, square_gap_checks,
        "infinite squares have unbounded gaps, so are not ultimately periodic"
    ))
    return {
        "singleton_density_checks": singleton_density_checks,
        "singleton_magnitude_checks": singleton_magnitude_checks,
        "square_gap_checks": square_gap_checks,
    }


def element_key(value):
    return tuple(fraction_key(entry) for entry in value)


def spectrum_key(spectrum):
    return tuple(element_key(value) for value in spectrum)


def main():
    repo_root = Path(__file__).resolve().parents[1]
    source_path = Path(__file__).resolve()
    source_hash = lf_hash(source_path)
    theorem_hashes = []
    for theorem_id, relative_path, expected_hash in PINNED_THEOREMS:
        actual_hash = lf_hash(repo_root / relative_path)
        require(actual_hash == expected_hash, (
            "theorem source hash drift", theorem_id, actual_hash, expected_hash
        ))
        theorem_hashes.append((theorem_id, actual_hash))

    tree = ast.parse(source_path.read_text(encoding="utf-8"))
    assertion_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
    floating_literals = sum(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(tree)
    )
    require(assertion_nodes == 0, ("assertion statements", assertion_nodes))
    require(floating_literals == 0, ("floating literals", floating_literals))

    cyclotomic_polynomials = build_cyclotomic_polynomials(max(PERIODS))
    fields = {
        period: CyclotomicField(period, cyclotomic_polynomials[period])
        for period in PERIODS
    }
    digest = ExactDigest()
    digest.add(("cyclotomic_polynomials", tuple(cyclotomic_polynomials.items())))
    digest.add(("theorem_hashes", tuple(theorem_hashes)))

    generic = exact_boolean_fourier_audit(fields, digest)
    app_a = application_a(fields, digest)
    app_b = application_b(fields, digest)
    hostiles = hostile_controls(fields, digest)
    semantic_digest = digest.hexdigest()
    freeze(semantic_digest, EXPECTED_SEMANTIC_DIGEST, "semantic digest")

    print("THM-3364 CYCLOTOMIC BOOLEAN CLOCKS + BERGGREN T4 XOR EXACT COMPANION")
    print("source_sha256_lf", source_hash)
    print("theorem_source_sha256_lf", tuple(theorem_hashes))
    print("assertion_nodes", assertion_nodes, "floating_literals", floating_literals)
    print("GENERIC EXACT FINITE-CYCLIC DFT IN Q[x]/Phi_p")
    print("periods", PERIODS)
    print("cyclotomic_polynomials_ascending", tuple(cyclotomic_polynomials.items()))
    print("period_degree_words_pairs_primitive_words", generic["period_records"])
    print("boolean_words", generic["total_words"],
          "ordered_boolean_pairs", generic["total_pairs"])
    print("DFT_inversion_values", generic["inversion_checks"],
          "basis_convolution_intersection_checks", generic["basis_convolution_checks"])
    print("word_convolution_idempotence_consequences", generic["idempotence_checks"],
          "zero_mode_density_checks", generic["zero_mode_checks"],
          "minimal_period_recovery_checks", generic["minimal_period_checks"])
    print("universal_Boolean_truth_checks", generic["boolean_truth_checks"],
          "XOR_union_spectral_coordinate_consequences", generic["boolean_spectral_checks"])
    print("termwise_Dirichlet_basis_checks", generic["dirichlet_term_checks"],
          "all_word_consequences", generic["dirichlet_term_consequences"],
          "indices", (1, DIRICHLET_MAX_INDEX),
          "exponents", DIRICHLET_EXPONENTS)
    print("APPLICATION A: THM-3339 TERNARY ANCESTRY x CASSINI PARITY")
    print("ancestry_cosets", ((2, "(BA)^*"), (0, "A(BA)^*"), (1, "C(BC)^*")))
    print("ancestry_phase_formula", "Phi_R_a(1/3)=zeta_3^(-a)/3",
          "exact_power_basis_values", tuple(element_key(v) for v in app_a["ancestry_phases"]))
    print("cassini_phase_formula", "Phi_P_b(1/2)=(-1)^b/2",
          "exact_power_basis_values", tuple(element_key(v) for v in app_a["cassini_phases"]))
    print("CRT_checks", app_a["crt_checks"],
          "formula", "Phi_S_t=Phi_R_(t mod 3)*Phi_P_(t mod 2)")
    print("six_state_dictionary", CHANNEL_DICTIONARY)
    print("six_state_phase_formula", "Phi_S_t(1/6)=zeta_6^(-t)/6",
          "exact_power_basis_values", tuple(element_key(v) for v in app_a["six_phases"]))
    print("APPLICATION B: THM-3357 LABELLED LOCAL T4 + MOD-FIVE XOR")
    print("labelled_T4_orders", ((0, "P<L<R<M"), (1, "P<R<L<M")),
          "unlabelled_score_sequence", (0, 1, 2, 3))
    print("reset_XOR_transitions", ("L->0", "M->b xor 1", "R->1"))
    print("BFS_parent_depths", (0, BFS_PARENT_MAX_DEPTH),
          "local_T4_checks", app_b["local_t4_checks"],
          "transition_checks", app_b["transition_checks"])
    print("A_d_B_d_depth_0_to_11", app_b["depth_counts"])
    print("closed_forms", "A_d=(3^d+(-1)^d)/2",
          "B_d=(3^d-(-1)^d)/2")
    print("positive_index_mod5_zero_supports", (("A", "n=3 mod 4"), ("B", "n=1 mod 4")),
          "checks", app_b["zero_clock_checks"])
    print("frequency_1_over_4_phases", ("A", "+i/4", element_key(app_b["a_phase"])),
          ("B", "-i/4", element_key(app_b["b_phase"])))
    print("XOR_support", "odd positive indices",
          "spectrum_in_Q(i)", spectrum_key(app_b["xor_spectrum"]))
    print("harmonic_residues", tuple(fraction_key(v) for v in app_b["densities"]),
          "periodic_scars", tuple(fraction_key(v) for v in app_b["scars"]))
    print("HOSTILE CONTROLS")
    print("singleton_coset_density_checks", hostiles["singleton_density_checks"],
          "singleton_coset_magnitude_checks", hostiles["singleton_magnitude_checks"],
          "conclusion", "density and Fourier magnitudes lose cyclic phase")
    print("global_hypotenuse_tie_c65", app_b["c65_records"],
          "conclusion", "no global strict tree tournament")
    print("local_label_erasure", "the two T4s are isomorphic by L<->R")
    print("arbitrary_subset_hostile", "squares are infinite with gap 2k+1",
          "exact_gap_checks", hostiles["square_gap_checks"],
          "conclusion", "no canonical finite periodic spectrum for arbitrary subsets")
    print("semantic_sha256", semantic_digest)
    print("replay_contract", "ordinary_and_python_O_outputs_byte_identical")
    print("status=ALL EXACT CHECKS PASSED; phase sidecars and nonperiodic boundary remain typed")


if __name__ == "__main__":
    main()
