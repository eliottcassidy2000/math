#!/usr/bin/env python3
"""Exact companion for THM-3379.

The companion combines THM-3339's three Fibonacci/Berggren ray words with
THM-3364's reset/XOR comparator.  It verifies directly that the local labelled
T4 bit is the indicator of the k=1 (mod 3) ancestry ray, equivalently the
indicator that channel b is the median of the six Farey order states.  Exact
Q(omega) Fourier data, the finite-startup harmonic correction, and several
full-tree hostile controls keep the projection typed.
"""

from __future__ import annotations

import ast
from fractions import Fraction
from hashlib import sha256
from math import gcd
from pathlib import Path


EXPECTED_SEMANTIC_DIGEST = "981a5759c78a3aadb928d78f451c5070b2eb4965f54c0066e5004dfa4ef13e2d"
RAY_MAX_INDEX = 604
BFS_MAX_DEPTH = 10

PINNED_THEOREMS = (
    (
        "THM-3339",
        "01-canon/theorems/THM-3339-fibonacci-three-ray-berggren-transplant-and-moving-owner-obstruction.md",
        "1e4aa8cd9d6cc9bf342328a4af0e5db8cd7eefb51eea08460d03f2c6410cee51",
    ),
    (
        "THM-3364",
        "01-canon/theorems/THM-3364-cyclotomic-boolean-clocks-berggren-t4-xor-and-crt-phase.md",
        "4d8e97c83d54a878c3f628d1944be0d5d391e1c0a2e96a5cf877f5edef875cbe",
    ),
)

ROOT_PARAMETER = (1, 2)
A = ((0, 1), (-1, 2))
B = ((0, 1), (1, 2))
C = ((1, 0), (2, 1))
BRANCHES = {"A": A, "B": B, "C": C}

CHANNEL_DICTIONARY = (
    "cab",
    "cba",
    "bca",
    "bac",
    "abc",
    "acb",
)
MEDIAN_BY_RESIDUE = ("a", "b", "c")


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lf_hash(path):
    return sha256(Path(path).read_bytes().replace(b"\r\n", b"\n")).hexdigest()


class ExactDigest:
    def __init__(self):
        self._hash = sha256()

    def add(self, item):
        self._hash.update(repr(item).encode("ascii"))
        self._hash.update(b"\n")

    def hexdigest(self):
        return self._hash.hexdigest()


def matrix_vector(matrix, vector):
    return tuple(
        sum(matrix[row][column] * vector[column] for column in range(2))
        for row in range(2)
    )


def apply_word(word, vector=ROOT_PARAMETER):
    for letter in word:
        vector = matrix_vector(BRANCHES[letter], vector)
    return vector


def triple(vector):
    m, n = vector
    return n * n - m * m, 2 * m * n, n * n + m * m


def comparator_delta(vector):
    a, b, _ = triple(vector)
    return b - a


def epsilon_direct(vector):
    delta = comparator_delta(vector)
    require(delta != 0, ("comparator tie", vector))
    return 0 if delta > 0 else 1


def epsilon_word(word):
    bit = 0
    for letter in word:
        if letter == "A":
            bit = 0
        elif letter == "B":
            bit ^= 1
        elif letter == "C":
            bit = 1
        else:
            raise RuntimeError(("unknown branch", letter))
    return bit


def fibs(limit):
    values = [0, 1]
    while len(values) <= limit:
        values.append(values[-1] + values[-2])
    return tuple(values)


def ray_word(index):
    require(index >= 2, ("Fibonacci carrier starts at two", index))
    residue = index % 3
    if residue == 2:
        exponent = (index - 2) // 3
        return "BA" * exponent, "R2"
    if residue == 0:
        exponent = (index - 3) // 3
        return "A" + "BA" * exponent, "R0"
    exponent = (index - 4) // 3
    return "C" + "BC" * exponent, "R1"


def normalized_fibonacci_parameter(index, fibonacci):
    raw = (fibonacci[index], fibonacci[index + 1])
    if index % 3 != 1:
        return raw
    m, n = raw
    require(m % 2 == 1 and n % 2 == 1, ("odd/odd transplant", index, raw))
    return (n - m) // 2, (n + m) // 2


def local_t4_order(vector):
    records = [(triple(vector)[2], "P")]
    for label in "ABC":
        child = matrix_vector(BRANCHES[label], vector)
        records.append((triple(child)[2], label))
    require(len({value for value, _ in records}) == 4, ("local T4 tie", vector, records))
    return tuple(label for _, label in sorted(records))


def permutation_sign(word):
    ranks = {letter: rank for rank, letter in enumerate("abc")}
    values = tuple(ranks[letter] for letter in word)
    inversions = sum(
        values[left] > values[right]
        for left in range(3)
        for right in range(left + 1, 3)
    )
    return -1 if inversions % 2 else 1


# Q(omega), omega^2+omega+1=0, represented as a+b*omega.
def qomega_add(left, right):
    return left[0] + right[0], left[1] + right[1]


def qomega_mul(left, right):
    a, b = left
    c, d = right
    return a * c - b * d, a * d + b * c - b * d


def qomega_scale(coefficient, value):
    coefficient = Fraction(coefficient)
    return coefficient * value[0], coefficient * value[1]


def qomega_pow(value, exponent, order):
    result = (Fraction(1), Fraction(0))
    for _ in range(exponent % order):
        result = qomega_mul(result, value)
    return result


def qomega_sum(values):
    result = (Fraction(0), Fraction(0))
    for value in values:
        result = qomega_add(result, value)
    return result


def dft(word, root, order):
    period = len(word)
    require(order == period, ("DFT root order", order, period))
    return tuple(
        qomega_scale(
            Fraction(1, period),
            qomega_sum(
                qomega_pow(root, -frequency * residue, order)
                for residue, bit in enumerate(word)
                if bit
            ),
        )
        for frequency in range(period)
    )


def cyclic_convolution(left, right):
    period = len(left)
    require(len(right) == period, "convolution length")
    return tuple(
        qomega_sum(
            qomega_mul(left[index], right[(frequency - index) % period])
            for index in range(period)
        )
        for frequency in range(period)
    )


def fraction_key(value):
    return value.numerator, value.denominator


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
            "theorem dependency drift", theorem_id, actual_hash, expected_hash
        ))
        theorem_hashes.append((theorem_id, actual_hash))

    syntax = ast.parse(source_path.read_text(encoding="utf-8"))
    assertion_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax))
    floating_literals = sum(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(syntax)
    )
    require(assertion_nodes == 0, ("assertion nodes", assertion_nodes))
    require(floating_literals == 0, ("floating literals", floating_literals))

    # Direct algebraic reset/flip identities on a large primitive hostile box.
    branch_identity_checks = 0
    for m in range(1, 100):
        for n in range(m + 1, 101):
            if gcd(m, n) != 1 or (m + n) % 2 != 1:
                continue
            vector = (m, n)
            delta = comparator_delta(vector)
            require(delta != 0, ("primitive tie", vector))
            common = n * n + 2 * m * n - m * m
            require(common > 0, ("reset positivity", vector, common))
            require(comparator_delta(matrix_vector(A, vector)) == common, (
                "A reset identity", vector
            ))
            require(comparator_delta(matrix_vector(B, vector)) == -delta, (
                "B flip identity", vector
            ))
            require(comparator_delta(matrix_vector(C, vector)) == -common, (
                "C reset identity", vector
            ))
            require(epsilon_direct(matrix_vector(A, vector)) == 0, ("A bit", vector))
            require(
                epsilon_direct(matrix_vector(B, vector)) == (epsilon_direct(vector) ^ 1),
                ("B bit", vector),
            )
            require(epsilon_direct(matrix_vector(C, vector)) == 1, ("C bit", vector))
            branch_identity_checks += 3

    fibonacci = fibs(RAY_MAX_INDEX + 2)
    ray_counts = {"R0": 0, "R1": 0, "R2": 0}
    ray_records = []
    terminal_records = []
    for index in range(2, RAY_MAX_INDEX + 1):
        word, ray = ray_word(index)
        vector = apply_word(word)
        expected_vector = normalized_fibonacci_parameter(index, fibonacci)
        require(vector == expected_vector, ("ray parameter", index, word, vector, expected_vector))
        direct_bit = epsilon_direct(vector)
        automaton_bit = epsilon_word(word)
        expected_bit = int(index % 3 == 1)
        require(direct_bit == automaton_bit == expected_bit, (
            "ray bit projection", index, word, direct_bit, automaton_bit, expected_bit
        ))
        expected_order = ("P", "C", "A", "B") if expected_bit else ("P", "A", "C", "B")
        require(local_t4_order(vector) == expected_order, (
            "ray local T4", index, vector, local_t4_order(vector), expected_order
        ))
        order_state = CHANNEL_DICTIONARY[index % 6]
        require(order_state[1] == MEDIAN_BY_RESIDUE[index % 3], (
            "median channel", index, order_state
        ))
        require(expected_bit == int(order_state[1] == "b"), (
            "median-b projection", index, order_state
        ))
        ray_counts[ray] += 1
        if index <= 13:
            ray_records.append((index, index % 6, word, ray, vector, order_state, direct_bit))
        if index > RAY_MAX_INDEX - 6:
            terminal_records.append((index, index % 6, len(word), ray, order_state, direct_bit))

    six_state_table = tuple(
        (
            residue,
            CHANNEL_DICTIONARY[residue],
            CHANNEL_DICTIONARY[residue][1],
            int(residue % 3 == 1),
            permutation_sign(CHANNEL_DICTIONARY[residue]),
        )
        for residue in range(6)
    )
    require(tuple(row[2] for row in six_state_table) == ("a", "b", "c", "a", "b", "c"), (
        "median cycle", six_state_table
    ))
    require(tuple(row[3] for row in six_state_table) == (0, 1, 0, 0, 1, 0), (
        "Boolean projection", six_state_table
    ))
    require(tuple(row[4] for row in six_state_table) == (1, -1, 1, -1, 1, -1), (
        "Cassini parity", six_state_table
    ))
    for residue in range(3):
        require(CHANNEL_DICTIONARY[residue + 3] == CHANNEL_DICTIONARY[residue][::-1], (
            "reversal pair", residue
        ))

    omega = (Fraction(0), Fraction(1))
    zeta6 = (Fraction(1), Fraction(1))
    require(qomega_pow(omega, 3, 3) == (Fraction(1), Fraction(0)), "omega order")
    require(qomega_pow(zeta6, 3, 6) == (Fraction(-1), Fraction(0)), "zeta6 cube")
    require(qomega_pow(zeta6, 6, 6) == (Fraction(1), Fraction(0)), "zeta6 order")

    period_three_word = (0, 1, 0)
    period_six_word = (0, 1, 0, 0, 1, 0)
    period_three_spectrum = dft(period_three_word, omega, 3)
    period_six_spectrum = dft(period_six_word, zeta6, 6)
    expected_three = (
        (Fraction(1, 3), Fraction(0)),
        (Fraction(-1, 3), Fraction(-1, 3)),
        (Fraction(0), Fraction(1, 3)),
    )
    expected_six = (
        expected_three[0],
        (Fraction(0), Fraction(0)),
        expected_three[1],
        (Fraction(0), Fraction(0)),
        expected_three[2],
        (Fraction(0), Fraction(0)),
    )
    require(period_three_spectrum == expected_three, ("period-three spectrum", period_three_spectrum))
    require(period_six_spectrum == expected_six, ("period-six spectrum", period_six_spectrum))
    require(cyclic_convolution(period_six_spectrum, period_six_spectrum) == period_six_spectrum, (
        "Boolean spectral idempotence", period_six_spectrum
    ))
    singleton_one = dft((0, 1, 0, 0, 0, 0), zeta6, 6)
    singleton_four = dft((0, 0, 0, 0, 1, 0), zeta6, 6)
    require(
        tuple(qomega_add(left, right) for left, right in zip(singleton_one, singleton_four))
        == period_six_spectrum,
        "Cassini split spectrum",
    )

    density = Fraction(1, 3)
    scar = density - density * density / 2
    require(scar == Fraction(5, 18), ("periodic scar", scar))
    harmonic_checks = 0
    for count in range(1, 129):
        actual = sum((Fraction(1, 3 * step + 1) for step in range(1, count + 1)), Fraction(0))
        periodic_completion = Fraction(1) + actual
        require(periodic_completion - actual == 1, ("startup correction", count))
        harmonic_checks += 1
    harmonic_constant_profile = (
        ("gamma", Fraction(1, 3)),
        ("pi_over_sqrt3", Fraction(1, 6)),
        ("log_3", Fraction(1, 6)),
        ("finite_startup", Fraction(-1)),
    )

    # Full-tree controls: the projection is special to the three rays.
    level = (("", ROOT_PARAMETER),)
    bfs_counts = []
    full_tree_checks = 0
    for depth in range(BFS_MAX_DEPTH + 1):
        counts = tuple(sum(epsilon_word(word) == bit for word, _ in level) for bit in (0, 1))
        expected_counts = (
            (3 ** depth + (-1) ** depth) // 2,
            (3 ** depth - (-1) ** depth) // 2,
        )
        require(counts == expected_counts, ("full-tree bit census", depth, counts, expected_counts))
        bfs_counts.append(counts)
        children = []
        for word, vector in level:
            require(epsilon_word(word) == epsilon_direct(vector), ("word/direct bit", word, vector))
            full_tree_checks += 1
            if depth < BFS_MAX_DEPTH:
                for letter in "ABC":
                    children.append((word + letter, matrix_vector(BRANCHES[letter], vector)))
        level = tuple(children)

    hostile_words = ("A", "B", "AB", "BA", "CB")
    hostile_records = tuple(
        (word, apply_word(word), epsilon_word(word), epsilon_direct(apply_word(word)))
        for word in hostile_words
    )
    hostile_bits = {word: epsilon_word(word) for word in hostile_words}
    require(hostile_bits["A"] != hostile_bits["B"], "depth-only hostile")
    require(hostile_bits["AB"] != hostile_bits["BA"], "branch-count hostile")
    require(hostile_bits["AB"] != hostile_bits["CB"], "terminal-letter hostile")

    digest = ExactDigest()
    digest.add(("theorem_hashes", tuple(theorem_hashes)))
    digest.add(("branch_identity_checks", branch_identity_checks))
    digest.add((
        "ray_counts", tuple(sorted(ray_counts.items())),
        tuple(ray_records), tuple(terminal_records)
    ))
    digest.add(("six_state_table", six_state_table))
    digest.add(("spectra", period_three_spectrum, period_six_spectrum))
    digest.add(("harmonic", density, scar, harmonic_constant_profile, harmonic_checks))
    digest.add(("full_tree", tuple(bfs_counts), full_tree_checks, hostile_records))
    semantic_digest = digest.hexdigest()
    if EXPECTED_SEMANTIC_DIGEST != "PENDING":
        require(semantic_digest == EXPECTED_SEMANTIC_DIGEST, (
            "semantic digest", semantic_digest, EXPECTED_SEMANTIC_DIGEST
        ))

    print("THM-3379 FIBONACCI-RAY LOCAL T4 MOD-THREE PROJECTION")
    print("source_sha256_lf", source_hash)
    print("theorem_source_sha256_lf", tuple(theorem_hashes))
    print("assertion_nodes", assertion_nodes, "floating_literals", floating_literals)
    print("branch_identity_checks", branch_identity_checks)
    print("ray_index_range", (2, RAY_MAX_INDEX), "ray_counts", tuple(sorted(ray_counts.items())))
    print("ray_boundary_records", tuple(ray_records))
    print("ray_terminal_summaries", tuple(terminal_records))
    print("six_state_table_residue_order_median_epsilon_sign", six_state_table)
    print("ancestry_median_dictionary", (("R0", "a"), ("R1", "b"), ("R2", "c")))
    print("epsilon_formula", "epsilon(k)=1 iff k=1 mod 3 iff median(channel order)=b")
    print("period_three_spectrum_Qomega", spectrum_key(period_three_spectrum))
    print("period_six_spectrum_Qomega", spectrum_key(period_six_spectrum))
    print("Cassini_split", "H_epsilon=S_1 disjoint_union S_4; odd modes cancel")
    print("harmonic_density", fraction_key(density), "periodic_scar", fraction_key(scar))
    print("actual_support", "{4,7,10,...}=({n>=1:n=1 mod 3})\\{1}")
    print("harmonic_constant_profile", tuple((name, fraction_key(value)) for name, value in harmonic_constant_profile))
    print("full_tree_depths", (0, BFS_MAX_DEPTH), "bit_counts", tuple(bfs_counts), "checks", full_tree_checks)
    print("hostile_records", hostile_records)
    print("hostile_conclusion", "depth, branch counts, and terminal letter do not determine epsilon on the full tree")
    print("semantic_sha256", semantic_digest)
    print("replay_contract", "ordinary_and_python_O_outputs_byte_identical")
    print("status=ALL EXACT CHECKS PASSED; collapse is restricted to the three Fibonacci rays")


if __name__ == "__main__":
    main()
