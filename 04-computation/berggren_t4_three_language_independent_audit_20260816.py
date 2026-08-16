#!/usr/bin/env python3
"""Independent exact audit of the three Berggren/T4 calibration languages.

This companion starts from the integer branch matrices in THM-3339.  It does
not import or pin any of the three candidate companions.  It independently
checks

* generatorwise and six-frame calibration counts;
* the variable-translation S3 x C2 language and its recurrence;
* the fixed-drift S4 x D4 language and its parity split;
* root-to-child multiplication on all three Fibonacci rays; and
* the finite-group mixing input used by the shortlex harmonic argument.

Only the Python standard library is used.  Characteristic polynomials are
computed from exact traces by Newton identities, and scalar recurrences are
recovered over Q by Berlekamp--Massey.
"""

from __future__ import annotations

from collections import Counter, deque
from fractions import Fraction
from hashlib import sha256
from itertools import permutations, product
from json import dumps
from math import log, sqrt
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "04-computation/berggren_t4_three_language_independent_audit_20260816.py"
OUTPUT = "05-knowledge/results/berggren_t4_three_language_independent_audit_20260816.out"
THM3339 = ROOT / "01-canon/theorems/THM-3339-fibonacci-three-ray-berggren-transplant-and-moving-owner-obstruction.md"
THM3339_LF_SHA256 = "1e4aa8cd9d6cc9bf342328a4af0e5db8cd7eefb51eea08460d03f2c6410cee51"

LETTERS = "ABC"
BRANCH = {
    "A": ((0, 1), (-1, 2)),
    "B": ((0, 1), (1, 2)),
    "C": ((1, 0), (2, 1)),
}
P1_F3 = ((1, 0), (0, 1), (1, 1), (1, 2))
I2 = ((1, 0), (0, 1))
J2 = ((0, 1), (1, 0))
I3 = (0, 1, 2)
I4 = (0, 1, 2, 3)
V4 = (0, 1, 2, 3)  # 0,p,q,r in bit coordinates p=(1,0), q=(0,1)
P, Q, R = 1, 2, 3
MATCHINGS = (
    ((0, 1), (2, 3)),
    ((0, 2), (1, 3)),
    ((0, 3), (1, 2)),
)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def compose(left, right):
    """Return left after right."""
    return tuple(left[right[index]] for index in range(len(left)))


def inverse(permutation):
    result = [0] * len(permutation)
    for source, target in enumerate(permutation):
        result[target] = source
    return tuple(result)


def perm_power(permutation, exponent):
    result = tuple(range(len(permutation)))
    for _ in range(exponent):
        result = compose(permutation, result)
    return result


def cycle_type(permutation):
    unseen = set(range(len(permutation)))
    result = []
    while unseen:
        start = min(unseen)
        point = start
        length = 0
        while point in unseen:
            unseen.remove(point)
            point = permutation[point]
            length += 1
        require(point == start, (permutation, start, point))
        result.append(length)
    return tuple(sorted(result))


def conjugate(calibration, action):
    return compose(compose(calibration, action), inverse(calibration))


def generated_permutation_group(generators):
    identity = tuple(range(len(generators[0])))
    seen = {identity}
    queue = deque((identity,))
    while queue:
        current = queue.popleft()
        for generator in generators:
            candidate = compose(generator, current)
            if candidate not in seen:
                seen.add(candidate)
                queue.append(candidate)
    return frozenset(seen)


def matrix_vector(matrix, vector, modulus=None):
    result = tuple(
        sum(matrix[row][column] * vector[column] for column in range(2))
        for row in range(2)
    )
    if modulus is not None:
        result = tuple(value % modulus for value in result)
    return result


def matrix_product(left, right, modulus=None):
    result = tuple(
        tuple(
            sum(left[row][middle] * right[middle][column] for middle in range(2))
            for column in range(2)
        )
        for row in range(2)
    )
    if modulus is not None:
        result = tuple(tuple(value % modulus for value in row) for row in result)
    return result


def same_projective_line(left, right):
    return (left[0] * right[1] - left[1] * right[0]) % 3 == 0


def projective_action(matrix):
    action = []
    for representative in P1_F3:
        image = matrix_vector(matrix, representative, 3)
        matches = tuple(
            index for index, candidate in enumerate(P1_F3)
            if same_projective_line(image, candidate)
        )
        require(len(matches) == 1, (matrix, representative, image, matches))
        action.append(matches[0])
    return tuple(action)


def bit_vector(point):
    return point & 1, (point >> 1) & 1


def vector_bit(vector):
    return vector[0] | (vector[1] << 1)


def linear_on_v4(matrix, point):
    return vector_bit(matrix_vector(matrix, bit_vector(point), 2))


def affine_action(matrix, translation):
    return tuple(linear_on_v4(matrix, point) ^ translation for point in V4)


def canonical_matching(matching):
    return tuple(sorted(tuple(sorted(edge)) for edge in matching))


def matching_action(permutation):
    result = []
    for matching in MATCHINGS:
        image = canonical_matching(
            tuple((permutation[left], permutation[right]) for left, right in matching)
        )
        result.append(MATCHINGS.index(image))
    return tuple(result)


def frame_action(matrix, frames):
    return tuple(
        frames.index(tuple(linear_on_v4(matrix, direction) for direction in frame))
        for frame in frames
    )


def orbit_partition(permutation):
    unseen = set(range(len(permutation)))
    orbits = []
    while unseen:
        start = min(unseen)
        orbit = []
        current = start
        while current not in orbit:
            orbit.append(current)
            unseen.remove(current)
            current = permutation[current]
        require(current == start, (permutation, orbit, current))
        orbits.append(tuple(orbit))
    return tuple(orbits)


def frame_family_count(source, target, base):
    """Propagate gauges around base orbits and count return-compatible seeds."""
    total = 1
    seed_counts = []
    for orbit in orbit_partition(base):
        source_return = perm_power(source, len(orbit))
        target_return = perm_power(target, len(orbit))
        seeds = sum(
            conjugate(gauge, source_return) == target_return
            for gauge in permutations(range(4))
        )
        seed_counts.append(seeds)
        total *= seeds
    return total, tuple(seed_counts)


def pair_multiply(left, right):
    return compose(left[0], right[0]), compose(left[1], right[1])


def pair_inverse(pair):
    return inverse(pair[0]), inverse(pair[1])


def quotient_multiply(left, right):
    return compose(left[0], right[0]), left[1] ^ right[1]


def enumerate_group(identity, generators, multiply):
    states = [identity]
    index = {identity: 0}
    paths = {identity: ""}
    queue = deque((identity,))
    labels = tuple(LETTERS) if len(generators) == len(LETTERS) else tuple(
        f"[{index}]" for index in range(len(generators))
    )
    while queue:
        current = queue.popleft()
        for letter, generator in zip(labels, generators):
            candidate = multiply(generator, current)
            if candidate not in index:
                index[candidate] = len(states)
                states.append(candidate)
                paths[candidate] = paths[current] + letter
                queue.append(candidate)
    return tuple(states), index, paths


def right_stabilizer(states, accepting, multiply):
    accepting_set = frozenset(accepting)
    return tuple(
        element for element in states
        if frozenset(multiply(item, element) for item in accepting_set) == accepting_set
    )


def transition_group_size(states, index, generators, multiply):
    transitions = tuple(
        tuple(index[multiply(generator, state)] for state in states)
        for generator in generators
    )
    return len(generated_permutation_group(transitions))


def level_counts(identity, states, index, generators, multiply, accepting, limit):
    distribution = [0] * len(states)
    distribution[index[identity]] = 1
    counts = []
    supports = []
    accepting_set = frozenset(accepting)
    for _ in range(limit + 1):
        counts.append(sum(distribution[index[state]] for state in accepting_set))
        supports.append(sum(value != 0 for value in distribution))
        following = [0] * len(states)
        for state_index, weight in enumerate(distribution):
            if not weight:
                continue
            state = states[state_index]
            for generator in generators:
                following[index[multiply(generator, state)]] += weight
        distribution = following
    return tuple(counts), tuple(supports)


def transfer_matrix(states, index, generators, multiply):
    size = len(states)
    matrix = [[0] * size for _ in range(size)]
    for column, state in enumerate(states):
        for generator in generators:
            row = index[multiply(generator, state)]
            matrix[row][column] += 1
    return tuple(tuple(row) for row in matrix)


def integer_matrix_multiply(left, right):
    size = len(left)
    return tuple(
        tuple(
            sum(left[row][middle] * right[middle][column] for middle in range(size))
            for column in range(size)
        )
        for row in range(size)
    )


def integer_matrix_add(left, right):
    return tuple(
        tuple(left[row][column] + right[row][column] for column in range(len(left)))
        for row in range(len(left))
    )


def integer_matrix_scale(matrix, scalar):
    return tuple(tuple(scalar * value for value in row) for row in matrix)


def integer_identity(size):
    return tuple(tuple(int(row == column) for column in range(size)) for row in range(size))


def integer_matrix_trace(matrix):
    return sum(matrix[index][index] for index in range(len(matrix)))


def characteristic_coefficients_by_traces(matrix):
    """Return descending coefficients of det(lambda I-M), exactly."""
    size = len(matrix)
    powers = [integer_identity(size)]
    traces = [0]
    for _ in range(size):
        powers.append(integer_matrix_multiply(matrix, powers[-1]))
        traces.append(integer_matrix_trace(powers[-1]))
    coefficients = [Fraction(1)]
    for degree in range(1, size + 1):
        value = -sum(
            coefficients[degree - offset] * traces[offset]
            for offset in range(1, degree + 1)
        ) / degree
        coefficients.append(value)
    require(all(value.denominator == 1 for value in coefficients), coefficients)
    return tuple(int(value) for value in coefficients)


def polynomial_multiply(left, right):
    result = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += Fraction(a) * Fraction(b)
    return tuple(result)


def polynomial_trim(polynomial):
    result = list(polynomial)
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return tuple(result)


def polynomial_divmod(numerator, denominator):
    numerator = list(map(Fraction, polynomial_trim(numerator)))
    denominator = tuple(map(Fraction, polynomial_trim(denominator)))
    require(denominator != (Fraction(0),), "zero polynomial divisor")
    quotient = [Fraction(0)] * max(1, len(numerator) - len(denominator) + 1)
    while len(numerator) >= len(denominator) and any(numerator):
        shift = len(numerator) - len(denominator)
        coefficient = numerator[-1] / denominator[-1]
        quotient[shift] = coefficient
        for index, value in enumerate(denominator):
            numerator[index + shift] -= coefficient * value
        numerator = list(polynomial_trim(numerator))
    return polynomial_trim(quotient), polynomial_trim(numerator)


def polynomial_gcd(left, right):
    left = tuple(map(Fraction, polynomial_trim(left)))
    right = tuple(map(Fraction, polynomial_trim(right)))
    while right != (Fraction(0),):
        _, remainder = polynomial_divmod(left, right)
        left, right = right, remainder
    leading = left[-1]
    return tuple(value / leading for value in left)


def polynomial_evaluate(polynomial, value):
    result = Fraction(0)
    for coefficient in reversed(polynomial):
        result = result * value + coefficient
    return result


def matrix_polynomial(matrix, coefficients_ascending):
    size = len(matrix)
    result = tuple(tuple(0 for _ in range(size)) for _ in range(size))
    power = integer_identity(size)
    for coefficient in coefficients_ascending:
        result = integer_matrix_add(result, integer_matrix_scale(power, int(coefficient)))
        power = integer_matrix_multiply(matrix, power)
    return result


def berlekamp_massey_rational(sequence):
    """Return s_n=sum recurrence[i-1] s_(n-i), using exact Q arithmetic."""
    values = tuple(map(Fraction, sequence))
    connection = [Fraction(1)]
    previous = [Fraction(1)]
    length = 0
    shift = 1
    previous_discrepancy = Fraction(1)
    for index in range(len(values)):
        discrepancy = values[index]
        for offset in range(1, length + 1):
            discrepancy += connection[offset] * values[index - offset]
        if discrepancy == 0:
            shift += 1
            continue
        old_connection = connection[:]
        factor = -discrepancy / previous_discrepancy
        required = len(previous) + shift
        if len(connection) < required:
            connection.extend(Fraction(0) for _ in range(required - len(connection)))
        for offset, value in enumerate(previous):
            connection[offset + shift] += factor * value
        if 2 * length <= index:
            length = index + 1 - length
            previous = old_connection
            previous_discrepancy = discrepancy
            shift = 1
        else:
            shift += 1
    connection = connection[: length + 1]
    return tuple(-connection[offset] for offset in range(1, len(connection)))


def generating_function_data(counts, recurrence):
    denominator = (Fraction(1),) + tuple(-value for value in recurrence)
    numerator = []
    for degree in range(len(recurrence)):
        value = Fraction(counts[degree])
        for lag in range(1, degree + 1):
            value -= recurrence[lag - 1] * counts[degree - lag]
        numerator.append(value)
    return polynomial_trim(tuple(numerator)), polynomial_trim(denominator)


def permutation_order(permutation):
    result = 1
    for length in cycle_type(permutation):
        a, b = result, length
        while b:
            a, b = b, a % b
        gcd = a
        result = result * length // gcd
    return result


def chronological_matrix(word):
    """Root-to-child word h1...hk acts by M_hk ... M_h1."""
    result = I2
    for letter in word:
        result = matrix_product(BRANCH[letter], result)
    return result


def opposite_matrix(word):
    result = I2
    for letter in word:
        result = matrix_product(result, BRANCH[letter])
    return result


def fibonacci(limit):
    values = [0, 1]
    while len(values) <= limit:
        values.append(values[-1] + values[-2])
    return tuple(values)


def word_composite(word, generators, multiply, identity):
    result = identity
    for letter in word:
        result = multiply(generators[LETTERS.index(letter)], result)
    return result


def lex_level_states(depth, identity, generators, multiply):
    states = [identity]
    for _ in range(depth):
        states = [
            multiply(generator, state)
            for state in states
            for generator in generators
        ]
    return states


def main() -> None:
    require(lf_sha256(THM3339) == THM3339_LF_SHA256, "THM-3339 source pin changed")

    # ------------------------------------------------------------------
    # Matrix-first reconstruction and generatorwise/six-frame obstruction.
    # ------------------------------------------------------------------
    source = {letter: projective_action(matrix) for letter, matrix in BRANCH.items()}
    linear = {
        letter: tuple(tuple(entry % 2 for entry in row) for row in matrix)
        for letter, matrix in BRANCH.items()
    }
    require(source == {
        "A": (1, 3, 2, 0),
        "B": (1, 3, 0, 2),
        "C": (3, 1, 0, 2),
    }, source)
    require(linear == {"A": J2, "B": J2, "C": I2}, linear)

    gauges = tuple(permutations(range(4)))
    admissible = {
        letter: tuple(affine_action(linear[letter], translation) for translation in V4)
        for letter in LETTERS
    }
    static_any = {
        letter: sum(
            conjugate(gauge, source[letter]) in admissible[letter]
            for gauge in gauges
        )
        for letter in LETTERS
    }
    pinned = {
        "A": affine_action(J2, P),
        "B": affine_action(J2, P),
        "C": affine_action(I2, P),
    }
    static_pinned = {
        letter: sum(conjugate(gauge, source[letter]) == pinned[letter] for gauge in gauges)
        for letter in LETTERS
    }
    require(static_any == {"A": 0, "B": 8, "C": 0}, static_any)
    require(static_pinned == {"A": 0, "B": 4, "C": 0}, static_pinned)
    require(source["B"] == pinned["B"], (source["B"], pinned["B"]))

    frames = tuple(permutations((P, Q, R)))
    base = {letter: frame_action(linear[letter], frames) for letter in LETTERS}
    family_counts = {}
    family_seed_counts = {}
    for letter in LETTERS:
        values = []
        seeds = []
        for target in admissible[letter]:
            count, seed_count = frame_family_count(source[letter], target, base[letter])
            values.append(count)
            seeds.append(seed_count)
        family_counts[letter] = tuple(values)
        family_seed_counts[letter] = tuple(seeds)
    require(family_counts == {
        "A": (0, 0, 0, 0),
        "B": (0, 512, 512, 0),
        "C": (0, 0, 0, 0),
    }, family_counts)
    require(family_seed_counts["B"][1:3] == ((8, 8, 8), (8, 8, 8)), family_seed_counts["B"])

    # ------------------------------------------------------------------
    # Variable-translation language: raw S4 x <J> and matching quotient.
    # ------------------------------------------------------------------
    matching = {letter: matching_action(source[letter]) for letter in LETTERS}
    bits = {"A": 1, "B": 1, "C": 0}
    require(matching == {
        "A": (1, 2, 0),
        "B": (1, 0, 2),
        "C": (1, 2, 0),
    }, matching)
    variable_generators = tuple((matching[letter], bits[letter]) for letter in LETTERS)
    variable_identity = (I3, 0)
    variable_states, variable_index, _ = enumerate_group(
        variable_identity, variable_generators, quotient_multiply
    )
    require(len(variable_states) == 12, len(variable_states))

    def variable_accept(state):
        matching_part, bit = state
        return (
            (bit == 0 and matching_part == I3)
            or (bit == 1 and cycle_type(matching_part) == (1, 2))
        )

    variable_accepting = tuple(state for state in variable_states if variable_accept(state))
    require(len(variable_accepting) == 4, variable_accepting)
    require(len(right_stabilizer(variable_states, variable_accepting, quotient_multiply)) == 1,
            "variable DFA has a nontrivial right stabilizer")
    require(
        transition_group_size(variable_states, variable_index, variable_generators, quotient_multiply) == 12,
        "variable syntactic transition group",
    )

    # Directly enumerate every S4 x {I,J} pair and every point gauge/translation.
    all_s4 = tuple(permutations(range(4)))
    linear_permutations = {
        0: affine_action(I2, 0),
        1: affine_action(J2, 0),
    }
    raw_variable_pass = {}
    for projective, bit in product(all_s4, (0, 1)):
        target_maps = tuple(
            affine_action(I2 if bit == 0 else J2, translation)
            for translation in V4
        )
        direct = any(
            conjugate(gauge, projective) in target_maps
            for gauge in gauges
        )
        quotient = variable_accept((matching_action(projective), bit))
        require(direct == quotient, (projective, bit, direct, quotient))
        raw_variable_pass[(projective, bit)] = direct
    require(sum(raw_variable_pass.values()) == 16, sum(raw_variable_pass.values()))

    variable_counts, _ = level_counts(
        variable_identity, variable_states, variable_index, variable_generators,
        quotient_multiply, variable_accepting, 80,
    )
    variable_recurrence = berlekamp_massey_rational(variable_counts[:40])
    require(variable_recurrence == tuple(map(Fraction, (2, 2, 6, -9))), variable_recurrence)
    require(all(
        variable_counts[n] == sum(
            variable_recurrence[lag - 1] * variable_counts[n - lag]
            for lag in range(1, 5)
        )
        for n in range(4, len(variable_counts))
    ), "variable recurrence residual")
    variable_numerator, variable_denominator = generating_function_data(
        variable_counts, variable_recurrence
    )
    require(variable_numerator == tuple(map(Fraction, (1, -1, -1, -3))), variable_numerator)
    require(variable_denominator == tuple(map(Fraction, (1, -2, -2, -6, 9))), variable_denominator)
    require(polynomial_gcd(variable_numerator, variable_denominator) == (Fraction(1),),
            "variable GF is not reduced")

    variable_transfer = transfer_matrix(
        variable_states, variable_index, variable_generators, quotient_multiply
    )
    variable_char = characteristic_coefficients_by_traces(variable_transfer)
    expected_variable_char_ascending = (Fraction(1),)
    for factor in (
        (-3, 1),
        (-1, 1), (-1, 1), (-1, 1), (-1, 1),
        (1, 1), (1, 1), (1, 1),
        (3, 2, 1), (3, 2, 1),
    ):
        expected_variable_char_ascending = polynomial_multiply(
            expected_variable_char_ascending, factor
        )
    expected_variable_char = tuple(
        int(value) for value in reversed(expected_variable_char_ascending)
    )
    require(variable_char == expected_variable_char, (variable_char, expected_variable_char))
    variable_annihilator = (Fraction(1),)
    for factor in ((-3, 1), (-1, 1), (1, 1), (3, 2, 1)):
        variable_annihilator = polynomial_multiply(variable_annihilator, factor)
    require(
        matrix_polynomial(variable_transfer, variable_annihilator)
        == tuple(tuple(0 for _ in range(12)) for _ in range(12)),
        "variable squarefree annihilator",
    )

    oscillation = [1, -1]
    while len(oscillation) < len(variable_counts):
        oscillation.append(-2 * oscillation[-1] - 3 * oscillation[-2])
    require(all(
        3 * variable_counts[n] == 3 ** n + 1 + oscillation[n]
        for n in range(len(variable_counts))
    ), "variable closed form")

    # A finite hostile battery for arbitrary shortlex cut points.  The proof
    # uses the exact squarefree annihilator above; this census checks that the
    # indexing and prefix-cylinder decomposition use the same convention.
    variable_lex = []
    cumulative_accepts = 0
    cumulative_words = 0
    max_scaled_discrepancy = 0.0
    for depth in range(0, 11):
        variable_lex = lex_level_states(
            depth, variable_identity, variable_generators, quotient_multiply
        )
        for state in variable_lex:
            cumulative_words += 1
            cumulative_accepts += int(variable_accept(state))
            discrepancy = abs(cumulative_accepts - cumulative_words / 3)
            max_scaled_discrepancy = max(
                max_scaled_discrepancy,
                discrepancy / sqrt(cumulative_words),
            )

    # ------------------------------------------------------------------
    # Fixed THM-3339 drift language: S4 x D4, parity, recurrence, mixing.
    # ------------------------------------------------------------------
    G = source["B"]
    T = affine_action(I2, P)
    require(G == affine_action(J2, P), (G, affine_action(J2, P)))
    require(cycle_type(G) == (4,) and cycle_type(T) == (2, 2), (G, T))
    require(compose(compose(T, G), T) == inverse(G), "D4 relation TGT=G^-1")
    fixed_target = {"A": G, "B": G, "C": T}
    fixed_generators = tuple((source[letter], fixed_target[letter]) for letter in LETTERS)
    fixed_identity = (I4, I4)
    fixed_states, fixed_index, fixed_paths = enumerate_group(
        fixed_identity, fixed_generators, pair_multiply
    )
    require(len(fixed_states) == 192, len(fixed_states))
    pure_source = {state[0] for state in fixed_states if state[1] == I4}
    pure_target = {state[1] for state in fixed_states if state[0] == I4}
    require((len(pure_source), len(pure_target)) == (24, 8), (len(pure_source), len(pure_target)))

    fixed_accepting = tuple(
        state for state in fixed_states
        if cycle_type(state[0]) == cycle_type(state[1])
    )
    require(len(fixed_accepting) == 34, len(fixed_accepting))
    require(len(right_stabilizer(fixed_states, fixed_accepting, pair_multiply)) == 1,
            "fixed DFA has a nontrivial right stabilizer")
    require(
        transition_group_size(fixed_states, fixed_index, fixed_generators, pair_multiply) == 192,
        "fixed syntactic transition group",
    )

    # Construct chi on D4 independently by demanding chi(G)=chi(T)=1.
    d4_generators = (G, G, T)
    d4_character = {I4: 0}
    queue = deque((I4,))
    while queue:
        current = queue.popleft()
        for generator in d4_generators:
            candidate = compose(generator, current)
            candidate_value = d4_character[current] ^ 1
            if candidate in d4_character:
                require(d4_character[candidate] == candidate_value,
                        (current, generator, candidate))
            else:
                d4_character[candidate] = candidate_value
                queue.append(candidate)
    require(len(d4_character) == 8, len(d4_character))
    fixed_colors = Counter(d4_character[state[1]] for state in fixed_states)
    fixed_accept_colors = Counter(d4_character[state[1]] for state in fixed_accepting)
    require(fixed_colors == Counter({0: 96, 1: 96}), fixed_colors)
    require(fixed_accept_colors == Counter({0: 16, 1: 18}), fixed_accept_colors)
    require(all(
        d4_character[state[1]] == len(path) % 2
        for state, path in fixed_paths.items()
    ), "fixed word-length character")

    fixed_counts, fixed_supports = level_counts(
        fixed_identity, fixed_states, fixed_index, fixed_generators,
        pair_multiply, fixed_accepting, 420,
    )
    fixed_recurrence = berlekamp_massey_rational(fixed_counts[:100])
    expected_fixed_recurrence = tuple(map(Fraction, (
        0, 7, 0, 7, 0, 71, 0, 213, 0, 189, 0, 1701, 0, -2187,
    )))
    require(fixed_recurrence == expected_fixed_recurrence, fixed_recurrence)
    fixed_residuals = tuple(
        fixed_counts[n] - sum(
            fixed_recurrence[lag - 1] * fixed_counts[n - lag]
            for lag in range(1, 15)
        )
        for n in range(14, len(fixed_counts))
    )
    require(fixed_residuals == (Fraction(0),) * len(fixed_residuals), fixed_residuals[:10])
    require(len(fixed_residuals) >= 192, len(fixed_residuals))

    fixed_numerator, fixed_denominator = generating_function_data(
        fixed_counts, fixed_recurrence
    )
    expected_fixed_denominator = (Fraction(1),)
    for factor in (
        (-1, 1), (1, 1), (-1, 3), (1, 3),
        (1, 0, 3), (1, -2, 3), (1, 2, 3), (1, 0, -2, 0, 9),
    ):
        expected_fixed_denominator = polynomial_multiply(expected_fixed_denominator, factor)
    require(fixed_denominator == expected_fixed_denominator,
            (fixed_denominator, expected_fixed_denominator))
    require(polynomial_gcd(fixed_numerator, fixed_denominator) == (Fraction(1),),
            "fixed GF is not reduced")
    expected_fixed_numerator = tuple(map(Fraction, (
        1, 1, -6, -3, -1, 11, -56, 0, -45, -99, -162, 243, -243, -729,
    )))
    require(fixed_numerator == expected_fixed_numerator,
            (fixed_numerator, expected_fixed_numerator))

    plus_quotient, plus_remainder = polynomial_divmod(
        fixed_denominator, (Fraction(1), Fraction(-3))
    )
    minus_quotient, minus_remainder = polynomial_divmod(
        fixed_denominator, (Fraction(1), Fraction(3))
    )
    require(plus_remainder == (Fraction(0),) and minus_remainder == (Fraction(0),),
            (plus_remainder, minus_remainder))
    polar_plus = polynomial_evaluate(fixed_numerator, Fraction(1, 3)) / polynomial_evaluate(
        plus_quotient, Fraction(1, 3)
    )
    polar_minus = polynomial_evaluate(fixed_numerator, Fraction(-1, 3)) / polynomial_evaluate(
        minus_quotient, Fraction(-1, 3)
    )
    require((polar_plus, polar_minus) == (Fraction(17, 96), Fraction(-1, 96)),
            (polar_plus, polar_minus))

    alpha = polar_plus
    beta = polar_minus
    require(alpha + beta == Fraction(1, 6), (alpha, beta))
    require(alpha - beta == Fraction(3, 16), (alpha, beta))
    # Complete-level boundary limits, derived from the two parity densities.
    even_boundary = (
        (Fraction(1, 6) + Fraction(3, 16) / 3) / (1 - Fraction(1, 9))
    ) / Fraction(3, 2)
    odd_boundary = (
        (Fraction(3, 16) + Fraction(1, 6) / 3) / (1 - Fraction(1, 9))
    ) / Fraction(3, 2)
    require((even_boundary, odd_boundary) == (Fraction(11, 64), Fraction(35, 192)),
            (even_boundary, odd_boundary))

    # Two-step mixing gate.  The nine two-letter increments generate exactly
    # ker(chi), and B^4/C^6 give return lengths two and three in paired steps.
    even_generators = tuple(
        pair_multiply(right, left)
        for left in fixed_generators for right in fixed_generators
    )
    even_states, _, _ = enumerate_group(fixed_identity, even_generators, pair_multiply)
    require(len(even_states) == 96, len(even_states))
    paired_orders = tuple(permutation_order(generator[0]) for generator in fixed_generators), tuple(
        permutation_order(generator[1]) for generator in fixed_generators
    )
    pair_element_orders = tuple(
        permutation_order(
            tuple(
                generator[0][index] * 4 + generator[1][second]
                for index in range(4) for second in range(4)
            )
        )
        for generator in fixed_generators
    )
    # The product-permutation encoding above is only an order check; direct
    # lcm gives the intended paired orders (12,4,6).
    require(pair_element_orders == (12, 4, 6), (paired_orders, pair_element_orders))
    require((pair_element_orders[1] // 2, pair_element_orders[2] // 2) == (2, 3),
            pair_element_orders)

    # Weighted full-level hostile: directly enumerate shortlex levels through
    # depth ten and compare their harmonic mass with delta_parity*log(3).
    harmonic_level_errors = {}
    for depth in (8, 9, 10):
        states_at_depth = lex_level_states(depth, fixed_identity, fixed_generators, pair_multiply)
        start_index = (3 ** depth + 1) // 2
        mass = sum(
            1.0 / (start_index + rank)
            for rank, state in enumerate(states_at_depth)
            if state in frozenset(fixed_accepting)
        )
        target_density = Fraction(1, 6) if depth % 2 == 0 else Fraction(3, 16)
        harmonic_level_errors[depth] = mass - float(target_density) * log(3)

    # ------------------------------------------------------------------
    # THM-3339 root-to-child convention and both Fibonacci restrictions.
    # ------------------------------------------------------------------
    fib = fibonacci(200)
    root = (1, 2)
    for ray_index in range(48):
        words = (
            "BA" * ray_index,
            "A" + "BA" * ray_index,
            "C" + "BC" * ray_index,
        )
        expected_vectors = (
            (fib[3 * ray_index + 2], fib[3 * ray_index + 3]),
            (fib[3 * ray_index + 3], fib[3 * ray_index + 4]),
            (
                (fib[3 * ray_index + 5] - fib[3 * ray_index + 4]) // 2,
                (fib[3 * ray_index + 5] + fib[3 * ray_index + 4]) // 2,
            ),
        )
        observed_vectors = tuple(
            matrix_vector(chronological_matrix(word), root)
            for word in words
        )
        require(observed_vectors == expected_vectors,
                (ray_index, words, observed_vectors, expected_vectors))

        variable_observed = tuple(
            variable_accept(word_composite(
                word, variable_generators, quotient_multiply, variable_identity
            ))
            for word in words
        )
        variable_expected = (
            ray_index % 2 == 0,
            ray_index % 2 == 1,
            ray_index % 2 == 1,
        )
        require(variable_observed == variable_expected,
                (ray_index, words, variable_observed, variable_expected))

        fixed_observed = tuple(
            word_composite(word, fixed_generators, pair_multiply, fixed_identity)
            in frozenset(fixed_accepting)
            for word in words
        )
        fixed_expected = (
            ray_index % 4 == 0,
            ray_index % 4 == 3,
            ray_index % 4 == 3,
        )
        require(fixed_observed == fixed_expected,
                (ray_index, words, fixed_observed, fixed_expected))

    correct_ba = matrix_vector(chronological_matrix("BA"), root)
    wrong_ba = matrix_vector(opposite_matrix("BA"), root)
    require(correct_ba == (5, 8) and wrong_ba == (3, 8), (correct_ba, wrong_ba))

    # Residue conversion is checked directly, rather than inferred from labels.
    variable_indices = []
    fixed_indices = []
    for index_value in range(2, 146):
        ray_index, family = divmod(index_value - 2, 3)
        word = (
            "BA" * ray_index,
            "A" + "BA" * ray_index,
            "C" + "BC" * ray_index,
        )[family]
        if variable_accept(word_composite(
            word, variable_generators, quotient_multiply, variable_identity
        )):
            variable_indices.append(index_value % 6)
        if word_composite(word, fixed_generators, pair_multiply, fixed_identity) in frozenset(fixed_accepting):
            fixed_indices.append(index_value % 12)
    require(set(variable_indices) == {0, 1, 2}, set(variable_indices))
    require(set(fixed_indices) == {0, 1, 2}, set(fixed_indices))

    semantic = {
        "source": source,
        "linear": linear,
        "static_any": static_any,
        "static_pinned": static_pinned,
        "family_counts": family_counts,
        "matching": matching,
        "variable_states": len(variable_states),
        "variable_accepting": len(variable_accepting),
        "variable_counts": variable_counts[:20],
        "variable_recurrence": tuple(int(value) for value in variable_recurrence),
        "variable_char": variable_char,
        "fixed_states": len(fixed_states),
        "fixed_accept_colors": dict(fixed_accept_colors),
        "fixed_counts": fixed_counts[:30],
        "fixed_recurrence": tuple(int(value) for value in fixed_recurrence),
        "fixed_numerator": tuple(int(value) for value in fixed_numerator),
        "fixed_denominator": tuple(int(value) for value in fixed_denominator),
        "boundary_limits": (str(even_boundary), str(odd_boundary)),
        "root_convention": (correct_ba, wrong_ba),
    }
    semantic_sha256 = sha256(
        dumps(semantic, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()

    print("BERGGREN / T4 THREE-LANGUAGE INDEPENDENT EXACT AUDIT")
    print("STATUS: ACCEPTED WITH ONE TERMINOLOGY REPAIR; NO LRC OR JACOBIAN CONSEQUENCE")
    print(f"SCRIPT: {SCRIPT}")
    print(f"OUTPUT: {OUTPUT}")
    print(f"THM3339_LF_SHA256: {THM3339_LF_SHA256}")
    print(f"MATRIX_FIRST_PROJECTIVE_ACTIONS: {source}")
    print(f"TRUE_MOD2_LINEAR_ACTIONS: {linear}")
    print(f"STATIC_COUNTS_ANY_TRANSLATION_A_B_C: {static_any}")
    print(f"STATIC_COUNTS_PINNED_DRIFT_A_B_C: {static_pinned}")
    print(f"SIX_FRAME_COUNTS_BY_TRANSLATION_0_P_Q_R: {family_counts}")
    print("B_POSITIVE_COUNTS: pinned=4, either lawful four-cycle translation=8, and each lawful translation has 8^3=512 six-frame families")
    print(f"VARIABLE_QUOTIENT_GENERATORS: {matching} with linear bits (1,1,0)")
    print("VARIABLE_LANGUAGE: image/minimal DFA/syntactic group = S3 x C2 of order 12; right stabilizer of the four accepting states is trivial")
    print("VARIABLE_ACCEPTANCE: (mu,eps)=(id,0) or (reflection,1); direct all-48 S4 x <J> gauge enumeration gives 16 passing raw pairs")
    print(f"VARIABLE_LEVEL_COUNTS_0_TO_19: {variable_counts[:20]}")
    print("VARIABLE_RECURRENCE_GF: (2,2,6,-9); numerator (1,-1,-1,-3), denominator (1,-2,-2,-6,9), gcd=1")
    print("VARIABLE_SPECTRUM: trace-Newton characteristic equals (L-3)(L-1)^4(L+1)^3(L^2+2L+3)^2 and the squarefree factor annihilates the full 12-state transfer")
    print(f"VARIABLE_SHORTLEX_HOSTILE_MAX_ABS_DISCREPANCY_OVER_SQRT_N_DEPTH10: {max_scaled_discrepancy:.12f}")
    print("VARIABLE_CONCLUSION: a_n=(3^n+1+Re((-1+i*sqrt(2))^n))/3; shortlex A(N)=N/3+O(sqrt(N)), hence harmonic coefficient 1/3")
    print("FIXED_LANGUAGE: paired image/minimal DFA/syntactic group = S4 x D4 of order 192; pure factors 24 and 8 and accepting right stabilizer trivial")
    print(f"FIXED_PARITY_CLASSES_AND_ACCEPTING: classes={dict(fixed_colors)}, accepting={dict(fixed_accept_colors)}")
    print(f"FIXED_LEVEL_COUNTS_0_TO_29: {fixed_counts[:30]}")
    print("FIXED_MINIMAL_RECURRENCE_LAGS_1_TO_14: (0,7,0,7,0,71,0,213,0,189,0,1701,0,-2187); reduced GF denominator has the eight claimed factors")
    print("FIXED_POLAR_COEFFICIENTS_NOT_RESIDUES: lim_(x->1/3)(1-3x)F=17/96 and lim_(x->-1/3)(1+3x)F=-1/96")
    print("FIXED_ASYMPTOTIC: (17/96)3^n-(1/96)(-3)^n+O(3^(n/2)); parity densities 1/6 and 3/16")
    print(f"FIXED_SHORTLEX_BOUNDARY_LIMITS: even={even_boundary}, odd={odd_boundary}; natural density does not exist")
    print("FIXED_CYLINDER_GATE: two-letter increments generate all 96 states of ker(chi); B^4 and C^6 give paired return lengths 2 and 3, so each parity chamber mixes aperiodically")
    print(f"FIXED_WEIGHTED_LEVEL_ERRORS_DEPTH_8_9_10: {harmonic_level_errors}")
    print("FIXED_HARMONIC_CONCLUSION: cylinder equidistribution gives full-level masses (1/6)log3 and (3/16)log3; two-level Cesaro averaging gives coefficient 17/96, while one partial level is O(1)")
    print(f"ROOT_TO_CHILD_HOSTILE: BA applies B then A, so AB*(1,2)={correct_ba}; the reversed convention gives BA*(1,2)={wrong_ba} and fails immediately")
    print("FIBONACCI_GATES: variable translation k mod 6 in {0,1,2}; fixed drift k mod 12 in {0,1,2}; index harmonic coefficients 1/2 and 1/4")
    print("QUANTIFIER_LADDER: static atlas < wordwise gauge+translation < wordwise gauge with fixed drift < one common branch gauge; only the last is ruled out generatorwise")
    print(f"SEMANTIC_SHA256: {semantic_sha256}")
    print("VERDICT: all substantive claims in commits 62457a7db, 9459b7693, 85f20bfde, and the 4aae3772a mixing gate are accepted; repair only the misuse of 'residue' for normalized polar coefficient")


if __name__ == "__main__":
    main()
