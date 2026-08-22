#!/usr/bin/env python3
"""Exact companion for promoted THM-3359 modular C-finite supports.

The certificate is deliberately finite and typed.  It checks modular output
supports of finite integer linear representations; it does not identify a
tournament with a ternary relation, nor a recurrence index with a Berggren
address.  Its exhaustive universes are:

* all 64 labelled tournaments of order four and all eight cut signs modulo
  global sign, for every modulus 2 through 12;
* all 64 complete deterministic transition tables on two states and a
  three-letter alphabet, with all four accepting masks, for the same moduli;
* Fibonacci parity, the three-ray Berggren depth count, ternary return modulo
  three with and without a tie letter, and a genuinely transient recurrence;
* the THM-3356 U-spine shell clock at modulus 1105, its antipodal K4 quotient,
  coprime orbit pullbacks, and the ramified charge-five hostile;
* four finite multiplicative-scar controls through cutoff 256.

Matrix powers and permutation-determinant characteristic recurrences are
independent paths.  All support densities and collision-corrected scar
coefficients are Fraction-exact.  Validation uses RuntimeError gates rather
than assertion statements, so ordinary and optimized Python have identical
semantics and deterministic output.
"""

from __future__ import annotations

import ast
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations
from math import gcd, isqrt
from pathlib import Path


MODULI = tuple(range(2, 13))
TOURNAMENT_PAIRS = tuple(combinations(range(4), 2))
EXPECTED_TOURNAMENT_DIGEST = "5a26a103a11b34945a7219a4ece74c06b4c534cdcd90457e96dda773a6853168"
EXPECTED_DFA_DIGEST = "bc4226375938d87c753dd34f9fc00086b536eca749a3a9f7df8e1d0bd239f2c8"
EXPECTED_EXAMPLE_DIGEST = "c174578f4f3db70696028ee8736429db0eb0b3f879a106f31384b3680a141899"
EXPECTED_SHELL_CLOCK_DIGEST = "6d0cde313928acc8058c9f75304df502d52621ce85784da6c2868412d8e2f6b3"
EXPECTED_SCAR_DIGEST = "8769884c09dac2fff5647f5bb68a2b382b3db76f084ba0e696190f51abf3bb98"
EXPECTED_SEMANTIC_DIGEST = "6cc242e4abee4d13d513c5f5016043aa0e77359c66798ca5c172867bf2daf3d1"


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


def poly_mul(left, right):
    out = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i + j] += a * b
    return tuple(out)


def permutation_sign(perm):
    inversions = sum(
        1
        for i in range(len(perm))
        for j in range(i + 1, len(perm))
        if perm[i] > perm[j]
    )
    return -1 if inversions % 2 else 1


def characteristic_polynomial(matrix):
    """Return det(lambda I-M), with coefficients in ascending order."""
    degree = len(matrix)
    require(degree >= 1, "empty matrix")
    require(all(len(row) == degree for row in matrix), ("nonsquare", matrix))
    total = [0] * (degree + 1)
    for perm in permutations(range(degree)):
        term = (1,)
        for i, j in enumerate(perm):
            factor = (-matrix[i][j], 1) if i == j else (-matrix[i][j],)
            term = poly_mul(term, factor)
        sign = permutation_sign(perm)
        for index, coefficient in enumerate(term):
            total[index] += sign * coefficient
    require(total[-1] == 1, ("nonmonic characteristic polynomial", matrix, total))
    return tuple(total)


def identity_matrix(size):
    return tuple(
        tuple(1 if i == j else 0 for j in range(size))
        for i in range(size)
    )


def zero_matrix(size):
    return tuple(tuple(0 for _ in range(size)) for _ in range(size))


def matrix_add(left, right):
    return tuple(
        tuple(left[i][j] + right[i][j] for j in range(len(left)))
        for i in range(len(left))
    )


def matrix_scale(coefficient, matrix):
    return tuple(tuple(coefficient * entry for entry in row) for row in matrix)


def matrix_multiply(left, right):
    size = len(left)
    return tuple(
        tuple(
            sum(left[i][k] * right[k][j] for k in range(size))
            for j in range(size)
        )
        for i in range(size)
    )


def cayley_hamilton_gate(matrix, polynomial):
    size = len(matrix)
    power = identity_matrix(size)
    value = zero_matrix(size)
    for coefficient in polynomial:
        value = matrix_add(value, matrix_scale(coefficient, power))
        power = matrix_multiply(power, matrix)
    require(value == zero_matrix(size), ("Cayley-Hamilton failure", matrix, polynomial, value))


def matrix_vector(matrix, vector, modulus=None):
    out = tuple(
        sum(matrix[i][j] * vector[j] for j in range(len(vector)))
        for i in range(len(matrix))
    )
    if modulus is not None:
        out = tuple(value % modulus for value in out)
    return out


def dot(left, right, modulus=None):
    value = sum(a * b for a, b in zip(left, right))
    return value if modulus is None else value % modulus


def direct_terms(matrix, initial, observer, count, modulus=None):
    state = tuple(value if modulus is None else value % modulus for value in initial)
    terms = []
    for _ in range(count):
        terms.append(dot(observer, state, modulus))
        state = matrix_vector(matrix, state, modulus)
    return tuple(terms)


def companion_step(state, recurrence, modulus):
    next_value = sum(
        recurrence[index] * state[index]
        for index in range(len(state))
    ) % modulus
    return state[1:] + (next_value,)


def linear_certificate(matrix, initial, observer, modulus):
    """Audit one modular output by direct dynamics and a companion state."""
    rank = len(matrix)
    polynomial = characteristic_polynomial(matrix)
    cayley_hamilton_gate(matrix, polynomial)
    recurrence = tuple(-polynomial[index] for index in range(rank))
    seed = direct_terms(matrix, initial, observer, rank, modulus)

    seen = {}
    states = []
    state = seed
    while state not in seen:
        seen[state] = len(states)
        states.append(state)
        require(len(states) <= modulus ** rank, (
            "finite-state bound exceeded", modulus, rank, len(states)
        ))
        state = companion_step(state, recurrence, modulus)
    transient = seen[state]
    period = len(states) - transient
    require(period >= 1, ("zero period", matrix, modulus))
    require(transient + period <= modulus ** rank, (
        "pigeonhole bound", transient, period, modulus, rank
    ))

    comparison_windows = transient + 2 * period + rank + 1
    direct = direct_terms(
        matrix, initial, observer, comparison_windows + rank, modulus
    )
    recurrence_state = seed
    for index in range(comparison_windows):
        require(recurrence_state == tuple(direct[index:index + rank]), (
            "direct/recurrence state mismatch", index, modulus,
            recurrence_state, direct[index:index + rank]
        ))
        recurrence_state = companion_step(recurrence_state, recurrence, modulus)
    for index in range(transient, transient + period + rank):
        require(direct[index] == direct[index + period], (
            "eventual output period failure", index, transient, period, modulus
        ))

    if gcd(abs(polynomial[0]), modulus) == 1:
        require(transient == 0, (
            "invertible companion has transient", polynomial, modulus, transient
        ))

    histogram = [0] * modulus
    for cycle_state in states[transient:transient + period]:
        histogram[cycle_state[0]] += 1
    require(sum(histogram) == period, (histogram, period))

    return {
        "polynomial": polynomial,
        "recurrence": recurrence,
        "transient": transient,
        "period": period,
        "histogram": tuple(histogram),
        "windows": comparison_windows,
    }


def add_support_record(digest, prefix, modulus, certificate):
    transient = certificate["transient"]
    period = certificate["period"]
    histogram = certificate["histogram"]
    digest.add((
        prefix, modulus, certificate["polynomial"], certificate["recurrence"],
        transient, period, histogram
    ))
    for residue, hits in enumerate(histogram):
        density = Fraction(hits, period)
        scar = density - density * density / 2
        digest.add((
            "singleton", prefix, modulus, residue,
            fraction_key(density), fraction_key(scar)
        ))


def tournament_matrix(code):
    matrix = [[0] * 4 for _ in range(4)]
    for bit, (i, j) in enumerate(TOURNAMENT_PAIRS):
        if (code >> bit) & 1:
            matrix[i][j] = 1
        else:
            matrix[j][i] = 1
    return tuple(tuple(row) for row in matrix)


def switched_tournament_code(code, normalized_cut):
    target = code
    for bit, (i, j) in enumerate(TOURNAMENT_PAIRS):
        side_i = 0 if i == 0 else (normalized_cut >> (i - 1)) & 1
        side_j = 0 if j == 0 else (normalized_cut >> (j - 1)) & 1
        if side_i != side_j:
            target ^= 1 << bit
    return target


def tournament_class(matrix):
    scores = tuple(sorted(sum(row) for row in matrix))
    classes = {
        (0, 1, 2, 3): "transitive",
        (1, 1, 1, 3): "source_C3",
        (0, 2, 2, 2): "C3_sink",
        (1, 1, 2, 2): "strong_2211",
    }
    require(scores in classes, ("unknown order-four score class", scores, matrix))
    return classes[scores]


def audit_tournaments():
    expected_base = (
        ("C3_sink", 8),
        ("source_C3", 8),
        ("strong_2211", 24),
        ("transitive", 24),
    )
    expected_switched = (
        ("C3_sink", 64),
        ("source_C3", 64),
        ("strong_2211", 192),
        ("transitive", 192),
    )
    expected_signatures = {
        "transitive": ((0, 0, 0, 0, 1), (4, 6, 4, 1, 0, 0, 0, 0)),
        "source_C3": ((0, -1, 0, 0, 1), (4, 6, 6, 6, 6, 6, 6, 6)),
        "C3_sink": ((0, -1, 0, 0, 1), (4, 6, 6, 6, 6, 6, 6, 6)),
        "strong_2211": ((-1, -2, 0, 0, 1), (4, 6, 8, 11, 16, 22, 30, 43)),
    }

    base_counts = Counter()
    signatures = {}
    for code in range(64):
        matrix = tournament_matrix(code)
        name = tournament_class(matrix)
        base_counts[name] += 1
        signature = (
            characteristic_polynomial(matrix),
            direct_terms(matrix, (1, 1, 1, 1), (1, 1, 1, 1), 8),
        )
        if name in signatures:
            require(signatures[name] == signature, (name, signatures[name], signature))
        else:
            signatures[name] = signature
    require(tuple(sorted(base_counts.items())) == expected_base, base_counts)
    require(signatures == expected_signatures, signatures)

    digest = ExactDigest()
    switched_counts = Counter()
    target_multiplicity = Counter()
    cycle_pairs = Counter()
    state_audits = 0
    singleton_filters = 0
    windows = 0
    maximum_transient = 0
    maximum_period = 0
    strong_mod_two = set()

    for base_code in range(64):
        for cut in range(8):
            target_code = switched_tournament_code(base_code, cut)
            target_multiplicity[target_code] += 1
            matrix = tournament_matrix(target_code)
            name = tournament_class(matrix)
            switched_counts[name] += 1
            for modulus in MODULI:
                certificate = linear_certificate(
                    matrix, (1, 1, 1, 1), (1, 1, 1, 1), modulus
                )
                prefix = ("tournament", base_code, cut, target_code, name)
                add_support_record(digest, prefix, modulus, certificate)
                state_audits += 1
                singleton_filters += modulus
                windows += certificate["windows"]
                maximum_transient = max(maximum_transient, certificate["transient"])
                maximum_period = max(maximum_period, certificate["period"])
                cycle_pairs[(certificate["transient"], certificate["period"])] += 1
                if name == "strong_2211" and modulus == 2:
                    strong_mod_two.add((
                        certificate["transient"], certificate["period"],
                        certificate["histogram"]
                    ))

    require(tuple(sorted(switched_counts.items())) == expected_switched, switched_counts)
    require(set(target_multiplicity.values()) == {8}, target_multiplicity)
    require(len(target_multiplicity) == 64, len(target_multiplicity))
    require(state_audits == 5_632, state_audits)
    require(singleton_filters == 39_424, singleton_filters)
    require(strong_mod_two == {(0, 4, (3, 1))}, strong_mod_two)

    strong_matrix = next(
        tournament_matrix(code)
        for code in range(64)
        if tournament_class(tournament_matrix(code)) == "strong_2211"
    )
    strong_parity = direct_terms(
        strong_matrix, (1, 1, 1, 1), (1, 1, 1, 1), 32, 2
    )
    require(strong_parity == (0, 0, 0, 1) * 8, strong_parity)

    result = {
        "base_counts": tuple(sorted(base_counts.items())),
        "switched_counts": tuple(sorted(switched_counts.items())),
        "signatures": tuple((key, signatures[key]) for key in sorted(signatures)),
        "state_audits": state_audits,
        "singleton_filters": singleton_filters,
        "windows": windows,
        "maximum_transient": maximum_transient,
        "maximum_period": maximum_period,
        "cycle_pair_types": len(cycle_pairs),
        "cycle_pair_digest": sha256(repr(tuple(sorted(cycle_pairs.items()))).encode("ascii")).hexdigest(),
        "strong_parity": strong_parity[:8],
        "digest": digest.hexdigest(),
    }
    freeze(result["digest"], EXPECTED_TOURNAMENT_DIGEST, "tournament digest")
    return result


def dfa_transition_matrix(table_code):
    row_matrix = [[0, 0], [0, 0]]
    for state in range(2):
        for letter in range(3):
            target = (table_code >> (3 * state + letter)) & 1
            row_matrix[state][target] += 1
    return (
        (row_matrix[0][0], row_matrix[1][0]),
        (row_matrix[0][1], row_matrix[1][1]),
    )


def audit_two_state_ternary_dfas():
    digest = ExactDigest()
    cycle_pairs = Counter()
    state_audits = 0
    singleton_filters = 0
    windows = 0
    maximum_transient = 0
    maximum_period = 0
    matrix_counts = Counter()

    for table_code in range(64):
        matrix = dfa_transition_matrix(table_code)
        matrix_counts[matrix] += 1
        require(all(sum(row) == 3 for row in zip(*matrix)), (
            "DFA transition is not complete", table_code, matrix
        ))
        for accepting_mask in range(4):
            observer = (
                accepting_mask & 1,
                (accepting_mask >> 1) & 1,
            )
            for modulus in MODULI:
                certificate = linear_certificate(
                    matrix, (1, 0), observer, modulus
                )
                prefix = ("dfa", table_code, accepting_mask, matrix, observer)
                add_support_record(digest, prefix, modulus, certificate)
                state_audits += 1
                singleton_filters += modulus
                windows += certificate["windows"]
                maximum_transient = max(maximum_transient, certificate["transient"])
                maximum_period = max(maximum_period, certificate["period"])
                cycle_pairs[(certificate["transient"], certificate["period"])] += 1

    require(len(matrix_counts) == 16, ("unexpected two-state count matrices", len(matrix_counts)))
    require(sum(matrix_counts.values()) == 64, matrix_counts)
    require(state_audits == 2_816, state_audits)
    require(singleton_filters == 19_712, singleton_filters)

    result = {
        "transition_tables": 64,
        "count_matrices": len(matrix_counts),
        "languages": 256,
        "state_audits": state_audits,
        "singleton_filters": singleton_filters,
        "windows": windows,
        "maximum_transient": maximum_transient,
        "maximum_period": maximum_period,
        "cycle_pair_types": len(cycle_pairs),
        "cycle_pair_digest": sha256(repr(tuple(sorted(cycle_pairs.items()))).encode("ascii")).hexdigest(),
        "matrix_multiplicity_digest": sha256(repr(tuple(sorted(matrix_counts.items()))).encode("ascii")).hexdigest(),
        "digest": digest.hexdigest(),
    }
    freeze(result["digest"], EXPECTED_DFA_DIGEST, "DFA digest")
    return result


def harmonic_number(n):
    total = Fraction(0, 1)
    for value in range(1, n + 1):
        total += Fraction(1, value)
    return total


def residue_return_matrix(increments):
    row_matrix = [[0, 0, 0] for _ in range(3)]
    for residue in range(3):
        for increment in increments:
            row_matrix[residue][(residue + increment) % 3] += 1
    return tuple(
        tuple(row_matrix[column][row] for column in range(3))
        for row in range(3)
    )


def audit_examples():
    digest = ExactDigest()

    fibonacci_matrix = ((0, 1), (1, 1))
    fibonacci_periods = []
    fibonacci_certificates = {}
    for modulus in MODULI:
        certificate = linear_certificate(
            fibonacci_matrix, (0, 1), (1, 0), modulus
        )
        fibonacci_certificates[modulus] = certificate
        add_support_record(digest, ("Fibonacci",), modulus, certificate)
        fibonacci_periods.append((
            modulus, certificate["transient"], certificate["period"]
        ))
    expected_fibonacci_periods = (
        (2, 0, 3), (3, 0, 8), (4, 0, 6), (5, 0, 20),
        (6, 0, 24), (7, 0, 16), (8, 0, 12), (9, 0, 24),
        (10, 0, 60), (11, 0, 10), (12, 0, 24),
    )
    require(tuple(fibonacci_periods) == expected_fibonacci_periods, fibonacci_periods)
    require(fibonacci_certificates[2]["histogram"] == (1, 2), fibonacci_certificates[2])

    fibonacci_values = [0, 1]
    while len(fibonacci_values) <= 512:
        fibonacci_values.append(fibonacci_values[-1] + fibonacci_values[-2])
    fibonacci_harmonic_checks = 0
    for cutoff in range(1, 513):
        direct = sum(
            (Fraction(1, index) for index in range(1, cutoff + 1)
             if fibonacci_values[index] % 2),
            Fraction(0, 1),
        )
        formula = harmonic_number(cutoff) - harmonic_number(cutoff // 3) / 3
        require(direct == formula, ("Fibonacci harmonic identity", cutoff, direct, formula))
        digest.add(("Fibonacci harmonic", cutoff, fraction_key(direct)))
        fibonacci_harmonic_checks += 1

    berggren_matrix = ((0, 1), (1, 0))
    berggren_periods = []
    berggren_certificates = {}
    for modulus in MODULI:
        certificate = linear_certificate(
            berggren_matrix, (1, 2), (1, 0), modulus
        )
        berggren_certificates[modulus] = certificate
        add_support_record(digest, ("Berggren depth",), modulus, certificate)
        berggren_periods.append((
            modulus, certificate["transient"], certificate["period"]
        ))
    require(tuple(berggren_periods) == tuple((m, 0, 2) for m in MODULI), berggren_periods)
    require(berggren_certificates[2]["histogram"] == (1, 1), berggren_certificates[2])

    berggren_level_counts = []
    berggren_harmonic_checks = 0
    for depth in range(512):
        first_ray = 1 if depth % 2 == 0 else 0
        second_ray = 1 if depth % 2 == 1 else 0
        third_ray = 1 if depth % 2 == 1 else 0
        direct_count = first_ray + second_ray + third_ray
        recurrence_count = 1 if depth % 2 == 0 else 2
        require(direct_count == recurrence_count, (
            "Berggren three-ray depth count", depth, direct_count, recurrence_count
        ))
        berggren_level_counts.append(direct_count)
        digest.add(("Berggren level", depth, direct_count))
    require(Counter(berggren_level_counts) == Counter({1: 256, 2: 256}), (
        "Berggren depth census", Counter(berggren_level_counts)
    ))
    for cutoff in range(1, 513):
        direct = sum(
            (Fraction(1, depth) for depth in range(1, cutoff + 1)
             if (1 if depth % 2 == 0 else 2) % 2),
            Fraction(0, 1),
        )
        formula = harmonic_number(cutoff // 2) / 2
        require(direct == formula, ("Berggren harmonic identity", cutoff, direct, formula))
        digest.add(("Berggren harmonic", cutoff, fraction_key(direct)))
        berggren_harmonic_checks += 1

    tie_matrix = residue_return_matrix((-1, 0, 1))
    no_tie_matrix = residue_return_matrix((-1, 1))
    tie_exact = direct_terms(tie_matrix, (1, 0, 0), (1, 0, 0), 257)
    no_tie_exact = direct_terms(no_tie_matrix, (1, 0, 0), (1, 0, 0), 257)
    for length in range(257):
        tie_formula = 1 if length == 0 else 3 ** (length - 1)
        no_tie_numerator = 2 ** length + 2 * ((-1) ** length)
        require(no_tie_numerator % 3 == 0, ("no-tie divisibility", length))
        no_tie_formula = no_tie_numerator // 3
        require(tie_exact[length] == tie_formula, (
            "tie return formula", length, tie_exact[length], tie_formula
        ))
        require(no_tie_exact[length] == no_tie_formula, (
            "no-tie return formula", length, no_tie_exact[length], no_tie_formula
        ))
        if length >= 1:
            require(tie_exact[length] % 2 == 1, ("tie parity", length))
            require(no_tie_exact[length] % 2 == 0, ("no-tie parity", length))
        digest.add(("ternary return", length, tie_exact[length], no_tie_exact[length]))
    tie_certificate = linear_certificate(tie_matrix, (1, 0, 0), (1, 0, 0), 2)
    no_tie_certificate = linear_certificate(no_tie_matrix, (1, 0, 0), (1, 0, 0), 2)
    require((tie_certificate["transient"], tie_certificate["period"], tie_certificate["histogram"]) == (0, 1, (0, 1)), tie_certificate)
    require((no_tie_certificate["transient"], no_tie_certificate["period"], no_tie_certificate["histogram"]) == (1, 1, (1, 0)), no_tie_certificate)
    add_support_record(digest, ("ternary tie return mod 3",), 2, tie_certificate)
    add_support_record(digest, ("binary no-tie return mod 3",), 2, no_tie_certificate)

    transient_matrix = ((2,),)
    transient_certificate = linear_certificate(transient_matrix, (1,), (1,), 4)
    require((
        transient_certificate["transient"], transient_certificate["period"],
        transient_certificate["histogram"]
    ) == (2, 1, (1, 0, 0, 0)), transient_certificate)
    transient_terms = direct_terms(transient_matrix, (1,), (1,), 16, 4)
    finite_support = tuple(
        index for index in range(1, len(transient_terms))
        if transient_terms[index] == 2
    )
    require(finite_support == (1,), finite_support)
    add_support_record(digest, ("transient powers of two",), 4, transient_certificate)

    densities = {
        "strong_odd": Fraction(1, 4),
        "Fibonacci_odd": Fraction(2, 3),
        "Berggren_depth_odd": Fraction(1, 2),
        "ternary_tie_odd": Fraction(1, 1),
        "ternary_no_tie_odd": Fraction(0, 1),
        "transient_residue_two": Fraction(0, 1),
    }
    scars = {
        key: value - value * value / 2
        for key, value in densities.items()
    }
    require(scars["strong_odd"] == Fraction(7, 32), scars)
    require(scars["Fibonacci_odd"] == Fraction(4, 9), scars)
    require(scars["Berggren_depth_odd"] == Fraction(3, 8), scars)
    require(scars["ternary_tie_odd"] == Fraction(1, 2), scars)
    for key in sorted(densities):
        digest.add((key, fraction_key(densities[key]), fraction_key(scars[key])))

    result = {
        "fibonacci_periods": tuple(fibonacci_periods),
        "fibonacci_mod2_histogram": fibonacci_certificates[2]["histogram"],
        "fibonacci_harmonic_checks": fibonacci_harmonic_checks,
        "berggren_periods": tuple(berggren_periods),
        "berggren_mod2_histogram": berggren_certificates[2]["histogram"],
        "berggren_level_census": tuple(sorted(Counter(berggren_level_counts).items())),
        "berggren_harmonic_checks": berggren_harmonic_checks,
        "tie_cycle": (
            tie_certificate["transient"], tie_certificate["period"],
            tie_certificate["histogram"]
        ),
        "no_tie_cycle": (
            no_tie_certificate["transient"], no_tie_certificate["period"],
            no_tie_certificate["histogram"]
        ),
        "transient_cycle": (
            transient_certificate["transient"], transient_certificate["period"],
            transient_certificate["histogram"]
        ),
        "finite_support": finite_support,
        "finite_harmonic_mass": fraction_key(Fraction(1, 1)),
        "densities": tuple((key, fraction_key(densities[key])) for key in sorted(densities)),
        "scars": tuple((key, fraction_key(scars[key])) for key in sorted(scars)),
        "digest": digest.hexdigest(),
    }
    freeze(result["digest"], EXPECTED_EXAMPLE_DIGEST, "example digest")
    return result


def shell_norm(charge, parameter):
    """Norm of charge*(1,0)+parameter*(1,1)."""
    return 2 * parameter * parameter + 2 * charge * parameter + charge * charge


def audit_shell_clock():
    modulus = 1_105
    base_roots = tuple(
        parameter
        for parameter in range(modulus)
        if shell_norm(1, parameter) % modulus == 0
    )
    expected_roots = (23, 231, 418, 431, 673, 686, 873, 1081)
    require(base_roots == expected_roots, ("U-spine roots", base_roots))

    local_roots = tuple(
        (prime, tuple(
            parameter
            for parameter in range(prime)
            if shell_norm(1, parameter) % prime == 0
        ))
        for prime in (5, 13, 17)
    )
    require(all(len(roots) == 2 for _, roots in local_roots), local_roots)

    unpaired = set(base_roots)
    vertices = []
    while unpaired:
        root = min(unpaired)
        antipode = (-1 - root) % modulus
        require(antipode in unpaired, ("missing antipode", root, antipode, unpaired))
        require(shell_norm(1, antipode) % modulus == 0, (root, antipode))
        vertex = tuple(sorted((root, antipode)))
        vertices.append(vertex)
        unpaired.remove(root)
        unpaired.remove(antipode)
    vertices = tuple(vertices)
    expected_vertices = (
        (23, 1081), (231, 873), (418, 686), (431, 673)
    )
    require(vertices == expected_vertices, ("antipodal K4 vertices", vertices))

    edges = []
    for left, right in combinations(range(4), 2):
        support = tuple(sorted(vertices[left] + vertices[right]))
        require(len(support) == 4 and len(set(support)) == 4, (
            "K4 edge support", left, right, support
        ))
        edges.append(((left, right), support))
    edges = tuple(edges)
    require(len({support for _, support in edges}) == 6, edges)

    matching_edge_pairs = (
        ((0, 1), (2, 3)),
        ((0, 2), (1, 3)),
        ((0, 3), (1, 2)),
    )
    matching_supports = []
    edge_dictionary = dict(edges)
    for matching in matching_edge_pairs:
        support = tuple(sorted(
            set(edge_dictionary[matching[0]]) | set(edge_dictionary[matching[1]])
        ))
        require(support == base_roots, ("perfect matching support", matching, support))
        matching_supports.append((matching, support))
    matching_supports = tuple(matching_supports)

    harmonic_poles = (
        ("vertex", fraction_key(Fraction(2, modulus))),
        ("edge", fraction_key(Fraction(4, modulus))),
        ("perfect_matching", fraction_key(Fraction(8, modulus))),
    )

    controls = ((1, 0), (2, 7), (3, 41), (7, 1104), (-3, 29))
    control_records = []
    period_checks = 0
    digest = ExactDigest()
    digest.add(("base", modulus, local_roots, base_roots, vertices, edges,
                matching_supports, harmonic_poles))
    for charge, orbit in controls:
        require(gcd(charge, modulus) == 1, ("noncoprime control", charge, modulus))
        charge_inverse = pow(charge % modulus, -1, modulus)
        charge_roots = tuple(
            parameter
            for parameter in range(modulus)
            if shell_norm(charge, parameter) % modulus == 0
        )
        scaled_roots = tuple(sorted((charge * root) % modulus for root in base_roots))
        require(charge_roots == scaled_roots, (
            "scaled shell roots", charge, charge_roots, scaled_roots
        ))
        pullback = tuple(sorted(
            (charge_inverse * (root - orbit)) % modulus
            for root in charge_roots
        ))
        direct = tuple(
            clock
            for clock in range(modulus)
            if shell_norm(charge, orbit + charge * clock) % modulus == 0
        )
        require(direct == pullback, (
            "orbit clock pullback", charge, orbit, direct, pullback
        ))
        require(len(direct) == 8, ("clock root count", charge, orbit, direct))
        for clock in range(modulus):
            present = shell_norm(charge, orbit + charge * clock) % modulus
            next_period = shell_norm(
                charge, orbit + charge * (clock + modulus)
            ) % modulus
            require(present == next_period, (
                "period-M shell clock", charge, orbit, clock,
                present, next_period
            ))
            period_checks += 1
        record = (charge, orbit, charge_roots, pullback)
        control_records.append(record)
        digest.add(("coprime control", record))
    control_records = tuple(control_records)
    require(period_checks == len(controls) * modulus, period_checks)
    digest.add(("period checks", period_checks))

    hostile_modulus = 5
    hostile_charge = 5
    hostile_orbits = []
    hostile_roots = tuple(
        parameter
        for parameter in range(hostile_modulus)
        if shell_norm(hostile_charge, parameter) % hostile_modulus == 0
    )
    require(hostile_roots == (0,), ("ramified roots", hostile_roots))
    for orbit in range(hostile_modulus):
        hits = tuple(
            clock
            for clock in range(hostile_modulus)
            if shell_norm(
                hostile_charge, orbit + hostile_charge * clock
            ) % hostile_modulus == 0
        )
        hostile_orbits.append((orbit, hits))
    hostile_orbits = tuple(hostile_orbits)
    require(hostile_orbits == (
        (0, (0, 1, 2, 3, 4)),
        (1, ()), (2, ()), (3, ()), (4, ()),
    ), hostile_orbits)
    digest.add(("ramified hostile", hostile_modulus, hostile_charge,
                hostile_roots, hostile_orbits))

    result = {
        "modulus": modulus,
        "local_roots": local_roots,
        "roots": base_roots,
        "vertices": vertices,
        "edges": edges,
        "matching_supports": matching_supports,
        "support_sizes": (2, 4, 8),
        "harmonic_poles": harmonic_poles,
        "controls": control_records,
        "period_checks": period_checks,
        "hostile_roots": hostile_roots,
        "hostile_orbits": hostile_orbits,
        "digest": digest.hexdigest(),
    }
    freeze(result["digest"], EXPECTED_SHELL_CLOCK_DIGEST, "shell-clock digest")
    return result


def audit_multiplicative_scars():
    controls = (
        ("strong_mod4", lambda value: value % 4 == 3),
        ("Fibonacci_mod3", lambda value: value % 3 != 0),
        ("Berggren_even_depth", lambda value: value % 2 == 0),
        ("finite_transient", lambda value: value == 1),
    )
    digest = ExactDigest()
    observable_checks = 0
    pointwise_checks = 0

    for name, predicate in controls:
        cumulative_incidence = 0
        cumulative_weak = 0
        cumulative_strict = 0
        for product in range(1, 257):
            incidence = sum(
                1 for divisor in range(1, product + 1)
                if product % divisor == 0 and predicate(divisor)
            )
            collision = 0
            diagonal = 0
            direct_weak = 0
            direct_strict = 0
            for left in range(1, isqrt(product) + 1):
                if product % left:
                    continue
                right = product // left
                if left > right:
                    continue
                left_hole = predicate(left)
                right_hole = predicate(right)
                if left_hole or right_hole:
                    direct_weak += 1
                    if left < right:
                        direct_strict += 1
                if left < right and left_hole and right_hole:
                    collision += 1
                if left == right and left_hole:
                    diagonal += 1
            require(direct_weak == incidence - collision, (
                "weak scar identity", name, product, direct_weak,
                incidence, collision
            ))
            require(direct_strict == incidence - collision - diagonal, (
                "strict scar identity", name, product, direct_strict,
                incidence, collision, diagonal
            ))
            pointwise_checks += 2
            cumulative_incidence += incidence
            cumulative_weak += direct_weak
            cumulative_strict += direct_strict
            floor_sum = sum(
                product // hole
                for hole in range(1, product + 1)
                if predicate(hole)
            )
            require(cumulative_incidence == floor_sum, (
                "incidence floor sum", name, product,
                cumulative_incidence, floor_sum
            ))
            digest.add((
                name, product, cumulative_incidence,
                cumulative_weak, cumulative_strict
            ))
            observable_checks += 3

    require(observable_checks == 3_072, observable_checks)
    require(pointwise_checks == 2_048, pointwise_checks)
    result = {
        "controls": tuple(name for name, _ in controls),
        "cutoffs": 256,
        "observable_checks": observable_checks,
        "pointwise_identity_checks": pointwise_checks,
        "digest": digest.hexdigest(),
    }
    freeze(result["digest"], EXPECTED_SCAR_DIGEST, "scar digest")
    return result


def source_audit():
    source_path = Path(__file__).resolve()
    text = source_path.read_text(encoding="utf-8")
    tree = ast.parse(text, filename=str(source_path))
    assertion_nodes = sum(
        1 for node in ast.walk(tree) if type(node).__name__ == "Assert"
    )
    floating_literals = sum(
        1 for node in ast.walk(tree)
        if type(node).__name__ == "Constant" and type(getattr(node, "value", None)) is float
    )
    require(assertion_nodes == 0, ("assertion nodes", assertion_nodes))
    require(floating_literals == 0, ("floating literals", floating_literals))
    return lf_hash(source_path), assertion_nodes, floating_literals


def main():
    source_digest, assertion_nodes, floating_literals = source_audit()
    tournaments = audit_tournaments()
    dfas = audit_two_state_ternary_dfas()
    examples = audit_examples()
    shell_clock = audit_shell_clock()
    scars = audit_multiplicative_scars()

    semantic_payload = (
        tournaments["digest"], dfas["digest"], examples["digest"],
        shell_clock["digest"], scars["digest"],
        tournaments["base_counts"], tournaments["switched_counts"],
        tournaments["state_audits"], tournaments["singleton_filters"],
        dfas["state_audits"], dfas["singleton_filters"],
        examples["densities"], examples["scars"],
    )
    semantic_digest = sha256(repr(semantic_payload).encode("ascii")).hexdigest()
    freeze(semantic_digest, EXPECTED_SEMANTIC_DIGEST, "semantic digest")

    print("THM-3359 MODULAR C-FINITE HARMONIC-SUPPORT EXACT COMPANION")
    print("source_sha256_lf", source_digest)
    print("assertion_nodes", assertion_nodes, "floating_literals", floating_literals)
    print("moduli", MODULI)
    print("TOURNAMENT ORDER FOUR + NORMALIZED CUT CUBE")
    print("labelled_class_counts", tournaments["base_counts"])
    print("normalized_cut_signs_per_tournament", 8, "switch_cases", 512)
    print("switched_target_class_counts", tournaments["switched_counts"])
    print("class_signatures_charpoly_then_w0_to_w7", tournaments["signatures"])
    print("switch_modulus_state_audits", tournaments["state_audits"],
          "direct_recurrence_windows", tournaments["windows"])
    print("tournament_singleton_residue_filters", tournaments["singleton_filters"])
    print("tournament_cycle_max_mu_p", (tournaments["maximum_transient"], tournaments["maximum_period"]),
          "cycle_pair_types", tournaments["cycle_pair_types"],
          "cycle_pair_digest", tournaments["cycle_pair_digest"])
    print("strong_class_mod2_cycle", (0, 4, (3, 1)),
          "parity_w0_to_w7", tournaments["strong_parity"],
          "density", (1, 4), "scar", (7, 32))
    print("tournament_semantic_sha256", tournaments["digest"])
    print("COMPLETE TWO-STATE TERNARY DFAS")
    print("transition_tables", dfas["transition_tables"],
          "count_matrices", dfas["count_matrices"],
          "accepting_masks", 4, "languages", dfas["languages"])
    print("dfa_modulus_state_audits", dfas["state_audits"],
          "direct_recurrence_windows", dfas["windows"])
    print("dfa_singleton_residue_filters", dfas["singleton_filters"])
    print("dfa_cycle_max_mu_p", (dfas["maximum_transient"], dfas["maximum_period"]),
          "cycle_pair_types", dfas["cycle_pair_types"],
          "cycle_pair_digest", dfas["cycle_pair_digest"])
    print("dfa_count_matrix_multiplicity_sha256", dfas["matrix_multiplicity_digest"])
    print("dfa_semantic_sha256", dfas["digest"])
    print("FIBONACCI + BERGGREN DEPTH")
    print("fibonacci_modulus_cycles", examples["fibonacci_periods"])
    print("fibonacci_mod2_histogram", examples["fibonacci_mod2_histogram"],
          "density", (2, 3), "scar", (4, 9),
          "exact_harmonic_checks", examples["fibonacci_harmonic_checks"])
    print("berggren_depth_modulus_cycles", examples["berggren_periods"])
    print("berggren_depth_mod2_histogram", examples["berggren_mod2_histogram"],
          "density", (1, 2), "scar", (3, 8),
          "level_census_d0_to_d511", examples["berggren_level_census"],
          "exact_harmonic_checks", examples["berggren_harmonic_checks"])
    print("TERNARY TIE / NO-TIE MOD-3 HOSTILE")
    print("tie_return_count", "1_at_n0_then_3^(n-1)",
          "mod2_cycle", examples["tie_cycle"],
          "positive_odd_support", "all", "density", (1, 1), "scar", (1, 2))
    print("no_tie_return_count", "(2^n+2*(-1)^n)/3",
          "mod2_cycle", examples["no_tie_cycle"],
          "positive_odd_support", "empty", "density", (0, 1), "scar", (0, 1))
    print("FINITE / TRANSIENT HOSTILE")
    print("a_n=2^n_mod4_cycle", examples["transient_cycle"],
          "residue2_positive_support", examples["finite_support"],
          "cycle_density", (0, 1),
          "finite_harmonic_mass", examples["finite_harmonic_mass"])
    print("example_densities", examples["densities"])
    print("collision_corrected_scar_coefficients", examples["scars"])
    print("example_semantic_sha256", examples["digest"])
    print("THM-3356 U-SPINE SHELL CLOCK")
    print("modulus", shell_clock["modulus"],
          "local_roots", shell_clock["local_roots"])
    print("Q(n)=2n^2+2n+1_roots", shell_clock["roots"])
    print("antipodal_n_to_-1-n_K4_vertices", shell_clock["vertices"])
    print("K4_edge_endpoint_supports", shell_clock["edges"])
    print("K4_perfect_matching_supports", shell_clock["matching_supports"])
    print("K4_support_sizes_vertex_edge_perfect_matching", shell_clock["support_sizes"])
    print("clock_harmonic_poles", shell_clock["harmonic_poles"])
    print("coprime_charge_orbit_controls_charge_r_roots_pullback",
          shell_clock["controls"], "period_M_checks", shell_clock["period_checks"])
    print("ramified_c=M=5_roots", shell_clock["hostile_roots"],
          "orbit_clock_hits_t0_to_t4", shell_clock["hostile_orbits"])
    print("shell_clock_semantic_sha256", shell_clock["digest"])
    print("MULTIPLICATIVE-SCAR POINTWISE CONTROLS")
    print("controls", scars["controls"], "cutoffs", scars["cutoffs"],
          "three_observable_checks", scars["observable_checks"],
          "pointwise_identity_checks", scars["pointwise_identity_checks"])
    print("scar_semantic_sha256", scars["digest"])
    print("semantic_sha256", semantic_digest)
    print("replay_contract ordinary_and_python_O_outputs_byte_identical")
    print("status=ALL EXACT CHECKS PASSED; index carriers and lost sidecars remain typed")


if __name__ == "__main__":
    main()
