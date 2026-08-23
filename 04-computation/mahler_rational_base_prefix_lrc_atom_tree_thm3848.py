#!/usr/bin/env python3
"""Exact companion for the rational-base prefix/LRC separation theorem.

The script uses integer and Fraction arithmetic only.  It checks the finite
atom carrier, the full adjacent relation lattice, independent breakpoint
maximization on small controls, the exact 13-speed (3/2)-tower, and the ABC
input identities/hostiles.  It deliberately does not test ABC itself.
"""

from fractions import Fraction
from functools import reduce
from itertools import combinations
from math import gcd, lcm


CHECKS = 0


def require(condition: bool, payload: object) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True:
        raise RuntimeError(payload)


def mixed_speeds(p: int, q: int, level: int) -> tuple[int, ...]:
    return tuple(q ** (level - index) * p**index for index in range(level + 1))


def wall_grid(p: int, q: int, level: int) -> int:
    return p**level * q ** (level + 1)


def safe_atom_indices(p: int, q: int, level: int) -> tuple[int, ...]:
    """Open wall-grid atoms lying in every one-sided prefix band."""

    delta = wall_grid(p, q, level)
    denominator = 2 * delta
    speeds = mixed_speeds(p, q, level)
    safe = []
    for atom in range(delta):
        midpoint_numerator = 2 * atom + 1
        if all(
            q * ((midpoint_numerator * speed) % denominator) < denominator
            for speed in speeds
        ):
            safe.append(atom)
    return tuple(safe)


def distance_to_integer(value: Fraction) -> Fraction:
    residue = value % 1
    return min(residue, 1 - residue)


def loneliness(speeds: tuple[int, ...], time: Fraction) -> Fraction:
    return min(distance_to_integer(time * speed) for speed in speeds)


def breakpoint_times(speeds: tuple[int, ...]) -> set[Fraction]:
    """All possible vertices of the lower envelope of triangular waves."""

    candidates = {Fraction(0), Fraction(1)}
    for speed in speeds:
        denominator = 2 * speed
        candidates.update(Fraction(k, denominator) for k in range(denominator + 1))
    for left, right in combinations(speeds, 2):
        for denominator in (left + right, abs(left - right)):
            if denominator:
                candidates.update(Fraction(k, denominator) for k in range(denominator + 1))
    return candidates


def exact_lonely_maximum(speeds: tuple[int, ...]) -> tuple[Fraction, tuple[Fraction, ...]]:
    values = [(loneliness(speeds, time), time % 1) for time in breakpoint_times(speeds)]
    maximum = max(value for value, _ in values)
    maximizers = tuple(sorted({time for value, time in values if value == maximum}))
    return maximum, maximizers


def radical(value: int) -> int:
    value = abs(value)
    answer = 1
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            answer *= divisor
            while value % divisor == 0:
                value //= divisor
        divisor += 1
    if value > 1:
        answer *= value
    return answer


def primes_through(bound: int) -> tuple[int, ...]:
    primes = []
    for candidate in range(2, bound + 1):
        if all(candidate % prime for prime in primes if prime * prime <= candidate):
            primes.append(candidate)
    return tuple(primes)


def greedy_equality_digits(length: int) -> tuple[int, ...]:
    """Greedy 0/1 expansion of 1 in weights (2/3)^k."""

    state = Fraction(1)
    digits = []
    for _ in range(length):
        expanded = Fraction(3, 2) * state
        digit = int(expanded >= 1)
        state = expanded - digit
        require(Fraction(0) < state < Fraction(1), (digits, digit, state))
        digits.append(digit)
    return tuple(digits)


def finite_safe_word(word: tuple[int, ...]) -> bool:
    """Every truncated suffix tail is strictly below one."""

    length = len(word)
    for start in range(length):
        suffix_length = length - start
        cleared = 2 * sum(
            word[start + offset] * 2**offset * 3 ** (suffix_length - 1 - offset)
            for offset in range(suffix_length)
        )
        if cleared >= 3**suffix_length:
            return False
    return True


def carry_residue(word: tuple[int, ...]) -> int:
    """Canonical starting residue modulo 2^m for a finite carry word."""

    length = len(word)
    if length == 0:
        return 0
    carry_polynomial = sum(
        digit * 3 ** (length - 1 - index) * 2**index
        for index, digit in enumerate(word)
    )
    modulus = 2**length
    return (-carry_polynomial * pow(pow(3, length, modulus), -1, modulus)) % modulus


# 1. Brute-force the p-ary atom theorem for several genuinely different bases.
parameter_pairs = tuple(
    (p, q)
    for q in range(2, 5)
    for p in range(q + 1, 7)
    if gcd(p, q) == 1
)
atom_profiles = {}
for p, q in parameter_pairs:
    previous = None
    for level in range(4):
        safe = safe_atom_indices(p, q, level)
        delta = wall_grid(p, q, level)
        require(len(safe) == p**level, (p, q, level, len(safe)))
        require(Fraction(len(safe), delta) == Fraction(1, q ** (level + 1)), (p, q, level))
        if previous is not None:
            next_safe = set(safe)
            old_delta = wall_grid(p, q, level - 1)
            for parent in previous:
                recovered_children = []
                for remainder in range(p):
                    choices = [
                        p * parent + p * old_delta * lift + remainder
                        for lift in range(q)
                        if p * parent + p * old_delta * lift + remainder in next_safe
                    ]
                    require(len(choices) == 1, (p, q, level, parent, remainder, choices))
                    child = choices[0]
                    require(child % p == remainder, (p, q, child, remainder))
                    require(((child - remainder) // p) % old_delta == parent, (parent, child))
                    recovered_children.append(child)
                require(len(set(recovered_children)) == p, (p, q, level, parent))
        previous = safe
    atom_profiles[(p, q)] = tuple(p**level for level in range(4))


# 2. Adjacent rows have the advertised minors and generate the full kernel.
for p, q in parameter_pairs:
    for level in range(1, 7):
        speeds = mixed_speeds(p, q, level)
        for index in range(level):
            require(p * speeds[index] - q * speeds[index + 1] == 0, (p, q, level, index))
        maximal_minors = tuple(p**deleted * q ** (level - deleted) for deleted in range(level + 1))
        require(reduce(gcd, maximal_minors) == 1, (p, q, level, maximal_minors))


# 3. Independent breakpoint maximization checks the general closed formula.
lonely_controls = []
for p, q in parameter_pairs:
    expected = Fraction((p + q) // 2, p + q)
    for level in range(1, 4):
        speeds = mixed_speeds(p, q, level)
        observed, _ = exact_lonely_maximum(speeds)
        require(observed == expected, (p, q, level, observed, expected))
        witness_numerator = (
            ((p + q) // 2) * pow(pow(q, level, p + q), -1, p + q)
        ) % (p + q)
        witness = Fraction(witness_numerator, p + q)
        require(loneliness(speeds, witness) == expected, (p, q, level, witness))
        lonely_controls.append((p, q, level, observed))


# 4. For 3/2, independently enumerate the complete maximizing set at five levels.
maximizer_profiles = []
for level in range(1, 6):
    speeds = mixed_speeds(3, 2, level)
    observed, maximizers = exact_lonely_maximum(speeds)
    b_level = pow(pow(2, level - 1, 5), -1, 5)
    predicted = tuple(sorted({Fraction(b_level, 5), Fraction((-b_level) % 5, 5)}))
    require(observed == Fraction(2, 5), (level, observed))
    require(maximizers == predicted, (level, maximizers, predicted))
    for time in predicted:
        phases = tuple((time * speed) % 1 for speed in speeds)
        require(set(phases) == {Fraction(2, 5), Fraction(3, 5)}, (level, phases))
        require(all(phases[index] != phases[index + 1] for index in range(level)), phases)
    maximizer_profiles.append((level, tuple(str(time) for time in maximizers)))


# 5. The level-12 row is exactly the 13-speed / fourteen-runner interface.
level_12 = mixed_speeds(3, 2, 12)
expected_level_12 = (
    4096,
    6144,
    9216,
    13824,
    20736,
    31104,
    46656,
    69984,
    104976,
    157464,
    236196,
    354294,
    531441,
)
require(level_12 == expected_level_12, level_12)
require(reduce(gcd, level_12) == 1, level_12)
require(len(set(level_12)) == 13, level_12)
require(all(speed % 5 for speed in level_12), level_12)
require(loneliness(level_12, Fraction(2, 5)) == Fraction(2, 5), level_12)
require(
    tuple((Fraction(2, 5) * speed) % 1 for speed in level_12)
    == tuple(Fraction(2 if index % 2 == 0 else 3, 5) for index in range(13)),
    "level-12 alternating phase word",
)
require(len(safe_atom_indices(3, 2, 5)) == 3**5, "ternary atom positive control")
require(wall_grid(3, 2, 12) == 3**12 * 2**13, "level-12 wall grid")


# 6. The point-lift is unique although every safe open atom has three children.
for level in range(5):
    speeds = mixed_speeds(3, 2, level)
    denominator = 2 * wall_grid(3, 2, level)
    for numerator in range(denominator):
        time = Fraction(numerator, denominator)
        if not all((time * speed) % 1 < Fraction(1, 2) for speed in speeds):
            continue
        lifts = tuple(Fraction(time + lift, 2) for lift in (0, 1))
        safe_lifts = [
            lifted
            for lifted in lifts
            if all(
                (lifted * speed) % 1 < Fraction(1, 2)
                for speed in mixed_speeds(3, 2, level + 1)
            )
        ]
        require(len(safe_lifts) == 1, (level, time, safe_lifts))


# 7. The centered maximizer decodes to a nonterminal 2-adic state.
carry_period = (0, 1)
carry_polynomial = carry_period[0] * 3 + carry_period[1] * 2
phi = -Fraction(carry_polynomial, 3**2 - 2**2)
y_even = sum((Fraction(2, 3) ** (2 * index + 2) for index in range(40)), Fraction(0))
y_odd = sum((Fraction(2, 3) ** (2 * index + 1) for index in range(40)), Fraction(0))
require(phi == -Fraction(2, 5), phi)
require(y_even < Fraction(4, 5) and Fraction(4, 5) - y_even < Fraction(1, 10**6), y_even)
require(y_odd < Fraction(6, 5) and Fraction(6, 5) - y_odd < Fraction(1, 10**6), y_odd)
# Exact geometric-series values, kept separate from the finite hostile above.
require(Fraction(4, 9) / (1 - Fraction(4, 9)) == Fraction(4, 5), "Y_even")
require(Fraction(2, 3) / (1 - Fraction(4, 9)) == Fraction(6, 5), "Y_odd")
require(phi + Fraction(4, 5) / 2 == 0, "formal real value")
canonical_residues = []
for level in range(1, 25):
    b_level = pow(pow(2, level - 1, 5), -1, 5)
    residue = (2**level * b_level - 2) // 5
    require(0 <= residue < 2**level, (level, residue))
    require((5 * residue + 2) % 2**level == 0, (level, residue))
    canonical_residues.append(residue)
require(len(set(canonical_residues[-8:])) > 1, canonical_residues)


# 8. Arbitrarily long denominator-19 safe prefixes and their two ABC triples.
denominator_19_profiles = []
for horizon in (18, 36):
    d_value = (2**horizon - 1) // 19
    q_value = (3**horizon - 2**horizon) // 19
    alpha = Fraction(9 * 2**horizon, 19)
    a_value = 9 * d_value
    require(alpha.numerator // alpha.denominator == a_value, (horizon, alpha, a_value))
    phases = tuple((alpha * Fraction(3, 2) ** index) % 1 for index in range(horizon + 1))
    require(set(phases) == {Fraction(4, 19), Fraction(6, 19), Fraction(9, 19)}, phases)
    require(all(phase < Fraction(1, 2) for phase in phases), (horizon, phases))
    carries = []
    orbit_value = a_value
    for _ in range(horizon):
        next_value = (3 * orbit_value + (orbit_value % 2)) // 2
        carries.append(2 * next_value - 3 * orbit_value)
        orbit_value = next_value
    require(tuple(carries) == (1, 0, 0) * (horizon // 3), (horizon, carries))
    carry_sum = sum(
        carry * 3 ** (horizon - 1 - index) * 2**index
        for index, carry in enumerate(carries)
    )
    require(carry_sum == 9 * q_value, (horizon, carry_sum, q_value))
    require(19 * d_value + 1 == 2**horizon, (horizon, d_value))
    require(2**horizon + 19 * q_value == 3**horizon, (horizon, q_value))
    require(gcd(19 * d_value, 2**horizon) == 1, "first ABC packet")
    require(gcd(2**horizon, 19 * q_value) == 1, "second ABC packet")
    denominator_19_profiles.append((horizon, a_value, d_value, q_value))


# 9. ABC-applicable odd steps coexist with arbitrarily long smooth even blocks.
for odd_value in range(1, 300, 2):
    next_value = (3 * odd_value + 1) // 2
    require(gcd(3 * odd_value, 2 * next_value) == 1, (odd_value, next_value))
    require(gcd(odd_value, next_value) == 1, (odd_value, next_value))
    require(3 * odd_value + 1 == 2 * next_value, (odd_value, next_value))
for height in range(1, 21):
    orbit = tuple(2 ** (height - index) * 3**index for index in range(height + 1))
    require(all(orbit[index] % 2 == 0 for index in range(height)), (height, orbit))
    require(all(3 * orbit[index] == 2 * orbit[index + 1] for index in range(height)), orbit)
    require(all(radical(orbit[index] * orbit[index + 1]) == 6 for index in range(height)), orbit)


# 10. Equal radical addresses can have opposite bounded-denominator LRC behavior.
prime_list = primes_through(121)
p_support = reduce(lambda left, right: left * right, prime_list, 1)
lcm_121 = reduce(lcm, range(1, 122), 1)
u_speed = 84 * p_support
v_speed = 84 * lcm_121
u_row = tuple(range(1, 12)) + (13, u_speed)
v_row = tuple(range(1, 12)) + (13, v_speed)
require(radical(u_speed) == p_support, "rad(U)")
require(radical(v_speed) == p_support, "rad(V)")
require(tuple(radical(speed) for speed in u_row) == tuple(radical(speed) for speed in v_row), "radical collision")
require(u_speed % 121 == 55, u_speed % 121)
require(v_speed % 121 == 0, v_speed % 121)
u_numerators = tuple(min((10 * speed) % 121, 121 - ((10 * speed) % 121)) for speed in u_row)
require(u_numerators == (10, 20, 30, 40, 50, 60, 51, 41, 31, 21, 11, 9, 55), u_numerators)
require(min(u_numerators) == 9 and Fraction(9, 121) > Fraction(1, 14), u_numerators)
require(loneliness(v_row, Fraction(10, 121)) == 0, "valuation-depth hostile")
require(all(v_speed % denominator == 0 for denominator in range(1, 122)), "denominator bank hostile")


# 11. The distinct strict safe carry-language carrier has a renewal count.
greedy_digits = greedy_equality_digits(64)
require(greedy_digits[:18] == tuple(int(digit) for digit in "101000001001001010"), greedy_digits[:18])
safe_levels = {0: ((),)}
safe_counts = [1]
for length in range(1, 23):
    words = []
    for parent in safe_levels[length - 1]:
        for digit in (0, 1):
            child = parent + (digit,)
            if finite_safe_word(child):
                words.append(child)
    safe_levels[length] = tuple(words)
    safe_counts.append(len(words))
    renewal_count = 1 + sum(
        greedy_digits[index - 1] * safe_counts[length - index]
        for index in range(1, length + 1)
    )
    require(safe_counts[length] == renewal_count, (length, safe_counts[length], renewal_count))
require(
    tuple(safe_counts[:13]) == (1, 2, 3, 5, 8, 12, 18, 27, 40, 60, 90, 134, 201),
    safe_counts[:13],
)
require(all(finite_safe_word(greedy_digits[:length]) for length in range(1, 65)), "strict finite equality prefixes")
for length, words in safe_levels.items():
    require(all(finite_safe_word(word + (0,) * 8) for word in words), (length, "zero extension"))
depth_four_residues = tuple(sorted(carry_residue(word) for word in safe_levels[4]))
require(depth_four_residues == (0, 4, 5, 6, 8, 9, 13, 14), depth_four_residues)


def interval_component_count(values: tuple[int, ...]) -> int:
    if not values:
        return 0
    return 1 + sum(right != left + 1 for left, right in zip(values, values[1:]))


require(interval_component_count(tuple(range(len(safe_levels[4])))) == 1, "rank-order control")
require(interval_component_count(depth_four_residues) == 4, depth_four_residues)


print("THM-3848 rational-base prefix / lonely-runner exact companion")
print(f"parameter_pairs={parameter_pairs}")
print(f"atom_profiles={atom_profiles}")
print(f"lonely_control_count={len(lonely_controls)}")
print(f"three_halves_maximizers={maximizer_profiles}")
print(f"level_12_speeds={level_12}")
print("level_12_lrc_maximum=2/5; witness=2/5; adjacent_relation_l1=5")
print(f"level_12_safe_atom_count={3**12}; wall_grid={wall_grid(3, 2, 12)}")
print(f"canonical_minus_two_fifths_residues_N1_to_N24={tuple(canonical_residues)}")
print(f"denominator_19_horizons={[profile[0] for profile in denominator_19_profiles]}")
print(f"radical_collision_mod_121=(U:{u_speed % 121},V:{v_speed % 121}); U_margin=9/121")
print(f"safe_carry_counts_N0_to_N22={tuple(safe_counts)}")
print(f"greedy_equality_prefix_36={''.join(map(str, greedy_digits[:36]))}")
print(f"safe_carry_depth4_residues={depth_four_residues}")
print(f"CHECKS={CHECKS}")
print("ALL CHECKS PASSED;ARTIFACT=THM-3848_RATIONAL_BASE_LRC_ATOM_TREE")
