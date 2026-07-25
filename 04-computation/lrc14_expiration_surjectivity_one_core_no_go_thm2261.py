#!/usr/bin/env python3
"""Exact checks for THM-2261.

This companion verifies the one-core Markov carrier, its comparison with the
THM-2255 expiration floor, and the strict-row hostile point.  The universal
surjectivity statement is proved analytically in the theorem.
"""

from fractions import Fraction
from itertools import product


def norm_circle(x: Fraction) -> Fraction:
    """Distance of x to the nearest integer."""
    r = x % 1
    return min(r, 1 - r)


def in_danger(speed: int, x: Fraction) -> bool:
    return norm_circle(speed * x) < Fraction(1, 14)


def in_guard(speed: int, x: Fraction) -> bool:
    return norm_circle(speed * x) > Fraction(1, 7)


P = (
    (Fraction(11, 13), Fraction(2, 13)),
    (Fraction(12, 13), Fraction(1, 13)),
)
PI = (Fraction(6, 7), Fraction(1, 7))

assert all(sum(row) == 1 for row in P)
assert all(
    sum(PI[i] * P[i][j] for i in (0, 1)) == PI[j]
    for j in (0, 1)
)
assert PI[0] * P[0][1] == PI[1] * P[1][0]


def word_probability(word: tuple[int, ...]) -> Fraction:
    ans = PI[word[0]]
    for left, right in zip(word, word[1:]):
        ans *= P[left][right]
    return ans


carrier_words = []
carrier_mass = Fraction(0)
total_mass = Fraction(0)
for word in product((0, 1), repeat=5):
    mass = word_probability(word)
    total_mass += mass
    in_carrier = any(word[0:3]) and any(word[2:5])
    if in_carrier:
        carrier_words.append(word)
        carrier_mass += mass

assert total_mass == 1
assert len(carrier_words) == 25
assert carrier_mass == Fraction(6055, 28561)

inclusion_exclusion_mass = (
    1
    - 2 * PI[0] * P[0][0] ** 2
    + PI[0] * P[0][0] ** 4
)
assert inclusion_exclusion_mass == carrier_mass

expiration_floor = Fraction(88159, 415800)
crossing_gap = expiration_floor - carrier_mass
assert crossing_gap == Fraction(240199, 11875663800)
assert crossing_gap > 0

x = Fraction(79, 338)
assert in_guard(1, x)

unit_speeds = tuple(1 + 338 * i for i in range(1, 6))
assert len(set(unit_speeds)) == 5
assert all(q % 13 for q in unit_speeds)
assert all(not in_danger(q, x) for q in unit_speeds)
assert all((q * x - x).denominator == 1 for q in unit_speeds)

strict_profiles = []
for c in range(5, 20):
    for b in range(2, c):
        c1, c2, c3 = 13, 13**b, 13**c
        assert in_danger(c1, x)
        assert not in_danger(c2, x)
        assert not in_danger(c3, x)
        strict_profiles.append((1, b, c))

assert len(strict_profiles) == 150

y = (13**2 * x) % 1
assert y == Fraction(1, 2)
future_word = tuple(
    int(in_danger(1, (13**k * y) % 1))
    for k in range(12)
)
assert future_word == (0,) * 12
assert not (any(future_word[0:3]) and any(future_word[2:5]))

print("THM-2261 exact verification")
print(f"transition={P}")
print(f"stationary={PI}")
print(f"carrier_word_count={len(carrier_words)}")
print(f"carrier_mass={carrier_mass}")
print(f"expiration_floor={expiration_floor}")
print(f"crossing_gap={crossing_gap}")
print(f"witness_x={x}")
print(f"unit_speeds={unit_speeds}")
print(f"strict_profiles_checked={len(strict_profiles)}")
print(f"witness_expiration_y={y}")
print("witness_future_word=" + "".join(map(str, future_word)))
print("ALL_CHECKS_PASSED")
