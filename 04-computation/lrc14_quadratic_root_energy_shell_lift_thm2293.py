#!/usr/bin/env python3
"""Exact arithmetic audit for THM-2293.

The proof is analytic; this companion freezes every finite census,
occupancy extremum, rational floor, tail comparison, smoothing bandwidth,
word-height bound, and boundary constant. It deliberately uses explicit
raising checks so ``python`` and ``python -O`` execute the same audit.
"""

from fractions import Fraction
from math import gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def strict_profiles() -> list[tuple[int, int]]:
    return [(b, c) for c in range(5, 20) for b in range(2, c)]


def root_grid_counts() -> set[int]:
    """Counts in a translated 13-grid inside the centered arc of length 1/7."""

    endpoints = {
        (Fraction(sign, 14) - Fraction(r, 13)) % 1
        for sign in (-1, 1)
        for r in range(13)
    }
    ordered = sorted(endpoints)
    cells: list[Fraction] = []
    for index, left in enumerate(ordered):
        right = ordered[(index + 1) % len(ordered)]
        if index + 1 == len(ordered):
            right += 1
        cells.append(((left + right) / 2) % 1)

    counts: set[int] = set()
    for phase in cells:
        count = 0
        for r in range(13):
            value = (phase + Fraction(r, 13)) % 1
            distance = min(value, 1 - value)
            if distance < Fraction(1, 14):
                count += 1
        counts.add(count)
    return counts


profiles = strict_profiles()
interior = [(b, c) for b, c in profiles if b >= 3 and c >= b + 2]
first_b_two = [(b, c) for b, c in profiles if b == 2]
adjacent = [(b, c) for b, c in profiles if c == b + 1]

require(len(profiles) == 150, "strict profile census")
require(len(interior) == 120, "interior profile census")
require(len(first_b_two) == 15, "b=2 boundary census")
require(len(adjacent) == 15, "adjacent boundary census")
require(sorted(interior + first_b_two + adjacent) == sorted(profiles),
        "strict branch partition")
require({c - b - 1 for b, c in profiles} == set(range(17)),
        "deepest-successor depth range")

require(root_grid_counts() == {1, 2}, "translated 13-grid occupancy")

proper_energy = {n: 13 * n - n * n for n in range(13)}
require(min(proper_energy[n] for n in range(1, 13)) == 12,
        "proper-mask energy minimum")
require(max(proper_energy[n] for n in range(13)) == 42,
        "first-label energy cap")
require(max(proper_energy[n] for n in range(3)) == 22,
        "second-label energy cap")

# The inscribed-dodecagon proof of pi > 31/10.
require(Fraction(1733, 1000) ** 2 > 3, "sqrt(3) rational upper bound")
require(36 * (Fraction(2) - Fraction(1733, 1000)) > Fraction(961, 100),
        "dodecagon squared margin")

Y0 = Fraction(5696989, 76962600)
ratio_a = 12 * Y0 / 64
tail_ratio_a = Fraction(33800, 961 * 2535)
margin_a = ratio_a - tail_ratio_a

require(ratio_a == Fraction(5696989, 410467200),
        "small-mode energy/cap ratio")
require(tail_ratio_a == Fraction(40, 2883),
        "small-mode tail ratio")
require(margin_a == Fraction(1910429, 394458979200),
        "small-mode exact margin")
require(margin_a > 0, "small-mode tail clears")
require(2534 % 7 == 0, "2534 seven boundary")
require(2535 % 13 == 0, "2535 thirteen boundary")

L0 = Fraction(5696989, 367580070)
a1 = Fraction(1183, 72)
a2 = Fraction(13013, 12)
require(a2 == 66 * a1, "weighted-root coefficient ratio")
effective = a1 * a2 / (a1 + a2)
require(effective == Fraction(13013, 804), "quadratic allocation constant")

nu0 = effective * L0 * L0
require(nu0 == Fraction(227189785662847, 58436012221844400),
        "weighted source energy floor")

ratio_b = nu0 / 64
tail_ratio_b = Fraction(33800, 961 * 578982)
margin_b = ratio_b - tail_ratio_b
require(
    margin_b
    == Fraction(
        296921577147599,
        346814897688821607884467200,
    ),
    "actual-source tail margin",
)
require(margin_b > 0, "actual-source tail clears")

# The finite Jackson realization keeps a uniform nonzero covariance head.
smoothing_tail_energy = Fraction(200, 961 * 578982)
smoothing_head_gap = Fraction(22, 26) * (
    ratio_b - 169 * smoothing_tail_energy
)
require(
    smoothing_tail_energy == Fraction(100, 278200851),
    "smoothing tail energy",
)
require(
    smoothing_head_gap
    == Fraction(
        296921577147599,
        409872151814061900227097600,
    ),
    "smoothing head gap",
)
require(smoothing_head_gap > 0, "smoothing head remains positive")

# From ||V_N-V||_2 < 78*sqrt(63/N), it is enough to make this error
# strictly smaller than half the preceding exact head gap.
smoothing_threshold = (
    Fraction(4 * 78**2 * 63) / smoothing_head_gap**2
)
smoothing_n = smoothing_threshold.numerator // smoothing_threshold.denominator + 1
smoothing_h = 2 * smoothing_n - 2
smoothing_word_difference_height = 2 * smoothing_h
smoothing_relation_height = smoothing_word_difference_height + 578982
smoothing_squared_margin = (
    smoothing_head_gap**2 / 4
    - Fraction(78**2 * 63, smoothing_n)
)

require(
    smoothing_n == 2921480906639115573490032947784,
    "smoothing bandwidth",
)
require(
    Fraction(smoothing_n - 1) <= smoothing_threshold < Fraction(smoothing_n),
    "smoothing certificate minimality",
)
require(
    smoothing_h == 5842961813278231146980065895566,
    "smoothing polynomial degree",
)
require(
    smoothing_word_difference_height
    == 11685923626556462293960131791132,
    "smoothing word-difference height",
)
require(
    smoothing_relation_height
    == 11685923626556462293960132370114,
    "smoothing relation-word height",
)
require(
    smoothing_squared_margin
    == Fraction(
        1512520856980458045994136291,
        81799118868347936026371657539647828685067698505281323958615849827440899472752640000,
    ),
    "smoothing squared margin",
)
require(smoothing_squared_margin > 0, "smoothing error clears head gap")

# Exact shell arithmetic for every strict profile.  Test every possible
# nonzero residue multiplier modulo 13 and every allowed profile.
for b, c in profiles:
    d = c - b - 1
    require(b + 1 + d == c, "shell exponent identity")
    for m in range(1, 13):
        if m % 13 == 0:
            continue
        # v_13(m * 13^c) is exactly c.
        value = m * 13**c
        quotient = value
        valuation = 0
        while quotient % 13 == 0:
            quotient //= 13
            valuation += 1
        require(valuation == c, "terminal exact valuation")

require(gcd(2533, 91) == 1, "small bound endpoint is a unit")
require(gcd(578982, 91) == 1, "source-pair bound endpoint is a unit")

delta3_bound = 54991358114
delta4_bound = 2533 * delta3_bound
require(delta4_bound == 139293110102762, "augmented determinant bound")

print("THM-2293 exact arithmetic audit")
print(f"strict_profiles={len(profiles)}")
print(f"interior_profiles={len(interior)}")
print(f"b2_boundary_profiles={len(first_b_two)}")
print(f"adjacent_boundary_profiles={len(adjacent)}")
print("deepest_successor_depths=0..16")
print("translated_13_grid_counts=1,2")
print("proper_root_energy_min=12")
print("unweighted_caps=42,22")
print("pi_lower=31/10")
print(f"small_mode_ratio={ratio_a}")
print(f"small_mode_tail_ratio={tail_ratio_a}")
print(f"small_mode_margin={margin_a}")
print("small_mode_multiplier_bound=2533")
print("excluded_boundary=2534:7,2535:13")
print(f"weighted_coefficients={a1},{a2}")
print(f"weighted_effective_coefficient={effective}")
print(f"weighted_source_energy_floor={nu0}")
print(f"source_pair_ratio={ratio_b}")
print(f"source_pair_tail_ratio={tail_ratio_b}")
print(f"source_pair_margin={margin_b}")
print("source_pair_multiplier_bound=578982")
print("fixed_character_multiplier_bound=UNIFORM_NOT_PROVED")
print(f"smoothing_tail_energy={smoothing_tail_energy}")
print(f"smoothing_head_gap={smoothing_head_gap}")
print(f"smoothing_N={smoothing_n}")
print(f"smoothing_H={smoothing_h}")
print(f"smoothing_word_difference_height={smoothing_word_difference_height}")
print(f"smoothing_relation_height={smoothing_relation_height}")
print(f"smoothing_squared_margin={smoothing_squared_margin}")
print("smoothing_genuine_relation=NOT_FORCED")
print(f"augmented_determinant_bound={delta4_bound}")
print("status=PASS")
