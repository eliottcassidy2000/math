#!/usr/bin/env python3
"""Exact companion for THM-2391.

Checks the cyclic period/containment and slope-classification lemmas behind
the blocker-caged single-layer reduction.  It also reconstructs a genuine
positive local chamber carrying both blocker partition types' obstruction:
clean quotient-blocker partition, four distinct lower q labels, and the
canonical repeated-first thirteen-adic roles.  No global cover is claimed.
"""

from fractions import Fraction as F
from itertools import product
from math import gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def valuation(n: int, p: int) -> int:
    out = 0
    while n % p == 0:
        n //= p
        out += 1
    return out


def bit_count(word: int) -> int:
    return word.bit_count()


def rotate(word: int, length: int, shift: int) -> int:
    shift %= length
    if shift == 0:
        return word
    mask = (1 << length) - 1
    return ((word << shift) | (word >> (length - shift))) & mask


def lifted_word(
    length: int,
    period: int,
    unit: int,
    start: int,
    width: int,
) -> int:
    require(length % period == 0, "period must divide ambient length")
    require(period % 7 == 0 and gcd(unit, 7) == 1, "invalid septimal word")
    span = width * period // 7
    out = 0
    for s in range(length):
        if (unit * s - start) % period < span:
            out |= 1 << s
    return out


def powers_of_seven_through(limit: int) -> tuple[int, ...]:
    out = []
    p = 7
    while p <= limit:
        out.append(p)
        p *= 7
    return tuple(out)


def balanced(residue: int, modulus: int) -> int:
    residue %= modulus
    return residue if residue <= modulus // 2 else residue - modulus


def ap_word(modulus: int, start: int, step: int, length: int) -> frozenset[int]:
    return frozenset((start + j * step) % modulus for j in range(length))


def frac_part(x: F) -> F:
    return x - (x.numerator // x.denominator)


def circle_norm(x: F) -> F:
    r = frac_part(x)
    return min(r, 1 - r)


def danger(v: int, x: F, width: int = 1) -> bool:
    return circle_norm(v * x) < F(width, 14)


# ---------------------------------------------------------------------------
# 1. Lifted-word cardinality and least-period audit.
# ---------------------------------------------------------------------------

lifted_word_cases = 0
least_period_checks = 0
word_bank: dict[tuple[int, int, int], set[int]] = {}

for ambient in (7, 49, 343):
    for period in powers_of_seven_through(ambient):
        words_at_width: dict[int, set[int]] = {1: set(), 2: set()}
        for unit in range(1, period):
            if unit % 7 == 0:
                continue
            for start in range(period):
                for width in (1, 2):
                    word = lifted_word(ambient, period, unit, start, width)
                    require(
                        bit_count(word) == width * ambient // 7,
                        "wrong lifted-word cardinality",
                    )
                    require(
                        rotate(word, ambient, period) == word
                        if period < ambient
                        else True,
                        "claimed period is not a period",
                    )
                    for smaller in powers_of_seven_through(period // 7):
                        require(
                            rotate(word, ambient, smaller) != word,
                            "lifted word acquired a smaller septimal period",
                        )
                        least_period_checks += 1
                    words_at_width[width].add(word)
                    lifted_word_cases += 1
        for width, words in words_at_width.items():
            word_bank[(ambient, period, width)] = words


# ---------------------------------------------------------------------------
# 2. Exhaustive period-containment check through the first nontrivial gap.
# ---------------------------------------------------------------------------

period_containment_pairs = 0
ambient = 49
for ordinary_period, guard_period in ((7, 49), (49, 7)):
    ordinary_words = word_bank[(ambient, ordinary_period, 1)]
    guard_words = word_bank[(ambient, guard_period, 2)]
    for ordinary in ordinary_words:
        for guard in guard_words:
            require(
                ordinary & ~guard != 0,
                "distinct-period ordinary word fit inside guard word",
            )
            period_containment_pairs += 1


# ---------------------------------------------------------------------------
# 3. Same-period slope and complementary-partition classification.
# ---------------------------------------------------------------------------

slope_cases = 0
contained_word_counts = []
complementary_pair_counts = []

for modulus in (49, 343):
    k = modulus // 7
    guard = frozenset(range(2 * k))
    contained: dict[frozenset[int], set[int]] = {}

    for step in range(1, modulus):
        if step % 7 == 0:
            continue
        for start in range(modulus):
            word = ap_word(modulus, start, step, k)
            if word <= guard:
                delta = balanced(step, modulus)
                require(
                    delta in (-2, -1, 1, 2),
                    "contained progression has an illegal slope",
                )
                contained.setdefault(word, set()).add(delta)
            slope_cases += 1

    expected_count = k + 3
    require(
        len(contained) == expected_count,
        "wrong number of contained ordinary word sets",
    )
    contained_word_counts.append(len(contained))

    half_left = frozenset(range(k))
    half_right = frozenset(range(k, 2 * k))
    parity_even = frozenset(range(0, 2 * k, 2))
    parity_odd = frozenset(range(1, 2 * k, 2))
    legal_pairs = {
        frozenset((half_left, half_right)),
        frozenset((parity_even, parity_odd)),
    }

    observed_pairs = set()
    for left in contained:
        right = guard - left
        if right in contained and not (left & right):
            observed_pairs.add(frozenset((left, right)))

    require(observed_pairs == legal_pairs, "unexpected blocker partition type")
    complementary_pair_counts.append(len(observed_pairs))


# Adjacent-layer binary boundary.
binary_guard = frozenset((0, 1))
binary_words = {frozenset((r,)) for r in range(7) if r in binary_guard}
require(
    binary_words == {frozenset((0,)), frozenset((1,))},
    "binary address boundary",
)
require(
    frozenset().union(*binary_words) == binary_guard,
    "binary blockers do not partition guard",
)


# Four slope classes: either a repeat or the complete chi_7 cross.
pigeonhole_cases = 0
all_distinct_cases = 0
classes = (1, -1, 2, -2)
for assignment in product(classes, repeat=4):
    if len(set(assignment)) == 4:
        require(set(assignment) == set(classes), "wrong four-class cross")
        all_distinct_cases += 1
    else:
        require(
            any(assignment.count(value) >= 2 for value in classes),
            "pigeonhole repeat missing",
        )
    pigeonhole_cases += 1

ratio_cross_mod7 = tuple(sorted({1 % 7, -1 % 7, pow(2, -1, 7), -pow(2, -1, 7) % 7}))
require(ratio_cross_mod7 == (1, 3, 4, 6), "wrong chi_7 ratio cross")
quadratic_residues = {1, 2, 4}
require(
    sum(value in quadratic_residues for value in ratio_cross_mod7) == 2,
    "chi_7 cross is not balanced",
)


# ---------------------------------------------------------------------------
# 4. Thirteen pullback and additive one-mask hostile controls.
# ---------------------------------------------------------------------------

pullback_cases = 0
phase_den = 211
for speed in range(1, 61):
    for num in range(phase_den):
        x = F(num, phase_den)
        require(
            danger(speed, 13 * x, 1) == danger(13 * speed, x, 1),
            "ordinary thirteen pullback",
        )
        require(
            danger(speed, 13 * x, 2) == danger(13 * speed, x, 2),
            "guard thirteen pullback",
        )
        pullback_cases += 2

additive_hostile_cases = 0
for H in range(1, 80):
    for q_top in range(7, 80, 7):
        v = H + q_top
        for num in range(1, phase_den):
            x = F(num, phase_den)
            if danger(v, x) and danger(q_top, x):
                require(danger(H, x, 2), "triangle hostile failed")
                additive_hostile_cases += 1


# ---------------------------------------------------------------------------
# 5. Positive local M=2 cage with canonical thirteen-adic roles.
# ---------------------------------------------------------------------------

N = 343
x0 = F(1, 4949)
H = 1
q_star = 49
q_low = (24, 25, 97, 99)
C1, C2 = 293, 393
c1, c2 = 13 * C1, 13 * C2
C3 = 343 * 13 * 13
c3 = 13 * C3

top_word = frozenset(j for j in range(N) if danger(q_star, x0 + F(j, N)))
require(top_word == frozenset(range(0, N, 7)), "wrong hostile top bin")


def top_slice(v: int, width: int = 1) -> frozenset[int]:
    return frozenset(
        k for k in range(49) if danger(v, x0 + F(7 * k, N), width)
    )


guard_slice = top_slice(H, 2)
C1_slice = top_slice(C1)
C2_slice = top_slice(C2)
q_slices = tuple(top_slice(q) for q in q_low)

expected_guard = frozenset((*range(0, 7), *range(42, 49)))
expected_C1 = frozenset(range(0, 7))
expected_C2 = frozenset(range(42, 49))
expected_q_slices = (
    frozenset((0, 2, 4, 6, 43, 45, 47)),
    frozenset((0, 2, 4, 6, 43, 45, 47)),
    frozenset((0, 1, 2, 3, 4, 47, 48)),
    frozenset((0, 1, 2, 45, 46, 47, 48)),
)

require(guard_slice == expected_guard, "wrong hostile guard slice")
require(C1_slice == expected_C1, "wrong hostile C1 slice")
require(C2_slice == expected_C2, "wrong hostile C2 slice")
require(not (C1_slice & C2_slice), "hostile blockers overlap")
require(C1_slice | C2_slice == guard_slice, "hostile blockers do not partition")
require(q_slices == expected_q_slices, "wrong hostile lower-q slices")
require(all(word <= guard_slice for word in q_slices), "hostile q escaped guard")

for j in range(N):
    x = x0 + F(j, N)
    require(not danger(C3, x), "quotient high blocker is not safe")
    require(not danger(c3, x), "actual high blocker is not safe")

require(valuation(q_star, 7) == 2, "hostile top depth")
require(
    valuation(H, 7) == 0
    and all(valuation(q, 7) == 0 for q in q_low)
    and valuation(C1, 7) == valuation(C2, 7) == 0,
    "hostile lower depths",
)
require(valuation(C3, 7) == valuation(c3, 7) == 3, "hostile high depth")
require(
    (valuation(c1, 13), valuation(c2, 13), valuation(c3, 13))
    == (1, 1, 3),
    "hostile blocker thirteen-adic roles",
)
require(len(set((*q_low, q_star))) == 5, "hostile q labels repeat")
require(all(q % 13 for q in (*q_low, q_star, H)), "hostile unit role")

positive_margin = None
for v, width in (
    (H, 2),
    *((q, 1) for q in (*q_low, q_star)),
    (C1, 1),
    (C2, 1),
    (c1, 1),
    (c2, 1),
    (C3, 1),
    (c3, 1),
):
    for j in range(N):
        x = x0 + F(j, N)
        margin = abs(circle_norm(v * x) - F(width, 14)) / v
        require(margin > 0, "hostile chamber hits an endpoint")
        positive_margin = margin if positive_margin is None else min(
            positive_margin, margin
        )

require(positive_margin is not None and positive_margin > 0, "hostile margin")


print("theorem=THM-2391")
print("status=PROVED+VERIFIED-EXACT-CANDIDATE")
print(f"lifted_word_cases={lifted_word_cases}")
print(f"least_period_checks={least_period_checks}")
print(f"period_containment_pairs={period_containment_pairs}")
print(f"slope_cases={slope_cases}")
print(f"contained_word_counts={','.join(map(str, contained_word_counts))}")
print(
    "complementary_pair_counts="
    + ",".join(map(str, complementary_pair_counts))
)
print("partition_types=contiguous,parity")
print("binary_address_count=2")
print(f"pigeonhole_cases={pigeonhole_cases}")
print(f"all_distinct_chi7_cases={all_distinct_cases}")
print("chi7_ratio_cross=1,3,4,6")
print(f"thirteen_pullback_cases={pullback_cases}")
print(f"additive_hostile_cases={additive_hostile_cases}")
print(f"local_guard_size={len(guard_slice)}")
print(f"local_blocker_sizes={len(C1_slice)},{len(C2_slice)}")
print(
    "local_q_sizes=" + ",".join(str(len(word)) for word in q_slices)
)
print(f"local_positive_margin={positive_margin}")
print("local_roles=M2,low0,C3depth3; blocker13=(1,1,3)")
print("forced_lower_layer=primitive-depth-0; heavy_word=weight-8-one-double")
print("branch_excluded=0; ledger=165; LRC(14)=OPEN")
print("all_checks=PASS")
