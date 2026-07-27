#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2577.

Builds the canonical typed-row endpoint sets on their exact rational grid and
checks the support-image mechanism behind the all-clock collision law.  The
packet clock never enters: every nonzero returned subpacket is supported in
its word set, so the displayed image separations/inclusions decide its first
owner-normalized collision uniformly.
"""

from fractions import Fraction
from math import gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


W = (1, 14, 27, 40, 53, 66, 13, 13**3, 2 * 13**5)
LCM_W = 1
for speed in W:
    LCM_W = LCM_W * speed // gcd(LCM_W, speed)
T_DEN = 182 * LCM_W
require(T_DEN == 297836897838480, "canonical grid changed")


def in_comb(index):
    """Intervals of {x: ||W[index] x|| < 1/14} on the exact grid."""
    speed = W[index]
    unit = T_DEN // (182 * speed)
    out = []
    for tooth in range(speed):
        start = (-13 + 182 * tooth) * unit % T_DEN
        stop = start + 26 * unit
        if stop <= T_DEN:
            out.append((start, stop))
        else:
            out.append((start, T_DEN))
            out.append((0, stop - T_DEN))
    return sorted(out)


def subtract_comb(intervals, speed, period, lo, hi):
    """Subtract all periodic windows [lo+period*n,hi+period*n)."""
    unit = T_DEN // (period * speed)
    step = period * unit
    width = (hi - lo) * unit
    base = (lo % period) * unit
    out = []
    for left, right in intervals:
        first = left - width + 1
        k = -((base - first) // step)
        window = base + k * step
        cursor = left
        while window < right:
            end = window + width
            if end > cursor:
                if window > cursor:
                    out.append((cursor, window))
                cursor = end
                if cursor >= right:
                    break
            window += step
        if cursor < right:
            out.append((cursor, right))
    return out


def build_set(pattern):
    """Build one Boolean comb atom; labels are in/out/gout."""
    inside = [index for index, mode in pattern.items() if mode == "in"]
    if inside:
        start = min(inside, key=lambda index: W[index])
        intervals = in_comb(start)
    else:
        start = None
        intervals = [(0, T_DEN)]

    for index, mode in pattern.items():
        if mode == "gout":
            intervals = subtract_comb(
                intervals, W[index], 91, -13, 13
            )

    remaining = sorted(
        (W[index], index)
        for index, mode in pattern.items()
        if index != start and mode in ("in", "out")
    )
    for _, index in remaining:
        if pattern[index] == "out":
            intervals = subtract_comb(
                intervals, W[index], 182, -13, 13
            )
        else:
            intervals = subtract_comb(
                intervals, W[index], 182, 13, 169
            )
    return intervals


def measure(intervals):
    return sum(right - left for left, right in intervals)


def check_intervals(intervals):
    previous = -1
    for left, right in intervals:
        require(
            0 <= left < right <= T_DEN and left >= previous,
            "interval list is not sorted/disjoint",
        )
        previous = right
    return measure(intervals)


def intersect(left_intervals, right_intervals):
    out = []
    i = j = 0
    while i < len(left_intervals) and j < len(right_intervals):
        left = max(left_intervals[i][0], right_intervals[j][0])
        right = min(left_intervals[i][1], right_intervals[j][1])
        if left < right:
            out.append((left, right))
        if left_intervals[i][1] < right_intervals[j][1]:
            i += 1
        else:
            j += 1
    return out


def image_support(intervals, multiplier):
    """Essential support of P_multiplier(1_intervals)."""
    pieces = []
    for left, right in intervals:
        width = multiplier * (right - left)
        if width >= T_DEN:
            return [(0, T_DEN)]
        start = multiplier * left % T_DEN
        stop = start + width
        if stop <= T_DEN:
            pieces.append((start, stop))
        else:
            pieces.append((start, T_DEN))
            pieces.append((0, stop - T_DEN))

    merged = []
    for left, right in sorted(pieces):
        if merged and left <= merged[-1][1]:
            if right > merged[-1][1]:
                merged[-1] = (merged[-1][0], right)
        else:
            merged.append((left, right))
    return merged


COMMON = {
    0: "gout",
    1: "out",
    2: "out",
    3: "out",
    4: "out",
    5: "out",
}
PATTERNS = {
    "E": {**COMMON, 6: "in", 7: "out", 8: "out"},
    "a": {**COMMON, 6: "out", 7: "in", 8: "out"},
    "b": {**COMMON, 6: "out", 7: "out", 8: "in"},
    "ab": {**COMMON, 6: "out", 7: "in", 8: "in"},
}


print("== THM-2577: all-clock word-depth support law ==")
sets = {name: build_set(pattern) for name, pattern in PATTERNS.items()}
for intervals in sets.values():
    check_intervals(intervals)

expected_anchors = {
    "E": (57072, Fraction(1882176, 28589561)),
    "a": (22478, Fraction(143103830843, 5727632650740)),
    "b": (131762, Fraction(4839079319, 190921088358)),
    "ab": (22018, Fraction(11009, 2599051)),
}
print("\n== exact typed-row atoms ==")
for name in ("E", "a", "b", "ab"):
    observed = (len(sets[name]), Fraction(measure(sets[name]), T_DEN))
    require(observed == expected_anchors[name], f"{name} anchor changed")
    print(f"  {name:>2}: intervals={observed[0]}, mass={observed[1]}")


specifications = {
    "a": 3,
    "b": 5,
    "ab": 5,
}
display_word = {"a": "a", "b": "b", "ab": "a,b"}
print("\n== exact support images ==")
for word, depth in specifications.items():
    word_set = sets[word]
    zero_depths = []
    for exponent in range(1, depth + 1):
        owner_image = image_support(sets["E"], 13**exponent)
        word_image = image_support(word_set, 13**exponent)
        overlap = measure(intersect(owner_image, word_image))
        require(overlap == 0, f"{word} overlaps before depth {depth}")
        zero_depths.append(exponent)

    exponent = depth + 1
    owner_image = image_support(sets["E"], 13**exponent)
    word_image = image_support(word_set, 13**exponent)
    overlap = measure(intersect(owner_image, word_image))
    require(
        overlap == measure(word_image),
        f"{word} next image is not contained in the owner image",
    )

    if word == "a":
        require(owner_image == word_image, "a-word terminal images differ")
        require(
            (len(owner_image), Fraction(measure(owner_image), T_DEN))
            == (26, Fraction(6, 7)),
            "a-word terminal image anchor changed",
        )
        terminal = "equal 26-interval images of mass 6/7"
    else:
        require(
            owner_image == [(0, T_DEN)]
            and word_image == [(0, T_DEN)],
            f"{word} terminal image is not the full circle",
        )
        terminal = "both images are the full circle"

    print(
        f"  word={{{display_word[word]}}}: disjoint at exponents {zero_depths}; "
        f"exponent {exponent}: {terminal}"
    )


print("\n== arbitrary returned-subpacket controls ==")
for word, depth in specifications.items():
    word_set = sets[word]
    indices = (0, len(word_set) // 2, len(word_set) - 1)
    for index in indices:
        left, right = word_set[index]
        midpoint = (left + right) // 2
        require(midpoint > left, "control interval too short")
        packet = [(left, midpoint)]
        require(measure(packet) > 0, "control packet has zero mass")
        before = image_support(packet, 13**depth)
        owner_before = image_support(sets["E"], 13**depth)
        require(
            measure(intersect(before, owner_before)) == 0,
            "positive subpacket collided too early",
        )
        after = image_support(packet, 13 ** (depth + 1))
        owner_after = image_support(sets["E"], 13 ** (depth + 1))
        require(
            measure(after) > 0
            and measure(intersect(after, owner_after)) == measure(after),
            "positive subpacket failed terminal inclusion",
        )
    print(
        f"  word={{{display_word[word]}}}: "
        "first/middle/last positive subintervals pass"
    )


print("\nabstract consequence:")
print("  f_K = 1_Q P^(K)1_E nonzero => supp(f_K) subset Q")
print("  image separation gives I_s=0 before the word depth")
print("  terminal image inclusion gives I_depth>0 for every such K")
print("  therefore r({a},K)=3 and r({b},K)=r({a,b},K)=5")
print("\ntyped-row scope only; no physical current or LRC(14) row exclusion")
print("\nall exact checks passed")
