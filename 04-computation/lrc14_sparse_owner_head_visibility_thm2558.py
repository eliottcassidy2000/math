#!/usr/bin/env python3
"""Exact finite referee for THM-2558.

The referee classifies the image of the THM-2531 lexicographic selected head
when all twelve nonzero slopes are used on a nonconstant 13-bit mask.  It also
checks the sharp middle-layer blind control and the elementary unit/guard
root-capacity bounds used at the LRC interface.
"""

from collections import Counter
from fractions import Fraction


P = 13


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def support(code):
    return {r for r in range(P) if (code >> r) & 1}


def cyclic_word(mask, tau, a):
    return tuple(int((a + j * tau) % P in mask) for j in range(P))


def selected_head(mask, tau):
    words = [(cyclic_word(mask, tau, a), a) for a in range(P)]
    alpha = max(words)[1]
    q = next(j for j in range(1, P) if (alpha + j * tau) % P not in mask)
    source = (alpha + (q - 1) * tau) % P
    head = (alpha + q * tau) % P
    require(source in mask and head not in mask, "selector is not a boundary")
    return head


def head_image(mask):
    return {selected_head(mask, tau) for tau in range(1, P)}


expected_histograms = {
    1: {12: 13},
    2: {11: 78},
    3: {6: 52, 8: 156, 10: 78},
    4: {6: 52, 7: 390, 8: 156, 9: 117},
    5: {6: 468, 7: 468, 8: 351},
    6: {5: 78, 6: 650, 7: 988},
    7: {5: 156, 6: 1560},
    8: {5: 1287},
    9: {4: 715},
    10: {3: 286},
    11: {2: 78},
    12: {1: 13},
}

histograms = {k: Counter() for k in range(1, P)}
good_by_weight = Counter()
good_codes = set()

for code in range(1, 2**P - 1):
    mask = support(code)
    heads = head_image(mask)
    require(heads <= set(range(P)) - mask, "a selected head is occupied")
    histograms[len(mask)][len(heads)] += 1
    if heads == set(range(P)) - mask:
        good_codes.add(code)
        good_by_weight[len(mask)] += 1

require(
    {k: dict(histograms[k]) for k in histograms} == expected_histograms,
    "head-diversity census changed",
)
require(len(good_codes) == 5564, "wrong all-empty-visible labelled count")
require((2**P - 2) - len(good_codes) == 2626, "wrong blind labelled count")

# Rotation does not change visibility; every nonconstant orbit has size 13.
def rotate_code(code, b):
    return sum(1 << ((r + b) % P) for r in support(code))


for code in range(1, 2**P - 1):
    require(
        all((rotate_code(code, b) in good_codes) == (code in good_codes) for b in range(P)),
        "visibility is not rotation invariant",
    )

rotation_representatives = {
    min(rotate_code(code, b) for b in range(P)) for code in range(1, 2**P - 1)
}
good_necklaces = sum(rep in good_codes for rep in rotation_representatives)
require(len(rotation_representatives) == 630, "wrong necklace count")
require(good_necklaces == 428, "wrong visible-necklace count")
require(len(rotation_representatives) - good_necklaces == 202, "wrong blind-necklace count")


print("== THM-2558: all-slope selected-head visibility on F_13 ==")
print("  weight : head-diversity histogram ; all-empty-visible count")
for k in range(1, P):
    hist = ", ".join(f"{h}:{histograms[k][h]}" for h in sorted(histograms[k]))
    print(f"  {k:>2} : {hist} ; {good_by_weight[k]}")
print(f"  labelled masks: visible {len(good_codes)}, blind {2626}")
print(f"  rotation necklaces: visible {good_necklaces}, blind {202}")


print("\n== sharp singleton and pair laws ==")
single = {0}
single_heads = [selected_head(single, tau) for tau in range(1, P)]
require(set(single_heads) == set(range(1, P)), "singleton misses an empty root")
require(len(set(single_heads)) == len(single_heads), "singleton repeats a head")

pair = {0, 1}
pair_heads = [selected_head(pair, tau) for tau in range(1, P)]
pair_midpoint = 7  # 2^{-1} modulo 13
require(set(pair_heads) == set(range(2, P)), "pair misses an empty root")
pair_counts = Counter(pair_heads)
require(pair_counts[pair_midpoint] == 2, "pair midpoint multiplicity changed")
require(all(pair_counts[t] == 1 for t in range(2, P) if t != pair_midpoint), "bad pair multiplicity")
print("  singleton: all 12 empty roots, once each")
print("  pair {0,1}: all 11 empty roots; midpoint 7 twice, all others once")
print("  any neutral union of size <=10 misses a singleton/pair selected head")


print("\n== exact middle-layer blind control ==")
hostile = {0, 1, 4}
hostile_heads = head_image(hostile)
hostile_blind = (set(range(P)) - hostile) - hostile_heads
target_active_failure = {3}
require(hostile_heads == {2, 7, 8, 9, 11, 12}, "hostile head image changed")
require(hostile_blind == {3, 5, 6, 10}, "hostile blind set changed")
require(target_active_failure <= hostile_blind, "target failure is not hidden")
require(not (target_active_failure & hostile_heads), "hostile has a target-active hit")
print("  mask {0,1,4}: selected heads {2,7,8,9,11,12}")
print("  blind empty roots: {3,5,6,10}")
print("  singleton target-active failure {3}: missed by all 12 slopes")


print("\n== unit and guard root-capacity controls ==")


def circle_distance(x):
    y = x - (x.numerator // x.denominator)
    return min(y, 1 - y)


unit_counts = set()
guard_counts = set()
for w in range(1, P):
    for numerator in range(1, 37):
        z = Fraction(numerator, 37)
        unit = sum(circle_distance(Fraction(w, P) * (z + r)) < Fraction(1, 14) for r in range(P))
        guard = sum(circle_distance(Fraction(w, P) * (z + r)) < Fraction(1, 7) for r in range(P))
        require(unit in {1, 2}, "unit danger capacity is not one or two")
        require(guard in {3, 4}, "guard failure capacity is not three or four")
        unit_counts.add(unit)
        guard_counts.add(guard)
require(unit_counts == {1, 2} and guard_counts == {3, 4}, "sharp capacities not attained")
print("  ordinary unit danger count: 1 or 2 (both attained)")
print("  guard failure count: 3 or 4 (both attained)")
print("  four neutral roles cover at most 4+2+2+2=10 roots")


print("\nall checks passed")
