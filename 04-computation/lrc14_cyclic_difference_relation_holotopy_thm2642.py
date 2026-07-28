#!/usr/bin/env python3
"""Exact referee for THM-2642.

All logical checks use ``require`` so optimized Python executes the same
certificate as ordinary Python.
"""

from collections import Counter
from itertools import combinations, product


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def is_prime(p):
    if p < 2:
        return False
    return all(p % d for d in range(2, int(p ** 0.5) + 1))


def members(mask, p):
    return tuple(i for i in range(p) if (mask >> i) & 1)


def convolution_count(left, right, c, p):
    right_set = set(right)
    return sum(1 for x in left if (c - x) % p in right_set)


def sumset(left, right, p):
    return {(x + y) % p for x in left for y in right}


def cycle_convolution(sets, p):
    counts = [0] * p
    counts[0] = 1
    for allowed in sets:
        nxt = [0] * p
        for partial, value in enumerate(counts):
            for step in allowed:
                nxt[(partial + step) % p] += value
        counts = nxt
    return tuple(counts)


def direct_two_edge_sections(left, right, c, p):
    total = 0
    for x0 in range(p):
        for s in left:
            x1 = (x0 + s) % p
            for t in right:
                x2 = (x1 + t) % p
                total += x2 == (x0 + c) % p
    return total


def exact_extremizer(p, a, b, c=0):
    """Sets attaining max(0,a+b-p) representations of c."""
    left = set(range(a))
    complement = [x for x in range(p) if x not in left]
    overlap = max(0, a + b - p)
    reflected_right = set(complement[: b - overlap])
    reflected_right.update(range(overlap))
    require(len(reflected_right) == b, "extremizer reflected size")
    right = {(-x + c) % p for x in reflected_right}
    return left, right


def main():
    p = 13
    require(is_prime(p), "p=13 must be prime")

    # Exhaust every nonempty two-edge relation on F_7.  This simultaneously
    # checks Cauchy--Davenport and the exact section/convolution identity.
    p_small = 7
    require(is_prime(p_small), "p=7 must be prime")
    subsets = [members(mask, p_small) for mask in range(1, 1 << p_small)]
    pair_count = 0
    for left in subsets:
        for right in subsets:
            pair_count += 1
            support = sumset(left, right, p_small)
            lower_support = min(p_small, len(left) + len(right) - 1)
            require(len(support) >= lower_support, "Cauchy-Davenport failure")
            for c in range(p_small):
                reps = convolution_count(left, right, c, p_small)
                sections = direct_two_edge_sections(left, right, c, p_small)
                require(sections == p_small * reps, "section convolution failure")
                require(reps >= max(0, len(left) + len(right) - p_small),
                        "two-edge representation floor failure")

    # The two-edge lower bound is sharp for every nonempty size profile.
    extremizer_profiles = 0
    for a in range(1, p + 1):
        for b in range(1, p + 1):
            left, right = exact_extremizer(p, a, b)
            reps = convolution_count(left, right, 0, p)
            require(reps == max(0, a + b - p), "sharp extremizer failure")
            extremizer_profiles += 1

    # One below the iterated saturation threshold can genuinely miss a clutch.
    near_left = tuple(range(6))
    near_right = tuple(range(7))
    near_support = sumset(near_left, near_right, p)
    require(len(near_support) == 12 and 12 not in near_support,
            "near-threshold hostile failure")
    require((len(near_left) - 1) + (len(near_right) - 1) == p - 2,
            "near-threshold excess")

    # At the threshold every clutch survives.
    threshold_left = tuple(range(6))
    threshold_right = tuple(range(8))
    threshold_support = sumset(threshold_left, threshold_right, p)
    require(len(threshold_support) == p, "threshold saturation failure")
    require((len(threshold_left) - 1) + (len(threshold_right) - 1) == p - 1,
            "threshold excess")

    # Exhaust all C(13,11)^2 pairs of eleven-sheet relations and all carries.
    eleven_sets = []
    for missing in combinations(range(p), 2):
        missing_set = set(missing)
        eleven_sets.append(tuple(x for x in range(p) if x not in missing_set))
    require(len(eleven_sets) == 78, "eleven-sheet set count")

    representation_census = Counter()
    for left in eleven_sets:
        for right in eleven_sets:
            counts = cycle_convolution((left, right), p)
            for c, reps in enumerate(counts):
                require(reps == convolution_count(left, right, c, p),
                        "dense convolution mismatch")
                require(reps >= 9, "eleven-sheet floor failure")
                representation_census[reps] += 1

    expected_census = Counter({9: 55770, 10: 22308, 11: 1014})
    require(representation_census == expected_census,
            "eleven-sheet representation census")

    # Thin transport retains one clutch class; thick transport saturates all.
    thin_counts = cycle_convolution(((2,), (5,), (9,)), p)
    require(sum(value > 0 for value in thin_counts) == 1, "thin clutch support")
    require(max(thin_counts) == 1, "thin representation multiplicity")

    print("THM2642 cyclic difference-relation exact referee")
    print(f"prime={p}")
    print(f"p7_nonempty_two_set_pairs={pair_count} cauchy_davenport=PASS convolution=PASS")
    print(f"sharp_two_edge_size_profiles={extremizer_profiles}")
    print("near_threshold_sizes=(6,7) excess=11 support=12 missing=(12,)")
    print("threshold_sizes=(6,8) excess=12 support=13")
    print(f"eleven_sheet_pairs={len(eleven_sets) ** 2} carry_triples={len(eleven_sets) ** 2 * p}")
    print("representation_census={9:55770,10:22308,11:1014}")
    print("section_census={117:55770,130:22308,143:1014}")
    print("eleven_sheet_min_per_base=9 min_sections=117")
    print("thin_cycle_supported_clutches=1 sections_on_supported_clutch=13")
    print("PASS")


if __name__ == "__main__":
    main()
