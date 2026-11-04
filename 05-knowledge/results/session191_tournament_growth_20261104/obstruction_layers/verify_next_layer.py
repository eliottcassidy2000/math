#!/usr/bin/env python3
"""Exact finite certificate for ADDENDUM.md; standard library only."""
from collections import Counter
from fractions import Fraction
from itertools import permutations
from math import comb
import json
import verify as v


def choose(n, k):
    return 0 if n < k or n < 0 else comb(n, k)


def expected(r, theta):
    total = 6*max(r-1, 0)*(theta >= Fraction(5, 9))
    total += choose(r-1, 2)*(Fraction(5, 9) <= theta < Fraction(25, 27))
    total += choose(r-2, 5)*(theta >= Fraction(625, 729))
    total += (2*choose(r-2, 2)+6*choose(r-2, 3)+8*choose(r-2, 4))*(theta >= Fraction(25, 27))
    return total


def main():
    profiles = {}
    class_counts = []
    strong_counts = []
    by_rows = {}
    for n in range(1, 8):
        reps = v.frame_census(n)
        class_counts.append(len(reps))
        profiles[n] = []
        for rows in reps:
            if n < 3 or not v.strong(rows):
                continue
            counts = v.subset_counts(rows)
            full = (1 << n)-1
            p = {"rows": rows, "h": counts[full],
                 "d": max(counts[full ^ (1 << i)] for i in range(n))}
            profiles[n].append(p)
            by_rows[rows] = p
        strong_counts.append(len(profiles[n]))
    v.require(class_counts == [1, 1, 2, 4, 12, 56, 456], "incomplete census")
    q = profiles[4][0]
    v.require(len(profiles[4]) == 1 and (q['h'], q['d']) == (5, 3), "Q data")
    v.require(sorted(p['h'] for p in profiles[5]) == [9, 9, 11, 13, 15, 15], "S5 counts")
    v.require(all(p['d'] == 5 for p in profiles[5]), "S5 deletion maxima")
    m6 = [p for p in profiles[6] if p['h'] <= 16]
    m7 = [p for p in profiles[7] if p['h'] <= 26]
    v.require(len(m6) == 2 and all((p['h'], p['d']) == (15, 9) for p in m6), "M6 classification")
    v.require(len(m7) == 1 and (m7[0]['h'], m7[0]['d']) == (25, 15), "M7 classification")
    core = m7[0]['rows']
    v.require(v.literal_count(core) == 25, "literal M7 H")
    v.require(max(v.literal_count(core, tuple(w for w in range(7) if w != i))
                  for i in range(7)) == 15, "literal M7 deletion maximum")

    ear_minimum = None
    ear_survivors = []
    strong_ears = 0
    for mask in range(128):
        rows = list(core)+[0]
        for i in range(7):
            if mask >> i & 1:
                rows[i] |= 1 << 7
            else:
                rows[7] |= 1 << i
        v.check_tournament(rows)
        if not v.strong(rows):
            continue
        strong_ears += 1
        counts = v.subset_counts(rows)
        d = max(counts[255 ^ (1 << i)] for i in range(8))
        ear_minimum = d if ear_minimum is None else min(d, ear_minimum)
        if d <= 26:
            ear_survivors.append(tuple(rows))
    v.require(not ear_survivors, "unexpected order-eight template")
    v.require(125**6 > 9*81**6, "E<=6 cutoff")

    words = set()
    def add_multiset(parts):
        rows = tuple(p['rows'] for p in parts)
        words.update(permutations(rows))
    for p in profiles[5]:
        add_multiset([p])
    add_multiset([q, q])
    for p in profiles[7]:
        add_multiset([p, q])
    for p in profiles[6]:
        for w in profiles[5]:
            add_multiset([p, w])
        add_multiset([p, q, q])
    for p in profiles[5]:
        for w in profiles[5]:
            add_multiset([p, w, q])
        add_multiset([p, q, q, q])
    add_multiset([q]*5)

    templates = []
    pattern_counts = Counter()
    for word in sorted(words):
        parts = [by_rows[rows] for rows in word]
        s = sum(map(len, word))
        v.require((s+1) % 3 == 0, "template congruence")
        t = (s+1)//3
        h = 1
        for p in parts:
            h *= p['h']
        exceptional_d = max(h//p['h']*p['d'] for p in parts)
        upper = min(Fraction(h, 3**t), Fraction(1))
        lower = max(Fraction(exceptional_d, 3**t), Fraction(1, 3))
        # In every survivor, exceptional deletion already dominates
        # triangle deletion. Check rather than assume the a=0 boundary.
        if lower < upper:
            v.require(exceptional_d*3 >= h, "a=0 boundary matters")
            pattern = tuple(sorted(map(len, word)))
            pattern_counts[pattern] += 1
            templates.append((t, len(word), lower, upper, pattern))
    expected_patterns = {(5,): 6, (4, 4): 1, (4, 7): 2,
                         (4, 4, 6): 6, (4, 4, 4, 5): 8,
                         (4, 4, 4, 4, 4): 1}
    v.require(dict(pattern_counts) == expected_patterns, "template table mismatch")

    boundaries = [Fraction(1, 3), Fraction(5, 9), Fraction(625, 729),
                  Fraction(25, 27), Fraction(1)]
    probes = boundaries[:-1]+[(x+y)/2 for x, y in zip(boundaries, boundaries[1:])]
    comparisons = 0
    for r in range(1, 41):
        for theta in probes:
            observed = sum(choose(r-t+k, k) for t, k, lo, hi, _ in templates
                           if r >= t and lo <= theta < hi)
            v.require(observed == expected(r, theta), f"formula mismatch r={r} theta={theta}")
            comparisons += 1
    v.require(expected(3, Fraction(21, 27)) == 13, "B21 control")
    v.require(expected(3, Fraction(24, 27)) == 13, "before Q-square disappears")
    v.require(expected(3, Fraction(25, 27)) == 12, "Q-square boundary hostile")
    v.require(expected(4, Fraction(74, 81)) == 21, "B74 control")
    v.require(expected(4, Fraction(75, 81)) == 20, "B75 replacement control")

    print(json.dumps({
        "status": "ALL CHECKS PASSED",
        "ordinary_classes_1_to_7": class_counts,
        "nontrivial_strong_classes_1_to_7": strong_counts,
        "M6_witnesses": m6,
        "M7_witness": m7[0],
        "order8_extensions_checked": 128,
        "strong_order8_extensions": strong_ears,
        "minimum_maximum_deletion_of_strong_extensions": ear_minimum,
        "strong_order8_extensions_with_D_at_most26": len(ear_survivors),
        "exceptional_ordered_class_lists_checked": len(words),
        "surviving_ordered_templates": [
            {"orders": list(pattern), "count": count}
            for pattern, count in sorted(pattern_counts.items())],
        "threshold_intervals": [
            {"orders": list(pattern), "minimum_r": t,
             "lower_closed": str(lo), "upper_open": str(hi)}
            for t, _, lo, hi, pattern in sorted(set(templates))],
        "endpoint_and_chamber_formula_comparisons": comparisons,
        "controls": {"A1(21)": 13, "A1(24)": 13, "A1(25)": 12,
                     "A1(74)": 21, "A1(75)": 20}
    }, indent=2))


if __name__ == '__main__':
    main()
