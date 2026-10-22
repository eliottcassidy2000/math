"""Exact number of labelled tournaments with a given score sequence (memoized last-vertex recurrence), used to certify
the completeness of the almost-regular class lists by labelled mass: almost-regular n-tournaments (scores n/2, n/2 - 1)
correspond bijectively to regular (n+1)-tournaments with a distinguished vertex (unique completion).
Session opus-2026-10-22-S54.  Usage: python count_regular.py"""
from functools import lru_cache
from itertools import combinations
from math import comb


@lru_cache(maxsize=None)
def F(d):
    """d: sorted tuple of required out-degrees on m labelled vertices."""
    m = len(d)
    if m == 0:
        return 1
    if sum(d) != m * (m - 1) // 2 or min(d) < 0 or max(d) > m - 1:
        return 0
    last, rest = d[-1], d[:-1]
    total = 0
    # the last vertex beats a set O of size `last`; every other vertex beats it and loses one required out-degree
    # count by multiset: group equal degrees
    from collections import Counter
    cnt = sorted(Counter(rest).items())
    def rec(i, need, acc_mult, new):
        nonlocal total
        if i == len(cnt):
            if need == 0:
                total += acc_mult * F(tuple(sorted(new)))
            return
        val, c = cnt[i]
        for k in range(0, min(c, need) + 1):          # k of these are beaten by `last` (keep degree), c - k beat it
            rec(i + 1, need - k, acc_mult * comb(c, k), new + [val] * k + [val - 1] * (c - k))
    rec(0, last, 1, [])
    return total


for m in (3, 5, 7, 9, 11):
    print(f'labelled regular {m}-tournaments: {F(tuple([(m - 1) // 2] * m))}')
