#!/usr/bin/env python3
"""alphabets.py -- value sets of ribbon words over restricted alphabets, by direct enumeration."""
import sys, itertools
from ribbon import transfer_system, letters, count_with_matrices, word_str

w = 2
sysw = transfer_system(w)
order, states, mats, final = sysw
F = frozenset
alphabets = {
    "{1},{2}": [F([1]), F([2])],
    "{1},{2},{1,2}": [F([1]), F([2]), F([1, 2])],
    "{},{2}": [F(), F([2])],
    "{2},{1,2}": [F([2]), F([1, 2])],
    "{1},{1,2}": [F([1]), F([1, 2])],
}
for name, alph in alphabets.items():
    vals = {}
    for n in range(1, 13):
        for word in itertools.product(alph, repeat=n):
            h = count_with_matrices(w, list(word), sysw)
            if h not in vals or len(vals[h]) > n:
                vals[h] = word
    mx = max(vals)
    # misses among odd numbers up to the largest value attained with words of length <= 11 (so that
    # length-12 words cannot change the picture below it): use max over n<=11 as the horizon
    horizon = 0
    for h, word in vals.items():
        if len(word) <= 11 and h > horizon:
            horizon = h
    miss = [h for h in range(1, horizon + 1, 2) if h not in vals]
    print(f"alphabet {name}: {len(vals)} distinct values with n<=12, max {mx}; horizon(n<=11)={horizon}; odd misses below horizon: {len(miss)}")
    print("   first misses:", miss[:40])
    print("   smallest values:", sorted(vals)[:25])
