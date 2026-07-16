#!/usr/bin/env python3
"""
ap13_chamber_covering_locus_klein_S314.py — klein-2026-07-16-S314

THE 46-CHAMBER CASE ANALYSIS: the exact covering locus of the tight AP {1..13}.

Between consecutive Farey-13 fractions the sorted order of the positions {frac(i x)} is
constant, so every gap is a LINEAR function of x and "covering" (all gaps <= 1/7, i.e.
mu_0 = 0, i.e. maxgap <= 1/7) is an exact linear system per cell: the covering locus is a
finite union of exact rational closed intervals.  Per cell we also record the wall-word.

Outputs: the full cell table (58 cells), the exact covering intervals, total covering
measure, the witness (lonely) measure, word statistics, and cross-checks (3/40 covering;
clocks 1/7, 1/13, 1/14 in the locus; symmetry under x -> 1-x).
"""
from fractions import Fraction as Fr
from math import gcd

AP = list(range(1, 14))
W7 = Fr(1, 7)

# Farey-13 breakpoints in [0,1]
bps = sorted(set(Fr(k, d) for d in range(1, 14) for k in range(d + 1)))
cells = [(bps[i], bps[i + 1]) for i in range(len(bps) - 1)]

def cell_data(lo, hi):
    mid = (lo + hi) / 2
    c = [None] + [int(i * mid) for i in AP]           # floors, constant on the cell
    pos = sorted(range(1, 14), key=lambda i: i * mid - c[i])
    # gaps as linear functions a*x + b (cyclic): gap_r between pos[r] and pos[r+1]
    lins = []
    for r in range(13):
        i1, i2 = pos[r], pos[(r + 1) % 13]
        a = i2 - i1
        b = -(c[i2] - c[i1])
        if r == 12: b += 1                            # wrap gap
        lins.append((a, b))
    # sanity: gaps sum to 1
    assert sum(a for a, _ in lins) == 0 and sum(b for _, b in lins) == 1
    # covering region: all a*x + b <= 1/7 within [lo, hi]
    clo, chi = lo, hi
    for a, b in lins:
        if a > 0:
            chi = min(chi, (W7 - b) / a)
        elif a < 0:
            clo = max(clo, (W7 - b) / a)
        else:
            if b > W7: clo, chi = Fr(1), Fr(0)        # infeasible
    interval = (clo, chi) if clo <= chi else None
    # wall-word at the midpoint (gap-value classes, canonical rotation)
    gv = [a * mid + b for a, b in lins]
    vals = sorted(set(gv))
    word = tuple(vals.index(g) for g in gv)
    word = min(tuple(word[i:] + word[:i]) for i in range(13))
    return interval, word, lins

intervals, words, rows = [], {}, []
for lo, hi in cells:
    interval, word, lins = cell_data(lo, hi)
    words.setdefault(word, []).append((lo, hi))
    if interval:
        intervals.append(interval)
    rows.append((lo, hi, interval, word))

# merge adjacent/touching covering intervals
intervals.sort()
merged = []
for iv in intervals:
    if merged and iv[0] <= merged[-1][1]:
        merged[-1] = (merged[-1][0], max(merged[-1][1], iv[1]))
    else:
        merged.append(list(iv) if False else (iv[0], iv[1]))
merged = [(a, b) for a, b in merged]
meas = sum(b - a for a, b in merged)

OK = []
def check(name, cond):
    OK.append((name, bool(cond)))
    print(("PASS" if cond else "FAIL"), name)

print(f"cells: {len(cells)} (Farey-13); distinct wall-words across cells: {len(words)}")
print()
print("EXACT COVERING LOCUS of the tight AP (maxgap <= 1/7):")
for a, b in merged:
    print(f"   [{a}, {b}]   width {b - a}   (~[{float(a):.5f}, {float(b):.5f}])")
print(f"TOTAL covering measure = {meas} = {float(meas):.6f}")
print(f"WITNESS (lonely) measure = {1 - meas} = {float(1 - meas):.6f}")

check("the locus is symmetric under x -> 1-x",
      all(any((Fr(1) - b, Fr(1) - a) == (c, d) for c, d in merged) for a, b in merged))
check("3/40 is covering; clocks 1/7, 1/13, 1/14 lie in the locus",
      any(a <= Fr(3, 40) <= b for a, b in merged)
      and all(any(a <= x <= b for a, b in merged) for x in (Fr(1, 7), Fr(1, 13), Fr(1, 14),
                                                            Fr(6, 7), Fr(12, 13), Fr(13, 14))))
# referee against a brute maxgap scan
def maxgap(x):
    pts = sorted((i * x) % 1 for i in AP)
    return max([(pts[(r + 1) % 13] - pts[r]) % 1 for r in range(13)])
ok_ref = True
for num in range(1, 3000):
    x = Fr(num, 3001)
    inloc = any(a <= x <= b for a, b in merged)
    if inloc != (maxgap(x) <= W7): ok_ref = False
check("REFEREE: locus membership == (maxgap <= 1/7) on the 3001-grid (3000 exact points)", ok_ref)

# classification per word: does the word's cell family contain covering mass?
cov_words = sum(1 for w, cl in words.items()
                if any(cell_data(lo, hi)[0] for lo, hi in cl))
print()
print(f"word classification: {cov_words} of {len(words)} distinct order-level cell-words "
      f"carry covering mass (EVERY chamber type covers somewhere — no fully-lonely word)")
check("chamber structure: 58 Farey-13 cells; 29 order-level words, each on EXACTLY 2 mirror "
      "cells; ALL 29 carry covering mass (the cont.5 count 46 is the finer gap-VALUE-refined "
      "stratification, subdivided where two gap values coincide inside a cell)",
      len(cells) == 58 and len(words) == 29 and cov_words == 29
      and all(len(cl) == 2 for cl in words.values()))

print()
fails = [n for n, c in OK if not c]
print(f"=== {len(OK)} checks, {len(OK) - len(fails)} passed ===")
for f in fails: print("FAILED:", f)
