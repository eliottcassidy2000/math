#!/usr/bin/env python3
"""Recompute the arithmetic of the note's section 7 (coincidence ledger)."""
from math import comb
L = {15, 16, 17, 23, 24, 25, 31, 32, 33, 34, 38, 39, 40, 41, 42}   # the note's landmark set
W = range(20, 51)
T = [20, 24, 32, 34, 41, 60, 94, 105]                              # Thomassen orders (MathWorld)
inW = [x for x in L if x in W]
Tw = [t for t in T if t in W]
hits = [t for t in Tw if t in L]
print("|L| =", len(L), " |L in [20,50]| =", len(inW), sorted(inW), " |[20,50]| =", len(W))
print("Thomassen orders in [20,50]:", Tw, " hits:", hits)
for marked, label in ((len(inW), "correct (12/31)"), (len(L), "as written in the note (15/31)")):
    exp = len(Tw) * marked / len(W)
    p4 = sum(comb(marked, k) * comb(len(W) - marked, len(Tw) - k) for k in range(len(hits), len(Tw) + 1)) / comb(len(W), len(Tw))
    print("%-32s expected hits %.2f   P(>= %d hits) hypergeometric = %.4f" % (label, exp, len(hits), p4))
