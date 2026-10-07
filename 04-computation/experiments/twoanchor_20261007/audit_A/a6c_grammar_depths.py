#!/usr/bin/env python3
"""Audit A: per-depth masses of the merge shapes (from a6b_grammar), to compare with the session's heads_plus1 union
(whose generator gen() only builds words of total <= SMAX - 1, i.e. sum(u) <= 21 for SMAX = 22)."""
import sys
import a6b_grammar as G
from collections import defaultdict
G.LMAX = 23
for r, pre in ((1, (1, 1)), (2, (1, 0, 0))):
    pats = G.dfs(pre)
    cum = defaultdict(float); nheads = defaultdict(int)
    for p in pats:
        sh, i, ok, uw, vw = G.classify(p)
        m = 2.0**(-(len(p) - len(pre)))
        for d in range(len(p), 24):
            cum[(d, sh, i)] += m
        if sh == 'child-ladder' and i == 1:
            nheads[len(p)] += 1
    for d in (12, 20, 21, 22, 23):
        i1 = cum[(d, 'child-ladder', 1)]
        tot = sum(v for (dd, s, i), v in cum.items() if dd == d)
        print(f"r={r} by post-run depth {d}: child-ladder i=1 mass {i1:.6f}; all shapes {tot:.6f}; i=1 share {i1/tot:.1%}")
    # number of distinct i=1 source heads u (absorption depth = sum(u)+1) with sum(u) <= 21
    print(f"r={r}: i=1 absorbing patterns with depth <= 22 (= sum(u) <= 21): {sum(c for d, c in nheads.items() if d <= 22)}")
