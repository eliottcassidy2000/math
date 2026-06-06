#!/usr/bin/env python3
"""
S630 part 3 — the REFRAME: min M  <->  max resonance (3-term fold count).
Connections (Explore): S550o resonance energy, S577 'the fold IS the delta-clock witness',
S579 2-adic tower, S617/S621 additive chains. Hypothesis: the extremal (min-M) multiple-of-19
config is the FOLD-RICHEST one (most a+b=c relations) -> search by maximizing folds, not min M.
"""
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
import itertools, random, sys
sys.path.insert(0, '04-computation')
from lrc_fast_s628 import gap_shells, gap_brute

N = 19

def folds(S):
    s = set(S); return sum(1 for a in S for b in S if a <= b and a+b in s)

def resonance_min_len(S):
    """minimal resonance length proxy: shortest |m1|+|m2|+|m3| with m.v=0; ~ presence of folds (len 3)."""
    s = set(S)
    if any(a+b in s for a in S for b in S if a < b): return 3
    return 99   # no 3-term relation

if __name__ == "__main__":
    print("[1] correlation: fold-count vs gap_shells over multiple-of-19 configs")
    rng = random.Random(2); data = []
    for _ in range(8000):
        S = {19}
        while len(S) < 18: S.add(rng.randint(1, 38))
        S = tuple(sorted(S))
        if len(S) != 18 or reduce(gcd, S) != 1: continue
        data.append((folds(S), float(gap_shells(list(S), N))))
    import statistics
    hi = [g for f, g in data if f >= 12]; lo = [g for f, g in data if f <= 3]
    print(f"   configs: {len(data)};  mean gap_shells when folds>=12: {statistics.mean(hi):.4f} (n={len(hi)})")
    print(f"                          when folds<=3 : {statistics.mean(lo):.4f} (n={len(lo)})")
    print("   => fold-rich configs are TIGHTER (smaller gap) — the reframe holds" if statistics.mean(hi) < statistics.mean(lo) else "   => no clear correlation")

    print("\n[2] FOLD-MAXIMIZING construction: greedy additive chain containing 19, 18 speeds")
    # build a dense additive chain: start {1,2}, repeatedly add the smallest new sum; ensure 19 in it
    best = None
    for seed in range(1, 6):
        S = [1, seed+1] if seed > 0 else [1, 2]
        S = sorted(set([1, 2]))
        while len(S) < 18:
            sums = sorted({a+b for a in S for b in S if a <= b} - set(S))
            S = sorted(set(S) | {sums[0]})
        if 19 not in S:
            # force a multiple of 19: replace the element nearest 19 with 19
            S = sorted(set(S) | {19}); S = S[:18] if len(S) > 18 else S
            while len(S) < 18:
                sums = sorted({a+b for a in S for b in S if a <= b} - set(S)); S = sorted(set(S)|{sums[0]})
        S = sorted(set(S))[:18]
        if len(S) == 18 and 19 in S and reduce(gcd, S) == 1:
            g = gap_brute(S)
            print(f"   chain {tuple(S)}: folds={folds(S)} M={g}={float(g):.5f}  ==2/37? {g==Fr(2,37)}  <1/18? {g<Fr(1,18)}")
            if best is None or g < best[0]: best = (g, S)

    print("\n[3] fold-max SEARCH among multiple-of-19 configs (hill-climb on fold-count), report M of best")
    cur = sorted(set([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,19]))  # AP-ish + 19
    curf = folds(cur)
    for it in range(8000):
        i = rng.randrange(18); newv = rng.randint(1, 40)
        cand = sorted(set(cur[:i]+cur[i+1:]+[newv]))
        if len(cand) != 18 or 19 not in cand and not any(v%19==0 for v in cand): continue
        if reduce(gcd, cand) != 1: continue
        f = folds(cand)
        if f > curf: cur, curf = cand, f
    g = gap_brute(cur)
    print(f"   fold-max config {tuple(cur)}: folds={curf}  M={g}={float(g):.5f}")
    print(f"   2/37={float(Fr(2,37)):.5f}  1/18={float(Fr(1,18)):.5f}  -> reframe finds min M ~ {float(g):.5f}")
