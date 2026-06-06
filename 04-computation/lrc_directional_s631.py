#!/usr/bin/env python3
"""
S631 Part 1 — LRC with reversible runners (opposite directions = negative speeds).
Single flip is INVISIBLE (||v t||=||-v t||); the content lives in PAIRWISE relative speeds:
runner i (dir e_i) vs j (dir e_j) has relative speed e_i v_i - e_j v_j -> SAME dir = difference,
OPPOSITE dir = sum. A directional assignment e in {+-1} is a TOURNAMENT orientation; the
"mutual loneliness" gap depends on e. Question: does opposite-direction (sums) push the natural
modulus n -> 2n-1 (the pair-sum sieve, THM-401)?
"""
from fractions import Fraction as Fr
from math import gcd
import itertools, sys
sys.path.insert(0, '04-computation')
from lrc_fast_s628 import gap_brute   # max_t min over a set of |speeds|

def norm(x):
    f = x - (x.numerator // x.denominator)
    return f if f <= Fr(1, 2) else 1-f

def mutual_gap(speeds, eps):
    """max_t min_{i<j} ||(e_i v_i - e_j v_j) t||  — the reversible mutual-loneliness gap."""
    rel = []
    for i in range(len(speeds)):
        for j in range(i+1, len(speeds)):
            r = abs(eps[i]*speeds[i] - eps[j]*speeds[j])
            if r != 0: rel.append(r)
    if not rel: return Fr(1, 2)
    return gap_brute(rel), sorted(set(rel))

if __name__ == "__main__":
    print("Single sign-flip invisibility (observer gap):")
    for S in [[1,2,3,4],[1,3,4,7]]:
        g0 = gap_brute(S)
        g1 = gap_brute([-S[0]]+S[1:])
        print(f"   {S}: M={g0}, M(flip one)={g1}, equal={g0==g1}")

    print("\nMUTUAL gap over all direction assignments e (tournament orientations), AP {1..n-1}:")
    for n in (4, 5, 6):
        S = list(range(1, n)); m = 2*n-1
        results = {}
        for eps in itertools.product([1, -1], repeat=len(S)):
            if eps[0] == -1: continue          # global sign symmetry
            g, rel = mutual_gap(S, eps)
            results[eps] = (g, tuple(rel))
        gaps = sorted(set(g for g, _ in results.values()))
        # the all-same (all differences) vs the ones producing the largest relative speed (sums)
        allsame = mutual_gap(S, [1]*len(S))
        maxrel = max(max(r) for g, r in results.values())
        print(f"  n={n} (2n-1={m}): #distinct mutual-gaps over orientations = {len(gaps)};"
              f" range [{min(gaps)},{max(gaps)}]")
        print(f"     all-same-direction (pure differences): gap={allsame[0]} rel-speeds={allsame[1]}")
        # the orientation maximizing the gap (easiest mutual loneliness) and its rel-speed set
        best = max(results.items(), key=lambda kv: kv[1][0])
        worst = min(results.items(), key=lambda kv: kv[1][0])
        print(f"     BEST orientation e={best[0]} gap={best[1][0]} rel-speeds={best[1][1]}")
        print(f"     WORST orientation e={worst[0]} gap={worst[1][0]} rel-speeds={worst[1][1]}")
        # do the witness denominators hit 2n-1 in the sum (opposite) orientations?
        dens = set()
        for eps, (g, rel) in results.items():
            dens.add(g.denominator)
        print(f"     witness-denominators across orientations: {sorted(dens)}  (contains 2n-1={m}? {m in dens})")
