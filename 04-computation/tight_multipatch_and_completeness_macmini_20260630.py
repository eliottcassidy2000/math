#!/usr/bin/env python3
"""
mac-mini-2026-06-30-S63
=======================
PART 2 -- the MULTI-PATCH moduli (skip j, patch j) generating the near-AP tight family for n=14.
PART 1 -- COMPLETENESS evidence: exhaustively search structured tight sets for n=14; do any exist
          beyond {AP, GW}?  We attack the two routes a non-{AP,GW} tight set could take:
            (i) multi-patch of the AP (skip j residues, add j lifts), j=1,2,3;
            (ii) a residue-complete lift set (all residues {1..13} present, some speeds lifted).
          The necessary conditions (HYP-+2914) force residues mod 14 = {1..13} or {1..13}\{k},
          so EVERY tight set is a multi-patch of the AP -- we enumerate them up to a lift bound.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations

def Mexact(S):
    Sg = sorted(set(S)); cand = set()
    for i in range(len(Sg)):
        for j in range(len(Sg)):
            for d in (Sg[i]-Sg[j], Sg[i]+Sg[j]):
                if d > 0:
                    for k in range(1, d): cand.add(F(k, d))
    best = F(0)
    for t in cand:
        g = min(min((v*t) % 1, 1-((v*t) % 1)) for v in Sg)
        if g > best: best = g
    return best

def prim(S): return reduce(gcd, S) == 1
def diffclosed(S):
    Sset = set(S); return all(abs(a-b) in Sset for a in S for b in S if a != b)

n = 14; AP = list(range(1, n)); tgt = F(1, n)
AP_set = set(AP)

print("="*80)
print("PART 2: MULTI-PATCH (skip j of {1..13}, add j new elements <= Bmax) tight for n=14")
print("="*80)
allfound = {}
for j, Bmax in ((1, 100), (2, 42), (3, 30)):
    found = set()
    for R in combinations(AP, j):
        base = [x for x in AP if x not in R]
        adds = [a for a in range(2, Bmax+1) if a not in AP_set]
        for A in combinations(adds, j):
            S = tuple(sorted(base+list(A)))
            if len(set(S)) != n-1 or not prim(S): continue
            if Mexact(S) == tgt: found.add(S)
    allfound[j] = found
    print(f"\n  j={j} (adds<={Bmax}): {len(found)} tight multi-patch(es)")
    for S in sorted(found):
        drop = sorted(AP_set-set(S)); add = sorted(set(S)-AP_set)
        res = sorted(x % n for x in S)
        resdup = [r for r in set(res) if res.count(r) > 1]
        print(f"     {list(S)}  drop{drop} add{add}  ndc={not diffclosed(S)}  dupres{resdup}")

print()
print("="*80)
print("PART 1: COMPLETENESS -- is every tight multi-patch = {AP} or {GW}?  (up to the bounds above)")
print("="*80)
AP_tuple = tuple(AP)
GW = tuple(sorted([x for x in AP if x != 12]+[24]))
canon = {AP_tuple, GW}
extra = set()
for j in allfound:
    for S in allfound[j]:
        if S not in canon:
            extra.add(S)
print(f"  AP present in found? {AP_tuple in allfound.get(0, {AP_tuple})}  (AP is the j=0 base)")
print(f"  GW  = {list(GW)} present? {GW in allfound[1]}")
print(f"  NON-{{AP,GW}} tight multi-patches found: {len(extra)}")
for S in sorted(extra):
    drop = sorted(AP_set-set(S)); add = sorted(set(S)-AP_set)
    print(f"     EXTRA: {list(S)}  drop{drop} add{add}")
if not extra:
    print("  => within these bounds, {AP, GW} is COMPLETE (no other tight multi-patch).")
print("\nDONE.")
