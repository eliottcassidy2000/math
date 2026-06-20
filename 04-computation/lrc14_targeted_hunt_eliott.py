#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""TARGETED adversarial hunt: can a TOWER-BROKEN core beat THR2?
Strategy: tower-broken means missing >=1 of {1,2,4,8}. The smallest-measure
cores keep small dense elements. So fix: drop exactly one tower bit, keep the
rest of [1,13] as much as possible, and let the 12th element roam to large tails.
Also do random sampling over wide range for tower-broken cores."""
import itertools, random
from fractions import Fraction
exec(open('04-computation/lrc14_adversarial_verify_eliott.py').read().split('# ---- DECISIVE')[0].split('if __name__')[0])

best = None
checked = 0
cex = []

# Direction 1: tower-broken, deterministic. Remove one tower bit b, then choose
# 12 elements from ([1,13]\{b}) plus tails up to T. Keep as many small as possible:
# i.e. take all of [1,13]\{b} (that's 12 elements) then swap k of the largest for tails.
T = 200
for b in [1,2,4,8]:
    base = [d for d in range(1,14) if d != b]  # 12 elements, tower-broken
    # swap up to 3 of these out for tails (to explore)
    for nsw in range(0,4):
        for out in itertools.combinations(base, nsw):
            keep = [d for d in base if d not in out]
            for tails in itertools.combinations(range(14,T+1), nsw):
                C = tuple(sorted(keep)+list(tails))
                if len(C)!=12: continue
                if not primitive(C): continue
                if has_tower(C): continue  # must stay tower-broken
                checked += 1
                L = lonely_measure(C)
                if best is None or L < best[0]:
                    best = (L, C)
                if L < THR2:
                    cex.append((L,C))
        if nsw>=2:  # limit big tail loops to nsw<=1 fully; nsw 2,3 too big with T=200
            break
print(f"[deterministic, tails<=200, swaps<=1 fully + base] checked {checked}")
if best: print(f"  min tower-broken meas = {best[0]} = {float(best[0]):.6f} at {best[1]}")
print(f"  sub-THR2 tower-broken counterexamples: {len(cex)}")
for L,C in sorted(cex)[:10]: print("   ***", L, C)

# Direction 2: random tower-broken cores, wide range
random.seed(1)
NR = 200000
rbest=None; rcex=0
for _ in range(NR):
    # build a tower-broken 12-set: drop one tower bit guaranteed missing
    miss = random.choice([1,2,4,8])
    pool = [x for x in range(1,40) if x!=miss]
    C = tuple(sorted(random.sample(pool, 12)))
    if has_tower(C): continue
    if not primitive(C): continue
    L = lonely_measure(C)
    if rbest is None or L<rbest[0]: rbest=(L,C)
    if L<THR2: rcex+=1; print("   RAND CEX", L, C)
print(f"\n[random {NR} tower-broken in 1..39] min meas={rbest[0]}={float(rbest[0]):.6f} at {rbest[1]}")
print(f"  sub-THR2 random counterexamples: {rcex}")
