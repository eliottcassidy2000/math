#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Lean targeted + random adversarial hunt for tower-broken sub-THR2 cores."""
import itertools, random, sys
from fractions import Fraction
exec(open('04-computation/lrc14_adversarial_verify_eliott.py').read().split('# ---- DECISIVE')[0].split('if __name__')[0])

cex=[]; best=None; checked=0
# Deterministic: drop one tower bit, keep rest of [1,13] (=12 elts tower-broken),
# then swap exactly ONE kept element for a tail up to T=300 (the densest tower-broken family).
T=300
for b in [1,2,4,8]:
    base=[d for d in range(1,14) if d!=b]   # size 12
    # nsw=0 (base itself)
    for nsw in (0,1):
        if nsw==0:
            cands=[tuple(base)]
        else:
            cands=[]
            for out in base:
                keep=[d for d in base if d!=out]
                for t in range(14,T+1):
                    cands.append(tuple(sorted(keep+[t])))
        for C in cands:
            if len(C)!=12: continue
            if not primitive(C): continue
            if has_tower(C): continue
            checked+=1
            L=lonely_measure(C)
            if best is None or L<best[0]: best=(L,C)
            if L<THR2: cex.append((L,C))
print(f"[det: drop tower bit, swap<=1 kept->tail<=300] checked={checked}")
print(f"  min tower-broken meas={best[0]}={float(best[0]):.6f} at {best[1]}")
print(f"  sub-THR2 counterexamples: {len(cex)}")
for L,C in sorted(cex)[:10]: print("   ***",L,C)
sys.stdout.flush()

# Random tower-broken, wide
random.seed(7); NR=300000; rbest=None; rcex=0
for _ in range(NR):
    miss=random.choice([1,2,4,8])
    pool=[x for x in range(1,50) if x!=miss]
    C=tuple(sorted(random.sample(pool,12)))
    if has_tower(C) or not primitive(C): continue
    L=lonely_measure(C)
    if rbest is None or L<rbest[0]: rbest=(L,C)
    if L<THR2: rcex+=1; print("   RAND CEX",L,C); sys.stdout.flush()
print(f"[random {NR} tower-broken in 1..49] min={rbest[0]}={float(rbest[0]):.6f} at {rbest[1]}")
print(f"  sub-THR2 random counterexamples: {rcex}")
print("VERDICT:", "COUNTEREXAMPLE" if (cex or rcex) else "NO tower-broken core beat THR2")
