# opus-2026-07-19-S393 -- HYP-7790: THE n=14 STABILITY GAP (branch-2 sharpened).
#
# ARCHAEOLOGY.  klein-S313d/e proved the n=13 analogue (12 speeds, floor 1/13):
#   THM-1004 Hamming-1 rigidity -- replacing ANY single element of {1,...,12}
#     forces M >= 2/25, with equality EXACTLY at {1,...,11,24};
#   THM-1005 Hamming-2 rigidity -- same bound at radius 2.  Radius >= 3 OPEN.
# So at n=13 there is a STABILITY GAP: M = 1/13 (the AP, tight) or M >= 2/25,
# and nothing in between.  Note {1,...,11,24} is the SAME 12->24 substitution
# that produced my n=14 second tight family {1,...,11,13,24} (THM-1120).
#
# THE n=14 ANALOGUE, which is exactly my THM-1215 branch-2 target sharpened:
#   is there a gap  M = 1/14  or  M >= 1/13,  with NOTHING strictly between?
# The near-extremal ladder of THM-1215 (1/14, 1/13, 1/12) is consistent with it.
# A family with M strictly in (1/14, 1/13) would break the gap; a family with
# M < 1/14 would break LRC(14) itself.  Adversarial search for both.
from fractions import Fraction as F
from functools import reduce
from math import gcd
from itertools import combinations
import random
LO=F(1,14); HI=F(1,13)
def fn(x):
    r=x-int(x)
    if r<0: r+=1
    return min(r,1-r)
def gap_at(V,t): return min(fn(t*v) for v in V)
def M(V):
    D=set()
    for a,b in combinations(V,2):
        D.add(a+b)
        if a!=b: D.add(abs(a-b))
    for v in V: D.add(2*v)
    best=F(0)
    for q in sorted(D):
        if q<2: continue
        for p in range(1,q):
            g=gap_at(V,F(p,q))
            if g>best: best=g
    return best

print(f"    LO = 1/14 = {float(LO):.8f}   HI = 1/13 = {float(HI):.8f}")
print()
print("(1) DOES ANY FAMILY LAND STRICTLY INSIDE (1/14, 1/13)?")
print("    -- adversarial search minimising distance into the open interval")
random.seed(393)
inside=[]; seen=set(); n=0
best_below=(F(1),None)
for trial in range(26):
    V=sorted(random.sample(range(1,42),13))
    if reduce(gcd,V)!=1: continue
    def score(W):
        m=M(W)
        if LO<m<HI: return F(0)          # inside the gap: jackpot
        return min(abs(m-LO),abs(m-HI))  # distance to the interval
    cur=score(V); stall=0
    for step in range(150):
        W=list(V); i=random.randrange(13)
        W[i]=max(1,W[i]+random.choice([-3,-2,-1,1,2,3]))
        W=sorted(set(W))
        if len(W)!=13 or reduce(gcd,W)!=1: continue
        s2=score(W)
        if s2<=cur:
            if s2<cur: stall=0
            V,cur=W,s2
        else: stall+=1
        if cur==0: break
        if stall>100: break
    n+=1
    m=M(V)
    if LO<m<HI and tuple(V) not in seen:
        seen.add(tuple(V)); inside.append((m,list(V)))
    if m<best_below[0]: best_below=(m,list(V))
print(f"    {n} adversarial runs")
print(f"    families found with M strictly in (1/14, 1/13): {len(inside)}")
for m,V in inside[:4]:
    print(f"      M = {m} = {float(m):.8f}   V = {V}")
print(f"    smallest M seen overall: {best_below[0]} = {float(best_below[0]):.8f}"
      f"  {'BELOW 1/14 -- LRC(14) COUNTEREXAMPLE!' if best_below[0]<LO else '(>= 1/14, no counterexample)'}")

print()
print("(2) THE KNOWN LADDER -- what values actually occur near the floor?")
from collections import Counter
C=Counter()
random.seed(3931)
for _ in range(60):
    V=sorted(random.sample(range(1,30),13))
    if reduce(gcd,V)!=1: continue
    m=M(V)
    if m<F(1,9): C[m]+=1
for m,c in sorted(C.items())[:8]:
    tag=""
    if m==LO: tag=" <-- 1/14 (tight)"
    if m==HI: tag=" <-- 1/13"
    print(f"      M = {str(m):9s} = {float(m):.8f}  count {c}{tag}")
print()
print("(3) THE 12->24 PARALLEL across n")
print("    n=13 (klein THM-1004): {1,...,12} tight at 1/13; equality case {1,...,11,24}")
print("    n=14 (my THM-1120):    {1,...,13} tight at 1/14; second tight  {1,...,11,13,24}")
print("    the SAME substitution 12 -> 24 produces the exceptional family at both n.")
for nm,V,fl in [("{1,...,12}",list(range(1,13)),F(1,13)),
                ("{1,...,11,24}",list(range(1,12))+[24],F(1,13)),
                ("{1,...,13}",list(range(1,14)),F(1,14)),
                ("{1..11,13,24}",[1,2,3,4,5,6,7,8,9,10,11,13,24],F(1,14))]:
    m=M(V)
    print(f"      {nm:16s} M = {str(m):8s} = {float(m):.8f}   floor {fl} -> {'AT FLOOR' if m==fl else 'above'}")
