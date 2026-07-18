# opus-2026-07-17-S377 -- HYP-7600: MAPPING THE TIGHT FAMILIES.
#
# A family is TIGHT if uncovered measure is EXACTLY 0 -- the danger arcs cover
# [0,1] except finitely many points, so the extremal gap is exactly 1/14 and is
# attained only on a null set.  These are the extremal objects of LRC(14).
# Known so far: {1,...,13} (classical), its dilates k*{1,...,13}, and
# {1,...,11,13,24} (S376, found by the swap 12 -> 24).
# That swap is suggestive: 24 = 2*12, and D_24 = {t : ||12t - k/2|| < 1/28},
# i.e. narrower at the D_12 centres but adding coverage at the half-points.
# So the natural first sweep is SINGLE-SPEED MULTIPLICATION.
from fractions import Fraction as F
from functools import reduce
from math import gcd
LAM=F(1,14)
def teeth01(x):
    w=LAM/x; out=[]
    for j in range(0,x+1):
        a,b=max(F(j,x)-w,F(0)), min(F(j,x)+w,F(1))
        if a<b: out.append((a,b))
    return out
def union(ivs):
    ivs=sorted(ivs); out=[]
    for a,b in ivs:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def uncovered(V):
    allv=[]
    for x in V: allv.extend(teeth01(x))
    return 1-sum(b-a for a,b in union(allv))
def tight(V): return uncovered(V)==0

base=list(range(1,14))
print("(1) SINGLE-SPEED MULTIPLICATION on {1,...,13}:  replace i by i*m")
print("    i    multipliers m (2..12) giving a TIGHT family")
found=[]
for i in base:
    good=[]
    for m in range(2,13):
        V=sorted([x for x in base if x!=i]+[i*m])
        if len(set(V))!=13: continue
        if tight(V):
            good.append(m); found.append(tuple(V))
    if good: print(f"    {i:2d}   {good}")
print(f"    total tight families from this sweep: {len(found)}")

print()
print("(2) SINGLE-SPEED REPLACEMENT by ANY value (not just multiples)")
print("    i -> r  giving TIGHT, r <= 120")
found2=[]
for i in base:
    good=[]
    for r in range(14,121):
        V=sorted([x for x in base if x!=i]+[r])
        if len(set(V))!=13: continue
        if tight(V):
            good.append(r); found2.append(tuple(V))
    if good: print(f"    {i:2d} -> {good}")
print(f"    total: {len(found2)}")

print()
print("(3) ARE THE NEW ONES DILATES OR GENUINELY NEW?")
allf=sorted(set(found)|set(found2))
for V in allf[:14]:
    g=reduce(gcd,V)
    red=tuple(x//g for x in V)
    isdil = (red==tuple(range(1,14)))
    print(f"    {str(list(V)):46s} gcd={g}  {'DILATE of {1..13}' if isdil else 'GENUINELY NEW'}")
print(f"    ... {len(allf)} total; genuinely new: "
      f"{sum(1 for V in allf if tuple(x//reduce(gcd,V) for x in V)!=tuple(range(1,14)))}")
