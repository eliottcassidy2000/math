#!/usr/bin/env python3
"""
lrc14_covering_min_hunt_klein_S267.py
=====================================
klein-2026-07-12-S267

THE QUESTION: pin the COVERING-MINIMUM = min M(S) over primitive COVERING (divisor-
complete) 13-sets. Three conflicting fleet claims — kps cont.51 "DC floor 1/12",
boxeph-S20 "1/13", klein-S266 deep well {1..12,182}=14/183 — are all box/search
artifacts of different ranges. The true object is the covering-min; if it exceeds 1/14
the LRC(14) covering case is proved.

M(S) = max_t min_i ||v_i t||, exact via THM-668 pair-sum ruler (M = p/q, q = v_i+v_j).
covering(S) := a multiple of every d in {2..14} is present.

Reference values: 1/14=0.07143 (tight, non-covering) < 2/27=0.07407 (H5 lower bracket)
< 14/183=0.07650 (deep well {1..12,182} = n/Phi6, Phi6(14)=183) < 1/13=0.07692
< 1/12=0.08333.  H5 reportedly brackets covering-min in [2/(2n-1), n/Phi6]=[2/27,14/183].

This script:
 (1) verifies the three claimed floor families + the deep well (exact M).
 (2) HUNTS for covering families with M < 14/183 (toward 2/27 or 1/14) using STRUCTURED
     shapes (deep-well variants, Farey rungs, spread sets, base+outlier), NOT hill-climb
     (which misses structured extremals — MISTAKE-101/138).
 (3) reports the minimum M found and where it sits vs the H5 bracket [2/27, 14/183].
"""
import math
from fractions import Fraction
from itertools import combinations

def dist0(r,q): return min(r,q-r)
def exact_M(v):
    n=len(v); best=Fraction(0); argq=None
    qs=set()
    for a in range(n):
        for b in range(a,n): qs.add(v[a]+v[b])
    for q in sorted(qs):
        mx=0
        for p in range(1,q):
            mn=q
            for vl in v:
                d=dist0((vl*p)%q,q)
                if d<mn:
                    mn=d
                    if mn<=mx: break
            if mn>mx: mx=mn
        val=Fraction(mx,q)
        if val>best: best=val; argq=q
    return best,argq
def primitive(v):
    g=0
    for x in v: g=math.gcd(g,x)
    return g==1
def covering(v):
    return len(set(v))==13 and all(any(x%d==0 for x in v) for d in range(2,15))

one14=Fraction(1,14); two27=Fraction(2,27); dw=Fraction(14,183); one13=Fraction(1,13)

print("="*72)
print("(1) THE CLAIMED FLOOR FAMILIES + THE DEEP WELL (exact M)")
print("="*72)
claims = [
    ("kps cont.51 floor {1,2,3,4,10..18}", [1,2,3,4,10,11,12,13,14,15,16,17,18]),
    ("boxeph-S20 floor 2*{1..12}u{13}", sorted([2*i for i in range(1,13)]+[13])),
    ("deep well {1..12,182} (=n/Phi6)", list(range(1,13))+[182]),
    ("deep well x2 {1..12,364}", list(range(1,13))+[364]),
    ("{1..11,13,182}", list(range(1,12))+[13,182]),
]
for nm,v in claims:
    if not covering(v):
        print(f"  {nm:38s} NOT covering (skip)"); continue
    M,q=exact_M(v)
    print(f"  {nm:38s} M={str(M):>8s}={float(M):.5f} (q={q}) prim={primitive(v)}")
print(f"  refs: 1/14={float(one14):.5f} 2/27={float(two27):.5f} 14/183={float(dw):.5f} 1/13={float(one13):.5f}")

print()
print("="*72)
print("(2) HUNT covering families with M < 14/183 (structured shapes)")
print("="*72)
best=Fraction(1); bestv=None; below_dw=[]
def consider(v):
    global best,bestv
    v=sorted(set(v))
    if len(v)!=13 or not primitive(v) or not covering(v): return
    M,q=exact_M(v)
    if M<best: best=M; bestv=(list(v),q)
    if M<dw: below_dw.append((M,q,list(v)))

# (a) deep-well shape {1..12, x} — x must supply 13 AND 14 => 182|x
for m in range(1,12):
    consider(list(range(1,13))+[182*m])
# (b) {1..11, a, b} with a supplying 13, b supplying 14 (a=13k, b=14j), spread
for a in [13,26,39,52,65,78,91,104,117,130,143,156,169,182]:
    for b in [14,28,42,56,70,84,98,112,126,140,154,168,182]:
        consider(list(range(1,12))+[a,b])
# (c) Farey-rung deep wells {1..n-2, n(n-1)} and neighbors {1..12, 182+-d}
for d in range(-6,7):
    consider(list(range(1,13))+[182+d])
# (d) spread beater search: base {1..k} truncated + spread covering completion
#     try {1,2,3, a,b,c,...} minimizing M over covering completions <= cap
import random
random.seed(267)
cap=260
# ensure covering: force witnesses for the primes/prime-powers 8,9,5,7,11,13, and 14
witness_pools = {
    8:[8,16,24,32,40,48,56,64,72,80,88,96,104,112,120,128,136,144,152,160,168,176,184,192,200,208,216,224,232,240,248,256],
    9:[9,18,27,36,45,54,63,72,81,90,99,108,117,126,135,144,153,162,171,180,189,198,207,216,225,234,243,252],
    5:[5,10,15,20,25,30,35,40,45,50,55,60],
    7:[7,14,21,28,35,42,49,56],
    11:[11,22,33,44,55,66,77,88,99,110,121,132],
    13:[13,26,39,52,65,78,91,104,117,130,143,156,169,182],
    14:[14,28,42,56,70,84,98,112,126,140,154,168,182,196,210,224,238,252],
    6:[6,12,18,24,30,36,42,48],
    10:[10,20,30,40,50,60],
    12:[12,24,36,48,60,72],
}
for _ in range(15000):
    v=set([1])
    for d in (8,9,5,7,11,13,14):     # the binding divisor demands
        v.add(random.choice(witness_pools[d]))
    while len(v)<13:
        v.add(random.randint(2,cap))
    v=sorted(v)[:13]
    if len(v)==13:
        consider(v)
# (e) targeted: covering families whose M sits at a low Farey rung r/q, q a pair-sum near 183
for base_top in range(150,260):
    consider(list(range(1,12))+[13,base_top] if (base_top%14==0) else list(range(1,13))+[base_top] if base_top%182==0 else list(range(1,12))+[13,14])

print(f"  minimum M over all covering families tried: {best} = {float(best):.5f}")
if bestv:
    print(f"    achiever: {bestv[0]} (ruler q={bestv[1]})")
print(f"  families found with M < 14/183 = {float(dw):.5f}: {len(below_dw)}")
seen=set()
for M,q,v in sorted(below_dw):
    if M in seen: continue
    seen.add(M)
    loc = "in H5 bracket (2/27,14/183)" if two27<M<dw else ("= 2/27 LOWER edge!" if M==two27 else ("< 2/27 !!" if M<two27 else ""))
    print(f"    M={str(M):>8s}={float(M):.5f} (q={q}) {loc}  e.g.{v}")
    if len(seen)>=15: break

print()
print("="*72)
print("(3) VERDICT vs H5 bracket [2/27, 14/183]")
print("="*72)
print(f"  covering-min (searched) = {best} = {float(best):.5f}")
print(f"  H5 bracket: [2/27={float(two27):.5f}, 14/183={float(dw):.5f}]")
if best >= dw:
    print(f"  => 14/183 (deep well) is the covering-min in this search; H5 UPPER bound tight.")
elif best > two27:
    print(f"  => a BEATER below the deep well exists, inside the H5 bracket. covering-min < 14/183.")
elif best == two27:
    print(f"  => the covering-min reaches the H5 LOWER edge 2/27.")
else:
    print(f"  => covering-min BELOW 2/27 -- would break H5's lower bracket! (recheck)")
print(f"  margin over 1/14: {best} - 1/14 = {best-one14} = {float(best-one14):.5f}  ({'comfortable' if best-one14>Fraction(1,400) else 'TIGHT'})")
print("\ndone.")
