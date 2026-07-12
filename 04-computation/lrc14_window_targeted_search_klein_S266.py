#!/usr/bin/env python3
"""
lrc14_window_targeted_search_klein_S266.py
==========================================
klein-2026-07-12-S266

The [1..20] scan (companion script) found the window (1/14, 2/27) EMPTY -- but only
covers pair-sum rulers q <= 40. Since M = d/q on a pair-sum ruler, the Farey fraction
3/41 = 0.07317 lies INSIDE (1/14, 2/27) and needs q=41 (Vmax >= 21). So the empty-window
claim must be tested at LARGER Vmax. This does a TARGETED search of the small-M region
(near-AP / near-dilate families, where all small-M families live) up to Vmax ~ 80, hunting
for ANY primitive 13-set with M strictly in (1/14, 2/27).

If none is found under real effort, that SUPPORTS opus-S246 (empty window => rigidity
[M<2/27 => interval]).  If one is found, it REFINES the claim.

Families searched (the small-M region):
  A. {1..12} u {x},           x = 13..80            (single far element)
  B. {1..11} u {x,y},         13 <= x<y <= 80        (two far elements)
  C. g*{1..12} u {delta},     g=2..4, delta swept    (compressed near-dilates)
  D. {1..13} with j in {1,2,3} elements shifted up by t (near-AP perturbations)
  E. targeted d=3,q=41 constructions: 13 residues in [3,38] mod 41, lifted to integers
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

one14=Fraction(1,14); two27=Fraction(2,27)
def inwin(M): return one14 < M < two27

hits=[]; nseen=0; nsmall=0
small=[]   # (M, family) for M <= 0.078 (just above the window top), to see what's near

def consider(v):
    global nseen,nsmall
    v=sorted(set(v))
    if len(v)!=13 or not primitive(v): return
    nseen+=1
    M,q=exact_M(v)
    if inwin(M):
        hits.append((M,q,list(v)))
    if M <= Fraction(78,1000):
        nsmall+=1
        small.append((M,q,list(v)))

# A. {1..12} u {x}
for x in range(13,81):
    consider(list(range(1,13))+[x])
# B. {1..11} u {x,y}
for x in range(12,81):
    for y in range(x+1,81):
        consider(list(range(1,12))+[x,y])
# C. g*{1..12} u {delta}
for g in (2,3,4):
    base=[g*i for i in range(1,13)]
    for delta in range(1, 81):
        if delta % g == 0: continue
        consider(base+[delta])
# D. {1..13} with 1..3 top elements shifted up by t
base=list(range(1,14))
for t in range(1,60):
    consider(base[:12]+[13+t])                      # shift top
    for t2 in range(1,40):
        consider(base[:11]+[12+t,13+t+t2])          # shift top two
# E. targeted d=3, q=41: pick 13 residues in [3,38] mod 41 that are "spread", lift
import random
random.seed(266)
q=41
band=list(range(3,39))
for trial in range(20000):
    res=random.sample(band,13)
    # lift residues to distinct positive integers (res, or res+41, ...) keeping small
    v=sorted(res)
    consider(v)
    # also a spread lift
    v2=sorted(set(res[:9]+[r+41 for r in res[9:]]))
    if len(v2)==13: consider(v2)

print("="*70)
print(f"TARGETED WINDOW SEARCH (1/14, 2/27) = ({float(one14):.5f}, {float(two27):.5f})")
print("="*70)
print(f"primitive 13-sets evaluated: {nseen}")
print(f"families with M in the OPEN window (1/14, 2/27): {len(hits)}")
if hits:
    for M,q,v in sorted(hits)[:20]:
        print(f"  !! M={M}={float(M):.5f} (ruler q={q})  {v}")
    print("  => WINDOW IS NOT EMPTY at larger Vmax -- opus-S246's threshold needs refinement.")
else:
    print("  NONE. The window stays EMPTY under targeted near-AP/near-dilate search up to Vmax~80.")
    print("  => supports opus-S246 [M<2/27 => interval]; the d=3,q=41 candidate 3/41 is NOT realized")
    print("     as a global max (such families clear at a higher d/q, pushing M >= 2/27).")

print(f"\nSmallest M-values found (M <= 0.078), the boundary of the window from above:")
seen=set()
for M,q,v in sorted(small):
    if M in seen: continue
    seen.add(M)
    tag = "  <-- IN WINDOW" if inwin(M) else (" = 1/14" if M==one14 else (" = 2/27 (window top)" if M==two27 else ""))
    print(f"  M={str(M):>8s}={float(M):.5f} (q={q}) e.g. {v}{tag}")
    if len(seen)>=12: break
print("\ndone.")
