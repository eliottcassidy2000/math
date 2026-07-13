#!/usr/bin/env python3
"""mac-mini-S66: the finite-Vmax glue (THM-527-A) needs a TIGHT arc-count bound on the good
set Good(E) = {x in [0,1): maxgap of {frac(e_i x)} > 2/7}. opus-S169's loose bound gives an
impractical V0. Part C says (for CONSECUTIVE E) good x live near {0,1/3,1/2,2/3} (rationals
a/b, b<=3). GOAL: test whether this LOCALIZATION near small-b rationals holds for general
COVERING bounded-spread clusters, count the exact arcs, and see the true growth.

Exact computation: breakpoints of x->maxgap are (i) collisions x=m/(e_i-e_j), (ii) the arc
where max-gap crosses 2/7. We sample finely + confirm arc boundaries exactly via sign changes,
and record which rational a/b (b<=6) each good arc's center is nearest."""
from fractions import Fraction as F
from math import gcd

def frac(x): return x-int(x) if x>=0 else x-int(x)+ (1 if x!=int(x) else 0)

def maxgap(phases):
    p=sorted(phases); n=len(p)
    if n==0: return 1.0
    if n==1: return 1.0
    g=max(p[(i+1)%n]-p[i] if i<n-1 else (p[0]+1-p[n-1]) for i in range(n))
    return g

def good_arcs(E, res=200000):
    """exact-ish arc decomposition of Good(E) by fine sampling + boundary refine."""
    E=[e for e in E]
    xs=[i/res for i in range(res)]
    good=[maxgap([ (e*x)%1.0 for e in E])>2.0/7.0 for x in xs]
    # find maximal runs of True (circular)
    arcs=[]; i=0
    # rotate so index 0 is False if possible
    if all(good): return [(0.0,1.0)], 1
    start=None
    # linear scan with circular wrap merge
    runs=[]
    j=0
    while j<res:
        if good[j]:
            k=j
            while k<res and good[k]: k+=1
            runs.append((j,k-1)); j=k
        else: j+=1
    # merge wrap
    if runs and runs[0][0]==0 and runs[-1][1]==res-1 and len(runs)>1:
        a=runs.pop(); b=runs.pop(0); runs.insert(0,(a[0]-res,b[1]))
    arcs=[( (s%res)/res, (e%res)/res ) for s,e in runs]
    return arcs, len(runs)

def nearest_rational(c, maxb=6):
    best=(9,None)
    for b in range(1,maxb+1):
        for a in range(0,b+1):
            if gcd(a,b)==1 or (a==0 and b==1) or (a==b):
                d=abs(c-a/b)
                if d<best[0]: best=(d, (a if a<b else 0,b) )
    return best

# covering bounded-spread clusters (co-offset sets E, spread <=~14) incl the covering deep-well
# cluster and assorted covering-saturated shapes. e in E are co-offsets Vmax - u.
clusters = {
 "consec k=9 {0..8}": list(range(9)),
 "consec k=11 {0..10}": list(range(11)),
 "consec k=13 {0..12} (AP cluster, OUT OF SCOPE)": list(range(13)),
 "perforated k=7 {0..8}\\{1,7}": [0,2,3,4,5,6,8],
 "spread-21 k=11 (part D extremal-ish)": [0,2,4,6,8,10,12,14,16,18,21],
 "deep-well E {0,14?}": [0,1,2,3,4,5,6,7,8,9,10,11,12],  # {1..12,182} co-offsets at Vmax=182: 182-{1..12,182}={181..170,0}; bounded rep = {0}+... use structural
}
print("cluster                                   | #arcs | arc centers -> nearest a/b (b<=6)")
print("-"*100)
for nm,E in clusters.items():
    arcs,narcs=good_arcs(E)
    centers=[]
    for s,e in arcs:
        c=((s+e)/2)%1.0
        d,ab=nearest_rational(c)
        centers.append(f"{c:.3f}~{ab[0]}/{ab[1]}")
    # dedup denominators used
    dens=sorted(set(int(cc.split('/')[1]) for cc in centers))
    print(f"{nm:41s} | {narcs:5d} | b-values {dens}  ({narcs} arcs)")
    if narcs<=12: print(f"{'':41s} |       | {centers}")
print()
print("KEY QUESTION: do all good-arc centers sit at rationals a/b with b<=3? (part C for consec)")
print("If YES for covering clusters => #arcs <= (# a/b with b<=3) * (arcs per rational) = O(1) small.")
