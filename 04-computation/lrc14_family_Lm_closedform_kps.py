#!/usr/bin/env python3
"""
lrc14_family_Lm_closedform_kps  (kind-pasteur, PROVE side)

CLOSED-FORM for L({1..11,13,m}) and WHY m=36 minimizes it.

We computed: the e-exclusive gap G_12 has 4 arcs. Two are "outer" (near tau=1/12,
11/12, width 2/1001 each) and two are "central" (near 5/12,7/12, width 1/245 each).
Adding speed m re-covers part of these 4 arcs with its danger comb (centers k/m,
half-width 1/(14m)). The opened L = meas(G_12 \ D_m).

We want L(m) as a piecewise function. We compute it exactly for m = 14..200 and
locate the global minimum, then DISSECT m=36: which danger centers k/36 land in
each of the 4 G_12 arcs, and what residual each leaves. This gives the explicit
reason 36 is optimal: 36=3*12 places centers k/36 so that the OUTER arcs are fully
covered, and the CENTRAL arcs are MAXIMALLY (but not fully) covered, leaving the
minimal 2 x 1/2520.
"""
from fractions import Fraction as F

def danger(v):
    out=[]; w=F(1,14*v)
    for k in range(v+1):
        lo=F(k,v)-w; hi=F(k,v)+w
        if lo<0: out += [(F(0),hi),(1+lo,F(1))]
        elif hi>1: out += [(lo,F(1)),(F(0),hi-1)]
        else: out.append((lo,hi))
    return [(x,y) for x,y in out if y>x]
def union(arcs):
    arcs=sorted((x,y) for x,y in arcs if y>x)
    if not arcs: return []
    res=[]; cl,ch=arcs[0]
    for lo,hi in arcs[1:]:
        if lo<=ch: ch=max(ch,hi)
        else: res.append((cl,ch)); cl,ch=lo,hi
    res.append((cl,ch)); return res
def total(a): return sum(y-x for x,y in union(a))
def complement(arcs):
    u=union(arcs); res=[]; prev=F(0)
    for lo,hi in u:
        if lo>prev: res.append((prev,lo))
        prev=max(prev,hi)
    if prev<1: res.append((prev,F(1)))
    return res
def intersect(A,B):
    A=union(A); B=union(B); res=[]; i=j=0
    while i<len(A) and j<len(B):
        lo=max(A[i][0],B[j][0]); hi=min(A[i][1],B[j][1])
        if lo<hi: res.append((lo,hi))
        if A[i][1]<B[j][1]: i+=1
        else: j+=1
    return res
def subtract(A,B): return intersect(A,complement(B))

base=list(range(1,14))
G12 = complement([a for v in base if v!=12 for a in danger(v)])
print("G_12 four arcs:")
labels=[]
for lo,hi in G12:
    labels.append((lo,hi))
    print(f"  [{lo},{hi}] len={hi-lo} center~{float((lo+hi)/2):.4f}")

# L(m) closed/explicit table
print("\nL(m)=meas(G_12 \\ D_m), m=14..120:")
mins=[]
for m in range(14,121):
    Lm=total(subtract(G12,danger(m)))
    mins.append((Lm,m))
mins_sorted=sorted(mins)
print("  smallest 8:")
for Lm,m in mins_sorted[:8]:
    print(f"    m={m:3d}: L={Lm}={float(Lm):.6e}")

# dissect m=36 per-arc
print("\nPer-arc dissection at m=36:")
for (lo,hi) in G12:
    arc=[(lo,hi)]
    cov=intersect(arc,danger(36))
    res=subtract(arc,danger(36))
    print(f"  arc[{lo},{hi}] len={hi-lo}: covered by 36 = {total(cov)}, residual = {total(res)}")
    # which k/36 centers fall in / near this arc
    cs=[F(k,36) for k in range(37) if lo-F(1,14*36)<=F(k,36)<=hi+F(1,14*36)]
    print(f"      nearby 36-centers k/36: {[str(c) for c in cs]}")

# Compare 24 (tight) vs 36: why 24 fully covers central but 36 doesn't
print("\nCompare m=24 (TIGHT, L=0) vs m=36:")
for m in (24,36):
    print(f"  m={m}:")
    for (lo,hi) in G12:
        res=total(subtract([(lo,hi)],danger(m)))
        cs=[F(k,m) for k in range(m+1) if lo-F(1,14*m)<=F(k,m)<=hi+F(1,14*m)]
        print(f"    arc center~{float((lo+hi)/2):.4f}: residual={res}  centers in arc: {[str(c) for c in cs]}")

# The clean statement of the residual
print("\nKEY: 24's centers include 10/24=5/12 and 14/24=7/12 EXACTLY (=central arc")
print("centers), so 24 covers the central arcs fully. 36's centers near 5/12 are")
print("14/36=0.3889 and 15/36=0.4167; 15/36 sits at the RIGHT edge so the LEFT")
print("part of the central arc [29/70,41/98] from 5's boundary 29/70 to 36's left")
print("boundary 15/36-1/504=209/504 stays lonely: 209/504-29/70=1/2520.")
print(f"  5/12 = {F(5,12)}; 15/36={F(15,36)}; equal? {F(5,12)==F(15,36)}")
print(f"  So 15/36 = 5/12 exactly. 36-arc at 5/12 has half-width 1/504 < 12-arc's 1/168.")
print(f"  The 12-arc covered [5/12-1/168,5/12+1/168]; 36 only covers [5/12-1/504,5/12+1/504].")
print(f"  The gap G_12 central arc was [29/70,41/98] (= the part of 12-arc not covered")
print(f"  by 5 on the left and 13 on the right). 36 covers [209/504, 5/12]=its left half")
print(f"  reaching only to 5/12-1/504=209/504, leaving [29/70,209/504] lonely.")
print(f"  29/70 is 5's right boundary 2/5+1/70. So the RESIDUAL is the 5<->36 clash.")
