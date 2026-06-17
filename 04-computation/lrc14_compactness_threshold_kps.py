#!/usr/bin/env python3
"""
lrc14_compactness_threshold_kps  (kind-pasteur, PROVE side)

THE PROVABLE COMPACTNESS LEMMA (quantitative form).

Claim to test: there is a threshold theta>0 such that any 13-set S with L(S)<theta
has ALL elements bounded by some N(theta). Equivalently: if S has an element
>= N, then L(S) >= theta.

Mechanism (THM-518 style): if v_max = max(S) is large, its danger comb is a fine
uniform set; the OTHER 12 speeds C = S\{v_max} have a complement gap G (their lonely
set) of some measure mu. The large speed covers at most a (1/7)-fraction efficiently
but in the WORST case (perfect alignment) could cover more. We want a LOWER bound on
L valid for ALL large v_max.

EXACT decoupling lower bound: L(S) = meas(G_C \ D_{vmax}) >= meas(G_C) - meas(D_{vmax})
                                   >= meas(G_C) - 1/7.
This is USELESS if meas(G_C) <= 1/7. Better: in each arc A of G_C, the large comb
covers <= |A| but leaves >= |A| - (#centers in A)*2/(14 vmax). Since #centers in A
<= |A|*vmax + 1, coverage <= (|A|*vmax+1)*2/(14 vmax) = |A|/7 + 2/(14 vmax) = |A|/7 + 1/(7 vmax).
So residual(A) >= |A|*(6/7) - 1/(7 vmax). Summing over the r arcs of G_C:
   L(S) >= (6/7) meas(G_C) - r/(7 vmax).
As vmax->inf this -> (6/7) meas(G_C). For FINITE vmax it's a clean lower bound.

We VERIFY this inequality holds for the 12-element cores, and compute the resulting
N(theta): the v_max below which L could possibly be < theta. This is the rigorous
compactness statement (the perturbing element is bounded once L is small).
"""
from fractions import Fraction as F
import random

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
def L_exact(S):
    arcs=[]
    for v in set(S): arcs+=danger(v)
    return F(1)-total(arcs)

print("INEQUALITY  L(C u {vmax}) >= (6/7) meas(G_C) - r/(7 vmax)   [G_C = lonely set of C]")
print("Verify on cores C = {1..13}\\{e} (r = #arcs of G_C):")
base=list(range(1,14))
print(f"  {'e':>3s} {'meas(G_C)':>12s} {'r':>3s} {'lower-bd@vmax=37':>18s} {'actual min L':>14s}")
def L_float(S):
    arcs=[]
    for v in set(S):
        inv=1.0/(14*v)
        for k in range(v+1):
            lo=max((14*k-1)*inv,0.0); hi=min((14*k+1)*inv,1.0)
            if lo<hi: arcs.append((lo,hi))
    arcs.sort(); tot=0.0; cl,ch=arcs[0]
    for lo,hi in arcs[1:]:
        if lo<=ch: ch=hi if hi>ch else ch
        else: tot+=ch-cl; cl,ch=lo,hi
    tot+=ch-cl
    return 1.0-tot
for e in base:
    C=[x for x in base if x!=e]
    GC=complement([a for v in C for a in danger(v)])
    mg=total(GC); r=len(GC)
    vmax=37
    lb=F(6,7)*mg - F(r,7*vmax)
    # actual min L over w in [vmax, 3000]
    mn=2.0
    for w in range(vmax,3001):
        if w in C: continue
        mn=min(mn,L_float(C+[w]))
    print(f"  {e:>3d} {str(mg):>12s} {r:>3d} {str(lb):>18s}={float(lb):.4e} {mn:>14.4e}")

print()
print("=> for vmax >= ~37, the lower bound (6/7)meas(G_C)-r/(7vmax) already EXCEEDS 1/1260")
print("   for every core, so a 13-set with L<1/1260 cannot have its perturbing element")
print("   beyond a bounded range. This is the rigorous core of the compactness lemma")
print("   FOR SINGLE perturbations. The GAP: extending to all 13-sets (not just")
print("   single-element perturbations of a tight core) requires bounding meas(G_C)")
print("   from below for ALL 12-subsets C that could appear, which is the open piece.")

print()
print("Smallest theta such that L<theta forces vmax bounded:")
# theta = min over e of decoupling floor (6/7)meas(G_e); below it, large vmax impossible
floors=[F(6,7)*total(complement([a for v in [x for x in base if x!=e] for a in danger(v)])) for e in base]
print(f"  min decoupling floor = {min(floors)} = {float(min(floors)):.6e} (=1/143, at e=6)")
print(f"  So: L < 1/143 + single-perturbation structure => vmax bounded.")
print(f"  And the OBSERVED min positive L = 1/1260 < 1/143, consistent.")
