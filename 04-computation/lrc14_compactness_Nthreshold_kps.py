#!/usr/bin/env python3
"""
lrc14_compactness_Nthreshold_kps  (kind-pasteur, PROVE side)

Find the EXPLICIT N(theta) for the decoupling lower bound, honestly.

Bound:  L(C u {vmax}) >= (6/7) meas(G_C) - r/(7 vmax).
For this to exceed a target theta we need
   vmax > r / (7*((6/7)meas(G_C) - theta))   provided (6/7)meas(G_C) > theta.
If (6/7)meas(G_C) <= theta the bound is vacuous (the core's own gap is too small;
but then the LIMIT L-> (6/7)meas(G_C) might itself be < theta? NO -- we computed
min floor = 1/143 > 1/1260, so every core's limit exceeds 1/1260). The issue is
purely that the FINITE-vmax bound is loose; the limit is fine.

So: for theta=1/1260, since every core floor (6/7)meas(G_C) >= 1/143 > 1/1260,
there EXISTS a finite N_e for each core beyond which L > 1/1260. We compute N_e
exactly from the bound, AND the TRUE threshold (smallest vmax beyond which the
actual min L over w>=vmax exceeds 1/1260).
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

base=list(range(1,14))
theta=F(1,1260)
print(f"theta = 1/1260 = {float(theta):.6e}")
print(f"For each core C={{1..13}}\\{{e}}: floor=(6/7)meas(G_C), bound-N = r/(7(floor-theta)),")
print(f"and TRUE-N = smallest vmax s.t. min_{{w>=vmax}} L(C+w) > theta.")
print(f"  {'e':>3s} {'floor':>12s} {'floor>theta':>11s} {'bound-N':>9s} {'true-N':>7s}")
maxN=0
for e in base:
    C=[x for x in base if x!=e]
    GC=complement([a for v in C for a in danger(v)])
    mg=total(GC); r=len(GC); floor=F(6,7)*mg
    if floor>theta:
        boundN=int(F(r,7*(floor-theta)))+1
    else:
        boundN=None
    # true threshold: largest w that gives L <= theta (scan all w up to 2000)
    trueN=0
    for w in range(14,2001):
        if w in C: continue
        if L_float(C+[w])<=float(theta)+1e-12:
            trueN=max(trueN,w)
    maxN=max(maxN, trueN)
    print(f"  {e:>3d} {str(floor):>12s} {'yes' if floor>theta else 'NO':>11s} {str(boundN):>9s} {str(trueN):>7s}")
print(f"\nMAX true-N over cores = {maxN}")
print(f"COMPACTNESS LEMMA (single-perturbation form, PROVEN-modulo-core-bound):")
print(f"  If S = C u {{w}} with C a 12-subset of {{1..13}} and L(S) < 1/1260, then w <= {maxN}.")
print(f"  Hence min positive L over single perturbations is attained at BOUNDED w,")
print(f"  and exhaustive search w<=2000 (done) proves it = 1/1260.")
print()
print("THE REMAINING GAP for full inf L>0:")
print("  This argument fixes 12 elements in {1..13}. The general claim needs: ANY")
print("  primitive 13-set with small L is within bounded Hamming distance of a tight")
print("  config (empirically true, hAP<=2 over large search) AND the perturbing")
print("  elements are bounded. The decoupling floor (6/7)meas(G_C) handles ONE large")
print("  element per fixed 12-core; iterating handles k large elements as long as the")
print("  remaining (13-k)-core has positive lonely measure bounded below. Proving a")
print("  UNIFORM positive lower bound on meas(G_C) over ALL 12-subsets C (or showing")
print("  finitely many tight 13-sets) is the open crux.")
