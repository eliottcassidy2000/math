#!/usr/bin/env python3
"""
lrc14_exact_floor_macmini_0618s1.py  (mac-mini-2026-06-18-S1)

EXACT (rational) computation of the LRC(14) S3 residual floor under the validated
reformulation:  good x  <=>  x in G_P  AND  maxgap{frac(e_i x)} > 2/7.
Worst cluster (verified) = consecutive offsets E={0,1,...,k-1} (three-distance set).

(I)  EXACT pure measure mu_k = meas{x: maxgap{0,x,...,(k-1)x} > 2/7}, k=3..13, via the
     collision-breakpoint method (cyclic order constant between x=m/d collisions; each gap
     affine; {gap>2/7} has rational endpoints). Fully exact.
(II) EXACT floor: c0 = min over k=3..13 and P subset {1..13}, |P|=13-k, of
     meas(GoodSet_k INTERSECT G_P)  (both exact rational arc-unions).
(III) EXHAUSTIVE consecutive-extremality check, k=4,5: over all E (0 in E, |E|=k, max<=11),
     is consecutive {0..k-1} the minimizer of mu? (uses exact pure measure.)
"""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
H=F(1,14); TWO7=F(2,7)

# ---------- exact safe set of P (subset of small speeds) ----------
def danger_arcs(u,h=H):
    iv=[]
    for j in range(u):
        c=F(j,u); a=(c-h/u)%1; b=(c+h/u)%1
        if a<b: iv.append((a,b))
        else: iv.append((a,F(1))); iv.append((F(0),b))
    return iv
def merge(iv):
    iv=sorted(iv); out=[]
    for a,b in iv:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def safe_set(A,h=H):
    if not A: return [(F(0),F(1))]
    dz=merge([iv for u in A for iv in danger_arcs(u,h)]); safe=[]; prev=F(0)
    for a,b in dz:
        if a>prev: safe.append((prev,a))
        prev=max(prev,b)
    if prev<1: safe.append((prev,F(1)))
    return safe
def meas(arcs): return sum((b-a for a,b in arcs), F(0))
def intersect(A,B):
    out=[]; i=j=0
    while i<len(A) and j<len(B):
        lo=max(A[i][0],B[j][0]); hi=min(A[i][1],B[j][1])
        if lo<hi: out.append((lo,hi))
        if A[i][1]<B[j][1]: i+=1
        else: j+=1
    return out

# ---------- EXACT good set for offset multiset E: {x: maxgap{frac(e x)} > 2/7} ----------
def good_set_exact(E):
    E=sorted(set(E)); k=len(E)
    # breakpoints: collisions x=m/d for differences d=e-e' (>0), m=0..d ; also where any
    # gap crosses 2/7 we handle analytically inside each cyclic-order cell, so cell endpoints
    # only need collision points (order changes) + 0,1.
    diffs=set()
    for a in range(k):
        for b in range(a+1,k):
            diffs.add(E[b]-E[a])
    bps=set([F(0),F(1)])
    for d in diffs:
        for m in range(0,d+1):
            bps.add(F(m,d))
    bps=sorted(x for x in bps if 0<=x<=1)
    good=[]
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        # cyclic order of points e*xm mod 1 (fixed on this cell); record which e at each slot
        pts=sorted(((E[t]*xm)%1, E[t]) for t in range(k))
        order=[e for _,e in pts]
        # gap between consecutive (incl wrap): gap_i(x)= frac(e_{i+1} x) - frac(e_i x) but with
        # fixed integer wrap chosen at xm. Represent each frac as e*x - floor(e*xm) on the cell.
        fl=[ (E[t]*xm).__floor__() if False else int((order_e*xm)//1) for order_e in order]  # placeholder
        # compute floors at xm for the ordered e's:
        floors=[int((e*xm)//1) for e in order]
        # value_i(x) = e_i x - floors_i  (the actual position in [0,1) on this cell, affine in x)
        # gaps: for i=0..k-1, g_i(x) = value_{i+1} - value_i, with last wrapping (+1).
        for idx in range(k):
            e_cur=order[idx]; f_cur=floors[idx]
            if idx<k-1:
                e_nx=order[idx+1]; f_nx=floors[idx+1]; wrap=F(0)
            else:
                e_nx=order[0];   f_nx=floors[0];   wrap=F(1)
            # g(x) = (e_nx x - f_nx + wrap) - (e_cur x - f_cur) = (e_nx-e_cur) x + (f_cur - f_nx + wrap)
            A=F(e_nx-e_cur); Cc=F(f_cur-f_nx)+wrap
            # solve g(x) > 2/7 on [x0,x1]
            # A x + Cc > 2/7
            if A==0:
                if Cc>TWO7: good.append((x0,x1))
                continue
            xb=(TWO7-Cc)/A
            if A>0:
                lo=max(x0,xb); hi=x1
            else:
                lo=x0; hi=min(x1,xb)
            if lo<hi: good.append((lo,hi))
    return merge(good)

print("="*80)
print("(I) EXACT pure measure mu_k (consecutive cluster {0..k-1})")
print("="*80)
mu={}
for k in range(3,14):
    g=good_set_exact(list(range(k))); mu[k]=meas(g)
    print(f"   k={k:2d}: mu_k = {mu[k]} = {float(mu[k]):.6f}   (#arcs={len(g)})")

print("\n"+"="*80)
print("(II) EXACT floor c0 = min over k, P subset {1..13} (|P|=13-k) of meas(Good_k ∩ G_P)")
print("="*80)
overall=(F(2),None,None)
for k in range(3,14):
    Gk=good_set_exact(list(range(k)))
    psz=13-k
    best=(F(2),None)
    for P in itertools.combinations(range(1,14),psz):
        val=meas(intersect(Gk, safe_set(list(P))))
        if val<best[0]: best=(val,P)
    print(f"   k={k:2d}: min rho* = {best[0]} = {float(best[0]):.6f}  at P={best[1]}")
    if best[0]<overall[0]: overall=(best[0],best[1],k)
print(f"\n   EXACT FLOOR c0 = {overall[0]} = {float(overall[0]):.6f}")
print(f"        at P={overall[1]}, consecutive cluster k={overall[2]}  (|P|={13-overall[2]})")

print("\n"+"="*80)
print("(III) EXHAUSTIVE consecutive-extremality: is {0..k-1} the min-mu shape? (k=4,5)")
print("="*80)
for k in (4,5):
    cons=mu[k]; viol=[]; checked=0
    for rest in itertools.combinations(range(1,12), k-1):
        E=[0]+list(rest); checked+=1
        m=meas(good_set_exact(E))
        if m<cons-F(1,10**9): viol.append((E,m))
    print(f"   k={k}: consecutive mu={float(cons):.5f}; checked {checked} shapes (max<=11); "
          f"violations (mu<consecutive): {len(viol)}")
    for E,m in viol[:5]: print(f"       E={E}: mu={float(m):.5f}")
print("\nDONE.")
