#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_Bk_gp_unionbound_kps-S5-wf.py  (kind-pasteur-2026-06-18-S5)
Companion to the gp-intersection-uniform angle: the k>=8 union bound and the 1/7 spread check.

  rho*_{1/7}(P,E) >= meas(G_P) + mu_{1/7}(E) - 1   (inclusion-exclusion / union bound).
Closes k>=8 IF  min_P meas(G_P) + min_E mu_{1/7}(E) > 1.   k<=7: pigeonhole mu_{1/7}=1.
"""
import sys, itertools, random
from fractions import Fraction as F
try: sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception: pass
random.seed(11)
H=F(1,14); ONE7=F(1,7)

def merge(iv):
    iv=sorted(iv); out=[]
    for a,b in iv:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def meas(arcs): return sum((b-a for a,b in arcs),F(0))
def complement(arcs):
    arcs=merge(arcs); out=[]; prev=F(0)
    for a,b in arcs:
        if a>prev: out.append((prev,a))
        prev=max(prev,b)
    if prev<1: out.append((prev,F(1)))
    return out
def danger_arcs(u,h=H):
    iv=[]
    for j in range(u):
        c=F(j,u); a=(c-h/u)%1; b=(c+h/u)%1
        if a<b: iv.append((a,b))
        else: iv.append((a,F(1))); iv.append((F(0),b))
    return iv
def safe_set(P,h=H):
    if not P: return [(F(0),F(1))]
    return complement(merge([iv for u in P for iv in danger_arcs(u,h)]))
def good_set_thr(E,thr):
    E=sorted(set(E)); k=len(E); diffs=set()
    for a in range(k):
        for b in range(a+1,k): diffs.add(E[b]-E[a])
    bps=set([F(0),F(1)])
    for d in diffs:
        for m in range(0,d+1): bps.add(F(m,d))
    bps=sorted(x for x in bps if 0<=x<=1); good=[]
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        pts=sorted(((E[t]*xm)%1,E[t]) for t in range(k))
        order=[e for _,e in pts]; floors=[int((e*xm)//1) for e in order]
        for idx in range(k):
            e_cur=order[idx]; f_cur=floors[idx]
            if idx<k-1: e_nx=order[idx+1]; f_nx=floors[idx+1]; wrap=F(0)
            else:       e_nx=order[0];     f_nx=floors[0];     wrap=F(1)
            A=F(e_nx-e_cur); Cc=F(f_cur-f_nx)+wrap
            if A==0:
                if Cc>thr: good.append((x0,x1))
                continue
            xb=(thr-Cc)/A
            if A>0: lo=max(x0,xb); hi=x1
            else:   lo=x0;        hi=min(x1,xb)
            if lo<hi: good.append((lo,hi))
    return merge(good)

print("="*84)
print("(1) Does consecutive E minimize mu_{1/7} at each k>=8? (structured + random search)")
print("="*84)
spread_min_consec=True
for k in range(8,14):
    muc=meas(good_set_thr(list(range(k)),ONE7)); below=None; minm=muc
    cand=[]
    base=list(range(k+1))
    for drop in range(1,k): cand.append(tuple(x for x in base if x!=drop))
    base=list(range(k+2))
    for d1,d2 in itertools.combinations(range(1,k+1),2):
        cand.append(tuple(x for x in base if x not in (d1,d2)))
    for _ in range(2500):
        sp=random.choice([k-1,k,k+1,k+2,k+3,k+5,2*k,3*k])
        e=tuple([0]+sorted(random.sample(range(1,sp+1),min(k-1,sp))))
        cand.append(e)
    for E in set(cand):
        if len(set(E))!=k or 0 not in E: continue
        m=meas(good_set_thr(list(E),ONE7))
        if m<minm: minm=m; below=E
    if below is None:
        print(f"  k={k:2d}: consecutive mu_1/7={muc}={float(muc):.4f} IS the min found")
    else:
        spread_min_consec=False
        print(f"  k={k:2d}: consec mu={float(muc):.4f}; LOWER E={below} mu={float(minm):.4f}")
print(f"  consecutive-minimizes-mu_1/7 (search): {spread_min_consec}")

print("\n"+"="*84)
print("(2) Per-pair union LB = meas(G_P)+mu_1/7(consecutive E)-1, min over admissible P, k>=8")
print("="*84)
gmin=(F(10),None,None)
for k in range(8,14):
    psz=13-k; muc=meas(good_set_thr(list(range(k)),ONE7)); mn=(F(10),None)
    for P in itertools.combinations(range(1,14),psz):
        lb=meas(safe_set(list(P)))+muc-1
        if lb<mn[0]: mn=(lb,P)
    print(f"  k={k:2d}: min union-LB = {str(mn[0]):>14s} = {float(mn[0]):.4f}  at P={mn[1]}")
    if mn[0]<gmin[0]: gmin=(mn[0],mn[1],k)
print(f"\n  GLOBAL min per-pair union-LB (k>=8) = {gmin[0]} = {float(gmin[0]):.4f}  (k={gmin[2]} P={gmin[1]})")
print(f"  => for k>=8, rho*_1/7 >= {gmin[0]} > 0 PROVIDED consecutive minimizes mu_1/7 (the 1/7 spread-bound).")

print("\n"+"="*84)
print("(3) The clean combined floor:  c0' = min( m_P[k<=7 pigeonhole=meas G_P],  union-LB[k>=8] )")
print("="*84)
# k<=7: rho*_1/7 = meas(G_P) >= m_P(|P|) ; the smallest such is at k=7 (|P|=6): 3029/10780
mP7=None
for P in itertools.combinations(range(1,14),6):
    m=meas(safe_set(list(P)))
    if mP7 is None or m<mP7: mP7=m
# but the global admissible m_P is at |P|=10 (k=3): smaller. for k<=7 the relevant min is at |P|=6.
# Actually rho*_1/7 = meas(G_P) and we take min over the k<=7 admissible P-sizes (|P|=6..10):
mP_k_le_7=min(min(meas(safe_set(list(P))) for P in itertools.combinations(range(1,14),psz)) for psz in range(6,11))
print(f"  k<=7 branch: min meas(G_P) over |P| in 6..10 = {mP_k_le_7} = {float(mP_k_le_7):.6f}")
print(f"  k>=8 branch: union-LB floor = {gmin[0]} = {float(gmin[0]):.4f}")
c0p=min(mP_k_le_7, gmin[0])
print(f"  => combined c0' = {c0p} = {float(c0p):.6f}  (the global-witness uniform floor candidate)")
print("\nDONE.")
