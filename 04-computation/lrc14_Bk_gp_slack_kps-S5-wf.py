#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_Bk_gp_slack_kps-S5-wf.py  (kind-pasteur-2026-06-18-S5)
How weak a mu_{1/7} lower bound suffices to close k>=8 via the union bound?
  Need:  mu_{1/7}(E) >= 1 - min_{|P|=13-k} meas(G_P) =: thr_k.
Report thr_k and the SLACK vs the (searched) min mu_{1/7}.  Also: a fully threshold-free
sufficient lemma -- does mu_{1/7}(E) >= 1 - meas(Bad over a single small speed) suffice trivially?
"""
import sys, itertools, random
from fractions import Fraction as F
try: sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception: pass
random.seed(101)
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
print("Required mu_{1/7} threshold thr_k = 1 - min_P meas(G_P), and slack vs searched min mu")
print("="*84)
for k in range(8,14):
    psz=13-k
    minGP=min(meas(safe_set(list(P))) for P in itertools.combinations(range(1,14),psz))
    thr_k=1-minGP
    # searched min mu_1/7 (consecutive is the min per part (1) of unionbound script; reconfirm quickly)
    minmu=meas(good_set_thr(list(range(k)),ONE7))
    cand=[tuple(range(k))]
    base=list(range(k+2))
    for d1,d2 in itertools.combinations(range(1,k+1),2):
        cand.append(tuple(x for x in base if x not in (d1,d2)))
    for _ in range(800):
        sp=random.choice([k-1,k,k+1,k+3,2*k])
        e=tuple([0]+sorted(random.sample(range(1,sp+1),min(k-1,sp))))
        cand.append(e)
    for E in set(cand):
        if len(set(E))==k and 0 in E:
            m=meas(good_set_thr(list(E),ONE7))
            if m<minmu: minmu=m
    slack=minmu-thr_k
    print(f"  k={k:2d}: thr_k=1-{minGP}={float(thr_k):.4f} | min mu_1/7={float(minmu):.4f} | SLACK={float(slack):.4f} ({'OK' if slack>0 else 'FAILS'})")
print()
print("INTERPRETATION: any spread-bound proving mu_{1/7}(E) >= thr_k closes k>=8.")
print("thr_8 = 1 - 2243/5880 = 3637/5880 ~ 0.6185 is the binding requirement; the searched")
print("min mu_1/7 ~ 0.940 clears it by ~0.32. A loose discrepancy bound suffices.")
print("\nDONE.")
