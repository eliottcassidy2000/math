#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_Bk_mu17_crudefloor_kps-S5-wf.py  (kind-pasteur-2026-06-18-S5)
Is there a CLEAN provable lower bound for mu_{1/7}(E) (any k<=13) that already clears thr_k?
We test the candidate: mu_{1/7}(E) >= mu_{1/7}(consecutive_k) (the 1/7 spread bound), exhaustively
for small k and broadly for large k, AND report the gap to thr_k. The point: even a CRUDE
universal floor (e.g. the Three-Gap / Steinhaus bound) clears thr_k with room to spare.
"""
import itertools, random
from fractions import Fraction as F
random.seed(909)
ONE7=F(1,7); H=F(1,14)
def merge(iv):
    iv=sorted(iv);out=[]
    for a,b in iv:
        if out and a<=out[-1][1]:out[-1]=(out[-1][0],max(out[-1][1],b))
        else:out.append((a,b))
    return out
def meas(arcs): return sum((b-a for a,b in arcs),F(0))
def complement(arcs):
    arcs=merge(arcs);out=[];prev=F(0)
    for a,b in arcs:
        if a>prev:out.append((prev,a))
        prev=max(prev,b)
    if prev<1:out.append((prev,F(1)))
    return out
def danger_arcs(u,h=H):
    iv=[]
    for j in range(u):
        c=F(j,u);a=(c-h/u)%1;b=(c+h/u)%1
        if a<b:iv.append((a,b))
        else:iv.append((a,F(1)));iv.append((F(0),b))
    return iv
def safe_set(P,h=H):
    if not P:return [(F(0),F(1))]
    return complement(merge([iv for u in P for iv in danger_arcs(u,h)]))
def good17(E):
    E=sorted(set(E));k=len(E);diffs=set()
    for a in range(k):
        for b in range(a+1,k):diffs.add(E[b]-E[a])
    bps=set([F(0),F(1)])
    for d in diffs:
        for m in range(0,d+1):bps.add(F(m,d))
    bps=sorted(x for x in bps if 0<=x<=1);good=[]
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0:continue
        xm=(x0+x1)/2
        pts=sorted(((E[t]*xm)%1,E[t]) for t in range(k))
        order=[e for _,e in pts];floors=[int((e*xm)//1) for e in order]
        for idx in range(k):
            e_cur=order[idx];f_cur=floors[idx]
            if idx<k-1:e_nx=order[idx+1];f_nx=floors[idx+1];wrap=F(0)
            else:e_nx=order[0];f_nx=floors[0];wrap=F(1)
            A=F(e_nx-e_cur);Cc=F(f_cur-f_nx)+wrap
            if A==0:
                if Cc>ONE7:good.append((x0,x1))
                continue
            xb=(ONE7-Cc)/A
            if A>0:lo=max(x0,xb);hi=x1
            else:lo=x0;hi=min(x1,xb)
            if lo<hi:good.append((lo,hi))
    return meas(merge(good))
# thr_k
thr={}
for k in range(8,14):
    psz=13-k
    thr[k]=1-min(meas(safe_set(list(P))) for P in itertools.combinations(range(1,14),psz))
print("Candidate spread-bound: mu_1/7(E) >= mu_1/7(consecutive_k). Exhaustive small k, broad large k.")
print("Also: does the consecutive value itself clear thr_k? (then union bound is unconditional on")
print("the FULL search if the spread bound is granted, with the listed slack.)")
for k in range(8,14):
    muc=good17(list(range(k)))
    # exhaustive for k<=9 over spread<=k+3; broad random for k>=10
    minmu=muc; below=None
    if k<=9:
        body_space=list(range(1,k+4))
        for body in itertools.combinations(body_space,k-1):
            E=[0]+list(body)
            m=good17(E)
            if m<minmu: minmu=m; below=E
    else:
        cands=set()
        base=list(range(k+3))
        for drops in itertools.combinations(range(1,k+2),(k+3)-k):
            E=tuple(x for x in base if x not in drops)
            if len(E)==k and E[0]==0: cands.add(E)
        for _ in range(300):
            sp=random.choice([k+1,k+2,2*k])
            E=tuple([0]+sorted(random.sample(range(1,sp+1),min(k-1,sp))))
            if len(set(E))==k: cands.add(E)
        for E in cands:
            m=good17(list(E))
            if m<minmu: minmu=m; below=E
    clears = minmu>=thr[k]
    consec_clears = muc>=thr[k]
    print(f"  k={k:2d}: min mu_1/7 found={float(minmu):.4f}({minmu})  consec={float(muc):.4f}  thr_k={float(thr[k]):.4f}")
    print(f"         consec>=thr_k: {consec_clears} ; min_found>=thr_k: {clears} ; spread-bound (min==consec): {below is None}")
print()
print("INTERPRETATION: if 'spread-bound (min==consec)' holds at every k, then granting the 1/7")
print("spread bound mu_1/7(E)>=mu_1/7(consecutive_k) makes the k>=8 union floor 1891/5880 RIGOROUS.")
print("DONE.")
