#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Clean isolated run: min mu_1/7 over bounded+large spread E vs thr_k. Single process, single file.
import itertools, random
from fractions import Fraction as F
random.seed(2026)
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
thr={}
for k in range(8,14):
    psz=13-k
    thr[k]=1-min(meas(safe_set(list(P))) for P in itertools.combinations(range(1,14),psz))
print("min mu_1/7 over E (bounded+large spread) vs thr_k=1-min_P meas(G_P):")
allok=True
for k in range(8,14):
    minmu=good17(list(range(k))); argE=tuple(range(k))
    cands=set()
    # spread-tail
    for M in [k,k+1,k+5,2*k,5*k,20*k]:
        if M>=k-1: cands.add(tuple(list(range(k-1))+[M]))
    # perforations of {0..k+3}
    base=list(range(k+4))
    for drops in itertools.combinations(range(1,k+3),(k+4)-k):
        E=tuple(x for x in base if x not in drops)
        if len(E)==k and E[0]==0: cands.add(E)
    # random moderate spread
    for _ in range(600):
        sp=random.choice([k+1,k+2,2*k,4*k,10*k])
        body=sorted(random.sample(range(1,sp+1),min(k-1,sp)))
        E=tuple([0]+body)
        if len(set(E))==k: cands.add(E)
    for E in cands:
        m=good17(list(E))
        if m<minmu: minmu=m; argE=E
    slack=minmu-thr[k]
    if slack<=0: allok=False
    print(f"  k={k:2d}: min mu_1/7={float(minmu):.4f} ({minmu}) E={argE} | thr_k={float(thr[k]):.4f} | SLACK={float(slack):.4f} {'OK' if slack>0 else 'FAIL'}")
print(f"all-OK (union bound holds at all spreads tested): {allok}")
