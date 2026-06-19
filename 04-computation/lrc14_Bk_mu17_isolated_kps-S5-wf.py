#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Isolated, authoritative mu_{1/7} for the consecutive E (the searched minimizer) at each k>=8,
# plus the specific E=(0,2,3,4,5,7,8,9,10) that a corrupted interleaved file showed as >1.
import itertools
from fractions import Fraction as F
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
print("consecutive mu_1/7 and the suspect E:")
for k in range(8,14):
    print(f"  k={k}: consecutive mu_1/7 = {good17(list(range(k)))}")
sus=[0,2,3,4,5,7,8,9,10]
print(f"  suspect E={sus}: mu_1/7 = {good17(sus)} (must be <=1)")
# Also the EXACT per-pair union floor over k>=8 using the (now verified <=1) consecutive minima:
print()
print("EXACT k>=8 union floor = min_{P,k} [ meas(G_P) + mu_1/7(consec_k) - 1 ]:")
gm=(F(10),None,None)
for k in range(8,14):
    psz=13-k; muc=good17(list(range(k)))
    for P in itertools.combinations(range(1,14),psz):
        lb=meas(safe_set(list(P)))+muc-1
        if lb<gm[0]: gm=(lb,P,k)
print(f"  = {gm[0]} = {float(gm[0]):.6f}  at k={gm[2]} P={gm[1]}")
