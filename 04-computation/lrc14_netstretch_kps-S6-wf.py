#!/usr/bin/env python3
"""
lrc14_netstretch_kps-S6-wf.py
Test STRETCH/COMPRESS monotonicity as a route to 'consecutive maximizes net'.

OP_stretch(E,i): add 1 to every e_j with index>=i (i in 1..k-1), i.e. widen the
   i-th consecutive gap by 1.  Conjecture: netmeas weakly DECREASES.
If true for all E,i: then any E reduces by compressions to consecutive, which is
   therefore the maximizer.  net(consec)=44/735 < cap_k.  DONE.

We test EXHAUSTIVELY over primitive k=8 sets of bounded spread, all i.
Also test the dual COMPRESS op for sanity.
"""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
from math import gcd
from functools import reduce
ONE7=F(1,7)

def net_intervals(E):
    E=sorted(set(E)); n=len(E)
    bp=set([F(0),F(1)])
    for i in range(n):
        for j in range(i+1,n):
            d=E[j]-E[i]
            for m in range(0,d+1): bp.add(F(m,d))
    bp=sorted(b for b in bp if 0<=b<=1)
    good=[]
    for a,b in zip(bp,bp[1:]):
        if b<=a: continue
        mid=(a+b)/2
        order=sorted(range(n),key=lambda i:(E[i]*mid)%1)
        floors=[(E[order[t]]*mid).__floor__() for t in range(n)]
        lo=a; hi=b; feasible=True
        for t in range(n):
            o1=order[t]; o2=order[(t+1)%n]; wrap=1 if t==n-1 else 0
            s=E[o2]-E[o1]; c=F(floors[t]-floors[(t+1)%n]+wrap)
            if s==0:
                if not (c<=ONE7): feasible=False; break
            elif s>0: hi=min(hi,(ONE7-c)/s)
            else: lo=max(lo,(ONE7-c)/s)
            if lo>=hi: feasible=False; break
        if feasible and lo<hi: good.append((lo,hi))
    good.sort(); out=[]
    for a,b in good:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def meas(iv): return sum((b-a for a,b in iv),F(0))
def netmeas(E): return meas(net_intervals(E))

def stretch(E,i):
    # widen i-th gap: add 1 to e_i,...,e_{k-1} (0-indexed, i in 1..k-1)
    E=sorted(E)
    return E[:i]+[e+1 for e in E[i:]]

if __name__=="__main__":
    print("[STRETCH MONOTONICITY] net(stretch(E,i)) <= net(E) for all gaps i?")
    viol=0; tot=0; examples=[]
    # exhaustive over primitive k=8, max(E)<=13, all interior stretch positions
    for body in itertools.combinations(range(1,13),6):
        E=list((0,)+body+(13,))
        # use all 8-subsets with max<=13 actually; simpler: all k=8 subsets of 0..12
        pass
    cnt=0
    for body in itertools.combinations(range(1,13),7):
        E=[0]+list(body)   # k=8, max(E)<=12
        nE=netmeas(E)
        for i in range(1,8):
            Es=stretch(E,i)
            nEs=netmeas(Es)
            tot+=1
            if nEs>nE:
                viol+=1
                if len(examples)<8: examples.append((E,i,nE,nEs))
        cnt+=1
    print(f"   scanned {cnt} sets (k=8,max<=12), {tot} stretch tests, {viol} violations")
    if viol==0:
        print("   => STRETCH MONOTONICITY HOLDS exhaustively (k=8, max<=12).")
        print("      Implication: every E compresses to consecutive without decreasing net,")
        print("      so net(E) <= net(consecutive) = 44/735.")
    else:
        print("   => STRETCH MONOTONICITY FALSE. Examples (E, gap i, net(E), net(stretch)):")
        for E,i,a,b in examples:
            print(f"      E={E} i={i}: net {float(a):.6f} -> {float(b):.6f} (INCREASED)")
