#!/usr/bin/env python3
# Decision-support scan: how close to cap_k can WIDE, STRUCTURED (no dissociated stranger)
# primitive k-sets get? And does consec genuinely dominate L_y at moderate spread with NO ties?
# Exact Fractions. k=8 (tightest finite check region) and k=9 (tightest cap margin).
from fractions import Fraction as F
from math import gcd
from functools import reduce
import itertools, random

def Nat(E,x):
    hit=set()
    for e in E:
        v=e*x; v=v-(v.numerator//v.denominator)
        hit.add((v.numerator*7)//v.denominator)
    return sum(1 for j in range(1,7) if j not in hit)

def dist(E):
    E=sorted(set(E)); bps={F(0),F(1)}
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(F(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1)
    p=[F(0)]*7
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi<=lo: continue
        t=Nat(E,(lo+hi)/2); p[t]+=hi-lo
    return p

def gpoly(k):
    g=[]
    for t in range(7):
        if k==8: g.append(F((t-1)*(t-2)*(t-4)*(t-5),40))
        elif k in (9,10): g.append(F(-(t-2)*(t-3)*(t-6),36))
        else: g.append(F((t-3)*(t-4),12))
    return g

def Ly(E,k):
    p=dist(E); g=gpoly(k); return sum(p[t]*g[t] for t in range(7))

caps={8:0.3815,9:0.49426,10:0.6044}
def primitive(E):
    return reduce(gcd,E)==1

# 1) consec value and cap margins
for k in (8,9,10):
    C=list(range(k)); print(f"k={k}: L_y(consec)={float(Ly(C,k)):.6f}  cap={caps[k]}  margin={caps[k]-float(Ly(C,k)):.6f}")

print("\n=== WIDE-STRUCTURED hunt (primitive, every elt in a height-2 relation, large spread) ===")
# build sets that are unions of short APs / have rich relation lattice but wide spread
random.seed(7)
for k in (8,9):
    cap=caps[k]; C=list(range(k)); LC=Ly(C,k)
    best=(None,F(0))
    tries=0
    for _ in range(4000):
        # union of two short APs -> rich relations, can be wide
        d1=random.randint(1,4); d2=random.randint(1,4)
        a=random.randint(8,40)
        s1=[i*d1 for i in range(k//2+1)]
        s2=[a+i*d2 for i in range(k-len(s1))]
        E=sorted(set([0]+s1+s2))
        if len(E)<k: continue
        E=E[:k]
        if not primitive(E): continue
        if max(E)<=16: continue  # wide only
        L=Ly(E,k); tries+=1
        if L>best[1]: best=(E,L)
    print(f"k={k}: wide-structured max L_y found = {float(best[1]):.6f} (set {best[0]}), cap={cap}, margin-to-cap={cap-float(best[1]):.6f}, vs consec {float(LC):.6f}")

print("\n=== single-outlier configs (consec_{k-1} + one far point) — the high wide-set ceiling ===")
for k in (8,9):
    cap=caps[k]; best=(None,F(0))
    for d in range(k, 60):
        E=list(range(k-1))+[d]
        if not primitive(E): continue
        L=Ly(E,k)
        if L>best[1]: best=(E,L)
    print(f"k={k}: max single-outlier L_y={float(best[1]):.6f} (set {best[0]}), cap={cap}, margin={cap-float(best[1]):.6f}")
