#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ADVERSARIAL part 4:
 (a) check the propagation-budget arithmetic: sum_{r>=1} |y_r| C(6,r), and
     the per-A budget eps needed so that |L_y - L_inf| stays below (cap - L_inf).
 (b) hunt WIDE primitive sets (spread > 16) trying to push L_y above the claimed
     0.179 ceiling (k=8). Includes adversarial near-AP / near-dissociated mixes.
"""
import sys, itertools, random
from fractions import Fraction as F
from functools import reduce
from math import gcd, comb

def N_at(E,x):
    hit=set()
    for e in E:
        v=e*x; v=v-(v.numerator//v.denominator)
        hit.add((v.numerator*7)//v.denominator)
    return sum(1 for j in range(1,7) if j not in hit)

def dist_p(E):
    E=sorted(set(E))
    bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(F(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1)
    p=[F(0)]*7
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        p[N_at(E,(lo+hi)/2)]+=(hi-lo)
    return p

def g_poly(k):
    g=[]
    for t in range(7):
        if k==8: g.append(F((t-1)*(t-2)*(t-4)*(t-5),40))
        elif k in (9,10): g.append(F(-(t-2)*(t-3)*(t-6),36))
        else: g.append(F((t-3)*(t-4),12))
    return g

def y_coeffs(k):
    g=g_poly(k); y=[]
    for r in range(7):
        y.append(sum((-1)**(r-j)*comb(r,j)*g[j] for j in range(r+1)))
    return y

def L_y(E,k):
    p=dist_p(E); g=g_poly(k)
    return sum(p[t]*g[t] for t in range(7))

caps={8:0.3815,9:0.49426,10:0.6044}

print("=== (a) propagation budget constants sum_{r>=1}|y_r| C(6,r) ===")
for k in [8,9,10]:
    y=y_coeffs(k)
    const=sum(abs(y[r])*comb(6,r) for r in range(1,7))
    Linf=sum(y[r]*comb(6,r)*F(7-r,7)**(k-1) for r in range(7))
    cap=caps[k]
    slack=cap-float(Linf)
    budget=slack/float(const)
    print(f"k={k}: sum|y_r|C(6,r) (r>=1) = {const} = {float(const):.4f}")
    print(f"      L_inf={float(Linf):.6f} cap={cap} slack={slack:.6f}  per-A eps budget={budget:.6f}")
print()
print("Claim states constants 48, 14.33, 14.33 and budgets 0.00692,0.0239,0.0269.")
print()

print("=== (b) WIDE primitive-set hunt (k=8, spread>16), maximize L_y ===")
k=8
best=F(0); bestE=None
random.seed(12345)
trials=0
# structured wide families
fams=[]
# arithmetic-progression-with-large-step but primitive via a +1 jitter
for D in range(2,9):
    for jit in range(D):
        E=[0]+[i*D + (1 if i==1 else 0)*jit for i in range(1,8)]
        E=sorted(set(E))
        if len(E)==8 and reduce(gcd,E)==1 and max(E)>16:
            fams.append(E)
# near-consec blocks pushed apart
for shift in range(17,40,3):
    E=[0,1,2,3]+[shift,shift+1,shift+2,shift+3]
    if reduce(gcd,E)==1 and max(E)>16: fams.append(E)
# two clusters
for gap in range(17,60,7):
    E=[0,1,2,3,4,5,6,gap]
    if reduce(gcd,E)==1: fams.append(E)
# random primitive wide sets, max in [17,60]
for _ in range(4000):
    top=random.randint(17,60)
    rest=random.sample(range(1,top),6)
    E=sorted(set([0]+rest+[top]))
    if len(E)!=8: continue
    if reduce(gcd,E)!=1: continue
    fams.append(E)
seen=set()
for E in fams:
    key=tuple(E)
    if key in seen: continue
    seen.add(key)
    trials+=1
    L=L_y(E,k)
    if L>best:
        best=L; bestE=E
print(f"wide primitive sets tested: {trials}")
print(f"max L_y over wide sets = {best} = {float(best):.6f}  at E={bestE}")
print(f"claimed wide ceiling ~0.179; cap={caps[8]}; margin to cap = {caps[8]-float(best):.4f}")
