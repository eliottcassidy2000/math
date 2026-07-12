#!/usr/bin/env python3
"""cont.48: quantify kps's looseness dichotomy (HYP-6120). For DC families, is M bounded
AWAY from 1/14 when diameter is large? Test M vs diameter and vs coverage-deficiency.
If large-diameter DC families are all LOOSE (M >> 1/14), the dichotomy's hard half closes."""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random
def M_pairsum(v, Qextra=60):
    # M via best clearance over pair-sum moduli (THM-668) + small q
    best=F(0)
    qs=set()
    for i in range(len(v)):
        for j in range(len(v)):
            if i!=j:
                qs.add(v[i]+v[j]); 
                if v[i]>v[j]: qs.add(v[i]-v[j])
    for q in list(qs)+list(range(2,Qextra)):
        if q<2: continue
        for p in range(1,q):
            if gcd(p,q)!=1: continue
            m=None
            for e in v:
                r=(e*p)%q; d=min(r,q-r); dd=F(d,q)
                if m is None or dd<m: m=dd
            if m>best: best=m
    return best
def deficiency(v):
    E=list(v)
    pts=set([F(0),F(1)])
    for e in E:
        for m in range(e):
            for s in range(8): pts.add(F(m,e)+F(s,7*e))
    pts=sorted(x for x in pts if 0<=x<=1); d=F(0)
    for a,b in zip(pts,pts[1:]):
        mid=(a+b)/2; hit=set()
        for e in E:
            fr=e*mid-(e*mid).__floor__(); hit.add(int(fr*7))
        if len(hit)<7: d+=b-a
    return float(d)
def make_dc(seed, scale):
    random.seed(seed)
    base=[2,3,4,6,12]; extra=sorted(random.sample(range(13,13+scale),8))
    v=sorted(set(base+extra))
    while len(v)<13: v.append(max(v)+random.randint(1,3)); v=sorted(set(v))
    v=v[:13]; g=reduce(gcd,v); return [x//g for x in v]
print(f"1/14 = {float(F(1,14)):.5f}; looseness dichotomy: large diameter => M bounded away from 1/14?")
print(f"{'family':30s} {'diam':>6s} {'defic':>7s} {'M(pairsum)':>11s}  M/(1/14)")
def show(nm,v):
    diam=max(v)-min(v); d=deficiency(v); M=M_pairsum(v)
    print(f"  {nm:28s} {diam:6d} {d:7.3f} {float(M):11.5f}  {float(M)/float(F(1,14)):.2f}x")
    return diam,d,float(M)
# kps's blocker example (M should be ~0.234)
show("kps blocker",[200,496,540,656,851,921,935,1122,1482,1680,1835,1849,1856])
# AP (the wall, M=1/14)
show("AP {1..13}",list(range(1,14)))
# DC families increasing diameter
rows=[]
for scale in [10,25,50,100,200]:
    Ms=[]
    for s in range(4):
        diam,d,M=show(f"DC scale{scale} s{s}",make_dc(s,scale))
        Ms.append((diam,M))
    rows.append((scale,min(m for _,m in Ms)))
print(f"\nmin M over DC families by scale: {[(s,round(m,4)) for s,m in rows]}")
print(f"=> if min M > 1/14 = {float(F(1,14)):.4f} and grows with diameter, large-diameter DC is LOOSE (dichotomy closes)")
