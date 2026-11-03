#!/usr/bin/env python3
"""Independent closed-interval and event-period checks. Does not import check.py."""
from fractions import Fraction as F
from itertools import combinations
from math import gcd
from functools import reduce
import json


def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def intersect(left,right):
    out=[]
    for a,b in left:
        for c,d in right:
            lo=max(a,c);hi=min(b,d)
            if lo<=hi:out.append((lo,hi))
    return out

def feasible(vs,delta=F(1,8)):
    out=[(F(0),F(1))]
    for v in vs:
        out=intersect(out,[((k+delta)/v,(k+1-delta)/v) for k in range(v)])
    return sorted(set(out))

def events(v,delta):
    return {((F(k)+p)/v)%1 for k in range(v) for p in (delta,F(1,2)-delta,F(1,2)+delta,1-delta)}

def detector(vs,t,delta):
    return tuple(int(((v*t-g)%1)<=F(1,2)) for v in vs for g in (delta,F(1,2)-delta))

def trace_period(vs,p,delta):
    points={F(0),F(1)}
    for v in vs:
        es=events(v,delta);points|=es;points|={(t-p)%1 for t in es}
    points=sorted(points)
    samples=points[:-1]+[(a+b)/2 for a,b in zip(points,points[1:])]
    return all(detector(vs,t,delta)==detector(vs,t+p,delta) for t in samples)

def clock_period(vs,delta):
    es=[events(v,delta) for v in vs]
    first=next(iter(es[0]))
    candidates=sorted({(e-first)%1 for e in es[0]}-{F(0)})+[F(1)]
    return next(p for p in candidates if all({(e+p)%1 for e in E}==E for E in es))

def main():
    examples=[(1,2,3,4,5,6,7),(7,13),(1,7,13,91),(1,2,4,8),(1,3,4,5,7,13,18)]
    certificates=[]
    for vs in examples:
        ints=feasible(vs)
        isolated=[str(a) for a,b in ints if a==b]
        certificates.append(dict(speeds=vs,closed_components=[[str(a),str(b)] for a,b in ints],isolated=isolated))
    need(certificates[0]['isolated']==['1/8','3/8','5/8','7/8'],'tight control')
    checks=0
    for n in range(1,5):
        for vs in combinations(range(1,13),n):
            intervals=feasible(vs)
            isolated=[a for a,b in intervals if a==b]
            for t in isolated:
                ent=[v for v in vs if v*t%1==F(1,8)]
                ext=[v for v in vs if v*t%1==F(7,8)]
                need(ent and ext,'isolated point has opposing active constraints')
                need(all((a+b)//gcd(a,b)%8==0 for a in ent for b in ext),'isolated point resonance')
            for a,b in intervals:
                t=(a+b)/2
                need(all(F(1,8)<=v*t%1<=F(7,8) for v in vs),'interval midpoint')
                if a<b:need(all(F(1,8)<v*t%1<F(7,8) for v in vs),'strict interval midpoint')
            checks+=1
    periods=[]
    for delta in (F(1,10),F(1,8),F(1,6),F(1,5)):
        for vs in ((1,),(2,),(3,7),(6,10),(2,4,6),(4,6,10),(7,13)):
            g=reduce(gcd,vs);m=4 if delta==F(1,8) else 2
            p=clock_period(vs,delta)
            need(p==F(1,m*g),'minimal labeled event period')
            for k in range(1,m+1):need(trace_period(vs,k*p,delta)==(k==m),'exact physical return multiplier')
            periods.append(dict(delta=str(delta),speeds=vs,clock_period=str(p),trace_period=str(m*p),laps=m))
    print(json.dumps(dict(status='PASS',interval_tuple_checks=checks,clock_cases=len(periods),examples=certificates,periods=periods),indent=2,sort_keys=True))

if __name__=='__main__':main()
