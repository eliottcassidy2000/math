#!/usr/bin/env python3
"""Exact finite controls for the proved, all-integer short-block classification."""
import json
from collections import Counter

def F(n):
    c=n & -n
    return n-n//c+2*c-1

def check(ok, message):
    if not ok:
        raise RuntimeError(message)

def branch(n):
    c=n & -n; p=n//c; q=p+1
    if q>2*c: return 'immediate',1
    if q==2*c: return 'fixed',0
    if q & (q-1): return 'nonpower_q',2
    if p==1: return 'power_start',3
    d=q
    if c==d: return 'first_resonance',4
    if c==d*d: return 'second_resonance',6 if d==4 else 5
    if c<d*d: return 'middle',3
    return 'high',4 if c==2*d*d else 3

exceptions={2,4,8,12,16}
counts=Counter(); cases=0

def audit(n):
    global cases
    kind,bound=branch(n); counts[kind]+=1; cases+=1
    if n in exceptions: return
    if bound==0:
        check(F(n)==n,('fixed',n)); return
    x=n
    for j in range(bound):
        x=F(x)
        if x<n: return
    raise RuntimeError(('descent failed',n,kind,bound,x))

for n in range(2,1000001,2): audit(n)
for a in range(1,181):
    c=1<<a
    for t in range(1,181):
        for k in (1,3,5,7,11):
            audit(c*((1<<t)*k-1))
    for p in (1,3,5,2*c-3,2*c-1,2*c+1):
        if p>0: audit(c*p)
for orbit in ([8,22,14,10,8],[16,46,26,16],[48,76,64,190,98,52,46],[2,4,10,8],[12,16]):
    check(all(F(a)==b for a,b in zip(orbit,orbit[1:])),('orbit',orbit))
check(F(120)==120,'formal fixed point')
print(json.dumps({'cases':cases,'branches':dict(sorted(counts.items())),'status':'PASS','scope':'Exact finite controls; universal proof is in the accompanying manuscript.'},sort_keys=True,indent=2))
