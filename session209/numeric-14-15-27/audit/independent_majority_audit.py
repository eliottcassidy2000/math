#!/usr/bin/env python3
"""Independent gap-composition and exact sign audit for the fixed-H27 majority."""
import json
from fractions import Fraction
from math import comb

def need(p,msg):
    if not p:raise RuntimeError(msg)

rows=[];gap_vectors=0
for n in range(6,61):
    reps=set();pals=0
    for a in range(max(0,n-9)+1):
        for b in range(max(0,n-9-a)+1):
            for c in range(max(0,n-9-a-b)+1):
                d=n-9-a-b-c
                if d<0:continue
                g=(a,b,c,d);gap_vectors+=1
                if g==g[::-1]:pals+=1
                else:reps.add(min(g,g[::-1]))
    # Count the ordered strong words without collapsing to the cubic first.
    N=2*(n-5)+3*max(0,n-6)
    if n>=8:N+=4*comb(n-6,2)
    if n>=9:N+=comb(n-6,3)
    S=(n-5)//2 if n%2 else 0
    merged=(N+S)//2
    sign=6*(4*len(reps)-N-S)
    poly=n**3-33*n*n+(263 if n%2 else 272)*n-(615 if n%2 else 672)
    need(sign==poly,('polynomial',n,sign,poly))
    need((2*len(reps)>merged)==(n>=23),('first majority',n))
    rows.append({'n':n,'leaves':len(reps),'merged':merged,
                 'fraction':str(Fraction(len(reps),merged)),'sign':sign})
need(all(x['sign']<0 for x in rows if x['n']<=22),'all initial signs negative')
need(6*(23-10)**2-198>0,'positive lower bound for both increasing parity tails')
print(json.dumps({'status':'PASS','gap_vectors_examined':gap_vectors,
                  'literal_max_order':60,'first_majority':23,
                  'selected':[r for r in rows if r['n'] in (14,21,22,23,24)]},indent=2))
