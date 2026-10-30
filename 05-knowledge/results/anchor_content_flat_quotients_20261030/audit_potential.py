#!/usr/bin/env python3
"""Read C++ potentials; independently check all needed inequalities in Python."""
from math import gcd
from pathlib import Path
import hashlib
import json
import sys

prefix=sys.argv[1] if len(sys.argv)>1 else '/tmp/five-potential'
rows=[]
for E in range(1,5):
    q=5**E
    states=[(m,b) for m in range(q) if gcd(m,q)==1 for b in range(q)]
    p=Path(f'{prefix}-{E}.txt');data=p.read_bytes();values=list(map(int,data.split()))
    if len(values)!=len(states):raise RuntimeError('potential state count')
    h=dict(zip(states,values));num=q-1;den=8*5**(E-1);g=gcd(num,den);num//=g;den//=g
    if not (min(values)==0 and h[1,0]==max(values)):raise RuntimeError('identity maximum')
    cap=2*E+4
    if num*(cap+1)<max(values)+den*E:raise RuntimeError('unbounded-letter tail')
    rewards={s:sum((s[0]-1)%5**e==0 and s[1]%5**e==0 for e in range(1,E+1)) for s in states}
    checks=0
    for a in range(1,cap+1):
        inverse=pow(2**a,-1,q)
        for s in states:
            t=(3*s[0]*inverse%q,(3*s[1]+1)*inverse%q)
            if h[t]<h[s]+den*rewards[t]-num*a:raise RuntimeError((E,a,s,t))
            checks+=1
    rows.append(dict(E=E,q=q,states=len(states),checks=checks,
                     maximum=h[1,0],sha256=hashlib.sha256(data).hexdigest()))
print(json.dumps(rows,indent=2,sort_keys=True))
