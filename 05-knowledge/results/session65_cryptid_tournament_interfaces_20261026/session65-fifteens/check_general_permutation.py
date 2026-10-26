#!/usr/bin/env python3
"""Root's general signed pair-orbit formula; independent bounded replay."""
from itertools import combinations
from collections import Counter
from math import comb,prod
from pathlib import Path
import json,hashlib

def need(c,m):
    if not c:raise RuntimeError(m)
def parts(n,m=1):
    if n==0:yield []
    for k in range(m,n+1):
        for rest in parts(n-k,k):yield [k]+rest

def check(part):
    n=sum(part);sigma=list(range(n));start=0
    for size in part:
        for j in range(size):sigma[start+j]=start+(j+1)%size
        start+=size
    pairs=list(combinations(range(n),2));idx={p:i for i,p in enumerate(pairs)}
    action=[]
    for a,b in pairs:
        x,y=sigma[a],sigma[b];action.append((idx[min(x,y),max(x,y)],int(x>y)))
    seen=set();orbits=[];pred=[1]
    for i in range(len(pairs)):
        if i in seen:continue
        cur=i;sign=0;members=[]
        while cur not in seen:
            seen.add(cur);members.append(cur);cur,twist=action[cur];sign^=twist
        need(cur==i,'orbit closure')
        L=len(members);epsilon=(-1)**sign
        factor=[comb(L,j)*(1+epsilon*(-1)**j) for j in range(L+1)]
        nxt=[0]*(len(pred)+L)
        for x,c in enumerate(pred):
            for y,d in enumerate(factor):nxt[x+y]+=c*d
        pred=nxt;orbits.append({'length':L,'epsilon':epsilon})
    got=Counter()
    for mask in range(1<<len(pairs)):
        defect=sum(((mask>>i)&1)!=(((mask>>j)&1)^twist) for i,(j,twist) in enumerate(action))
        got[defect]+=1
    want={j:x for j,x in enumerate(pred) if x}
    need(dict(got)==want,('polynomial',part,got,want))
    even=sum(x%2==0 for x in part)
    neg=[o['length'] for o in orbits if o['epsilon']==-1]
    need(len(neg)==even and sorted(neg)==sorted(x//2 for x in part if x%2==0),'negative orbit classification')
    need(min(got)==even and got[even]==2**len(orbits)*prod(neg),'minimum and count')
    need(all(d%2==even%2 for d in got),'parity')
    return {'vertex_cycle_type':part,'unordered_pair_orbits':orbits,'minimum':even,'minimizers':got[even],'defect_distribution':dict(sorted(got.items()))}
rows=[check(part) for n in range(7) for part in parts(n)]
print(json.dumps({'status':'FINITE-EXACT all vertex-permutation conjugacy types and all tournament orientations for n<=6','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'cases':len(rows),'tournaments_tested_with_multiplicity':sum(sum(r['defect_distribution'].values()) for r in rows),'rows':rows},sort_keys=True,indent=2))
