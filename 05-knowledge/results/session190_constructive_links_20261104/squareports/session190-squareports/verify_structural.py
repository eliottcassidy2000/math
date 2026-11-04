#!/usr/bin/env python3
"""Adjacency and witness checks for the human deletion-support-loss proof."""
from math import isqrt
import json

def need(ok,message):
    if not ok:raise RuntimeError(message)

def graph(n):
    return {v:{w for w in range(1,n+1) if w!=v and w!=3 and isqrt(v+w)**2==v+w}
            for v in range(1,n+1) if v!=3}

G=graph(26)
expected={18:{7},22:{14},9:{7,16},13:{12,23},26:{10,23},
          2:{7,14,23},7:{2,9,18},23:{2,13,26}}
for v,nb in expected.items():need(G[v]==nb,('adjacency',v))
P=[18,7,9,16,20,5,11,25,24,1,8,17,19,6,10,15,21,4,12,13,23,2,14,22]
H=graph(25)
need(set(P)==set(H) and len(P)==len(H),'positive witness vertices')
need(all(b in H[a] for a,b in zip(P,P[1:])),'positive witness edges')
print(json.dumps({'status':'local adjacency checks and positive witness passed',
    'Q26_minus3_adjacency':{v:sorted(nb) for v,nb in expected.items()},
    'Q25_minus3_path':P},sort_keys=True,indent=2))
