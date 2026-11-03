#!/usr/bin/env python3
"""Independent vertex-extension count, compared with the C++ edge-submask count."""
from itertools import combinations
from pathlib import Path
import json
n=5; edges=list(combinations(range(n),2));edgeid={p:i for i,p in enumerate(edges)};all_edges=(1<<len(edges))-1
clique=[sum(1<<edgeid[p] for p in combinations([i for i in range(n) if s>>i&1],2)) for s in range(1<<n)]
triangles=[clique[s] for s in range(32) if s.bit_count()==3]
trianglefree=[all(x&t!=t for t in triangles) for x in range(1<<len(edges))]
stars=[]
for r in range(32):
    available=31^r;b=available
    while True:
        stars.append((clique[r],clique[b],available^b))
        if not b:break
        b=(b-1)&available
count5=0;counts=[0,0,0]
for red in range(1<<len(edges)):
    if not trianglefree[red]:continue
    available=all_edges^red;blue=available
    while True:
        if trianglefree[blue]:
            count5+=1;green=available^blue;core5=(green==all_edges)
            g4=[s for s in range(32) if s.bit_count()==4 and green&clique[s]==clique[s]]
            bad5=[core5 or any(s&t==t for t in g4) for s in range(32)]
            for req_red,req_blue,gstar in stars:
                if red&req_red or blue&req_blue:continue
                counts[0]+=1;counts[1]+=not bad5[gstar];counts[2]+=not(core5 and gstar==31)
        if not blue:break
        blue=(blue-1)&available
expected=json.loads((Path(__file__).parent/'ramsey_counts.jsonl').read_text().splitlines()[-1])
actual={'n':6,'edge_slots':15,'red_blue_triangle_free':counts[0],'valid_335':counts[1],'valid_336':counts[2]}
if count5!=29197 or actual!=expected:raise RuntimeError((count5,actual,expected))
print(json.dumps({'method':'independent extension of every labelled five-vertex coloring','cores':count5,**actual}))
print('ALL CHECKS PASSED')
