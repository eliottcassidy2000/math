#!/usr/bin/env python3
"""Exact small4-cycle obstruction to a single internal segment reversal."""
from itertools import combinations
from math import isqrt
import json
records=[]
for n in (42,45,46,48):
    neighbors={u:{v for v in range(1,n+1) if u!=v and isqrt(u+v)**2==u+v} for u in range(1,n+1)}
    cycles=set()
    for u,v in combinations(range(1,n+1),2):
        for a,b in combinations(sorted(neighbors[u]&neighbors[v]),2):
            cycles.add(tuple(sorted(tuple(sorted(e)) for e in ((u,a),(a,v),(v,b),(b,u)))))
    records.append({'N':n,'C4_count':len(cycles),'cycles':sorted(cycles)})
if [r['C4_count'] for r in records]!=[0,0,1,2]:raise RuntimeError('small graph boundary')
print(json.dumps({'status':'finite-exact','records':records,
                  'scope':'A nontrivial single-segment reversal needs a4-cycle on its boundary. This does not prove path uniqueness or exclude other moves.'},indent=2,sort_keys=True))
