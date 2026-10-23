#!/usr/bin/env python3
"""Complete 409-path atlas, retaining actual endpoints and pointed partner type."""
from collections import Counter,defaultdict
from itertools import permutations,combinations
import json
import native_cube as nc

def canonical(K,first,last):
    remaining=[u for u in range(8) if u not in (first,last)]
    best=None;order=None
    for middle in permutations(remaining):
        p=(first,last)+middle
        code=0
        for i,j in combinations(range(8),2):code=(code<<1)|K[p[i]][p[j]]
        if best is None or code<best:best,order=code,p
    return best,order

def main():
    paths=nc.core();fibres=defaultdict(list)
    for p in paths:
        K=nc.partner(nc.A,p);code,order=canonical(K,p[0],p[-1])
        fibres[p[0],p[-1],code].append((p,K,order))
    hist=Counter(map(len,fibres.values()))
    for rows in fibres.values():
        p,K,order=rows[0]
        for q,L,neworder in rows:
            g=dict(zip(order,neworder))
            nc.need(g[p[0]]==p[0] and g[p[-1]]==p[-1], 'literal endpoints fixed')
            nc.need(all(K[u][v]==L[g[u]][g[v]] for u in range(8) for v in range(8)), 'canonical witness authenticated as actual partner map')
    moving=sum(len(rows) for rows in fibres.values() if len(rows)>1)
    nc.need(sum(size*count for size,count in hist.items())==409,'complete path coverage without quotient weights')
    nc.need(any({nc.P[0],nc.P[1]}<=set(p for p,K,o in rows) for rows in fibres.values()),'seed pair belongs to one pointed chart')
    return dict(status='FINITE-EXACT complete original-source atlas',actual_paths=409,
                charts=len(fibres),chart_size_histogram=dict(sorted(hist.items())),
                paths_in_nontrivial_charts=moving,checks=nc.CHECKS,
                charts_data=[dict(first=k[0],last=k[1],pointed_code=k[2],paths=[list(p) for p,K,o in rows])
                             for k,rows in sorted(fibres.items())])

if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
