#!/usr/bin/env python3
"""Independent finite labelled-port control; no inherited helper imports."""
from itertools import permutations
import json


def need(ok,message):
    if not ok:raise RuntimeError(message)


def local(n,edge):
    rows=[((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]
    a,b=edge;rows[a]^=1<<b;rows[b]^=1<<a
    return tuple(rows)


def canonical(rows,port=None):
    n=len(rows);best=None
    for p in permutations(range(n)):
        if port is not None and p[0]!=port:continue
        key=tuple(sum((rows[a]>>b&1)<<j for j,b in enumerate(p)) for a in p)
        if best is None or key<best:best=key
    return best


def main():
    full=[]
    for name,edge,port,degree in [('12',(0,2),4,0),('13',(0,4),2,2),('23',(2,4),0,4)]:
        rows=local(5,edge)
        need(rows[port].bit_count()==degree,'five-cell marked degree')
        full.append(dict(pair=name,rows=rows,port=port,degree=degree,
                         rooted_key=canonical(rows,port)))
    need(len({tuple(r['rooted_key']) for r in full})==3,'three full-cell records remain distinct')
    mixed=[]
    for name,edge,port,degree in [('X13',(0,3),1,2),('X23',(1,3),0,3),
                                 ('Y12',(0,2),3,0),('Y13',(0,3),2,1)]:
        rows=local(4,edge)
        need(rows[port].bit_count()==degree,'four-cell marked degree')
        mixed.append(dict(name=name,rows=rows,port=port,degree=degree,
                          rooted_key=canonical(rows,port),bare_key=canonical(rows)))
    need(len({tuple(r['rooted_key']) for r in mixed})==4,'all four mixed-cell owners distinct')
    need(mixed[0]['bare_key']==mixed[3]['bare_key'],'forgetting the port gives a genuine collision')
    need(mixed[0]['rooted_key']!=mixed[3]['rooted_key'],'the actual port repairs that collision')
    print(json.dumps(dict(status='PASS exact finite port controls',
                         full_two_unit_cells=full,missing_plus_unit_cells=mixed),indent=2,sort_keys=True))


if __name__=='__main__':main()
