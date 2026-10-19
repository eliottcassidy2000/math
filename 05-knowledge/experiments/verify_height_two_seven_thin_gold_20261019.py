#!/usr/bin/env python3
"""Exact independent small-family controls; no external dependencies."""
from itertools import combinations, permutations
from collections import Counter
import json

def require(ok, msg):
    if not ok: raise RuntimeError(msg)

def near(x,y):
    return (x & ~y).bit_count()<=1 and (y & ~x).bit_count()<=1

def orbit(rows,b):
    return {tuple(sorted(sum(1<<p[i] for i in range(b) if r>>i&1) for r in rows)) for p in permutations(range(b))}

def main():
    records=[]; survivors={}
    for a in range(2,7):
        for b in range(1,a+1):
            masks=[m for m in range(1<<b) if m.bit_count()<=8-a]
            families=[]
            def visit(rs,start):
                if len(rs)==a:
                    cols=[sum((bool(r>>j&1)<<i) for i,r in enumerate(rs)) for j in range(b)]
                    require(any(near(c,d) for c,d in combinations(cols,2)), 'no near column pair')
                    families.append(tuple(rs));return
                for k in range(start,len(masks)):
                    r=masks[k]
                    if all(not near(r,s) for s in rs):visit(rs+[r],k+1)
            visit([],0)
            records.append({'a':a,'b':b,'separated_families':len(families)})
            if families:survivors[a,b]=set(families)
    require(set(survivors)=={(2,2),(4,4),(5,5)},'unexpected dimensions')
    require(survivors[2,2]=={(0,3)},'two rows')
    four=orbit((0,3,12,15),4)
    five_empty=orbit((0,3,12,21,26),5)
    five_single=orbit((16,3,12,21,26),5)
    require(survivors[4,4]==four and len(four)==3,'four-row classification')
    require(survivors[5,5]==five_empty|five_single,'five-row classification')
    require(len(five_empty)==len(five_single)==30 and not five_empty&five_single,'five-row orbits')
    print(json.dumps({'status':'ALL_SMALL_FAMILY_CONTROLS_PASSED','records':records,'orbit_sizes':{'four':[3],'five':[30,30]},'surviving_families':64},indent=2))
if __name__=='__main__':main()
