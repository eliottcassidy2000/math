from collections import Counter
from pathlib import Path
import csv
import itertools
import sys
from selector import profile,dp,solve

def direct(rows):
    n=len(rows);loads=Counter(); h=0
    for p in itertools.permutations(range(n)):
        hops=list(zip(p,p[1:]))
        bad=[(u,v) for u,v in hops if not rows[u]>>v&1]
        if not bad:h+=1;loads.update(hops)
        if len(bad)==1:loads.update(bad)
    return h,loads

hist, examples={},{}
first_inert=None
for order in range(1,7):
    for target in range(1,46,2):
        W=solve(target,True,order)
        if not W or len(W[0])!=order:continue
        part=hist.setdefault(order,Counter())
        for w in W:
            P=profile(w)
            h,loads=direct(w)
            assert h==target
            for u,v,a,b in P['arcs_uv_load_reverse_load']:
                assert (loads[u,v],loads[v,u])==(a,b)
                changed=list(w);changed[u]^=1<<v;changed[v]^=1<<u
                assert sum(dp(changed)[-1])==h+b-a
            part[P['unused'],P['inert']]+=1
print('birth-class unused and inert histogram through six (all classes)')
for n,A in hist.items():print(n,sorted(A.items()))
if len(sys.argv)>1:
    atlas={}
    with open(sys.argv[1]) as f:
        for r in csv.DictReader(f,delimiter='\t'):
            n,h=int(r['order']),int(r['count'])
            rows=tuple(map(int,r['outgoing_masks'].split(',')))
            P=profile(rows);assert P['H']==h
            part=atlas.setdefault(n,Counter())
            part[P['unused'],P['inert']]+=1
            examples.setdefault(n,P)
            if P['inert'] and first_inert is None:first_inert=(n,h,rows,P['zero_arcs'])
    print('input birth atlas (all classes if first_classes.cpp; one per value if inherited witnesses.tsv)')
    for n,A in atlas.items():print(n,sorted(A.items()))
    print('first representative with an inert arc',first_inert)
P=profile((10,28,25,16,1))
assert P['H']==9 and P['zero_arcs']==[[1,4,4]]
print('H9 birth example',P['rows'],'zero arc and reverse load',P['zero_arcs'])
