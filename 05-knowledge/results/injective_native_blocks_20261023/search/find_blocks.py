"""Discover four-source balanced charges on an actual shared target pair."""
import json
from itertools import combinations
from collections import defaultdict
import explore as e
HERE=e.HERE
data=json.loads((HERE/'cache_all.json').read_text());by={r['source']:r for r in data}
menus=defaultdict(dict)
for s,r in by.items():
    for u,v,w in r['menu']:
        if s not in (u,v):menus[u,v][s]=w

def code(s):return (e.X[s]%128,)+e.SIG[s]
mods=(128,4,2,2,2,2)
def add(a,b):return tuple((x+y)%m for x,y,m in zip(a,b,mods))
found=[]
for (u,v),ss in menus.items():
    if (e.X[u]-e.X[v])%128:continue
    buckets=defaultdict(list)
    for a,b in combinations(sorted(ss),2):
        key=add(code(a),code(b))
        for c,d in buckets[key]:
            if len({a,b,c,d})<4:continue
            # No cancellation by a legal two-slot proper subblock.
            if any(code(i)==code(j) for i in (a,b) for j in (c,d)):continue
            slots=[]
            for s,sgn in [(a,1),(b,1),(c,-1),(d,-1)]:
                w=ss[s];cost,pins,i,j=w
                if sgn<0:i,j=j,i
                slots.append(dict(source=s,sign=sgn,cost=cost,pins=pins,
                                  before=by[s]['rows'][i],after=by[s]['rows'][j]))
            found.append(dict(targets=[u,v],target_x=[e.X[u],e.X[v]],
                              sources=[a,b,c,d],source_x=[e.X[s] for s in [a,b,c,d]],
                              source_codes=[code(s) for s in [a,b,c,d]],slots=slots,
                              cost=sum(z['cost'] for z in slots),max_cost=max(z['cost'] for z in slots)))
        buckets[key].append((a,b))
found.sort(key=lambda r:(r['cost'],r['max_cost'],r['targets'],r['sources']))
(HERE/'candidate_blocks.json').write_text(json.dumps(found,indent=2))
print('candidate blocks',len(found))
print(json.dumps(found[:3],indent=2))
