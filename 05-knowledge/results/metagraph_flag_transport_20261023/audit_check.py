"""Independent direct-adjacency audit; imports none of the producer's code."""
from itertools import combinations, permutations
from pathlib import Path
import json

folder=Path(__file__).resolve().parent
r=json.loads((folder/'receipt.json').read_text())
pairs=list(combinations(range(5),2)); perms=list(permutations(range(5)))
ordinary={}; rooted={};checks=0
for f in range(16):
    arcs={(i,i+1) for i in range(4) if f>>i&1}
    for code in range(1024):
        adjacency={(i,j) if code>>k&1 else (j,i) for k,(i,j) in enumerate(pairs)}
        if not arcs<=adjacency:continue
        target=(adjacency&arcs)|{(j,i) for i,j in adjacency-arcs}
        canon=[None]*6;canon2=[None]*6
        for p in perms:
            mark=sum((tuple(sorted((p[i],p[j]))) in arcs)<<k for k,(i,j) in enumerate(pairs))
            a=mark*1024+sum(((p[i],p[j]) in adjacency)<<k for k,(i,j) in enumerate(pairs))
            b=mark*1024+sum(((p[i],p[j]) in target)<<k for k,(i,j) in enumerate(pairs))
            for ix in (0,p[0]+1):
                canon[ix]=a if canon[ix] is None else min(canon[ix],a)
                canon2[ix]=b if canon2[ix] is None else min(canon2[ix],b)
        for ix in range(6):
            mapping=ordinary if ix==0 else rooted
            if canon[ix] in mapping and mapping[canon[ix]]!=canon2[ix]:
                raise ValueError('inconsistent direct quotient')
            mapping[canon[ix]]=canon2[ix];checks+=1
if sorted(ordinary.items())!=[tuple(x) for x in r['flag_map']]:
    raise ValueError('ordinary table mismatch')
if sorted(rooted.items())!=[tuple(x) for x in r['rooted_flag_map']]:
    raise ValueError('rooted table mismatch')
print(json.dumps({'direct_quotient_checks':checks,'ordinary_flags':len(ordinary),
 'rooted_flags':len(rooted),'ordinary_fixed':sum(a==b for a,b in ordinary.items()),
 'rooted_fixed':sum(a==b for a,b in rooted.items()),'complete_tables_match':True},sort_keys=True))
