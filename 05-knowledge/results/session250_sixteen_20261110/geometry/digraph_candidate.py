from itertools import permutations
from collections import Counter
import json
n=3; edges=[(i,j)for i in range(n)for j in range(n)if i!=j];pos={e:k for k,e in enumerate(edges)}
def transport(x,p):return sum(1<<pos[p[i],p[j]]for k,(i,j)in enumerate(edges)if x>>k&1)
def canon(x):return min(transport(x,p)for p in permutations(range(n)))
def tr(x):return sum(1<<pos[j,i]for k,(i,j)in enumerate(edges)if x>>k&1)
def h(x):return sum(all(x>>pos[i,j]&1 for i,j in zip(p,p[1:]))for p in permutations(range(n)))
reps=sorted({canon(x)for x in range(64)})
rows=[]
for x in reps:
 rows.append(dict(mask=x,arcs=x.bit_count(),H=h(x),transpose=canon(tr(x)),complement=canon(x^63),complement_transpose=canon(tr(x)^63),aut=sum(transport(x,p)==x for p in permutations(range(n)))))
summary={'arc_grades':dict(sorted(Counter(x.bit_count()for x in reps).items())), 'fixed':{name:sum(row[name]==row['mask']for row in rows)for name in ['transpose','complement','complement_transpose']},'klein_orbits':len({min(x,canon(tr(x)),canon(x^63),canon(tr(x)^63))for x in reps}),'rows':rows}
assert [summary['arc_grades'][i]for i in range(7)]==[1,1,4,4,4,1,1]
assert all(h(tr(x))==h(x)for x in range(64))
print(json.dumps(summary,indent=2))
