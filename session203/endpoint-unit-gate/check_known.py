from itertools import permutations
from collections import Counter
import json

def tensor(rows,mask,q):
 n=len(rows); t=Counter()
 for p in permutations(range(n)):
  d={(0,0,0):1}
  for g in range(n+1):
   l=p[g-1] if g else None;r=p[g] if g<n else None
   e=l is None or r is None or (rows[l]>>r&1)
   f=(l is None or l==q or not(mask>>l&1)) and (r is None or r==q or(mask>>r&1))
   nd=Counter()
   for (a,i,j),v in d.items():
    if e:nd[a,i,j]+=v
    if f:nd[a+1,i+int(r==q),j+int(l==q)]+=v
   d=nd
  t.update(d)
 return t
P=[2,4,41,3,15,91,31];T=[0,5,41,131,15,219,159,23]
p=tensor(P,27,5);c=tensor(T,155,5)
d={key:c[key]-2*p[key] for key in set(p)|set(c) if c[key]!=2*p[key]}
print(json.dumps({'P':{str(k):v for k,v in sorted(p.items()) if v},'C':{str(k):v for k,v in sorted(c.items()) if v},'difference':{str(k):v for k,v in sorted(d.items())}},indent=2))
