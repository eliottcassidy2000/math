"""Produce a finite NAE refutation tree. No external solver needed."""
from pathlib import Path
import json, hashlib
from collections import defaultdict
from probe_finite_triangles import field
OUT=Path(__file__).resolve().parents[1]/'results'
p=17; d,zeta,pts,U,edges=field(p); n=p*p
incid=defaultdict(list)
for i,e in enumerate(edges):
 for v in e:incid[v].append(i)
nodes=0

def dfs(a):
 global nodes
 nodes+=1; forces=[]
 while True:
  changed=False
  for i,e in enumerate(edges):
   assigned=[v for v in e if a[v]>=0]
   if len(assigned)==3:
    if len({a[v] for v in e})==1:return {'forces':forces,'conflict':i}
   elif len(assigned)==2 and a[assigned[0]]==a[assigned[1]]:
    v=next(v for v in e if a[v]<0);c=1-a[assigned[0]]
    a[v]=c; forces.append([v,c,i]);changed=True
  if not changed:break
 rem=[v for v in range(n) if a[v]<0]
 if not rem:raise AssertionError('Found coloring')
 def score(v):return sum(10 if sum(a[w]>=0 for w in edges[i])==1 else 1 if all(a[w]<0 for w in edges[i]) else 0 for i in incid[v])
 v=max(rem,key=score);children=[]
 for c in [0,1]:
  b=a[:];b[v]=c;children.append(dfs(b))
 return {'forces':forces,'branch':v,'children':children}

a=[-1]*n;a[0]=0
cert={'p':p,'d':d,'zeta':zeta,'vertices':pts,'edges':edges,'root_fixed':[0,0],'proof':dfs(a)}
raw=json.dumps(cert,separators=(',',':')).encode();(OUT/'q17_nae_certificate.json').write_bytes(raw)
print(json.dumps({'p':p,'vertices':n,'edges':len(edges),'proof_nodes':nodes,'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()},sort_keys=True))
