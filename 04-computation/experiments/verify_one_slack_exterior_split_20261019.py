#!/usr/bin/env python3
"""Independent standard-library verification of a fixed-core exterior split."""
from itertools import combinations
from collections import defaultdict
import json
P=(0,28,0,0,12,4,13,93,125)
O=(2,3,0,4,5,6,1,7,8)
CERTS=(((6,1,2),(4,5,1),(5,6,2)),((0,4,2),(4,5,2),(5,6,1)))
def require(test,msg):
 if not test:raise RuntimeError(msg)
N=len(P)
require(all(not(P[i]>>i&1) and all(P[j]&~P[i]==0 for j in range(N) if P[i]>>j&1) for i in range(N)), 'transitive strict core')
I=[{j for j in range(N) if i!=j and not(P[i]>>j&1 or P[j]>>i&1)} for i in range(N)]
require(len(I[5])==6,'one-slack vertex')
words=[]
def walk(prefix,remaining):
 if not remaining:words.append(prefix);return
 used=((1<<N)-1)^remaining
 for i in range(N):
  if remaining>>i&1 and not(P[i]&~used):walk(prefix+(i,),remaining^(1<<i))
walk((),(1<<N)-1)
require(len(words)==540,'extension count')
require(O in words,'order extends core')
require(O.index(1)-O.index(0)==4,'root span')
require(max(abs(O.index(i)-O.index(j)) for i in range(N) for j in I[i])==4,'maximum span')
require(I[0]&I[1]=={5},'common root vertex')
summary=[]
for side,cert in enumerate(CERTS):
 for a,b,w in cert:require(b in I[a] and O.index(a)<O.index(b) and w>0,'certificate row')
 counts={}
 for order in words:
  r={a:i for i,a in enumerate(order)}
  first=min(r[0],r[1]);last=max(r[0],r[1])
  for left in range(first+1):
   for right in range(last+1,N+1):
    key=order[:left],order[right:]
    if 5 in key[side]:continue
    if key not in counts:counts[key]=[0]+[0]*len(cert)
    counts[key][0]+=1
    for k,(a,b,w) in enumerate(cert):counts[key][k+1]+=r[a]<r[b]
 margins=[sum(w*(3*c[k+1]-2*c[0]) for k,(a,b,w) in enumerate(cert)) for c in counts.values()]
 require(len(counts)==112,'generator count')
 require(max(margins)==-2,'strict generator upper bound')
 # Verify the stronger probability envelope, using only cross-multiplication.
 require(all(3*sum(w*c[k+1] for k,(a,b,w) in enumerate(cert))<=10*c[0] for c in counts.values()),'weighted bound')
 summary.append({'case':'all A below5' if side==0 else 'all B above5','certificate':cert,'generators':len(counts),'minimum_integer_margin':min(margins),'maximum_integer_margin':max(margins)})
print(json.dumps({'status':'PASS','core':P,'order':O,'extensions':len(words),'common_vertex':5,'core_incomparability_degree':6,'cases':summary},indent=2))
