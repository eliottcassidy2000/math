"""Literal depth-first enumeration of all directed Hamilton paths in G_25."""
from math import isqrt
from collections import Counter
import json
n=25;a=[[] for _ in range(n)]
for u in range(1,n+1):
 for v in range(u+1,n+1):
  if isqrt(u+v)**2==u+v:a[u-1].append(v-1);a[v-1].append(u-1)
counts=Counter(); endpoints=Counter();full=(1<<n)-1
def walk(v,visited,start,oo,ee):
 if visited==full:
  counts[oo,ee]+=1;endpoints[start+1,v+1]+=1;return
 for w in a[v]:
  if not visited>>w&1:walk(w,visited|(1<<w),start,oo+int((v+1)%2==1 and (w+1)%2==1),ee+int((v+1)%2==0 and (w+1)%2==0))
for v in range(n):walk(v,1<<v,v,0,0)
def need(b,msg):
 if not b:raise RuntimeError(msg)
need(sum(counts.values())==20,'paths');need([i+1 for i,row in enumerate(a) if len(row)==1]==[18],'degree one');need(all(18 in e for e in endpoints),'endpoint 18')
print(json.dumps({'status':'PASS','paths':sum(counts.values()),'odd_odd_even_even_histogram':[[u,v,c] for (u,v),c in sorted(counts.items())],'endpoint_pairs':[[u,v,c] for (u,v),c in sorted(endpoints.items())]},indent=2))
