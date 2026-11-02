from explore import *
import json,time
for n in [5,6]:
 t=time.time();N=1<<(n*(n-1)//2);rows={};pairs=[(i,j) for i in range(n) for j in range(i+1,n)];idx={p:k for k,p in enumerate(pairs)}
 for code in range(N):
  a=tournament(n,code)
  if strong(a):rows[code]=(H(a),a)
 p={c:c for c in rows}
 def get(c):
  while p[c]!=c:p[c]=p[p[c]];c=p[c]
  return c
 def union(c,d):
  if d in rows and d>c and abs(rows[d][0]-rows[c][0])<=2:p[get(c)]=get(d)
 def delete(a,r):
  vs=[i for i in range(n) if i!=r]
  return tuple(sum(1<<j for j,v in enumerate(vs) if a[u]>>v&1) for u in vs)
 for c,(h,a) in rows.items():
  for r in range(n):
   if not strong(delete(a,r)):continue
   es=[(v,idx[tuple(sorted([r,v]))]) for v in range(n) if v!=r]
   for v,k in es:union(c,c^(1<<k))
   for v,k in es:
    for w,l in es:
     if k<l and ((a[r]>>v)^(a[r]>>w))&1:union(c,c^(1<<k)^(1<<l))
 comp={}
 for c,(h,a) in rows.items():comp.setdefault(get(c),set()).add(h)
 intervals={tuple(sorted(v)) for v in comp.values()};maximal=sorted(v for v in intervals if not any(set(v)<set(w) for w in intervals))
 print(json.dumps({'n':n,'strong_labelled':len(rows),'components':len(comp),'distinct_images':len(intervals),'maximal_images':maximal,'seconds':time.time()-t}))
