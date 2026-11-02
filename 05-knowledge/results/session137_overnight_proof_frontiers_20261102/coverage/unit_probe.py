from explore import *
import json,time
for n in [5,6]:
 t=time.time();N=1<<(n*(n-1)//2);rows={};pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
 for code in range(N):
  a=tournament(n,code)
  if strong(a):rows[code]=H(a)
 p={c:c for c in rows}
 def get(c):
  while p[c]!=c:p[c]=p[p[c]];c=p[c]
  return c
 for c,h in rows.items():
  for k,ij in enumerate(pairs):
   d=c^(1<<k)
   if d in rows and d>c and abs(rows[d]-h)<=2:p[get(c)]=get(d)
 comp={}
 for c,h in rows.items():comp.setdefault(get(c),set()).add(h)
 intervals={tuple(sorted(v)) for v in comp.values()}
 maximal=sorted(v for v in intervals if not any(set(v)<set(w) for w in intervals))
 print(json.dumps({'n':n,'strong_labelled':len(rows),'components':len(comp),'distinct_images':len(intervals),'maximal_images':maximal,'seconds':time.time()-t}))
