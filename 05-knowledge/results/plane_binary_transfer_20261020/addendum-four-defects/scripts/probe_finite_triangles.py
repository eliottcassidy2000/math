"""Finite exact NAE solver for anisotropic norm equilateral triangle hypergraphs."""
from collections import defaultdict
import json,sys,time

def field(p):
 d=next(d for d in range(2,p) if pow(d,(p-1)//2,p)==p-1)
 def add(z,w):return ((z[0]+w[0])%p,(z[1]+w[1])%p)
 def mul(z,w):return ((z[0]*w[0]+d*z[1]*w[1])%p,(z[0]*w[1]+z[1]*w[0])%p)
 norm=lambda z:(z[0]*z[0]-d*z[1]*z[1])%p
 pts=[(x,y) for x in range(p) for y in range(p)]
 U=[z for z in pts if norm(z)==1]
 zeta=next(z for z in U if add(mul(z,z),(1,0))==z)
 ix={z:i for i,z in enumerate(pts)}
 edges=sorted(set(tuple(sorted([ix[a],ix[add(a,u)],ix[add(a,mul(u,zeta))]])) for a in pts for u in U))
 return d,zeta,pts,U,edges

def nae_solve(n,edges,limit=1000000):
 incid=defaultdict(list)
 for e in edges:
  for v in e:incid[v].append(e)
 nodes=0
 def dfs(a):
  nonlocal nodes
  nodes+=1
  if nodes>limit: raise TimeoutError
  while True:
   changes=[]
   for e in edges:
    colors=[a[v] for v in e if a[v]>=0]
    if len(colors)==3 and len(set(colors))==1:return None
    if len(colors)==2 and colors[0]==colors[1]:
     v=next(v for v in e if a[v]<0); changes.append((v,1-colors[0]))
   if not changes:break
   for v,c in changes:
    if a[v]>=0 and a[v]!=c:return None
    a[v]=c
  rem=[v for v in range(n) if a[v]<0]
  if not rem:return a
  def score(v):
   return sum(10 if sum(a[w]>=0 for w in e)==1 else 1 if all(a[w]<0 for w in e) else 0 for e in incid[v])
  v=max(rem,key=score)
  for c in [0,1]:
   b=a[:];b[v]=c
   ans=dfs(b)
   if ans is not None:return ans
  return None
 a=[-1]*n;a[0]=0
 try:ans=dfs(a);status='SAT' if ans else 'UNSAT'
 except TimeoutError:ans=None;status='LIMIT'
 return status,ans,nodes

if __name__=='__main__':
 out=[]
 for p in map(int,sys.argv[1:] or ['5','11']):
  d,zeta,pts,U,edges=field(p);t=time.time()
  status,ans,nodes=nae_solve(p*p,edges)
  row=dict(p=p,d=d,zeta=zeta,vertices=p*p,unit_count=len(U),triangle_count=len(edges),status=status,nodes=nodes,elapsed=round(time.time()-t,4),coloring=ans)
  if ans is not None:assert all(len({ans[v] for v in e})==2 for e in edges)
  out.append(row);print(json.dumps(row),flush=True)
