"""Exact small controls of the seam-fan response law; no spectrum census."""
import itertools,json
from functools import lru_cache

def rows(n,mask):
 r=[0]*n
 for bit,(i,j) in enumerate(itertools.combinations(range(n),2)):
  u,v=(i,j) if mask>>bit&1 else (j,i);r[u]|=1<<v
 return r

def profile(r):
 n=len(r)
 @lru_cache(None)
 def end(m,v):
  if m==1<<v:return 1
  return sum(end(m^(1<<v),u) for u in range(n) if m>>u&1 and u!=v and r[u]>>v&1)
 @lru_cache(None)
 def h(m):return sum(end(m,v) for v in range(n) if m>>v&1) if m else 1
 @lru_cache(None)
 def start(m,a):
  if m==1<<a:return 1
  return sum(start(m^(1<<a),v) for v in range(n) if v!=a and m>>v&1 and r[a]>>v&1)
 def l(m,a):
  z=m;out=0
  while z:
   if z>>a&1:out+=start(z,a)*h(m^z)
   z=(z-1)&m
  return out
 return h,end,l

def signature(r,a,t):
 n=len(r);h,e,l=profile(r);Q=E=L=D=0
 for colors in itertools.product(range(t),repeat=n):
  ms=[0]*t
  for v,c in enumerate(colors):ms[c]|=1<<v
  root=colors[a];factor=1
  for c in range(t):
   if c!=root:factor*=h(ms[c])
  Q+=factor*h(ms[root]);E+=factor*e(ms[root],a);L+=factor*l(ms[root],a)
  D+=factor*h(ms[root]^(1<<a))
 return Q,E,L,D

def fan(r,a,k):
 r=r[:];n=len(r)
 for j in range(k):
  v=n+j;r.append(1<<a)
  for u in range(v):
   if u!=a:r[u]|=1<<v
 return r

def run():
 checks=0
 for n in range(1,5):
  for mask in range(1<<(n*(n-1)//2)):
   r=rows(n,mask)
   for a in range(n):
    for t in range(1,4):
     Q,E,L,D=signature(r,a,t)
     for k in range(4):
      out=signature(fan(r,a,k),a,t)[0]
      pred=(L-D)*(t+1)**k+(Q-E-L+2*D)*t**k+(E-D)*(t-1)**k
      assert out==pred,(n,mask,a,t,k,out,pred)
      checks+=1
 return {'status':'PASS','checks':checks,'universe':'all labelled tournaments n=1..4, every root, t=1..3, fan k=0..3','triangle_H':[signature(fan([2,4,1],0,k),0,1)[0] for k in range(7)]}
if __name__=='__main__': print(json.dumps(run(),indent=2))
