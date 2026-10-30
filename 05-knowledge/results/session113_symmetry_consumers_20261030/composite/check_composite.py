"""Exact generated groups, parity constraints, and tuple moments. Stdlib only."""
from math import gcd,factorial
from itertools import product
from collections import deque
import json
CHECKS=0
def need(ok,msg):
 global CHECKS
 CHECKS+=1
 if not ok:raise RuntimeError(msg)
def compose(a,b):return tuple(a[b[i]] for i in range(len(a)))
def generated(gs):
 n=len(gs[0]);identity=tuple(range(n));seen={identity};queue=deque([identity])
 while queue:
  x=queue.popleft()
  for g in gs:
   y=compose(g,x)
   if y not in seen:seen.add(y);queue.append(y)
 return seen
def generators(r,m):
 n=r*m;a=tuple((x+1)%n for x in range(n));v=tuple((x+r)%n if x%r==0 else x for x in range(n));vi=[0]*n
 for x,y in enumerate(v):vi[y]=x
 b=compose(v,compose(a,tuple(vi)))
 return a,v,b

def bank_dimensions(r,m):
 a,v,b=generators(r,m);n=r*m;edges=[(x,y) for x in range(n) for y in range(x+1,n)];ix={e:i for i,e in enumerate(edges)}
 pa=list(range(len(edges)));par=[0]*len(edges)
 def find(x):
  if pa[x]!=x:
   p=pa[x];pa[x]=find(p);par[x]^=par[p]
  return pa[x]
 def add(i,j,s):
  ri,rj=find(i),find(j);d=par[i]^par[j]^s
  if ri==rj:need(not d,'signed equations inconsistent')
  else:pa[ri]=rj;par[ri]=d
 for perm in (a,b):
  for (x,y),i in ix.items():
   u,w=perm[x],perm[y];add(i,ix[tuple(sorted((u,w)))],int(u>w))
 roots={find(i) for i in range(len(edges))};dim=len(roots);pins={}
 for x in range(n-1):
  i=ix[(x,x+1)];root=find(i);value=1^par[i]
  need(root not in pins or pins[root]==value,'shared path illegal');pins[root]=value
 need(dim==(r+m-2)//2,'common direction dimension')
 need(dim-len(pins)==(r+m-4)//2,'native chart intersection dimension')
 # All source points for small common bank are exactly circulant products.
 if n<=25:
  for bits in product(range(2),repeat=dim):
   vals=dict(zip(sorted(roots),bits));edge=[vals[find(i)]^par[i] for i in range(len(edges))]
   def arc(x,y):return edge[ix[(x,y)]] if x<y else 1-edge[ix[(y,x)]]
   for x,y in edges:
    i,j=x%r,y%r;u,w=x//r,y//r
    if i==j:
     d=(w-u)%m;want=arc(0,r*d)
    else:want=arc(0,(j-i)%r)
    need(arc(x,y)==want,'bank not lexicographic circulant product')
 return {'r':r,'m':m,'common_dimension':dim,'native_dimension':dim-len(pins)}
def stir(t,j):
 s=[1]+[0]*j
 for _ in range(t):s=[0]+[i*s[i]+s[i-1] for i in range(1,j+1)]
 return s[j]
def falling(n,j):return factorial(n)//factorial(n-j)
def tuple_formula(r,m,t):
 if not t:return 1,1
 w=sum(falling(r,j)*stir(t,j)*m**(t-j) for j in range(1,min(r,t)+1))//r
 extra=0 if t<r else (gcd(r,m)-1)*factorial(r-1)*stir(t,r)*m**(t-r)
 return w,w+extra

def group_case(r,m):
 a,v,b=generators(r,m);G=generated((a,b));W=generated((a,v));g=gcd(r,m);n=r*m
 need(len(G)==r*m**r//g,'generated G order');need(len(W)==r*m**r,'generated W order');need(a!=b,'parents identical')
 kernel=set()
 for perm in W:
  t=perm[0]%r;z=[]
  for i in range(r):
   zi=perm[i]//r;z.append(zi)
   need(all(perm[i+r*j]==(i+t)%r+r*((j+zi)%m) for j in range(m)),'wreath coordinates')
  if (sum(z)-t)%g==0:kernel.add(perm)
 need(kernel==G,'kernel differs from generated group')
 fw=[sum(i==y for i,y in enumerate(p)) for p in W];fg=[sum(i==y for i,y in enumerate(p)) for p in G]
 counts=[]
 for t in range(r+3):
  w=sum(x**t for x in fw)//len(W);gg=sum(x**t for x in fg)//len(G)
  need((w,gg)==tuple_formula(r,m,t),'tuple moment formula')
  counts.append([w,gg])
 # A named transversal and one local shift detect the missing phase.
 if g>1:
  trans=tuple(range(r));target=tuple(v[i] for i in trans)
  need(any(tuple(p[i] for i in trans)==target for p in W),'W phase transport absent')
  need(not any(tuple(p[i] for i in trans)==target for p in G),'G sees erased phase')
  for omit in range(r):
   sub=tuple(i for i in trans if i!=omit);tar=tuple(v[i] for i in sub)
   need(any(tuple(p[i] for i in sub)==tar for p in G),'proper subtransversal distinguishes phase')
 return {'r':r,'m':m,'G_order':len(G),'W_order':len(W),'tuple_counts_W_G':counts}

def main():
 banks=[bank_dimensions(r,m) for r,m in [(3,3),(3,5),(3,7),(5,3),(5,5),(5,7),(7,7),(9,3),(9,9)]]
 groups=[group_case(r,m) for r,m in [(3,3),(3,5),(5,3),(5,5)]]
 print(json.dumps({'status':'PASS','checks':CHECKS,'banks':banks,'groups':groups},indent=2,sort_keys=True))
if __name__=='__main__':main()
