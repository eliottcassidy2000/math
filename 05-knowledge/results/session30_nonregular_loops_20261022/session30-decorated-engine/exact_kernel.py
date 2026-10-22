from itertools import product
from math import factorial
from fractions import Fraction

def matchings(items):
 if not items: yield ();return
 a,*rest=items
 for m in matchings(rest):yield m
 for j,b in enumerate(rest):
  for m in matchings(rest[:j]+rest[j+1:]):yield ((a,b),)+m

def kernel(edges):
 inc={}
 for e,(a,b) in enumerate(edges):
  inc.setdefault(a,[]).append(2*e);inc.setdefault(b,[]).append(2*e+1)
 opts=[]
 for hs in inc.values():
  oo=[]
  for pairs in matchings(hs):
   s=len(hs)-len(pairs);mu=(-1)**(s-1)*factorial(s-1)
   oo.append((pairs,mu))
  opts.append(oo)
 result=0
 for choices in product(*opts):
  nxt={};w=1
  for pairs,mu in choices:
   w*=mu
   for a,b in pairs:nxt[a]=b;nxt[b]=a
  endpoints=[h for h in range(2*len(edges)) if h not in nxt]
  visited=set();bad=False
  for start in endpoints:
   if start//2 in visited:continue
   h=start;length=0;back=0
   while True:
    visited.add(h//2);length+=1;back+=h%2
    h=h^1
    if h not in nxt:break
    h=nxt[h]
   if length%2:bad=True;break
   w*=2*(-1)**back
  if not bad and len(visited)==len(edges):result+=w
 return result

def lollipop(a,b):
 edges=[(i,(i+1)%a) for i in range(a)]
 last=0
 for j in range(b):
  new=a+j;edges.append((last,new));last=new
 return edges
if __name__=='__main__':
 for a in [2,4,6,8]:
  for b in [1,2,3,4,5,6]:
   print(a,b,Fraction(kernel(lollipop(a,b)),2))
