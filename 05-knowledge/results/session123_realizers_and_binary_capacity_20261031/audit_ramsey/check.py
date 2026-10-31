"""Independent Session123 audit: maximal-port enumeration, no producer imports."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
import json
import numpy as np

def require(t,m):
 if not t: raise RuntimeError(m)
@lru_cache(None)
def weight(a):
 return sum((Q(3,2)**i.bit_count() for i in range(a)),Q(0))
def kernel(a,b):
 return sum(2**(i&j).bit_count() for i in range(a) for j in range(b))
def L(a): return 1<<(a-1).bit_length()
def H(a,b):
 l=L(a)
 return int(weight(a)*l)*(b//l)+kernel(a,b%l)
@lru_cache(None)
def maximal(b):
 out=set()
 for u in range(1,b//2+1):
  for v in range(u,b-u+1):
   p=(u,b-u,u,v,b-v)
   for k in range(5):out.add(p[k:]+p[:k])
 return tuple(sorted(out))
def orbit(p):
 return min(q[k:]+q[:k] for q in (p,p[::-1]) for k in range(5))
def exhaustive(b):
 return set(p for p in product(range(1,b),repeat=5) if all(p[i]+p[(i+1)%5]<=b and (p[i]+p[(i+1)%5]==b or p[i]+p[(i-1)%5]==b) for i in range(5)))
def ceil(q):return -(-q.numerator//q.denominator)
@lru_cache(None)
def root(a,b):
 if b<2:return (0,None)
 pp=np.array(maximal(b),dtype=int)
 best=0;witness=None
 for p in sorted(set(map(orbit,maximal(a)))):
  reordered=[p[i] for i in (0,2,4,1,3)]
  vals=np.zeros(len(pp),dtype=np.int64)
  for j,x in enumerate(reordered):
   lookup=np.array([H(x,z) for z in range(b)],dtype=np.int64)
   vals+=lookup[pp[:,j]]
  k=int(np.argmax(vals))
  if int(vals[k])>best:best=int(vals[k]);witness=[list(p),[int(x) for x in pp[k]]]
 return best,witness

def audit(a):
 l=L(a);s=weight(a);d=max(s*r-H(a,r) for r in range(l));cut=l;records=[]
 require(d>0,('deficit',a))
 for p in sorted(set(map(orbit,maximal(a)))):
  c=[weight(p[i]) for i in (0,2,4,1,3)]
  half=sum(c)/2;m=max([half]+c+[c[i]+c[(i+2)%5] for i in range(5)])
  if half==s:
   # Solve the dual using a forward rather than backward alternating sum.
   y=[sum((-1)**k*c[(i-k)%5] for k in range(5))/2 for i in range(5)]
   require(all(y[i]+y[(i-1)%5]==c[i] for i in range(5)),('dual',a,p))
   require(min(y)>0,('positive',a,p))
   require(all((l//2)%L(x)==0 for x in p),('child periods',a,p))
   k=ceil(d/min(y))+l
   typ='sharp'
  else:
   require(m<s,('gap',a,p,m,s))
   k=ceil(d/(s-m))-1;typ='nonsharp'
  cut=max(cut,k);records.append({'p':p,'kind':typ,'cutoff':k})
 for x in range(1,a):
  gap=s-weight(x)-weight(a-x);require(gap>0,('join gap',a,x))
  cut=max(cut,ceil(d/gap)-1)
 seeds=[]
 for b in range(1,cut+1):
  n,wit=root(a,b)
  joins=[H(x,b)+H(a-x,b) for x in range(1,a)]
  require(max([n]+joins)<=H(a,b),('upper',a,b,n,joins,H(a,b)))
  if b<=l:
   unions=[H(a,x)+H(a,b-x) for x in range(1,b)]
   v=max([1,n]+joins+unions)
   require(v==H(a,b),('seed',a,b,v,H(a,b)))
   seeds.append({'b':b,'order':v,'c5_witness':wit if n==v else None,'join':next((x for x in range(1,a) if H(x,b)+H(a-x,b)==v),None),'union':next((x for x in range(1,b) if H(a,x)+H(a,b-x)==v),None)})
 for r in range(l):
  for t in range(l):require(H(a,r)+H(a,t)<=H(a,r+t),('union',a,r,t))
 irr=[b for b in range(1,l+1) if all(H(a,b)>H(a,x)+H(a,b-x) for x in range(1,b))]
 require(irr==list(range(1,l+1)),('irreducibles',a,irr))
 return {'A':a,'L':l,'prefix':[H(a,b) for b in range(l+1)],'cutoff':cut,'full_patterns':len(maximal(a)),'orbit_patterns':len(records),'records':records,'seeds':seeds,'irreducibles':irr}
if __name__=='__main__':
 for b in range(2,11):require(set(maximal(b))==exhaustive(b),('exhaustive maximal',b))
 # Orbit reduction preserves feasible clique and complementary independence constraints.
 for a in range(2,17):
  for p in maximal(a):require(orbit(p) in maximal(a),('orbit',a,p))
 results=[]
 for a in range(2,17):
  result=audit(a);results.append(result)
 print(json.dumps({'status':'PASS','algorithm':'tight-edge-cover maximal-port enumeration','exhaustive_controls':'budgets2..10','rows':results},sort_keys=True))
