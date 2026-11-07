#!/usr/bin/env python3
"""Two explicit order-six controls; no repository or session-worker imports."""
from itertools import permutations,combinations,product
from collections import Counter
from fractions import Fraction as F
from math import factorial
import json

def check(ok,msg):
 if not ok:raise RuntimeError(msg)

def paths(a,mask):
 vs=[i for i in range(len(a)) if mask>>i&1]
 return sum(all(a[u]>>v&1 for u,v in zip(p,p[1:])) for p in permutations(vs))

def odd_cycles(a):
 out=[]
 for k in (3,5):
  for vs in combinations(range(6),k):
   for tail in permutations(vs[1:]):
    p=(vs[0],)+tail
    if all(a[u]>>v&1 for u,v in zip(p,p[1:]+p[:1])):out.append(p)
 return out

def poly_from_hist(hist,n):
 out=[F(0)]*(n+1)
 for f,count in hist.items():
  p=[F(1)]
  for j in range(n):
   q=[F(0)]*(len(p)+1)
   for i,c in enumerate(p):q[i]+=(f-j)*c;q[i+1]+=c
   p=q
  for i,c in enumerate(p):out[i]+=F(count,factorial(n))*c
 check(all(c.denominator==1 for c in out),'integral coefficients')
 return [int(c) for c in out]

def direct_cover(a,q):
 n=len(a); h=[paths(a,m) for m in range(1<<n)]
 total=0
 for bins in product(range(q),repeat=n):
  masks=[0]*q
  for v,b in enumerate(bins):masks[b]|=1<<v
  z=1
  for m in masks:z*=h[m]
  total+=z
 return total

def analyze(a):
 n=len(a)
 check(all(not(a[i]>>i&1) for i in range(n)),'no loops')
 check(all(((a[i]>>j)&1)+((a[j]>>i)&1)==1 for i in range(n) for j in range(i+1,n)),'tournament')
 hist=Counter(sum(a[u]>>v&1 for u,v in zip(p,p[1:])) for p in permutations(range(n)))
 polynomial=poly_from_hist(hist,n)
 cyc=odd_cycles(a);pairs=[(i,j) for i,j in combinations(range(len(cyc)),2) if not set(cyc[i])&set(cyc[j])]
 # At order six every packing has at most two cycles, so this is complete.
 profile=[1,len(cyc),len(pairs)]
 while profile[-1]==0:profile.pop()
 parent=list(range(len(cyc)))
 def find(i):
  while parent[i]!=i:i=parent[i]
  return i
 for i,j in pairs:parent[find(i)]=find(j)
 components=len({find(i) for i in range(len(cyc))})
 vertices=len(cyc);edges=2*len(pairs)
 b1=edges-vertices+components
 for q in range(5):
  check(direct_cover(a,q)==sum(c*q**i for i,c in enumerate(polynomial)),('direct bins',a,q))
 return {'rows':a,'H':paths(a,63),'forward_adjacency_histogram':[hist[f] for f in range(6)],'P_coefficients_ascending':polynomial,'odd_cycles':cyc,'disjoint_cycle_pairs':pairs,'packing_coefficients':profile,'cell_counts':[vertices]+([edges] if edges else []),'betti_Q_and_F2':[components,b1]}

def main():
 out=[analyze(a) for a in ([58,4,1,6,14,30],[38,8,2,5,15,30])]
 check(all(x['H']==17 for x in out),'H')
 check(all(x['P_coefficients_ascending']==[0,0,8,0,8,0,1] for x in out),'P')
 check(all(x['forward_adjacency_histogram']==[17,105,238,238,105,17] for x in out),'histogram')
 check([x['packing_coefficients'] for x in out]==[[1,8],[1,6,1]],'packing profiles')
 check([x['betti_Q_and_F2'] for x in out]==[[8,0],[5,1]],'topology')
 print(json.dumps({'status':'PASS','scope':'Two explicit order-six tournaments; exact permutations, odd cycles, bin assignments q=0..4 and polynomial expansion.','controls':out},indent=2))
if __name__=='__main__':main()
