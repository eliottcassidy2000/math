"""Exact transitive-core charge fibers; no unbounded coverage inference."""
from itertools import product
from collections import defaultdict,Counter
import json

def f(d): return 1 if d==1 else 1<<(d-2)
def count(n,a):
 return 1+2*sum(f(j-i) for i in range(n) for j in range(i+1,n) if a>>i&1 and not(a>>j&1))
def charge(n,a):return sum(((1<<(n-1-i))-(1<<i)) for i in range(n) if a>>i&1)
def signature(n,a):return tuple(((a>>i)&1)-((a>>(n-1-i))&1) for i in range(n//2))
def graph(n):
 groups=defaultdict(list)
 for a in range(1<<n):groups[charge(n,a)].append(a)
 best=[];hist=Counter();maxwidth=0;example=None
 for c,aa in groups.items():
  if len({signature(n,a) for a in aa})!=1:raise AssertionError('charge loses ternary digits')
  aset=set(aa);hs={a:count(n,a) for a in aa}; seen=set()
  moves=[(1<<i)|(1<<(n-1-i)) for i in range(n//2)]
  if n%2:moves+=[1<<(n//2)]
  for a in aa:
   if a in seen:continue
   seen.add(a);q=[a];states=[]
   while q:
    v=q.pop();states.append(v)
    for m in moves:
     u=v^m
     if u in aset and u not in seen and abs(hs[v]-hs[u])<=2:seen.add(u);q.append(u)
   vv=sorted({hs[v] for v in states});w=vv[-1]-vv[0]
   assert vv==list(range(vv[0],vv[-1]+1,2))
   hist[w]+=1
   if w>maxwidth:maxwidth=w;example={'charge':c,'cuts':states,'counts':vv,'signature':signature(n,a)}
 return {'core_order':n,'fibers':len(groups),'fiber_size_histogram':dict(Counter(map(len,groups.values()))),'component_width_histogram':dict(sorted(hist.items())),'max_width':maxwidth,'example':example}
if __name__=='__main__':
 print(json.dumps([graph(n) for n in range(2,15)],indent=2))
