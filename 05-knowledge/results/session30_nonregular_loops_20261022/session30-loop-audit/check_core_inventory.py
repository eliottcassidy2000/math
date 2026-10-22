"""Independent degree-sequence-first multigraph census and automorphism audit."""
from itertools import combinations_with_replacement,permutations
from math import factorial,prod
from collections import Counter
from pathlib import Path
import json,sys
HERE=Path(__file__).resolve().parent

def canonical(V,edges):
 return min(tuple(sorted(tuple(sorted((p[a],p[b])))for a,b in edges))for p in permutations(range(V)))
def enumerate_cores(w):
 out={}
 for V in range(1,2*w+1):
  slots=[(a,b)for a in range(V)for b in range(a,V)]
  for E in range(1,V+w+1):
   for ds in combinations_with_replacement([1]+list(range(3,2*E+1)),V):
    if sum(ds)!=2*E or E-V+ds.count(1)!=w:continue
    rem=list(ds);edges=[]
    def rec(j):
     if j==len(slots):
      if any(rem):return
      reach={0}
      while True:
       r=reach|{a for a,b in edges if b in reach}|{b for a,b in edges if a in reach}
       if r==reach:break
       reach=r
      if len(reach)!=V:return
      c=canonical(V,edges);mult=Counter(c)
      vertex_aut=sum(Counter(tuple(sorted((p[a],p[b])))for a,b in c)==mult for p in permutations(range(V)))
      full=vertex_aut*prod(factorial(m)for m in mult.values())*2**sum(a==b for a,b in c)
      out[V,c]=(E-V,ds.count(1),full);return
     a,b=slots[j];mx=rem[a]//2 if a==b else min(rem[a],rem[b])
     for k in range(mx+1):
      rem[a]-=k;rem[b]-=k;edges.extend([(a,b)]*k)
      rec(j+1)
      if k:del edges[-k:]
      rem[a]+=k;rem[b]+=k
    rec(0)
 return out
sys.path.insert(0,str(HERE.parent/'session30-decorated-engine'))
from decorated import decorated_cores
report=[]
for w in (1,2):
 ref=enumerate_cores(w);engine={}
 for d in decorated_cores(w):
  key=d['V'],canonical(d['V'],d['edges']);engine[key]=d['excess'],d['leaves'],d['automorphisms']
 if ref!=engine:raise RuntimeError((w,ref.keys()-engine.keys(),engine.keys()-ref.keys(),[(k,ref[k],engine[k])for k in ref.keys()&engine.keys()if ref[k]!=engine[k]]))
 row={'weight':w,'count':len(ref),'by_type':{str(k):v for k,v in Counter((e,l)for e,l,a in ref.values()).items()},'automorphisms':'ALL_PASS'};report.append(row);print(row)
(HERE/'core-inventory-results.json').write_text(json.dumps(report,indent=2)+'\n')
