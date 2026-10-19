#!/usr/bin/env python3
"""Independent reconstruction: order-ideal level profiles + core-word interleavings.
Does not import any author code, profile arrays, counts, or numerical certificates.
"""
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import product
import hashlib,json
import numpy as np

Q=(0,0,0,6,0,19,2,18,86,151)
ROOTS=(0,3)
TREES=(((3,7),(0,2)),((5,8),(0,6)),((6,7),(5,8)))
SAFE=((0,0),(0,1),(1,1),(1,2),(1,3),(2,2),(2,3))
CERTS=(((1,0,2),(1,2,0)),((1,0,2),(1,1,0),(2,2,1)),((1,0,3),(1,2,0)),((1,0,3),(1,1,0),(2,2,1)))

def req(test,message):
 if not test:raise ValueError(message)

def lt(a,b):return bool(Q[b] & (1<<a))
SLACK=tuple(7-sum(a!=b and not lt(a,b) and not lt(b,a) for b in range(10)) for a in range(10))
# Independently enumerate profiles by nested order ideals or order filters.
def profiles(upper):
 levels=[]
 for i in range(1,1+max(SLACK)):
  choices=[]
  for bits in range(1<<10):
   support=[q for q in range(10) if bits>>q&1]
   if len(support)>7 or any(q in ROOTS or SLACK[q]<i for q in support):continue
   if upper:
    valid=all(not lt(a,b) or not(bits>>a&1) or (bits>>b&1) for a in range(10) for b in range(10))
   else:
    valid=all(not lt(a,b) or not(bits>>b&1) or (bits>>a&1) for a in range(10) for b in range(10))
   if valid:choices.append(bits)
  levels.append(choices)
 out=[]
 def visit(i,last,p):
  if i==len(levels):out.append(tuple(p));return
  for bits in levels[i]:
   if bits & ~last:continue
   visit(i+1,bits,[p[q]+bool(bits>>q&1) for q in range(10)])
 visit(0,(1<<10)-1,[0]*10)
 return sorted(out)
L,U=profiles(False),profiles(True)
req((len(L),len(U))==(46,48),'profile counts')
req(len(set(L))==len(L) and len(set(U))==len(U),'unique profiles')
PAIRS=[(i,j) for i,a in enumerate(L) for j,b in enumerate(U) if all(a[q]+b[q]<=SLACK[q] for q in range(10))]
req(len(PAIRS)==1764,'compatible profile count')

# Independently construct every ambient model as a set relation and close it.
# Check that closure does not silently create any extra core or exterior edge.
ambient_orders=[]
for li,ui in PAIRS:
 a,b=L[li],U[ui]
 low=max(a);up=max(b);order=10+low+up
 edges={(x,y) for x in range(10) for y in range(10) if lt(x,y)}
 lower=list(range(10,10+low));upper=list(range(10+low,order))
 edges.update((x,y) for chain in (lower,upper) for x in chain for y in chain if x<y)
 edges.update((x,y) for x in lower for y in upper)
 edges.update((x,q) for x in lower for q in range(10) if x-10<low-a[q])
 edges.update((q,y) for y in upper for q in range(10) if y-10-low>=b[q])
 intended=edges.copy()
 for mid in range(order):
  incoming=[x for x in range(order) if (x,mid) in edges]
  outgoing=[y for y in range(order) if (mid,y) in edges]
  edges.update((x,y) for x in incoming for y in outgoing)
 req(edges==intended,'intended ambient relation is already transitive')
 req(all((x,x) not in edges for x in range(order)),'ambient acyclic')
 incomparable=lambda x,y:x!=y and (x,y) not in edges and (y,x) not in edges
 req(all(sum(incomparable(x,y) for y in range(order))<=7 for x in range(order)),'ambient degree bound')
 fullcore={x for x in range(order) if x in ROOTS or any(incomparable(x,r) for r in ROOTS)}
 req(fullcore==set(range(10)),'ambient full core')
 pred=lambda y:{x for x in range(order) if (x,y) in edges}
 succ=lambda x:{y for y in range(order) if (x,y) in edges}
 req(pred(0)<=pred(3) and succ(3)<=succ(0),'ambient critical roots')
 ambient_orders.append(order)

# Group all core extensions only by boundary words. Each group records its
# entire joint answer-cell counts, not merely individual probabilities.
GROUPS={}
core_extensions=0
def visit(word,placed):
 global core_extensions
 if len(word)==10:
  core_extensions+=1
  pos={q:i for i,q in enumerate(word)}
  first=min(pos[r] for r in ROOTS);last=max(pos[r] for r in ROOTS)
  key=(word[:first],tuple(reversed(word[last+1:])))
  row=GROUPS.setdefault(key,[0]*(1+4*len(TREES)))
  row[0]+=1
  for t,(ab,cd) in enumerate(TREES):
   yes1=pos[ab[0]]<pos[ab[1]];yes2=pos[cd[0]]<pos[cd[1]]
   cell=2*(not yes1)+(not yes2)
   row[1+4*t+cell]+=1
  return
 for q in range(10):
  if not(placed>>q&1) and not(Q[q]&~placed):visit(word+(q,),placed|1<<q)
visit((),0)

# Count nonincreasing nonnegative integer gap words by a direct transfer table.
# This is intentionally separate from the recursive G formula in discovery.
@lru_cache(None)
def gaps(caps):
 if not caps:return 1
 row=[1]*(caps[0]+1)
 for cap in caps[1:]:
  row=[sum(row[v:]) for v in range(cap+1)]
 return sum(row)
keys=sorted(GROUPS)
LC=np.array([[gaps(tuple(a[q] for q in key[0])) for key in keys] for a in L],dtype=np.int64)
UC=np.array([[gaps(tuple(b[q] for q in key[1])) for key in keys] for b in U],dtype=np.int64)
EVENTS=np.array([GROUPS[k] for k in keys],dtype=np.int64)
WEIGHTS=np.stack([LC[i]*UC[j] for i,j in PAIRS])
COUNTS=WEIGHTS@EVENTS
# An explicit crude arithmetic bound proves every dot product is exact int64.
req(int(WEIGHTS.max())*core_extensions<2**63,'dot product overflow guard')
E=COUNTS[:,0];CELLS=COUNTS[:,1:].reshape(-1,len(TREES),4)
req(bool(np.all(CELLS>=0)) and bool(np.all(CELLS.sum(axis=2)==E[:,None])),'cell partitions')
BRANCHES=np.repeat(CELLS.reshape(-1,len(TREES),2,2).sum(axis=3),2,axis=2)
M=E[:,None,None]-BRANCHES-CELLS
mins_safe=[];mins_cert=[]
for t,c in SAFE:
 req(bool(np.all(200*M[:,t,c]>E)),'safe inequality')
 mins_safe.append(min(Fraction(int(v),int(e)) for v,e in zip(M[:,t,c],E)))
for cert in CERTS:
 total=sum(co for co,t,c in cert)
 vals=sum(co*M[:,t,c] for co,t,c in cert)
 req(bool(np.all(200*vals>total*E)),'positive certificate inequality')
 mins_cert.append(min(Fraction(int(v),total*int(e)) for v,e in zip(vals,E)))
patterns=0
for cells in product(range(4),repeat=len(TREES)):
 if any(cells[t]==c for t,c in SAFE):continue
 patterns+=1
 req(any(all(cells[t]==c for co,t,c in cert) for cert in CERTS),'uncovered failure pattern')
req(patterns==4,'residual pattern count')

# Independently exercise the stated deterministic selection rule on every atom.
selection_counts=[0,0,0]
for z,e in enumerate(E):
 if 200*M[z,2,0]<=e:choice=0
 elif 200*M[z,2,1]>e:choice=2
 elif 200*M[z,1,0]>e:choice=1
 else:choice=0
 req(bool(np.all(200*M[z,choice,:]>e)),'deterministic selection rule')
 selection_counts[choice]+=1
 # The selected first question is strictly inside the stated balanced interval.
 first=int(CELLS[z,choice,0]+CELLS[z,choice,1])
 req(101*int(e)<300*first<199*int(e),'balanced first-comparison corollary')

# Reconstruct the author's publicly specified table format to compare independent
# integer counts across EVERY profile and every joint cell, using one digest.
digest=hashlib.sha256()
for z,(i,j) in enumerate(PAIRS):
 line=json.dumps([list(L[i]),list(U[j]),int(E[z]),CELLS[z].tolist()],separators=(',',':'))+'\n'
 digest.update(line.encode())
expected='2f83a7b75f377b5160c37b83b30f64da1ec3207531dcd575983a5d48da9c767d'
req(digest.hexdigest()==expected,'independent all-count table matches direct ideal DP')
result={'status':'PASS','method':'nested ideals/filters; integer core-word gap transfer counts','core_extensions':core_extensions,'boundary_groups':len(keys),'profile_counts':[len(L),len(U)],'compatible_atoms':len(PAIRS),'ambient_orders':[min(ambient_orders),max(ambient_orders)],'count_table_sha256':digest.hexdigest(),'safe_minima':list(map(str,mins_safe)),'certificate_minima':list(map(str,mins_cert)),'residual_failure_patterns':patterns,'selection_counts':selection_counts,'max_gap_weight':int(WEIGHTS.max()),'max_atom_extensions':int(E.max()),'universal_strict_margin':'1/200'}
print(json.dumps(result,indent=2))
