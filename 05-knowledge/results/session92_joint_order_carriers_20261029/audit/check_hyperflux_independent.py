from check_flux_independent import rank,require
from itertools import combinations,product
import json
cases=0
for m in range(2,5):
 candidates=[frozenset(c) for k in range(2,m+1) for c in combinations(range(m),k)]
 families=([set(candidates[j] for j in range(len(candidates)) if mask>>j&1) for mask in range(1<<len(candidates))] if m<=3 else [{F} for F in candidates]+[set(candidates)])
 states=list(product((-1,0,1),repeat=m));idx={s:i for i,s in enumerate(states)}
 for family in families:
  L=[]
  for i in range(m):
   for others in product((-1,0,1),repeat=m-1):
    s=list(others);s.insert(i,0);t=s.copy();t[i]=1
    parity=sum(all(s[j]==1 for j in F if j!=i) for F in family if max(F)==i)%2
    row=[0]*len(states);row[idx[tuple(s)]]=1;row[idx[tuple(t)]]=(-1)**parity;L.append(row)
  inds=[S for k in range(m+1) for c in combinations(range(m),k) for S in [set(c)] if not any(F<=S for F in family)]
  for S in inds:
   f=[(-1)**sum(x==1 for x in s) if {i for i,x in enumerate(s) if x>=0}==S else 0 for s in states]
   require(all(sum(a*b for a,b in zip(row,f))==0 for row in L),'integer kernel')
  require(rank(L)==3**m-len(inds),'hypergraph rank');cases+=1
print(json.dumps({'hypergraph_cases':cases,'status':'pass'}))
