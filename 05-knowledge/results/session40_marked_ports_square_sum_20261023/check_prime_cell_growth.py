from itertools import permutations, product
from math import factorial
import json
def involutions(n):
 return sum(factorial(n)//(2**j*factorial(j)*factorial(n-2*j)) for j in range(n//2+1))
def formula(p):
 N=(factorial(p-1)+(p-1)**2)//p; F=(involutions(p)+p-1)//p
 return dict(ordinary=N,self_converse=F,merged=(N+F)//2,involutions=involutions(p))
rows=[]
for p in [3,5,7]:
 k=(p-1)//2
 def orbit(q):return {tuple((q[(i-a)%p]+b)%p for i in range(p)) for a in range(p) for b in range(p)}
 seen=set();rs=[]
 def necklace(q):
  d=tuple((q[(i+1)%p]-q[i])%p for i in range(p))
  return min(d[i:]+d[:i] for i in range(p))
 for q in permutations(range(p)):
  if q in seen:continue
  o=orbit(q);seen.update(o);q=min(o);inv=[q.index(i) for i in range(p)]
  cq=min(orbit(tuple(-inv[-i%p]%p for i in range(p))))
  J=sum(1 if (q[j]-q[i])%p in range(1,k+1) else-1 for i in range(p) for j in range(p) if(j-i)%p in range(1,k+1))
  K=sum(1 if(q[(i+1)%p]-q[i])%p in range(1,k+1) else-1 for i in range(p))
  d=tuple((q[(i+1)%p]-q[i])%p for i in range(p))
  partial=[0]
  for x in d[:-1]:partial.append((partial[-1]+x)%p)
  if sum(d)%p or len(set(partial))!=p or necklace(tuple(partial))!=necklace(q):raise ValueError('necklace reconstruction')
  if any(necklace(x)!=necklace(q) for x in o):raise ValueError('necklace gauge')
  rs.append(dict(permutation=q,orbit_size=len(o),converse=cq,J=J,K=K))
 N=(factorial(p-1)+(p-1)**2)//p
 if len({necklace(tuple(r['permutation'])) for r in rs})!=N:raise ValueError('necklace completeness')
 if len(rs)!=N or len(seen)!=factorial(p):raise ValueError('Burnside')
 merged={min(tuple(r['permutation']),tuple(r['converse'])) for r in rs}
 f=formula(p)
 if len(merged)!=f['merged'] or sum(tuple(r['permutation'])==tuple(r['converse']) for r in rs)!=f['self_converse']:raise ValueError('converse Burnside')
 buckets={}
 for r in rs:buckets.setdefault(r['J'],set()).add(min(tuple(r['permutation']),tuple(r['converse'])))
 if p==7 and max(map(len,buckets.values()))<2:raise ValueError('missing scalar hostile')
 rows.append(dict(p=p,ordinary=N,merged=len(merged),self_converse=sum(tuple(r['permutation'])==tuple(r['converse']) for r in rs),exact_J_values=len(buckets),exact_JK_values=len({(r['J'],r['K']) for r in rs}),mod8_J_mod4_K_values=len({(r['J']%8,r['K']%4) for r in rs}),classes=rs))
# Exhaust every candidate displacement word for the small prime controls.
word_controls=[]
for p in [3,5,7]:
 valid=0
 for d in product(range(1,p),repeat=p-1):
  last=-sum(d)%p
  if last==0:continue
  partial=[0]
  for x in d:partial.append((partial[-1]+x)%p)
  if len(set(partial))==p:valid+=1
 if valid!=factorial(p-1):raise ValueError('word realizability count')
 word_controls.append(dict(p=p,normalized_valid_words=valid))
print(json.dumps(dict(status='PASS',word_controls=word_controls,exhaustive=rows,formula_only={str(p):formula(p) for p in [11,13,17]}),indent=2))
