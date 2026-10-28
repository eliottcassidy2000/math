"""Separate verification formulas, by the same agent (not an independent reviewer)."""
from pathlib import Path
from math import comb, factorial, isqrt, prod
from fractions import Fraction as F
from itertools import combinations, permutations, product
from collections import Counter
import json
P=Path(__file__).resolve().parent
R=json.loads((P/'receipts.json').read_text());count=0

def check(v):
 global count
 count+=1
 if not v:raise RuntimeError('audit failure '+str(count))

def parity_integer(n):
 p=2;bits=set()
 while p*p<=n:
  parity=0
  while n%p==0:n//=p;parity^=1
  if parity:bits.add(p)
  p+=1
 if n>1:bits.add(n)
 return bits

for name,B in R['banks'].items():
 rows=B['rows'];N=len(rows)
 for row in rows:
  p=row['prime'];check(p>=2 and all(p%d for d in range(2,isqrt(p)+1)))
  check(p+1==row['squarefree_part']*row['multiplier']**2)
 enum=Counter()
 for w in range(1<<N):
  val=prod(rows[i]['prime']+1 for i in range(N) if w>>i&1)
  parity=set()
  for i in range(N):
   if w>>i&1:parity.symmetric_difference_update(parity_integer(rows[i]['prime']+1))
  check((isqrt(val)**2==val)==(not parity))
  if not parity:enum[w.bit_count()]+=1
 check(dict(enum)=={int(k):v for k,v in B['enumerator'].items()})

# All scalar moments via the distribution of h(c)-Z, not the multivariate expansion.
def uniform_sum_moments(n,d):
 mom=[F(1)]+[F(0)]*d
 for _ in range(n):mom=[sum(F(comb(k,j))*mom[j]/(k-j+1) for j in range(k+1)) for k in range(d+1)]
 return mom
momdata={}
for rank in [2,3,4]:
 n=2**rank-1;w=2**(rank-1);W=[Counter(),Counter()]
 for bits in range(1<<n):
  s=0
  for i in range(n):
   if bits>>i&1:s^=i+1
  if s in [0,1]:W[s][n-bits.bit_count()]+=1
 zm=uniform_sum_moments(n,w)
 moments=[]
 for q in [0,1]:
  size=sum(W[q].values())
  moments.append([sum(cnt*sum(F(comb(k,j))*(-1)**j*h**(k-j)*zm[j] for j in range(k+1)) for h,cnt in W[q].items())/size for k in range(w+1)])
 check(moments[0][:-1]==moments[1][:-1]);check(moments[0][-1]-moments[1][-1]==F(2**rank,2**w)*factorial(w))
 momdata[str(rank)]=[str(x) for x in moments[0]]+[str(moments[1][-1])]

# Reconstruct matching -> oriented triangle without reading the producer's matrix logic.
def all_matchings(S):
 if not S:yield ();return
 a=min(S)
 for b in sorted(S-{a}):
  for rest in all_matchings(S-{a,b}):yield tuple(sorted(((a,b),)+rest))
M4=[];inverses=[]
for M in all_matchings(set(range(8))):
 D={a^b for a,b in M}
 if len(D)!=4:continue
 K={a^b for a in D for b in D};f={}
 for a,b in M:f[a]=b;f[b]=a
 fixed=f[0];T={x:f[x]^fixed for x in K};T.update({x^fixed:T[x]^fixed for x in list(K)})
 check(all(T[x^y]==T[x]^T[y] for x in range(8) for y in range(8)))
 outside=set(range(8))-K-{fixed};a=min(outside);tri=(a,T[a],T[T[a]])
 check(set(tri)==outside and T[tri[-1]]==a)
 aa,bb,cc=tri
 expected=tuple(sorted(tuple(sorted(e)) for e in [(0,aa^bb^cc),(aa,aa^bb),(bb,bb^cc),(cc,cc^aa)]))
 check(expected==M);M4.append(M);inverses.append(tri)
check(len(M4)==len(set(inverses))==56)

# Source endpoint replay by a different affine numerator recurrence.
for row in R['collatz_cosets']:
 endpoints=[];descents=[]
 for n in row['sources']:
  A=1;B=0;word=0;values=[]
  for j in range(7):
   value=(A*n+B)//2**j
   if value%2:A*=3;B=3*B+2**j;word|=1<<j
   check((A*n+B)%2**(j+1)==0)
   values.append((A*n+B)//2**(j+1))
  sy=0
  for j in range(7):
   if word>>j&1:sy^=j+1
  check(sy==row['syndrome']);endpoints.append(values[-1]);descents.append(min(values)<n)
 check(sum(endpoints)==row['target_sum'])
 check(sum(t<n for t,n in zip(endpoints,row['sources']))==row['endpoint_below_source'])
 check(sum(descents)==row['descend_by7'])

# Inherited square-clock witnesses and complete finite lower-order obstruction gate.
shifts=json.loads((P/'shift_search.json').read_text())
for row in shifts:
 if row['path']:
  x,y=row['path'],row['other_path'];n=row['n'];c=row['shift']
  check(sorted(x)==sorted(y)==list(range(1,n+1)))
  check(all(isqrt(a+b)**2==a+b for a,b in zip(x,x[1:])))
  check([a+b for a,b in zip(x,x[1:])]==[a+b for a,b in zip(y,y[1:])])
  check(all(abs(a-b)==c for a,b in zip(x,y)))
small=[]
for n in range(2,36,2):
 for c in range(1,n//2+1):
  if n%(2*c):continue
  A={x for x in range(1,n+1) if (x-1)//c%2==0}
  adj={x:{y for y in range(1,n+1) if (x in A)!=(y in A) and isqrt(x+y)**2==x+y} for x in range(1,n+1)}
  seen={1};todo=[1]
  while todo:
   x=todo.pop()
   new=adj[x]-seen;seen|=new;todo+=list(new)
  reason='isolated' if any(not adj[x] for x in adj) else ('too_many_leaves' if sum(len(adj[x])==1 for x in adj)>2 else ('disconnected' if len(seen)<n else None))
  check(reason is not None);small.append((n,c,reason))
# Exact observability: equal widths lose all but cardinality, perturbed widths recover all.
heights=[];eps=F(1,512)
for bits in range(128):heights.append(sum(1+eps*2**i for i in range(7) if bits>>i&1))
check(len(set(heights))==128)
check(min(b-a for a,b in zip(sorted(heights),sorted(heights)[1:]))>=eps)
result={'status':'PASS; separate formulas by same agent, not independent-agent audit','checks':count,'scalar_moment_controls':momdata,'matching_inverses':56,'square_clock_lower_order_exclusions':len(small),'lower_order_exclusions':small,'generic_width_distinct_heights':128}
(P/'audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['lower_order_exclusions']},indent=2))
