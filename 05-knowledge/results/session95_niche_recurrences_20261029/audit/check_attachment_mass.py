"""Independent permutation-action checks of recursive orbit mass moments."""
from collections import defaultdict
from fractions import Fraction as Q
from itertools import product
from functools import lru_cache
import json
checks=0

def check(a,b,tag):
 global checks
 checks+=1
 if a!=b:raise RuntimeError((tag,a,b))

def perms2():
 for r in range(3):
  for shifts in product(range(3),repeat=3):
   yield tuple(3*((b+r)%3)+(i+shifts[b])%3 for b in range(3) for i in range(3))
G=list(perms2());alpha=tuple(3*((-b)%3)+(-i)%3 for b in range(3) for i in range(3))
C=[g for g in G if all(g[alpha[i]]==alpha[g[i]] for i in range(9))]
check(len(G),81,'group');check(len(C),3,'centralizer')
def orbit_mass(group,p):
 out=defaultdict(Q)
 for f in product(range(2),repeat=9):
  canonical=min(tuple(f[g[i]] for i in range(9)) for g in group)
  weight=p**sum(f)*(1-p)**(9-sum(f))
  out[canonical]+=weight
 return out
# Use Fraction power for r0, avoiding float in3^(-1).
@lru_cache(None)
def exact_recurse(k,r,palette):
 if k==0:return sum((v**r for v in palette),Q(0))
 f=Q(3)**(r-1)
 return f*exact_recurse(k-1,r,palette)**3+(1-f)*exact_recurse(k-1,3*r,palette)
for p in [Q(1,2),Q(1,3),Q(1,5)]:
 A=orbit_mass(G,p);D=orbit_mass(C,p);palette=(p,1-p);paired=tuple(x*y for x in palette for y in palette)
 for r in range(5):
  check(sum(v**r for v in A.values()),exact_recurse(2,r,palette),('A',p,r))
  pred=sum(v**r for v in palette)*exact_recurse(0,r,paired)*exact_recurse(1,r,paired)
  check(sum(v**r for v in D.values()),pred,('C',p,r))
 # Global converse quotient exact source mass, not halfmass per orbit.
 folded=defaultdict(Q)
 for f,w in D.items():
  rev=tuple(1-f[alpha[i]] for i in range(9))
  crev=min(tuple(rev[g[i]] for i in range(9)) for g in C)
  folded[min(f,crev)]+=w
 check(sum(folded.values()),1,('mass',p))
 check(len(folded)*2,len(D),('free converse',p))
 if p!=Q(1,2):
  check(any(D[f]!=D[min(tuple(1-f[alpha[g[i]]] for i in range(9)) for g in C)] for f in D),True,('unequal partner masses',p))
# Heterogeneous source collision needs the six-way overlap, not merely pair data.
def cyclic_key(t):return min(t,t[1:]+t[:1],t[2:]+t[:2])
def pushed(ds):
 out=defaultdict(Q)
 for t in product(*(tuple(d) for d in ds)):
  out[cyclic_key(t)]+=ds[0][t[0]]*ds[1][t[1]]*ds[2][t[2]]
 return out
answers=[]
for parity in (0,1):
 sets=[{x:Q(1,16) for x in range(32) if not(x>>j&1)} for j in range(5)]
 sets.append({x:Q(1,16) for x in range(32) if x.bit_count()%2==parity})
 A,B=sets[:3],sets[3:]
 overlaps=[[sum((a.get(x,0)*b.get(x,0) for x in range(32)),Q(0)) for b in B] for a in A]
 check(overlaps,[[Q(1,32)]*3 for _ in range(3)],('same marginal overlaps',parity))
 six=sum((__import__('math').prod(d.get(x,0) for d in sets) for x in range(32)),Q(0))
 pa,pb=pushed(A),pushed(B)
 collision=sum((v*pb.get(k,0) for k,v in pa.items()),Q(0))
 rhs=sum((__import__('math').prod(overlaps[i][(i+j)%3] for i in range(3)) for j in range(3)),Q(0))-2*six
 check(collision,rhs,('exact orbit collision',parity));answers.append(collision)
check(answers,[Q(1534,16**6),Q(1536,16**6)],'higher-overlap hostile')
print(json.dumps({'status':'PASS','checks':checks,'group_size':len(G),'centralizer_size':len(C),'producer_imports':False},sort_keys=True))
