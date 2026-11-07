from itertools import product
from math import gcd
from fractions import Fraction as F
from collections import defaultdict,Counter

def require(v,msg):
 if not v:raise RuntimeError(msg)
def move(s,c,q):
 D,E,U,d=s
 vals=(D,0,2*U,d) if c==0 else (D-E+U,d,2*U-d,d) if c==1 else (D-E+U,D-E,D+U,D)
 return tuple(x%q for x in vals)
def key(s,q):
 D,E,U,d=s;g=gcd(gcd(U,d),q)
 return d,g,D%g

def audit(q):
 groups=defaultdict(list)
 for d,D,U in product(range(q),repeat=3):
  for E in {0,d}:
   s=(D,E,U,d);groups[key(s,q)].append(s)
 for k,states in groups.items():
  st=set(states);rev=defaultdict(list);incoming=Counter()
  for s in states:
   for c in [0,1]:
    v=move(s,c,q);require(v in st,('key invariant',s,c));rev[v].append(s);incoming[v]+=1
  require(all(incoming[s]==2 for s in states),'regularity')
  for backwards in [False,True]:
   seen={states[0]};todo=[states[0]]
   for s in todo:
    for v in rev[s] if backwards else [move(s,c,q) for c in [0,1]]:
     if v not in seen:seen.add(v);todo.append(v)
   require(seen==st,('class connectivity',q,k))
 # Push arbitrary distributions by explicit first letter then normalized gcd class.
 def project(dist):
  weights=defaultdict(F)
  for s,p in dist.items():
   for c in [0,1]:weights[key(move(s,c,q),q)]+=p/2
  out={}
  for k,p in weights.items():
   for s in groups[k]:out[s]=p/len(groups[k])
  return out
 dist=project({(1,0,2,1):F(1)})
 for j in range(3):
  probs=[F(0)]*q
  for (D,E,U,d),p in dist.items():probs[(D-E+U)%q]+=p
  print('q',q,'J',j,'probabilities',','.join(map(str,probs)))
  require(sum(probs)==1,'mass')
  nxt=defaultdict(F)
  for s,p in dist.items():nxt[move(s,2,q)]+=p
  dist=project(nxt)
 print('PASS gcd classes',q,len(groups),sum(map(len,groups.values())))
for q in [3,9,27]:audit(q)
