from itertools import permutations,combinations,product
from fractions import Fraction as F
from math import comb
import random,json
r=random.Random(27801)
def paths(A):
 return [p for p in permutations(range(len(A))) if all(A[u][v] for u,v in zip(p,p[1:]))]
def contraction(A,E):
 n=len(A); nxt={}; prv={}
 for u,v in E:
  if u in nxt or v in prv:return 0
  nxt[u]=v;prv[v]=u
 blocks=[];seen=set()
 for u in range(n):
  if u in prv:continue
  seq=[u];seen.add(u)
  while seq[-1] in nxt:
   v=nxt[seq[-1]]
   if v in seen:return 0
   seq.append(v);seen.add(v)
  blocks.append(seq)
 if len(seen)<n:return 0
 C=[[int(i!=j and A[a[-1]][b[0]]) for j,b in enumerate(blocks)] for i,a in enumerate(blocks)]
 return len(paths(C))
cc=0
for n in range(1,7):
 for _ in range(12):
  A=[[int(i!=j and r.randrange(2)) for j in range(n)] for i in range(n)]
  ps=paths(A); arcs=[(i,j) for i in range(n) for j in range(n) if A[i][j]]
  for k in range(min(len(arcs),4)+1):
   for E in r.sample(list(combinations(arcs,k)),min(8,comb(len(arcs),k))):
    brute=sum(all(e in set(zip(p,p[1:])) for e in E) for p in ps)
    assert brute==contraction(A,E),(A,E,brute)
    cc+=1
# Independent Bernoulli derivation of the enumeration remainder, all rational q.
bc=0
for m in range(8):
 for _ in range(15):
  qs=[F(r.randrange(5),4) for i in range(m)]; e=[F(1)]+[F(0)]*m
  for q in qs:
   for j in range(m,0,-1):e[j]+=q*e[j-1]
  exact=F(1)
  for q in qs:exact*=1-q
  for k in range(m):
   partial=sum((-1)**j*e[j] for j in range(k+1)); rem=(-1)**k*(partial-exact)
   expectation=F(0)
   for bits in product((0,1),repeat=m):
    mass=F(1)
    for bit,q in zip(bits,qs):mass*=q if bit else 1-q
    K=sum(bits)
    if K>k:expectation+=mass*comb(K-1,k)
   assert rem==expectation
   assert F(k+1,m)*e[k+1]<=rem<=e[k+1]
   bc+=1
print(json.dumps({'forced_contraction_checks':cc,'rational_Bernoulli_remainder_checks':bc,'status':'PASS'},indent=2))
