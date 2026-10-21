from itertools import combinations,permutations,product
from fractions import Fraction
from pathlib import Path
import json,hashlib
count=0
def need(v,label):
 global count
 count+=1
 if not v:raise RuntimeError(label)
R={1,2,4}
def edge(i,j):return (j-i)%7 in R
def transitive(S):
 return any(all(edge(p[i],p[j]) for i in range(len(p)) for j in range(i+1,len(p))) for p in permutations(S))
need(sum(transitive(s) for s in combinations(range(7),3))==21,'transitive triples')
need(not any(transitive(s) for s in combinations(range(7),4)),'no transitive four')
base={0:[],1:[(0,1,2),(3,4),(5,6)],2:[(0,1,2),(0,1,2),(3,4,5),(4,5,6),(3,6)],3:[tuple((j+t)%7 for j in (0,1,2)) for t in range(7)]}
for k in range(121):
 q,r=divmod(k,3);pal=base[3]*q+base[r]
 need(len(pal)==(7*k+2)//3,'exact palette size')
 need(all(sum(i in s for s in pal)==k for i in range(7)),'exact vertex multiplicities')
 need(all(transitive(s) for s in pal),'legal supports')
a=[(7*k+2)//3 for k in range(121)]
for k in range(4,121):need(a[k]==a[k-1]+a[k-3]-a[k-4],'order four')
# Exact denominator multiplication of the proposed GF.
for k in range(121):
 c=a[k]-(a[k-1] if k>=1 else 0)-(a[k-3] if k>=3 else 0)+(a[k-4] if k>=4 else 0)
 need(c==({1:3,2:2,3:2}.get(k,0)),'rational numerator')
# Fano line hypergraph: triples are cyclic in this P7, each pair occurs once.
lines=[tuple(sorted((j+t)%7 for j in (0,1,3))) for t in range(7)]
for i,j in combinations(range(7),2):need(sum(i in L and j in L for L in lines)==1,'Fano pair partition')
need(all(not transitive(L) for L in lines),'Fano lines directed cycles')
for mask in range(128):need(any(len({(mask>>v)&1 for v in L})==1 for L in lines),'not Property B')
# Structural substitution test at every vertex subset of C3[C3] (nine vertices).
blocks=3;size=3
edge3=lambda i,j:(j-i)%3==1
def e9(i,j):
 b,x=divmod(i,3);c,y=divmod(j,3)
 return edge3(x,y) if b==c else edge3(b,c)
def acyclic(V,e):
 V=set(V)
 while V:
  src=[v for v in V if all(v==w or e(v,w) for w in V)]
  if not src:return False
  V.remove(src[0])
 return True
for mask in range(512):
 V=[i for i in range(9) if mask>>i&1];occupied={i//3 for i in V}
 rhs=acyclic(occupied,edge3) and all(acyclic([i%3 for i in V if i//3==b],edge3) for b in occupied)
 need(acyclic(V,e9)==rhs,'substitution exact support criterion')
# Finite recurrence controls, never a finite proof of nonrationality.
seq=[1];carry=[]
for i in range(250):
 nxt=(7*seq[-1]+2)//3;carry.append(3*nxt-7*seq[-1]);seq.append(nxt)
 need(carry[-1] in (0,1,2),'carry range')
 need(Fraction(seq[-1]*3**(i+1),7**(i+1))<=Fraction(3,2),'normalized bound')
need(seq[:6]==[1,3,7,17,40,94],'depth prefix')
for A,B in [(3,2),(5,2),(7,3),(8,3),(11,7)]:
 x=1
 for n in range(100):
  y=(A*x+B-1)//B;e=B*y-A*x
  need(0<=e<B and y>x,'general rounded recurrence');x=y
root=Path(__file__).resolve().parent
print(json.dumps({'status':'PASS','checks':count,'scope':'exact finite palette, hypergraph, substitution and recurrence controls; nonrationality follows from the proof, not these finite prefixes','proof_sha256':hashlib.sha256((root/'PALEY_MULTICOLORING_AND_CARRY.md').read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'depth_prefix':seq[:15],'first_carries':carry[:30]},indent=2))
