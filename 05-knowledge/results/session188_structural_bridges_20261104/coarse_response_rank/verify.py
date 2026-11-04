"""Exact universal ten-dimensional block-uniform cyclic response contrast."""
from fractions import Fraction as F
from itertools import product
from math import comb,factorial
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parent
KEYS=list(product(range(8),repeat=2))
BASIS=[(0,1),(0,2),(1,0),(1,1),(1,2),(1,3),(1,4),(2,0),(2,1),(2,2)]
WITNESSES=[(1,1,1),(1,1,2),(1,2,1),(2,1,1),(1,1,3),(1,2,2),(1,3,1),(2,1,2),(1,2,3),(2,1,3)]
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def add(*ps):
 out={}
 for p in ps:
  for e,c in p.items():out[e]=out.get(e,0)+c
 return {e:c for e,c in out.items() if c}
def scale(p,c):return {e:c*v for e,v in p.items() if c*v}
def mul(*ps):
 out={(0,0,0):1}
 for p in ps:
  nxt={}
  for e,c in out.items():
   for f,d in p.items():
    g=tuple(x+y for x,y in zip(e,f));nxt[g]=nxt.get(g,0)+c*d
  out={e:c for e,c in nxt.items() if c}
 return out
ONE={(0,0,0):1};X={(1,0,0):1};Y={(0,1,0):1};Z={(0,0,1):1};P=mul(X,Y,Z);R=add(ONE,scale(P,-1))
G=[[X,mul(X,Y),P],[P,Y,mul(Y,Z)],[mul(X,Z),P,Z]]
def numerator(G,a,c):
 A=[i for i in range(3) if a>>i&1];O=[i for i in range(3) if not(a>>i&1)];C=[i for i in range(3) if c>>i&1];B=[i for i in range(3) if not(c>>i&1)]
 def left(I):return add(R,*[G[u][w] for u in range(3) for w in I])
 def right(I):return add(R,*[G[u][w] for u in I for w in range(3)])
 def middle(I,J):return add(*[G[u][w] for u in I for w in J])
 return add(mul(left(A),add(R,middle(O,C)),right(B)),mul(left(C),middle(B,A),right(O)))
def coordinates(polys,basis):
 monomials=sorted(set().union(*(set(p) for p in polys)))
 A=[[F(p.get(m,0)) for p in basis]+[F(p.get(m,0)) for p in polys] for m in monomials]
 row=0;n=len(basis)
 for col in range(n):
  k=next((i for i in range(row,len(A)) if A[i][col]),None);need(k is not None,'rank deficient claimed basis');A[row],A[k]=A[k],A[row]
  d=A[row][col];A[row]=[x/d for x in A[row]]
  for i in range(len(A)):
   if i!=row and A[i][col]:
    d=A[i][col];A[i]=[x-d*y for x,y in zip(A[i],A[row])]
  row+=1
 need(all(not any(r) for r in A[n:]),'a numerator is outside the basis span')
 return [[A[i][n+j] for i in range(n)] for j in range(len(polys))],monomials
def determinant(A):
 A=[[F(x) for x in row] for row in A];z=F(1);n=len(A)
 for j in range(n):
  i=next((i for i in range(j,n) if A[i][j]),None)
  if i is None:return 0
  if i!=j:A[i],A[j]=A[j],A[i];z=-z
  d=A[j][j];z*=d
  for i in range(j+1,n):
   t=A[i][j]/d
   for k in range(j,n):A[i][k]-=t*A[j][k]
 need(z.denominator==1,'integer determinant');return z.numerator
def functional(p,qs):
 total=0
 for k in range(max(map(len,qs))):
  for e,c in p.items():
   term=c*comb(k+2,2)
   for q,i in zip(qs,e):term*=q[k+i] if k+i<len(q) else 0
   total+=term
 return total
def tt(n):return [((1<<n)-1)^((1<<(u+1))-1) for u in range(n)]
def tt_profile(n):
 q=[0,1]
 for m in range(2,n+1):q=[0]+[k*((q[k] if k<len(q) else 0)+q[k-1]) for k in range(1,m+1)]
 return q
def build(blocks,reverse=False):
 offsets=[0,len(blocks[0]),len(blocks[0])+len(blocks[1])];out=[0]*sum(map(len,blocks))
 for i,B in enumerate(blocks):
  j=(i+(2 if reverse else 1))%3;target=((1<<len(blocks[j]))-1)<<offsets[j]
  for u,r in enumerate(B):out[offsets[i]+u]=(r<<offsets[i])|target
 return out,offsets
def physical_mask(m,blocks,offsets):return sum(((1<<len(B))-1)<<o for i,(B,o) in enumerate(zip(blocks,offsets)) if m>>i&1)
def extend(T,a,c):
 n=len(T);out=list(T)+[1<<(n+1),0]
 for u in range(n):
  if a>>u&1:out[u]|=1<<n
  else:out[n]|=1<<u
  if c>>u&1:out[u]|=1<<(n+1)
  else:out[n+1]|=1<<u
 return out
def hp(T):
 n=len(T);d=[[0]*n for _ in range(1<<n)]
 for S in range(1,1<<n):
  rem=S
  while rem:
   bit=rem&-rem;rem-=bit;u=bit.bit_length()-1;rest=S^bit
   if not rest:d[S][u]=1;continue
   destinations=T[u]&rest;z=0
   while destinations:
    b=destinations&-destinations;destinations-=b;z+=d[rest][b.bit_length()-1]
   d[S][u]=z
 return sum(d[-1])
def main():
 Gt=list(map(list,zip(*G)));nums=[numerator(G,*key) for key in KEYS];opnums=[numerator(Gt,*key) for key in KEYS];D=[add(a,scale(b,-1)) for a,b in zip(nums,opnums)]
 basis=[D[KEYS.index(key)] for key in BASIS];coords,mons=coordinates(D,basis);need(len(mons)==48,'complete numerator support');twice=[]
 for d,cs in zip(D,coords):
  need(all((2*x).denominator==1 for x in cs),'coefficient denominator exceeds two');ints=[int(2*x) for x in cs];need(max(map(abs,ints))<=4,'coefficient magnitude bound');twice.append(ints)
  need(scale(d,2)==add(*[scale(p,c) for p,c in zip(basis,ints)]),'coefficientwise doubled numerator identity')
 literal=[];rank_rows=[]
 cases=[([tt(n) for n in ns],[tt_profile(n) for n in ns],'rank_witness',ns) for ns in WITNESSES]
 cases += [([[2,4,1]]*3,[[0,3,6,6]]*3,'pair_free_equal',(3,3,3)),([tt(3),tt(3),[2,4,1]],[[0,1,6,6],[0,1,6,6],[0,3,6,6]],'one_probe_cancellation',(3,3,3))]
 for blocks,qs,kind,sizes in cases:
  T,offs=build(blocks);U,_=build(blocks,True);values=[];gaps=[]
  for i,(a,c) in enumerate(KEYS):
   aa=physical_mask(a,blocks,offs);cc=physical_mask(c,blocks,offs);f=hp(extend(T,aa,cc));g=hp(extend(U,aa,cc));values.append([f,g]);gaps.append(f-g)
   need(f==functional(nums[i],qs) and g==functional(opnums[i],qs),'literal count versus full physical-word functional')
  observed=[gaps[KEYS.index(k)] for k in BASIS]
  for gap,co in zip(gaps,twice):need(2*gap==sum(a*b for a,b in zip(co,observed)),'literal reconstruction')
  if kind=='rank_witness':rank_rows.append(observed)
  literal.append({'kind':kind,'child_sizes':sizes,'child_rows':blocks,'child_profiles':qs,'basis_contrasts':observed,'full64_response_pairs':values})
 det=determinant(rank_rows);need(det==-3057647616,'literal witness determinant')
 return {'status':'PASS','scope':'All-order rank-ten compression of all64 block-uniform cyclic quotient-converse contrasts; not universal separation or full named reconstruction.','basis_block_queries':BASIS,'twice_reconstruction_coefficients_in_lexicographic_query_order':twice,'kernel_monomials':mons,'kernel_rank':10,'physical_witness_determinant':det,'literal_physical_counts':len(cases)*128,'full64_literal_controls':literal,'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'PROOF.md',ROOT/'verify.py']}}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
