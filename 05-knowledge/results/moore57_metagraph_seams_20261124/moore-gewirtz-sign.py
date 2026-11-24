from itertools import product, combinations
from collections import Counter
from math import lcm

def mul(x,y):
 z=0
 while y:
  if y&1:z^=x
  x<<=1
  if x&4:x^=7
  y>>=1
 return z

def norm(p):
 x=next(t for t in p if t); inv=next(i for i in range(1,4) if mul(i,x)==1)
 return tuple(mul(t,inv) for t in p)
V=sorted({norm(x) for x in product(range(4),repeat=3) if any(x)})
def dot(x,y):return mul(x[0],y[0])^mul(x[1],y[1])^mul(x[2],y[2])
L=[sum(1<<i for i,p in enumerate(V) if dot(p,l)==0) for l in V]
def actvec(p,i,j,c):
 p=list(p);p[i]^=mul(c,p[j]);return norm(p)
def pmask(p,h):return sum(1<<p[i] for i in range(21) if h>>i&1)
GG=[[V.index(actvec(p,i,j,c)) for p in V] for i,j in [(0,1),(1,2),(2,0)] for c in [1,2]]
start=next(h for p in combinations(range(21),6) if all((m&(h:=sum(1<<i for i in p))).bit_count() in [0,2] for m in L))
H=[start];hset={start}
for h in H:
 for g in GG:
  z=pmask(g,h)
  if z not in hset:hset.add(z);H.append(z)
assert len(H)==56
H.sort()
def dual(h):return sum(1<<i for i,m in enumerate(L) if not m&h)
frob=[V.index(tuple(mul(x,x) for x in p)) for p in V]
GG=[tuple(H.index(pmask(g,h)) for h in H) for g in GG]+[tuple(H.index(pmask(frob,h)) for h in H),tuple(H.index(dual(h)) for h in H)]
A=[[not (j&k) for j in H] for k in H]
assert all(sum(x)==10 for x in A)
assert all(A[g[i]][g[j]]==A[i][j] for g in GG for i in range(56) for j in range(56))
I=tuple(range(56))
K=[I];seen={I}
for g in K:
 for h in GG:
  z=tuple(g[i] for i in h)
  if z not in seen:K.append(z);seen.add(z)
assert len(K)==80640
def ds(g):
 seen=set();parts=[]
 for i in I:
  if i in seen:continue
  x=i;k=0
  while x not in seen:
   seen.add(x);k+=1;x=g[x]
  parts.append(k)
 m=lcm(*parts)
 s=(-1)**(56-len(parts))
 if m%2:return s,1,1
 x=I
 alt=[0,0]
 for t in range(m):
  f=sum(i==v for i,v in enumerate(x))
  a=sum(A[i][v] for i,v in enumerate(x))
  # projector traces multiplied by 6
  alt[0]+=(-1)**t*(a+4*f-14)
  alt[1]+=(-1)**t*(2*f-a+8)
  x=tuple(g[i] for i in x)
 assert all(v%(6*m)==0 for v in alt)
 d2,d4=((-1)**(v//(6*m)) for v in alt)
 assert d2*d4==s
 return s,d2,d4
print("generator signs", [ds(g) for g in GG],flush=True)
# The identity and the six elementary transformations generate PSL3(4), so
# the two projective point/line transformations represent the outer cosets.
print("positional derangement signs",Counter(ds(g) for g in K if all(i!=x for i,x in enumerate(g))),flush=True)
B = (25,27,33,35,37,39,45,47,49,50,53,54)
C = frozenset((9,11,14,15,25,26,29,31,32,34,37,39,41,43,44,47,48,51,53,54))
def parity(g):
 unseen=set(I); cyc=0
 while unseen:
  x=min(unseen);cyc+=1
  while x in unseen:
   unseen.remove(x);x=g[x]
 return cyc%2
assert len(B)*len(C)==240
images={0:set(),1:set()}
hist={0:Counter(),1:Counter()}
for g in K:
 p=parity(g)
 img=frozenset(g[x] for x in B)
 q=len(img&C)
 assert q%2==p
 images[p].add(img)
 hist[p][q]+=1
assert not(images[0]&images[1])
print("sign rectangle:",len(B),len(C), tuple(len(images[p]) for p in range(2)), hist,flush=True)
print("internal degrees",sorted(sum(A[x][y] for y in B) for x in B),sorted(sum(A[x][y] for y in C) for x in C),flush=True)
