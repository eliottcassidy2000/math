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
K=[I];seen={I:0}
for g in K:
 for h,tag in zip(GG,[0]*6+[1,2]):
  z=tuple(g[i] for i in h)
  new=seen[g]^tag
  if z not in seen:K.append(z);seen[z]=new
  else:assert seen[z]==new
assert len(K)==80640


def span3(bank,title):
 B={};N=56*56
 def add(P,Q):
  a,b=P;c,d=Q;u=~(a|b);v=~(c|d)
  return ((a&v)|(u&c)|(b&d),(b&v)|(u&d)|(a&c))
 for number,g in enumerate(bank,1):
  a=sum(1<<(56*i+x+1) for i,x in enumerate(g))|(seen[g]!=0);b=0
  while (a|b)>=2:
   k=(a|b).bit_length()-1;amt=(a>>k)&1
   if k not in B:
    B[k]=(a,b) if amt else (b,a);break
   c,d=B[k]
   a,b=add((a,b),(d,c) if amt else (c,d))
  else:
   if a|b:
    print(title,'tag fails at row',number,'rank',len(B),flush=True);return
 a=((1<<N)-1)<<1;b=0
 while (a|b)>=2:
  k=(a|b).bit_length()-1;amt=(a>>k)&1
  if k not in B:break
  c,d=B[k];a,b=add((a,b),(d,c) if amt else (c,d))
 print(title,'rank',len(B),'target reduction',a,b,flush=True)
span3([g for g in K if seen[g] in (0,3)],'all even')

