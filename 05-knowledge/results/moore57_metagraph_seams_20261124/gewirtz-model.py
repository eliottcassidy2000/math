from itertools import product,combinations
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
G=[[V.index(actvec(p,i,j,c)) for p in V] for i,j in [(0,1),(1,2),(2,0)] for c in [1,2]]
start=next(h for p in combinations(range(21),6) if all((m&(h:=sum(1<<i for i in p))).bit_count() in [0,2] for m in L))
H=[start];hset={start}
for h in H:
 for g in G:
  z=pmask(g,h)
  if z not in hset:hset.add(z);H.append(z)
assert len(H)==56
H.sort()
def dual(h):return sum(1<<i for i,m in enumerate(L) if not m&h)
frob=[V.index(tuple(mul(x,x) for x in p)) for p in V]
G=[tuple(H.index(pmask(g,h)) for h in H) for g in G]+[tuple(H.index(pmask(frob,h)) for h in H),tuple(H.index(dual(h)) for h in H)]
def group(G):
 S=[tuple(range(56))];seen=set(S)
 for g in S:
  for h in G:
   z=tuple(g[i] for i in h)
   if z not in seen:S.append(z);seen.add(z)
 return S
n=56
def check(cond,msg):
 if not cond:raise RuntimeError(msg)
def inverse(g):
 return tuple(sorted(range(n),key=g.__getitem__))
def derange(g):
 return all(i!=j for i,j in enumerate(g))

