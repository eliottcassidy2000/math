from itertools import combinations
from collections import Counter
n=10
E=[(i,j) for i,j in combinations(range(n),2) if j>i+1];M=(1<<len(E))-1
def allm(t):
 A=[0]*n
 for i in range(n-1):A[i]|=1<<(i+1)
 for k,(i,j) in enumerate(E):
  if t>>k&1:A[j]|=1<<i
  else:A[i]|=1<<j
 out={}
 def run(q,S):
  if not S:
   f=sum((not(A[q[i]]>>q[j]&1))<<k for k,(i,j) in enumerate(E))
   out.setdefault(f,[]).append(q);return
  Z=A[q[-1]]&S if q else S
  while Z:
   a=Z&-Z;Z^=a;run(q+(a.bit_length()-1,),S^a)
 run((),(1<<n)-1)
 return A,out
refl=[E.index((n-1-j,n-1-i)) for i,j in E]
sym=lambda t:all(not ((t>>i)^(t>>j))&1 for i,j in enumerate(refl))
x,y=4304406596,4311892004
A,L=allm(x); B,R=allm(y)
assert len(L)==len(R)==5213 and y not in L and not any(map(sym,L)) and not any(map(sym,R))
assert not any((t^M) in R for t in L)
print('scores',sorted(a.bit_count() for a in A),'H',sum(map(len,L.values())),'classsize',len(L),'opposite',y in L)
print('symmetric reroots',[(t,L[t]) for t in L if sym(t)])
print('first complement H,blue',sum(map(len,allm(x^M)[1].values())),sym(x))
print('adjacent opposite?',[(t,t^M) for t in L if (t^M) in R])
from itertools import permutations,product
def rows(t):
 A=[0]*n
 for i in range(n-1):A[i]|=1<<(i+1)
 for k,(i,j) in enumerate(E):
  if t>>k&1:A[j]|=1<<i
  else:A[i]|=1<<j
 return A
def canon(A):
 K=[a.bit_count() for a in A]
 while True:
  sizes={key:sum(1<<i for i,k in enumerate(K) if k==key) for key in sorted(set(K))}
  new=[(K[i],)+tuple((A[i]&mask).bit_count() for mask in sizes.values()) for i in range(n)]
  codes={k:i for i,k in enumerate(sorted(set(new)))}
  K=[codes[c] for c in new]
  if len(codes)==len(sizes):break
 parts=[[i for i in range(n) if K[i]==k] for k in sorted(set(K))]
 return min(sum((A[p[i]]>>p[j]&1)<<b for b,(i,j) in enumerate(combinations(range(n),2))) for oo in product(*[permutations(p) for p in parts]) for p in [sum(oo,())])
def book(h):
 D={}
 for t in h:D.setdefault(canon(rows(t^M)),t)
 return D
P,Q=book(L),book(R)
common=P.keys()&Q.keys()
assert len(common)==220
print('neighbor counts',len(P),len(Q),'2-path options',len(common))
if common:
 k=min(common); a,b=P[k],Q[k]
 print('two moves',a,a^M,b^M,b)
 print('reroot first from x',L[a][0])
 print('middle reroot', allm(a^M)[1][b^M][0])
 print('right reroot',allm(b)[1][y][0])

def literal(t):
 A=rows(t)
 values=[]
 for z in ('5613479280','1376042859','8597436012','0124367598'):
  p=tuple(map(int,z))
  assert all(A[u]>>v&1 for u,v in zip(p,p[1:]))
  q=sum((not(A[p[i]]>>p[j]&1))<<k for k,(i,j) in enumerate(E))
  for i,j in E:
   a,b=p[i],p[j];A[a]^=1<<b;A[b]^=1<<a
  r=sum((not(A[p[i]]>>p[j]&1))<<k for k,(i,j) in enumerate(E))
  assert q^r==M;values.append((z,q,r))
 assert sorted(x.bit_count() for x in A)==list(range(n))
 print('literal 4-edge witness from direct chord table',values)
literal(x)
