"""Independent recursive implementation; no producer imports; explicit gates."""
import random,json
from fractions import Fraction
rng=random.Random(116931)
gates=0
def eq(a,b):
 global gates
 gates+=1
 if a!=b: raise RuntimeError((gates,a,b))
def add(a,b):return tuple(x+y for x,y in zip(a,b))
def neg(a):return tuple(-x for x in a)
def conj(a):return (a[0],)+tuple(-x for x in a[1:])
def mul(x,y):
 if len(x)==1:return (x[0]*y[0],)
 h=len(x)//2;a,b=x[:h],x[h:];c,d=y[:h],y[h:]
 return add(mul(a,c),neg(mul(conj(d),b)))+add(mul(d,a),mul(b,conj(c)))
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def norm(a):return dot(a,a)
def basis(n,j):return tuple(int(i==j) for i in range(n))
def mm(a,b):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def det(a):
 a=[list(map(Fraction,r)) for r in a];v=Fraction(1)
 for k in range(len(a)):
  p=next((i for i in range(k,len(a)) if a[i][k]),None)
  if p is None:return 0
  if p!=k:a[k],a[p]=a[p],a[k];v=-v
  z=a[k][k];v*=z
  for i in range(k+1,len(a)):
   t=a[i][k]/z
   for j in range(k+1,len(a)):a[i][j]-=t*a[k][j]
 return v

def rankmod(a,p):
 a=[[x%p for x in r] for r in a];r=0
 for c in range(len(a[0])):
  k=next((i for i in range(r,len(a)) if a[i][c]),None)
  if k is None:continue
  a[r],a[k]=a[k],a[r];v=pow(a[r][c],-1,p);a[r]=[x*v%p for x in a[r]]
  for i in range(r+1,len(a)):
   v=a[i][c];a[i]=[(x-v*y)%p for x,y in zip(a[i],a[r])]
  r+=1
 return r
bank=[(0,)*16,basis(16,0),add(basis(16,0),basis(16,9))]
# Exhaustive signed orthogonal basis zero-divisor cone.
for i in range(1,8):
 for j in range(1,8):
  if i!=j:bank.append(basis(8,i)+basis(8,j))
# Collinear/commuting imaginary parts, arbitrary real parts.
for t in range(-2,3):bank.append((t,1,2,0,0,0,0,0)+(1,2,4,0,0,0,0,0))
bank += [tuple(rng.randrange(-2,3) for _ in range(16)) for z in range(25)]
for x in bank:
 a,b=x[:8],x[8:];n=norm(x);delta=norm(a[1:])*norm(b[1:])-dot(a[1:],b[1:])**2
 cols=[mul(x,basis(16,j)) for j in range(16)];L=list(map(list,zip(*cols)));G=mm(list(map(list,zip(*L))),L)
 D=[[G[i][j]-n*(i==j) for j in range(16)] for i in range(16)];D2=mm(D,D);D3=mm(D2,D)
 eq(D3,[[4*delta*v for v in row] for row in D]);eq(sum(D[i][i] for i in range(16)),0);eq(sum(D2[i][i] for i in range(16)),32*delta)
 eq(det(L),n**4*(n*n-4*delta)**2)
 for p in [2,3,5,7,11,13]:eq(rankmod(L,p)<16,(n*(n*n-4*delta))%p==0)
 if n and n*n==4*delta:eq(rankmod(L,1000003),12)
# Associator norm defect, independently sampled arbitrary four octonions.
for trial in range(200):
 a,b,c,d=[tuple(rng.randrange(-2,3) for i in range(8)) for _ in range(4)]
 assoc=add(mul(mul(d,a),c),neg(mul(d,mul(a,c))))
 eq(norm(mul(a+b,c+d))-norm(a+b)*norm(c+d),2*dot(assoc,b))
# All pure basis singular pairs and all H-perpendicular basis c: kernel formula.
for i in range(1,8):
 for j in range(1,8):
  if i==j:continue
  a=basis(8,i);b=basis(8,j);ab=mul(a,b)
  for k in range(8):
   c=basis(8,k)
   if c[0] or dot(c,a) or dot(c,b) or dot(c,ab):continue
   d=mul(b,mul(c,a));eq(mul(a+b,c+d),(0,)*16)
# Every stage: quadratic relation and inverse conjugation norm.
for n in [1,2,4,8,16,32,64]:
 for z in range(5):
  x=tuple(rng.randrange(-2,3) for _ in range(n));unit=basis(n,0)
  eq(mul(x,conj(x)),tuple(norm(x)*u for u in unit));eq(mul(conj(x),x),tuple(norm(x)*u for u in unit))
  eq(mul(x,x),tuple(2*x[0]*v-norm(x)*u for v,u in zip(x,unit)))
# Exhaustive independent W zero-product relation, including degenerate factors.
from itertools import product
qb=[(0,)*4]+[tuple(s*v for v in basis(4,i)) for i in range(4) for s in [-1,1]]
for a,s,c,t in product(qb,repeat=4):
 x=a+(0,)*4+(0,)*4+s;y=c+(0,)*4+(0,)*4+t
 reduced=add(mul(a,c),neg(mul(conj(s),t)))+(0,)*8+add(mul(t,conj(a)),mul(s,c))
 eq(mul(x,y),reduced)
 if norm(x) and norm(y):
  cond=norm(a)>0 and norm(c)>0 and a[0]==0 and c[0]==0 and dot(a,c)==0 and norm(s)==norm(a) and tuple(norm(a)*v for v in t)==mul(mul(s,a),c)
  eq(norm(mul(x,y))==0,cond)
print(json.dumps({'status':'PASS','gates':gates,'spectrum_cases':len(bank),'seed':116931},sort_keys=True))
