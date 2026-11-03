from fractions import Fraction as F
from itertools import permutations, combinations
from math import factorial
from functools import lru_cache
# Tournament on 4 blocks; 0->1->2->0, 3->0, 1->3, 3->2.
B=((0,1,0,0),(0,0,1,1),(1,0,0,0),(1,0,1,0))
N=6; weights=(1,2,3,4)
def mul(a,b):return [sum(a[j]*b[i-j] for j in range(i+1)) for i in range(N+1)]
def det(A):
 out=[F(0)]*(N+1)
 for perm in permutations(range(4)):
  inversions=sum(perm[i]>perm[j] for i in range(4) for j in range(i+1,4))
  p=[F(1)]+[F(0)]*N
  for i in range(4):p=mul(p,A[i][perm[i]])
  out=[x+(-1)**inversions*y for x,y in zip(out,p)]
 return out
A=[]
for i in range(4):
 row=[]
 for j in range(4):
  row.append(([F(1)]+[F(0)]*N) if i==j else [F(0)]+[-B[i][j]*F(weights[i]**k,factorial(k)) for k in range(1,N+1)])
 A.append(row)
D=det(A); Dm=[(-1)**n*d for n,d in enumerate(D)]
R=[F(0)]*(N+1)
for n in range(N+1):R[n]=Dm[n]-sum(D[k]*R[n-k] for k in range(1,n+1))
exp_s=[F(sum(weights)**n,factorial(n)) for n in range(N+1)]
def predicted(t):
 alpha=F(t,2); power=[F(1)]+[F(0)]*N
 for n in range(1,N+1):
  power[n]=sum(((alpha+1)*k-n)*R[k]*power[n-k] for k in range(1,n+1))/n
 return mul(exp_s,power)
def direct(ns,t):
 n=sum(ns); blocks=sum(([i]*ns[i] for i in range(4)),[])
 def edge(i,j):return i<j if blocks[i]==blocks[j] else bool(B[blocks[i]][blocks[j]])
 cycles=[]
 for size in range(3,n+1,2):
  for S in combinations(range(n),size):
   for p in permutations(S[1:]):
    v=(S[0],)+p
    if all(edge(v[i],v[(i+1)%size]) for i in range(size)):
     cycles.append(sum(1<<k for k in S))
 @lru_cache(None)
 def pack(mask):
  if not mask:return 1
  bit=mask&-mask; result=pack(mask^bit)
  for cyc in cycles:
   if cyc&bit and cyc&mask==cyc:result+=t*pack(mask^cyc)
  return result
 return pack((1<<n)-1)
def compositions(n,d):
 if d==1:
  yield (n,);return
 for i in range(n+1):
  for rest in compositions(n-i,d-1):yield (i,)+rest
for t in (0,1,2,3,4):
 got=[F(0)]*(N+1)
 for n in range(N+1):
  for ns in compositions(n,4):
   weight=F(1)
   for w,k in zip(weights,ns):weight*=F(w**k,factorial(k))
   got[n]+=weight*direct(ns,t)
 if got!=predicted(t):raise ValueError((t,got,predicted(t)))
print('Verified 4-block determinant gas EGF through total order 6, with t=0,1,2,3,4.')
