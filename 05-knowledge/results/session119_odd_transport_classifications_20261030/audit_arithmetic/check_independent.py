"""Independent literal determinant, modular image, and cyclic-trajectory audit."""
from itertools import product,permutations,combinations
from math import gcd,prod
from fractions import Fraction
from collections import Counter
import json

def check(ok,message):
 if not ok:raise RuntimeError(message)
def determinant(M):
 n=len(M)
 if n==0:return 1
 total=0
 for q in permutations(range(n)):
  inv=sum(q[i]>q[j] for i in range(n) for j in range(i+1,n))
  total+=(-1)**inv*prod(M[i][q[i]] for i in range(n))
 return total
def matrix(p,u,a):
 k=len(a);M=[[0]*k for _ in a]
 for i,v in enumerate(a):M[i][i]-=u;M[i][(i+1)%k]+=p**v
 return M
def divisors(M):
 k=len(M);answer=[]
 for r in range(1,k+1):
  d=0
  for rows in combinations(range(k),r):
   for cols in combinations(range(k),r):d=gcd(d,determinant([[M[i][j] for j in cols] for i in rows]))
  answer.append(d)
 return answer
def numerator(p,u,a,b,i=0):
 k=len(a);v=0;prefix=1
 for j in range(k):
  v+=u**(k-1-j)*prefix*b[(i+j)%k]
  prefix*=p**a[(i+j)%k]
 return v
counts={'minor_matrices':0,'modular_systems':0,'modular_rhs':0,'fraction_solutions':0,'collatz_words':0,'integral_collatz_words':0}
for k in range(1,5):
 for a in product([1,2],repeat=k):
  for p,u in product(range(1,5),range(-3,4)):
   if gcd(p,u)!=1:continue
   M=matrix(p,u,a);D=p**sum(a)-u**k
   check(determinant(M)==(-1)**(k-1)*D,('det',p,u,a))
   check(divisors(M)==[1]*(k-1)+[abs(D)],('minors',p,u,a))
   counts['minor_matrices']+=1
   for b in [tuple(range(k)),tuple((-1)**i for i in range(k))]:
    if D:
     x=[Fraction(numerator(p,u,a,b,i),D) for i in range(k)]
     check(all(sum(M[i][j]*x[j] for j in range(k))==b[i] for i in range(k)),('inverse',p,u,a,b))
     check(all(v.denominator==1 for v in x)==(numerator(p,u,a,b)%D==0),'integer obstruction')
     counts['fraction_solutions']+=1
   if k<=3 and p<=3:
    for m in range(1,6):
     image=Counter(tuple(sum(M[i][j]*x[j] for j in range(k))%m for i in range(k)) for x in product(range(m),repeat=k))
     g=gcd(abs(D),m)
     for b in product(range(m),repeat=k):
      predicted=g if numerator(p,u,a,b)%g==0 else 0
      check(image[b]==predicted,('modular fiber',p,u,a,m,b,image[b],predicted))
      counts['modular_rhs']+=1
     counts['modular_systems']+=1
# Actual closed Collatz words, not merely final affine endpoints.
for k in range(1,7):
 for a in product([1,2,3],repeat=k):
  D=2**sum(a)-3**k;b=(1,)*k;N=[numerator(2,3,a,b,i) for i in range(k)]
  counts['collatz_words']+=1
  if N[0]%D:continue
  x=[n//D for n in N];counts['integral_collatz_words']+=1
  check(all(v and v%2 for v in x),'nonzero odd coordinates')
  check(all((v>0)==(D>0) for v in x),'positivity sign')
  for i,v in enumerate(x):
   z=3*v+1;e=0
   while z%2==0:z//=2;e+=1
   check(e==a[i] and z==x[(i+1)%k],'exact accelerated step')
w=(1,1,1,2,1,1,4);v=(1,1,1,1,1,2,4)
check(sorted(w)==sorted(v),'same multiset')
D=2**sum(w)-3**len(w)
check(D==-139,'hostile determinant')
check(numerator(2,3,w,(1,)*7)==2363 and numerator(2,3,v,(1,)*7)==2123,'hostile numerators')
cycle=[numerator(2,3,w,(1,)*7,i)//D for i in range(7)]
check(cycle==[-17,-25,-37,-55,-41,-61,-91],'negative cycle')
# Signed graph compatibility: triangle and square with exact prescribed RHS.
triangle=[[1,1,0],[0,1,1],[1,0,1]];square=[[1,1,0,0],[0,1,1,0],[0,0,1,1],[1,0,0,1]]
for m in range(2,8):
 for M,b,expected in [(triangle,(1,0,0),1 if m%2 else 0),(square,(1,0,0,0),0)]:
  k=len(M)
  actual=sum(all(sum(M[i][j]*x[j] for j in range(k))%m==b[i] for i in range(k)) for x in product(range(m),repeat=k))
  check(actual==expected,('signed RHS',m))
print(json.dumps({'status':'PASS','counts':counts,'hostile_cycle':cycle,'producer_imported':False,'scope':'exact finite checks plus independent proofs; no Collatz convergence result'},indent=2,sort_keys=True))
