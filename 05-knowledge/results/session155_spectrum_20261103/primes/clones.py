from itertools import permutations
from math import comb

def paths(A):
 n=len(A); f=[[0]*n for _ in range(1<<n)]
 for x in range(n):f[1<<x][x]=1
 for m in range(1,1<<n):
  for x in range(n):
   if f[m][x]:
    for y in range(n):
     if not(m>>y&1) and A[x][y]:f[m|1<<y][y]+=f[m][x]
 return sum(f[-1])

def expansion(A,r,k):
 old=[x for x in range(len(A)) if x!=r];n=len(old)+k; B=[[0]*n for _ in range(n)]
 for i,x in enumerate(old):
  for j,y in enumerate(old):B[i][j]=A[x][y]
  for j in range(len(old),n):B[i][j]=A[x][r];B[j][i]=A[r][x]
 for i in range(len(old),n):
  for j in range(i+1,n):B[i][j]=1
 return B

def slot_count(A,r,k):
 old=[x for x in range(len(A)) if x!=r];total=0
 if not old:return 1
 for pi in permutations(old):
  t=A[r][pi[0]]+A[pi[-1]][r];q=0;bad=False
  for x,y in zip(pi,pi[1:]):
   valid=A[x][r] and A[r][y]
   t+=valid
   if not A[x][y]:q+=1;bad|=not valid
  if not bad:total+=sum((-1)**j*comb(q,j)*(t-j)**k for j in range(q+1))
 return total

def tournament(n,c):
 A=[[0]*n for _ in range(n)];p=0
 for i in range(n):
  for j in range(i+1,n):
   if c>>p&1:A[i][j]=1
   else:A[j][i]=1
   p+=1
 return A

N=0
for n in range(1,5):
 for c in range(1<<(n*(n-1)//2)):
  A=tournament(n,c)
  for r in range(n):
   for k in range(1,5):
    x,y=slot_count(A,r,k),paths(expansion(A,r,k));N+=1
    if x!=y:raise SystemExit(f'FAIL {n} {c} {r} {k} {x} {y}')
print('slot formula checks:',N,'PASS')
