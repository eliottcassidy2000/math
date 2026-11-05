from itertools import product
import json

def req(x,m):
 if not x:raise RuntimeError(m)
def g(n):
 if n==1:return 1
 if n%2==0:return g(n//2)
 l=(n+1)//2;return g(l)+l-l%2

def O(v):D,E,U,d=v;return D,0,2*U,d
def Z(v):D,E,U,d=v;return D-E+U,d,2*U-d,d
def J(v):D,E,U,d=v;return D-E+U,D-E,D+U,D

def direct(v,k,A):
 for i in range(k):v=(Z if A>>i&1 else O)(v)
 return v

def reverse(A,k):return sum(((A>>i)&1)<<(k-i-1)for i in range(k))
def block(v,k,A):
 D,E,U,d=v;B=1<<k;R=reverse(A,k);p=A%2;t=A>>(k-1)
 return D+(U-d)*A+d*(g(A+1)-1)+(d-E)*p,d*t,B*U-d*R,d

def readout(v,f):
 D,E,U,d=v
 return 3*D-2*E+6*U-2*d if f=='N' else 5*D-3*E+13*U-2*d

def scalar(v,k,A,f):
 D,E,U,d=v;B=1<<k;R=reverse(A,k);p=A%2;t=A>>(k-1);G=g(A+1)
 if f=='N':return 3*D+3*(U-d)*A+6*B*U+d*(3*G-6*R-2*t-5)+3*(d-E)*p
 return 5*D+5*(U-d)*A+13*B*U+d*(5*G-13*R-3*t-7)+5*(d-E)*p

basis=[tuple(int(i==j)for i in range(4))for j in range(4)]
physical={(1,0,2,1)}
level={(1,0,2,1)}
for _ in range(5):
 level={X(v)for v in level for X in(O,Z,J)};physical.update(level)
checks=0;complements=0
for k in range(1,13):
 for A in range(1<<k):
  B=1<<k;R=reverse(A,k);p=A%2;t=A>>(k-1)
  req(g(A+1)-g(B-A)==A-R+t-p,'complement')
  complements+=1
  for v in basis:
   req(direct(v,k,A)==block(v,k,A),'linear block basis')
   for f in('N','M'):req(readout(direct(v,k,A),f)==scalar(v,k,A,f),'scalar basis')
   checks+=1
physical_checks=0
for v in physical:
 for k in range(1,7):
  for A in range(1<<k):
   req(direct(v,k,A)==block(v,k,A),'physical block')
   physical_checks+=1
print(json.dumps(dict(status='PASS',complement_tests=complements,linear_basis_block_tests=checks,physical_prefixes=len(physical),physical_block_tests=physical_checks),sort_keys=True,indent=2))
