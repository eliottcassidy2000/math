from itertools import product
from math import isqrt
import json

def require(p,msg):
 if not p:raise ValueError(msg)
def polygon(n,bits):
 A=[1<<((i+1)%n) for i in range(n)]
 for i,s in enumerate(bits):
  j=(i+2)%n
  if s:A[i]|=1<<j
  else:A[j]|=1<<i
 return A
def outside_paths(n,bits):
 if all(bits) or not any(bits):
  d=2 if bits[0] else -2
  return [[(i+d*k)%n for k in range(n)] for i in range(n)]
 paths=[]
 for i,s in enumerate(bits):
  if s and not bits[(i+1)%n]:
   paths.append([(i+k)%n for k in range(4,n+1)]+[(i+2)%n,(i+3)%n,(i+1)%n])
 return paths
for n in (5,7,9,11):
 for z in product((0,1),repeat=n):
  A=polygon(n,z); P=outside_paths(n,z)
  require(len(P)==len(set(map(tuple,P)))>0,(n,z))
  for q in P:
   require(sorted(q)==list(range(n)),(n,q))
   require(all(A[i]>>j&1 for i,j in zip(q,q[1:])),(n,z,q))
   require(any((j-i)%n!=1 for i,j in zip(q,q[1:])),(n,q))

def F(u):
 g=1
 while u:
  x=u//2
  if not u&1:g+=x+(x&1)
  u=x
 return g
def independentF(u):
 if u==0:return 1
 n=u.bit_length();value=1+int(u%2==0)
 for a in range(n):
  if not u>>a&1:continue
  for b in range(n):
   if u>>b&1:continue
   if a>b:value+=1<<(a-b-1)
   elif b==a+1:value+=1
 return value

def C(v):
 n=v.bit_length()
 return [sum(1<<j for j in range(i+1,n))+(((v>>i)&1)<<n) for i in range(n)]+[((1<<n)-1)^v]
def join(b,e):
 A,B=C(b-1),C(e-1)
 l,r=len(A),len(B);x=l-1;y=l+r-1
 R=A+[(sum(((B[j]>>i)&1)<<j for j in range(r))<<l)|((1<<l)-1) for i in range(r)]
 R[y]^=1<<x;R[x]|=1<<y
 return R
def h(A):
 n=len(A);D=[[0]*n for _ in range(1<<n)]
 for i in range(n):D[1<<i][i]=1
 for S,row in enumerate(D):
  for i,d in enumerate(row):
   if not d:continue
   w=A[i]&~S
   while w:
    b=w&-w;w^=b;D[S|b][b.bit_length()-1]+=d
 return sum(D[-1])
def val(b,e):
 return b*e+F(b-1)*F(e-1)-(b%2)*(e%2)
for b,e in product(range(1,17),repeat=2):
 require(h(join(b,e))==val(b,e),(b,e))
for p,b,e in [(1000000007,3,191311747),(1000000009,7,86246185),
(1000000000039,1613,583862840),(1000000000061,7,88583919141),
(10000000000037,19,404043168698)]:
 require(all(p%d for d in range(2,isqrt(p)+1)),p)
 require(F(b-1)==independentF(b-1) and F(e-1)==independentF(e-1),(b,e))
 R=join(b,e);n=len(R)
 require(all((R[i]>>j&1)+(R[j]>>i&1)==1 for i in range(n) for j in range(i+1,n)),R)
 require(val(b,e)==p,(p,b,e))
 print(json.dumps(dict(prime=p,addresses=[b,e],order=n,rows=R)))
