from collections import deque
from itertools import product
from math import factorial
def poly(A):
 n=len(A)
 F=(1<<n)-1
 d=[[0]*n for _ in range(F+1)]
 for i in range(n):d[1<<i][i]=1
 for u in range(1,F+1):
  for j,a in enumerate(d[u]):
   for k in range(n):
    if (A[j]>>k)&1 and not(u>>k&1):d[u|1<<k][k]+=a
 H=[1]+list(map(sum,d[1:]))
 prev=[1]+[0]*F
 fd=[0]*(n+1); f=[1]
 for j in range(1,n+1):
  cur=[0]*(F+1)
  for s in range(1,F+1):
   u=s
   while u:
    cur[s]+=H[u]*prev[s^u]
    u=(u-1)&s
  a=f+[0]
  f=[(a[i-1] if i else 0)-(j-1)*a[i] for i in range(j+1)]
  for i,q in enumerate(f):fd[i]+=q*cur[F]//factorial(j)
  prev=cur
 return fd
def val(f,x):return sum(a*x**i for i,a in enumerate(f))
S3=[(10,12,9,16,39,15),(62,4,8,18,38,14),(62,36,24,18,34,12),(58,12,9,16,38,14),(6,4,8,19,39,15)]
S5=[(62,36,24,18,34,12),(234,220,25,144,225,206,140,4),(170,244,201,18,229,76,137,40)]
def test(p, seeds):
 d=p//2
 q=deque(); seen=set()
 def put(x):
  if x not in seen:seen.add(x);q.append(x)
 for A in seeds:
  f=poly(A);a=tuple(val(f,i) for i in range(1,d+1))
  assert all(x%p==0 for x in a)
  s=tuple(x//p%p for x in a)+tuple(val([j*f[j] for j in range(1,len(f))],i)%p for i in range(1,d+1))
  for v in product(range(p), repeat=d):
   put(tuple(x*v[j%d]%p for j,x in enumerate(s)))
  while q:
   z=q.popleft();a=[0]+list(z[:d])+[(z[d-1]-z[-1])%p];b=[0]+list(z[d:])+[-z[-1]]
   put(tuple(i*i*a[i]%p for i in range(1,d+1))+tuple(i*i*b[i]%p for i in range(1,d+1)))
   put(tuple((i*i*a[i]+i*(a[i+1]-a[i-1]))%p for i in range(1,d+1))+tuple((i*i*b[i]+i*(b[i+1]-b[i-1]))%p for i in range(1,d+1)))
  print(p,s,len(seen))
 assert len(seen)==p**(2*d)
if __name__=='__main__':
 test(3,S3)
 test(5,S5)
