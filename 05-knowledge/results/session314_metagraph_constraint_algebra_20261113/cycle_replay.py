from itertools import product
from collections import Counter
from math import gcd

def brute(n):
 m=n//2; D=list(range(2,m+1)) if n%2 else list(range(2,m))
 b=len(D)*n+(m if n%2==0 else 0);A=Counter();Q=Counter()
 for bits in product(range(2),repeat=b):
  X={d:bits[i*n:(i+1)*n] for i,d in enumerate(D)}
  if n%2==0:
   Y=bits[-m:]
   X[m]=Y+tuple(1-z for z in Y)
  g=next(g for g in range(1,n+1) if n%g==0 and all(all(w[i]==w[(i+g)%n] for i in range(n)) for w in X.values()))
  A[g]+=1
  for k in range(n):
   if all(all(w[i]==w[(k-i-d)%n] for i in range(n)) for d,w in X.items()):
    Q[k,g]+=1
 return A,Q

def theory(n):
 m=n//2; D=list(range(2,m+1)) if n%2 else list(range(2,m))
 A=Counter();Q=Counter()
 for g in (x for x in range(1,n+1) if n%x==0):
  ex0=g*len(D)
  okay=(n%2 or m%g)
  A[g]=(1<<(ex0+(0 if n%2 else g//2))) if okay else 0
  A[g]-=sum(e for h,e in A.items() if h<g and g%h==0)
  for k in range(n):
   f=lambda d:sum(2*i%g==(k-d)%g for i in range(g))
   ex=sum((g+f(d))//2 for d in D)
   if n%2==0:
    if not (g%2==0 or not okay):raise ValueError
    if okay and k%2:ex+=(g+f(m))//4
   Q[k,g]=(1<<ex) if okay and (n%2 or k%2) else 0
   Q[k,g]-=sum(e for (a,h),e in Q.items() if a==k and h<g and g%h==0)
 return A,Q

for n in range(3,8):
 A,Q=theory(n)
 assert (A,Q)==brute(n)
 B=Counter({g:sum(Q[k,g] for k in range(n))//n for g in A})
 assert all(sum(Q[k,g] for k in range(n))%n==0 and A[g]%g==0 and (A[g]//g-B[g])%2==0 for g in A)
 a=sum(A.values());cl=sum(A[g]//g for g in A);b=sum(B.values())
 assert b*(1 if n%2 else 2)==2**(((n-1)**2//4)-1)
 print(n,a,cl,b,tuple((g,A[g]//g,B[g]) for g in A if A[g]))
