from arcscan import stats, strong, rows
from itertools import permutations

def gate(T,u,v):
 n=len(T)
 return [r+((i==u)<<n)+((i!=v)<<(n+1)) for i,r in enumerate(T)] + [((1<<(n+2))-1)^((1<<u)|(1<<n)), 1<<v]

def verify(T,u,v):
 h,M,*_=stats(T)
 internal=0
 for p in permutations(range(len(T))):
  if not all(T[a]>>b&1 for a,b in zip(p,p[1:])):continue
  if any((a,b)==(u,v) for a,b in zip(p[1:-1],p[2:-1])):internal+=1
 G=gate(T,u,v);D=gate(G,u,v)
 if stats(G)[1][u][v]!=internal or stats(D)[1][u][v]!=internal:raise ValueError((T,u,v))
 if strong(T) and not (strong(G) and strong(D)):raise ValueError("strongness")
 return M[u][v],internal
for n in range(2,5):
 for x in range(1<<(n*(n-1)//2)):
  T=rows(n,x)
  for u in range(n):
   for v in range(n):
    if T[u]>>v&1:verify(T,u,v)
for m in range(1,7):
 n=m+2;N=(1<<n)-1
 T=[2,N&~3]+[1+sum(1<<j for j in range(i+1,n)) for i in range(2,n)]
 for _ in range(2):
  T=gate(T,0,1)
  if not strong(T) or stats(T)[1][0][1]!=2**m-2:raise ValueError(m)
print("all marked arcs through order four, two gate applications; wrap controls m=1..6")
