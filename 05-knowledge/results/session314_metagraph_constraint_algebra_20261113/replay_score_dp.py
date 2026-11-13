from itertools import combinations
from collections import defaultdict
def build(n):
 m=n//2;W=[n**i for i in range(m)]
 s=sum((n-1-i)*a for i,a in enumerate(W));A={s:1}
 def mul(k):
  nonlocal A
  B=A.copy()
  for a,b in A.items():B[a+k]=B.get(a+k,0)+b
  A=B
 for i in range(m):
  for j in range(i+1,m):
   if j>=i+2:mul(W[j]-W[i])
   mul(-W[j]-W[i])
 for i in range(m):
  u=int( i<m-1 or n%2)+int(n%2 and i<m-1)
  for _ in range(u):mul(-W[i])
 return A
def calc(n):
 m=n//2;A=build(n);P=defaultdict(int);V={}
 def sc(S):return tuple(sorted(abs(2*x-n+1) for x in S))
 def isstrong(z):
  S=sorted([(n-1-x)//2 for x in z]+[(n-1+x)//2 for x in z]+([m] if n%2 else []))
  return all(sum(S[:i])>i*(i-1)//2 for i in range(1,n))
 for q,w in A.items():
  S=[q//n**i%n for i in range(m)]
  z=sc(S);V[z]=isstrong(z)
  S[0]-=1;t=sc(S);P[tuple(sorted((z,t)))]+=w
 for k,v in P.items():assert v%2==0
 P=set(P)
 D={tuple(sorted((u,v))) for u,v in combinations(V,2) if sum(abs(a-b) for a,b in zip(u,v))==2 and (V[u] or V[v])}
 print(n,len(A),len(V),len(P),len(D),sorted(D-P)[:10],flush=True)
for n in range(3,14,2):calc(n)
