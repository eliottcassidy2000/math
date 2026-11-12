from itertools import combinations,product,permutations
from functools import lru_cache
import json

def tours(k):
 e=list(combinations(range(k),2))
 for b in range(1<<len(e)):
  a=[0]*k
  for i,(u,v) in enumerate(e):
   if b>>i&1:a[u]|=1<<v
   else:a[v]|=1<<u
  yield tuple(a)
def wordgraph(A,w):
 k=len(A);n=k+len(w);t=list(A)+[0]*len(w)
 for i,l in enumerate(w,k):
  for a in range(k):
   if l>>a&1:t[a]|=1<<i
   else:t[i]|=1<<a
  t[i]|=sum(1<<j for j in range(i+1,n))
 return tuple(t)
def H(T):
 n=len(T);d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for S in range(1,1<<n):
  for i,x in enumerate(d[S]):
   R=T[i]&~S
   while R:
    b=R&-R;R-=b;d[S|b][b.bit_length()-1]+=x
 return sum(d[-1])
def strong(T):
 n=len(T)
 for i in range(n):
  seen=0;todo=1<<i
  while todo:
   b=todo&-todo;todo-=b;j=b.bit_length()-1
   if seen&b:continue
   seen|=b;todo|=T[j]&~seen
  if seen!=(1<<n)-1:return False
 return True
def bound(N):
 b=0;p=1
 while p*5<=N**3:b+=1;p*=5
 return b
@lru_cache(None)
def cores(A,N):
 k=len(A);out={}
 for m in range(max(-1,bound(N)+1-k)+1):
  for w in product(range(1<<k),repeat=m):
   T=wordgraph(A,w)
   if strong(T):
    h=H(T)
    if N%h==0:out[h,m]=out.get((h,m),0)+1
 return out

def partitions(k):
 for s in range(1,k+1):
  for a in product(range(s),repeat=k):
   if set(a)==set(range(s)):yield tuple(tuple(i for i in range(k) if a[i]==j) for j in range(s))
def induced(A,B):return tuple(sum(1<<j for j,b in enumerate(B) if A[a]>>b&1) for a in B)
def coeff(A,N,m):
 from math import comb
 z=0;k=len(A)
 for P in partitions(k):
  if any(not(A[a]>>b&1) for i,B in enumerate(P) for C in P[i+1:] for a in B for b in C):continue
  d={(1,0):1}
  for B in P:
   nxt={}
   for (h,l),x in d.items():
    for (hh,ll),y in cores(induced(A,B),N).items():
     if N%(h*hh)==0:nxt[h*hh,l+ll]=nxt.get((h*hh,l+ll),0)+x*y
   d=nxt
  for (h,l),x in d.items():
   if h==N and l<=m:
    assert l<=bound(N)+len(P)-k
    z+=x*comb(m-l+len(P),len(P))
 return z
checks=0
for k in range(1,4):
 for A in tours(k):
  for m in range(4):
   actual={}
   for w in product(range(1<<k),repeat=m):
    h=H(wordgraph(A,w));actual[h]=actual.get(h,0)+1
   for N in (1,3,5,7,9):
    want=coeff(A,N,m);assert want==actual.get(N,0),(A,N,m,want,actual.get(N,0));checks+=1
print(json.dumps(dict(status='PASS',coefficient_checks=checks,apex_orders=[1,2,3],targets=[1,3,5,7,9],spine_lengths=[0,1,2,3]),indent=2))
