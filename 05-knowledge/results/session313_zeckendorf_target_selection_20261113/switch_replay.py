from itertools import product
from switch_profiles import stats
def formula(L):
 H=0
 for i,(h,q,s,ss,e,ee,b) in enumerate(L):
  others=1
  for j in range(len(L)):
   if j!=i:others*=L[j][0]
  H+=others*(b+(s if i==0 else 0)+(e if i==len(L)-1 else 0))
  for j in range(i+1):
   if j==i and i+1<len(L):
    out=1
    for k in range(len(L)):
     if k not in (i,i+1):out*=L[k][0]
    H+=out*e*L[i+1][2]
   elif j<i:
    out=1
    for k in range(len(L)):
     if k<j or k>i:out*=L[k][0]
     if j<k<i:out*=L[k][1]
    H+=out*ee*L[j][3]
 return H
def raw(M):
 G=[1 if cut&1<<i else 0 for A,cut in [] for i in range(len(A))];G=[0];k=1
 for A,cut in M:
  for i in range(len(A)):
   if cut>>i&1:G[0]|=1<<(k+i)
   G.append((A[i]<<k)| (0 if cut>>i&1 else 1))
  for i in range(1,k):
   G[i]|=((1<<len(A))-1)<<k
  k+=len(A)
 n=len(G);N=1<<n;D=[0]*(N*n)
 for j in range(n):D[(1<<j)*n+j]=1
 for S in range(1,N):
  for j in range(n):
   if D[S*n+j]:
    Q=G[j]&~S
    while Q:
     L=Q&-Q;Q^=L;D[(S|L)*n+L.bit_length()-1]+=D[S*n+j]
 return sum(D[(N-1)*n:])
for A in ((20,1,2,7,14),(24,17,3,6,12),(40,1,3,6,15,30)):
 a,b={5:((1,4) if A[0]==20 else (3,5)),6:(17,40)}[len(A)]
 h,q,d=stats(A)
 def prof(c):
  (s,ss),(e,ee,k)=d[c]
  return (h,q,s,ss,e,ee,k)
 f={}
 for C in product((a,b),repeat=2):
  r=raw([(A,c) for c in C]);f[C]=formula([prof(c) for c in C])
  assert f[C]==r,(A,C,f[C],r)
 assert f[a,a]+f[b,b]==f[a,b]+f[b,a]
 print(A,h,q,d[a],d[b],f)
