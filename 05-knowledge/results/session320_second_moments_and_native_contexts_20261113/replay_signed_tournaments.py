def info(A):
 n=len(A);N=(1<<n)-1
 dp=[[0]*n for _ in range(N+1)];h=[1]+[0]*N
 for j in range(n):dp[1<<j][j]=1
 for s in range(1,N+1):
  h[s]=sum(dp[s])
  for j,d in enumerate(dp[s]):
   if not d:continue
   a=A[j]&~s
   while a:
    b=a&-a;a-=b;dp[s|b][b.bit_length()-1]+=d
 def v(s):
  p=s; c=0
  while True:
   c+=h[p]*h[s^p]
   if p==0:return c
   p=(p-1)&s
 V=[v(s)for s in range(N+1)]
 return h[N],V[N],sum(h[s]*V[N^s]for s in range(N+1))
def strong(A):
 n=len(A)
 for j in range(n):
  s=1<<j
  for _ in A:
   for i,r in enumerate(A):
    if s>>i&1:s|=r
  if s!=(1<<n)-1:return False
 return True

A=(452,349,136,385,141,351,412,34,148)
B=(16,413,209,277,160,15,315,105,181)
def wrap(T):
 n=len(T)
 return tuple(a+(1<<n) for a in T)+(1<<(n+1),(1<<n)-1)
for T, counts in ((A,(931,19556,189849)),(B,(931,19388,191097))):
 n=len(T)
 if any((T[i]>>j&1)==(T[j]>>i&1) for i in range(n) for j in range(i+1,n)):
  raise ValueError('not a tournament')
 if not strong(T) or info(T)!=counts:raise ValueError(T)
w=info(wrap(A))[0]-info(wrap(B))[0]
z=info(wrap(wrap(A)))[0]-info(wrap(wrap(B)))[0]
if (w,z)!=(168,-1656):raise ValueError((w,z))
print('literal counts and both signed strong continuations verified')

from itertools import permutations
C=(24,9,3,4,14)
X=(388,449,202,195,335,351,129,48,204)
Y=(284,33,506,418,138,17,59,355,114)
def inflate(R,U):
    m=len(U); q=(0,)*m+tuple(range(1,len(R)))
    return tuple(sum(1<<j for j,b in enumerate(q) if i!=j and
                     ((U[i]>>j&1) if a==b else (R[a]>>b&1)))
                 for i,a in enumerate(q))
def hq2(A):
    n=len(A);N=(1<<n)-1
    dp=[[0]*n for _ in range(N+1)];h=[1]+[0]*N
    for j in range(n):dp[1<<j][j]=1
    for s in range(1,N+1):
        h[s]=sum(dp[s])
        for j,d in enumerate(dp[s]):
            if not d:continue
            a=A[j]&~s
            while a:
                b=a&-a;a-=b;dp[s|b][b.bit_length()-1]+=d
    return h[N],sum(h[s]*h[N^s]for s in range(N+1))
def coeffs(R):
    n=len(R); c=[]
    for j in (1,2,3):
        total=0
        for p in permutations(range(n+j-1)):
            if [a for a in p if a==0 or a>=n]!=[0]+list(range(n,n+j-1)):continue
            w=[a if a<n else 0 for a in p]
            total+=all(R[a]>>b&1 for a,b in zip(w,w[1:]))
        c.append(total)
    return c
if not (strong(C) and strong(X) and strong(Y)) or coeffs(C)!=[13,12,2]:raise ValueError('context')
if info(X)!=(1047,21492,204453) or info(Y)!=(1047,21252,205173):raise ValueError('source')
P,Q=inflate(C,X),inflate(C,Y)
px,py=hq2(P),hq2(Q)
if not (strong(P) and strong(Q)) or px[0]!=532623 or py[0]!=532623 or px[1]==py[1]:raise ValueError((px,py))
print('exact common strong return, output H and Q2:',px,py)
