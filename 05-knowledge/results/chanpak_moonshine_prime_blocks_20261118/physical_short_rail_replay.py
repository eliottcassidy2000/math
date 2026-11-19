def add(A,bit):
 n=len(A)
 return [r|((1<<n) if (i!=0 or bit==0) else 0) for i,r in enumerate(A)]+[bit]
def profile(A):
 n=len(A);L=1<<n;P=[[0]*n for _ in range(L)]
 for j in range(n):P[1<<j][j]=1
 for s,row in enumerate(P):
  for j,v in enumerate(row):
   if v:
    m=A[j]&~s
    while m:
     b=m&-m;m^=b;P[s|b][b.bit_length()-1]+=v
 back=[sum(((A[j]>>i)&1)<<j for j in range(n)) for i in range(n)]
 Q=[[0]*n for _ in range(L)];Q[1][0]=1
 for s,row in enumerate(Q):
  for j,v in enumerate(row):
   if v:
    m=A[j]&~s
    while m:
     b=m&-m;m^=b;Q[s|b][b.bit_length()-1]+=v
 h=[1]+list(map(sum,P[1:]))
 return h[-1],P[-1][0],sum(sum(Q[s])*h[(L-1)^s] for s in range(1,L,2)),h[L-2]
def block(A):
 for b in (0,1,1):A=add(A,b)
 return A
A=block(block([2,0]));p=101
for i in range(3):
 h,e,t,d=profile(A)
 if (e,t,d)!=(1,p,1) or 7*h!=6*p+3+4*i:raise ValueError((i,h,e,t,d))
 B=add(A,1)
 n=len(A)
 C=B.copy();C[n]^=1
 if profile(B)[0]!=h+t-e or profile(B)[0]-profile(C)[0]!=t:raise ValueError(i)
 print(i,len(A),h,t,h+t-e)
 A=block(A);p=8*p-3

def ends(A):
 n=len(A);L=1<<n;D=[[0]*n for _ in range(L)]
 for j in range(n):D[1<<j][j]=1
 for s,row in enumerate(D):
  for j,v in enumerate(row):
   if not v:continue
   for k in range(n):
    if A[j]>>k&1 and not s>>k&1:D[s|1<<k][k]+=v
 return D
A=block(block([2,0]));p=101
for step in range(3):
 n=len(A);F=(1<<n)-1
 E=ends(A); R=ends([sum(((A[j]>>i)&1)<<j for j in range(n)) for i in range(n)])
 W=[[0]*n for _ in range(n)]
 for s in range(1,F):
  for i,v in enumerate(E[s]):
   if v:
    for j,w in enumerate(R[F^s]):
     W[i][j]+=v*w
 hits=[]
 for s in range(F+1):
  d=sum(R[F][i] if s>>i&1 else E[F][i] for i in range(n))
  d+=sum(W[i][j] for i in range(n) if not s>>i&1 for j in range(n) if s>>j&1)
  if d==p:hits.append(s)
 if hits:raise ValueError((step,hits))
 print(n,p,len(hits),'one-vertex scalar readers',flush=True)
 A=block(A);p=8*p-3
