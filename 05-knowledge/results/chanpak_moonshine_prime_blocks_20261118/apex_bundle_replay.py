from itertools import combinations
def dp(A):
 n=len(A);d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for S,row in enumerate(d):
  for i,v in enumerate(row):
   if not v:continue
   rest=A[i]&~S
   while rest:
    b=rest&-rest;rest^=b
    d[S|b][b.bit_length()-1]+=v
 return d
def profile(A,v):
 I=[i for i in range(len(A)) if i!=v]; n=len(I);F=(1<<n)-1
 R=[sum((A[i]>>j&1)<<k for k,j in enumerate(I)) for i in I]
 C=[sum((R[j]>>i&1)<<j for j in range(n)) for i in range(n)]
 E=dp(R); B=dp(C)
 old=sum((A[v]>>j&1)<<k for k,j in enumerate(I))
 W=[[sum(E[S][i]*B[F^S][j] for S in range(1,F)) for j in range(n)] for i in range(n)]
 return [sum((B[-1][i] if m>>i&1 else E[-1][i]) for i in range(n))+sum(W[i][j] for i in range(n) for j in range(n) if not m>>i&1 and m>>j&1) for m in (old^b for b in range(F+1))]

def flip(A,v,M):
 R=A.copy(); I=[i for i in range(len(A)) if i!=v]
 for k,i in enumerate(I):
  if M>>k&1:
   R[v]^=1<<i; R[i]^=1<<v
 return R
def H(A):return sum(dp(A)[-1])
def test(A,B,v,U,V,expect):
 if U&V:raise ValueError("overlap")
 corners=(0,U,V,U^V)
 vals=[(H(flip(A,v,m)),H(flip(B,v,m))) for m in corners]
 if vals!=expect:raise ValueError(vals)
 for Z in (A,B):
  for m,h in enumerate(profile(Z,v)):
   if H(flip(Z,v,m))!=h:raise ValueError(m)
 print(vals)
 print("mixed",*[vals[3][i]-vals[1][i]-vals[2][i]+vals[0][i] for i in range(2)])

test([0,1,3,7,15],[10,0,3,6,15],4,3,4,[(1,3),(11,13),(3,5),(9,9)])
A=[38,64,34,103,111,2,37]
test(flip(A,6,30),flip(A,6,31),6,3,36,[(9,11),(37,39),(17,19),(45,47)])
A=[50,224,123,179,130,16,57,101]
test(flip(A,7,3),flip(A,7,112),7,48,78,[(173,175),(171,173),(71,73),(69,71)])
