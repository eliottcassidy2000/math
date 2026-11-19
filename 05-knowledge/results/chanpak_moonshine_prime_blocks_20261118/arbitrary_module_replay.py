from math import isqrt

U=[16,37,9,179,38,69,31,375,127]
mark=8
deleted=4

def require(cond,msg):
 if not cond:raise ValueError(msg)
def induced(A,I):
 return [sum(((A[i]>>j)&1)<<k for k,j in enumerate(I)) for i in I]
I=[i for i in range(len(U)) if i!=deleted]
V=induced(U,I)
vmark=I.index(mark)

def profiles(A):
 n=len(A);L=1<<n
 D=[[0]*n for _ in range(L)]
 for i in range(n):D[1<<i][i]=1
 for s,row in enumerate(D):
  for i,z in enumerate(row):
   if not z:continue
   remaining=A[i]&~s
   while remaining:
    j=remaining&-remaining
    remaining^=j;D[s|j][j.bit_length()-1]+=z
 hs=[1]+list(map(sum,D[1:]))
 return hs[-1],sum(y*hs[L-1^i] for i,y in enumerate(hs))
def covers(A,v):
 n=len(A);R=(n+1)//2;L=1<<n;start=n
 D=[[[0]*(n+1) for _ in range(R+1)] for _ in range(L)]
 D[0][0][start]=1
 for s in range(L):
  if s>>v&1:continue
  for r in range(R+1):
   for i,count in enumerate(D[s][r]):
    if not count:continue
    if r<R and (i==start or A[i]>>v&1):D[s][r+1][v]+=count
    for j in range(n):
     if j!=v and not s>>j&1 and (i==start or A[i]>>j&1):
      D[s|1<<j][r][j]+=count
 return [sum(a) for a in D[L-1^(1<<v)][1:]]
def substitute(A,v,T):
 L=[(i,j) for i in range(len(A)) for j in range(len(T) if i==v else 1)]
 return [sum(1<<q for q,(j,b) in enumerate(L) if (T[a]>>b&1 if i==j==v else (A[i]>>j&1))) for i,a in L]
def tt(m):return [(1<<i)-1 for i in range(m)]
def isprime(n):return n>1 and all(n%d for d in range(2,isqrt(n)+1))
def strong(A):
 for i in range(len(A)):
  seen=1<<i
  for _ in A:
   for j,y in enumerate(A):
    if seen>>j&1:seen|=y
  if seen!=(1<<len(A))-1:return False
 return True

require(covers(U,mark)==[343,126,0,0,0],"child coefficients")
require(covers(V,vmark)==[49,18,0,0],"source coefficients")
require(strong(U),"child quotient strong")
print('Marked coefficient vectors',covers(V,vmark),covers(U,mark))
for T in [tt(m) for m in range(1,6)]+[[8,1,3,6],[32,1,3,7,15,30],[48,41,3,5,14,28]]:
 h,p=profiles(T);N=13*h+18*p
 F=substitute(U,mark,T);G=substitute(V,vmark,T)
 require(profiles(F)[0]==7*N and profiles(G)[0]==N,(T,N))
 J=[i for i in range(len(F)) if i!=deleted]  # mark is after deleted
 require(induced(F,J)==G,"not literal deletion")
 require(strong(F),"inflated child not strong")
 if len(T)>=2:require(h%72==(61*N)%72,(T,h,N))
 print((len(T),h,p),N,7*N,'prime input' if isprime(N) else 'composite input')
