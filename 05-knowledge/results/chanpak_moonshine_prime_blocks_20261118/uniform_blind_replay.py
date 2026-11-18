from itertools import combinations,permutations
from functools import lru_cache
from math import comb
L=(20,1,2,7,14);R=(16,1,3,7,14)
def dp(A):
 n=len(A);out=[[0]*n for _ in range(1<<n)]
 for i in range(n):out[1<<i][i]=1
 for S in range(1,1<<n):
  for i,v in enumerate(out[S]):
   todo=A[i]&~S
   while v and todo:
    b=todo&-todo;todo-=b;out[S|b][b.bit_length()-1]+=v
 return out
def H(A):return sum(dp(A)[-1])
def reverse(A):return tuple(((1<<len(A))-1)^a^(1<<i) for i,a in enumerate(A))
def profile(A,r):
 E=dp(A);F=dp(reverse(A));M=(1<<len(A))-1;hs=[1]+[sum(x) for x in E[1:]]
 return hs[M],E[M][r],sum(F[S][r]*hs[M^S] for S in range(1,M+1)),hs[M^(1<<r)]
def covers(A,q):
 tab=dp(A);hs=[1]+[sum(x) for x in tab[1:]];M=len(tab)-1
 @lru_cache(None)
 def go(S,k):
  if not k:return int(S==0)
  ans=0;U=S
  while True:
   ans+=hs[U]*go(S^U,k-1)
   if not U:return ans
   U=(U-1)&S
 return go(M,q)
def cycles(A,k):
 result=[]
 for S in combinations(range(len(A)),k):
  for tail in permutations(S[1:]):
   c=(S[0],)+tail
   if all(A[c[i]]>>c[(i+1)%k]&1 for i in range(k)):result.append(c)
 return result
for A in(L,R):
 assert len(cycles(A,3))==3 and len(cycles(A,5))==1
 for q in range(7):assert covers(A,q)==q**5+6*q**3+2*q
for B in(R,reverse(R)):
 assert not any(all(((L[i]>>j)&1)==((B[p[i]]>>p[j])&1) for i in range(5) for j in range(5)) for p in permutations(range(5)))
assert profile(L,3)==(9,1,23,5) and profile(R,3)==(9,1,25,5)
def substitute(host,mark,A):
 groups=[];at=0
 for i in range(len(host)):
  size=len(A) if i==mark else 1;groups.append(list(range(at,at+size)));at+=size
 out=[0]*at
 for i in range(len(host)):
  for j in range(len(host)):
   if host[i]>>j&1:
    for u in groups[i]:
     for v in groups[j]:out[u]|=1<<v
 for i,u in enumerate(groups[mark]):
  for j,v in enumerate(groups[mark]):
   if A[i]>>j&1:out[u]|=1<<v
 return tuple(out)
for host,mark,expected in[((10,12,9,0),0,93),((2,4,9,3),0,177),((2,4,9,3),1,243)]:
 assert H(substitute(host,mark,L))==H(substitute(host,mark,R))==expected
# Every labeled four-vertex host with its vertex0 replaced:64 hosts.
E=list(combinations(range(4),2))
for code in range(64):
 host=[0]*4
 for k,(i,j) in enumerate(E):
  if code>>k&1:host[i]|=1<<j
  else:host[j]|=1<<i
 assert H(substitute(host,0,L))==H(substitute(host,0,R))
def append(A,r,W):
 A=list(A)
 for k in range(W.bit_length()):
  bit=W>>k&1;n=len(A)
  A=[(a|1<<n)^((bit if i==r else 0)<<n) for i,a in enumerate(A)]+[bit<<r]
 return tuple(A)
def K(W):
 s=0
 while W:
  v=W//2
  if W&1:s+=v+(v&1)
  W=v
 return s
for W in(0,1,2,3,5,17,31):
 a=append(L,3,W);b=append(R,3,W)
 hl=H(a);hr=H(b)
 assert hl==9+23*W-5*K(W)-(W&1)
 assert hr-hl==2*W
 assert [(i,j) for i in range(len(a)) for j in range(i+1,len(a)) if (a[i]>>j&1)!=(b[i]>>j&1)]==[(0,2)]
 print('marked word',W,hl,hr)
print('full slot equality, uniform contexts, nonconverse control, marked2W: PASS')
