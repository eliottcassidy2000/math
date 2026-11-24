S={5:(6,8,10,1),17:(30,60,56,48,32,1),37:(30,44,48,52,34,1),77:(222,252,248,240,96,129,160,16),129:(254,508,504,496,480,448,384,256,1),751:(150,504,506,257,40,73,153,312,113)}
def H(A,prec=False):
 n=len(A);F=(1<<n)-1
 D=[[0]*n for _ in range(1<<n)]; C=[0]*n
 for i in range(n):D[1<<i][i]=1
 for s,row in enumerate(D):
  for i,a in enumerate(row):
   if not a:continue
   u=A[i]&~s
   while u:
    b=u&-u;u-=b;D[s|b][b.bit_length()-1]+=a
 h=sum(D[-1])
 if not prec:return h
 # all oriented Hamilton words
 def dfs(s,v):
  if s==F:return 1
  ans=0;u=A[v]&~s
  while u:
   b=u&-u;u-=b;j=b.bit_length()-1
   x=dfs(s|b,j);ans+=x
   if x&1:
    a=s
    while a:
     c=a&-a;a-=c;C[c.bit_length()-1]^=b
  return ans
 for i in range(n):dfs(1<<i,i)
 return h,C,sorted((a.bit_count() for a in C))==list(range(n))
def flip(A,i,j):
 R=list(A);R[i]^=1<<j;R[j]^=1<<i;return tuple(R)
def apex(A,mask):
 n=len(A);return tuple(a|((not(mask>>i&1))<<n) for i,a in enumerate(A))+(mask,)
from itertools import combinations
found = []
for h,A in S.items():
 for u,v in combinations(range(len(A)),2):
  B=flip(A,u,v); q=H(B)
  if abs(q-h)==2:
   _, D, good=H(B,True)
   tri=next((i,j,k) for i,j,k in combinations(range(len(A)),3) if (D[i]>>j&1)==(D[j]>>k&1)!=(D[i]>>k&1))
   assert not good
   found.append((h,q,(u,v),tri))
assert found == [(5,3,(0,2),(0,1,3)), (5,3,(1,3),(0,2,3)),
                 (17,15,(0,2),(0,1,2)), (17,19,(0,3),(0,1,4)),
                 (17,19,(2,5),(0,1,2)), (17,15,(3,5),(3,4,5))]

def products(factors, bound):
    if not factors:
        yield 1
        return
    f, *rest = factors
    power = 1
    while power <= bound:
        for n in products(rest, bound // power):
            yield n * power
        power *= f

values = set(products(tuple(S), 10**40))
assert len(values) == 755771
assert not any(x+2 in values for x in values)
print(len(values), "six-bank values through 10^40; no gap two")
