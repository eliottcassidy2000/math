from itertools import product
from functools import lru_cache
from math import comb
A=[14,92,8,48,101,71,13]
B=[88,17,91,98,40,71,18]
def plus(c,r):
 return c[:r]+(c[r]+1,)+c[r+1:]
def coeff(T):
 n=len(T);nr=(1<<(n-3))-1
 # The source marks are transitive in the displayed order (or a permutation).
 if sum(((T[i]>>j)&1 and (T[j]>>k)&1 and (T[k]>>i)&1)
         for i in range(3) for j in range(3) for k in range(3)):
  raise ValueError("cyclic marks")
 @lru_cache(None)
 def f(s,c,last):
  W=T[last] if last<n else (1<<n)-1
  acc={c:1} if s==nr and all(c) else {}
  for j in range(n):
   if not W>>j&1:continue
   if j<3:
    if c[j]>=(n-2):continue
    more=f(s,plus(c,j),j)
   else:
    b=1<<(j-3)
    if s&b:continue
    more=f(s|b,c,j)
   for k,a in more.items():acc[k]=acc.get(k,0)+a
  return acc
 return f(0,(0,0,0),n)
def surj(m,r):
 return sum((-1)**(r-i)*comb(r,i)*i**m for i in range(r+1))
def expand(T,m):
 L=sum(([i]*x for i,x in enumerate(m)),[])+list(range(3,len(T)))
 return [sum(1<<j for j,b in enumerate(L) if (i<j if a==b else (T[a]>>b)&1)) for i,a in enumerate(L)]
def h(T):
 n=len(T);D=[[0]*n for _ in range(1<<n)]
 for i in range(n):D[1<<i][i]=1
 for s,row in enumerate(D):
  for i,d in enumerate(row):
   if d:
    w=T[i]&~s
    while w:
     b=w&-w;w^=b;D[s|b][b.bit_length()-1]+=d
 return sum(D[-1])
def val(c,m):
 return sum(d*surj(m[0],r[0])*surj(m[1],r[1])*surj(m[2],r[2]) for r,d in c.items())
c,d=coeff(A),coeff(B)
if c[(1,1,1)]!=89 or d!={**c,(1,1,1):91}:raise ValueError((c,d))
print("Exact lower certificate:",sorted(c.items()))
for m in product(range(1,4),repeat=3):
 x=h(expand(A,m));y=h(expand(B,m))
 if x!=val(c,m) or y!=val(d,m) or y-x!=2:raise ValueError((m,x,y))
print("One-clock first steps:",[(m,val(c,m),val(d,m)) for m in ((1,1,1),(2,1,1),(1,2,1),(1,1,2))])
print("Passed 54 independent literal inflation counts")

from collections import defaultdict
E=defaultdict(int)
for r,t in c.items():
 for q in product(*(range(1,a+1) for a in r)):
  d=t
  for a,b in zip(r,q):d*=(-1)**(a-b)*comb(a,b)
  E[q]+=d
E[(1,1,1)]+=1
V={(1,1):8,(1,2):-7,(1,3):2,(2,1):-9,(2,2):11,(2,3):-4,(3,1):2,(3,2):-4,(3,3):2}
W={(1,1):24,(1,2):-28,(2,1):-28,(1,3):8,(3,1):8,(2,2):66,(2,3):-40,(3,2):-40,(3,3):44,(2,4):6,(4,2):6,(3,4):-12,(4,3):-12,(4,4):6}
if {q:d for q,d in E.items() if d}!={**{(a,b,1):d for (a,b),d in V.items()},**{(a,b,2):d for (a,b),d in W.items()}}:raise ValueError('exponential expansion')
