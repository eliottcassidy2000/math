def require(x):
 if not x: raise AssertionError("bridge check failed")
from itertools import permutations
from collections import Counter
from math import comb
ROWS=(62,60,120,16,96,8,43)
def arc(a,b):return (ROWS[a]>>b)&1
def cert():
 c=[0]*8; hist=Counter()
 for p in permutations(range(1,7)):
  slots=[arc(0,p[0])]+[arc(a,0)*arc(0,b) for a,b in zip(p,p[1:])]+[arc(p[-1],0)]
  defects=[i for i,(a,b) in enumerate(zip(p,p[1:]),1) if not arc(a,b)]
  if any(not slots[i] for i in defects):continue
  t,q=sum(slots),len(defects); hist[t,q]+=1
  for j in range(q+1):c[t-j]+=(-1)**j*comb(q,j)
 print('slot histogram', sorted(hist.items()), 'coefficient vector',c)
 require(c==[0,-1,18,0,0,0,0,0])

def dp_count(rows):
 n=len(rows); N=1<<n; d=[[0]*n for _ in range(N)]
 for i in range(n):d[1<<i][i]=1
 for s in range(1,N):
  for i in range(n):
   if not s>>i&1:continue
   x=d[s][i]
   for j in range(n):
    if not s>>j&1 and rows[i]>>j&1:d[s|1<<j][j]+=x
 return sum(d[-1])
def blown(k):
 # root clones 0..k-1, old others shifted by k-1
 n=k+6; rows=[0]*n
 old=lambda v: 0 if v<k else v-k+1
 for i in range(n):
  for j in range(n):
   if i==j:continue
   if i<k and j<k: beat=i<j
   else:beat=arc(old(i),old(j))
   if beat:rows[i]|=1<<j
 return rows
cert()
for k in range(1,9):
 h=dp_count(blown(k)); expected=18*2**k-1
 if h!=expected: raise AssertionError((k,h,expected))
 print('k',k,'H',h,'v7_exact_one', h%7==0 and h%49!=0)
for k in range(1,300):
 h=18*2**k-1
 require((h%7==0)==(k%3==1))
 if k%3==1:
  require((h%49==0)==(k%21==7))
print('independent subset DP k=1..8 and congruence diagnostics k=1..299 passed')
