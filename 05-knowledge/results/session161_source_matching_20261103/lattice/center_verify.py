from itertools import product
from math import gcd
from collections import Counter

def test(m,N):
 S=[[0]*m for _ in range(m)]
 for i in range(m):
  for j in range(i+2,m):S[i][j]=2**(j-i-2);S[j][i]=-S[i][j]
 b=[2**max(0,m-2-j) for j in range(m)]
 cols=[[S[i][j]%N for i in range(m)] for j in range(m)]
 counts=Counter()
 for x in product(range(N),repeat=m):
  y=tuple(sum(cols[j][i]*x[j] for j in range(m))%N for i in range(m))
  counts[y]+=1
 checks=[]
 for t in range(N):
  actual=counts[tuple((t*v)%N for v in b)]
  expected=(gcd((m//2)-1,N)**2 if m%2==0 and m>=4 and t%gcd((m//2)-1,N)==0 else (N if m%2==1 and t==0 else N**m if m==2 and t==0 else 0))
  if actual != expected: raise AssertionError((m,N,t,actual,expected))
  checks.append(actual)
 return checks
for m,Ns in ((2,(2,3,4)),(3,(2,3,4)),(4,(2,3,4,5)),(5,(2,3,4)),(6,(2,3,4)),(7,(2,3)),(8,(2,3))):
 for N in Ns:print(f'm={m} N={N} counts(t=0..N-1)={test(m,N)}')
