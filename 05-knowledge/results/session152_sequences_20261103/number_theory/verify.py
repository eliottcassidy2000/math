def require(x):
 if not x: raise AssertionError("verification failed")
#!/usr/bin/env python3
"""Exact arithmetic checks; no external packages."""
from heapq import heappush, heappop
from itertools import combinations_with_replacement
A=(15,17,19,23,25,27,29,33,37,41)
B=(2,3,5,7,11,13,17,19,23,27)
expected=(0,46,17,33,19,50,36,37,23,54,25,41,27,58,29)
def apery(gens):
 m=min(gens); d=[10**9]*m; d[0]=0; q=[(0,0)]
 while q:
  n,r=heappop(q)
  if n != d[r]: continue
  for a in gens:
   s=n+a; t=s%m
   if s<d[t]: d[t]=s; heappush(q,(s,t))
 return tuple(d)
w=apery(A); require(w==expected)
require(max(w)-15==43)
require(sum((x-r)//15 for r,x in enumerate(w))==26)
# Independent exhaustive certificate: all relevant sums need at most three terms.
small={0}
for k in range(1,4):
 small.update(map(sum,combinations_with_replacement(A,k)))
for r,v in enumerate(w):
 require(v in small and not any(n%15==r and n<v for n in small))
# A elements are all odd and < 3*min(A), so a nontrivial representation is impossible.
require(all(a%2 and a<3*min(A) for a in A))
require(set(A)&set(B)=={17,19,23,27})
require(27==3**3)
print('Apéry:',w)
print('Frobenius 43; conductor 44; genus 26; all assertions passed')
