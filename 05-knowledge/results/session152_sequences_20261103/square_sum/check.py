def require(x):
 if not x: raise AssertionError("square-sum check failed")
from math import isqrt
from itertools import combinations
A={15,17,19,23,25,27,29,33,37,41}
B={2,3,5,7,11,13,17,19,23,27}
def square(x): return isqrt(x)**2==x
def edges(S): return {(x,y) for x,y in combinations(sorted(S),2) if square(x+y)}
require(edges(A)=={(17,19),(23,41),(27,37)})
require(edges(B)=={(2,7),(2,23),(3,13),(5,11),(13,23),(17,19)})
for n in range(1,1001):
 b=isqrt(2*n-1)-isqrt(n)
 require(b==sum(square(n+x) for x in range(1,n)))
 require(b==sum(k*k//2+1<=n<=k*k-1 for k in range(2,isqrt(2*n)+2)))
print('PASS induced edges and birth identities through n=1000')
