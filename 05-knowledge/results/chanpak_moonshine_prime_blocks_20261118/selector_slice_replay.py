from itertools import permutations
A=[46,28,72,0,109,78,11]
def H(m):
 T=A.copy()
 for j in range(6):
  if m>>j&1:
   T[j]^=64;T[6]^=1<<j
 return sum(all(T[i]>>j&1 for i,j in zip(p,p[1:])) for p in permutations(range(7)))
h={m:H(m) for m in range(64)}
if [h[m] for m in (36,38,20,22)]!=[9,11,13,15]:raise ValueError("wrong window")
def g(m):return h[m^2]-h[m]
board=[m for m in h if not m&2]
coef={1<<j:g(1<<j)-g(0) for j in range(6) if j!=1}
if not all(g(m)==g(0)+sum(c for b,c in coef.items() if m&b) for m in board):
 raise ValueError("not affine")
if coef[16]!=coef[32] or 36^20!=16^32:raise ValueError("not an exchange")
print("Affine gap coefficients:",coef)
print("Literal +2 windows:",[(h[m],h[m^2]) for m in (36,20)])
from itertools import combinations
def hcode(n,m):
 E=list(combinations(range(n),2))
 arc={e if m>>i&1 else e[::-1] for i,e in enumerate(E)}
 return sum(all((i,j) in arc for i,j in zip(p,p[1:])) for p in permutations(range(n)))
if [hcode(5,i) for i in (71,70,77,76)]!=[9,11,13,15]:raise ValueError("5-vertex window")
print("Five-vertex edge-exchange replay passed")
