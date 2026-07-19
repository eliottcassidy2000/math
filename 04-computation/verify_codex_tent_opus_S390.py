from fractions import Fraction as F
# Independent check of codex THM-1201 on the cell it names.
V=[1,2,3,4,5,6,7,8,9,10,11,13,24]
L,R=F(1,24),F(1,13)
def fn(x):
    r=x-int(x)
    if r<0: r+=1
    return min(r,1-r)
def g(t): return min(fn(t*v) for v in V)
# 1. no runner has a zero strictly inside the cell?
zeros=[F(k,v) for v in V for k in range(0,v+1)]
inside=[z for z in zeros if L<z<R]
print(f"(1) zeros strictly inside [1/24,1/13]: {inside}  (empty => genuine zero-cell)")
print(f"    g(L) = {g(L)}, g(R) = {g(R)}   (both must be 0)")
# 2. codex's claimed piecewise form
def gc(t):
    if t<=F(1,23): return 24*t-1
    if t<=F(1,14): return t
    return 1-13*t
bad=0
for i in range(1,400):
    t=L+(R-L)*F(i,400)
    if g(t)!=gc(t): bad+=1
print(f"(2) codex's piecewise form vs direct min: mismatches {bad}/399")
# 3. the max H and the two competing area claims
H=max(g(L+(R-L)*F(i,4000)) for i in range(0,4001))
print(f"(3) cell max H = {H} = {float(H):.8f}   (1/14 = {float(F(1,14)):.8f})")
# exact area by integrating the three affine pieces
def area_affine(a,b,f):
    return (f(a)+f(b))*(b-a)/2
A = area_affine(L,F(1,23),gc)+area_affine(F(1,23),F(1,14),gc)+area_affine(F(1,14),R,gc)
tent = H*(R-L)/2
print(f"    exact cell area = {A} = {float(A):.10f}")
print(f"    tent H*L/2      = {tent} = {float(tent):.10f}")
print(f"    => area {'>=' if A>=tent else '<'} tent   "
      f"{'CONCAVE direction (codex right, my THM-1195 backwards)' if A>=tent else 'my direction'}")
# 4. is g concave on the cell? check midpoint inequality on random triples
import random
random.seed(1)
viol=0
for _ in range(300):
    a,b=sorted(L+(R-L)*F(random.randint(1,399),400) for _ in range(2))
    m=(a+b)/2
    if g(m) < (g(a)+g(b))/2 - F(1,10**12): viol+=1
print(f"(4) concavity midpoint test on the cell: violations {viol}/300")
