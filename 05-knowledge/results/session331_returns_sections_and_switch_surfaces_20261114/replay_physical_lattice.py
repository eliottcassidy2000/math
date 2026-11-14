from itertools import permutations
from math import comb, gcd
from sympy import Poly, Symbol, expand
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent /
                      "session328_fermat_gluing_and_phase_repair_20261114"))
from profile_tools import strong
from replay_matching_sections import coefficients, cyclic, replace, trans

z=Symbol('z')
def L(f,p):
    return expand(z**(p+1)*f+z*((z+1)**p-z**p)*(f.subs(z,z+1)-f)
              +z*((z-1)**p-z**p)*(f.subs(z,z-1)-f))

def residual(f,d):
    while (F:=Poly(f,z)).degree()>d:
        f=expand(f-F.LC()*L(z**(F.degree()-2),1))
    return int(f.coeff(z,d))

def q(A,k):
    return sum(comb(k,i)*x for i,x in enumerate(coefficients(A,k),1))

for U in ((0,),(2,4,1),(2,4,9,3),(6,12,24,17,3)):
    n=len(U)
    # Newton basis interpolation on optional slots.
    term=1
    f=Poly(0,z).as_expr()
    from math import factorial
    B=coefficients(U,n)
    for i in range(1,n+1):
        term*=z-(i-1)
        f+=term*B[i-1]/factorial(i)
    f=expand(f)
    for p in (1,2):
        G=cyclic(U,trans(p),(0,))
        assert strong(G)
        g=L(f,p)
        assert all(g.subs(z,k)==q(G,k) for k in range(len(G)+1))

assert residual(L(z**2,2),1)==-4
assert [residual(L(z**i,2),2) for i in (3,5)]==[-14,114]
assert residual((L(z**5,2)-L(z**6,1))/2,2)==57
C5=(6,12,24,17,3)
def ham_cycles(A):
    n=len(A)
    return sum(all(A[x]>>y&1 for x,y in zip((0,)+p,p+(0,))) for p in permutations(range(1,n)))
assert ham_cycles(C5)==2
assert [ham_cycles(replace(C5,0,trans(n-4)))%2 for n in (5,7,9)]==[0,0,0]
for n in (5,6):
    R=[(2,4,9,19,7),(2,20,9,3,13),(2,12,9,17,7)] if n==5 else [
        (2,4,17,7,43,15),(2,16,3,7,45,15),(2,4,41,19,7,27)]
    assert all(strong(A) for A in R)
    (h,p),(i,o),(j,r)=[(q(U,1),q(U,2)) for U in R]
    assert abs((i-h)*(r-p)-(j-h)*(o-p))==(24 if n==5 else 48)
print('literal cycle-wrap polynomials, residuals, parity controls and primitive bases: PASS')
