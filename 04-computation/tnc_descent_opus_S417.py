import sympy as sp, numpy as np
from itertools import product
u,v=sp.symbols('u v')
def saddles(R,N):
    R=sp.expand(R); S=sp.expand(u*sp.diff(R,u)-N*R)
    if sp.Poly(S,u).degree()<1: return []
    out=[]
    for r in sp.nroots(sp.Poly(S,u)):
        rv=complex(r)
        if abs(complex(R.subs(u,r)))<1e-12 or abs(rv)<1e-12: continue
        out.append((rv, complex(R.subs(u,r))/rv**N))
    return out
def kfold(R):
    """largest k with R(u) = S(u^k)"""
    P=sp.Poly(sp.expand(R),u); ex=[e[0] for e,c in zip(P.monoms(),P.coeffs()) if c!=0]
    from math import gcd
    g=0
    for e in ex: g=gcd(g,e)
    return g if g>0 else 1
print("ARE ALL DOMINANT-VALUE COLLISIONS CAUSED BY A SYMMETRY R(u) = S(u^k)?")
coll=[]; N=2
for c0,c1,c2,c3 in product([-2,-1,1,2],[-2,-1,0,1,2],[-2,-1,0,1,2],[-2,-1,0,1,2]):
    R=c0+c1*u+c2*u**2+c3*u**3+u**4
    d=saddles(R,N)
    if not d: continue
    rho=max(abs(w) for _,w in d)
    vals=[w for _,w in d if abs(w)>rho*(1-1e-9)]
    if any(abs(vals[i]-vals[j])<1e-8 for i in range(len(vals)) for j in range(i+1,len(vals))):
        coll.append(R)
ks=[kfold(R) for R in coll]
print(f"  collision cases: {len(coll)}   their k-fold symmetries: {sorted(set(ks))}")
print(f"  ALL have k >= 2 (a genuine symmetry)? {all(k>=2 for k in ks)}")
print(f"  examples: {[str(R) for R in coll[:5]]}")
print()
print("THE DESCENT.  If R(u) = S(u^k) then R^m is a polynomial in u^k, so")
print("  [u^{Nm}] R^m = 0 unless k | Nm; and when k | Nm it EQUALS a constant term of a")
print("  SMALLER instance in v = u^k.  Verify on R = u^4 - 2u^2 - 2, N=2, k=2, S(v)=v^2-2v-2:")
R=u**4-2*u**2-2; S=v**2-2*v-2
for m in range(1,7):
    lhs=sp.Poly(sp.expand(R**m),u).coeff_monomial(u**(2*m))
    rhs=sp.Poly(sp.expand(S**m),v).coeff_monomial(v**m)
    print(f"   m={m}: [u^{{2m}}]R^m = {lhs}   [v^m]S^m = {rhs}   {'MATCH' if lhs==rhs else 'differ'}")
print()
print("  => the N=2 instance with k=2 symmetry IS the N'=1 instance Lambda' = v^{-1}S(v),")
print("     and N'=1 (min exponent -1) is PROVED by THM-1530(B) via Lagrange-Burmann.")
print("     So symmetric collisions DESCEND to settled cases.")
