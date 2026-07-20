# {-1,0,1} stratum, NON-CONSTANT coeffs a(r),b(r),c(r) in C[r].  P = a(r)W + b(r) + c(r)Z.
# E[P^m] = E_r[ L_m(alpha,beta) ],  alpha = a*c*r, beta = b,  E_r[poly] = sum_j coeff_j * j!,
#   L_m(alpha,beta) = sum_k m!/(k!^2 (m-2k)!) alpha^k beta^{m-2k}   (= CT_w(P^m)).
# NC2 says: E[P^m]=0 for all m  <=>  b==0 AND a*c==0 (one-sided).  Search for any violation.
from fractions import Fraction as Fr
from math import factorial
from itertools import product
def pmul(p,q):
    r={}
    for i,u in enumerate(p):
        if u==0: continue
        for j,v in enumerate(q):
            if v==0: continue
            r[i+j]=r.get(i+j,Fr(0))+u*v
    d=max(r)+1 if r else 1
    out=[Fr(0)]*d
    for k,v in r.items(): out[k]=v
    return out
def padd(p,q):
    d=max(len(p),len(q)); out=[Fr(0)]*d
    for i,u in enumerate(p): out[i]+=u
    for i,u in enumerate(q): out[i]+=u
    return out
def pscale(p,s): return [u*s for u in p]
def Er(poly): return sum(c*factorial(j) for j,c in enumerate(poly))  # Gamma(1) moments
def poly_is_zero(p): return all(u==0 for u in p)
def Lm_Er(a,b,c,M):
    # alpha = a*c*r  (shift a*c up by one r-power), beta=b
    ac=pmul(a,c); alpha=[Fr(0)]+ac   # multiply by r
    beta=b
    # precompute powers
    res=[]
    for m in range(1,M+1):
        tot=[Fr(0)]
        for k in range(0,m//2+1):
            # coeff m!/(k!^2 (m-2k)!)
            coef=Fr(factorial(m),factorial(k)**2*factorial(m-2*k))
            term=[Fr(1)]
            for _ in range(k): term=pmul(term,alpha)
            for _ in range(m-2*k): term=pmul(term,beta)
            tot=padd(tot,pscale(term,coef))
        res.append(Er(tot))
    return res
# sanity: constant a=c=1,b=1 -> E[P^m] = s^m He_m(1/s), s=sqrt(-2). Check nonzero, matches earlier.
print("sanity const a=b=c=1, E[P^m] m=1..6:", [str(x) for x in Lm_Er([Fr(1)],[Fr(1)],[Fr(1)],6)])
# SEARCH: a,b,c linear (deg<=1), integer coeffs in [-2,2], find E[P^m]=0 for m=1..M0 with NOT one-sided
vals=[Fr(x) for x in range(-2,3)]
M0=10
checked=0; nullcone_nononesided=[]
polys1=[list(t) for t in product(vals,repeat=2)]  # a=a0+a1 r
for a in polys1:
    for c in polys1:
        for b in polys1:
            if poly_is_zero(a) and poly_is_zero(b) and poly_is_zero(c): continue
            checked+=1
            E=Lm_Er(a,b,c,M0)
            if all(x==0 for x in E):
                ac=pmul(a,c)
                onesided = poly_is_zero(b) and poly_is_zero(ac)
                if not onesided:
                    nullcone_nononesided.append((a,b,c,E))
print(f"\nsearched {checked} triples (a,b,c linear, coeffs in [-2,2]); m=1..{M0}")
print(f"nullcone members that are NOT one-sided (NC2 COUNTEREXAMPLES): {len(nullcone_nononesided)}")
for a,b,c,E in nullcone_nononesided[:5]:
    print("   a,b,c=",[str(x) for x in a],[str(x) for x in b],[str(x) for x in c])
print("=> NC2 holds on this box:", len(nullcone_nononesided)==0)
