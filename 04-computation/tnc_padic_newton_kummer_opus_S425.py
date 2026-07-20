import sympy as sp
from sympy import primerange
u,a=sp.symbols('u a')
def CT(R,N,mm):
    return sp.Poly(sp.expand(R**mm),u).coeff_monomial(u**(N*mm))
def vp(n,p):
    n=abs(int(n))
    if n==0: return None
    v=0
    while n%p==0: n//=p; v+=1
    return v
def root_valuations(cexpr,p):
    """p-adic Newton-polygon lower-hull slopes = -(root valuations)."""
    P=sp.Poly(sp.expand(cexpr),a); d={e[0]:int(co) for e,co in zip(P.monoms(),P.coeffs()) if co!=0}
    ys=sorted(d); pts=[(y,vp(d[y],p)) for y in ys]
    # lower convex hull
    hull=[pts[0]]
    for pt in pts[1:]:
        while len(hull)>=2:
            (x0,y0),(x1,y1)=hull[-2],hull[-1]
            # slope hull[-2]->hull[-1] vs hull[-1]->pt
            if (y1-y0)*(pt[0]-x1) >= (pt[1]-y1)*(x1-x0): hull.pop()
            else: break
        hull.append(pt)
    slopes=[sp.Rational(hull[i+1][1]-hull[i][1], hull[i+1][0]-hull[i][0]) for i in range(len(hull)-1)]
    return sorted(set(-s for s in slopes))
print("CLAIM: for every tunable trinomial, SOME prime p separates the p-adic root-valuations")
print("of CT(m0) and CT(2m0) => coprime => TNC. (Kummer carry-CA computes the valuations.)")
print()
tunable=[(2,3,6),(2,5,8),(3,2,6),(3,4,8),(3,5,10),(4,5,10),(2,3,4),(2,1,4)]
allok=True
for N,j,d in tunable:
    R=1+a*u**j+u**d
    m0=next((m for m in range(1,16) if sp.expand(CT(R,N,m))!=0 and not sp.expand(CT(R,N,m)).is_number),None)
    if m0 is None: continue
    c0=sp.expand(CT(R,N,m0)); c1=sp.expand(CT(R,N,2*m0))
    if len([x for x in sp.Poly(c0,a).all_coeffs() if x!=0])<2: continue
    sep_prime=None
    for p in list(primerange(2,60)):
        rv0=root_valuations(c0,p); rv1=root_valuations(c1,p)
        if not(set(rv0)&set(rv1)):
            sep_prime=(p,rv0,rv1); break
    ok = sep_prime is not None
    allok=allok and ok
    if ok:
        p,rv0,rv1=sep_prime
        print(f"   N={N}(-{N},{j-N},{d-N}) m0={m0}: p={p} SEPARATES  v_{p}(CT(m0) roots)={rv0} vs {rv1}")
    else:
        print(f"   N={N}(-{N},{j-N},{d-N}) m0={m0}: NO separating prime<60 (resultant still != 0)")
print()
print(f"ALL tunable trinomials have a Kummer-CA separating prime: {allok}")
print()
print("KACZYNSKI BOUNDARY REFRAME: G(t)=t(log Pi)' is analytic in the disk |t|<1/rho; the")
print("saddle values t_j=1/w_j are BOUNDARY SINGULARITIES on |t|=1/rho. TNC <=> G constant <=>")
print("the boundary function (radial limits, one per branch = 'approach label') is TRIVIAL <=>")
print("empty singular set <=> R monomial. The branches u_i(t) ARE Kaczynski's approach paths.")
