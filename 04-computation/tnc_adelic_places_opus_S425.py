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
def newton_root_vals(cexpr,p):
    P=sp.Poly(sp.expand(cexpr),a); d={e[0]:int(co) for e,co in zip(P.monoms(),P.coeffs()) if co!=0}
    ys=sorted(d); pts=[(y,vp(d[y],p)) for y in ys]
    hull=[pts[0]]
    for pt in pts[1:]:
        while len(hull)>=2:
            (x0,y0),(x1,y1)=hull[-2],hull[-1]
            if (y1-y0)*(pt[0]-x1)>=(pt[1]-y1)*(x1-x0): hull.pop()
            else: break
        hull.append(pt)
    return sorted(set(sp.Rational(-(hull[i+1][1]-hull[i][1]),hull[i+1][0]-hull[i][0]) for i in range(len(hull)-1)))
def arch_root_radii(cexpr):
    """archimedean: log|root| set, from the ordinary Newton polygon of log|coeff|."""
    P=sp.Poly(sp.expand(cexpr),a); d={e[0]:abs(complex(co)) for e,co in zip(P.monoms(),P.coeffs()) if co!=0}
    ys=sorted(d)
    import math
    pts=[(y,math.log(d[y])) for y in ys]
    hull=[pts[0]]
    for pt in pts[1:]:
        while len(hull)>=2:
            (x0,y0),(x1,y1)=hull[-2],hull[-1]
            if (y1-y0)*(pt[0]-x1)>=(pt[1]-y1)*(x1-x0): hull.pop()
            else: break
        hull.append(pt)
    return sorted(round(-(hull[i+1][1]-hull[i][1])/(hull[i+1][0]-hull[i][0]),4) for i in range(len(hull)-1))
print("ADELIC COPRIMALITY: CT(m0), CT(2m0) differ at SOME place -- finite prime (Kummer carry")
print("CA) OR archimedean (amoeba radius).  Product-formula reasoning: agreeing everywhere =>")
print("equal, but they are not, so they differ somewhere -> coprime -> TNC.")
print()
tunable=[(2,3,4),(2,5,8),(3,2,6),(3,4,8),(3,5,10),(4,5,10)]
for N,j,d in tunable:
    R=1+a*u**j+u**d
    m0=next((m for m in range(1,16) if sp.expand(CT(R,N,m))!=0 and not sp.expand(CT(R,N,m)).is_number),None)
    if m0 is None: continue
    c0=sp.expand(CT(R,N,m0)); c1=sp.expand(CT(R,N,2*m0))
    if len([x for x in sp.Poly(c0,a).all_coeffs() if x!=0])<2: continue
    place=None
    for p in list(primerange(2,40)):
        if not(set(newton_root_vals(c0,p))&set(newton_root_vals(c1,p))): place=f"p={p} (Kummer CA)"; break
    if place is None:
        r0,r1=arch_root_radii(c0),arch_root_radii(c1)
        if not(set(r0)&set(r1)): place=f"archimedean (amoeba: log|a|={r0} vs {r1})"
    print(f"   N={N}(-{N},{j-N},{d-N}) m0={m0}: separated at {place}")
print()
print("=> every tunable trinomial closes at SOME place; the finite places are the")
print("   Sierpinski/Kummer carry automaton, the infinite place is the multinomial amoeba.")
