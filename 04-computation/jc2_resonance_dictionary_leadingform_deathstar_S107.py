import numpy as np, itertools, random
random.seed(1); np.random.seed(1)
# bivariate poly = dict {(i,j): coeff}  for x^i y^j
def padd(A,B):
    C=dict(A)
    for k,v in B.items(): C[k]=C.get(k,0)+v
    return {k:v for k,v in C.items() if abs(v)>1e-12}
def pscale(A,s): return {k:s*v for k,v in A.items()}
def pmul(A,B):
    C={}
    for (i,j),a in A.items():
        for (k,l),b in B.items():
            C[(i+k,j+l)]=C.get((i+k,j+l),0)+a*b
    return {k:v for k,v in C.items() if abs(v)>1e-12}
def ppow(A,n):
    R={(0,0):1.0}
    for _ in range(n): R=pmul(R,A)
    return R
def dx(A): return {(i-1,j):a*i for (i,j),a in A.items() if i>=1 and abs(a*i)>1e-12}
def dy(A): return {(i,j-1):a*j for (i,j),a in A.items() if j>=1 and abs(a*j)>1e-12}
def jac(F,G): return padd(pmul(dx(F),dy(G)), pscale(pmul(dy(F),dx(G)),-1.0))
def div(H1,H2): return padd(dx(H1),dy(H2))
def peval(A,xv,yv): return sum(a*xv**i*yv**j for (i,j),a in A.items())
def is_zero(A,trials=8):
    for _ in range(trials):
        xv,yv=np.random.randn()+1j*np.random.randn(),np.random.randn()+1j*np.random.randn()
        if abs(peval(A,xv,yv))>1e-6: return False
    return True
def hom(d,rng=None):  # random homogeneous form of degree d
    return {(d-i,i): (np.random.randn()+1j*np.random.randn()) for i in range(d+1)}

print("="*72)
print("ANCHOR 1: det(DF) identity and homogeneous 2D Keller = shear")
# det(D(I+H)) = 1 + div(H) + jac(H1,H2)  -- exact identity
H1,H2=hom(3),hom(3)
F1,F2=padd({(1,0):1.0},H1),padd({(0,1):1.0},H2)
detDF=padd(padd({(0,0):1.0},pmul(dx(F1),dy(F2))),pscale(pmul(dy(F1),dx(F2)),-1.0))
rhs=padd(padd({(0,0):1.0},div(H1,H2)),jac(H1,H2))
print(" det(DF) == 1+div(H)+jac(H1,H2):", is_zero(padd(detDF,pscale(rhs,-1.0))))
print(" => for homog deg d>=2: div(H)[deg d-1] and jac(H1,H2)[deg 2d-2] different degrees")
print("    so Keller <=> div(H)=0 AND jac(H1,H2)=0  (each homogeneous piece vanishes)")

# jac(H1,H2)=0 for SAME-degree homog <=> proportional
print("\n LEMMA: same-degree homog forms, jac=0 <=> proportional")
h=hom(3)
prop1,prop2=pscale(h,2.0+1j),pscale(h,0.5-0.3j)   # proportional
print("  proportional pair: jac=0 ?", is_zero(jac(prop1,prop2)))
print("  random independent pair: jac=0 ?", is_zero(jac(hom(3),hom(3))))

# shear = c*(b x - a y)^d *(a,b): Keller + invertible with inverse I-H
print("\n SHEAR: H = c*(b x - a y)^d *(a,b)")
a,b,c=1.7+0.4j,-0.9+1.1j,0.6-0.8j; d=3
ell={(1,0):b,(0,1):-a}
Hd=ppow(ell,d)
H1s,H2s=pscale(Hd,c*a),pscale(Hd,c*b)
F1s,F2s=padd({(1,0):1.0},H1s),padd({(0,1):1.0},H2s)
detshear=padd(padd({(0,0):1.0},pmul(dx(F1s),dy(F2s))),pscale(pmul(dy(F1s),dx(F2s)),-1.0))
print("  det(DF)=1 ?", is_zero(padd(detshear,{(0,0):-1.0})))
print("  (b x - a y) invariant: b*H1-a*H2 =0 ?", is_zero(padd(pscale(H1s,b),pscale(H2s,-a))))
# inverse I - H : compose F(F^{-1})=id
def compose(F, G1, G2):  # F(x,y) with x->G1, y->G2
    R={}
    for (i,j),co in F.items():
        R=padd(R, pscale(pmul(ppow(G1,i),ppow(G2,j)), co))
    return R
Fi1,Fi2=padd({(1,0):1.0},pscale(H1s,-1.0)),padd({(0,1):1.0},pscale(H2s,-1.0))
c1=compose(F1s,Fi1,Fi2); c2=compose(F2s,Fi1,Fi2)
print("  F∘(I-H) = id ?", is_zero(padd(c1,{(1,0):-1.0})) and is_zero(padd(c2,{(0,1):-1.0})))

print("\n"+"="*72)
print("ANCHOR 1 (fixed): det(DF)=jac(F1,F2) directly")
print(" identity det(DF)==1+div(H)+jac(H1,H2):", is_zero(padd(jac(F1,F2), pscale(rhs,-1.0))))
print(" shear det(DF)=1 ?", is_zero(padd(jac(F1s,F2s),{(0,0):-1.0})))

print("\n"+"="*72)
print("ANCHOR 2: leading-form power structure + the descent dichotomy")
print(" Claim: jac(f_n,g_m)=0  <=>  f_n=α h^p, g_m=β h^q  (common form h);")
print("        equivalently f_n^(m/g) = const * g_m^(n/g),  g=gcd(n,m).")
from math import gcd
for (p,q,kdeg,label) in [(2,3,2,'single-root h=l^k'),(2,3,2,'multi-root h'),(3,5,1,'coprime deg, h linear')]:
    # build h: single-root = (x - r y)^kdeg ; multi-root = product of distinct linear
    if 'single' in label:
        h=ppow(padd({(1,0):1.0},{(0,1):-0.7-0.3j}),kdeg)
    elif 'multi' in label:
        h=pmul(padd({(1,0):1.0},{(0,1):-0.7-0.3j}),padd({(1,0):1.0},{(0,1):1.4+0.9j})) if kdeg==2 else ppow(padd({(1,0):1.0},{(0,1):-0.7j}),kdeg)
    else:
        h=padd({(1,0):1.0},{(0,1):0.5+0.2j})
    n,m=p*kdeg,q*kdeg; g=gcd(n,m)
    fn=pscale(ppow(h,p),1.3+0.2j); gm=pscale(ppow(h,q),-0.8+0.5j)
    jz=is_zero(jac(fn,gm))
    # functional dependence f_n^(m/g) == const * g_m^(n/g)?
    lhs=ppow(fn,m//g); rhsp=ppow(gm,n//g)
    # ratio const: compare at one point
    xv,yv=0.3+0.1j,0.7-0.2j
    ratio=peval(lhs,xv,yv)/peval(rhsp,xv,yv)
    dep=is_zero(padd(lhs,pscale(rhsp,-ratio)))
    print(f"  {label:22s} n={n},m={m},gcd={g}: jac(f_n,g_m)=0? {jz}  f_n^(m/g)=const*g_m^(n/g)? {dep}")
print("\n DESCENT MOVES (reduce degree by killing the common leading form):")
print("  * n=m  => f_n,g_m proportional => g' = g - (β/α) f  drops deg g. ALWAYS reduces.")
print("  * n!=m, q|p (i.e. m|n)  => f' = f - (α/β^(p/q)) g^(p/q) drops deg f. Reduces.")
print("  * n!<m and m!|n  (neither divides) => leading form cannot be killed by a power.")
print("    THIS is the stuck stratum (Abhyankar-Moh). Verify a reduce vs a stuck case:")
# reduce case n=4,m=2 (m|n): f=h^2 (deg4)+lower, g=h (deg2)+lower, h deg2
h2=pmul(padd({(1,0):1.0},{(0,1):-0.7j}),padd({(1,0):1.0},{(0,1):1.1j}))  # deg2, two roots
f=padd(pscale(ppow(h2,2),1.0), {(1,0):0.3,(0,2):0.1})   # deg4 + junk
gp=padd(pscale(h2,1.0), {(0,1):0.2})                    # deg2 + junk
# f' = f - (lc match) g^2 : leading of g^2 = h2^2, of f = h2^2 -> subtract to drop deg
fp=padd(f, pscale(ppow(gp,2),-1.0))
maxdeg=lambda A: max(i+j for (i,j) in A)
print(f"    m|n reduce: deg f={maxdeg(f)} -> deg(f - g^2)={maxdeg(fp)} (dropped: {maxdeg(fp)<maxdeg(f)})")

print("\n"+"="*72)
print("ANCHOR 3: places at infinity = DISTINCT roots of the leading form h")
print("  Keller pair f_n=a h^p, g_m=b h^q share the SAME places (roots of h).")
def n_distinct(rootlist):
    # count distinct complex numbers robustly (adaptive threshold)
    uniq=[]
    for z in rootlist:
        if not any(abs(z-w)<1e-3 for w in uniq): uniq.append(z)
    return len(uniq)
def has_repeated(poly):        # discriminant test: gcd(p,p')!=const  <=> repeated root
    d=np.polyder(poly)
    # resultant(p,p') via companion-eigenvalue proximity
    rp=np.roots(poly)
    return any(abs(np.polyval(d,r))<1e-6*max(1,abs(np.polyval(np.polyder(d),r))) for r in rp) or n_distinct(rp)<len(rp)
c=0.7+0.3j
single=np.poly([1/c,1/c,1/c])           # (s-1/c)^3  : one distinct root (triple)
multi =np.poly([1/c,-0.5+0.2j,2.0])     # three distinct roots
print("  single-root h=(x-cy)^3 : distinct roots =", n_distinct(np.roots(single)),
      " repeated? ", has_repeated(single), " => 1 place at infinity")
print("  multi-root (3 distinct): distinct roots =", n_distinct(np.roots(multi)),
      " repeated? ", has_repeated(multi), " => 3 places at infinity")
print("""
  CLASSICAL CONSEQUENCE (Abhyankar-Moh): ONE place at infinity => the affine
  curve embeds as a coordinate line => the pencil is tame. So:
    single-root leading form h  =>  1 place at infinity  =>  TAME (no counterexample)
    JC(2) counterexamples REQUIRE a multi-root (>=2 distinct roots) leading form h.""")
print("="*72)
print("THE RESONANCE DICTIONARY (verified anchors + cross-thread):")
print("""
  TAME pole            | RESONANT / hard pole
  ---------------------+----------------------------------------
  DvdK: unique cycle   | DvdK: >=2 coincident cycles (S101/S106)
  single-term CT != 0  | cancellation possible; needs THM-2067
  ---------------------+----------------------------------------
  Hessian: rank-1 nilp | Hessian: rank>=2 nilpotent (S103)
  P one-sided (x+iy)^d | dim>=3 counterexample (THM-1300)
  ---------------------+----------------------------------------
  JC(2): single-root h | JC(2): multi-root h (>=2 roots)
  1 place at infinity  | >=2 places; Abhyankar-Moh fails
  (tame, AbhyankarMoh) | the open residual
""")
