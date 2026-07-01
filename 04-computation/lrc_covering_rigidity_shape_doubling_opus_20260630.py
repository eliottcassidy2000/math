"""
Chase the covering-rigidity SHAPE. (1) patchability criterion: k doubling-patchable requires 2k not in AP
(k>=n/2, so 2k is NEW) THEN Jacobsthal-gated. (2) is the AP special in int m^2 (smoothest covering)?
"""
from fractions import Fraction
def frac(x): x=x%1.0; return min(x,1-x)
def mult(S,t,n): return sum(1 for v in S if frac(v*t)<=1.0/n+1e-12)
def m2(S,n,G): return sum(mult(S,i/G,n)**2 for i in range(G))/G
def Mexact(S,Q):
    b=Fraction(0)
    for q in range(1,Q+1):
        for a in range(1,q):
            m=min(min((Fraction(s*a,q))%1,1-(Fraction(s*a,q))%1) for s in S)
            if m>b: b=m
    return b
def Mf(S,Q):
    b=0.0
    for q in range(1,Q+1):
        for a in range(1,q):
            m=min(frac(s*a/q) for s in S)
            if m>b: b=m
    return b
print("(1) single-swap doubling patch k->2k: needs 2k NOT in AP (k>=n/2) AND tight. Scan k, check g=2k:")
for n in [8,10,12,14,16]:
    base=list(range(1,n)); res=[]
    for k in base:
        g=2*k
        if g in base: continue  # 2k already a runner -> can't be a NEW patch
        S=[v for v in base if v!=k]+[g]
        tight = (Mf(S,5*n)<=1.0/n+1e-9) and (Mexact(S,10*n)==Fraction(1,n))
        res.append((k,g,tight))
    patchable=[k for (k,g,t) in res if t]
    print(f"  n={n}: k with 2k-new: {[k for k,g,t in res]}; DOUBLING-tight k = {patchable}  (all k>=n/2={n/2}? {all(k>=n/2 for k in patchable)})")
print("   => 2k-new (k>=n/2) is NECESSARY; tight is the further JACOBSTHAL gate (opens rarely: n=8,14).")
print()
print("(2) int m^2 (covering smoothness): AP vs tight retiling vs a NON-tight near-AP:")
for n in [8,14]:
    G=200*n; AP=list(range(1,n))
    ap2=m2(AP,n,G)
    print(f"  n={n}: AP int m^2 = {ap2:.3f}")
    if n==14:
        GW=list(range(1,12))+[13,24]; print(f"    GW (tight)  int m^2 = {m2(GW,n,G):.3f}")
        nt=[v for v in AP if v!=12]+[26]; print(f"    swap 12->26 (non-tight?) int m^2 = {m2(nt,n,G):.3f}, M={float(Mf(nt,10*n)):.4f} vs 1/14={1/14:.4f}")
    if n==8:
        GW=[1,2,3,4,5,7,12]; print(f"    n8 sporadic (tight) int m^2 = {m2(GW,n,G):.3f}")
print("   => compare: is the AP the MIN int m^2 (smoothest)? tight retilings close; non-tight differ.")
