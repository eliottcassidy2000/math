# Confirm the NEAR-EQUAL COLLAPSE: bad-kick-set union is ~one killer for near-equal, ~j for spread.
from fractions import Fraction as F
from math import gcd
def norm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return min(r,1-r)
def minnorm(fam,t): return min(norm(F(v)*t) for v in fam)
def bad_fraction(P,K,N=4000):
    # fraction of kick range [-smax,smax] where SOME killer is dangerous (||k(t0+s)||<1/14)
    # t0 = P-maximizer
    best=F(0); t0=None
    for q in range(2,2*max(P)+2):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            v=minnorm(P,F(a,q))
            if v>best: best=v; t0=F(a,q)
    smax=(best-F(1,14))/max(P)
    bad=0; badper=[0]*len(K)
    for i in range(N):
        s=smax*F(2*i-N,N); t=t0+s
        anybad=False
        for idx,k in enumerate(K):
            if norm(F(k)*t)<F(1,14): badper[idx]+=1; anybad=True
        if anybad: bad+=1
    return best,bad/N,[b/N for b in badper]

P=list(range(1,11))
for K in [[257,258,263],[300,301,302],[200,300,450]]:  # near-equal, tight, spread
    mp,uf,per=bad_fraction(P,K)
    print(f"K={K}: M(P)={float(mp):.4f} bad-union frac={uf:.3f}, per-killer≈{[round(x,3) for x in per]}, sum-per={sum(per):.3f}")
    print(f"   => union {uf:.3f} vs sum {sum(per):.3f}: {'COLLAPSED (overlap)' if uf<0.8*sum(per) else 'additive'}; good s exists: {uf<1.0}")
