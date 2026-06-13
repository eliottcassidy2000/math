"""The FTA n+1↔n map applied to the repo: the covering-depth distribution {p_k} IS a polynomial
P(z)=Σ p_k z^k (the PGF). n+1 coefficients (p_0..p_m, the constant p_0 = the lonely measure) ↔ m
roots (with multiplicity). WORRY-SET ⟺ p_0=0 ⟺ z=0 is a root (the constant term vanishes, n+1→n
collapse). Independent/'free' dangers ⟹ all roots coincide (multiplicity m); correlation (worry-set)
⟹ roots SPREAD. Compute roots & multiplicities for AP (worry-set) vs loose. opus-2026-06-06-S699l."""
from fractions import Fraction as F
try:
    import numpy as np; HAVE=True
except Exception: HAVE=False
def dist(x): x%=1; return min(x,1-x)
def depth_dist(V,delta):
    eps=set([F(0)])
    for v in V:
        for k in range(v+1):
            for s in(1,-1): eps.add((F(k)+s*delta)/v % 1)
    pts=sorted(eps); p={}; L=len(pts)
    for i in range(L):
        a=pts[i]; b=pts[(i+1)%L]; ln=(b-a) if b>a else (b-a+1); mid=(a+ln/2)%1
        d=sum(1 for v in V if dist(v*mid)<delta); p[d]=p.get(d,F(0))+ln
    return p
def main():
    print("Covering-depth POLYNOMIAL P(z)=Σ p_k z^k : n+1 coeffs (p_0=lonely measure) ↔ m roots.")
    print("WORRY-SET ⟺ p_0=0 ⟺ z=0 is a root (constant term vanishes).\n")
    configs=[("AP n=5 (worry-set)",(1,2,3,4),5),("chain (worry)",(1,3,4,7),5),
             ("loose",(1,4,6,9),5),
             ("AP n=6 (worry)",(1,2,3,4,5),6),("loose n=6",(2,3,7,8,12),6)]
    for tag,V,n in configs:
        d=F(1,n); p=depth_dist(V,d); m=max(p)
        coeffs=[float(p.get(k,F(0))) for k in range(m+1)]  # p_0..p_m (ascending)
        p0=p.get(0,F(0))
        print(f"{tag}: V={V}, gap 1/{n}; p_0={p0}={'WORRY (z=0 root)' if p0==0 else 'loose'}")
        print(f"   coeffs (p_0..p_{m}) = {[round(c,4) for c in coeffs]}; Σp_k={round(sum(coeffs),4)}")
        if HAVE:
            # numpy.roots wants highest-degree first
            hi=coeffs[::-1]
            # strip leading zeros
            while len(hi)>1 and abs(hi[0])<1e-12: hi=hi[1:]
            r=np.roots(hi)
            rr=sorted([complex(round(x.real,3),round(x.imag,3)) for x in r], key=lambda z:(z.real,z.imag))
            mods=sorted(round(abs(x),3) for x in r)
            print(f"   roots = {rr}")
            print(f"   |roots| = {mods}  (all-equal ⟺ 'free/independent'; spread ⟺ correlated/worry)")
        print()
    print("So: the constant term p_0 (the lonely measure / the OBSERVER's safe mass) is the '+1';")
    print("its vanishing = the worry-set = z=0 root = the n+1→n degree collapse. Root multiplicity")
    print("= degeneracy: free dangers ⟹ a single repeated root; correlation spreads the roots.")
if __name__=='__main__': main()
