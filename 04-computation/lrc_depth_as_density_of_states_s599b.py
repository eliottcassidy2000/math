"""Is the covering-depth distribution p_k the 'interacting density of states', whose
'free'/independent version is Poisson(lambda=2n*delta)? Generic configs ~ Poisson(2);
additive chains = the anti-Poisson p_0=0 extreme (arithmetic correlation). Grounds the
master-object synthesis. opus-2026-06-03-S599b. Convention: n runners, gap 1/(n+1)."""
from fractions import Fraction as F
from math import exp, factorial
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
def poisson(lam,k): return exp(-lam)*lam**k/factorial(k)
def main():
    print("DEPTH DISTRIBUTION vs POISSON(lambda=2n*delta) — the 'free'/independent baseline")
    print("(generic = near-Poisson [arcs ~ independent]; additive chain = p_0=0 anti-Poisson [correlated])\n")
    cases=[
      ("generic",(1,4,6,9,11),5),("generic",(2,3,7,8,12),5),
      ("generic",(1,5,8,11,13,17),6),("generic",(3,4,9,10,16,19),6),
      ("AP-chain",(1,2,3,4,5),5),("sporadic-chain",(1,3,4,5,9),5),
      ("AP-chain",(1,2,3,4,5,6),6),
    ]
    for tag,V,n in cases:
        d=F(1,n+1); p=depth_dist(V,d); lam=float(2*n*d)
        pk={k:float(p.get(k,F(0))) for k in range(0,5)}
        pois={k:poisson(lam,k) for k in range(0,5)}
        p0=float(p.get(0,F(0)))
        # total variation distance to Poisson over k=0..n
        kk=range(0,n+1); tv=0.5*sum(abs(float(p.get(k,F(0)))-poisson(lam,k)) for k in kk)
        print(f"{tag:16s} V={V} lam={lam:.3f}: p_0={p0:.3f} (Pois e^-lam={exp(-lam):.3f}); TV(p,Pois)={tv:.3f}")
        print(f"   p_k   ={[round(pk[k],3) for k in range(5)]}")
        print(f"   Pois_k={[round(pois[k],3) for k in range(5)]}")
    print("\nTakeaway: generic configs sit near Poisson(2) with p_0~e^-2~0.135>0 (INDEPENDENCE => lonely);")
    print("additive chains force p_0=0 (CORRELATION kills the ground state). LRC = correlation can")
    print("empty the bulk p_0 but never the measure-0 witness set. The worry-set = the large-deviation residual.")
if __name__=='__main__': main()
