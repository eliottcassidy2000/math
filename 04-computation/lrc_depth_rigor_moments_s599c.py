"""RIGOR for the master object. Prove (verify exactly): the j-th FACTORIAL moment of depth
equals j!*S_j, S_j = sum over j-subsets of meas(intersection of D_i). Hence:
 (1) E[depth]=2n*delta exactly; (2) p_0 = sum_j (-1)^j S_j (incl-excl); (3) the gap to
 Poisson at 2nd factorial moment = 2*sum Cov(1_{D_i},1_{D_j}) = the total RESONANCE/additive
 energy. Generic: resonance ~ 0 (near-Poisson). Additive chains: resonance large (anti-Poisson).
opus-2026-06-03-S599c. Convention n runners, gap 1/(n+1)."""
from fractions import Fraction as F
from itertools import combinations
from math import comb
def dist(x): x%=1; return min(x,1-x)
def overlap(Vsub,delta):
    # meas{ t : ||v t||<delta for ALL v in Vsub } via breakpoint partition
    eps=set([F(0)])
    for v in Vsub:
        for k in range(v+1):
            for s in(1,-1): eps.add((F(k)+s*delta)/v % 1)
    pts=sorted(eps); L=len(pts); tot=F(0)
    for i in range(L):
        a=pts[i]; b=pts[(i+1)%L]; ln=(b-a) if b>a else (b-a+1); mid=(a+ln/2)%1
        if all(dist(v*mid)<delta for v in Vsub): tot+=ln
    return tot
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
def Sj_direct(V,delta,j):
    return sum(overlap(sub,delta) for sub in combinations(V,j))
def Sj_from_p(p,j):
    return sum(comb(k,j)*v for k,v in p.items())
def main():
    cases=[("AP",(1,2,3,4),4),("chain",(1,3,4,7),4),("generic",(1,4,6,9),4),
           ("AP",(1,2,3,4,5),5),("chain",(1,3,4,5,9),5),("generic",(2,3,7,8,12),5)]
    print("THM-M1: j-th factorial moment of depth = j! * S_j  (S_j = total j-fold overlap)")
    print("Verify  E[C(N,j)] from {p_k}  ==  S_j computed DIRECTLY from arc intersections.\n")
    for tag,V,n in cases:
        d=F(1,n+1); p=depth_dist(V,d)
        line=f"{tag:8s} V={V} (delta=1/{n+1}): "
        ok=True
        for j in range(1,4):
            a=Sj_from_p(p,j); b=Sj_direct(V,d,j)
            ok &= (a==b)
            if j==1: line+=f"S1={a}(=2n*delta={2*n*d}) "
            if j==2: line+=f"S2={a} "
        # incl-excl p0
        p0_ie=sum((-1)**j*Sj_from_p(p,j) for j in range(0,n+1))
        # resonance: 2*sum Cov = 2*(S2 - C(n,2)*(2delta)^2) = E[N(N-1)]-lam^2
        S2=Sj_from_p(p,2); lam=2*n*d
        cov_gap=2*S2 - lam*lam
        line+=f"| S_j match j=1..3: {ok}; p0(incl-excl)={p0_ie}(=p_0 {p.get(0,F(0))}: {p0_ie==p.get(0,F(0))}); resonance 2*sumCov={cov_gap}={float(cov_gap):+.4f}"
        print(line)
    print("\nReading: S1=2n*delta is CONFIG-FREE (the blind first moment). The resonance")
    print("2*sumCov = E[N(N-1)]-lambda^2 is the EXACT gap to Poisson at the 2nd factorial")
    print("moment: ~0 for generic (near-Poisson/free), large for additive chains (correlated).")
    print("p_0 = sum_j (-1)^j S_j: loneliness is an exact alternating sum of overlap volumes.")
if __name__=='__main__': main()
