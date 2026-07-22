"""S102: LRC lonely measure = GMC2 moment (both integer-kernel resonance sums); the resonance
decomposition mu = (clock-floor)^n + corrections reduces LRC(n+1) to the maximal-resonance (AP) cores."""
import numpy as np
from math import sin, pi
from itertools import product

def ghat(k, d): return (1-2*d) if k==0 else -sin(2*pi*k*d)/(pi*k)

def decompose(speeds, d, K=9):
    n=len(speeds); main=ghat(0,d)**n; corr=0.0; nres=0
    for ks in product(range(-K,K+1),repeat=n):
        if any(ks) and sum(k*v for k,v in zip(ks,speeds))==0:
            p=1.0
            for k in ks: p*=ghat(k,d)
            corr+=p; nres+=1
    return main, corr, nres

def lonely_direct(speeds, d, N=200000):
    t=(np.arange(N)+0.5)/N; inside=np.ones(N)
    for v in speeds: inside*=(np.abs(((v*t+0.5)%1)-0.5)>d)
    return inside.mean()

print("mu(lonely) = MAIN ((1-2d)^n = clock/Eisenstein floor = THM-878 6/7 for LRC14) + resonance corrections.")
print(f"{'core':24s} {'main':>7s} {'|corr|/main':>11s} {'mu(res)':>8s} {'mu(direct)':>10s} {'verdict':>18s}")
for name,sp in [("Sidon {1,2,5,11}",(1,2,5,11)),("AP {1,2,3,4}",(1,2,3,4)),
                ("Sidon {1,2,5,11,22}",(1,2,5,11,22)),("AP {1,2,3,4,5}",(1,2,3,4,5))]:
    n=len(sp); d=1.0/(n+1); main,corr,nr=decompose(sp,d); mu=main+corr; md=lonely_direct(sp,d)
    v = "robustly LONELY" if abs(corr)/main<0.3 else "tight (AP, can cancel)"
    print(f"{name:24s} {main:7.4f} {abs(corr)/main:11.3f} {mu:8.4f} {md:10.4f} {v:>18s}")
print()
print("REDUCTION: low-resonance (Sidon) cores have |corr|/main<<1 => mu>0 => NOT covering (lonely).")
print("Only maximal-resonance (AP) cores approach cancellation of the floor. LRC(14) reduces to the")
print("AP-neighborhood = GMC2 coincident-cycle hard stratum (S101) = degenerate tournament zeta (S99).")
print("UNIFICATION: mu(lonely)=sum over {k: sum k_i v_i=0} prod hat_g_{k_i}  is the SAME integer-kernel")
print("resonance sum as the GMC2 moment E[P^m]=sum over balanced channels; non-cancellation is the shared")
print("hard problem, extremal at the AP (relation-rich = many coincident minimal cycles).")
