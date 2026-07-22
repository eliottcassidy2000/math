"""S102 finite-truncation illustration, corrected by MISTAKE-235.

The exact LRC relation-lattice expansion is valid with proper convergence
control. This finite K=9 experiment neither identifies it with a GMC moment nor
reduces LRC(n+1) to AP cores.
"""
import numpy as np
from math import sin, pi
from itertools import product

print("CORRECTED BY MISTAKE-235: K=9 samples illustrate an LRC kernel expansion; they do not identify GMC or reduce LRC(14) to AP cores.")

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

print("FINITE K=9 BOX SUM versus a midpoint-grid estimate; no Fourier-tail bound is asserted.")
print(f"{'core':24s} {'main':>7s} {'|box|/main':>11s} {'main+box':>8s} {'grid est.':>10s} {'scope':>18s}")
for name,sp in [("Sidon {1,2,5,11}",(1,2,5,11)),("AP {1,2,3,4}",(1,2,3,4)),
                ("Sidon {1,2,5,11,22}",(1,2,5,11,22)),("AP {1,2,3,4,5}",(1,2,3,4,5))]:
    n=len(sp); d=1.0/(n+1); main,corr,nr=decompose(sp,d); mu=main+corr; md=lonely_direct(sp,d)
    v = "finite sample only"
    print(f"{name:24s} {main:7.4f} {abs(corr)/main:11.3f} {mu:8.4f} {md:10.4f} {v:>18s}")
print()
print("NO REDUCTION: four n=4/5 rows and a K=9 box do not control the full Fourier tail or any LRC(14) row.")
print("SURVIVOR: the LRC measure has a Fejer-regularized weighted relation-lattice expansion.")
print("TYPE WARNING: a fixed GMC moment is a different finite mass-graded fiber with different weights;")
print("no map to GMC or tournament zeta and no AP-neighborhood theorem is supplied. See MISTAKE-235.")
