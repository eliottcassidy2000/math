#!/usr/bin/env python3
"""mac-mini-S81: pursue the tournament chi=2 orbit-rigidity via the FUGACITY POLYNOMIAL.
Q_S(w)=E[w^X]=Sum_x p_x w^x, p_x=meas{danger count X=x} (nonneg coeffs, prob gen fn).
L(S)=Q_S(0)=p_0; tournament H ~ Q_S(2)=E[2^X]. AP: p_0=0 => w=0 is a ROOT of Q_AP (the chi=2
equioscillation signature). Covering: p_0>0 => w=0 NOT a root. TEST: (1) the p_x distribution
(AP high-concentrated, covering has p_0>0); (2) roots of Q_S(w) -- does the AP have a special
root pattern (w=0 + a group/tournament structure) that covering BREAKS non-locally?"""
import numpy as np
c=1.0/14
def pdist(S,res=600000):
    S=sorted(set(S)); k=len(S); p=np.zeros(k+1)
    ts=(np.arange(res)+0.5)/res
    X=np.zeros(res,dtype=int)
    for v in S:
        r=(v*ts)%1.0; d=np.minimum(r,1-r); X+=(d<c).astype(int)
    for x in range(k+1): p[x]=(X==x).mean()
    return p
def Q_roots(p):
    # Q(w)=sum p_x w^x; roots (numpy: highest degree first)
    coeffs=p[::-1]
    # strip leading zeros (top p_x may be 0)
    nz=np.nonzero(coeffs)[0]
    coeffs=coeffs[nz[0]:]
    return np.roots(coeffs)
print("danger-count distribution p_x and Q_S(w)=E[w^X] structure:\n")
for nm,S in [("AP {1..13}",list(range(1,14))),("{1..11,13,84}",[*range(1,12),13,84]),
             ("deep well {1..12,182}",[*range(1,13),182]),("{2..14}",list(range(2,15)))]:
    S=sorted(set(S)); p=pdist(S)
    p0=p[0]; H=sum(p[x]*(2**x) for x in range(len(p)))  # Q(2)=E[2^X]
    print(f"{nm}: p0=L={p0:.5f}, E[2^X]=Q(2)={H:.2f}")
    print("   p_x = " + " ".join(f"{x}:{p[x]:.3f}" for x in range(len(p)) if p[x]>0.001))
    roots=Q_roots(p)
    real_roots=sorted(r.real for r in roots if abs(r.imag)<1e-6)
    near0=[r for r in real_roots if abs(r)<0.15]
    print(f"   Q(w) real roots: {[round(r,3) for r in real_roots]}")
    print(f"   root near w=0? {[round(r,4) for r in near0] if near0 else 'NONE (p0>0, no root at 0)'}")
print("\nTOURNAMENT/ORBIT READ: AP p0=0 => Q has w=0 root (chi=2 equioscillation). Covering p0>0")
print("=> NO w=0 root. If the AP's root pattern were group-forced (R_13 orbit) and covering broke")
print("it NON-locally, that'd bypass the resummation. HONEST CHECK: is the covering p0>0 a")
print("group-rigidity consequence, or just = M>1/14 (the metric, order-forgets-metric)?")
