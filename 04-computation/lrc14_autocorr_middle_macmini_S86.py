#!/usr/bin/env python3
"""mac-mini-S86: BUILD the autocorrelation bound on meas(G' cap middle). THM-729's device is a
1-D autocorrelation discrepancy (density route, sqrt-cancellation). Test whether meas(G' cap
middle) is 1-D-tractable or genuinely MULTI-LINEAR (Gowers, the covering asymmetry THM-729 names).
meas(G' cap middle) = INT_middle PROD_w(1-1_{D_w}). If independent, = meas(middle)*PROD(1-2c).
The GAP (independent - actual) = the correlation reduction; is it PAIRWISE (2nd moment, THM-729-
tractable) or MULTI-LINEAR (needs the full product, Gowers)? Decompose by truncation order."""
import numpy as np
from itertools import combinations
c=1.0/14; P=[2,3,5,7,11,13]
def cop(v): return all(v%p!=0 for p in P)
def build(S,N=300000):
    x=(np.arange(N)+0.5)/N
    nc=[v for v in S if not cop(v)]
    mid=((x>=1/14)&(x<=13/14)).astype(float)
    def dang(v): r=(v*x)%1.0; return (np.minimum(r,1-r)<c).astype(float)
    D=[dang(w) for w in nc]
    return x,nc,mid,D
for nm,S in [("{1..11,13,84}",[*range(1,12),13,84]),("deep well {1..12,182}",[*range(1,13),182])]:
    S=sorted(set(S)); x,nc,mid,D=build(S); N=len(x)
    Gp=np.ones(N)
    for d in D: Gp*= (1-d)
    actual=(Gp*mid).mean()
    indep=(mid.mean())*np.prod([1-d.mean() for d in D])   # if D_w independent
    # truncated inclusion-exclusion in the MIDDLE: order-k = sum over |T|=k of INT_mid prod_{w in T} D_w
    # meas(G' cap mid)=sum_k (-1)^k S_k, S_k = sum_{|T|=k} INT_mid prod D_w (co-danger overlaps in middle)
    Sk=[mid.mean()]
    m=len(nc)
    for k in range(1,min(m,8)+1):
        tot=0.0
        for T in combinations(range(m),k):
            pr=mid.copy()
            for i in T: pr=pr*D[i]
            tot+=pr.mean()
        Sk.append(tot)
    partial=[sum((-1)**k*Sk[k] for k in range(0,J+1)) for J in range(len(Sk))]
    print(f"{nm}: actual meas(G' cap mid)={actual:.5f}, INDEPENDENT est={indep:.5f} (ratio {actual/indep:.3f})")
    print(f"   S_k (co-danger overlaps in middle, |T|=k): "+" ".join(f"{k}:{Sk[k]:.3f}" for k in range(min(6,len(Sk)))))
    print(f"   partial IE sum_(0..J): "+" ".join(f"{p:+.3f}" for p in partial[:8]))
print("\nREAD: if actual << independent (ratio<<1), the smooth runners STRONGLY correlate (near-blanket)")
print("in the middle -- the reduction is the covering object. If the partial IE sums settle only at")
print("HIGH order (not order 2), the reduction is MULTI-LINEAR (Gowers), NOT the 1-D THM-729 discrepancy")
print("=> THM-729's density device does NOT transfer; the covering middle-order needs Gowers (the")
print("asymmetry THM-729 names). Honest test of whether the autocorrelation bound can be built here.")
