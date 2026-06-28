"""
S73c: ANGLE A -- a spectral proof route for the EVEN/provable half of the k=8 node.
kps S31ai: LRC(14) bounded core reduces to "consec maximizes total empty-sector covariance
  Sigma_{i<j} Cov(X_i,X_j)"  where X_j = 1[inner sector j empty], over the 6 inner sectors.
CREATIVE ANGLE: the AP orbit {frac(i x)} is a COHERENT rotation -> its 6x6 covariance matrix C
should be CIRCULANT, and the total covariance 1^T C 1 should be governed by C's PERRON eigenvalue
with the all-ones vector 1 as (near) Perron eigenvector. If so, "AP maximizes 1^T C 1" becomes a
spectral (Perron-Frobenius / circulant) statement -- the cleanest possible form of the even half.
TEST: is C circulant for consec? is 1 its top eigenvector? does consec maximize 1^T C 1 vs others?
"""
from fractions import Fraction as F
import numpy as np

def sector_of(p): return int((p%1)*7)
def cov_matrix(E):
    """6x6 covariance of empty-indicators X_1..X_6 (inner sectors 1..6; sector 0 = observer, always hit)."""
    E=sorted(set(E)); b=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): b.add(F(m,7*e))
    b=sorted(b)
    p=np.zeros(7); P=np.zeros((7,7)); tot=F(0)
    for i in range(len(b)-1):
        x0,x1=b[i],b[i+1]
        if x1<=x0: continue
        w=x1-x0; tot+=w
        hit=set(sector_of(e*((x0+x1)/2)) for e in E)
        empty=[s for s in range(7) if s not in hit]   # which sectors empty on this interval
        fw=float(w)
        for s in empty: p[s]+=fw
        for s in empty:
            for t in empty: P[s][t]+=fw
    # restrict to inner sectors 1..6
    idx=list(range(1,7))
    pe=p[idx]; Pe=P[np.ix_(idx,idx)]
    C=Pe-np.outer(pe,pe)
    return C, pe

def is_circulant(C, tol=1e-6):
    n=C.shape[0]
    # circulant: C[i,j] depends only on (i-j) mod n
    for d in range(n):
        vals=[C[i,(i+d)%n] for i in range(n)]
        if max(vals)-min(vals)>tol: return False, d, max(vals)-min(vals)
    return True, None, 0.0

def total_cov(C): return float(C.sum()-np.trace(C))/1.0   # sum of off-diagonal = 2*Sigma_{i<j}

sets={
 "consec {0..6}": tuple(range(7)),
 "consec {0..7}": tuple(range(8)),
 "consec {0..12}": tuple(range(13)),
 "2*consec {0,2..12}": tuple(range(0,13,2)),
 "skip {0,1,2,3,4,5,7}": (0,1,2,3,4,5,7),
 "dissoc {0,1,2,4,8,16,32}": (0,1,2,4,8,16,32),
 "primes {0,2,3,5,7,11,13}": (0,2,3,5,7,11,13),
 "random {0,1,5,11,17,23,30}": (0,1,5,11,17,23,30),
}
print("="*92)
print(" ANGLE A: empty-sector covariance C (6x6, inner sectors) -- circulant? Perron? consec = max 1^T C 1 ?")
print("="*92)
print(f"{'set':<30}{'Sigma_off C':>13}{'circulant?':>12}{'1 align Perron (cos)':>22}{'lambda_max':>12}")
base=None
for name,E in sets.items():
    C,pe=cov_matrix(E)
    sig=total_cov(C)
    circ,bd,dev=is_circulant(C)
    evals,evecs=np.linalg.eigh(C)
    lmax=evals[-1]; vmax=evecs[:,-1]
    ones=np.ones(6)/np.sqrt(6)
    cosang=abs(float(vmax@ones))                         # |cos| of 1 vs Perron eigvec
    tag="YES" if circ else f"no(d{bd},{dev:.0e})"
    print(f"{name:<30}{sig:>13.5f}{tag:>12}{cosang:>22.4f}{lmax:>12.5f}")
print("-"*92)
# detailed look at consec {0..12} (the LRC(14) size): the C row structure
C,pe=cov_matrix(tuple(range(13)))
print(" consec {0..12} covariance matrix C (6x6, inner sectors 1..6):")
for r in C: print("   ["+" ".join(f"{x:+.4f}" for x in r)+"]")
evals,evecs=np.linalg.eigh(C)
print(f" eigenvalues: {[f'{x:+.4f}' for x in evals]}")
print(f" all-ones vector: 1^T C 1 = {float(np.ones(6)@C@np.ones(6)):+.5f}; lambda_max*6 = {evals[-1]*6:+.5f}")
print(" => if C circulant & 1=Perron, 1^T C1 = 6*lambda_max and consec maximizes it = SPECTRAL proof of even half.")
