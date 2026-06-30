#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THE LEAST-EIGENVALUE CERTIFICATE for the apex sigma-EVEN lonely measure (klein-S19).

mac-mini-S6 (THM-581/582): the LRC floor is entirely sigma-EVEN (lonely = even; Redei = odd does NOT
apply -- a lonely tournament is not self-converse). So the floor's certificate is a Bochner-positive
(SOS = sum of squares = even) object: the cyclotomic Gram. THM-590's g(O) = min_{k!=0}|sum_{x in O}
zeta^{kx}|^2 is EXACTLY the smallest nonzero eigenvalue of the core's autocorrelation circulant C(O).
This script makes that explicit (matrix + eigenvalues + SOS eigenvector + the sigma=negation structure),
verifies lambda_min = 4cos^2(3pi/7) at the doublet, and tests the lift to the full lonely measure.
"""
import math, cmath, itertools
import numpy as np

P = 7
W = cmath.exp(2j*math.pi/P)

def autocorr_circulant(O):
    """C_{ij} = a((i-j) mod 7), a(d) = #{x in O: (x+d) mod 7 in O}. Real symmetric circulant."""
    Oset = set(x % P for x in O)
    a = [sum(1 for x in Oset if (x+d) % P in Oset) for d in range(P)]
    C = np.array([[a[(i-j) % P] for j in range(P)] for i in range(P)], dtype=float)
    return C, a

def g_formula(O):
    Oset = set(x % P for x in O)
    return min(abs(sum(W**((k*x) % P) for x in Oset))**2 for k in range(1,P))

print("="*84)
print(" THE LEAST-EIGENVALUE CERTIFICATE: g(O) = smallest nonzero eigenvalue of the cyclotomic Gram")
print("="*84)
doublet = 4*math.cos(3*math.pi/7)**2
print(f" target doublet value 4cos^2(3pi/7) = {doublet:.6f}\n")

# (1) verify lambda_min(nonzero) of C(O) == g(O) for ALL cores, matching THM-590's 5 values
vals=set()
mism=0
for r in range(1,P+1):
    for O in itertools.combinations(range(P), r):
        C,_ = autocorr_circulant(O)
        eig = sorted(np.linalg.eigvalsh(C))           # real symmetric => real eigenvalues
        # Perron = |O|^2 at k=0; nonzero spectral gap = smallest eigenvalue that is not the k=0 mode
        # eigenvalues are |Ohat(k)|^2; the k=0 one equals |O|^2 (the largest). gap = min over k!=0.
        # numerically: drop ONE copy of the Perron |O|^2, take min of the rest
        eig_nz = eig[:]  # all; the min eigenvalue IS min_k |Ohat(k)|^2 over k!=0 unless |O| in {0,7}
        lam_min = min(eig)            # smallest eigenvalue overall
        g = g_formula(O)
        if abs(lam_min - g) > 1e-6 and not (len(O) in (0,7)):
            mism += 1
        vals.add(round(g,6))
print(f"[1] For every proper nonempty core O: smallest eigenvalue of C(O) == g(O) (THM-590).  mismatches={mism}")
print(f"    distinct gap values = {sorted(vals)}")
print(f"    (these are THM-590's five values; min nonzero = {min(v for v in vals if v>1e-9):.6f} = 4cos^2(3pi/7))")

# (2) the DOUBLET as an explicit least-eigenvalue certificate: matrix, eigenvalues, SOS eigenvector
print("\n[2] The binding DOUBLET O={0,1}: the explicit certificate")
C,a = autocorr_circulant({0,1})
print(f"    autocorrelation a(d), d=0..6 = {a}")
w_,V = np.linalg.eigh(C)
print(f"    Gram C eigenvalues (sorted) = {[round(x,4) for x in sorted(w_)]}")
print(f"    smallest eigenvalue lambda_min = {min(w_):.6f}  vs 4cos^2(3pi/7) = {doublet:.6f}")
imin = int(np.argmin(w_))
vec = V[:,imin]
print(f"    SOS/least eigenvector v (the dual certificate) = {[round(x,3) for x in vec]}")
print(f"    check C v = lambda v: ||Cv - lam v|| = {np.linalg.norm(C@vec - w_[imin]*vec):.2e}")
print(f"    Bochner: ALL eigenvalues >= 0 ({all(x>=-1e-9 for x in w_)}) -- C is PSD = an SOS (sigma-even).")

# (3) the sigma = NEGATION (k <-> -k) structure: eigenvalues pair, the gap lives in Q(cos 2pi/7)
print("\n[3] sigma = negation (k <-> 7-k): C is real-symmetric circulant, so lambda_k = lambda_{7-k}")
print("    (the de Moivre cos-pairing). The 6 nonzero modes pair into 3 sigma-orbits -> a CUBIC")
print("    over Q(cos 2pi/7) {cos2pi/7, cos4pi/7, cos6pi/7}. The gap is a sigma-EVEN (cosine) value:")
for O in [{0,1},{0,1,2},{1,2,4}]:
    Oset=set(O); 
    lam=[round(abs(sum(W**((k*x)%P) for x in Oset))**2,4) for k in range(P)]
    print(f"    O={sorted(O)}: |Ohat(k)|^2 k=0..6 = {lam}  (k<->7-k symmetric: {lam[1]==lam[6] and lam[2]==lam[5] and lam[3]==lam[4]})")

# (4) the LIFT test: does the apex gap lower-bound the full lonely-measure spectral object?
print("\n[4] LIFT to the full lonely measure (the non-factoring test):")
import sys; sys.path.insert(0,'04-computation')
M=__import__('lrc14_floor_CV_sheetcount_bound_macmini_20260629'); lonely_set, measure = M.lonely_set, M.measure
def full_min_spectral(speeds, N=14*30):
    # discretize lonely set on Z_N, autocorrelation power spectrum |Lhat(k)|^2 over k!=0, smallest
    L = lonely_set(speeds)
    ind = np.zeros(N)
    for (lo,hi) in L:
        a0=int(math.ceil(float(lo)*N)); b0=int(math.floor(float(hi)*N))
        for t in range(a0,b0+1): ind[t%N]=1.0
    F = np.fft.fft(ind)
    pk = (np.abs(F)**2)/N
    nz = pk[1:]
    return float(min(nz)), float(pk[0]/N)  # min nonzero mode (could be ~0), and the DC
for name,O in [("doublet-ish {1,8}",[1,8]),("{1,3,5}",[1,3,5]),("{1,2,4}",[1,2,4])]:
    fm,dc = full_min_spectral(O)
    print(f"    {name}: full min |Lhat(k)|^2 = {fm:.4f} (DC=m^2 region {dc:.3f}); apex g(O mod7)={g_formula(O):.4f}")
print("    NOTE: the full power spectrum HAS near-zeros (L̂ vanishes at many k) -- the full 'least")
print("    eigenvalue' is ~0 generically. The certificate is NOT the naive full-spectrum min; it is")
print("    the APEX (mod-7) cyclotomic gap (THM-590), reached after the 2-adic descent factors out the")
print("    2-part. The 'danger does not factor' => the certificate lives at the apex, not the full grid.")
