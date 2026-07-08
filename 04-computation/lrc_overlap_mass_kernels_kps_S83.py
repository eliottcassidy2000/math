#!/usr/bin/env python3
r"""
lrc_overlap_mass_kernels_kps_S83.py  (kind-pasteur-2026-07-08-S83, HYP-5397)

THE TRIPLE / QUAD (and general) OVERLAP MASS FOURIER KERNELS -- derived + verified.

For a subset S of the k arcs (A_i = 1/7-arc at frac(e_i x)), the overlap L_S(x) =
meas(cap_{i in S} A_i).  Its mean over x (the "S-overlap mass") is, in Fourier form:

    E[L_S] = sum_{m in Lambda_S} prod_{a in S} c_{m_a},    c_m = int_0^{1/7} e(-2pi i m t) dt,

    Lambda_S = { m in Z^S : sum_a m_a = 0  AND  sum_a m_a e_a = 0 }   (verified below).

Lambda_S is the BALANCED ADDITIVE-RELATION LATTICE of S, of rank |S| - 2:
  * |S|=2: rank 0, Lambda={0}, E[L_ij] = c_0^2 = 1/49 (the pair-overlap mean).
  * |S|=3: rank 1 -- Lambda = Z * (d1,d2,d3)/gcd, the PRIMITIVE TRIANGLE of differences
           (d1=e_j-e_k, d2=e_k-e_i, d3=e_i-e_j; d1+d2+d3=0).  gcd = the dilation.
           E[L_ijk] = (1/7)^3 + sum_{t!=0} prod c_{t*(d/gcd)}.
  * |S|=4: rank 2 -- Lambda = triangle vectors + Sidon-violation vectors (additive relations
           of the quad).  E[L_ijkl] = sum over the rank-2 lattice.
  * general |S|: rank |S|-2 lattice of additive relations.
APEX-7: |c_m|^2 = sin^2(pi m/7)/(pi^2 m^2) = 0 at m = 0 mod 7, so relations with a component
divisible by 7 lose those harmonics (the THM-637/638 apex-7 invisibility, generalized).
DILATION-INVARIANT: E[L_S] depends only on the primitive relation structure (2*AP, 3*AP, AP
all give the same mass), consistent with mu_{1/7} dilation-invariance.

This is the |S|>=3 extension of THM-641 (the |S|=2 pair-mass law).  It is the building block
of Var(W) = sum_{|S|,|T|>=2} (-1)^{|S|+|T|} Cov(L_S,L_T); the full resonance needs the
2-window joint masses (a further lattice extension) and is NON-perturbative (kps-S82).
"""
import numpy as np, itertools
from math import gcd
from functools import reduce

def c(m):
    return 1/7 if m == 0 else (1 - np.exp(-2j*np.pi*m/7)) / (2j*np.pi*m)
def overlap(phs):
    yg = (np.arange(6000)+0.5)/6000; cov = np.ones(6000, bool)
    for p in phs: cov &= (np.mod(yg - p, 1.0) < 1/7)
    return cov.mean()
def EL_direct(elts, res=50000):
    xs = (np.arange(res)+0.5)/res
    return float(np.mean([overlap([(e*x) % 1.0 for e in elts]) for x in xs]))
def EL_lattice(elts, M=26):
    n = len(elts); s = 0j; rng = range(-M, M+1)
    for m in itertools.product(rng, repeat=n-1):
        mn = -sum(m)
        if abs(mn) > M: continue
        mm = m + (mn,)
        if sum(mm[a]*elts[a] for a in range(n)) != 0: continue
        p = 1.0
        for mi in mm: p *= c(mi)
        s += p
    return s.real
def EL_triple(elts, T=400):  # rank-1 closed form
    i, j, k = elts; d = (j-k, k-i, i-j); g = reduce(gcd, [abs(x) for x in d]) or 1
    dp = tuple(x//g for x in d)
    return (1/7)**3 + sum((c(t*dp[0])*c(t*dp[1])*c(t*dp[2])).real for t in range(-T, T+1) if t)

print("TRIPLE overlap mass kernel (rank-1 primitive triangle):")
for elts in [[0,1,3],[0,1,7],[0,3,12],[0,7,20]]:
    print(f"  {str(elts):>12}: direct={EL_direct(elts):.6f}  triangle-kernel={EL_triple(elts):.6f}")
print("QUAD overlap mass kernel (rank-2 additive-relation lattice):")
for elts in [[0,1,2,3],[0,2,4,6],[0,3,6,9],[0,2,3,5]]:
    print(f"  {str(elts):>14}: direct={EL_direct(elts):.6f}  lattice-kernel={EL_lattice(elts):.6f}")
print("=> E[L_S] = sum_{m in Lambda_S} prod c_{m_a}, Lambda_S = balanced additive-relation")
print("   lattice (rank |S|-2). |S|=2 pair = THM-641; |S|=3,4 = this. apex-7 + dilation-invariant.")
