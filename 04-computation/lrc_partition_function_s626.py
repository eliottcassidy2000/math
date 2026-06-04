#!/usr/bin/env python3
"""
S626 — Partition functions everywhere: the covering-depth as a statistical-mechanics
partition function, its transfer recursion, and the n→n+2 (±-pair / even-fold) structure.
Inspired by infinite-Go transfinite game values (Hamkins/Evans): game value = ordinal built
by recursion (= the iterated-log 'altitude', HYP-2180); ordinal natural sum = product of
generating functions; the 'n+2' skip-recursion = stepping by the ±-pair 2-orbit (S625 conj.).

Convention (canon THM-398): n runners, speeds {1,..,n-1} for the AP, gap 1/n.
COVERING-DEPTH PARTITION FUNCTION:
    Z_S(z) = ∫_0^1 z^{depth_S(t)} dt = Σ_k p_k z^k,   depth_S(t)=#{v∈S: ||v t||<1/n}.
  z = fugacity, runners = sites, depth = occupation number.  Then (THM-410):
    Z_S(1)=1,  Z_S(0)=p_0 (free measure / 'ground state'),  Z_S'(1)=mean depth=2(n-1)/n,
    Z_S(z) = Σ_k S_k (z-1)^k   (S_k = inclusion-exclusion / sieve moments).
TRANSFER RECURSION (adding a runner v):
    Z_{S∪{v}}(z) = Z_S(z) + (z-1)·∫_{bad_v} z^{depth_S(t)} dt.
"""
from fractions import Fraction as Fr
from math import gcd

def norm(x):
    f = x - (x.numerator // x.denominator)
    return f if f <= Fr(1, 2) else 1-f

def depth_dist(speeds, delta):
    bps = {Fr(0), Fr(1)}
    for v in speeds:
        for k in range(v):
            for s in (delta, -delta):
                t = (Fr(k)+s)/v; t -= t.numerator//t.denominator
                if 0 <= t < 1: bps.add(t)
    bps = sorted(bps); dist = {}
    for a, b in zip(bps, bps[1:]):
        mid = (a+b)/2
        d = sum(norm(v*mid) < delta for v in speeds)
        dist[d] = dist.get(d, Fr(0)) + (b-a)
    return dist

def Z_poly(dist):
    """partition function coefficients p_k (as Fractions), index = depth."""
    K = max(dist)
    return [dist.get(k, Fr(0)) for k in range(K+1)]

def eval_Z(p, z):
    return sum(c*z**k for k, c in enumerate(p))

def sieve_moments(p):
    """S_k = Σ_j C(j,k) p_j  = the inclusion-exclusion terms (THM-410)."""
    from math import comb
    K = len(p)
    return [sum(comb(j, k)*p[j] for j in range(K)) for k in range(K)]

if __name__ == "__main__":
    print("COVERING-DEPTH PARTITION FUNCTION Z_S(z)=Σ p_k z^k for the AP {1,..,n-1}, gap 1/n")
    print("="*74)
    Zs = {}
    for n in range(3, 10):
        delta = Fr(1, n); S = list(range(1, n))
        dist = depth_dist(S, delta); p = Z_poly(dist); Zs[n] = p
        Sk = sieve_moments(p)
        mean = sum(k*float(c) for k, c in enumerate(p))
        p0 = float(p[0]) if p else 0.0
        # verify Z(z)=Σ S_k (z-1)^k at a test z
        z = Fr(7, 3)
        lhs = eval_Z(p, z); rhs = sum(Sk[k]*(z-1)**k for k in range(len(Sk)))
        print(f" n={n}: p_k={[str(c) for c in p]}")
        print(f"      Z(0)=p0={p0:.4f}  mean=Z'(1)={mean:.4f} (=2(n-1)/n={2*(n-1)/n:.4f})"
              f"  S_k-identity holds={lhs==rhs}  S1={Sk[1] if len(Sk)>1 else '-'}")

    print("\n n→n+2 (±-pair / even-fold) structure: support length K(n)=max depth, and parity split")
    for n in range(3, 10):
        K = len(Zs[n])-1
        print(f"   n={n}: deg Z = {K}  (#depth levels {K+1})   p0={float(Zs[n][0]):.4f}")
    print("\n   even-n vs odd-n p0 sequences (the 2-adic/±-pair recursion lives here):")
    print("   even n:", [(n, round(float(Zs[n][0]), 4)) for n in range(4, 10, 2)])
    print("   odd  n:", [(n, round(float(Zs[n][0]), 4)) for n in range(3, 10, 2)])
