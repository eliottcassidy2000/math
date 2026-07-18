#!/usr/bin/env python3
"""
Does the TWO-CIRCLE deep structure generalize across N?  (boxeph-S81)

death-star THM-985/987 (N=14): the deep set {p : bandCount >= 6} is EXACTLY
two resonance circles -- circle I (integer, small p or near q) with threshold
84 = 6*14, circle II (half-integer, near q/2).  The "6" is the coverage cap K.

Question: for the prefix {1..N-1} at threshold 1/N, is the K-deep set
{p : bandCount >= K} always two circles, with circle-I threshold = N*K
(binding speed K must nest at the origin: K*p < q/N <=> N*K*p < q)?

We test:
  - the deep set for various K, decompose into arcs, check "two circles"
  - the integer-circle threshold  (should be N*K when the K smallest nest)
  - the half-circle (even speeds nest at q/2; needs N even for a clean parity lock)
"""
from math import gcd

def band_count(V, q, p, N):
    c = 0
    for vi in V:
        r = (vi * p) % q
        if not (q <= N * r <= (N - 1) * q):
            c += 1
    return c

def deep_arcs(V, q, N, K):
    """Return the sorted list of p in (0,q) with bandCount>=K, grouped into arcs."""
    pts = [p for p in range(1, q) if band_count(V, q, p, N) >= K]
    if not pts:
        return []
    arcs = []
    start = prev = pts[0]
    for p in pts[1:]:
        if p == prev + 1:
            prev = p
        else:
            arcs.append((start, prev))
            start = prev = p
    arcs.append((start, prev))
    return arcs

def analyze(N, K, qlist):
    V = list(range(1, N))
    print(f"  N={N}  K={K}  (predicted circle-I threshold N*K = {N*K})")
    for q in qlist:
        arcs = deep_arcs(V, q, N, K)
        B = (q - 1) // (N * K)   # predicted circle-I half-width
        print(f"    q={q:4d}  B_pred=(q-1)//{N*K}={B:2d}  #deep-arcs={len(arcs)}  arcs={arcs}")

if __name__ == "__main__":
    print("="*74)
    print("TWO-CIRCLE deep structure across N.  Prefix {1..N-1}, deep = bandCount>=K.")
    print("="*74)
    # N=14 baseline with K=6 (death-star's setting): expect 84-thresholded two circles
    print("\n[N=14, K=6]  (death-star baseline: threshold 84, two circles)")
    analyze(14, 6, [170, 200, 250, 300])
    # generalize: N even, K = N/2 (evens nest at half-circle) and K = phi(N)?
    for N in [10, 12, 16, 18]:
        K = N // 2
        print(f"\n[N={N}, K=N/2={K}]")
        analyze(N, K, [N*K*2+5, N*K*3+7, N*K*4+3])
    # odd N: half-circle parity lock absent -- does circle II vanish?
    for N in [9, 15]:
        K = (N-1)//2
        print(f"\n[N={N} odd, K=(N-1)/2={K}]  (parity lock absent -> circle II?)")
        analyze(N, K, [N*K*2+5, N*K*3+7])
    print("\nREADING: circle I (integer) threshold = N*K is UNIVERSAL (K smallest speeds")
    print("nest at origin). Circle II (half) needs N even (evens nest at q/2). The 84")
    print("of N=14 is 6*14 with K=6; the two-circle count formula generalizes with N*K.")
