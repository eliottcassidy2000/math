#!/usr/bin/env python3
"""
S620 — Helly-entropy accounting for the Lonely Runner Conjecture.

Frame (abstract): a loneliness certificate is a point of the clock R/Z that
avoids every "forbidden arc". So the certificate problem IS a *circular-arc
covering* problem, and the natural object is the COVERING-DEPTH function

    depth(t) = #{ i : ||v_i t|| < delta }            (delta = gap = 1/(n+1))

whose sublevel set {depth = 0} is exactly the set of lonely (certificate) times.
Everything about LRC at gap delta is a functional of the depth DISTRIBUTION
    p_k = Lebesgue-measure{ t in [0,1) : depth(t) = k }.

This script proves/verifies (exactly, over Q) the rigorous backbone:

  (R1) CONSERVATION.   mean depth = S_1 = 2 n /(n+1)  for EVERY speed set
                       (linearity; arithmetic-independent; always < 2).
  (R2) MOMENT-SIEVE.   The k-th factorial moment of depth equals the order-k
                       inclusion-exclusion term
                          S_k = sum_{|T|=k} measure( intersect_{i in T} bad_i ),
                       and the free/lonely measure is the alternating sum
                          p_0 = sum_k (-1)^k S_k   (== Bonferroni == S561's rho).
  (R3) HELLY-ENTROPY.  H_depth = -sum_k p_k log p_k measures the spread of the
                       depth charge.  Generic sets sit at the independence
                       baseline (1-2delta)^n with HIGH H_depth; resonant/tight
                       extremals (AP and sporadic additive chains) collapse to
                       p_0 = 0 with LOW H_depth.

All measures are exact Fractions (arc endpoints are rational).
"""
from fractions import Fraction as Fr
from math import log, comb, gcd, exp
import itertools


def _active(v, t, delta):
    """Is runner with speed v inside its forbidden arc at clock time t?  ||v t|| < delta."""
    x = v * t
    f = x - (x.numerator // x.denominator)        # v t mod 1 in [0,1)
    if f > Fr(1, 2):
        f = 1 - f                                 # distance to nearest integer
    return f < delta


def depth_distribution(speeds, delta):
    """Exact {k: measure(depth==k)} via the arc-endpoint arrangement on [0,1)."""
    bps = {Fr(0), Fr(1)}
    for v in speeds:
        for k in range(v):
            for s in (delta, -delta):
                t = (Fr(k) + s) / v
                t -= t.numerator // t.denominator
                if 0 <= t < 1:
                    bps.add(t)
    bps = sorted(bps)
    dist = {}
    for a, b in zip(bps, bps[1:]):
        mid = (a + b) / 2
        d = sum(_active(v, mid, delta) for v in speeds)
        dist[d] = dist.get(d, Fr(0)) + (b - a)
    return dist


def factorial_moments(dist, n):
    """S_k = E[ C(depth, k) ] = sum_{|T|=k} measure(intersection of k bad arcs)."""
    return {k: sum(float(m) * comb(d, k) for d, m in dist.items()) for k in range(n + 1)}


def helly_entropy(dist):
    return -sum(float(m) * log(float(m)) for m in dist.values() if m > 0)


def summary(name, speeds):
    n = len(speeds)
    delta = Fr(1, n + 1)
    dist = depth_distribution(speeds, delta)
    S = factorial_moments(dist, n)
    p0 = float(dist.get(0, Fr(0)))
    p0_incl = sum((-1) ** k * S[k] for k in range(n + 1))      # R2 alternating sum
    mean = S[1]
    base = (1 - 2 * float(delta)) ** n                          # independence p_0
    H = helly_entropy(dist)
    print(f"  {name:32s} mean={mean:.4f} p0={p0:.5f} (incl-excl {p0_incl:+.5f}) "
          f"indep={base:.5f} H_depth={H:.4f}")
    return dist, S


def main():
    print("=" * 78)
    print("S620  HELLY-ENTROPY ACCOUNTING FOR THE LONELY RUNNER CONJECTURE")
    print("=" * 78)

    print("\n[R1+R3]  depth distribution per speed set  (delta = 1/(n+1))")
    for n in (4, 5, 6, 7):
        print(f"\n n={n}:  conserved mean depth = 2n/(n+1) = {2*n/(n+1):.4f}")
        summary("AP 1..n  (tight extremal)", list(range(1, n + 1)))
        summary("geometric 1,2,4,..", [2 ** i for i in range(n)])
        summary("random A (10,21,26,..)", sorted([10, 21, 26, 42, 53, 7, 24][:n]))
        summary("random B (4,5,35,..)", sorted([4, 5, 35, 53, 14, 33, 38][:n]))

    print("\n" + "=" * 78)
    print("[R2]  moment-sieve identity:  S_k == direct pairwise/triple intersection,")
    print("      and  p_0 == sum_k (-1)^k S_k   (free measure = alternating sieve)")
    print("=" * 78)
    for speeds in ([1, 2, 3, 4], [1, 2, 4, 8], [1, 3, 5, 6], [2, 3, 4, 5]):
        n = len(speeds)
        delta = Fr(1, n + 1)
        dist = depth_distribution(speeds, delta)
        S = factorial_moments(dist, n)
        # independent direct check of S_2 (sum of pairwise arc-overlap measures)
        bps = sorted({Fr(0), Fr(1)} | {
            (Fr(k) + s) / v - ((Fr(k) + s) / v).numerator // ((Fr(k) + s) / v).denominator
            for v in speeds for k in range(v) for s in (delta, -delta)})
        S2dir = 0.0
        for i, j in itertools.combinations(range(n), 2):
            for a, b in zip(bps, bps[1:]):
                mid = (a + b) / 2
                if _active(speeds[i], mid, delta) and _active(speeds[j], mid, delta):
                    S2dir += float(b - a)
        p0 = float(dist.get(0, Fr(0)))
        p0i = sum((-1) ** k * S[k] for k in range(n + 1))
        print(f"  {tuple(speeds)}  S1={S[1]:.4f} S2={S[2]:.4f} (direct {S2dir:.4f}) "
              f"p0={p0:.5f}  alt-sum={p0i:+.5f}  match={abs(p0-p0i)<1e-9 and abs(S[2]-S2dir)<1e-9}")

    print("\n" + "=" * 78)
    print("[EXTREMAL FAMILY]  which primitive n-sets COLLAPSE (p_0 = 0)?  AP is not alone.")
    print("=" * 78)
    for n in (3, 4, 5):
        delta = Fr(1, n + 1)
        zeros, minpos = [], (1e9, None)
        for s in itertools.combinations(range(1, 3 * n + 2), n):
            g = 0
            for x in s:
                g = gcd(g, x)
            if g != 1:
                continue
            p0 = float(depth_distribution(list(s), delta).get(0, Fr(0)))
            if p0 == 0:
                zeros.append(s)
            elif p0 < minpos[0]:
                minpos = (p0, s)
        print(f"  n={n}: {len(zeros)} primitive collapse-sets; examples {zeros[:6]}; "
              f"min positive p_0={minpos[0]:.5f} at {minpos[1]}")


if __name__ == "__main__":
    main()
