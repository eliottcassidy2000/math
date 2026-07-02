#!/usr/bin/env python3
"""
klein-2026-07-02-S112 (HYP-4017) - band margin analysis: the epsilon-schedule base data.

For every covering 13-subset of [1..W] (the censused band), compute the EXACT max-min
optimum  M(S) = max_t min_i ||s_i t||  over the exhaustive breakpoint grid, and the
margin M(S) - 1/14.

WHY: kps-S17's design note — the analytic peel (far-element transport) needs the
induction invariant to carry an INTERVAL of lonely times, not a point. A row with
margin mu > 0 certifies an interval of length >= 2*mu/max_speed around its witness
(each runner's distance moves at speed <= max_speed). A TIGHT row (mu = 0) certifies
only isolated points — no interval. So the epsilon-schedule base = the non-tight rows;
the tight rows must be classified and handled by their own (finite) neighborhood
analysis. This script produces that classification.
"""
from itertools import combinations
from fractions import Fraction
import time

def covering(S):
    return all(any(s % q == 0 for s in S) for q in range(2, 15))

def exact_opt(S):
    """Exact max-min over the exhaustive breakpoint grid; returns (num, den) reduced."""
    qs = set()
    for i, u in enumerate(S):
        qs.add(2 * u)
        for v in S[i+1:]:
            qs.add(u + v)
            qs.add(v - u)
    qs.discard(0)
    bn, bd = 0, 1   # best value num/den
    arg = None
    for q in qs:
        for a in range(1, q):
            k = q  # min over runners of min(m, q-m)
            for s in S:
                m = (s * a) % q
                mm = m if m <= q - m else q - m
                if mm < k:
                    k = mm
                    if 14 * k < q and k * bd <= bn * q:
                        break
            if k * bd > bn * q:
                bn, bd, arg = k, q, (a, q)
    return bn, bd, arg

for W in (20, 22):
    t0 = time.time()
    tight, results = [], []
    n = 0
    for S in combinations(range(1, W + 1), 13):
        if not covering(S):
            continue
        n += 1
        bn, bd, arg = exact_opt(S)
        M = Fraction(bn, bd)
        mu = M - Fraction(1, 14)
        assert mu >= 0, (S, M)
        results.append((mu, S, arg, M))
        if mu == 0:
            tight.append((S, arg))
    results.sort()
    print(f"\n=== W = {W}: {n} covering rows ({time.time()-t0:.0f}s) ===")
    print(f"TIGHT rows (optimum exactly 1/14, zero-length interval): {len(tight)}")
    for S, arg in tight:
        print(f"   {S}   witness {arg[0]}/{arg[1]}")
    pos = [r for r in results if r[0] > 0]
    if pos:
        mu_min = pos[0]
        print(f"smallest POSITIVE margin: {mu_min[0]} = {float(mu_min[0]):.3e}  at {mu_min[1]} (opt {mu_min[3]})")
        print(f"  -> uniform interval floor for non-tight rows: 2*mu/W = {float(2*mu_min[0]/W):.3e}")
        qs_ = sorted(set(float(r[0]) for r in pos))
        import statistics
        print(f"margin stats over {len(pos)} non-tight rows: min {float(pos[0][0]):.3e}, "
          f"median {float(pos[len(pos)//2][0]):.3e}, max {float(pos[-1][0]):.3e}")
    # margin histogram by decade
    from collections import Counter
    import math
    c = Counter()
    for mu, S, arg, M in pos:
        c[math.floor(math.log10(float(mu)))] += 1
    print("margin histogram (log10 decade -> count):", dict(sorted(c.items())))
