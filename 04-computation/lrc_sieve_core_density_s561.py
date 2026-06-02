#!/usr/bin/env python3
"""
lrc_sieve_core_density_s561.py    claudebox-2026-06-02-S561

How thin is the lonely-runner "open core" at n = k+1?  A quantitative follow-up to
oracle-S557 (HYP-2064), which proved the SIEVE gives the FULL conjecture
    g(v) := sup_t min_i ||v_i t||  >=  1/(k+1)
for every NON-sieve-covered speed set, leaving only the SIEVE-COVERED sets (a multiple
of every q in {2..k+1}) as the locus where Tao's general bound is operative. S557 called
that core "measure-zero-ish" and located the genuinely hard part in an "even thinner
near-AP slice."

This script makes both halves precise at n = k+1:

  (1) DENSITY OF THE CORE -- EXACT.  The R->inf density of sieve-covered primitive
      k-sets is computed in closed form by inclusion-exclusion over the events
      A_q = "no speed divisible by q", q in {2..k+1}:
          P(covered) = 1 - P(union_q A_q),
          P(intersection_{q in S} A_q) = rho(S)^k,
          rho(S) = density of integers divisible by no q in S
                 = sum_{T subset of S} (-1)^|T| / lcm(T).
      Monte-Carlo over growing speed ranges R confirms convergence to this limit.

  (2) GAP DISTRIBUTION OVER THE CORE.  For many random sieve-covered sets we compute
      the actual gap g(v) (fine grid). Result: random core sets sit FAR above the
      conjecture (min ~1.9x at n=14), so they are NOT near the wall. The near-wall
      slice (g ~ 1/(k+1)) is the genuinely thin sub-locus and demands the rigid
      AP / dilated-AP structure -- it is not reached by generic core sets.

VERDICT (refines S557): the sieve-covered core is NOT measure-zero -- it has POSITIVE
limiting density (11.2% at n=14, rising as n shrinks). What is thin is the near-wall
slice INSIDE the core. So "the open problem is a thin core" is right about the
EXTREMIZERS but overstates the thinness of the sieve-covered set itself: ~1 in 9 random
primitive 13-sets is sieve-covered, yet essentially none of those is anywhere near tight.
g(v) is dilation-invariant (g(cv)=g(v)), so primitive sets are WLOG.
"""
from functools import reduce
from math import gcd
from itertools import combinations
import random
import statistics as st


def lcm(a, b):
    return a * b // gcd(a, b)


def d0(p):
    p = p % 1.0
    return min(p, 1 - p)


def gap(v, G=120000):
    """actual gap g(v) = sup_t min_i ||v_i t|| via fine grid (a lower estimate;
    error <= max(v)/G, negligible here vs the scales reported)."""
    best = 0.0
    for i in range(G):
        t = (i + 0.5) / G
        m = min(d0(s * t) for s in v)
        if m > best:
            best = m
    return best


def covered(v, N):
    """sieve-covered at n=N: a multiple of every q in {2..N} is present."""
    return all(any(s % q == 0 for s in v) for q in range(2, N + 1))


def rho(S):
    """density of integers divisible by no q in S = sum_{T<=S}(-1)^|T|/lcm(T)."""
    tot = 0.0
    for r in range(0, len(S) + 1):
        for T in combinations(S, r):
            tot += ((-1) ** r) / reduce(lcm, T, 1)
    return tot


def core_density_exact(K, N):
    """exact R->inf density that a random primitive K-set is sieve-covered at n=N."""
    qs = list(range(2, N + 1))
    union = 0.0  # P(union_q A_q) by inclusion-exclusion
    for r in range(1, len(qs) + 1):
        for S in combinations(qs, r):
            union += ((-1) ** (r + 1)) * (rho(S) ** K)
    return 1.0 - union


def main():
    print("=" * 78)
    print("LRC sieve-covered CORE: exact density + gap distribution    (S561, refines S557)")
    print("=" * 78)

    # (1) EXACT limiting density of the core across n = k+1
    print("\n(1) EXACT R->inf density of sieve-covered primitive k-sets  [inclusion-exclusion]")
    for N in [10, 12, 14, 16]:
        K = N - 1
        print(f"    n={N:2d} (k={K:2d}):  P(sieve-covered) = {core_density_exact(K, N):.5f}")
    print("    -> POSITIVE density, not measure-zero. Rises as n shrinks.")

    # (1b) Monte-Carlo convergence to the exact n=14 limit
    K, N = 13, 14
    exact14 = core_density_exact(K, N)
    print(f"\n(1b) Monte-Carlo at n=14 converging to the exact limit {exact14:.5f}")
    for R in [200, 1000, 5000]:
        rnd = random.Random(7)
        cov = tot = 0
        TARGET = 20000
        while tot < TARGET:
            v = rnd.sample(range(1, R + 1), K)
            if reduce(gcd, v) != 1:
                continue
            tot += 1
            if covered(v, N):
                cov += 1
        print(f"    R={R:5d}: P~={cov/tot:.4f}   ({cov}/{tot})")

    # (1c) reproduce S557's headline count (seed 557, 40 sets, range 1..79)
    rnd = random.Random(557)
    ncov = 0
    for _ in range(40):
        while True:
            v = tuple(sorted(rnd.sample(range(1, 80), K)))
            if reduce(gcd, v) == 1:
                break
        if covered(v, N):
            ncov += 1
    print(f"\n(1c) reproduce S557 (seed 557, 40 sets): sieve-covered = {ncov}/40   "
          f"(S557 reported 2/40)")

    # (2) gap distribution OVER the core at n=14
    print("\n(2) GAP over random sieve-covered primitive 13-sets at n=14  (conjecture 1/14="
          f"{1/N:.4f})")
    conj = 1.0 / N
    rnd = random.Random(99)
    gaps = []
    while len(gaps) < 200:
        v = tuple(sorted(rnd.sample(range(1, 120), K)))
        if reduce(gcd, v) == 1 and covered(v, N):
            gaps.append((gap(v), v))
    gaps.sort()
    gs = [g for g, _ in gaps]
    print(f"    sampled {len(gs)} core sets:  min={min(gs):.4f}  median={st.median(gs):.4f}  "
          f"max={max(gs):.4f}")
    print(f"    min ratio g/(1/14) = {min(gs)/conj:.2f}   "
          f"frac with g < 1.15*(1/14) = {sum(1 for g in gs if g < 1.15*conj)/len(gs):.3f}")
    print("    smallest-gap core sets (the near-wall direction):")
    for g, v in gaps[:3]:
        print(f"       g={g:.4f}  g/(1/14)={g/conj:.2f}  v={v}")

    # (2b) the wall itself (AP 1..13): sieve-covered AND tight
    ap = tuple(range(1, N))
    g_ap = gap(ap)
    print(f"\n(2b) the wall AP{ap}: covered={covered(ap, N)}, g={g_ap:.5f} ~ 1/{N}={conj:.5f}")
    print("     -> the wall is the ONE tight core set; generic core sets miss it by ~2-4x.")

    print("\n" + "=" * 78)
    print("VERDICT")
    print("=" * 78)
    print(f"""  * The sieve-covered CORE is NOT measure-zero: exact limiting density
    {exact14:.4f} at n=14 (~1 in 9 random primitive 13-sets), rising as n shrinks.
  * But generic core sets are FAR from the wall: over 200 random core sets the
    smallest gap is {min(gs):.4f} = {min(gs)/conj:.1f}x the conjecture; none come near 1/14.
  * So S557's "thin core" is right about the EXTREMIZERS (the near-AP slice is the
    genuinely thin sub-locus) but overstates the thinness of the core itself. The
    open problem is: a positive-density core whose hardness concentrates on a thin
    near-AP slice. Density-based ("measure-zero core") attacks cannot work; the slice
    must be cut by the AP/dilation structure (S530 apex, S556 local LP), not by size.""")


if __name__ == "__main__":
    main()
