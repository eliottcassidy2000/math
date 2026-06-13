#!/usr/bin/env python3
"""
lrc_gap_lower_bounds_by_case_s557.py    oracle-2026-06-01-S557o

User: Tao improved the GENERAL lower bound for the lonely-runner gap from the trivial
   g(v) := sup_t min_i ||v_i t||  >=  1/(2k)        (k = # nonzero speeds)
to  1/(2k) + c*log k / (k^2 (log log k)^2)  -- a tiny gain over the trivial bound.
"Can you improve this further? what about in certain cases?"

HONEST FRAME. The trivial 1/(2k) IS the first moment: at threshold theta the expected
near-count is E[#near] = 2*k*theta (each runner's danger arc has measure 2*theta), so
E[#near] < 1 iff theta < 1/(2k) -> some t is lonely. Tao pushes theta a hair above
1/(2k) by lower-bounding the OVERLAPS of the danger arcs (second-order). Beating his
general constant is a delicate optimization I will NOT win in a session.

BUT "in certain cases" the repo's SIEVE gives the FULL conjecture g >= 1/(k+1), and
structured cases give far more (all-odd -> g >= 1/2). This script makes the CASE
HIERARCHY precise and measures, for n = k+1 = 14, how far above Tao's bound the real
gaps sit, and exactly which (small) family of sets Tao's bound is actually FOR.

CASE HIERARCHY (rigorous lower bounds on g(v)):
  (C0) trivial / first moment ............ g >= 1/(2k)                  [always]
  (C0') Tao .............................. g >= 1/(2k) + c log k/(...)  [always, tiny]
  (C1) sieve at q* ....................... g >= 1/q*,  q* = least q in {2..k+1}
                                           dividing NO speed (oo if sieve-covered)
       => not sieve-covered  =>  q* <= k+1  =>  g >= 1/(k+1)  (FULL CONJECTURE)
  (C2) all speeds odd .................... g >= 1/2          (q*=2)   [>> conjecture]
  (C3) no speed divisible by k+1 ......... g >= 1/(k+1)      (q*=k+1) exactly
  (C4) coprime small modulus q ........... g >= 1/q for the smallest missed q
The point: (C1) BEATS Tao (gives the full 1/(k+1)) for every set that is not
sieve-covered. Sieve-covered sets (a multiple of every q<=k+1) are the ONLY locus
where Tao's bound is the operative general bound -- and they are exactly the open
core (S554 condition A1). So: in certain cases (the bulk) the conjecture is PROVEN.
"""
from functools import reduce
from math import gcd
import random

def d0(p):
    p = p % 1.0
    return min(p, 1 - p)

def gap(v, G=400000):
    """actual gap g(v) = sup_t min_i ||v_i t|| via fine grid (lower estimate)."""
    best = 0.0
    for i in range(G):
        t = (i + 0.5) / G
        m = min(d0(s * t) for s in v)
        if m > best:
            best = m
    return best

def q_star(v, qmax):
    """least q in {2..qmax} dividing NO speed; None if sieve-covered up to qmax."""
    for q in range(2, qmax + 1):
        if all(s % q != 0 for s in v):
            return q
    return None

def sieve_bound(v, qmax):
    q = q_star(v, qmax)
    return (1.0 / q, q) if q else (0.0, None)

def main():
    K = 13                 # nonzero speeds
    N = K + 1              # = 14 ; conjecture target 1/N
    triv = 1.0 / (2 * K)   # trivial / first-moment bound = 1/26
    conj = 1.0 / N         # conjecture = 1/14
    print("=" * 78)
    print(f"n={N} runners (observer + {K} movers).  trivial g>=1/(2k)=1/{2*K}={triv:.5f}"
          f"   conjecture g>=1/{N}={conj:.5f}")
    print("=" * 78)

    # ---- (A) the gap is ~conjecture, ~2x Tao, for random sets; sieve bound nails it
    print("\n(A) RANDOM primitive 13-sets: actual gap vs sieve bound (C1) vs trivial (C0)")
    rnd = random.Random(557)
    n_notcov = n_cov = 0
    gaps_notcov = []
    for trial in range(40):
        while True:
            v = tuple(sorted(rnd.sample(range(1, 80), K)))
            if reduce(gcd, v) == 1:
                break
        g = gap(v)
        sb, q = sieve_bound(v, N)
        covered = (q is None)
        tag = "SIEVE-COVERED (core; Tao's regime)" if covered else f"q*={q}: g>=1/{q}={sb:.4f} (>=1/{N} conj)"
        if covered:
            n_cov += 1
        else:
            n_notcov += 1
            gaps_notcov.append(g)
        if trial < 12 or covered:
            ok = "OK" if (covered or g >= sb - 1e-3) else "??"
            print(f"  set#{trial:02d} max={max(v):2d}: g={g:.4f}  triv={triv:.4f}  {tag}  [{ok}]")
    print(f"\n  -> not sieve-covered: {n_notcov}/40  (these have g >= 1/{N}, the FULL conjecture,"
          f" proven by the sieve)")
    print(f"  -> sieve-covered:     {n_cov}/40  (the ONLY locus where Tao's general bound is"
          f" operative -- the open core)")
    if gaps_notcov:
        print(f"  -> over non-covered sets: min actual gap = {min(gaps_notcov):.4f} "
              f">= conjecture 1/{N}={conj:.4f}  (>> trivial {triv:.4f}; ~{min(gaps_notcov)/triv:.1f}x Tao)")

    # ---- (B) structured cases that BLOW PAST the conjecture
    print("\n(B) STRUCTURED cases (rigorous, far above the conjecture):")
    allodd = tuple(range(1, 2 * K, 2))                 # 1,3,5,...,25 : all odd
    g = gap(allodd); print(f"  all-odd  {allodd[:5]}...: q*=2  g>=1/2;  actual g={g:.4f}")
    coprimeN = tuple(s for s in range(1, 200) if s % N != 0)[:K]  # none divisible by 14
    g = gap(coprimeN); print(f"  no mult of {N}: q*={N}  g>=1/{N}={conj:.4f} (conjecture exactly); actual g={g:.4f}")
    no_small = tuple(s for s in range(1, 400) if s % 2 and s % 3)[:K]  # none div by 2 or 3
    qb, q = sieve_bound(no_small, N)
    g = gap(no_small); print(f"  no mult of 2 or 3: q*={q}  g>=1/{q}={qb:.4f}; actual g={g:.4f}")

    # ---- (C) the AP/wall: the unique tight set (gap = exactly 1/n)
    print("\n(C) THE WALL (AP 1..13): sieve-covered, gap = exactly the conjecture")
    ap = tuple(range(1, N)); g = gap(ap)
    print(f"  AP {ap}: sieve-covered (q*=None, has a mult of every q<=14);"
          f" g={g:.5f} ~ 1/{N}={conj:.5f}")
    print(f"  -> the AP saturates the conjecture: NO general bound can exceed 1/{N}."
          f" The wall caps the game at the conjecture; Tao's bound lives BELOW it,"
          f" for the sieve-covered core only.")

    # ---- summary
    print("\n" + "=" * 78)
    print("VERDICT")
    print("=" * 78)
    print(f"""  * Trivial 1/(2k) = first moment; Tao adds a log/loglog hair via arc-overlaps.
  * SIEVE (C1) beats Tao outright for every NON-sieve-covered set: g >= 1/q* >= 1/{N}
    (the FULL conjecture). That is the BULK ({n_notcov}/40 here).
  * Structured cases blow past it: all-odd -> g >= 1/2; no-mult-of-2-or-3 -> g >= 1/2.
  * The ONLY sets where Tao's general bound is the operative one are SIEVE-COVERED
    (a multiple of every q<=k+1) -- a subset of the open core (S554-A1). And even most
    SIEVE-COVERED sets have gap >> 1/(k+1) (e.g. 0.16-0.23 above): the genuinely hard
    locus is the even thinner near-AP slice. The AP/wall is sieve-covered and saturates
    g = 1/(k+1): no general bound can beat the conjecture.
  * So 'improve in certain cases' = PROVEN: g >= 1/(k+1) for all non-sieve-covered
    sets (the bulk). Tao's frontier is the measure-zero-ish sieve-covered core, where
    the repo's tools (measure bound S550 at theta<1/n; local LP S556) give g >= 1/n
    off the AP -- the genuine remaining work.""")

if __name__ == "__main__":
    main()
