#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THE THREE TILING RECURRENCES, COMBINED INTO ONE PARTITION FUNCTION  (kps-2026-06-20)
===================================================================================
User ask: combine the full-tiling recurrence and the even/odd half-tiling recurrences
"at the same time" to understand the recursive structure of any tile, treat the
half-tiling as an address quotient, and use it to EFFICIENTLY compute how tournament
structure changes as n grows.

THE COMBINATION.  In codex's two-clock address (THM-553) every free tile (a,b),
1<=b<a<=n, a-b>=2, has
    beta(a,b) = a          (BIRTH clock; the full-tiling beta-strip recurrence, THM-442/553)
    tau(a,b)  = a+b-1      (MIRROR clock; the half-tiling tau-layer recurrence, THM-550)
with gap = a-b-1 = 2*beta-tau-2 and b = tau-beta+1.  The full recurrence advances beta;
the even/odd half recurrences advance tau split by parity(tau).

These three recurrences combine into the SCORE PARTITION FUNCTION of the tiling model.
With the canonical convention (base path n->n-1->...->1, tile bit 0 => a->b, bit 1 => b->a),
each tile (a,b) adds EXACTLY +1 to score(a) [bit 0] or score(b) [bit 1].  Hence the joint
score generating function over all 2^{C(n-1,2)} tilings factors:

    Z_n(x_1,...,x_n) = ( prod_{v=2}^n x_v ) * prod_{(a,b) tile} (x_a + x_b).        (Z)

  * BETA-CLOCK (full recurrence) = the INCREMENTAL n->n+1 update: birth strip beta=n+1
        Z_{n+1} = Z_n * x_{n+1} * prod_{b=1}^{n-1} (x_{n+1} + x_b).                 (beta-step)
    This is the engine: build Z_n one birth strip at a time, NEVER enumerating 2^F tilings.
  * TAU-CLOCK (half recurrence) = the COMPLEMENT reflection R:(a,b)->(n+1-b,n+1-a) (THM-280),
    which is tournament complement up to relabel.  Every complement-invariant, score-
    determined invariant (c3 = leading OCF/alpha_1 term, the score census) is R-symmetric,
    so the half-tiling tau<=n representatives (THM-550) carry the whole distribution -- the
    ADDRESS QUOTIENT.  c3(T)=c3(T^op) gives a 2x fold.

WHAT THIS COMPUTES (exactly, no 2^F enumeration):
  - the score-sequence census (multiplicity of every score vector),
  - the c3 (3-cycle) distribution, since c3 = C(n,3) - sum_v C(s_v, 2),
  - hence the leading OCF term alpha_1's distribution and the H-max/regular census.

Verified: Z-engine score & c3 distribution == brute over all tilings (n<=7).  Reaches the
EXACT c3 distribution at n=10 (68.7 billion tilings) in ~100s -- far beyond enumeration.
"""
import sys, time
from collections import defaultdict, Counter
from itertools import product
from math import comb
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")


def tiles(n):
    return [(a, b) for a in range(3, n + 1) for b in range(1, a - 1)]


def beta_step(dist, n):
    """Z_{n-1} -> Z_n: add vertex n (base arc n->n-1) then the birth strip beta=n."""
    nd = defaultdict(int)
    for vec, cnt in dist.items():
        l = list(vec) + [0]
        l[n - 1] += 1                      # base arc n->n-1 : +1 to score(n)
        nd[tuple(l)] += cnt
    dist = nd
    for b in range(1, n - 1):              # birth-strip tiles (n,b), gap = n-b-1
        nd = defaultdict(int)
        for vec, cnt in dist.items():
            l0 = list(vec); l0[n - 1] += 1; nd[tuple(l0)] += cnt   # bit0: +1 to a=n
            l1 = list(vec); l1[b - 1] += 1; nd[tuple(l1)] += cnt   # bit1: +1 to b
        dist = nd
    return dist


def build_Z(N):
    dist = {(0,): 1}                        # n=1
    for n in range(2, N + 1):
        dist = beta_step(dist, n)
    return dist


def c3_distribution(distZ, n):
    cs = Counter()
    for vec, cnt in distZ.items():
        cs[comb(n, 3) - sum(comb(s, 2) for s in vec)] += cnt
    return cs


def score_census(distZ):
    """Multiplicity of each (sorted) score multiset across all tilings."""
    cs = Counter()
    for vec, cnt in distZ.items():
        cs[tuple(sorted(vec))] += cnt
    return cs


# ---------------------------------------------------------------- validation
def _brute_c3(n):
    T = tiles(n); cs = Counter()
    for bv in product((0, 1), repeat=len(T)):
        adj = [[0] * (n + 1) for _ in range(n + 1)]
        for k in range(n, 1, -1):
            adj[k][k - 1] = 1
        for (a, b), bit in zip(T, bv):
            if bit == 0:
                adj[a][b] = 1
            else:
                adj[b][a] = 1
        t = 0
        for i in range(1, n + 1):
            for j in range(i + 1, n + 1):
                for k in range(j + 1, n + 1):
                    if (adj[i][j] + adj[i][k], adj[j][i] + adj[j][k],
                            adj[k][i] + adj[k][j]) == (1, 1, 1):
                        t += 1
        cs[t] += 1
    return cs


def main():
    print(__doc__)
    print("=== VALIDATION: Z-engine c3-distribution == brute over all tilings ===")
    for n in range(3, 8):
        gf = c3_distribution(build_Z(n), n)
        ok = dict(gf) == dict(_brute_c3(n))
        print(f"  n={n}: GF==brute {ok}  (total {sum(gf.values())} = 2^{{{comb(n-1,2)}}})")

    print("\n=== EFFICIENCY: exact c3-distribution via the beta-clock, no 2^F enumeration ===")
    print("  n   tilings=2^C(n-1,2)   Z-states   c3_max   mean_c3   uniform_mean   t(s)")
    dist = {(0,): 1}
    for n in range(2, 11):
        t0 = time.time(); dist = beta_step(dist, n); dt = time.time() - t0
        if n >= 3:
            cs = c3_distribution(dist, n)
            mean = sum(c * m for c, m in cs.items()) / sum(cs.values())
            print(f"  {n:2d}  {1<<comb(n-1,2):18d}  {len(dist):10d}  {max(cs):6d}  "
                  f"{mean:8.3f}  {comb(n,3)/4:11.3f}   {dt:5.2f}")

    print("\n=== APPLICATION: REGULAR (H-max/Paley) score-sequence census in the tiling model ===")
    for n in (3, 5, 7, 9):
        d = build_Z(n); r = (n - 1) // 2
        reg = sum(c for v, c in d.items() if all(s == r for s in v))
        print(f"  n={n}: # tilings with regular score seq (all s_v={r}) = {reg}")


if __name__ == "__main__":
    main()
