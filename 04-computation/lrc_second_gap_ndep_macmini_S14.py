#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S14 -- HYP-4422: THE n-DEPENDENCE OF THE SECOND GAP.

codex-S573 found the second gap (1/n, 2/(2n-1)) is NONEMPTY at n=7 (M=5/33) and
n=8 (M=3/23).  So (G) -- gap-emptiness at n=13 -- is n-SPECIFIC, not universal.
This sweep finds WHICH n have an empty second gap and correlates with arithmetic
(n prime? 2n-1 prime/composite? factorization?), to isolate what makes n=13
special -- the (G) proof mechanism.

For LRC(n): n-1 nonzero speeds, floor 1/n, second value 2/(2n-1); gap =
(1/n, 2/(2n-1)).  Search (n-1)-speed primitive families (bounded height) for a
member with M in the gap.
"""
import itertools, random
from fractions import Fraction as F
from math import gcd
from sympy import isprime, factorint

def exact_M(W):
    dens = set()
    for v, w in itertools.combinations(W, 2):
        dens.add(v + w)
        if v != w:
            dens.add(abs(v - w))
    for v in W:
        dens.add(2 * v)
    best = F(0); seen = set()
    for s in dens:
        if s == 0:
            continue
        for j in range(1, s):
            t = F(j, s)
            if t in seen:
                continue
            seen.add(t)
            mv = min(abs(v * t - round(v * t)) for v in W)
            if mv > best:
                best = mv
    return best

def float_M(W):
    dens = set()
    for v, w in itertools.combinations(W, 2):
        dens.add(v + w)
        if v != w:
            dens.add(abs(v - w))
    for v in W:
        dens.add(2 * v)
    best = 0.0
    for s in dens:
        if s == 0:
            continue
        for j in range(1, s):
            t = j / s
            mv = min(abs(v * t - round(v * t)) for v in W)
            if mv > best:
                best = mv
    return best

def primitive(W):
    g = 0
    for v in W:
        g = gcd(g, v)
    return tuple(sorted(v // g for v in W)) if g > 1 else tuple(sorted(W))

def search_gap(n, maxspeed, n_random=200000, seed=0):
    """find (n-1)-speed primitive families with M in (1/n, 2/(2n-1))."""
    lo = F(1, n); hi = F(2, 2 * n - 1)
    flo, fhi = float(lo), float(hi)
    k = n - 1
    members = []
    minM = None
    rng = random.Random(seed)
    # exhaustive for small k, random for large
    if k <= 5 and maxspeed <= 24:
        it = itertools.combinations(range(1, maxspeed + 1), k)
    else:
        def gen():
            for _ in range(n_random):
                yield tuple(sorted(rng.sample(range(1, maxspeed + 1), k)))
        it = gen()
    seen = set()
    for W in it:
        Wp = primitive(W)
        if Wp in seen or len(set(Wp)) != k:
            continue
        seen.add(Wp)
        fm = float_M(Wp)
        if flo - 1e-6 < fm < fhi + 1e-6:      # near the gap
            M = exact_M(Wp)
            if lo < M < hi:
                members.append((M, Wp))
                if minM is None or M < minM:
                    minM = M
    return members, minM, lo, hi

def main():
    print("=" * 90)
    print("THE n-DEPENDENCE OF THE SECOND GAP (1/n, 2/(2n-1))")
    print("=" * 90)
    print(f"  {'n':>3} {'floor':>6} {'2nd':>7} {'gap':>18} {'n prime':>8} {'2n-1':>6} {'2n-1 fact':>12} "
          f"{'#members':>9} {'min M':>10}")
    rows = []
    # bounds: small n exhaustive-ish, large n random with modest height
    config = {4: 20, 5: 20, 6: 22, 7: 20, 8: 22, 9: 22, 10: 24, 11: 26, 12: 28, 13: 30}
    for n in range(4, 14):
        ms = config[n]
        nr = 300000 if n >= 8 else 0
        members, minM, lo, hi = search_gap(n, ms, n_random=nr, seed=n)
        f2n1 = factorint(2 * n - 1)
        factstr = "*".join(f"{p}^{e}" if e > 1 else str(p) for p, e in f2n1.items())
        rows.append((n, len(members), minM, members))
        print(f"  {n:>3} {'1/'+str(n):>6} {'2/'+str(2*n-1):>7} "
              f"{f'({float(lo):.4f},{float(hi):.4f})':>18} {str(isprime(n)):>8} {2*n-1:>6} "
              f"{factstr:>12} {len(members):>9} {str(minM) if minM else '--':>10}")
    print()
    print("  MEMBERS by n (up to 3 each):")
    for n, cnt, minM, members in rows:
        if members:
            for M, W in sorted(members)[:3]:
                print(f"    n={n}: M={M}={float(M):.5f}  W={list(W)}")
    print()
    print("  PATTERN: empty-gap n vs nonempty-gap n -- correlate with 2n-1 factorization")
    empty = [n for n, c, m, mem in rows if c == 0]
    nonempty = [n for n, c, m, mem in rows if c > 0]
    print(f"    EMPTY second gap (bounded search): n = {empty}")
    print(f"    NONEMPTY second gap: n = {nonempty}")
    print(f"    (n=13: 2n-1=25=5^2; is it empty? the (G) target)")

if __name__ == "__main__":
    main()
