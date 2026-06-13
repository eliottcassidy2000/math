#!/usr/bin/env python3
"""
monad-explorer-2026-06-06-S702
Characterize the LOWER-HALF deletion ladder M(AP_n \ {a}) for a < n/2.

S701 proved the UPPER half exactly: M(AP_n\{a}) = 1/a for n/2 <= a <= n-1.
S701 only *verified* the lower half a < n/2 and (wrongly) called it a "2/(odd)"
ladder -- its own data has n=12,a=3 -> 1/8 = 2/16 (even denom).  This script
gathers exhaustive exact data for the lower half (both even AND odd n), records
M=p/d reduced, the optimal time t, the BINDING runners (argmin), and the gap
1/M, and hunts for the true closed form.

Exact M via the standard breakpoint set (half-integer maxima + pairwise pinches).
"""
from fractions import Fraction
from math import gcd


def frac(x):
    r = x % 1
    return min(r, 1 - r)


def candidate_times(V):
    c = set()
    for v in V:
        v = abs(v)
        for k in range(2 * v):
            c.add(Fraction(2 * k + 1, 2 * v) % 1)
    for i in range(len(V)):
        for j in range(len(V)):
            for s in (1, -1):
                d = V[i] + s * V[j]
                if d:
                    d = abs(d)
                    for k in range(d + 1):
                        c.add(Fraction(k, d) % 1)
    c.discard(Fraction(0))
    return c


def M_exact(V):
    V = list(V)
    best = Fraction(0)
    bt = None
    for t in candidate_times(V):
        mn = min(frac(v * t) for v in V)
        if mn > best:
            best, bt = mn, t
    return best, bt


def binders(V, t, M):
    return [v for v in V if frac(v * t) == M]


if __name__ == "__main__":
    print("=" * 78)
    print("S702: lower-half deletion ladder M(AP_n \\ {a}), a < n/2")
    print("=" * 78)
    # store delta = den(M) - n  where M reduces to p/d
    for n in range(4, 31):
        AP = list(range(1, n))
        half = n / 2
        rows = []
        for a in range(1, n):
            if a >= half:  # only lower half (strict) of interest, but compute middle too
                continue
            R = [v for v in AP if v != a]
            M, t = M_exact(R)
            d = M.denominator
            p = M.numerator
            bnd = binders(R, t, M)
            rows.append((a, M, t, p, d, bnd))
        if not rows:
            continue
        print(f"\n--- n={n}  (1/n=1/{n}, 2/n=2/{n}) ---")
        for (a, M, t, p, d, bnd) in rows:
            # express as 2/x if numerator 2
            asg = f"{p}/{d}"
            twoform = f"  =2/{Fraction(2,1)/M}" if p == 2 else ""
            inv = Fraction(1,1)/M
            print(f"  a={a:>2}: M={asg:>8}={float(M):.5f}  1/M={inv}  t={t}  binders={bnd}")
