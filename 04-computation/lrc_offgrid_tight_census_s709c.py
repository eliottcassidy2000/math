"""
monad-explorer-2026-06-06-S709c — Are there PRIMITIVE tight configs with a multiple of n?

s709b established: grid-witness exists <=> no speed ≡0 (mod n); grid-witness => M>=1/n
(via t=1/n: ||v_i/n|| = min(r,n-r)/n >= 1/n when r=v_i mod n != 0). Hence:

   LRC(n) is TRIVIAL for any primitive speed set with NO multiple of n.
   => every potential LRC counterexample MUST contain a multiple of n  (an "off-grid" speed).

If a runner v_j ≡ 0 (mod n), the whole n-grid {m/n} fails to witness loneliness
(v_j sits ON the observer there), so M must be attained OFF the n-grid. The worry-set
question becomes: do PRIMITIVE TIGHT configs (M=1/n exactly) with a multiple of n exist?
If yes, they are the genuinely hard 'off-grid' boundary the signed/2n-1 theory could touch.
If no (in a window), tight configs are all on-grid and the floor is a pure mod-n phenomenon.

Census: exhaustive over (n-1)-subsets of [1..W], primitive (gcd=1). Classify tight configs
by (has multiple of n) and (M attained at denominator n).
"""
import sys
from itertools import combinations
from math import gcd
from fractions import Fraction as F


def pinch_denoms(V):
    ds = set(V)
    L = len(V)
    for i in range(L):
        for j in range(i + 1, L):
            ds.add(V[i] + V[j]); ds.add(abs(V[i] - V[j]))
    ds.discard(0)
    return ds


def M_fast(V):
    best_n, best_d = -1, 1
    arg = (1, 0)
    for d in pinch_denoms(V):
        for m in range(1, d):
            mn = d
            for v in V:
                r = (v * m) % d
                rr = d - r
                x = r if r < rr else rr
                if x < mn:
                    mn = x
                    if mn == 0:
                        break
            if mn * best_d > best_n * d:
                best_n, best_d = mn, d
                arg = (d, m)
    g = gcd(best_n, best_d)
    return (best_n // g, best_d // g), arg


def gcd_list(V):
    g = 0
    for v in V:
        g = gcd(g, v)
    return g


print = print
def flush(): sys.stdout.flush()

windows = {4: 6*4, 5: 5*5, 6: 4*6, 7: 3*7, 8: 2*8+8}
for n in sorted(windows):
    W = windows[n]
    thr = F(1, n)
    tot_tight = 0
    with_mult = []     # primitive tight configs containing a multiple of n
    offgrid = []       # tight configs whose M is NOT attained at denominator dividing n
    for V in combinations(range(1, W + 1), n - 1):
        if gcd_list(V) != 1:
            continue
        (Mn, Md), (d, m) = M_fast(V)
        if F(Mn, Md) != thr:
            continue
        tot_tight += 1
        # reduced denominator of t* = m/d
        td = d // gcd(d, m)
        has_mult = any(v % n == 0 for v in V)
        on_grid_n = (n % td == 0)   # t* denominator divides n
        if has_mult:
            with_mult.append((V, td))
        if not on_grid_n:
            offgrid.append((V, td))
    print(f" n={n} (W={W}): {tot_tight} primitive tight configs")
    print(f"    with a multiple of n: {len(with_mult)}"
          + ("" if not with_mult else f"  e.g. {with_mult[:6]}"))
    print(f"    M NOT on n-grid (t* denom ∤ n): {len(offgrid)}"
          + ("" if not offgrid else f"  e.g. {offgrid[:6]}"))
    flush()

print()
print("Lemma spot-check: grid-witness(S) exists  <=>  no v in S with v % n == 0")
import random
random.seed(1)
bad = 0
for _ in range(200000):
    n = random.randint(3, 12)
    k = random.randint(1, min(n - 1, 6))
    V = random.sample(range(1, 6 * n), k)
    gw = any(all((m * v) % n != 0 for v in V) for m in range(1, n))
    nomult = all(v % n != 0 for v in V)
    if gw != nomult:
        bad += 1
print(f"   200000 random (n,S): grid-witness == no-multiple-of-n  mismatches = {bad}")
