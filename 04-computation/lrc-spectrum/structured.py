"""
Structured exhaustive over minimizer-shaped sets:
  S = (a near-prefix of {1,2,...} with up to G gaps) + (r killers from a bounded range)
total size k. This matches observed sigma_2 minimizers:
  k=4 {1,2,3,8}; k=5 {1,2,3,4,10}; k=7 {1,2,3,4,5,7,18}; k=9 {1..8,18}.
We exhaust gaps (positions removed from a length-(k-r+G) prefix) and killers (from a
window), compute EXACT M, track smallest M (largest level a). Polynomial in k.
"""
import sys
from fractions import Fraction
from math import gcd
from functools import reduce
from itertools import combinations
sys.path.insert(0, ".")
from fast_M import M_exact_fast
from primorial_family import level_a


def setgcd(S):
    return reduce(gcd, S)


def Mval(S):
    return M_exact_fast(sorted(set(S)))[0]


def search(k, r_max=2, G_max=2, killer_hi_mult=4, verbose=False):
    """r killers (1..r_max), G gaps (0..G_max). Block = first (k-r) survivors of a
    prefix of length (k-r+G); killers from [k, killer_hi_mult*k]."""
    floor = Fraction(1, k + 1)
    best = None
    for r in range(1, r_max + 1):
        block_size = k - r
        killer_lo = block_size + 1
        killer_hi = killer_hi_mult * k
        killer_range = list(range(killer_lo, killer_hi + 1))
        for G in range(0, G_max + 1):
            prefix_len = block_size + G
            if prefix_len < block_size:
                continue
            prefix = list(range(1, prefix_len + 1))
            # choose G positions to remove (gaps), but never remove 1 (keep gcd handle)
            removable = prefix[1:]  # keep element 1
            for gaps in combinations(removable, G):
                block = [x for x in prefix if x not in gaps]
                if len(block) != block_size:
                    continue
                for killers in combinations(killer_range, r):
                    S = block + list(killers)
                    if len(set(S)) != k:
                        continue
                    if setgcd(S) != 1:
                        continue
                    M = Mval(S)
                    if M <= floor:
                        continue
                    if best is None or M < best[0]:
                        best = (M, sorted(S), (r, G, gaps, killers))
                        if verbose:
                            print(f"   new best M={M} (a={float(level_a(M,k)):.3f}) S={sorted(S)}", flush=True)
    return best, floor


if __name__ == "__main__":
    import sympy
    ks = [int(x) for x in sys.argv[1:]] or [10, 11, 12, 13, 14]
    for k in ks:
        # bound killer range to keep it tractable; r_max=2 for larger k
        rmax = 2 if k > 11 else 2
        gmax = 2
        best, floor = search(k, r_max=rmax, G_max=gmax,
                             killer_hi_mult=4 if k <= 14 else 3)
        om = len(sympy.factorint(k - 1))
        if best:
            M, S, info = best
            g = M - floor
            print(f"k={k} sqrt(k)={k**0.5:.3f} omega(k-1)={om}: M={M} (~{float(M):.8f}) "
                  f"g*k^2={float(g*k*k):.5f} a={float(level_a(M,k)):.4f}")
            print(f"   S={S}  info={info}", flush=True)
        else:
            print(f"k={k}: none", flush=True)
