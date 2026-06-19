"""
Structured search tuned for larger k. Minimizer shape that produced a=3 at k=7,13:
   block = {1,...,k} with ONE interior element removed (gap at position g),
   plus ONE killer X (single killer), size k.  (r=1, G in {0,1,2})
Also r=2 killers with G in {0,1} and a bounded killer window.
EXACT M. Polynomial, fast.
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


def search(k, G_max=2, r_max=2, killer_hi_mult=4):
    floor = Fraction(1, k + 1)
    best = None
    for r in range(1, r_max + 1):
        block_size = k - r
        for G in range(0, G_max + 1):
            prefix_len = block_size + G
            prefix = list(range(1, prefix_len + 1))
            removable = prefix[1:]  # keep 1
            for gaps in combinations(removable, G):
                block = [x for x in prefix if x not in gaps]
                if len(block) != block_size:
                    continue
                klo = max(block) + 1
                khi = killer_hi_mult * k
                krange = list(range(klo, khi + 1))
                if r == 1:
                    for X in krange:
                        S = block + [X]
                        if setgcd(S) != 1:
                            continue
                        M = Mval(S)
                        if M <= floor:
                            continue
                        if best is None or M < best[0]:
                            best = (M, sorted(S), (r, G, gaps, (X,)))
                else:
                    # r==2: limit killer window to keep tractable
                    krange2 = krange[: min(len(krange), 3 * k)]
                    for killers in combinations(krange2, 2):
                        S = block + list(killers)
                        if setgcd(S) != 1:
                            continue
                        M = Mval(S)
                        if M <= floor:
                            continue
                        if best is None or M < best[0]:
                            best = (M, sorted(S), (r, G, gaps, killers))
    return best, floor


if __name__ == "__main__":
    import sympy
    ks = [int(x) for x in sys.argv[1:]] or [30, 31, 42, 43, 60, 61]
    for k in ks:
        # for larger k, single killer dominates; keep r_max=2 only if k small enough
        rmax = 2 if k <= 22 else 1
        gmax = 3 if k <= 22 else 2
        best, floor = search(k, G_max=gmax, r_max=rmax,
                             killer_hi_mult=4 if k <= 35 else 3)
        om = len(sympy.factorint(k - 1))
        M, S, info = best
        g = M - floor
        print(f"k={k} sqrt(k)={k**0.5:.3f} omega(k-1)={om}: M={M} (~{float(M):.9f}) "
              f"g*k^2={float(g*k*k):.6f} a={float(level_a(M,k)):.4f}")
        print(f"   S={S}  info={info}", flush=True)
