"""
Exhaustive sigma_2 for small k: find the non-tight gcd-1 k-set with smallest M
(closest to floor), over speeds bounded by V. Report the minimizer + structure.
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


def run(k, V):
    floor = Fraction(1, k + 1)
    best = None
    for S in combinations(range(1, V + 1), k):
        if setgcd(S) != 1:
            continue
        M, t = M_exact_fast(S)
        if M <= floor:
            continue  # tight (==) or (shouldn't be <)
        if best is None or M < best[0]:
            best = (M, S, t)
    return best, floor


if __name__ == "__main__":
    for k, V in [(4, 14), (5, 16), (6, 18), (7, 18), (8, 20)]:
        best, floor = run(k, V)
        if best:
            M, S, t = best
            g = M - floor
            print(f"k={k} V={V}: sigma_2 minimizer S={S}")
            print(f"   M={M} (~{float(M):.8f}) floor={floor} g*k^2={float(g*k*k):.5f} "
                  f"a={float(level_a(M,k)):.4f} t*={t}", flush=True)
        else:
            print(f"k={k}: none", flush=True)
