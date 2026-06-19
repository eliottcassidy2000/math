"""
Theory-guided search.

Lower bound: M = p/q (lowest terms), M > 1/(k+1) => g >= 1/(q(k+1)),
and q <= 2 max(S). So g*k^2 >= k^2/(q(k+1)) ~ k/q ~ k/(2 max S).
=> to get g*k^2 small (level a large) we MUST use large max speed.

Construction idea ("killer" speeds): Take AP {1,...,k-r} as a "dense block" that
forces the optimal t near 1/(k+1)-ish, then add r large killer speeds that are
multiples of a big modulus designed so ||killer * t*|| stays >= M.

Cleaner known idea: the "Goddyn-style" / dilation construction.
We test:
 (A) {1,2,...,k} with the LAST element replaced by a large L (single killer).
 (B) {1,...,k-r} + r killers each ~ j*(k+1) (multiples of k+1, shifted).
 (C) dilation: take a small extremal set, multiply by d, add 1 to restore gcd.

We measure M exactly and report g*k^2 and level a.
"""
import sys
from fractions import Fraction
from math import gcd
from functools import reduce
sys.path.insert(0, ".")
from lrc_maxmin import M_exact, lcm_list
from primorial_family import level_a


def setgcd(S):
    return reduce(gcd, S)


def evalset(S, k=None):
    S = sorted(set(int(x) for x in S))
    if k is None:
        k = len(S)
    if len(S) != k or setgcd(S) != 1 or min(S) <= 0:
        return None
    M, t = M_exact(S)
    floor = Fraction(1, k + 1)
    g = M - floor
    return dict(S=S, M=M, t=t, floor=floor, g=g, tight=(g == 0),
                gk2=g * k * k, a=(level_a(M, k) if g > 0 else None))


def fmt(r):
    if r is None:
        return "INVALID"
    a = r['a']
    return (f"M={r['M']} (~{float(r['M']):.10f}) g*k^2={float(r['gk2']):.6f} "
            f"a={float(a):.4f} tight={r['tight']} maxS={max(r['S'])}")


def best_over(gen, k):
    best = None
    for S in gen:
        r = evalset(S, k)
        if r is None or r['tight']:
            continue
        if r['g'] <= 0:
            continue
        if best is None or r['gk2'] < best['gk2']:
            best = r
    return best


if __name__ == "__main__":
    import argparse
    for k in [30, 31, 36, 40, 48, 60]:
        print(f"==== k={k}  (sqrt(k)={k**0.5:.3f}) ====")
        # (A) single killer: {1..k-1} + L, L large
        def genA():
            base = list(range(1, k))
            for L in range(k, 40 * k):
                yield base + [L]
        rA = best_over(genA(), k)
        print("  (A) single killer {1..k-1}+L:", fmt(rA))
