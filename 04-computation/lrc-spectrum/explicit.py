"""
Explicit family tests, exact M, fast. Battery per k:

(1) SINGLE-KILLER baseline: {1,...,k-1, X} with X = 2k+1 ... gives M=2/(2k+1), a=2.
    (We confirm and scan X.)

(2) PRIMORIAL/level-a construction via t* = m/q, q = a(k+1)-1:
    The set = the residues r in [a, q-a] hit by the k smallest integers v with
    (v*m mod q) in [a,q-a], BUT we then locally verify M == a/q. We scan a up to a cap
    and ALL coprime m, taking the k smallest safe speeds AND a 'gap-broken' variant where
    we drop the longest initial run (which is what creates competing t) and replace by
    the next safe speeds. Cheap, no full beam search.

We report, per k, the smallest exact M found (largest level a).
"""
import sys
from fractions import Fraction
from math import gcd
from functools import reduce
sys.path.insert(0, ".")
from fast_M import M_exact_fast
from primorial_family import level_a


def setgcd(S):
    return reduce(gcd, S)


def Mval(S):
    return M_exact_fast(sorted(set(S)))[0]


def safe_pool(m, q, a, hi):
    return [v for v in range(1, hi + 1)
            if (lambda r: (r if r <= q - r else q - r) >= a)((v * m) % q)]


def single_killer(k):
    floor = Fraction(1, k + 1)
    best = None
    base = list(range(1, k))
    for X in range(k, 6 * k):
        S = base + [X]
        if setgcd(S) != 1:
            continue
        M = Mval(S)
        if M <= floor:
            continue
        if best is None or M < best[0]:
            best = (M, sorted(S))
    return best


def level_a_variants(k, a, hi_mult=3):
    """Return list of (M, S) for various subset choices from the safe pool at level a."""
    q = a * (k + 1) - 1
    floor = Fraction(1, k + 1)
    target = Fraction(a, q)
    out = []
    for m in range(1, q):
        if gcd(m, q) != 1:
            continue
        pool = safe_pool(m, q, a, hi_mult * q)
        if len(pool) < k:
            continue
        # variant A: k smallest
        SA = pool[:k]
        # variant B: drop the longest initial consecutive run, shift to later safe speeds
        # find initial run length
        run = 1
        while run < len(pool) and pool[run] == pool[run - 1] + 1:
            run += 1
        # keep pool[run-? ...]; simplest: take pool but skip first (run-1) small ones if enough
        if len(pool) >= k + (run - 1):
            SB = pool[(run - 1):(run - 1) + k]
        else:
            SB = None
        for S in [SA, SB]:
            if S is None or len(set(S)) != k:
                continue
            if setgcd(S) != 1:
                continue
            M = Mval(S)
            if M <= floor:
                continue
            out.append((M, sorted(S), m))
    return out, target


if __name__ == "__main__":
    import sympy
    ks = [int(x) for x in sys.argv[1:]] or [30, 31, 42, 43, 60, 61]
    for k in ks:
        floor = Fraction(1, k + 1)
        om = len(sympy.factorint(k - 1))
        # baseline single killer
        sk = single_killer(k)
        bestM, bestS, why = sk[0], sk[1], "single-killer"
        amax = int(k**0.5) + 4
        for a in range(2, amax + 1):
            variants, tgt = level_a_variants(k, a)
            for M, S, m in variants:
                if M < bestM:
                    bestM, bestS, why = M, S, f"level-a a={a} q={a*(k+1)-1} m={m}"
        g = bestM - floor
        print(f"==== k={k} sqrt(k)={k**0.5:.3f} omega(k-1)={om} ====")
        print(f"  BEST M={bestM} (~{float(bestM):.10f}) g*k^2={float(g*k*k):.6f} "
              f"a={float(level_a(bestM,k)):.4f} via {why}")
        print(f"  S={bestS}", flush=True)
