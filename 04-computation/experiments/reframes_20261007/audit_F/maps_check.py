#!/usr/bin/env python3
"""audit_F (B): independent check of the maps in HYP-9244 / coalescence_phase.py.

For each map T(x) = (m_i x + r_i)/d on x = i mod d check, by brute force on residues:
  * gcd(m_i, d) = 1 and d | m_i i + r_i;
  * the branch maps its class i + d Z onto Z_d: for every K <= 4 the map z -> T(i + d z) mod d^K is a bijection of
    Z/d^K (exhaustive);
  * the rank of Gamma = <m_i/m_j> (rank over Q of the matrix of prime-exponent vectors of the ratios m_i/m_0);
  * Lambda = (1/d) sum ln(m_i/d) and its sign.
Standard library only.  Usage: python3 maps_check.py
"""
import math
from fractions import Fraction

MAPS = {  # name: (d, [(m_i, r_i)], claimed rank, claimed sign of Lambda)
    'x+1':      (2, [(1, 0), (1, 1)], 0, -1),
    '3x+1':     (2, [(1, 0), (3, 1)], 1, -1),
    'Z3_124':   (3, [(1, 0), (2, 1), (4, 1)], 1, -1),
    'Z3_125':   (3, [(1, 0), (2, 1), (5, 2)], 2, -1),
    'Z5_12311': (5, [(1, 0), (2, 3), (3, 4), (1, 2), (1, 1)], 2, -1),
    'Z5_12471': (5, [(1, 0), (2, 3), (4, 2), (7, 4), (1, 1)], 2, -1),
    'Z5_12371': (5, [(1, 0), (2, 3), (3, 4), (7, 4), (1, 1)], 3, -1),
    '5x+1':     (2, [(1, 0), (5, 1)], 1, +1),
    'Z3_1416':  (3, [(1, 0), (4, 2), (16, 1)], 1, +1),
    'Z3_157':   (3, [(1, 0), (5, 1), (7, 1)], 2, +1),
}


def factor(n):
    f, p = {}, 2
    while p * p <= n:
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def rank_of(vectors):
    """rank over Q of a list of integer vectors (Gaussian elimination with Fractions)"""
    rows = [[Fraction(x) for x in v] for v in vectors if any(v)]
    if not rows:
        return 0
    ncol = len(rows[0])
    r = 0
    for c in range(ncol):
        piv = next((i for i in range(r, len(rows)) if rows[i][c] != 0), None)
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i][c] != 0:
                f = rows[i][c] / rows[r][c]
                rows[i] = [a - f * b for a, b in zip(rows[i], rows[r])]
        r += 1
    return r


def debt_rank(ms):
    primes = sorted({p for m in ms for p in factor(m)})
    def vec(q):  # exponent vector of the rational q = a/b
        fa, fb = factor(q.numerator), factor(q.denominator)
        return [fa.get(p, 0) - fb.get(p, 0) for p in primes]
    vecs = [vec(Fraction(m, ms[0])) for m in ms]
    return rank_of(vecs) if primes else 0, primes


def branch_onto(d, i, m, r, K):
    seen = set()
    for z in range(d ** K):
        x = i + d * z
        assert (m * x + r) % d == 0
        seen.add(((m * x + r) // d) % d ** K)
    return len(seen) == d ** K


if __name__ == '__main__':
    allok = True
    for name, (d, br, rk_claim, sgn_claim) in MAPS.items():
        ms = [m for m, r in br]
        ok_div = all(math.gcd(m, d) == 1 and (m * i + r) % d == 0 for i, (m, r) in enumerate(br))
        ok_onto = all(branch_onto(d, i, m, r, K) for i, (m, r) in enumerate(br) for K in (1, 2, 3, 4) if d ** K <= 4000)
        rk, primes = debt_rank(ms)
        lam = sum(math.log(m / d) for m in ms) / d
        sgn = -1 if lam < 0 else (1 if lam > 0 else 0)
        ok = ok_div and ok_onto and rk == rk_claim and sgn == sgn_claim
        allok &= ok
        print(f"{name:9s} d={d} m={ms} r={[r for m, r in br]}: d|m_i i+r_i {ok_div}, onto Z_d (mod d^K, K<=4) {ok_onto}, "
              f"rank {rk} (primes {primes}; claimed {rk_claim}), Lambda = {lam:+.4f} (claimed sign {sgn_claim:+d}) -> {'OK' if ok else 'MISMATCH'}")
    print("ALL MAP CHECKS PASSED" if allok else "SOME MAP CHECK FAILED")
