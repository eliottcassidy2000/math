#!/usr/bin/env python3
"""collatz_cycle_shadows_20260926.py -- growth families as 2-adic shadows of negative (rational) cycles
(session fibonacci-two-copies-20260926, opus, 2026-09-26).

Syracuse map U(n) = (3n+1)/2^v on odd integers, both signs; for a valuation word w = (v_1..v_p), A = sum v,
the affine map is n -> (3^p n + S_w)/2^A, whose fixed point x_w = S_w/(2^A - 3^p) is the unique 2-adic point
with parity word w^infinity (Terras/Lagarias); it is a negative rational when 3^p > 2^A (a no-descent word:
growth factor 3^p/2^A per period) and positive when 3^p < 2^A.
 (1) Negative odd cycles of U with |n| <= 10^6: expected exactly -1 (word (1)), -5 (word (1,2)),
     -17 (word (1,1,1,2,1,1,4)); their words, factors, and the shadow rule n = x mod 2^(mA+1) => the orbit follows
     w^m and U^(mp)(n) = (3^p/2^A)^m (n - x) + x exactly (a Syracuse word with total valuation A is a residue
     class modulo 2^(A+1): the last valuation must be exact); checked on random members.
 (2) Rational cycles: for every primitive no-descent word of length p <= 8 (up to rotation), the fixed point
     x_w, whether it is an integer, and a numerical check of the shadow rule; count of integer ones.
 (3) The lane's families: n_(k+1) = 8 n_k + 35 (27, 251, 2043, ...) converges 2-adically to -5;
     n_(k+1) = 4 n_k - 17 ((4^k+17)/3: 7, 11, 27, 91, 347, ...) converges to 17/3, a positive rational whose orbit
     is 17/3 -> 9 -> 7 -> 11 -> 17 -> 13 -> 5 -> 1: the family's large members shadow a descending point.
 (4) At precision K bits, the density covered by the three integer-cycle shadow rules versus the whole
     no-descent residual (from collatz_precision_residual_20260926.py's D(k)): the residual is dominated by
     aperiodic words.
Usage: python3 collatz_cycle_shadows_20260926.py
"""
import math, random
from fractions import Fraction

LOG23 = math.log2(3)


def U(n):
    m = 3 * n + 1
    v = 0
    while m % 2 == 0:
        m //= 2; v += 1
    return m, v


def word_of(n, p):
    w = []
    for _ in range(p):
        n, v = U(n); w.append(v)
    return tuple(w), n


def affine(w):
    p = len(w); A = sum(w)
    S = 0
    At = 0
    for t in range(p):
        S = S + 3 ** (p - 1 - t) * 2 ** At
        At += w[t]
    return p, A, S


def negative_cycles(limit):
    seen = set(); cycles = []
    for n in range(-1, -limit, -2):
        if n in seen:
            continue
        path = []; m = n; pathset = set()
        while m not in pathset and abs(m) < 50 * limit:
            pathset.add(m); path.append(m); m, _ = U(m)
        seen |= pathset
        if m in pathset:
            cyc = path[path.index(m):]
            mn = min(cyc)
            if mn not in [min(c) for c in cycles]:
                cycles.append(cyc)
    return cycles


def twoadic(x, K):
    """x = a/b with b odd: the residue of x mod 2^K"""
    return (x.numerator * pow(x.denominator, -1, 2 ** K)) % 2 ** K


def main():
    print("== (1) negative odd cycles of U with |n| <= 10^6 ==")
    cycles = negative_cycles(10 ** 6)
    rng = random.Random(1)
    for cyc in cycles:
        x = cyc[0]
        w, back = word_of(x, len(cyc))
        assert back == x
        p, A, S = affine(w)
        xw = Fraction(S, 2 ** A - 3 ** p)
        assert xw == x, (xw, x)
        print(" cycle through %d: length %d, word %s, A = %d, 3^p/2^A = %.5f, fixed point S/(2^A-3^p) = %s" % (x, len(cyc), w, A, 3 ** p / 2 ** A, xw))
        # shadow rule
        for m in (1, 3, 10):
            K = m * A + 1
            r = twoadic(Fraction(x), K)
            ok = True
            for _ in range(20):
                n = r + 2 ** K * rng.randrange(1, 1000)
                ww, end = word_of(n, m * p)
                ok &= (ww == w * m) and (Fraction(end) == Fraction(3 ** p, 2 ** A) ** m * (n - x) + x)
            print("   n = %d mod 2^%d follows w^%d and U^(%d)(n) = (3^p/2^A)^%d (n - x) + x: %s" % (r, K, m, m * p, m, ok))
    print("== (2) rational cycles of primitive no-descent words, p <= 8 ==")
    integer_ones = []
    count = 0
    for p in range(1, 9):
        # enumerate words with 3^p > 2^A and no descent at any prefix: A_j < j log2 3 for all j
        def rec(prefix, At):
            nonlocal count
            j = len(prefix)
            if j == p:
                w = tuple(prefix)
                # primitive: not a power of a shorter word; canonical under rotation: minimal rotation
                rots = [w[i:] + w[:i] for i in range(p)]
                if w != min(rots):
                    return
                if any(w == rots[i] for i in range(1, p)):
                    return  # periodic (a power)
                pp, A, S = affine(w)
                x = Fraction(S, 2 ** A - 3 ** pp)
                count += 1
                if x.denominator == 1:
                    integer_ones.append((w, x))
                # numerical shadow check with m = 2
                K = 2 * A + 1
                r = twoadic(x, K)
                n = r + 2 ** K * 7
                ww, end = word_of(n, 2 * pp)
                assert ww == w * 2 and Fraction(end) == Fraction(3 ** pp, 2 ** A) ** 2 * (n - x) + x, (w, x)
                return
            vmax = int(math.floor((j + 1) * LOG23 - At - 1e-12))
            for v in range(1, vmax + 1):
                rec(prefix + [v], At + v)
        rec([], 0)
    print(" primitive no-descent words up to rotation with p <= 8: %d; all shadow rules verified (m = 2); integer fixed points: %s" % (count, integer_ones))
    print("== (3) the lane's families ==")
    fam = [27]
    for _ in range(6):
        fam.append(8 * fam[-1] + 35)
    print(" 8n+35 family:", fam, "; fixed point of n -> 8n+35 is -5; n_k = -5 mod 2^(3k):", [(n + 5) % 2 ** (3 * (k + 1)) == 0 for k, n in enumerate(fam)])
    fam2 = [(4 ** k + 17) // 3 for k in range(1, 12)]
    print(" (4^k+17)/3 family:", fam2[:8], "; fixed point of n -> 4n-17 is 17/3; 2-adic residues of 17/3 mod 2^(2k) match:", [(n - twoadic(Fraction(17, 3), 2 * k)) % 2 ** (2 * k) == 0 for k, n in zip(range(1, 12), fam2)])
    x = Fraction(17, 3); orb = [x]
    for _ in range(8):
        m = 3 * x + 1
        while twoadic(m, 1) == 0:
            m /= 2
        x = m; orb.append(x)
    print(" orbit of 17/3 under U (2-adic parity):", [str(o) for o in orb])
    sig = []
    for n in fam2[2:9]:
        m = n; j = 0
        while True:
            m, _ = U(m); j += 1
            if m < n:
                break
        sig.append(j)
    print(" first-descent times of (4^k+17)/3 for k = 3..9:", sig)
    print("== (4) density covered by the three integer-cycle rules at precision K ==")
    for K in (12, 24, 40, 64):
        d = sum(2.0 ** (-K) for _ in range(3))  # each rule is a single residue class mod 2^K
        print(" K = %d bits: three integer rules cover density 3 * 2^-%d = %.2e of all integers; the no-descent residual at ~K/1.58 Syracuse steps is of order 10^-2..10^-3 (D(k) table): aperiodic words dominate" % (K, K, d))


if __name__ == '__main__':
    main()
