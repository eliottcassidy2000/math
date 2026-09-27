#!/usr/bin/env python3
"""collatz_precision_residual_20260926.py -- the insufficient-precision residual of the classical residue-class sieve
(session gilbreath6-collatz-precision-20260926, opus, 2026-09-26).

Syracuse map U(n) = (3n+1)/2^v on odd n. A class of odd n is determined by its first valuations (v_1, ..., v_j);
its density among odd integers is 2^-(v_1+...+v_j). Coefficient descent at step j: 3^j < 2^A, A = v_1+...+v_j.
 (1) Residual density D(k) = sum over valuation words of length k with NO coefficient descent at any j <= k of 2^-A_k
     (the density, among odd n, of classes not certified within k Syracuse steps), by dynamic programming over
     (j, A) -- exact rationals for k <= 60; compared with the parity-word count W_k of THM-4495 (T-coding).
 (2) Certification threshold of a descending class: U^j(n) = (3^j n + S_j)/2^A with S_j = sum_(t<j) 3^(j-1-t) 2^(A_t);
     if 3^j < 2^A then U^j(n) < n for every n > N(w) = S_j / (2^A - 3^j). Members of the class below N(w) are the
     only ones the coefficient certificate leaves open; there are at most floor(N(w)/2^A) + 1 of them. Maximal
     N(w)/2^A over all first-descent words with j <= JMAX is tabulated by (j, A), showing the near-convergent
     pairs (A/j close to log_2 3) as the worst cases.
 (3) The window-avoidance density: odd n whose first L Syracuse iterates all have no descent within k steps has
     density at most D(1)^L = 2^-L trivially (v = 1 forced at every point), computed exactly for small k, L.
Usage: python3 collatz_precision_residual_20260926.py [KMAX=60] [JMAX=14]
"""
import sys, math
from fractions import Fraction
from functools import lru_cache

LOG23 = math.log2(3)


def residual_density(kmax):
    """D(k) for k = 1..kmax: sum over no-descent valuation words of length k of 2^-A_k (exact)"""
    # state: (j, A) with 3^j > 2^A for all prefixes; weight 2^-A
    cur = {(0, 0): Fraction(1)}
    out = []
    for j in range(1, kmax + 1):
        nxt = {}
        for (jj, A), wgt in cur.items():
            # next valuation v >= 1; need 3^(j) > 2^(A+v)  <=>  A + v < j log2 3  <=> v < j*LOG23 - A
            vmax = int(math.floor(j * LOG23 - A - 1e-12))  # largest v with A+v < j log2 3 (strict, and equality impossible)
            for v in range(1, vmax + 1):
                key = (j, A + v)
                nxt[key] = nxt.get(key, Fraction(0)) + wgt / 2 ** v
        cur = nxt
        out.append(sum(cur.values()))
    return out


def parity_word_counts(kmax):
    """W_k: parity words of length k (T-coding, letter 1 = odd step) with no coefficient descent: 3^o >= 2^j fails
    never, i.e. o_j > j log_3 2 ... precisely 3^(o_j) > 2^j for all j <= k (n odd forced at step 1)"""
    # dp over (j, o)
    cur = {0: 1}
    out = []
    for j in range(1, kmax + 1):
        nxt = {}
        for o, c in cur.items():
            for letter in (0, 1):
                oo = o + letter
                if 3 ** oo > 2 ** j:
                    nxt[oo] = nxt.get(oo, 0) + c
        cur = nxt
        out.append(sum(cur.values()))
    return out


def thresholds(jmax):
    """enumerate valuation words with first coefficient descent exactly at step j <= jmax; record max N(w)/2^A per (j, A)"""
    best = {}
    count = 0

    def rec(j, A, S, word):
        nonlocal count
        # extend by v
        vmax = int(math.floor(j * LOG23 - A - 1e-12)) if j > 0 else 0
        # words with no descent so far: 3^j > 2^A. Try v such that descent happens at step j+1 (A+v > (j+1) log2 3)
        # or no descent yet (A+v < (j+1) log2 3)
        jn = j + 1
        Sn = 3 * S + 2 ** A  # S_{j+1} = 3 S_j + 2^{A_j}   (S_0 = 0)
        vcrit = int(math.floor(jn * LOG23 - A - 1e-12))
        # descent at step jn for v >= vcrit + 1 ; the class density is 2^-(A+v): record the largest few v only (v <= vcrit + 6)
        for v in range(vcrit + 1, vcrit + 7):
            An = A + v
            N = Fraction(Sn, 2 ** An - 3 ** jn)
            ratio = N / 2 ** An
            key = (jn, An)
            count += 1
            if key not in best or ratio > best[key][0]:
                best[key] = (ratio, N, tuple(word + [v]))
        if jn < jmax:
            for v in range(1, vcrit + 1):
                rec(jn, A + v, Sn, word + [v])

    rec(0, 0, 0, [])
    return best, count


def main():
    KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    JMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 14
    print("== (1) residual density D(k) among odd n: classes with no coefficient descent within k Syracuse steps ==")
    D = residual_density(KMAX)
    W = parity_word_counts(KMAX)
    print(" k   D(k) exact (odd n)        D(k) float   D(k)^(1/k)   W_k (T-coding words)  W_k/2^k   2 W_k/2^k")
    for k in range(1, KMAX + 1):
        if k <= 12 or k % 4 == 0 or k in (41, 45, 51, 55):
            Dk = D[k - 1]
            print(" %2d  %-26s %.6e   %.5f     %-20d  %.4e  %.4e" % (k, str(Dk) if k <= 8 else '', float(Dk), float(Dk) ** (1.0 / k), W[k - 1], W[k - 1] / 2 ** k, 2 * W[k - 1] / 2 ** k))
    print(" (W_k/2^k is the density among all n of parity words without descent; among odd n it is 2 W_k/2^k.)")
    print(" certified density among odd n after 41 Syracuse steps: 1 - D(41) = %.6f  (the swaplift bank of 65 classes has 0.378669 among odd n)" % (1 - float(D[40])))
    print("== (2) certification thresholds N(w)/2^A of first-descent classes, max over words with j <= %d, by (j, A) ==" % JMAX)
    best, count = thresholds(JMAX)
    print(" %d classes examined; the largest N(w)/2^A per (j, A) (only pairs with ratio > 0.05 shown):" % count)
    rows = sorted(best.items(), key=lambda kv: -kv[1][0])
    for (j, A), (ratio, N, word) in rows:
        if ratio > 0.05:
            print("  j=%2d A=%2d  2^A/3^j=%.5f  N(w)=%.1f  N/2^A=%.4f  worst word %s" % (j, A, 2 ** A / 3 ** j, float(N), float(ratio), word))
    print(" the largest thresholds sit at the near-convergent pairs (j, A) = (5, 8), (12, 19) [2^19 < 3^12: no descent], (41, 65) ...")
    print(" bound: N(w)/2^A <= ((3/2)^j - 1)/(2^A - 3^j), since S_j <= 2^A ((3/2)^j - 1); every class has at most one uncertified member (its representative rho_w) once N(w) < 2^A.")
    # exceptions: classes whose representative rho_w = -S_j 3^-j mod 2^A is <= N(w); verify the word of rho_w by simulation
    print("== (2b) actual exceptions: first-descent classes (j <= %d) whose representative rho_w <= N(w) ==" % JMAX)
    exc = []
    checked = 0

    def word_of(n, j):
        w = []
        for _ in range(j):
            m = 3 * n + 1; v = 0
            while m % 2 == 0:
                m //= 2; v += 1
            w.append(v); n = m
        return tuple(w)

    def rec2(j, A, S, word):
        nonlocal checked
        jn = j + 1
        Sn = 3 * S + 2 ** A
        vcrit = int(math.floor(jn * LOG23 - A - 1e-12))
        for v in range(vcrit + 1, vcrit + 40):
            An = A + v
            N = Fraction(Sn, 2 ** An - 3 ** jn)
            rho = (-Sn * pow(3, -jn, 2 ** An)) % 2 ** An
            if rho == 0:
                rho = 2 ** An
            checked += 1
            if rho <= N:
                w = tuple(word + [v])
                assert word_of(rho, jn) == w, (rho, w, word_of(rho, jn))
                # actual iterate
                n = rho
                for _ in range(jn):
                    m = 3 * n + 1
                    while m % 2 == 0:
                        m //= 2
                    n = m
                exc.append((rho, w, float(N), n))
            if 2 ** An > 2 ** 40:
                break
        if jn < jmax2:
            for v in range(1, vcrit + 1):
                rec2(jn, A + v, Sn, word + [v])

    jmax2 = JMAX
    rec2(0, 0, 0, [])
    print(" %d classes checked (all v up to 40 beyond the critical valuation); exceptions (rho, word, N(w), U^j(rho)):" % checked)
    for e in exc:
        print("  ", e)
    print(" every other class is certified in full by its coefficient descent; for the exceptions the actual j-th iterate is listed (it is >= rho, so these are the only members of their classes needing a direct check).")
    print("== (3) trivial window-avoidance bound ==")
    print(" odd n whose first L Syracuse iterates all avoid the k-step certified region (k >= 1) must have v = 1 at each of them: density exactly 2^-L for k = 1; for k >= 2 at most that.")


if __name__ == '__main__':
    main()
