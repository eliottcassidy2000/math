#!/usr/bin/env python3
"""collatz_zeckendorf_20260926_interlock_probe.py -- do the Zeckendorf three-colouring and the odd-prime
decomposition interlock with the Collatz step? (session collatz-landing-20260926, opus, 2026-09-26)

Owner's seed: numbers Zeckendorf-decomposed and split into three colours along diagonals; odd numbers arranged
similarly and decomposed into copies of odd primes (composites); parity flows: addition mixes even/odd, odd+odd
returns to even; multiplication keeps any even factor until it is divided out. The Collatz shortcut is exactly
one multiplicative step (3n) followed by one additive step (+1, a parity flip) and the removal of the 2-adic
part. The probe measures, for n <= N, how much the Collatz-relevant quantities are predicted by the colourings.

Colourings (three colours each):
  cZ_low   = index of the smallest Fibonacci term in the Zeckendorf representation, mod 3 (F_2 = 1, F_3 = 2, ...);
  cZ_diag  = (number of terms + index of the largest term) mod 3   ("diagonal" of the (row = #terms, column =
             top index) arrangement);
  cZ_len   = number of Zeckendorf terms mod 3;
  cP_omega = Omega(n) mod 3 for odd n (number of odd prime factors with multiplicity), the additive count of
             "copies of odd primes";
  cP_mod3  = n mod 3 (the colour that the Collatz map actually sees: multiples of 3 have no odd preimage).
Targets: t_par = n mod 2; t_v2 = v_2(3n+1) mod 3 for odd n (how many 2s are factored out); t_desc = [first
descent below n within 20 shortcut steps]; t_stop = total stopping time mod 3 (n <= N reach 1).
We print the mutual information I(colour; target) in bits (uniform three-colour max is log2(3) = 1.585, target
entropies shown) and the colour transition matrix under T for cZ_diag, against the product-of-marginals
baseline. A colouring 'interlocks' if the mutual information is far above the random-colouring level ~ 1/N.
Usage: python3 collatz_zeckendorf_20260926_interlock_probe.py [N=200000]
"""
import math, sys
from collections import Counter


def fib_upto(N):
    F = [1, 2]
    while F[-1] <= N:
        F.append(F[-1] + F[-2])
    return F           # F[i] = F_(i+2) in the usual indexing: 1, 2, 3, 5, 8, ...


def zeck(n, F):
    idx = []
    i = len(F) - 1
    while n > 0:
        while F[i] > n:
            i -= 1
        idx.append(i + 2)
        n -= F[i]
        i -= 2
    return idx


def v2(n):
    c = 0
    while n % 2 == 0:
        n //= 2
        c += 1
    return c


def omega_odd(n, spf):
    c = 0
    while n > 1:
        p = spf[n]
        while n % p == 0:
            n //= p
            c += 1
    return c


def T(n):
    return n // 2 if n % 2 == 0 else (3 * n + 1) // 2


def mutual_info(pairs):
    N = len(pairs)
    cx, cy, cxy = Counter(), Counter(), Counter()
    for x, y in pairs:
        cx[x] += 1; cy[y] += 1; cxy[(x, y)] += 1
    I = 0.0
    for (x, y), c in cxy.items():
        I += (c / N) * math.log2(c * N / (cx[x] * cy[y]))
    Hy = -sum((c / N) * math.log2(c / N) for c in cy.values())
    return I, Hy


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
    F = fib_upto(N)
    spf = list(range(N + 1))
    for i in range(2, int(N ** 0.5) + 1):
        if spf[i] == i:
            for j in range(i * i, N + 1, i):
                if spf[j] == j:
                    spf[j] = i
    stop = {1: 0}
    def stopping(n):
        path = []
        m = n
        while m not in stop:
            path.append(m); m = T(m)
            if m > 10 ** 12:
                return None
        s = stop[m]
        for x in reversed(path):
            s += 1; stop[x] = s
        return stop[n]
    rows = []
    for n in range(2, N + 1):
        z = zeck(n, F)
        cZ_low = z[-1] % 3
        cZ_diag = (len(z) + z[0]) % 3
        cZ_len = len(z) % 3
        par = n % 2
        # first descent within 20 steps
        m = n; desc = 0
        for _ in range(20):
            m = T(m)
            if m < n:
                desc = 1; break
        st = stopping(n)
        tstop = (st % 3) if st is not None else -1
        if n % 2 == 1:
            tv2 = v2(3 * n + 1) % 3
            cP_omega = omega_odd(n, spf) % 3
        else:
            tv2 = -1; cP_omega = -1
        rows.append((n, cZ_low, cZ_diag, cZ_len, cP_omega, n % 3, par, tv2, desc, tstop))
    print("N = %d; mutual information I(colour; target) in bits [target entropy in brackets]" % N)
    cols = {"cZ_low": 1, "cZ_diag": 2, "cZ_len": 3, "cP_omega(odd)": 4, "n mod 3": 5}
    tg = {"parity": 6, "v2(3n+1) mod 3 (odd n)": 7, "descent<=20": 8, "stopping time mod 3": 9}
    print("%-16s" % "" + "".join("%26s" % t for t in tg))
    for cname, ci in cols.items():
        line = "%-16s" % cname
        for tname, ti in tg.items():
            pairs = [(r[ci], r[ti]) for r in rows if r[ci] >= 0 and r[ti] >= 0]
            I, Hy = mutual_info(pairs)
            line += "%18.5f [%5.3f]" % (I, Hy)
        print(line)
    # random-colouring baseline: colour = hash(n) mod 3
    pairs = [((r[0] * 2654435761) % 2 ** 32 % 3, r[8]) for r in rows]
    I, Hy = mutual_info(pairs)
    print("baseline pseudo-random colour vs descent<=20: I = %.6f  (noise level ~ %.1e)" % (I, 1.0 / N))
    # transition matrix of cZ_diag under T
    M = [[0] * 3 for _ in range(3)]
    for r in rows:
        n = r[0]; c = r[2]; m = T(n)
        if m <= N and m >= 2:
            z = zeck(m, F); c2 = (len(z) + z[0]) % 3
            M[c][c2] += 1
    print("cZ_diag transition matrix under T (rows: colour of n; cols: colour of T(n)); row-normalised:")
    for i in range(3):
        s = sum(M[i]) or 1
        print("   colour %d: " % i + "  ".join("%.3f" % (M[i][j] / s) for j in range(3)) + "   (n=%d)" % s)
    marg = [sum(M[i]) for i in range(3)]
    print("   marginal colour frequencies: " + "  ".join("%.3f" % (m / sum(marg)) for m in marg))
    # the multiplicative/additive parity flow: colour of 3n+1 (before halving) for odd n, in terms of cZ_diag(n)
    M2 = [[0] * 3 for _ in range(3)]
    for r in rows:
        n = r[0]
        if n % 2 == 1 and 3 * n + 1 <= N:
            z = zeck(3 * n + 1, F); c2 = (len(z) + z[0]) % 3
            M2[r[2]][c2] += 1
    print("cZ_diag(n) -> cZ_diag(3n+1) for odd n, row-normalised:")
    for i in range(3):
        s = sum(M2[i]) or 1
        print("   colour %d: " % i + "  ".join("%.3f" % (M2[i][j] / s) for j in range(3)) + "   (n=%d)" % s)


if __name__ == '__main__':
    main()
