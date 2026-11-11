#!/usr/bin/env python3
"""S91 blind proof helper: numerical checks of the lemmas and constants used in the independent proof of the
by-value local densities of F1.

Notation (CONTEXT_S91.md): q(2h) = q(h), q(2h+1) = q(h) + ceil(h/2); y = 2h + 1, x = y - 2q, rho = x/y,
g = gcd(x, y); M = x1 x2 + 2^s y1 y2, N = M - [s = 0]; sheet 0: s = 0, sheet 1: s >= 1; class c' = M mod p
(c' = 0 is the distinguished class on both sheets).

Statements checked (status of each check in brackets):
  1. [FINITE-EXACT, h < 2^NL] Lemma 1: 0 <= q(h) <= h/2 (so y/2 < x <= y), x = y iff h = 0 or a power of two, and
     |rho - (1 - phi/2)| < n/(2h+1) in the exact integer form |rev_n(h) (2h+1) - 2^(n+2) q(h)| < 2 n 2^n.
  2. [VERIFIED float] constants: I_s = int_{1/2}^1 rho^-1 ln((rho+2^s)/(rho/2+2^s)) = l_s - 2 l_(s+1) + l_(s+2),
     l_j = -Li_2(-2^-j); kappa_0 = I_0, kappa_1 = sum_{s>=1} I_s = l_1 - l_2 = int rho^-1 ln(1+rho/2);
     kappa_0 + kappa_1 = (pi^2/6 + 2 Li_2(-1/2))/2; shares; pair normalisation sum (s+1) I_s = pi^2/12.
  3. [VERIFIED float] Lemma 2 (L1), the norm argument: for every cell C (n digits, top J and bottom J prescribed),
     every k != 0 and l:  |sum_{h in C} zeta^(k q + l h)| <= || T_k^(n-2J) ||_(inf->inf)  (max row sum norm),
     and sup_r || T_k^r || / (2 rho')^r is finite for rho' slightly above the spectral rate.
  3b.[VERIFIED float, random regions] Lemma 4: a region {a x + b y <= X}, b >= a > 0, cuts at most 4 top words per
     bottom word among the cells (top j digits, bottom j digits) of a dyadic block.
  4. [FINITE-EXACT enumeration, N <= 2^K] Theorem A/B: T(X) against kappa X ln X, and the exact remainder
        R(c') = T(X; c') - T(X)/p - ([c' = 0] - 1/p) D(X),   D(X) = #triples whose SMALLER state is p-divisible,
     which the proof bounds by O(X^(1 - theta_p/2)); also D(X) - T(X)/p^2 against 2 A_p X.
  5. [NUMERICAL] the first-order constants A_p = sum_h ([p | g(h)] - 1/p^2) w(h) by dyadic blocks, h < 2^NA,
     w^(0) = x^-1 ln((x+y)/(x/2+y)), w^(1) = x^-1 ln(1 + x/(2y)).

Usage: python3 check_blind_proof.py [K] [NA]      (defaults K = 24, NA = 24)
"""
import math
import os
import sys

os.environ.setdefault("NUMBA_NUM_THREADS", "3")

import numpy as np
from numba import njit, prange
import mpmath as mp

PRIMES = np.array([3, 5, 7, 11, 13, 17, 31, 43], dtype=np.int64)
PMAX = int(PRIMES.max())


@njit
def q_table(nbits):
    size = 1 << nbits
    q = np.zeros(size, dtype=np.int32)
    for h in range(1, size):
        half = h >> 1
        if h & 1:
            q[h] = q[half] + ((half + 1) >> 1)
        else:
            q[h] = q[half]
    return q


@njit
def lemma1_check(q, nbits):
    """exact integer checks for h < 2^nbits; returns (ok, worst ratio |rho-(1-phi/2)|(2h+1)/n as float)."""
    ok = True
    worst = 0.0
    for h in range(1, 1 << nbits):
        qq = np.int64(q[h])
        if 2 * qq > h:
            ok = False
        # bit length and bit reversal
        n = 0
        t = h
        rev = np.int64(0)
        while t > 0:
            rev = (rev << 1) | (t & 1)
            t >>= 1
            n += 1
        # phi = rev / 2^n ; check |rev (2h+1) - 4 q 2^n| < 2 n 2^n
        d = rev * (2 * h + 1) - ((4 * qq) << n)
        if d < 0:
            d = -d
        if d >= ((2 * n) << n):
            ok = False
        r = d / float((2 * n) << n)
        if r > worst:
            worst = r
        ispow = (h & (h - 1)) == 0
        if (qq == 0) != ispow:
            ok = False
    return ok, worst


def twisted_matrix(m, k):
    inv2 = pow(2, -1, m)
    zeta = np.exp(2j * np.pi / m)
    T = np.zeros((2 * m, 2 * m), dtype=np.complex128)
    for eta in range(m):
        for beta in range(2):
            i = 2 * eta + beta
            T[i, 2 * ((2 * eta) % m) + 0] += 1.0
            e = (k * ((eta + beta) * inv2 % m)) % m
            T[i, 2 * ((2 * eta + 1) % m) + 1] += zeta ** e
    return T


def norm_argument_check(q, moduli, n, J):
    """max over cells, k != 0, l of |cell sum| against the max-row-sum norm of T_k^(n-2J)."""
    ok = True
    r = n - 2 * J
    h = np.arange(1 << (n - 1), 1 << n, dtype=np.int64)
    qh = q[1 << (n - 1): 1 << n].astype(np.int64)
    rows = []
    for m in moduli:
        zeta = np.exp(2j * np.pi / m)
        worst_ratio = 0.0
        worst_rel = 0.0
        supC = 0.0
        rate = 0.0
        for k in range(1, m):
            T = twisted_matrix(m, k)
            rho = float(np.max(np.abs(np.linalg.eigvals(T)))) / 2.0
            rate = max(rate, rho)
            P = np.eye(2 * m, dtype=np.complex128)
            norms = []
            for _ in range(81):
                norms.append(float(np.max(np.sum(np.abs(P), axis=1))))
                P = P @ T
            nr = norms[r]
            rp = min(0.999, rho * 1.02)
            supC = max(supC, max(norms[i] / (2.0 * rp) ** i for i in range(81)))
            base = zeta ** ((k * qh) % m)
            for l in range(m):
                z = base * zeta ** ((l * h) % m)
                cells = z.reshape(1 << (J - 1), 1 << r, 1 << J).sum(axis=1)
                mx = float(np.max(np.abs(cells)))
                worst_ratio = max(worst_ratio, mx / nr)
                worst_rel = max(worst_rel, mx / float(1 << r))
                if mx > nr * (1 + 1e-9):
                    ok = False
        rows.append((m, rate, worst_ratio, worst_rel, supC))
    return ok, rows


@njit(parallel=True)
def census(q, primes, K, ks, nchunk):
    """Exact enumeration of all triples with N <= 2^K.  Groups: 0: s = 0; 1: s = 1; 2: s = 2; 3: all s >= 1.
    Returns arrays indexed [group, kbucket, prime, ...] (buckets are NOT cumulative)."""
    npr = primes.shape[0]
    nk = ks.shape[0]
    pmax = 0
    for i in range(npr):
        if primes[i] > pmax:
            pmax = primes[i]
    cnt = np.zeros((nchunk, 4, nk, npr, pmax), dtype=np.int64)
    dsm = np.zeros((nchunk, 4, nk, npr), dtype=np.int64)     # smaller state divisible (h1 != h2), ordered pairs
    dfirst = np.zeros((nchunk, 4, nk, npr), dtype=np.int64)  # first state divisible, ordered pairs
    tot = np.zeros((nchunk, 4, nk), dtype=np.int64)
    diag = np.zeros((nchunk, 4, nk), dtype=np.int64)
    thr = np.zeros(nk, dtype=np.int64)
    for i in range(nk):
        thr[i] = np.int64(1) << ks[i]
    for c in prange(nchunk):
        div1 = np.zeros(npr, dtype=np.int64)
        for s in range(0, K + 1):
            two_s = np.int64(1) << s
            xmax = (np.int64(1) << K) + (1 if s == 0 else 0)
            coef = 4 * two_s + 1
            h1 = np.int64(c)
            while True:
                y1 = 2 * h1 + 1
                if coef * y1 * y1 >= 4 * xmax:
                    break
                x1 = y1 - 2 * np.int64(q[h1])
                for i in range(npr):
                    p = primes[i]
                    div1[i] = 1 if (x1 % p == 0 and y1 % p == 0) else 0
                h2 = h1
                while True:
                    y2 = 2 * h2 + 1
                    if coef * y1 * y2 >= 4 * xmax:
                        break
                    x2 = y2 - 2 * np.int64(q[h2])
                    M = x1 * x2 + two_s * y1 * y2
                    if M <= xmax:
                        N = M - (1 if s == 0 else 0)
                        b = 0
                        while thr[b] < N:
                            b += 1
                        w = 2 if h2 > h1 else 1
                        g0 = 0 if s == 0 else 3
                        g1 = s if (s == 1 or s == 2) else -1
                        tot[c, g0, b] += w
                        if g1 >= 0:
                            tot[c, g1, b] += w
                        if h2 == h1:
                            diag[c, g0, b] += 1
                            if g1 >= 0:
                                diag[c, g1, b] += 1
                        for i in range(npr):
                            p = primes[i]
                            cl = M % p
                            d2 = 1 if (x2 % p == 0 and y2 % p == 0) else 0
                            cnt[c, g0, b, i, cl] += w
                            if h2 > h1:
                                dsm[c, g0, b, i] += 2 * div1[i]
                                dfirst[c, g0, b, i] += div1[i] + d2
                            else:
                                dfirst[c, g0, b, i] += div1[i]
                            if g1 >= 0:
                                cnt[c, g1, b, i, cl] += w
                                if h2 > h1:
                                    dsm[c, g1, b, i] += 2 * div1[i]
                                    dfirst[c, g1, b, i] += div1[i] + d2
                                else:
                                    dfirst[c, g1, b, i] += div1[i]
                    h2 += 1
                h1 += nchunk
    return cnt.sum(axis=0), dsm.sum(axis=0), dfirst.sum(axis=0), tot.sum(axis=0), diag.sum(axis=0)


@njit(parallel=True)
def a_blocks(q, primes, NA, nchunk):
    """block sums over B_n (n = 0..NA) of w^(0), w^(1), w_1, w_2 and of the same restricted to p-divisible states.
    q is a table for h < 2^(NA-1); q(h) for the top block is computed from it."""
    npr = primes.shape[0]
    allw = np.zeros((nchunk, NA + 1, 4), dtype=np.float64)
    divw = np.zeros((nchunk, NA + 1, npr, 4), dtype=np.float64)
    size_q = q.shape[0]
    for c in prange(nchunk):
        wv = np.zeros(4, dtype=np.float64)
        h = np.int64(c)
        top = np.int64(1) << NA
        while h < top:
            if h < size_q:
                qq = np.int64(q[h])
            else:
                half = h >> 1
                qq = np.int64(q[half]) + (((half + 1) >> 1) if (h & 1) else 0)
            y = 2 * h + 1
            x = y - 2 * qq
            n = 0
            t = h
            while t > 0:
                t >>= 1
                n += 1
            fx = float(x)
            fy = float(y)
            wv[0] = math.log((fx + fy) / (0.5 * fx + fy)) / fx          # sheet 0
            wv[1] = math.log(1.0 + fx / (2.0 * fy)) / fx                # sheet 1 (all s >= 1)
            wv[2] = math.log((fx + 2.0 * fy) / (0.5 * fx + 2.0 * fy)) / fx   # s = 1
            wv[3] = math.log((fx + 4.0 * fy) / (0.5 * fx + 4.0 * fy)) / fx   # s = 2
            for j in range(4):
                allw[c, n, j] += wv[j]
            for i in range(npr):
                p = primes[i]
                if x % p == 0 and y % p == 0:
                    for j in range(4):
                        divw[c, n, i, j] += wv[j]
            h += nchunk
    return allw.sum(axis=0), divw.sum(axis=0)


def boundary_cells_check(q):
    """Lemma 4: for b >= a > 0 and any X, among the cells (top j digits, bottom j digits) of B_n at most 4 top words
    per bottom word are cut by {a x + b y <= X} (so at most 4 * 2^(n-j) points lie in cut cells)."""
    import random
    random.seed(7)
    worst, worst_pts, tests = 0, 0.0, 0
    for n in (8, 11, 14, 17, 18):
        hb = np.arange(1 << (n - 1), 1 << n, dtype=np.int64)
        yb = (2 * hb + 1).astype(np.float64)
        xb = yb - 2.0 * q[1 << (n - 1): 1 << n].astype(np.float64)
        for j in range(1, n // 2 + 1):
            r = n - 2 * j
            for _ in range(40):
                a = random.uniform(0.01, 1.0)
                b = a * random.choice([1.0, 1.0, 1.5, 2.0, 8.0, 64.0]) * random.uniform(1.0, 1.3)
                X = (2 ** n) * random.uniform(0.9, 2.2) * (a * random.uniform(0.5, 1.0) + b)
                inside = (a * xb + b * yb <= X)
                cnt = inside.reshape(1 << (j - 1), 1 << r, 1 << j).sum(axis=1)      # [top word, bottom word]
                mixed = (cnt > 0) & (cnt < (1 << r))
                worst = max(worst, int(mixed.sum(axis=0).max()))
                worst_pts = max(worst_pts, float(mixed.sum()) * (1 << r) / float(1 << (n - j)))
                tests += 1
    return (worst <= 4 and worst_pts <= 4.0), worst, worst_pts, tests


@njit
def smooth_sums(q, primes, K, ks):
    """Area predictions (no enumeration of pairs).  For X = 2^k and each s:
        Gamma_X(h1) = int_0^1 ( X'/(2 y1 (rho1 (1 - v/2) + 2^s)) - max(h1, 1) )^+ dv,      X' = X + [s = 0],
    the predicted number of h2 > h1 with M <= X'.  Returns sall[group, kbucket] = sum_h1 Gamma and
    sdiv[group, kbucket, prime] = the same over p-divisible h1 (groups as in census)."""
    npr = primes.shape[0]
    nk = ks.shape[0]
    sall = np.zeros((4, nk), dtype=np.float64)
    sdiv = np.zeros((4, nk, npr), dtype=np.float64)
    for bi in range(nk):
        k = ks[bi]
        for s in range(0, k + 1):
            two_s = float(np.int64(1) << s)
            X = float(np.int64(1) << k) + (1.0 if s == 0 else 0.0)
            g0 = 0 if s == 0 else 3
            g1 = s if (s == 1 or s == 2) else -1
            h1 = np.int64(0)
            while True:
                y1 = 2 * h1 + 1
                fy = float(y1)
                if fy * fy * (two_s + 0.25) >= X:
                    break
                x1 = y1 - 2 * np.int64(q[h1])
                fx = float(x1)
                rho = fx / fy
                cut = float(h1) if h1 >= 1 else 1.0
                tau = (X / (2.0 * fy * cut) - two_s) / rho
                rmax = 1.0 if tau > 1.0 else tau
                gam = 0.0
                if rmax > 0.5:
                    gam = (X / fx) * math.log((rho * rmax + two_s) / (0.5 * rho + two_s)) - 2.0 * cut * (rmax - 0.5)
                sall[g0, bi] += gam
                if g1 >= 0:
                    sall[g1, bi] += gam
                for i in range(npr):
                    p = primes[i]
                    if x1 % p == 0 and y1 % p == 0:
                        sdiv[g0, bi, i] += gam
                        if g1 >= 0:
                            sdiv[g1, bi, i] += gam
                h1 += 1
    return sall, sdiv


def main():
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    NA = int(sys.argv[2]) if len(sys.argv) > 2 else 24
    ok = True
    mp.mp.dps = 30

    # ---------------------------------------------------------------- 2. constants
    lam = [-mp.polylog(2, -mp.mpf(2) ** (-j)) for j in range(0, 70)]
    I = [lam[s] - 2 * lam[s + 1] + lam[s + 2] for s in range(0, 60)]
    quad_err = mp.mpf(0)
    for s in range(0, 8):
        val = mp.quad(lambda r: mp.log((r + 2 ** s) / (r / 2 + 2 ** s)) / r, [mp.mpf(1) / 2, 1])
        quad_err = max(quad_err, abs(val - I[s]))
    kappa0 = I[0]
    kappa1 = lam[1] - lam[2]
    k1_quad = mp.quad(lambda r: mp.log(1 + r / 2) / r, [mp.mpf(1) / 2, 1])
    k0_quad = mp.quad(lambda r: mp.log((1 + r) / (1 + r / 2)) / r, [mp.mpf(1) / 2, 1])
    tot_closed = (mp.pi ** 2 / 6 + 2 * mp.polylog(2, -mp.mpf(1) / 2)) / 2
    pair_sum = mp.fsum((s + 1) * I[s] for s in range(60))
    c_ok = (quad_err < mp.mpf(10) ** -20 and abs(kappa1 - k1_quad) < mp.mpf(10) ** -20
            and abs(kappa0 - k0_quad) < mp.mpf(10) ** -20
            and abs(mp.fsum(I[1:60]) - kappa1) < mp.mpf(10) ** -15
            and abs(kappa0 + kappa1 - tot_closed) < mp.mpf(10) ** -25
            and abs(kappa0 - (mp.pi ** 2 / 12 + 2 * mp.polylog(2, -mp.mpf(1) / 2) - mp.polylog(2, -mp.mpf(1) / 4))) < mp.mpf(10) ** -25
            and abs(pair_sum - mp.pi ** 2 / 12) < mp.mpf(10) ** -15)
    ok &= bool(c_ok)
    print("2. constants (mpmath, 30 digits):", bool(c_ok))
    print("   kappa_0 = pi^2/12 + 2 Li2(-1/2) - Li2(-1/4) = %s" % mp.nstr(kappa0, 15))
    print("   kappa_1 = Li2(-1/4) - Li2(-1/2)             = %s" % mp.nstr(kappa1, 15))
    print("   kappa_0 + kappa_1 = (pi^2/6 + 2 Li2(-1/2))/2 = %s   g_bar = %s" % (mp.nstr(kappa0 + kappa1, 15), mp.nstr(2 * (kappa0 + kappa1), 15)))
    print("   shares (triples) w_0 = %s  w_1 = %s" % (mp.nstr(kappa0 / (kappa0 + kappa1), 8), mp.nstr(kappa1 / (kappa0 + kappa1), 8)))
    print("   pair normalisation: sum (s+1) I_s = %s (pi^2/12 = %s), shares %s / %s"
          % (mp.nstr(pair_sum, 12), mp.nstr(mp.pi ** 2 / 12, 12), mp.nstr(kappa0 / pair_sum, 6), mp.nstr(1 - kappa0 / pair_sum, 6)))
    print("   I_s, s = 0..5: " + " ".join(mp.nstr(I[s], 10) for s in range(6)))

    # ---------------------------------------------------------------- tables
    nb = max(K - 1, NA - 1, 22)
    q = q_table(nb)

    # ---------------------------------------------------------------- 1. Lemma 1
    NL = min(nb, 24)
    l1_ok, worst = lemma1_check(q, NL)
    ok &= bool(l1_ok)
    print("1. Lemma 1 for 1 <= h < 2^%d (exact integers): q <= h/2, q = 0 iff h = 2^j, |rho - (1 - phi/2)| < n/(2h+1): %s"
          "   (max of |.|(2h+1)/n = %.4f)" % (NL, bool(l1_ok), worst))

    # ---------------------------------------------------------------- 3. norm argument
    n_cell, J = 16, 3
    na_ok, rows = norm_argument_check(q, [3, 5, 7, 9, 11, 13], n_cell, J)
    ok &= bool(na_ok)
    print("3. cell sums against the max-row-sum norm of T_k^(n-2J), n = %d, J = %d, all cells, k != 0, all l: %s" % (n_cell, J, bool(na_ok)))
    print("   %4s %10s %18s %18s %22s" % ("m", "max rate", "max|cell|/||T^r||", "max|cell|/2^r", "sup_r||T^r||/(2 rho')^r"))
    for (m, rate, wr, wrel, supC) in rows:
        print("   %4d %10.6f %18.4f %18.5f %22.3f" % (m, rate, wr, wrel, supC))
        ok &= bool(rate < 1.0 - 1e-9)

    bc_ok, bc_worst, bc_pts, bc_tests = boundary_cells_check(q)
    ok &= bool(bc_ok)
    print("3b. Lemma 4 boundary cells (%d random regions a x + b y <= X, b >= a, n <= 18, all j <= n/2): at most %d cut top"
          " words per bottom word (proved <= 4), cut points / 2^(n-j) <= %.3f: %s" % (bc_tests, bc_worst, bc_pts, bool(bc_ok)))

    # ---------------------------------------------------------------- 5. first-order constants A_p
    allw, divw = a_blocks(q, PRIMES, NA, 64)
    cum_all = np.cumsum(allw, axis=0)          # [n, j]
    cum_div = np.cumsum(divw, axis=0)          # [n, prime, j]
    wname = ["sheet 0", "sheet 1", "s = 1", "s = 2"]          # columns j of the weights
    wI = [float(kappa0), float(kappa1), float(I[1]), float(I[2])]
    print("5. block sums of w over B_n against I ln 2 (uniform law of rho; n = %d): " % NA
          + ", ".join("%s %.7f / %.7f" % (wname[j], allw[NA, j], wI[j] * math.log(2)) for j in range(4)))
    for j in range(4):
        ok &= bool(abs(allw[NA, j] - wI[j] * math.log(2)) < 1e-4)
        ok &= bool(abs(allw[NA - 4, j] - wI[j] * math.log(2)) < 1e-3)
    A = np.zeros((len(PRIMES), 4))
    nlist = list(range(NA - 10, NA + 1, 2))
    print("   partial sums A_p(n) = sum_{h < 2^n} ([p | g] - 1/p^2) w(h)   [sheet 0 | sheet 1]")
    print("   %4s " % "p" + " ".join("n=%-2d      " % n for n in nlist) + "| same for sheet 1")
    for i, p in enumerate(PRIMES):
        for j in range(4):
            A[i, j] = cum_div[NA, i, j] - cum_all[NA, j] / p ** 2
        s0 = " ".join("%+.6f" % (cum_div[n, i, 0] - cum_all[n, 0] / p ** 2) for n in nlist)
        s1 = " ".join("%+.6f" % (cum_div[n, i, 1] - cum_all[n, 1] / p ** 2) for n in nlist)
        print("   %4d %s | %s" % (p, s0, s1))

    # ---------------------------------------------------------------- 4. census
    ks = np.array([k for k in range(K - 8, K + 1, 2)], dtype=np.int64)
    cnt, dsm, dfirst, tot, diag = census(q, PRIMES, K, ks, 96)
    cnt = np.cumsum(cnt, axis=1)
    dsm = np.cumsum(dsm, axis=1)
    dfirst = np.cumsum(dfirst, axis=1)
    tot = np.cumsum(tot, axis=1)
    diag = np.cumsum(diag, axis=1)
    gname = {0: "sheet 0 (s = 0)", 3: "sheet 1 (s >= 1)", 1: "s = 1", 2: "s = 2"}
    gI = {0: float(kappa0), 3: float(kappa1), 1: float(I[1]), 2: float(I[2])}
    acol = {0: 0, 3: 1, 1: 2, 2: 3}
    print("4. census of all triples with N <= 2^k (exact enumeration)")
    print("   Theorem A:  per k:  T/(I X ln X),  (T - I X ln X)/X (bounded),  slope [T(4X) - 4 T(X)]/(4X ln 4) -> I")
    for g in (0, 3, 1, 2):
        line = []
        for bi, k in enumerate(ks):
            X = float(1 << int(k))
            sl = ""
            if bi > 0:
                slope = (tot[g, bi] - 4 * tot[g, bi - 1]) / (X * math.log(4))
                sl = " slope %.5f" % slope
                if bi == len(ks) - 1 and K >= 22:
                    ok &= bool(abs(slope - gI[g]) < 0.01 * gI[g] + 2e-4)
            line.append("k=%d: %.4f %+.4f%s" % (k, tot[g, bi] / (gI[g] * X * math.log(X)), (tot[g, bi] - gI[g] * X * math.log(X)) / X, sl))
        print("   %-16s I = %.6f | " % (gname[g], gI[g]) + " | ".join(line))
    sall, sdiv = smooth_sums(q, PRIMES, K, ks)
    print("   area formula (no enumeration): (T - 2 sum_h1 Gamma_X(h1) - diag) / X^(3/4)  (proved O(1)):")
    for g in (0, 3, 1, 2):
        vals = []
        for bi, k in enumerate(ks):
            X = float(1 << int(k))
            v = (tot[g, bi] - 2.0 * sall[g, bi] - diag[g, bi]) / X ** 0.75
            vals.append("k=%d: %+.4f" % (k, v))
            ok &= bool(abs(v) < 3.0)
        print("   %-16s " % gname[g] + " | ".join(vals))
    print("   Theorem B: dens = T(c'=0)/T, limit mu = 1/p + (1-1/p)/p^2, F = D/T (D = #triples whose smaller state is")
    print("              p-divisible), Fsm = 2 sum_{h1 div} Gamma_X(h1)/T (area prediction of F), FA = 1/p^2 + 2 A_p X/T,")
    print("              2f-1/p^2 = the session's two-sided Delta (f = share of triples with FIRST state divisible),")
    print("              rem = max_c |R(c')|/T, R = T(c') - T/p - ([c'=0]-1/p) D     (proved O(X^(-theta_p/2)/ln X)),")
    print("              rem2 = the same with D replaced by 2 Dfirst - T/p^2,  ripple = (max - min of T(c'), c' != 0)/T,")
    print("              dD = (D - 2 sum_{h1 div} Gamma_X(h1)) / X^(3/4)             (proved O(1), uniformly in p)")
    for g in (0, 3, 1, 2):
        print("   --- %s" % gname[g])
        print("   %4s %3s %9s %9s %9s %9s %9s %9s %9s %9s %9s %8s"
              % ("p", "k", "dens", "mu", "F", "Fsm", "FA", "2f-1/p^2", "rem", "rem2", "ripple", "dD"))
        for i, p in enumerate(PRIMES):
            p = int(p)
            for bi, k in enumerate(ks):
                X = float(1 << int(k))
                T = int(tot[g, bi])
                D = int(dsm[g, bi, i])
                Df = int(dfirst[g, bi, i])
                cl = [int(v) for v in cnt[g, bi, i, :p]]
                ok &= (sum(cl) == T)
                Rp = [p * cl[c] - T - ((p - 1) if c == 0 else -1) * D for c in range(p)]     # p R(c'), exact
                ok &= (sum(Rp) == 0)
                rem = max(abs(v) for v in Rp) / p / T
                D2 = 2 * Df - T / p ** 2
                rem2 = max(abs(cl[c] - T / p - ((1 if c == 0 else 0) - 1 / p) * D2) for c in range(p)) / T
                rip = (max(cl[1:]) - min(cl[1:])) / T
                mu = 1 / p + (1 - 1 / p) / p ** 2
                FA = 1 / p ** 2 + 2 * A[i, acol[g]] * X / T
                Fsm = 2.0 * sdiv[g, bi, i] / T
                dD = (D - 2.0 * sdiv[g, bi, i]) / X ** 0.75
                ok &= bool(abs(dD) < 3.0)
                if bi in (0, len(ks) // 2, len(ks) - 1):
                    print("   %4d %3d %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.2e %9.2e %9.2e %+8.4f"
                          % (p, k, cl[0] / T, mu, D / T, Fsm, FA, D2 / T, rem, rem2, rip, dD))
                if bi == len(ks) - 1 and K >= 22:
                    ok &= bool(rem < 0.03)
                    ok &= bool(abs(D / T - FA) < 0.03)
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
