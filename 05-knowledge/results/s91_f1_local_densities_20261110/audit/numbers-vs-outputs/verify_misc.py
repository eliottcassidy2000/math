#!/usr/bin/env python3
"""Auditor numbers-vs-outputs: small independent recomputations.

1. dilogarithm constants of THM-5128 (D1) to 15 digits (mpmath);
2. generic limits p/(p^2 - 1), 8/pi^2, theta_p/2, the candidate limits of HYP-9660 (A);
3. rates max_k rho(T_k)/2 from the full 2m x 2m operator of THM-5128 (B) at selected moduli;
4. F1 on windows of odd N: triples per odd N, misses, and the relative variance var(r)/mean(r)^2 on each sheet
   (the constant 0.75 of PREREGISTRATION.md / window_split.py / README section 5 / THM-5128 Remark 3).
"""
import math
import sys

import numpy as np
from numba import njit
import mpmath as mp


def constants():
    mp.mp.dps = 30
    L2 = mp.polylog(2, mp.mpf(-1) / 2)
    L4 = mp.polylog(2, mp.mpf(-1) / 4)
    k0 = mp.pi ** 2 / 12 + 2 * L2 - L4
    k1 = L4 - L2
    g = 2 * (k0 + k1)
    print("1. constants (mpmath, 30 digits)")
    print("   Li2(-1/2) =", mp.nstr(L2, 15), "  Li2(-1/4) =", mp.nstr(L4, 15))
    print("   kappa_0 =", mp.nstr(k0, 13), "  kappa_1 =", mp.nstr(k1, 13), "  g_bar =", mp.nstr(g, 13), "  pi^2/6 + 2 Li2(-1/2) =", mp.nstr(mp.pi ** 2 / 6 + 2 * L2, 13))
    print("   triple shares w_0 =", mp.nstr(k0 / (k0 + k1), 11), " w_1 =", mp.nstr(k1 / (k0 + k1), 11))
    w0 = k0 / (mp.pi ** 2 / 12)
    print("   word-pair shares w_0 =", mp.nstr(w0, 11), " w_1 =", mp.nstr(1 - w0, 11), "  sigma = 2 w_1 =", mp.nstr(2 * (1 - w0), 9))
    print("   8/pi^2 =", mp.nstr(8 / mp.pi ** 2, 10))
    # check sum_{s>=0} (s+1)(l_s - 2 l_{s+1} + l_{s+2}) = l_0
    l = [-mp.polylog(2, -mp.mpf(2) ** (-j)) for j in range(0, 120)]
    tot = sum((s + 1) * (l[s] - 2 * l[s + 1] + l[s + 2]) for s in range(0, 110))
    print("   sum_{s>=0} (s+1)(l_s - 2 l_(s+1) + l_(s+2)) =", mp.nstr(tot, 15), "  l_0 = pi^2/12 =", mp.nstr(mp.pi ** 2 / 12, 15))
    print("   sum_{s>=1} (l_s - 2 l_(s+1) + l_(s+2)) =", mp.nstr(sum(l[s] - 2 * l[s + 1] + l[s + 2] for s in range(1, 110)), 15), "  l_1 - l_2 =", mp.nstr(l[1] - l[2], 15))
    print("2. generic limits p/(p^2 - 1):", {p: round(p / (p * p - 1), 5) for p in (3, 5, 7, 11, 13, 17, 19, 29, 31, 43, 257, 65537)})
    for p, rate in ((3, 0.8660254), (5, 0.741760), (7, 0.932039), (11, 0.858710)):
        th = math.log(1 / rate) / math.log(2 / rate ** 2)
        rho = rate + 0.005
        th2 = math.log(1 / rho) / math.log(2 / rho ** 2)
        print(f"   p = {p}: theta_p at the spectral rate = {th:.4f} (theta_p/2 = {th / 2:.4f}); at rate + 0.005 (regions_check) = {th2:.4f}")
    sig = float(2 * (1 - w0))
    print("   HYP-9660 (A): rho(p) = L(0)/L(c) for sigma = 2 w_1 and sigma = 3/2")
    for p in (3, 5, 7, 11, 13, 17, 23, 31):
        out = []
        for s in (sig, 1.5):
            a = (s * p - 1) / p ** 3
            out.append((1 + (p - 1) * a) / (1 - a))
        print(f"     p = {p:2d}: {out[0]:.4f}  {out[1]:.4f}")


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


def rates():
    print("3. max_k rho(T_k)/2 from the full 2m x 2m operator (numpy eigvals), k = 1..(m-1)/2")
    for m in (65, 71, 79, 83, 89, 63, 93, 51, 85):
        best = 0.0
        arg = 0
        bestprim = 0.0
        for k in range(1, (m - 1) // 2 + 1):
            r = float(np.max(np.abs(np.linalg.eigvals(twisted_matrix(m, k))))) / 2
            if r > best:
                best, arg = r, k
            if math.gcd(k, m) == 1 and r > bestprim:
                bestprim = r
        print(f"   m = {m:3d}: {best:.6f} at k = +-{arg};  over k prime to m: {bestprim:.6f}")


@njit(cache=False)
def build_q(size):
    q = np.zeros(size, np.int64)
    for h in range(1, size):
        half = h >> 1
        q[h] = q[half] + (((half + 1) >> 1) if (h & 1) else 0)
    return q


@njit(cache=False)
def window(q, A, B, K, r0, r1):
    """r0[i], r1[i] = number of ordered triples (h1, h2, s) on sheet 0 / 1 with value the odd N = A + 1 + 2 i."""
    H = q.shape[0]
    tot0 = 0
    tot1 = 0
    for s in range(K + 1):
        d0 = 1 if s == 0 else 0
        for h1 in range(H):
            y1 = 2 * h1 + 1
            Y = y1 << s
            if Y >= B + d0:
                break
            x1 = y1 - 2 * q[h1]
            # x2 > y2 / 2, so N + d0 > y2 (x1 / 2 + Y): y2 < 2 (B + d0) / (x1 + 2 Y)
            y2max = (2 * (B + d0)) // (x1 + 2 * Y) + 1
            h2max = y2max // 2 + 1
            if h2max > H - 1:
                h2max = H - 1
            for h2 in range(0, h2max + 1):
                y2 = 2 * h2 + 1
                x2 = y2 - 2 * q[h2]
                N = x1 * x2 + Y * y2 - d0
                if N >= A and N < B:
                    i = (N - A - 1) >> 1
                    if s == 0:
                        r0[i] += 1
                        tot0 += 1
                    else:
                        r1[i] += 1
                        tot1 += 1
    return tot0, tot1


def f1_windows(wins):
    print("4. F1 windows: triples per odd N, misses, relative variance var(r)/mean(r)^2")
    for LO, K in wins:
        A, B = 1 << LO, 1 << K
        q = build_q(B // 2 + 2)
        nodd = (B - A) // 2
        r0 = np.zeros(nodd, np.int32)
        r1 = np.zeros(nodd, np.int32)
        t0, t1 = window(q, A, B, K, r0, r1)
        r = r0 + r1
        def rv(a):
            a = a.astype(np.float64)
            return a.var() / a.mean() ** 2
        def fano(a):
            a = a.astype(np.float64)
            return a.var() / a.mean()
        print(f"   [2^{LO}, 2^{K}): triples per odd N: sheet 0 {t0 / nodd:.4f}, sheet 1 {t1 / nodd:.4f}, total {(t0 + t1) / nodd:.4f}; "
              f"misses {int((r == 0).sum())} ({100.0 * (r == 0).sum() / nodd:.3f} %)")
        print(f"       var/mean^2: sheet 0 {rv(r0):.4f}, sheet 1 {rv(r1):.4f}, both sheets {rv(r):.4f};   var/mean: sheet 0 {fano(r0):.3f}, sheet 1 {fano(r1):.3f}, both {fano(r):.3f}")
        # misses by class modulo 3
        Nodd = np.arange(A + 1, B, 2)
        miss = r == 0
        print("       miss rate by N mod 3 (classes 0, 1, 2): " + ", ".join(f"{100.0 * miss[Nodd % 3 == c].mean():.3f} %" for c in range(3)))


if __name__ == "__main__":
    constants()
    rates()
    wins = [(18, 21), (21, 23), (23, 24)]
    if len(sys.argv) > 1 and sys.argv[1] == "big":
        wins.append((24, 27))
    f1_windows(wins)
