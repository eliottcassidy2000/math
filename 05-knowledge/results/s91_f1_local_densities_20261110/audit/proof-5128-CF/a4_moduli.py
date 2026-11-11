#!/usr/bin/env python3
"""Auditor proof-5128-CF, check A4: composite and prime-power moduli (THM-5128 (C) squarefree case, (D3) product law
and the "singular series").

usage: python3 a4_moduli.py K     (default 24)

(1) m = 15 and m = 21 (squarefree): the product law of (C),
      C(c) = sum_{d | m} D^(d) [c = c_sigma mod d] d/m + R(c),   D^(d) = #{ordered triples, smaller state with gcd(g, m) = d}
    and the limit prod_{p | m} mu_c^(p).
(2) m = 9, 25, 27 (prime powers, NOT covered by the theorem): is the density flat on the residue class c_sigma mod p,
    as a singular series with the single factor 1 - 1/p^2 + [p | N + [sigma = 0]]/p would say?
    First-order law at p^2 (same mechanism as (C)):   with v = v_p(g(smaller state)),
      C(c) ~ D_0/p^2 + D_1 [p | c']/p + D_2 [p^2 | c'],  c' = c - c_sigma,  D_v = #{triples with v_p(g_small) = v (v = 2: >= 2)}
    Limit p-adic factor (if the states are uniform in Z_p^2):  (1 + 1/p)(1 - p^(-v - 1)),  v = v_p(N - c_sigma).
"""
import sys
import time
import numpy as np
from a_lib import psi_table, census_mod

K = int(sys.argv[1]) if len(sys.argv) > 1 else 24
KMIN = 16
t0 = time.time()
HN = (1 << (K - 1)) + 2
psi = psi_table(HN)
x64 = 2 * psi + 1
x = x64.astype(np.int32)
h = np.arange(HN, dtype=np.int64)
g = np.gcd(x64, 2 * h + 1)


def squarefree(m, primes):
    divs = [d for d in range(1, m + 1) if m % d == 0]
    gm = np.gcd(g, m)
    cls = np.zeros(HN, dtype=np.uint8)
    for i, d in enumerate(divs):
        cls[gm == d] = i
    hist, diag = census_mod(K, KMIN, x, m, cls, len(divs))
    hist = np.cumsum(hist, axis=1)
    diag = np.cumsum(diag, axis=1)
    print("\nm = %d (squarefree), divisors %s" % (m, divs))
    for sheet in (0, 1):
        csig = (m - 1) if sheet == 0 else 0
        print("  sheet %d:  k |          T | m max|R|/T | forced-class violations | shares F^(d), d = %s | limits %s"
              % (sheet, divs, ["%.4f" % np.prod([(1 / p ** 2) if d % p == 0 else (1 - 1 / p ** 2) for p in primes])
                               for d in divs]))
        for k in range(KMIN, K + 1, 2):
            hk = hist[sheet, k - KMIN]
            dk = diag[sheet, k - KMIN]
            T = 2 * hk.sum() + dk.sum()
            C = 2 * hk.sum(axis=(0, 1)) + dk.sum(axis=0)
            pred = np.zeros(m)
            viol = 0
            Fd = []
            for i, d in enumerate(divs):
                Dd = 2 * hk[i].sum()
                Fd.append(Dd / T)
                comp = np.array([(c - csig) % d == 0 for c in range(m)])
                pred += comp * Dd * d / m
                viol += int(hk[i][:, ~comp].sum())
            R = C - pred
            print("            %2d | %10d |  %.5f   | %d | %s" % (k, T, m * np.abs(R).max() / T, viol,
                                                              " ".join("%.4f" % f for f in Fd)))
        # class shares at the largest X against the limit product
        hk = hist[sheet, K - KMIN]
        dk = diag[sheet, K - KMIN]
        T = 2 * hk.sum() + dk.sum()
        C = 2 * hk.sum(axis=(0, 1)) + dk.sum(axis=0)
        lim = np.ones(m)
        for p in primes:
            for c in range(m):
                lim[c] *= (1 - 1 / p ** 2 + (1 / p if (c - csig) % p == 0 else 0))
        print("    m * share(c) at X = 2^%d, grouped by the pattern of c - c_sigma:" % K)
        for d in divs:
            sel = [c for c in range(m) if np.gcd((c - csig) % m, m) == d or (d == m and (c - csig) % m == 0)]
            if not sel:
                continue
            print("      gcd(c - c_sigma, m) = %2d: measured mean %.4f (min %.4f, max %.4f), limit %.4f"
                  % (d, m * C[sel].mean() / T, m * C[sel].min() / T, m * C[sel].max() / T, lim[sel[0]]))


def prime_power(p, e):
    m = p ** e
    v = np.zeros(HN, dtype=np.uint8)
    for j in range(1, e + 1):
        v[g % (p ** j) == 0] = j
    hist, diag = census_mod(K, KMIN, x, m, v, e + 1)
    hist = np.cumsum(hist, axis=1)
    diag = np.cumsum(diag, axis=1)
    print("\nm = %d = %d^%d (prime power: outside the theorem)" % (m, p, e))
    first = [int(np.argmax(v >= j)) for j in range(1, e + 1)]
    print("  first states with p^j | g, j = 1..%d: %s" % (e, first))
    for sheet in (0, 1):
        csig = (m - 1) if sheet == 0 else 0
        print("  sheet %d: m * share(c) by v = v_p(c - c_sigma)   [theorem's single factor at p | c - c_sigma: "
              "1 - 1/p^2 + 1/p = %.4f for every v >= 1]" % (sheet, 1 - 1 / p ** 2 + 1 / p))
        hdr = "      k |" + "".join(" v=%d meas / 1st-order |" % j for j in range(e + 1)) + " mean over v>=1 | F_v (v=1..)"
        print(hdr)
        for k in range(KMIN, K + 1, 2):
            hk = hist[sheet, k - KMIN]
            dk = diag[sheet, k - KMIN]
            T = 2 * hk.sum() + dk.sum()
            C = 2 * hk.sum(axis=(0, 1)) + dk.sum(axis=0)
            Dv = [2 * hk[j].sum() for j in range(e + 1)]
            Dv[0] += dk.sum()          # diagonal triples put with v = 0 (negligible)
            out = []
            tot1 = 0.0
            n1 = 0
            for j in range(e + 1):
                sel = [c for c in range(m) if ((c - csig) % m == 0 and j == e) or
                       ((c - csig) % m != 0 and (c - csig) % (p ** j) == 0 and (c - csig) % (p ** (j + 1)) != 0)]
                meas = m * C[sel].mean() / T
                # first order: a triple with v_small = i puts its mass uniformly on the classes with p^i | c'
                fo = sum(Dv[i] * p ** i / m for i in range(0, j + 1)) * m / T
                out.append(" %.4f / %.4f       |" % (meas, fo))
                if j >= 1:
                    tot1 += C[sel].sum()
                    n1 += len(sel)
            print("     %2d |" % k + "".join(out) + "   %.4f       | " % (m * tot1 / n1 / T)
                  + " ".join("%.5f" % (Dv[j] / T) for j in range(1, e + 1)))
        print("      limit (uniform states in Z_p^2): " + "  ".join(
            "v=%d: %.4f" % (j, (1 + 1 / p) * (1 - p ** (-j - 1.0))) for j in range(e))
              + "  v>=%d: %.4f" % (e, (1 + 1 / p) - p ** (-e - 1.0) * (1 + 1 / p) * (1 - 1 / p) / (1 - 1 / p ** 2)))


squarefree(15, [3, 5])
squarefree(21, [3, 7])
prime_power(3, 2)
prime_power(5, 2)
prime_power(3, 3)
print("\ntotal time %.1f s" % (time.time() - t0))
