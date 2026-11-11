#!/usr/bin/env python3
"""Audit proof-5128-AB: THM-5128 (B), setting and (B1).

Independent code (pure recursion for psi; matrices rebuilt from the statement).
  1. psi recursion against the canon recursion g(0) = 1, g(2v+1) = g(v), g(2v) = g(v) + v + (v mod 2) (THM-4875 (D)):
     x(h) = 2 psi(h) + 1 = g(2h);  q(h) = h - psi(h) satisfies q(2h) = q(h), q(2h+1) = q(h) + ceil(h/2);
     gcd(x, y) = gcd(2h + 1, q).
  2. the matrix identity sum_{h < 2^n} zeta^(k q(h) + l h) = e_(0,0)^T T_k^n v_l (all k, l; n = 1..NB), brute force
     with direct complex sums (no FFT).
  3. T_0/2 is primitive: (T_0)^t > 0 entrywise as soon as 2^t >= 2m (the proof's claim), and for no smaller t in general.
  4. spectral radii rho(T_k)/2 for all k != 0: by eigenvalues, and a rigorous-in-exact-arithmetic style upper bound
     ||T_k^N||_inf^(1/N) (Gelfand: rho <= ||T^N||^(1/N) for every N), N = 64, 256.
  5. Jordan structure at the dominant modulus (is the top eigenvalue simple?).
  6. the listed rates of the statement.
"""
import math
import sys

import numpy as np


def psi_table(size):
    psi = [0] * size
    for h in range(1, size):
        half = h >> 1
        if h & 1:
            psi[h] = psi[half] + half // 2 + 1
        else:
            psi[h] = psi[half] + half
    return psi


def g_canon(limit):
    g = [0] * limit
    g[0] = 1
    for u in range(1, limit):
        v = u >> 1
        if u & 1:
            g[u] = g[v]
        else:
            g[u] = g[v] + v + (v & 1)
    return g


def T_matrix(m, k):
    inv2 = (m + 1) // 2
    T = np.zeros((2 * m, 2 * m), dtype=np.complex128)
    for eta in range(m):
        for beta in range(2):
            src = 2 * eta + beta
            T[src, 2 * ((2 * eta) % m)] = 1.0
            expo = (k * (((eta + beta) * inv2) % m)) % m
            T[src, 2 * ((2 * eta + 1) % m) + 1] = np.exp(2j * np.pi * expo / m)
    return T


def main():
    ok = True
    NB = 12
    size = 1 << 20
    psi = psi_table(size)
    g = g_canon(2 * size)
    for h in range(size):
        assert 2 * psi[h] + 1 == g[2 * h]
    q = [h - psi[h] for h in range(size)]
    for h in range(size // 2):
        assert q[2 * h] == q[h] and q[2 * h + 1] == q[h] + (h + 1) // 2
    for h in range(1 << 14):
        assert math.gcd(2 * psi[h] + 1, 2 * h + 1) == math.gcd(2 * h + 1, q[h])
    print(f"setting: x(h) = 2 psi(h) + 1 = g(2h) (canon), q recursion, gcd(x, y) = gcd(2h+1, q): OK for h < 2^20")

    # 2. matrix identity, direct sums
    qn = np.array(q[: 1 << NB], dtype=np.int64)
    hn = np.arange(1 << NB, dtype=np.int64)
    worst = 0.0
    for m in (3, 5, 7, 9, 15, 21):
        for k in range(m):
            T = T_matrix(m, k)
            row = np.zeros(2 * m, dtype=np.complex128)
            row[0] = 1.0
            for n in range(1, NB + 1):
                row = row @ T
                for l in range(m):
                    v = np.array([np.exp(2j * np.pi * ((l * eta) % m) / m) for eta in range(m) for _ in range(2)])
                    lhs = np.sum(np.exp(2j * np.pi * ((k * qn[: 1 << n] + l * hn[: 1 << n]) % m) / m))
                    worst = max(worst, abs(lhs - row @ v))
    print(f"matrix identity sum zeta^(k q + l h) = e^T T_k^n v_l, m in (3,5,7,9,15,21), n <= {NB}: max error {worst:.2e}")
    ok &= worst < 1e-7

    # 3. primitivity
    for m in (3, 5, 7, 9, 15, 31, 127):
        T0 = np.abs(T_matrix(m, 0))
        t_need = math.ceil(math.log2(2 * m))
        Pw = np.eye(2 * m)
        first = None
        for t in range(1, t_need + 1):
            Pw = Pw @ T0
            if first is None and np.all(Pw > 0):
                first = t
        ok &= first is not None and first <= t_need
        # doubly stochastic
        ok &= bool(np.all(T0.sum(axis=0) == 2) and np.all(T0.sum(axis=1) == 2))
        print(f"m = {m:3d}: T_0^t > 0 first at t = {first} (claim: for 2^t >= 2m, i.e. t >= {t_need}); row and column sums 2")

    # 4.-6. spectral radii
    listed = {3: 0.8660, 5: 0.7418, 7: 0.9320, 11: 0.8587, 13: 0.8272, 17: 0.8289, 19: 0.8053, 23: 0.7931,
              31: 0.9002, 43: 0.8342, 73: 0.8516, 127: 0.8754}
    primes = [p for p in range(3, 128) if all(p % d for d in range(2, int(p ** 0.5) + 1))]
    others = []
    print(f"{'m':>4} {'max rho/2':>10} {'arg k':>6} {'min rho/2':>10} {'Gelfand N=64':>13} {'N=256':>9} "
          f"{'top simple':>10} {'gap to 2nd modulus':>18} {'listed':>8}")
    for m in [3, 5, 7, 9, 15, 21, 25, 27] + [p for p in primes if p > 7]:
        best, bestk, least = 0.0, None, 9.0
        gel64 = gel256 = 0.0
        simple = True
        gap = None
        for k in range(1, m):
            T = T_matrix(m, k)
            ev = np.linalg.eigvals(T)
            mods = np.sort(np.abs(ev))[::-1]
            rho = mods[0]
            if rho / 2 > best:
                best, bestk = rho / 2, k
                # multiplicity of the top modulus
                ntop = int(np.sum(mods > rho - 1e-7))
                topvals = ev[np.abs(ev) > rho - 1e-7]
                # distinct eigenvalues at the top modulus?
                distinct = len(set(np.round(topvals, 6)))
                simple = (distinct == ntop)
                gap = float(mods[ntop] / rho) if ntop < len(mods) else 0.0
            least = min(least, rho / 2)
            A = np.linalg.matrix_power(T, 64)
            gel64 = max(gel64, float(np.max(np.sum(np.abs(A), axis=1))) ** (1 / 64) / 2)
            A = np.linalg.matrix_power(A, 4)
            gel256 = max(gel256, float(np.max(np.sum(np.abs(A), axis=1))) ** (1 / 256) / 2)
            ok &= rho < 2 - 1e-9
        tag = ""
        if m in listed:
            tag = f"{listed[m]:.4f}"
            if abs(listed[m] - best) > 5.1e-5:
                tag += " MISMATCH"
                ok = False
        elif m in primes and m <= 89:
            others.append((m, best))
        print(f"{m:4d} {best:10.6f} {bestk:6d} {least:10.6f} {gel64:13.6f} {gel256:9.6f} {str(simple):>10} "
              f"{gap:18.6f} {tag:>8}")
    lo = min(v for _, v in others)
    hi = max(v for _, v in others)
    print(f"other primes up to 89 (not listed): {[m for m, _ in others]}")
    print(f"   their max rates lie in [{lo:.4f}, {hi:.4f}]   (statement: between 0.776 and 0.806)")
    for m, v in others:
        if not (0.776 <= v <= 0.806):
            print(f"   OUTSIDE the stated band: m = {m}, rate {v:.6f}")
    # theta values
    for p in (3, 5, 7):
        T_rates = max(np.max(np.abs(np.linalg.eigvals(T_matrix(p, k)))) for k in range(1, p)) / 2
        L = math.log(1 / T_rates)
        theta = L / (math.log(2) + 2 * L)
        print(f"p = {p}: sup of admissible theta_p = {theta:.6f}, theta_p/2 < {theta / 2:.6f}")
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
