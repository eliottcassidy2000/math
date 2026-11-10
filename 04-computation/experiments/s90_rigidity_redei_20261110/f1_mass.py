#!/usr/bin/env python3
"""HYP-9657: representation mass of the F1 family (one-arc joins of two one-apex tournaments at their apexes).

In the notation psi(2h) = psi(h) + h, psi(2h+1) = psi(h) + floor(h/2) + 1, psi(0) = 0 (one-apex count = 2 psi(h) + 1,
see psi_sbf.py), the F1 values are

        N = (2 psi(h1) + 1)(2 psi(h2) + 1) + 2^s (2 h1 + 1)(2 h2 + 1) - [s = 0],      h1, h2, s >= 0

(THM-4875 (D) / THM-4885 in the coordinates g(u) g(v) + (1 + u)(1 + v) - e(u) e(v); s is the padding).
r(N) = number of triples (h1, h2, s).   This script tabulates r on the odd numbers of a window [2^lo, 2^K):
  * class means of r modulo small primes p (the pattern "N = 0 or -1 mod p is favoured");
  * mean r and miss rate on Rabinowitsch-rough N (N and (N+1)/4 both prime) and on other strata;
  * r at the Heegner numbers.
Usage:  python3 f1_mass.py K lo      (in-sample: 21 18;   pre-registered out-of-sample: 23 21)
"""
import sys

import numpy as np


def main():
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 21
    lo_exp = int(sys.argv[2]) if len(sys.argv) > 2 else 18
    X = 1 << K
    Hmax = X // 2 + 2
    psi = np.zeros(Hmax, dtype=np.int64)
    # vectorised by dyadic blocks
    psi[1] = 1
    k = 1
    while (1 << k) < Hmax:
        lo, hi = 1 << k, min(1 << (k + 1), Hmax)
        idx = np.arange(lo, hi, dtype=np.int64)
        half = idx >> 1
        psi[lo:hi] = psi[half] + np.where(idx % 2 == 0, half, half // 2 + 1)
        k += 1
    g = 2 * psi + 1
    y = 2 * np.arange(Hmax, dtype=np.int64) + 1
    r = np.zeros(X + 1, dtype=np.int32)
    r0 = np.zeros(X + 1, dtype=np.int32)
    for s in range(0, K):
        for h1 in range(0, Hmax):
            y1 = int(y[h1]) << s
            if y1 > X:
                break
            lim = X // y1
            h2max = (lim - 1) // 2
            if h2max < 0:
                break
            N = int(g[h1]) * g[: h2max + 1] + y1 * y[: h2max + 1] - (1 if s == 0 else 0)
            N = N[N <= X]
            np.add.at(r, N, 1)
            if s == 0:
                np.add.at(r0, N, 1)
    sieve = np.ones(X + 2, dtype=bool)
    sieve[:2] = False
    for i in range(2, int((X + 2) ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i:: i] = False
    odd = np.arange(1, X + 1, 2)
    win = odd[odd >= (1 << lo_exp)]
    rw = r[win]
    overall_mean = float(rw.mean())
    overall_miss = float((rw == 0).mean())
    print(f"F1 family, odd N in [2^{lo_exp}, 2^{K}): n = {len(win)}, mean r = {overall_mean:.4f}, "
          f"(s = 0 part {r0[win].mean():.4f}), miss rate = {overall_miss:.5f}, misses = {int((rw == 0).sum())}")
    print("first misses overall:", odd[r[odd] == 0][:40].tolist())
    q1_ok = True
    for p in (3, 5, 7, 11, 13, 17, 19, 23):
        res = win % p
        means = np.array([rw[res == c].mean() for c in range(p)])
        order = np.argsort(-means)
        top2 = set(order[:2].tolist())
        favoured = top2 == {0, p - 1}
        others = np.delete(means, [0, p - 1])
        print(f"p = {p:2d}: class 0: {means[0]:.4f}  class -1: {means[p-1]:.4f}  other classes: mean {others.mean():.4f}, "
              f"max {others.max():.4f}, min {others.min():.4f}   top two classes = {{0, -1}}: {favoured}")
        if p in (3, 5, 7, 11, 17):
            q1_ok &= bool(favoured)

    def stats(mask, name):
        sel = win[mask]
        if len(sel) == 0:
            print(f"{name}: empty")
            return None, None
        m = float(r[sel].mean())
        z = float((r[sel] == 0).mean())
        print(f"{name:62s} n = {len(sel):8d}  mean r = {m:8.4f}  ratio {m/overall_mean:.3f}   miss rate = {z:.5f}  ratio "
              f"{(z/overall_miss if overall_miss > 0 else float('nan')):.3f}")
        return m, z

    isp = sieve[win]
    m4 = (win % 4 == 3)
    qq = (win + 1) // 4
    stats(np.ones(len(win), dtype=bool), "all odd N")
    stats(isp, "N prime")
    mr, zr = stats(m4 & isp & sieve[qq], "N prime and (N+1)/4 prime [Rabinowitsch-rough]")
    stats(m4 & isp & ~sieve[qq], "N prime, N = 3 mod 4, (N+1)/4 composite")
    rough = np.ones(len(win), dtype=bool)
    for p in (3, 5, 7, 11, 13, 17, 19, 23, 29):
        rough &= (win % p != 0) & ((win + 1) % p != 0)
    stats(rough, "N(N+1) free of odd primes <= 29")
    m15 = ((win * (win + 1)) % 15 == 0)
    ms, zs = stats(m15, "15 | N(N+1)")
    print("r at the Heegner numbers:", {d: int(r[d]) for d in (3, 7, 11, 19, 43, 67, 163) if d <= X})
    cand = [n for n in range(7, min(X, 3000), 4) if sieve[n] and sieve[(n + 1) // 4]]
    print("N = 4q - 1 < 3000 with N, q prime, and r(N):", [(n, int(r[n])) for n in cand])
    print("   of which r = 0:", [n for n in cand if r[n] == 0])
    if len(sys.argv) > 3 and sys.argv[3] == "prereg":
        q2_ok = (mr <= 0.80 * overall_mean) and (zr >= 1.8 * overall_miss)
        q3_ok = (ms >= 1.12 * overall_mean) and (zs <= 0.6 * overall_miss)
        print("PRE-REGISTERED (HYP-9657) Q1 [classes 0 and -1 are the top two for p = 3, 5, 7, 11, 17]:", "CONFIRMED" if q1_ok else "REFUTED")
        print("PRE-REGISTERED (HYP-9657) Q2 [rough: mean ratio <= 0.80 and miss-rate ratio >= 1.8]:", "CONFIRMED" if q2_ok else "REFUTED")
        print("PRE-REGISTERED (HYP-9657) Q3 [15 | N(N+1): mean ratio >= 1.12 and miss-rate ratio <= 0.6]:", "CONFIRMED" if q3_ok else "REFUTED")


if __name__ == "__main__":
    main()
