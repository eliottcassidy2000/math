#!/usr/bin/env python3
"""Auditor proof-5134-FH.  THM-5134 (H)(iii): own by-value census of the base-b join family.

N_b(h1, h2, s) = X1 X2 + b (2^s y1 y2 - [s = 0]) / 2,   X_i = 1 + b psi(h_i),  y_i = 2 h_i + 1,
over ORDERED triples (h1, h2, s) with N_b <= X = 2^K; sheet 0: s = 0, sheet 1: s >= 1.

For an odd modulus m prime to b:
  type of a triple with h1 != h2  = gcd(X(h), y(h), m) of its state of smaller word h = min(h1, h2);
  D^(d) = number of ordered triples (h1 != h2) of type d;
  e(N) = N on sheet 1, 2 N + b on sheet 0   (so the distinguished class is N = 0, resp. N = -b/2);
  exact structure: every triple of type d has d | e(N)                          (violations counted);
  first-order law: #{N = c} ~ sum_{d | gcd(m, e(c))} (d/m) D^(d);
  remainder      : m * max_c | #{N = c} - sum_{d | gcd(m, e(c))} (d/m) D^(d) | / #T.
The diagonal triples h1 = h2 are counted in #T and in the classes, not in any D^(d) (as in THM-5128 (C2)); they are
reported separately.

The enumeration uses only lower bounds that are valid for every larger word:  psi(h) >= (h + 1) // 2, so for fixed
h1 the scan over h2 stops when X1 (1 + b ((h2 + 1) // 2)) + seam(h1, h2) > X, and the scan over h1 stops when the
same bound at h2 = h1 exceeds X.  No state is skipped because its own diagonal value is large.

usage: h3_base_b_census.py K sheet b1,b2,...        (sheet = 0 or 1)
"""
import sys
import time
from math import gcd
import numpy as np
from numba import njit


def q_table(nbits):
    q = np.zeros(1 << nbits, dtype=np.int64)
    for k in range(nbits):
        lo, hi = 1 << k, 1 << (k + 1)
        idx = np.arange(lo, hi, dtype=np.int64)
        half = idx >> 1
        q[lo:hi] = q[half] + np.where((idx & 1) == 1, (half + 1) >> 1, 0)
    return q


@njit(cache=True)
def gcd_nb(a, b):
    while b:
        a, b = b, a % b
    return a


@njit(cache=True)
def census(psi, b, X, sheet, mods, cnt, dtype_cnt, viol, diag, hmax_used):
    """cnt[j, c] : ordered triples with N = c mod mods[j]
       dtype_cnt[j, d] : ordered triples (h1 != h2) of type d
       viol[j] : triples of type d with d not dividing e(N)
       diag[j, 0] : diagonal triples;  diag[j, 1]: diagonal triples of type > 1 ; diag[j, 2]: of them with type | e(N)
       returns (#T ordered, #unordered pairs, bad congruence count N != 1 mod b)"""
    H = psi.shape[0]
    total = 0
    bad = 0
    nm = mods.shape[0]
    smax = 0 if sheet == 0 else 62
    s0 = 0 if sheet == 0 else 1
    for s in range(s0, smax + 1):
        if sheet == 1 and (b << (s - 1)) > X:
            break
        any_s = False
        for h1 in range(H):
            y1 = 2 * h1 + 1
            lb1 = 1 + b * ((h1 + 1) // 2)
            if sheet == 0:
                seam_lb = b * (y1 * y1 - 1) // 2
            else:
                seam_lb = (b << (s - 1)) * y1 * y1
            if lb1 * lb1 + seam_lb > X:
                break
            X1 = 1 + b * psi[h1]
            for h2 in range(h1, H):
                y2 = 2 * h2 + 1
                if sheet == 0:
                    seam = b * (y1 * y2 - 1) // 2
                else:
                    seam = (b << (s - 1)) * y1 * y2
                if X1 * (1 + b * ((h2 + 1) // 2)) + seam > X:
                    break
                X2 = 1 + b * psi[h2]
                N = X1 * X2 + seam
                if N > X:
                    continue
                any_s = True
                if h2 > hmax_used[0]:
                    hmax_used[0] = h2
                w = 1 if h2 == h1 else 2
                total += w
                if N % b != 1:
                    bad += 1
                e = N if sheet == 1 else 2 * N + b
                for j in range(nm):
                    m = mods[j]
                    cnt[j, N % m] += w
                    d = gcd_nb(gcd_nb(X1 % m, y1 % m), m)
                    if h2 != h1:
                        dtype_cnt[j, d] += w
                        if e % d != 0:
                            viol[j] += w
                    else:
                        diag[j, 0] += 1
                        if d > 1:
                            diag[j, 1] += 1
                            if e % d == 0:
                                diag[j, 2] += 1
        if not any_s and s > s0 + 2:
            pass
    return total, bad


def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def main():
    K = int(sys.argv[1])
    sheet = int(sys.argv[2])
    bases = [int(t) for t in sys.argv[3].split(",")]
    X = 1 << K
    nbits = min(K, 24)
    q = q_table(nbits)
    h = np.arange(1 << nbits, dtype=np.int64)
    psi = h - q
    assert bool(np.all(2 * psi[1:] > h[1:]))            # psi(h) > h/2 for h >= 1 (used in the lower bounds)
    assert bool(np.all(psi >= (h + 1) // 2))
    print(f"own census, auditor proof-5134-FH: N_b <= 2^{K}, sheet {'s = 0' if sheet == 0 else 's >= 1'}")
    print("class e(N) = N (sheet 1), 2N + b (sheet 0); remainder = m max_c |#{N = c} - sum_{d | gcd(m, e(c))} (d/m) D^(d)| / #T")
    print("columns per modulus: factor of the distinguished class [min, max of the other classes]; F m^2 = (share of triples "
          "whose smaller state has gcd = m) * m^2; rem; noise = sqrt(m/#T)")
    allmods = [3, 5, 7, 9, 11, 13, 25, 35, 49, 55, 77, 91]
    worst_prime = {}
    for b in bases:
        mods = np.array([m for m in allmods if gcd(m, 2 * b) == 1], dtype=np.int64)
        mmax = int(mods.max())
        cnt = np.zeros((len(mods), mmax), dtype=np.int64)
        dt = np.zeros((len(mods), mmax + 1), dtype=np.int64)
        viol = np.zeros(len(mods), dtype=np.int64)
        diag = np.zeros((len(mods), 3), dtype=np.int64)
        hmax = np.zeros(1, dtype=np.int64)
        t0 = time.time()
        total, bad = census(psi, b, X, sheet, mods, cnt, dt, viol, diag, hmax)
        assert int(hmax[0]) < psi.shape[0] - 1, "psi table too short"
        print(f"\nb = {b}: #T = {total} ordered triples ({diag[0, 0]} diagonal); N_b = 1 (mod b) for all: {bad == 0}; "
              f"largest word used {int(hmax[0])}; {time.time() - t0:.1f} s")
        for j, m in enumerate(mods.tolist()):
            m = int(m)
            assert int(cnt[j, :m].sum()) == total
            # distinguished class
            if sheet == 1:
                c0 = 0
            else:
                c0 = (-b * pow(2, -1, m)) % m
            # prediction
            rem = 0.0
            fac = []
            for c in range(m):
                e = c if sheet == 1 else (2 * c + b)
                g = gcd(m, e % m) if e % m else m
                pred = sum(d * int(dt[j, d]) for d in divisors(g)) / m
                rem = max(rem, abs(int(cnt[j, c]) - pred))
                fac.append(m * int(cnt[j, c]) / total)
            rem = m * rem / total
            others = [fac[c] for c in range(m) if c != c0]
            Fm = int(dt[j, m]) / total
            tag = "prime" if all(m % t for t in range(2, m)) else "composite"
            print(f"  m = {m:3d} ({tag}): class {c0:2d}: {fac[c0]:.4f} [{min(others):.4f}, {max(others):.4f}]  "
                  f"F m^2 = {Fm * m * m:.3f}  rem = {rem:.4f}  noise = {(m / total) ** 0.5:.4f}  "
                  f"type violations = {int(viol[j])}  types " +
                  ", ".join(f"{d}: {int(dt[j, d])}" for d in divisors(m)) +
                  f"  diagonal of type > 1: {int(diag[j, 1])} (all in a class with type | e(N): {int(diag[j, 1]) == int(diag[j, 2])})")
            if tag == "prime" and m in (3, 5, 7, 11, 13):
                worst_prime[(b, m)] = rem
    print("\nlargest remainder over the primes 3..13 not dividing b, per base:")
    for b in bases:
        items = [(m, r) for (bb, m), r in worst_prime.items() if bb == b]
        mm, rr = max(items, key=lambda t: t[1])
        print(f"  b = {b}: {rr:.4f} at p = {mm}")
    print(f"overall: {max(worst_prime.values()):.4f}")


if __name__ == "__main__":
    main()
