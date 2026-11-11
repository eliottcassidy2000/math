#!/usr/bin/env python3
"""S91 (THM-5128 (E), remark): the mod-4 law of the one-apex family and the symbol (1 + z)^2 (1 + z + z^2) over F_2.

q(h) = sum_i h_i (h_(i+1) + h_(i+2)) (mod 2) in the binary digits h_0 (lowest), h_1, ..., h_(n-1) of h < 2^n, and
psi(h) = h - q(h), so (H - 1)/2 = psi(h) = h_0 + q(h) (mod 2) for the one-apex count H = 2 psi(h) + 1 (the chain case
of THM-5127 (C)).  The polar form of q is the n x n band matrix B_n over F_2 with ones at distances 1 and 2 from the
diagonal: a Toeplitz matrix with symbol z + z^2 + z^-1 + z^-2 = z^-2 (1 + z)^2 (1 + z + z^2).

Computed exactly:
  * rank of B_n over F_2 for n <= NMAX and its period in n;
  * S_n = sum_{h < 2^n} (-1)^q(h) and D_n = sum_{h < 2^n} (-1)^psi(h) = #{H = 1 mod 4} - #{H = 3 mod 4} by the digit
    recursion (exact integers), against 2^(n - rank/2);
  * brute force for n <= 20.
"""
import sys

import numpy as np

NMAX = 48


def rank_f2(rows, n):
    rows = list(rows)
    r = 0
    for col in range(n):
        piv = None
        for i in range(r, len(rows)):
            if (rows[i] >> col) & 1:
                piv = i
                break
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for i in range(len(rows)):
            if i != r and (rows[i] >> col) & 1:
                rows[i] ^= rows[r]
        r += 1
    return r


def band(n):
    rows = []
    for i in range(n):
        v = 0
        for d in (1, 2):
            if i + d < n:
                v |= 1 << (i + d)
            if i - d >= 0:
                v |= 1 << (i - d)
        rows.append(v)
    return rows


def sums(nmax):
    """exact S_n, D_n for n <= nmax by a transfer over the top-down digits: state = last two digits read."""
    # reading digits from the top: appending digit d below the current digits (a = lowest so far, b = second lowest)
    # adds the terms d * (a + b) to the form (the new lowest digit pairs with the two digits above it).
    out = {}
    cur = {(0, 0): 1}          # value = signed count; state (a, b) = (lowest, second lowest) digits read so far
    for n in range(0, nmax + 1):
        s_n = sum(cur.values())
        d_n = sum(v * (-1 if a else 1) for (a, b), v in cur.items())   # (-1)^(h_0 + q)
        out[n] = (s_n, d_n)
        nxt = {}
        for (a, b), v in cur.items():
            for d in (0, 1):
                sign = -1 if (d and ((a + b) & 1)) else 1
                key = (d, a)
                nxt[key] = nxt.get(key, 0) + sign * v
        cur = nxt
    return out


def main():
    ok = True
    ex = sums(NMAX)
    # brute force
    nb = 20
    size = 1 << nb
    q = np.zeros(size, dtype=np.int64)
    for k in range(nb):
        lo, hi = 1 << k, 1 << (k + 1)
        idx = np.arange(lo, hi, dtype=np.int64)
        half = idx >> 1
        q[lo:hi] = q[half] + np.where((idx & 1) == 1, (half + 1) >> 1, 0)
    h = np.arange(size, dtype=np.int64)
    for n in range(1, nb + 1):
        s_b = int(np.where(q[: 1 << n] % 2 == 0, 1, -1).sum())
        d_b = int(np.where((h[: 1 << n] - q[: 1 << n]) % 2 == 0, 1, -1).sum())
        ok &= (s_b, d_b) == ex[n]
    print("digit recursion for S_n and D_n agrees with brute force for n <= 20:", ok)
    print(f"{'n':>3} {'rank B_n':>8} {'n-rank':>6} {'S_n':>14} {'log2|S_n|':>9} {'D_n':>14} {'log2|D_n|':>9}")
    corank = []
    for n in range(1, NMAX + 1):
        r = rank_f2(band(n), n)
        corank.append(n - r)
        s_n, d_n = ex[n]
        ls = "-inf" if s_n == 0 else f"{np.log2(abs(s_n)):.1f}"
        ld = "-inf" if d_n == 0 else f"{np.log2(abs(d_n)):.1f}"
        if n <= 30:
            print(f"{n:3d} {r:8d} {n - r:6d} {s_n:14d} {ls:>9} {d_n:14d} {ld:>9}")
        # a quadratic form of rank r on F_2^n has |sum| in {0, 2^(n - r/2)}
        for val in (s_n, d_n):
            ok &= val == 0 or abs(val) == 2 ** (n - r // 2)
    # period of the corank
    per = next(p for p in range(1, 25) if all(corank[i] == corank[i + p] for i in range(len(corank) - p)))
    print("corank of B_n for n = 1..24:", corank[:24], " period", per)
    ok &= per == 6
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
