#!/usr/bin/env python3
"""Audit proof-5128-AB: THM-5128 (B3), exact integer check for every h < 2^NB.

Claims checked (exact, int64; all quantities below 2^62):
  (a) 0 <= q(h) <= h/2 (so y/2 < x <= y);
  (b) q(h) = sum over one-bits i of ceil(H_i / 2), H_i = floor(h / 2^(i+1))        (for h < 2^16, pure Python);
  (c) |rho(h) - (1 - phi(h)/2)| < n/(2h + 1) on B_n, i.e. with rev = 2^n phi(h) (bit reversal on n digits):
          | (2h + 1) rev - 4 q 2^n |  <  2 n 2^n ;
      the proof's finer statement: S(h) := (2h+1) phi - 4q lies in (-n, 2n), every bracket in (-1, 2);
  (d) the radical inverse takes each value (2i + 1)/2^n once on B_n.
Also reported: the extreme values of S(h)/n on each block (how sharp the bound is), and the edge case h = 0.
"""
import sys

import numpy as np


def q_table(nb):
    size = 1 << nb
    q = np.zeros(size, dtype=np.int64)
    for k in range(nb):
        lo, hi = 1 << k, 1 << (k + 1)
        idx = np.arange(lo, hi, dtype=np.int64)
        half = idx >> 1
        q[lo:hi] = q[half] + (idx & 1) * ((half + 1) >> 1)
    return q


def bitrev(h, n):
    r = np.zeros_like(h)
    t = h.copy()
    for _ in range(n):
        r = (r << 1) | (t & 1)
        t >>= 1
    return r


def main():
    nb = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    q = q_table(nb)
    h = np.arange(1 << nb, dtype=np.int64)
    ok = True
    # (a)
    ok &= bool(np.all(q >= 0) and np.all(2 * q <= h))
    x = 2 * h + 1 - 2 * q
    y = 2 * h + 1
    ok &= bool(np.all(2 * x > y) and np.all(x <= y))
    print(f"(a) 0 <= q <= h/2 and y/2 < x <= y for all h < 2^{nb}: {ok}")
    # (b) pure python
    def q_bits(hh):
        s, i = 0, 0
        while (hh >> i):
            if (hh >> i) & 1:
                H = hh >> (i + 1)
                s += (H + 1) // 2
            i += 1
        return s
    okb = all(q_bits(hh) == int(q[hh]) for hh in range(1 << 16))
    ok &= okb
    print(f"(b) q(h) = sum over one-bits of ceil(H_i/2), h < 2^16: {okb}")
    # brackets, pure python, h < 2^16: every bracket in (-1, 2), exactly
    from fractions import Fraction
    okbr = True
    for hh in range(1, 1 << 12):
        i = 0
        while (hh >> i):
            if (hh >> i) & 1:
                H = hh >> (i + 1)
                br = Fraction(2 * hh + 1, 1 << (i + 1)) - 4 * ((H + 1) // 2)
                first = Fraction(2 * hh + 1, 1 << (i + 1))
                okbr &= (2 * H + 1 < first < 2 * H + 2) and (-1 < br < 2)
                okbr &= (1 < br < 2) if H % 2 == 0 else (-1 < br < 0)
            i += 1
    ok &= okbr
    print(f"    every bracket in (-1, 2) (in (1, 2) for H_i even, in (-1, 0) for H_i odd), h < 2^12, exact: {okbr}")
    # (c), (d)
    print(f"{'n':>3} {'max |S|/(2n)':>13} {'min S/n':>9} {'max S/n':>9} {'max (2h+1)|rho-1+phi/2|/n':>26}  bound holds")
    for n in range(1, nb + 1):
        lo, hi = 1 << (n - 1), 1 << n
        hh = h[lo:hi]
        rev = bitrev(hh, n)
        # (d)
        if n <= 20:
            srt = np.sort(rev)
            ok &= bool(np.array_equal(srt, np.arange(1, 1 << n, 2, dtype=np.int64)))
        Sn = (2 * hh + 1) * rev - 4 * q[lo:hi] * (1 << n)          # = 2^n * S(h)
        holds = bool(np.all(np.abs(Sn) < 2 * n * (1 << n)))
        finer = bool(np.all(Sn > -n * (1 << n)) and np.all(Sn < 2 * n * (1 << n)))
        ok &= holds and finer
        smin = float(Sn.min()) / (1 << n) / n
        smax = float(Sn.max()) / (1 << n) / n
        worst = float(np.max(np.abs(Sn))) / (1 << n) / (2 * n)
        if n <= 12 or n % 4 == 0 or n == nb:
            print(f"{n:3d} {worst:13.6f} {smin:9.4f} {smax:9.4f} {worst:26.6f}  {holds and finer}")
    # edge case
    print("edge case h = 0 (zero digits, n = 0): rho = 1, phi = 0, |rho - (1 - phi/2)| = 0, bound n/(2h+1) = 0: the strict"
          " inequality fails there (the statement needs n >= 1)")
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
