#!/usr/bin/env python3
"""Auditor proof-5134-FH.  THM-5134 (H)(ii): does statement (C) transfer to base b with the weights obtained by the
substitution a = b X1/2, b' = b 2^(s-1) y1 in THM-5128 (D1), (D2)?

Weights (derived by the auditor; for b = 2 they are the w^0, w^1 of THM-5128 (D2)):
    w^0_b(h) = (2/(b X)) ln(2 (X + y)/(X + 2 y)),      w^1_b(h) = (2/(b X)) ln(1 + X/(2 y)),     X = 1 + b psi(h), y = 2h + 1.
Constant:  A_p^sigma(b) = lim_H sum_{h < H} chi(h) w^sigma_b(h),   chi(h) = [p | X(h) and p | y(h)] - 1/p^2.
Prediction of (C) in base b:   D^(p)(X) = #T(X)/p^2 + 2 A_p^sigma(b) X + lower order,
so the census value  (D^(p) - #T/p^2)/(2X)  should approach the partial sums of the series (slowly: at X the census
sees the smaller states with y1 < (X/2^s)^(1/2) only, less a boundary term - THM-5128 (D4)).

Control: b = 2 against the table of THM-5128 (D4).
Census values: read from the auditor's complete censuses h3_census_26_sheet{0,1}.out (b = 3, 10; X = 2^26).
"""
import math
import os
import re
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def psi_table(nbits):
    q = np.zeros(1 << nbits, dtype=np.int64)
    for k in range(nbits):
        lo, hi = 1 << k, 1 << (k + 1)
        idx = np.arange(lo, hi, dtype=np.int64)
        half = idx >> 1
        q[lo:hi] = q[half] + np.where((idx & 1) == 1, (half + 1) >> 1, 0)
    return np.arange(1 << nbits, dtype=np.int64) - q


def partial_sums(psi, b, p, cuts):
    h = np.arange(psi.shape[0], dtype=np.int64)
    X = 1 + b * psi
    y = 2 * h + 1
    chi = ((X % p == 0) & (y % p == 0)).astype(np.float64) - 1.0 / (p * p)
    Xf = X.astype(np.float64)
    yf = y.astype(np.float64)
    w0 = (2.0 / (b * Xf)) * np.log(2.0 * (Xf + yf) / (Xf + 2.0 * yf))
    w1 = (2.0 / (b * Xf)) * np.log1p(Xf / (2.0 * yf))
    c0 = np.cumsum(chi * w0)
    c1 = np.cumsum(chi * w1)
    return [(float(c0[H - 1]), float(c1[H - 1])) for H in cuts]


def read_census(fname):
    """returns {(b, m): (T, D_m)} from an h3 census output."""
    out = {}
    b = None
    T = None
    with open(os.path.join(HERE, fname)) as f:
        for line in f:
            m1 = re.match(r"b = (\d+): #T = (\d+)", line)
            if m1:
                b, T = int(m1.group(1)), int(m1.group(2))
                continue
            m2 = re.match(r"\s+m =\s+(\d+) \(prime\).*types 1: (\d+), (\d+): (\d+)", line)
            if m2 and b is not None:
                m = int(m2.group(1))
                assert int(m2.group(3)) == m
                out[(b, m)] = (T, int(m2.group(4)))
    return out


def main():
    NB = 24
    psi = psi_table(NB)
    cuts = [1 << k for k in (8, 10, 12, 14, 16, 20, 24)]
    print("partial sums of A_p^0(b), A_p^1(b) over h < 2^k, k = 8, 10, 12, 14, 16, 20, 24")
    print("control b = 2 (THM-5128 (D4), states below 2^26: p=3: 0.0556, 0.0789; p=5: 0.02943, 0.04250; p=11: 0.01096, 0.01449; "
          "p=13: -0.00377, -0.00488; p=17: 0.01355, 0.01931)")
    for p in (3, 5, 11, 13, 17):
        ps = partial_sums(psi, 2, p, cuts)
        print(f"  b = 2, p = {p:2d}: sheet 0: " + ", ".join(f"{a:.5f}" for a, _ in ps) + " | sheet 1: " + ", ".join(f"{c:.5f}" for _, c in ps))
    X = float(1 << 26)
    cens = {0: read_census("h3_census_26_sheet0.out"), 1: read_census("h3_census_26_sheet1.out")}
    print("\nbases 3 and 10: partial sums against the complete census at X = 2^26, (D^(p) - #T/p^2)/(2X)")
    for b in (3, 10):
        for p in (3, 5, 7, 11, 13):
            if b % p == 0:
                continue
            ps = partial_sums(psi, b, p, cuts)
            cv = []
            for sheet in (0, 1):
                T, D = cens[sheet][(b, p)]
                cv.append((D - T / (p * p)) / (2 * X))
            print(f"  b = {b:2d}, p = {p:2d}: sheet 0: series " + ", ".join(f"{a:.5f}" for a, _ in ps) + f" | census {cv[0]:+.5f}")
            print(f"                  sheet 1: series " + ", ".join(f"{c:.5f}" for _, c in ps) + f" | census {cv[1]:+.5f}")


if __name__ == "__main__":
    main()
