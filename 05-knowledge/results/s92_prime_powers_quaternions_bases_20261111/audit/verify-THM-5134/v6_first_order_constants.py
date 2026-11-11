#!/usr/bin/env python3
"""Verifier of THM-5134 (G), Numerical support, sentence "Constants": partial sums of the series
A_p^sigma(b) = sum_h chi(h) w^sigma_b(h), chi(h) = [p | gcd(x_b(h), y(h))] - 1/p^2, with the weights of (G),
w^0_b = (2/(b x_b)) ln(2 (x_b + y)/(x_b + 2y)), w^1_b = (2/(b x_b)) ln(1 + x_b/(2y)); and the census value
(D^(p) - #T/p^2)/(2X) at X = 2^26 (own complete census by value).  Also the constants A_(m,d) of (C) at 9, 25, 49
(base 2) for the sentence "The factors".  Own code.
"""
import math
import sys

import numpy as np
from numba import njit

NB = 27


@njit(cache=False)
def psi_table(n):
    ps = np.zeros(n, dtype=np.int64)
    for h in range(1, n):
        k = h >> 1
        if h % 2 == 0:
            ps[h] = ps[k] + k
        else:
            ps[h] = ps[k] + (k >> 1) + 1
    return ps


@njit(cache=False)
def census_types(ps, b, X, sheet, p):
    """ordered triples with N_b <= X on the sheet: (total, number with h1 != h2 whose smaller state is p-divisible)."""
    H = ps.shape[0]
    total = 0
    D = 0
    smin = 0 if sheet == 0 else 1
    smax = 0 if sheet == 0 else 62
    for s in range(smin, smax + 1):
        twoc = b << s
        e0 = b if s == 0 else 0
        if 2 + twoc - e0 > 2 * X:
            break
        for h1 in range(H):
            y1 = 2 * h1 + 1
            lb1 = 1 + b * (h1 >> 1)
            if 2 * lb1 * lb1 + twoc * y1 * y1 - e0 > 2 * X:
                break
            x1 = 1 + b * ps[h1]
            div = (x1 % p == 0) and (y1 % p == 0)
            for h2 in range(h1, H):
                y2 = 2 * h2 + 1
                if 2 * x1 * (1 + b * (h2 >> 1)) + twoc * y1 * y2 - e0 > 2 * X:
                    break
                x2 = 1 + b * ps[h2]
                if 2 * x1 * x2 + twoc * y1 * y2 - e0 > 2 * X:
                    continue
                if h2 == h1:
                    total += 1
                else:
                    total += 2
                    if div:
                        D += 2
    return total, D


def main():
    ok = True
    ps = psi_table(1 << NB)
    h = np.arange(1 << NB, dtype=np.int64)
    y = (2 * h + 1).astype(np.float64)
    print("partial sums of A_p^sigma(b) over h < 2^k and the census value (D^(p) - #T/p^2)/(2X) at X = 2^26")
    want = {(3, 5, 0): ("-0.01338", "-0.01328"), (3, 5, 1): ("-0.01949", "-0.01889"),
            (3, 7, 0): ("0.03517", "0.03548"), (10, 7, 0): ("0.00496", "0.00508")}
    drift = None
    for (b, p) in ((3, 5), (3, 7), (10, 7)):
        xb = (1 + b * ps).astype(np.float64)
        chi = (((1 + b * ps) % p == 0) & ((2 * h + 1) % p == 0)).astype(np.float64) - 1.0 / (p * p)
        w0 = (2.0 / (b * xb)) * np.log(2.0 * (xb + y) / (xb + 2.0 * y))
        w1 = (2.0 / (b * xb)) * np.log(1.0 + xb / (2.0 * y))
        for sheet, w in ((0, w0), (1, w1)):
            cs = np.cumsum(chi * w)
            parts = {k: cs[(1 << k) - 1] for k in (8, 10, 12, 14, 16, 20, 24)}
            tot, D = census_types(ps[: (1 << 25)], b, 1 << 26, sheet, p)
            cen = (D - tot / float(p * p)) / (2.0 * (1 << 26))
            print("  b = %2d, p = %d, sheet %d: series %s | census %+.5f  (#T = %d)"
                  % (b, p, sheet, ", ".join("2^%d: %+.5f" % (k, v) for k, v in parts.items()), cen, tot))
            key = (b, p, sheet)
            if key in want:
                if "%.5f" % parts[12] != want[key][0] or "%.5f" % cen != want[key][1]:
                    ok = False
            if key == (3, 7, 0):
                drift = parts[24]
    print("  the four pairs quoted in the file (series over h < 2^12 against the census) are reproduced: %s" % ok)
    print("  b = 3, p = 7, sheet 0: partial sum over h < 2^24 = %.4f (file: 0.0276)" % drift)
    ok &= "%.4f" % drift == "0.0276"

    print("constants A_(m,d)^sigma of (C), base 2, partial sums over h < 2^27 (file: sentence 'The factors')")
    xs = (2 * ps + 1)
    xf = xs.astype(np.float64)
    w0 = (1.0 / xf) * np.log(2.0 * (xf + y) / (xf + 2.0 * y))
    w1 = (1.0 / xf) * np.log(1.0 + xf / (2.0 * y))
    kap = (0.1615389173, 0.2125139092)
    lnX = 28 * math.log(2.0)

    def cdv(m, d):
        r = 1.0 / (d * d)
        md = m // d
        for q in (3, 5, 7):
            if md % q == 0:
                r *= 1 - 1.0 / (q * q)
        return r

    g = np.gcd(xs, 2 * h + 1)
    res = {}
    for m, divs in ((9, (1, 3, 9)), (25, (1, 5, 25)), (49, (1, 7, 49))):
        gm = np.gcd(g, m)
        for sheet, w in ((0, w0), (1, w1)):
            A = {}
            for d in divs:
                chi = (gm == d).astype(np.float64) - cdv(m, d)
                cs = np.cumsum(chi * w)
                A[d] = cs[-1]
                if m == 49 and d == 7 and sheet == 0:
                    a49 = [cs[(1 << k) - 1] for k in (8, 14, 20, 27)]
            fac = []
            for g0 in divs:
                L = sum(d * cdv(m, d) for d in divs if g0 % d == 0)
                corr = (2.0 / kap[sheet]) * sum(d * A[d] for d in divs if g0 % d == 0)
                fac.append(L + corr / lnX)
            res[(m, sheet)] = fac
            print("  m = %2d sheet %d: A = %s; (D) at X = 2^28: %s"
                  % (m, sheet, ", ".join("d=%d: %+.5f" % (d, A[d]) for d in divs), " / ".join("%.4f" % f for f in fac)))
    c1 = ["%.3f" % f for f in res[(9, 0)]] == ["0.853", "1.220", "1.441"] and \
         ["%.3f" % f for f in res[(9, 1)]] == ["0.850", "1.222", "1.454"]
    c2 = ["%.4f" % f for f in res[(25, 0)]] == ["0.9412", "1.2339", "1.2399"] and \
         ["%.4f" % f for f in res[(25, 1)]] == ["0.9394", "1.2414", "1.2467"]
    c3 = ["%+.4f" % v for v in a49] == ["-0.0151", "+0.0017", "+0.0080", "+0.0054"]
    print("  at 9: 0.853 / 1.220 / 1.441 and 0.850 / 1.222 / 1.454: %s" % c1)
    print("  at 25: 0.9412 / 1.2339 / 1.2399 and 0.9394 / 1.2414 / 1.2467: %s" % c2)
    print("  partial sums of A_(49,7)^0 at H = 2^8, 2^14, 2^20, 2^27: %s : %s" % (", ".join("%+.4f" % v for v in a49), c3))
    ok &= c1 and c2 and c3
    print("ALL CHECKS PASSED" if ok else "FAILURE")
    sys.exit(0 if ok else 1)


main()
