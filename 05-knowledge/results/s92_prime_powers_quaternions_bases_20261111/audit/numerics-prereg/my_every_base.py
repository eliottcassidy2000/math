# Auditor numerics-prereg: independent base-b census of F1 (THM-5134 (H), HYP-9665 (C)); written before
# reading every_base.py / section S10 of bridges_exact.py.
#   N_b(h1,h2,s) = X1 X2 + b (2^s y1 y2 - [s=0]) / 2,   X = 1 + b psi(h),  y = 2h+1
#   alpha_1 = psi1 + psi2 + (2^s y1 y2 - [s=0])/2,  alpha_2 = psi1 psi2
# sheet 1 (s >= 1): distinguished class N_b = 0 (mod p); sheet 0: 2 N_b + b = 0 (mod p).
# First-order law at a prime p not dividing b:  #{class c} = D1/p + [c distinguished] Dp, where D1, Dp count the
# ordered triples with h1 != h2 whose smaller state is not / is p-divisible (p | gcd(X, y)).
# usage: python3 -B my_every_base.py K1 K0     (N_b <= 2^K1 on sheet 1, N_b <= 2^K0 on sheet 0)
import sys
import time
import numpy as np
import numba as nb
from aud_common import psi_table, gcd_nb

PRIMES = np.array([3, 5, 7, 11, 13], dtype=np.int64)


@nb.njit
def census_b(psi, b, X, sheet, primes):
    npr = primes.shape[0]
    cls = np.zeros((npr, 16), np.int64)      # class of N_b mod p
    Dp = np.zeros(npr, np.int64)             # off-diagonal ordered triples with p-divisible smaller state
    viol = np.zeros(npr, np.int64)           # such triples outside the distinguished class
    dig = np.zeros(b, np.int64)              # alpha_1 mod b
    tot = 0
    ndiag = 0
    n1bad = 0
    s = 0 if sheet == 0 else 1
    while True:
        e = 1 if s == 0 else 0
        # smallest value: h1 = h2 = 0: 1 + b (2^s - e)/2
        if 1 + (b * ((1 << s) - e)) // 2 > X:
            break
        h1 = 0
        while True:
            y1 = 2 * h1 + 1
            # trivial bound only: N_b >= b (2^s y1 y2 - e)/2 >= b (2^s y1 y1 - e)/2 for h2 >= h1
            if (b * ((y1 << s) * y1 - e)) // 2 > X:
                break
            X1 = 1 + b * psi[h1]
            g1 = gcd_nb(X1, y1)
            h2 = h1
            while True:
                y2 = 2 * h2 + 1
                c3 = ((y1 << s) * y2 - e) // 2
                if b * c3 > X:
                    break
                X2 = 1 + b * psi[h2]
                N = X1 * X2 + b * c3
                if N <= X:
                    w = 2 if h2 > h1 else 1
                    tot += w
                    if h2 == h1:
                        ndiag += 1
                    if N % b != 1:
                        n1bad += 1
                    a1 = psi[h1] + psi[h2] + c3
                    dig[a1 % b] += w
                    for i in range(npr):
                        p = primes[i]
                        cls[i, N % p] += w
                        if h2 > h1 and g1 % p == 0:
                            Dp[i] += 2
                            if sheet == 1:
                                if N % p != 0:
                                    viol[i] += 1
                            else:
                                if (2 * N + b) % p != 0:
                                    viol[i] += 1
                h2 += 1
            h1 += 1
        if sheet == 0:
            break
        s += 1
    return cls, Dp, viol, dig, tot, ndiag, n1bad


def main():
    K1 = int(sys.argv[1])
    K0 = int(sys.argv[2])
    bases1 = [2, 3, 4, 5, 6, 10]
    bases0 = [2, 3, 5, 10]
    if len(sys.argv) > 3:
        bases1 = [int(v) for v in sys.argv[3].split(",")]
        bases0 = []
    # largest h2 that can occur: sheet 1: N_b >= b y2, so h2 <= X/(2b); sheet 0: N_b >= b (y2 - 1)/2, so h2 <= X/b
    need = (1 << K1) // (2 * min(bases1)) + 8
    if bases0:
        need = max(need, (1 << K0) // min(bases0) + 8)
    psi = psi_table(need)
    for sheet, K, bases in ((1, K1, bases1), (0, K0, bases0)):
        X = 1 << K
        print("sheet %s, N_b <= 2^%d" % ("s >= 1" if sheet == 1 else "s = 0", K))
        for b in bases:
            t0 = time.time()
            cls, Dp, viol, dig, tot, ndiag, n1bad = census_b(psi, b, X, sheet, PRIMES)
            flat = max(abs(b * dig[c] / tot - 1.0) for c in range(b))
            cells = []
            worst = 0.0
            for i, p in enumerate(PRIMES):
                p = int(p)
                if b % p == 0:
                    cells.append("p=%d: N = 1 (mod p)" % p)
                    continue
                cdist = 0 if sheet == 1 else (-b * pow(2, -1, p)) % p
                D1 = tot - ndiag - Dp[i]
                rem = 0.0
                for c in range(p):
                    pred = D1 / p + (Dp[i] if c == cdist else 0)
                    rem = max(rem, p * abs(cls[i, c] - pred) / tot)
                worst = max(worst, rem)
                others = [p * cls[i, c] / tot for c in range(p) if c != cdist]
                cells.append("p=%d: class %d factor %.3f [%.3f, %.3f] F p^2=%.2f rem=%.3f viol=%d" % (
                    p, cdist, p * cls[i, cdist] / tot, min(others), max(others), Dp[i] / tot * p * p, rem, viol[i]))
            print("    b = %2d: %10d triples (diagonal %d); N = 1 (mod b) fails %d times; second digit: max_c |b share(alpha_1 = c) - 1| = %.3f; largest rem %.3f; %s  (%.0f s)" % (
                b, tot, ndiag, n1bad, flat, worst, "; ".join(cells), time.time() - t0))


if __name__ == "__main__":
    main()
