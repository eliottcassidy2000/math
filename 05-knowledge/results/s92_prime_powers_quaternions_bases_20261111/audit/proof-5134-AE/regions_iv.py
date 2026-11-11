# regions_iv.py -- auditor proof-5134-AE.
# Direct test of THM-5128 (B4)(ii) in the form (iv) at prime-power and composite odd moduli, which THM-5134 (A) uses:
#   #{h in B_n cap I cap Omega : (x, y) = w mod m} = m^-2 #(B_n cap I cap Omega) + O_m(2^(n (1 - theta_m)))
# for EVERY residue vector w modulo m, with Omega = {h : a x(h) + b y(h) <= X}, a = x1, b = 2^s y1, I = (h1, inf).
# Printed per block n: the block count, the largest deviation over the m^2 residue vectors, and the exponent
#   e_n = log2(max deviation) / n   (the lemma asserts e_n <= 1 - theta_m + o(1); trivial bound e_n = 1).
import math
import numpy as np
import numba as nb


@nb.njit(cache=True)
def build_psi(H):
    psi = np.zeros(H + 1, np.int32)
    for h in range(1, H + 1):
        k = h >> 1
        if (h & 1) == 0:
            psi[h] = psi[k] + k
        else:
            psi[h] = psi[k] + (k >> 1) + 1
    return psi


@nb.njit(cache=True)
def block_counts(psi, n, h1, a, b, X, m):
    cnt = np.zeros(m * m, np.int64)
    lo = 1 << (n - 1)
    hi = 1 << n
    if lo <= h1:
        lo = h1 + 1
    for h in range(lo, hi):
        x = np.int64(2 * psi[h] + 1)
        y = np.int64(2 * h + 1)
        if a * x + b * y <= X:
            cnt[(x % m) * m + (y % m)] += 1
    return cnt


def main():
    NMAX = 26
    psi = build_psi(1 << NMAX)
    rates = {9: 0.866025, 27: 0.866025, 25: 0.815263, 45: 0.892191, 49: 0.932039}
    for m in (9, 27, 25, 45, 49):
        th = math.log(1 / rates[m]) / math.log(2 / rates[m] ** 2)
        print("modulus %d  (1 - theta_m > %.4f)" % (m, 1 - th))
        for (h1, s) in ((0, 0), (1, 1), (3, 0), (5, 2)):
            x1 = 2 * int(psi[h1]) + 1
            y1 = 2 * h1 + 1
            a = x1
            b = (1 << s) * y1
            row = []
            for n in (14, 18, 22, 26):
                # choose X so that the region cuts the block B_n roughly in the middle: a x + b y ~ y (a rho + b)
                X = int((a * 0.75 + b) * (1.5 * 2 ** n))
                cnt = block_counts(psi, n, h1, a, b, X, m)
                tot = int(cnt.sum())
                dev = float(np.abs(cnt - tot / (m * m)).max())
                row.append("n=%d: #=%d maxdev=%.1f (rel %.3f) e=%.3f" % (n, tot, dev, dev / (tot / (m * m)), math.log2(max(dev, 1e-9)) / n))
            print("   (h1, s) = (%d, %d), a = %d, b = %d | " % (h1, s, a, b) + " | ".join(row))


if __name__ == "__main__":
    main()
