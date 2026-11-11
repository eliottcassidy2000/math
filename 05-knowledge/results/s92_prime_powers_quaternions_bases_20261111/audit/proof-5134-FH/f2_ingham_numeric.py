#!/usr/bin/env python3
"""Auditor proof-5134-FH.  THM-5134 (F): numerical illustration of Ingham's theorem (CITED, not proved here)
    D(n) = sum_{m=1}^{n-1} d(m) d(n-m) = #{(a1,a2,b1,b2) positive : a1 a2 + b1 b2 = n} ~ (6/pi^2) sigma(n) ln^2 n
and of the local reading: (6/pi^2) sigma(n)/n = prod_p (1 - p^-2) sigma(p^v)/p^v, the product of the p-adic densities
of the determinant form a1 a2 + b1 b2 at n (all primes, 2 included: 3/4 at odd n).
Own code.  (i) D(n) counted directly for small n by brute force over (a1, a2, b1, b2) to check the identity with
the divisor convolution; (ii) the ratio D(n) / ((6/pi^2) sigma(n) ln^2 n) on windows near 2*10^4, 2*10^5, 10^6.
"""
import math
import numpy as np


def main():
    # (i) the identity, brute force
    ok = True
    for n in range(2, 121):
        cnt = 0
        for a1 in range(1, n):
            for a2 in range(1, (n - 1) // a1 + 1):
                r = n - a1 * a2
                if r <= 0:
                    break
                # number of (b1, b2) positive with b1 b2 = r
                cnt += sum(1 for b1 in range(1, r + 1) if r % b1 == 0)
        dv = [0] * (n + 1)
        for d in range(1, n + 1):
            for k in range(d, n + 1, d):
                dv[k] += 1
        conv = sum(dv[m] * dv[n - m] for m in range(1, n))
        ok &= (cnt == conv)
    print(f"(i) #{{a1 a2 + b1 b2 = n, positive}} = sum_{{m<n}} d(m) d(n-m) for 2 <= n <= 120 (brute force): {ok}")

    # (ii) windows
    TOP = 1000400
    dv = np.zeros(TOP + 1, dtype=np.int64)
    sg = np.zeros(TOP + 1, dtype=np.int64)
    for d in range(1, TOP + 1):
        dv[d::d] += 1
        sg[d::d] += d
    c = 6 / math.pi ** 2
    print("(ii) ratio D(n) / ((6/pi^2) sigma(n) ln^2 n) on windows of 400 consecutive n; Ingham's error term is O(sigma(n) ln n ln ln n),")
    print("     i.e. the ratio is 1 + O(ln ln n / ln n)")
    for N0 in (20000, 200000, 1000000):
        ns = np.arange(N0, N0 + 400)
        D = np.array([int(np.dot(dv[1:n], dv[n - 1:0:-1])) for n in ns], dtype=float)
        ratio = D / (c * sg[ns] * np.log(ns) ** 2)
        corr = np.corrcoef(D / ns, sg[ns] / ns)[0, 1]
        odd = (ns % 2 == 1)
        print(f"   n in [{N0}, {N0 + 400}): mean {ratio.mean():.4f}, min {ratio.min():.4f}, max {ratio.max():.4f}, sd {ratio.std():.4f}; "
              f"ln ln n / ln n = {math.log(math.log(N0)) / math.log(N0):.3f}; corr(D(n)/n, sigma(n)/n) = {corr:.4f}; "
              f"odd n only: mean {ratio[odd].mean():.4f}, sd {ratio[odd].std():.4f}")
    print("constants: 6/pi^2 = %.10f = 1/zeta(2);  (6/pi^2)/(1 - 1/4) = %.10f;  8/pi^2 = %.10f" % (c, c / 0.75, 8 / math.pi ** 2))


if __name__ == "__main__":
    main()
