# Auditor recompute-5128, task 5: the constants A_p^sigma = sum_{h >= 0} chi(h) w^sigma(h) of (D2),
#   chi(h) = [p | gcd(x(h), y(h))] - 1/p^2,
#   w^0(h) = (1/x) ln(2 (x + y)/(x + 2y)),   w^1(h) = (1/x) ln(1 + x/(2y)),
# summed by blocks B_n (h with exactly n binary digits; h = 0 is the block n = 0).
# Also the per-padding weights w_s(h) = (1/x) ln((x + 2^s y)/(x/2 + 2^s y)) are checked to sum to w^1 over s >= 1.
import sys
import math
import numpy as np
from numba import njit
from a0_defs import psi_table32


@njit(cache=True)
def block_sums(psi, nb, primes):
    npr = primes.shape[0]
    SW = np.zeros((2, nb + 1), dtype=np.float64)        # sum over block of w^sigma
    DV = np.zeros((2, npr, nb + 1), dtype=np.float64)   # sum over p-divisible h in the block of w^sigma
    CN = np.zeros((npr, nb + 1), dtype=np.int64)        # number of p-divisible h in the block
    for n in range(0, nb + 1):
        lo = 0 if n == 0 else (1 << (n - 1))
        hi = 1 if n == 0 else (1 << n)
        s0 = 0.0
        s1 = 0.0
        c0 = 0.0
        c1 = 0.0
        for h in range(lo, hi):
            x = 2.0 * psi[h] + 1.0
            y = 2.0 * h + 1.0
            w0 = math.log(2.0 * (x + y) / (x + 2.0 * y)) / x
            w1 = math.log1p(x / (2.0 * y)) / x
            # Kahan summation
            t = w0 - c0
            u = s0 + t
            c0 = (u - s0) - t
            s0 = u
            t = w1 - c1
            u = s1 + t
            c1 = (u - s1) - t
            s1 = u
            xi = 2 * np.int64(psi[h]) + 1
            yi = 2 * np.int64(h) + 1
            for j in range(npr):
                p = primes[j]
                if yi % p == 0:
                    if xi % p == 0:
                        DV[0, j, n] += w0
                        DV[1, j, n] += w1
                        CN[j, n] += 1
        SW[0, n] = s0
        SW[1, n] = s1
    return SW, DV, CN


def main():
    nb = int(sys.argv[1]) if len(sys.argv) > 1 else 26
    psi = psi_table32(1 << nb)
    primes = np.array([3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 43, 257], dtype=np.int64)
    claimed = {3: (0.0557, 0.0789), 5: (0.02943, 0.04250), 7: (0.006, 0.007), 11: (0.01096, 0.01449),
               13: (-0.00377, -0.00488), 17: (0.01355, 0.01931), 31: (-0.00194, -0.00261), 43: (0.00235, 0.00323)}
    SW, DV, CN = block_sums(psi, nb, primes)
    # sanity: block sums of w^sigma tend to ln 2 * kappa_sigma (D1)
    import mpmath as mp
    l = [-mp.polylog(2, -mp.mpf(2) ** (-j)) for j in range(0, 4)]
    k0 = float(l[0] - 2 * l[1] + l[2])
    k1 = float(l[1] - l[2])
    print("== block sums of w^sigma over B_n (should tend to ln 2 * kappa_sigma) ==")
    for n in [8, 12, 16, 20, 24, nb]:
        print(f"  n = {n:2d}: sum w^0 = {SW[0, n]:.8f} (ln2*kappa_0 = {math.log(2) * k0:.8f}), sum w^1 = {SW[1, n]:.8f} (ln2*kappa_1 = {math.log(2) * k1:.8f})")
    # check w^1 = sum_{s >= 1} w_s on a few h
    worst = 0.0
    for h in [0, 1, 2, 3, 5, 38, 1000, 123457]:
        x = 2.0 * psi[h] + 1
        y = 2.0 * h + 1
        tot = sum(math.log((x + 2.0 ** s * y) / (x / 2 + 2.0 ** s * y)) / x for s in range(1, 80))
        worst = max(worst, abs(tot - math.log1p(x / (2 * y)) / x))
    print(f"  w^1(h) = sum_(s>=1) w_s(h): max abs difference on 8 test states = {worst:.2e}")
    print()
    print("== partial sums A_p^sigma(n) = sum_{h < 2^n} chi(h) w^sigma(h) ==")
    for j, p in enumerate(primes):
        p = int(p)
        A0 = np.cumsum(DV[0, j] - SW[0] / p ** 2)
        A1 = np.cumsum(DV[1, j] - SW[1] / p ** 2)
        cl = claimed.get(p, None)
        print(f" p = {p}: claimed {cl}")
        print("    n:      " + "  ".join(f"{n:8d}" for n in range(max(10, nb - 12), nb + 1)))
        print("    A^0(n): " + "  ".join(f"{A0[n]:+.5f}" for n in range(max(10, nb - 12), nb + 1)))
        print("    A^1(n): " + "  ".join(f"{A1[n]:+.5f}" for n in range(max(10, nb - 12), nb + 1)))
        last = range(nb - 5, nb + 1)
        print(f"    last 6 partial sums: A^0 in [{min(A0[n] for n in last):+.5f}, {max(A0[n] for n in last):+.5f}], A^1 in [{min(A1[n] for n in last):+.5f}, {max(A1[n] for n in last):+.5f}]; "
              f"mean of last 6: {np.mean([A0[n] for n in last]):+.5f}, {np.mean([A1[n] for n in last]):+.5f}")
        # block terms (the tail behaviour)
        b0 = DV[0, j] - SW[0] / p ** 2
        print("    block terms of A^0, last 8 blocks: " + "  ".join(f"{b0[n]:+.2e}" for n in range(nb - 7, nb + 1)))
        print(f"    p^2 * (density of p-divisible states in block n), last 8 blocks: " + "  ".join(f"{CN[j, n] * p * p / 2 ** (n - 1):.3f}" for n in range(nb - 7, nb + 1)))


if __name__ == "__main__":
    main()
