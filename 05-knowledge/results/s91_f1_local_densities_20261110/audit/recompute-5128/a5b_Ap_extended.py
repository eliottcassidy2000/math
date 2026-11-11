# Auditor recompute-5128, task 5 (extension): partial sums of A_p^sigma far beyond 2^26, for slowly mixing primes.
# Only the p-divisible states are visited (h = (p-1)/2 mod p), with q(h) computed without a table of size 2^n:
#   for h = a 2^16 + b (0 <= b < 2^16):  q(h) = q(a) + q(b) + a rev15(b mod 2^15) + [b >= 2^15] ceil(a/2),
# rev15 = reversal of the low 15 bits (from q(h) = sum over one-bits i of ceil(floor(h/2^(i+1))/2)).
# The block sums of w^sigma over ALL h are computed directly up to 2^26 and replaced by their limits ln 2 * kappa_sigma
# beyond (they agree to 2e-8 at n = 24).
# Usage: python3 a5b_Ap_extended.py p nmax [threads]
import sys
import math
import time
import numpy as np
from numba import njit
from concurrent.futures import ThreadPoolExecutor
from a0_defs import psi_table32


@njit(cache=True)
def rev15_table():
    r = np.zeros(1 << 15, dtype=np.int64)
    for b in range(1 << 15):
        v = 0
        for i in range(15):
            if (b >> i) & 1:
                v |= 1 << (14 - i)
        r[b] = v
    return r


@njit(cache=True)
def check_split(psi, rev15, nmax_check):
    bad = 0
    for h in range(1 << 16, 1 << nmax_check, 7):
        a = h >> 16
        b = h & 65535
        q = (a - psi[a]) + (b - psi[b]) + a * rev15[b & 32767]
        if b >= 32768:
            q += (a + 1) >> 1
        if q != h - psi[h]:
            bad += 1
    return bad


@njit(cache=True, nogil=True)
def div_sums(qlow, qa, rev15, p, a_lo, a_hi, nmax):
    DV = np.zeros((2, nmax + 1), dtype=np.float64)
    CN = np.zeros(nmax + 1, dtype=np.int64)
    r0 = (p - 1) // 2
    two16 = 65536 % p
    for a in range(a_lo, a_hi):
        n = 16
        t = a
        while t > 0:
            n += 1
            t >>= 1
        # b = r0 - a 2^16 (mod p)
        b0 = (r0 - (a % p) * two16) % p
        qa_ = qa[a]
        ca = (a + 1) >> 1
        base = np.int64(a) << 16
        s0 = 0.0
        s1 = 0.0
        c = 0
        b = b0
        while b < 65536:
            q = qa_ + qlow[b] + a * rev15[b & 32767]
            if b >= 32768:
                q += ca
            h = base + b
            y = 2 * h + 1
            x = y - 2 * q
            if x % p == 0:
                xf = float(x)
                yf = float(y)
                s0 += math.log(2.0 * (xf + yf) / (xf + 2.0 * yf)) / xf
                s1 += math.log1p(xf / (2.0 * yf)) / xf
                c += 1
            b += p
        DV[0, n] += s0
        DV[1, n] += s1
        CN[n] += c
    return DV, CN


@njit(cache=True)
def low_blocks(psi, p, nb):
    """blocks n <= nb by direct summation: SW[sigma, n], DV[sigma, n], CN[n]"""
    SW = np.zeros((2, nb + 1), dtype=np.float64)
    DV = np.zeros((2, nb + 1), dtype=np.float64)
    CN = np.zeros(nb + 1, dtype=np.int64)
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
            t = w0 - c0
            u = s0 + t
            c0 = (u - s0) - t
            s0 = u
            t = w1 - c1
            u = s1 + t
            c1 = (u - s1) - t
            s1 = u
            yi = 2 * np.int64(h) + 1
            if yi % p == 0:
                xi = 2 * np.int64(psi[h]) + 1
                if xi % p == 0:
                    DV[0, n] += w0
                    DV[1, n] += w1
                    CN[n] += 1
        SW[0, n] = s0
        SW[1, n] = s1
    return SW, DV, CN


def main():
    p = int(sys.argv[1])
    nmax = int(sys.argv[2])
    nthreads = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    nlow = 26
    t0 = time.time()
    psi = psi_table32(1 << max(nlow, nmax - 16))
    rev15 = rev15_table()
    bad = check_split(psi, rev15, 24)
    SW, DVl, CNl = low_blocks(psi, p, nlow)
    import mpmath as mp
    l = [-mp.polylog(2, -mp.mpf(2) ** (-j)) for j in range(0, 4)]
    lim = (math.log(2) * float(l[0] - 2 * l[1] + l[2]), math.log(2) * float(l[1] - l[2]))
    print(f"# p = {p}: split formula for q checked against the table on a sample of h < 2^24: {bad} mismatches; low blocks done in {time.time() - t0:.1f} s")
    print(f"# block sums of w at n = {nlow}: {SW[0, nlow]:.10f}, {SW[1, nlow]:.10f}; limits {lim[0]:.10f}, {lim[1]:.10f}")
    # high blocks through the split
    qlow = (np.arange(1 << 16, dtype=np.int64) - psi[: 1 << 16].astype(np.int64))
    na = 1 << (nmax - 16)
    qa = (np.arange(na, dtype=np.int64) - psi[:na].astype(np.int64))
    DV = np.zeros((2, nmax + 1))
    CN = np.zeros(nmax + 1, dtype=np.int64)
    # a from 1 .. na-1 covers h in [2^16, 2^nmax); chunk by ranges, interleaved for balance
    nchunks = 48
    edges = np.linspace(1, na, nchunks + 1).astype(np.int64)
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=nthreads) as ex:
        futs = [ex.submit(div_sums, qlow, qa, rev15, p, int(edges[i]), int(edges[i + 1]), nmax) for i in range(nchunks)]
        for f in futs:
            d, c = f.result()
            DV += d
            CN += c
    print(f"# high blocks done in {time.time() - t0:.1f} s")
    # consistency between the two methods on blocks 17..nlow
    worst = 0.0
    for n in range(17, nlow + 1):
        worst = max(worst, abs(DV[0, n] - DVl[0, n]), abs(DV[1, n] - DVl[1, n]))
        assert CN[n] == CNl[n], (n, CN[n], CNl[n])
    print(f"# divisible-state sums by the split formula agree with the direct ones on blocks 17..{nlow}: max abs difference {worst:.2e}; counts identical")
    A0 = 0.0
    A1 = 0.0
    print("  n    A^0(n)     A^1(n)     block term A^0   block term A^1   p^2 * density of divisible states in the block")
    for n in range(0, nmax + 1):
        if n <= nlow:
            b0 = DVl[0, n] - SW[0, n] / p ** 2
            b1 = DVl[1, n] - SW[1, n] / p ** 2
            cn = CNl[n]
        else:
            b0 = DV[0, n] - lim[0] / p ** 2
            b1 = DV[1, n] - lim[1] / p ** 2
            cn = CN[n]
        A0 += b0
        A1 += b1
        if n >= 14:
            dens = cn * p * p / (2 ** (n - 1)) if n >= 1 else 0.0
            print(f"  {n:2d}  {A0:+.6f}  {A1:+.6f}   {b0:+.3e}       {b1:+.3e}       {dens:.4f}")


if __name__ == "__main__":
    main()
