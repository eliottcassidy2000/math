# Auditor recompute-5128, task 4: p-divisible states (p | x(h) and p | y(h)), first occurrences, densities, g = 1 share,
# and the characterisation x = y  <=>  h = 0 or h = 2^j   (theorem part (F)).
import sys
import math
import numpy as np
from numba import njit
from a0_defs import psi_table32, psi_py


@njit(cache=True)
def divisible_list(psi, p, nmax, keep):
    """first `keep` p-divisible states below nmax and counts below 2^k for all k"""
    first = np.full(keep, -1, dtype=np.int64)
    nf = 0
    cnt = np.zeros(64, dtype=np.int64)  # cnt[k] = # divisible h with exactly k binary digits (h = 0 -> k = 0)
    h = (p - 1) // 2
    while h < nmax:
        x = 2 * np.int64(psi[h]) + 1
        if x % p == 0:
            if nf < keep:
                first[nf] = h
                nf += 1
            k = 0
            t = h
            while t > 0:
                k += 1
                t >>= 1
            cnt[k] += 1
        h += p
    return first, cnt


@njit(cache=True)
def gcd_stats(psi, nmax):
    """# h < 2^k with gcd(x, y) = 1, for all k; and the list of h with x == y"""
    cnt = np.zeros(64, dtype=np.int64)
    eq = np.zeros(64, dtype=np.int64)
    neq = 0
    for h in range(nmax):
        x = 2 * np.int64(psi[h]) + 1
        y = 2 * np.int64(h) + 1
        if x == y:
            if neq < 64:
                eq[neq] = h
            neq += 1
        a = x
        b = y
        while b != 0:
            a, b = b, a % b
        if a == 1:
            k = 0
            t = h
            while t > 0:
                k += 1
                t >>= 1
            cnt[k] += 1
    return cnt, eq, neq


def main():
    nb = int(sys.argv[1]) if len(sys.argv) > 1 else 26
    nmax = 1 << nb
    psi = psi_table32(nmax)
    claimed_first = {3: 1, 5: 2, 7: 38, 11: 16, 13: 32, 17: 8, 19: 218, 23: 103, 29: 1435, 31: 8664, 37: 1350, 41: 266,
                     43: 64, 127: 109283, 257: 128, 1039: 405729, 65537: 32768}
    primes = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 61, 73, 127, 257, 1039, 2731, 65537]
    print(f"== first p-divisible states and densities (h < 2^{nb}) ==")
    print("  p   first states (mine)                          claimed first   # below 2^24   p^2 * density below 2^24   # below 2^26   p^2 * density below 2^26")
    for p in primes:
        first, cnt = divisible_list(psi, p, nmax, 6)
        c24 = int(cnt[:25].sum())
        c26 = int(cnt[:27].sum()) if nb >= 26 else -1
        fl = [int(v) for v in first if v >= 0]
        # independent confirmation of the first state by the pure-python recursion (exhaustive scan below it)
        if fl and fl[0] < 2_000_000 and p < 3000:
            h0 = None
            for h in range(fl[0] + 1):
                if (2 * h + 1) % p == 0 and (2 * psi_py(h) + 1) % p == 0:
                    h0 = h
                    break
            assert h0 == fl[0], (p, h0, fl[0])
        cl = claimed_first.get(p, None)
        flag = "" if cl is None else ("OK" if cl == fl[0] else "MISMATCH")
        print(f" {p:6d}  {str(fl):44s}  {str(cl):>8s} {flag:8s}  {c24:9d}   {c24 * p * p / 2 ** 24:.5f}    {c26:9d}   {c26 * p * p / 2 ** 26:.5f}")
    print()
    print("== the states h = 2^(m-1) with x = y = 2^m + 1 and the primes dividing 2^m + 1 (part (F)) ==")
    for m in range(1, 17):
        h = 1 << (m - 1)
        x = 2 * psi_py(h) + 1
        y = 2 * h + 1
        from sympy import factorint
        print(f"  m = {m:2d}: h = {h:6d}, x = {x}, y = {y}, x == y == 2^m + 1: {x == y == 2 ** m + 1}; 2^m + 1 = {dict(factorint(2 ** m + 1))}")
    cnt, eq, neq = gcd_stats(psi, nmax)
    eql = [int(v) for v in eq[:min(neq, 64)]]
    print(f"  h < 2^{nb} with x(h) = y(h): {eql}")
    print(f"  all of them are 0 or powers of two: {all(h == 0 or (h & (h - 1)) == 0 for h in eql)}; count {neq} (expected {nb + 1})")
    print()
    print("== share of h with g(h) = gcd(x, y) = 1 ==")
    for k in [16, 20, 22, 24, 26]:
        if k <= nb:
            c = int(cnt[:k + 1].sum())
            print(f"  h < 2^{k}: {c} / {2 ** k} = {c / 2 ** k:.6f}   (8/pi^2 = {8 / math.pi ** 2:.6f}; claimed 0.81176 at 2^24)")


if __name__ == "__main__":
    main()
