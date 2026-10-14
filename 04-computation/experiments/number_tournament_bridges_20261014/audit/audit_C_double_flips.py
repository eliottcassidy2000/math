#!/usr/bin/env python3
"""Independent audit of the 'double flips' (odd a with s(a^2) a square or twice a square), referee script.

sigma(a^2) is computed by a SEGMENTED trial-division sieve over odd a (the author used a full smallest-prime-factor
table); exact integer square tests; then decade counts, prime/composite split, the amicable closure test
sigma(s(a^2)) = sigma(a^2) with s(a^2) = y^2 or 2 y^2 factored through y, the max ratio, the p = 71 conic, the count
of primes 2 b^2 - 1, and the int64 bound.
"""
import sys, time, math
import numpy as np
import numba as nb

A = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10 ** 9


def small_primes(B):
    sieve = bytearray([1]) * (B + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, int(B ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    return np.array([i for i in range(B + 1) if sieve[i]], dtype=np.int64)


@nb.njit(cache=True)
def isqrt64(m):
    if m < 0:
        return -1
    r = np.int64(math.sqrt(float(m)))
    while r * r > m:
        r -= 1
    while (r + 1) * (r + 1) <= m:
        r += 1
    return r


@nb.njit(parallel=True, cache=True)
def scan(A, oddprimes, seg):
    """returns a flag array over odd a (index i <-> a = 2 i + 1) marking s(a^2) in squares or twice squares,
    plus the value s(a^2) for flagged entries (0 elsewhere) and the max of sigma(a^2)/a^2 seen"""
    n_odd = (A + 1) // 2                       # a = 1, 3, ..., <= A
    flag = np.zeros(n_odd, dtype=np.uint8)
    nseg = (n_odd + seg - 1) // seg
    maxrat = np.zeros(nseg, dtype=np.float64)
    for si in nb.prange(nseg):
        i0 = si * seg
        i1 = min(i0 + seg, n_odd)
        L = 2 * i0 + 1
        R = 2 * (i1 - 1) + 1                   # last a in segment
        rem = np.empty(i1 - i0, dtype=np.int64)
        sg = np.ones(i1 - i0, dtype=np.int64)
        for i in range(i0, i1):
            rem[i - i0] = 2 * i + 1
        lim = isqrt64(R)
        for t in range(oddprimes.shape[0]):
            p = oddprimes[t]
            if p > lim:
                break
            # first odd multiple of p >= L
            a0 = ((L + p - 1) // p) * p
            if a0 % 2 == 0:
                a0 += p
            for a in range(a0, R + 1, 2 * p):
                j = (a - 1) // 2 - i0
                e = 0
                x = rem[j]
                while x % p == 0:
                    x //= p
                    e += 1
                rem[j] = x
                tsum = np.int64(1)
                pw = np.int64(1)
                for _ in range(2 * e):
                    pw *= p
                    tsum += pw
                sg[j] *= tsum
        lmax = 0.0
        for j in range(i1 - i0):
            a = np.int64(2 * (i0 + j) + 1)
            x = rem[j]
            if x > 1:
                sg[j] *= 1 + x + x * x
            a2 = a * a
            m = sg[j] - a2
            rat = sg[j] / float(a2)
            if rat > lmax:
                lmax = rat
            if a >= 3:
                r = isqrt64(m)
                ok = r * r == m
                if not ok and m % 2 == 0:
                    r2 = isqrt64(m // 2)
                    ok = r2 * r2 == m // 2
                if ok:
                    flag[i0 + j] = 1
        maxrat[si] = lmax
    return flag, maxrat


def factor_small(n, primes):
    f = {}
    m = n
    for p in primes:
        p = int(p)
        if p * p > m:
            break
        while m % p == 0:
            f[p] = f.get(p, 0) + 1
            m //= p
    if m > 1:
        f[m] = f.get(m, 0) + 1
    return f


def sigma_from(f):
    t = 1
    for p, e in f.items():
        t *= (p ** (e + 1) - 1) // (p - 1)
    return t


def is_prime_small(n, primes):
    if n < 2:
        return False
    for p in primes:
        p = int(p)
        if p * p > n:
            return True
        if n % p == 0:
            return n == p
    return True


def main():
    out = []
    t0 = time.time()
    primes = small_primes(40000)
    oddprimes = primes[1:]
    flag, maxrat = scan(A, oddprimes, 1 << 20)
    out.append(f'scan of odd a <= {A:.0e}: {time.time() - t0:.0f} s; max sigma(a^2)/a^2 over odd a <= A: {maxrat.max():.5f}')
    bound = math.prod(p / (p - 1) for p in (3, 5, 7, 11, 13, 17, 19, 23))
    out.append(f'prod_(3<=p<=23) p/(p-1) = {bound:.5f};  3*5*7*11*13*17*19*23 = {3*5*7*11*13*17*19*23}, *29 = {3*5*7*11*13*17*19*23*29}; '
               f'bound * 1e18 = {bound * 1e18:.3e} < 2^63 = {2 ** 63:.3e}')
    idx = np.nonzero(flag)[0]
    hits = []
    for i in idx:
        a = 2 * int(i) + 1
        fa0 = factor_small(a, primes)
        hits.append((a, sigma_from({p: 2 * e for p, e in fa0.items()}) - a * a))
    out.append(f'double flips: {len(hits)}')
    rows = []
    for a, m in hits:
        pa = is_prime_small(a, primes)
        r = math.isqrt(m)
        if r * r == m:
            fy = factor_small(r, primes)
            fm = {p: 2 * e for p, e in fy.items()}
            kind = 'square'
        else:
            r = math.isqrt(m // 2)
            fy = factor_small(r, primes)
            fm = {p: 2 * e for p, e in fy.items()}
            fm[2] = fm.get(2, 0) + 1
            kind = 'twice a square'
        fa = factor_small(a, primes)
        siga2 = sigma_from({p: 2 * e for p, e in fa.items()})
        assert siga2 - a * a == m
        closes = sigma_from(fm) == siga2
        rows.append((a, m, kind, pa, fa, closes))
    out.append(f'any amicable closure sigma(m) = sigma(a^2): {[r for r in rows if r[5]]}')
    for e in range(1, 10):
        X = 10 ** e
        if X > A:
            break
        tot = sum(1 for r in rows if r[0] <= X)
        comp = sum(1 for r in rows if r[0] <= X and not r[3])
        out.append(f'   a <= 10^{e}: {tot} double flips, {comp} composite')
    comps = [r for r in rows if not r[3]]
    primes_rows = [r for r in rows if r[3]]
    out.append(f'prime double flips: {len(primes_rows)}; all of the form a = 2b^2 - 1 or a = 3: '
               f'{all(r[0] == 3 or in_twice_sq_plus(r[0]) for r in primes_rows)}; all with s(a^2) = a + 1: {all(r[1] == r[0] + 1 for r in primes_rows)}')
    # independent count of primes 2 b^2 - 1 <= A
    cnt = 1   # a = 3
    b = 2
    while 2 * b * b - 1 <= A:
        if is_prime_small(2 * b * b - 1, primes):
            cnt += 1
        b += 1
    out.append(f'independent count: primes 2b^2 - 1 <= A plus a = 3: {cnt}')
    out.append(f'composite double flips ({len(comps)}):')
    special = []
    for a, m, kind, pa, fa, closes in comps:
        sp = [p for p in fa if is_twice_sq_minus1_prime(p)]
        if sp:
            special.append(a)
        out.append(f'   a = {a} = {fa}: s(a^2) = {m} ({kind}), ratio {m / (a * a):.4f}, primes p | a with p + 1 = 2b^2: {sp}')
    out.append(f'composite double flips having a prime factor p with p + 1 = 2 b^2: {len(special)} of {len(comps)}')
    best = max(rows, key=lambda r: r[1] / (r[0] * r[0]))
    out.append(f'max ratio s(a^2)/a^2 over all double flips: {best[1] / best[0] ** 2:.4f} at a = {best[0]} = {best[4]}')
    # p = 71 conic
    sols = []
    for t in range(1, 1000):
        num = 2 * t * t - 5113
        den = 5113 - 24 * t
        if den != 0 and num * den > 0 and num % den == 0:
            sols.append((t, num // den))
    out.append(f'p = 71 twice-square conic: (t, q) with q = (2t^2 - 5113)/(5113 - 24t) a positive integer: {sols};'
               f' primes among q: {[q for t, q in sols if is_prime_small(q, primes)]}; 85625 = {factor_small(85625, primes)}')
    # direct check: primes q <= 10^6, q != 71: s(71^2 q^2) twice a square
    direct = []
    for q in primes:
        q = int(q)
        if q in (2, 71):
            continue
        sv = 72 * q * q + 5113 * (q + 1)
        assert sv == (1 + 71 + 71 ** 2) * (1 + q + q * q) - 71 ** 2 * q * q
        if sv % 2 == 0 and math.isqrt(sv // 2) ** 2 == sv // 2:
            direct.append(q)
    out.append(f'direct check over primes q < 40000: s(71^2 q^2) twice a square for q in {direct}')
    print('\n'.join(out))


def in_twice_sq_plus(a):
    v = (a + 1) // 2
    return (a + 1) % 2 == 0 and math.isqrt(v) ** 2 == v


def is_twice_sq_minus1_prime(p):
    return in_twice_sq_plus(p)


if __name__ == '__main__':
    main()
