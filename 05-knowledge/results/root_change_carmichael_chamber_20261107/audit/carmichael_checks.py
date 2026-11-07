"""Auditor's checks of the elementary totient statements in CARMICHAEL_MAP.md (audit of kind-pasteur-2026-11-07-S14).
Inverse totient by recursion over primes p with (p - 1) | m; standard library only.
usage: python3 carmichael_checks.py
"""
from functools import lru_cache


def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    f = 3
    while f * f <= n:
        if n % f == 0:
            return False
        f += 2
    return True


def divisors(m):
    small, large = [], []
    d = 1
    while d * d <= m:
        if m % d == 0:
            small.append(d)
            if d * d != m:
                large.append(m // d)
        d += 1
    return small + large[::-1]


def inv_phi(m):
    """All x with phi(x) = m."""
    P = sorted(d + 1 for d in divisors(m) if is_prime(d + 1))

    @lru_cache(maxsize=None)
    def sols(mm, i):
        res = [1] if mm == 1 else []
        for j in range(i, len(P)):
            p = P[j]
            if mm % (p - 1):
                continue
            r = mm // (p - 1)
            pk = p
            while True:
                res.extend(y * pk for y in sols(r, j + 1))
                if r % p:
                    break
                r //= p
                pk *= p
        return tuple(res)
    return sorted(sols(m, 0))


def phi(n):
    r, m, p = n, n, 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            r -= r // p
        p += 1
    if m > 1:
        r -= r // m
    return r


def main():
    # self-test of inv_phi against brute force
    from collections import defaultdict
    brute = defaultdict(list)
    for x in range(1, 20001):
        brute[phi(x)].append(x)
    for m in range(1, 2001):
        got = inv_phi(m)
        exp = [x for x in brute.get(m, []) if x <= 20000]
        if m <= 2000 and got != exp and max(got, default=0) <= 20000:
            raise RuntimeError('inv_phi mismatch at %d' % m)
    print('inv_phi self-test against brute force (m <= 2000): OK')
    nont = [m for m in range(2, 301, 2) if not inv_phi(m)]
    print('even nontotients <= 300 (A005277 start):', nont)
    ev = [m for m in range(14, 167, 2)]
    nn = [m for m in ev if not inv_phi(m)]
    print('even numbers in [14,166]: %d, nontotients %d (%.1f%%)' % (len(ev), len(nn), 100.0 * len(nn) / len(ev)))
    doubles = [14, 26, 42, 82, 86, 110, 134, 154, 166]
    print('doubles of the first nine F_<=2 misses, nontotient?', [(m, not inv_phi(m)) for m in doubles])
    print('2*7^k nontotient for k = 1..9:', all(not inv_phi(2 * 7 ** k) for k in range(1, 10)),
          '; k = 0 (m = 2) totient:', bool(inv_phi(2)))
    bad = []
    for p in range(5, 20000):
        if is_prime(p):
            if (not inv_phi(2 * p)) != (not is_prime(2 * p + 1)):
                bad.append(p)
    print('for primes 3 < p < 20000: 2p nontotient iff 2p+1 composite; exceptions:', bad)
    print('A_sigma(3) = #{n : sigma(n) = 3} =', sum(1 for n in range(1, 10) if sum(d for d in range(1, n + 1) if n % d == 0) == 3))


if __name__ == '__main__':
    main()
