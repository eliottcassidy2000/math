"""Independent Held-Karp, exact rational partition, and high-precision checks."""
from fractions import Fraction
from itertools import permutations
from math import prod
import mpmath as mp

def path_count(out):
    n = len(out)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for end in range(n):
            if not mask & 1 << end:
                continue
            q = dp[mask][end]
            for nxt in range(n):
                if not mask & 1 << nxt and out[end] & 1 << nxt:
                    dp[mask | 1 << nxt][nxt] += q
    return dp[-1]

def cyclic_wrap(n):
    out = [0] * (n + 2)
    b, c = (n, n + 1)
    for i in range(n):
        out[i] |= 1 << b
        out[c] |= 1 << i
        for j in range(i + 1, n):
            out[i] |= 1 << j
    out[b] |= 1 << c
    return out
for n in range(1, 7):
    out = cyclic_wrap(n)
    if not sum(path_count(out)) == 2 ** n + 1:
        raise AssertionError()
    w = [Fraction((-1) ** i * (i + 1), 2 * i + 7) for i in range(n)]
    seen = 0
    partition = Fraction()
    for p in permutations(range(n + 2)):
        if all((out[a] >> b & 1 for a, b in zip(p, p[1:]))):
            seen += 1
            partition += prod((w[a] for a in p[:p.index(n)] if a < n))
    if not partition == prod((1 + x for x in w)) + prod(w):
        raise AssertionError()
    print('wrap n=', n, 'paths=', seen, 'partition=', partition)
flips = ((1, 2), (2, 3), (5, 6))
expected = [189, 159, 159, 159, 159, 141, 141, 145]
for state in range(8):
    out = [sum((1 << j for j in range(7) if (j - i) % 7 in (1, 2, 4))) for i in range(7)]
    for k, (i, j) in enumerate(flips):
        if state >> k & 1:
            out[i] ^= 1 << j
            out[j] ^= 1 << i
    if not all((bool(out[i] >> j & 1) != bool(out[j] >> i & 1) for i in range(7) for j in range(i + 1, 7))):
        raise AssertionError()
    ends = path_count(out)
    total = sum(ends)
    x, y, z = (state >> i & 1 for i in range(3))
    poly = 189 - 30 * (x + y + z) + 30 * x * y + 12 * x * z + 12 * y * z - 8 * x * y * z
    if not total == expected[state] == poly:
        raise AssertionError()
    if state == 7:
        if not out == [22, 40, 82, 53, 98, 5, 43]:
            raise AssertionError()
        if not ends == [21, 24, 20, 17, 18, 28, 17]:
            raise AssertionError()
        cyc = (0, 1, 3, 2, 4, 6, 5, 0)
        if not all((out[a] >> b & 1 for a, b in zip(cyc, cyc[1:]))):
            raise AssertionError()
    print('flip state=', state, 'total=', total, 'endpoints=', ends)
mp.mp.dps = 50
for N in [4, 10, 100, 1000]:
    w = [(-1) ** i * mp.expm1(mp.log(i) / i) for i in range(2, N + 1)]
    C = sum(w, mp.mpf(0))
    prod1 = mp.fprod((1 + v for v in w))
    exc = mp.fprod(w)
    Z = prod1 + exc
    correction = sum((v - mp.log1p(v) for v in w), mp.mpf(0))
    error = abs(C - (mp.log(Z - exc) + correction))
    if not error < mp.mpf('1e-47'):
        raise AssertionError()
    print('N=', N, 'identity error=', mp.nstr(error, 3), 'product=', mp.nstr(prod1, 12), 'exceptional=', mp.nstr(exc, 4))
