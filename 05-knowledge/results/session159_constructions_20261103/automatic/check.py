#!/usr/bin/env python3
from collections import deque

def A(n):
    if n == 0:
        return 1
    if n == 1:
        return 2
    m, b = divmod(n, 2)
    return A(m) + m if b == 0 else A(m) + 1 - (m & 1)

def H_word(w):
    N = len(w)
    v = N + 1
    out = [0] * v
    for i in range(N):
        for j in range(i + 1, N):
            out[i] |= 1 << j
        if w[i] == '1':
            out[N] |= 1 << i
        else:
            out[i] |= 1 << N
    dp = [[0] * v for _ in range(1 << v)]
    for i in range(v):
        dp[1 << i][i] = 1
    for mask in range(1, 1 << v):
        for i in range(v):
            c = dp[mask][i]
            if c:
                avail = out[i] & ~mask
                while avail:
                    low = avail & -avail
                    j = low.bit_length() - 1
                    dp[mask | low][j] += c
                    avail -= low
    return sum(dp[-1])
for k in range(1, 13):
    M = 1 << k
    vals = [A(n) for n in range(M)]
    if not sum(vals) == M * M // 4 + M:
        raise AssertionError()
    if not sum((n + vals[n] for n in range(M))) == 3 * M * M // 4 + M // 2:
        raise AssertionError()
    if not all(((n + vals[n]) % 2 for n in range(M))):
        raise AssertionError()
for n in range(1, 128):
    w = bin(n)[2:] + '0'
    if not H_word(w) == n + A(n):
        raise AssertionError((w, n, A(n)))
for q in range(2, 501):
    r = q
    v = 0
    while r % 2 == 0:
        r //= 2
        v += 1
    j = max(1, v)
    while pow(2, j, r) != 1 % r:
        j += 1
    n, a = (1 % q, 2 % q)
    values = []
    previous_n = None
    D = None
    for t in range(q + 2):
        oldn, olda = (n, a)
        n = (pow(2, 2 * j, q) * n + pow(2, j, q) - 1) % q
        a = (a + (pow(2, j, q) - 1) * oldn + 1) % q
        if previous_n is not None:
            if not n == previous_n:
                raise AssertionError()
            delta = (a - olda) % q
            if D is None:
                D = delta
            if not D == delta:
                raise AssertionError()
            if not (D % r == 1 % r and (v == 0 or D % (1 << v) == 2 % (1 << v))):
                raise AssertionError()
            values.append((n + a) % q)
        previous_n = n
    target = {x for x in range(q) if q % 2 == 1 or x % 2 == 1}
    if not set(values) == target:
        raise AssertionError((q, j, D, sorted(target - set(values))))
print('PASS: 127 independent tournament DP counts; exact dyadic moments k<=12; modular block coverage q<=500')
N = 2048
g = [(n + A(n) - 1) // 2 & 1 for n in range(N)]
for n in range(N):
    out = 0
    for k in range(8):
        if k & ~7 == 0 and n >= k and ((n - k) % 2 == 0):
            out ^= g[(n - k) // 2]
    for k in (0, 2, 4, 6):
        if n >= k:
            out ^= g[n - k]
    if n in (1, 2, 4):
        out ^= 1
    if not out == 0:
        raise AssertionError((n, out))
from fractions import Fraction
for j in range(1, 5):
    Q = 1 << j
    n = 1
    a = 2
    for t in range(8):
        closed = 3 + Fraction((Q + 2) ** 2 * (Q ** (2 * t) - 1), (Q + 1) ** 2) + Fraction(2 * t, Q + 1)
        if not closed == n + a:
            raise AssertionError((j, t, closed, n + a))
        n, a = (Q * Q * n + Q - 1, a + (Q - 1) * n + 1)
print('PASS: 2048 characteristic-two series coefficients; closed form j<=4, t<=7')
