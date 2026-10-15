#!/usr/bin/env python3
"""Audit (f/g) misc numerical facts quoted in the S35 notes.
  1. Collatz stopping-time sieve counts (classes mod 2^k not yet dropped within k steps) vs OEIS A076227; doublings.
  2. continued fraction of log_2 3: partial quotient 55 and the bracketing convergents.
  3. 'octahedron minus two edges at a vertex' = q + a_q also for q = 1 (mod 4)?  (the session tested q = 3 mod 4 only)
  4. sphere-packing identities and small arithmetic quoted in the atlas note.
"""
from fractions import Fraction
import mpmath
import numpy as np
import numba as nb
import sympy as sp

# 1. sieve: parity vectors of length k (T(n) = n/2 or (3n+1)/2) with 3^(o_j) > 2^j for all j <= k
counts = []
alive = [()]          # parity prefixes that have not dropped
for k in range(1, 21):
    new = []
    for w in alive:
        for b in (0, 1):
            ww = w + (b,)
            o = sum(ww)
            if 3 ** o > 2 ** k:
                new.append(ww)
    alive = new
    counts.append(len(alive))
print('1. sieve counts k = 1..20:', counts)
print('   A076227 (offset 0): 1, 1, 1, 2, 3, 4, 8, 13, 19, 38, 64, 128, 226, 367, 734, 1295, 2114, 4228, 7495, 14990')
print('   exact doublings at k =', [k + 1 for k in range(1, 20) if counts[k] == 2 * counts[k - 1]])

# 2. continued fraction of log2(3)
mpmath.mp.dps = 80
x = mpmath.log(3, 2)
cf = []
y = x
for _ in range(20):
    a = int(mpmath.floor(y))
    cf.append(a)
    y = 1 / (y - a)
print('2. log2 3 = ', cf)
h0, h1, k0, k1 = 1, cf[0], 0, 1
conv = [Fraction(cf[0], 1)]
for a in cf[1:]:
    h0, h1 = h1, a * h1 + h0
    k0, k1 = k1, a * k1 + k0
    conv.append(Fraction(h1, k1))
i55 = cf.index(55)
print(f'   partial quotient 55 at index {i55}; convergents before/after: {conv[i55 - 1]}, {conv[i55]}')


# 3. octahedron minus {35, 45} at q = 1 mod 4
def beukers(N):
    P = [0] * N
    P[0] = 1
    for m in range(1, N):
        for step in (2 * m, 4 * m):
            if step >= N:
                continue
            for _ in range(4):
                for i in range(N - 1, step - 1, -1):
                    P[i] -= P[i - step]
    return [0] + P[:N - 1]


@nb.njit(cache=True)
def U6(q, chi, eu, ev):
    tot = 0
    v = np.zeros(6, dtype=np.int64)
    v[1] = 1
    for code in range(q ** 4):
        c = code
        for i in range(4):
            v[2 + i] = c % q
            c //= q
        p = 1
        for k in range(eu.shape[0]):
            p *= chi[(v[ev[k]] - v[eu[k]]) % q]
            if p == 0:
                break
        tot += p
    return tot


a = beukers(120)
E = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 5), (2, 4), (2, 5), (3, 4)]
eu = np.array([e[0] for e in E], dtype=np.int64)
ev = np.array([e[1] for e in E], dtype=np.int64)
res = []
for q in sp.primerange(5, 80):
    chi = np.full(q, -1, dtype=np.int64)
    chi[0] = 0
    for t in range(1, q):
        chi[(t * t) % q] = 1
    val = int(U6(q, chi, eu, ev))
    res.append((q, val, q + a[q], val == q + a[q]))
print('3. octahedron - {35,45}: (q, U/(q(q-1)), q + a_q, equal) =', res)

# 4. arithmetic
print('4. 3^5 = 1 + 2*11 + 4*55:', 3 ** 5 == 1 + 2 * 11 + 4 * 55, '; 2^11 = 1 + 23 + C(23,2) + C(23,3):', 2 ** 11 == 1 + 23 + 253 + 1771,
      '; 2^7 = 11^2 + 7:', 2 ** 7 == 121 + 7, '; 2^11 - 3^7 =', 2 ** 11 - 3 ** 7, '; 2^8 - 3^5 =', 2 ** 8 - 3 ** 5,
      '; (5 - sqrt5)/10 =', (5 - 5 ** 0.5) / 10, '= 1/(1+phi^2):', 1 / (1 + ((1 + 5 ** 0.5) / 2) ** 2))
print('   Bell(8), Bell(12), Bell(14) =', sp.bell(8), sp.bell(12), sp.bell(14))
