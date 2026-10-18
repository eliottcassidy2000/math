#!/usr/bin/env python3
"""Audit F3, item (f) support: check that the source weights used in f3_twins.py reproduce the ROOT current
R(N) = H(4N) - H(N), H(N) = sum_(n <= N) w_n t(n), w_n = 3[n = 2 mod 3] - 1, t(n) = (-1)^(T-steps to 1),
i.e. R(N) = sum over seeds x <= 4N of c_0(x) t(x), for every N <= NMAX (brute-force Collatz parities).
Usage: python3 f3_root_invariant.py NMAX"""
import sys
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
M = 4 * NMAX + 10
tau = {1: 0}
def T(n):
    return n // 2 if n % 2 == 0 else (3 * n + 1) // 2
def steps(n):
    path = []
    while n not in tau:
        path.append(n); n = T(n)
    s = tau[n]
    for m in reversed(path):
        s += 1; tau[m] = s
    return tau[path[0]] if path else s
t = [0] * (M + 1)
for n in range(1, M + 1):
    t[n] = -1 if steps(n) % 2 else 1
w = [0] + [3 * (n % 3 == 2) - 1 for n in range(1, M + 1)]
H = [0] * (M + 1)
for n in range(1, M + 1):
    H[n] = H[n - 1] + w[n] * t[n]
def seed(m):
    while m % 8 == 5:
        m = (m - 1) // 4
    return m
bad = 0
for N in range(1, NMAX + 1):
    c0 = {}
    for m in range(1, 4 * N + 1, 2):
        r = m % 3
        K = (2 if r == 2 else -1) + ((1, -2, 1)[r] if m < 2 * N else 0)
        x = seed(m)
        c0[x] = c0.get(x, 0) + K
    R = sum(c * t[x] for x, c in c0.items())
    if R != H[4 * N] - H[N]:
        bad += 1
print('N = 1..%d: R(N) = sum_x c_0(x) t(x) equals H(4N) - H(N) at every N: %s (failures %d)' % (NMAX, bad == 0, bad))
