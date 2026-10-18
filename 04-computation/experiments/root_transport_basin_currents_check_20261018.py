#!/usr/bin/env python3
"""Basin form of the frontier lemma: c_Y(u) = f(u) R_u(N), R_u(N) = sum over the first-entry basin of u (S42).

For a source m (odd, m <= 4N) let u(m) be the first seed <= Y on the odd-transport path of seed(m). Along every path the
transport sign is sigma(m -> u) = f(m) f(u) (because f(P(s)) = (-1)^a f(s) and stripping B-chains preserves f), so the
frontier coefficient is c_Y(u) = f(u) R_u(N) with the basin current
    R_u(N) = sum_{m odd <= 4N, u(m) = u} f(m) k_N(m),   k_N(m) = w(m) (m > 2N), -3 chi(m) (m < 2N).
Hence E_Y(N) = sum_u R_u(N)^2 and R(N) = sum_u R_u(N), and Lemma F is Cauchy-Schwarz over the basins. This script checks
c_Y(u) = f(u) R_u(N) seed by seed, by brute force (orbit parities computed directly), for several N and Y.
Usage: python root_transport_basin_currents_check_20261018.py
Session opus-2026-10-17-S42.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def fpar(n):
    c = 0
    while n != 1:
        n = n // 2 if n % 2 == 0 else (3 * n + 1) // 2
        c += 1
    return 1 - 2 * (c & 1)


def seed(u):
    while u % 8 == 5:
        u = (u - 1) // 4
    return u


def step(s):
    u = 3 * s + 1; a = 0
    while u % 2 == 0:
        u //= 2; a += 1
    return seed(u), (1 if a == 2 else -1)


def k(m, N):
    if m > 2 * N:
        return 2 if m % 3 == 2 else -1
    r = m % 3
    return -3 if r == 1 else (3 if r == 2 else 0)


def main():
    total = 0; bad = 0
    for N in (1, 7, 100, 1000, 7910, 30001):
        for Y in (3, 30, 100, 1000):
            c = {}; Rb = {}
            for m in range(1, 4 * N, 2):
                x = seed(m); sg = 1
                while x > Y:
                    x, g = step(x); sg *= g
                c[x] = c.get(x, 0) + sg * k(m, N)
                Rb[x] = Rb.get(x, 0) + fpar(m) * k(m, N)
            for u in set(c) | set(Rb):
                total += 1
                if c.get(u, 0) != fpar(u) * Rb.get(u, 0):
                    bad += 1
            R = sum(Rb.values()); E = sum(v * v for v in c.values())
            Rdirect = sum((2 if n % 3 == 2 else -1) * fpar(n) for n in range(N + 1, 4 * N + 1))
            assert R == Rdirect, (N, Y, R, Rdirect)
            assert E == sum(v * v for v in Rb.values())
    print('basin identity c_Y(u) = f(u) R_u(N): %d seed checks, %d violations; R = sum_u R_u and E_Y = sum_u R_u^2 checked' % (total, bad))


if __name__ == '__main__':
    main()
