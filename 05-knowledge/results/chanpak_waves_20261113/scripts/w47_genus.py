"""Genera of X(N), X_0(N), X_1(N); where 2g+1 lands on the Chan–Pak holes {7, 21}."""
from fractions import Fraction
from math import gcd
def primes(n):
    ps, d = [], 2
    while d * d <= n:
        if n % d == 0:
            ps.append(d)
            while n % d == 0: n //= d
        d += 1
    if n > 1: ps.append(n)
    return ps
def phi(n):
    r = n
    for p in primes(n): r = r // p * (p - 1)
    return r
def legendre(a, p):
    a %= p
    if a == 0: return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1
def mu_full(N):          # |PSL(2, Z/N)|
    if N == 2: return 6
    m = Fraction(N ** 3, 2)
    for p in primes(N): m *= Fraction(p * p - 1, p * p)
    return int(m)
def g_full(N):
    if N <= 5: return 0
    return 1 + mu_full(N) * (N - 6) // (12 * N)
def g_0(N):
    mu = N
    for p in primes(N): mu = mu * (p + 1) // p
    nu2 = 0 if N % 4 == 0 else eval('*'.join(['1'] + [str(1 + (legendre(-1, p) if p != 2 else 0)) for p in primes(N)]))
    nu3 = 0 if N % 9 == 0 else eval('*'.join(['1'] + [str(1 + (legendre(-3, p) if p != 3 else 0)) for p in primes(N)]))
    cusps = sum(phi(gcd(d, N // d)) for d in range(1, N + 1) if N % d == 0)
    return int(1 + Fraction(mu, 12) - Fraction(nu2, 4) - Fraction(nu3, 3) - Fraction(cusps, 2))
def g_1(N):
    if N <= 4: return 0
    m = Fraction(N * N, 24)
    for p in primes(N): m *= Fraction(p * p - 1, p * p)
    s = sum(phi(d) * phi(N // d) for d in range(1, N + 1) if N % d == 0)
    return int(1 + m - Fraction(s, 4))
print("N:      ", list(range(2, 17)))
print("g(X(N)):", [g_full(N) for N in range(2, 17)])
print("2g+1:   ", [2 * g_full(N) + 1 for N in range(2, 17)])
print("X(N) with 2g+1 in {7,21}, N <= 2000:", [N for N in range(2, 2001) if 2 * g_full(N) + 1 in (7, 21)])
print("X_0(N) with 2g+1 = 7:", [N for N in range(1, 1001) if 2 * g_0(N) + 1 == 7])
print("X_0(N) with 2g+1 = 21:", [N for N in range(1, 1001) if 2 * g_0(N) + 1 == 21])
print("X_1(N) with 2g+1 in {7,21}:", [(N, 2 * g_1(N) + 1) for N in range(1, 400) if 2 * g_1(N) + 1 in (7, 21)])
print("checks: g(X(9)) = 2 g(X(8)):", g_full(9) == 2 * g_full(8), "; g(X(9)) = 3 g(X(7)) + 1:", g_full(9) == 3 * g_full(7) + 1,
      "; |PSL(2,Z/6)| =", mu_full(6), "; |PSL(2,7)| =", mu_full(7), "; |PSL(2,Z/8)| =", mu_full(8), "; |PSL(2,Z/9)| =", mu_full(9))
# sanity: known small values of g_0: g(X_0(11)) = 1, g(X_0(37)) = 2, g(X_0(N)) = 3 first at N = 30? ; g_1(13) = 2
print("sanity g0(11), g0(37), g0(23), g1(13), g1(11):", g_0(11), g_0(37), g_0(23), g_1(13), g_1(11))
