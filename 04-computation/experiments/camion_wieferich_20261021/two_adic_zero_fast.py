# Fast 2-adic analysis of f(n) = s(n)/n! = sum over strong classes of 1/|Aut| (in Z_2).
# 1/T(x) mod 2^K equals 1/Q(x) for the polynomial Q = T truncated at the degree where v2(2^C(n,2)/n!) >= K, so
# f(n) = -[x^n] 1/Q(x) (n >= 1) satisfies an order-d linear recurrence; f(n) mod 2^K for huge n by companion-matrix powers.
# Determines the 2-adic zero n* of the class n = 2 (mod 3) bit by bit (law: v2(f(n)) = 2 + v2(n - n*)), then tests whether
# n* is a small rational (rational reconstruction).  Also checks the quasi-period f(n + 3*2^k) = f(n) mod 2^(k+c).
# mac-mini-2026-10-21-fanopsp.
import sys
from math import comb
K = int(sys.argv[1]) if len(sys.argv) > 1 else 72
MOD = 1 << K
def v2(x): return (x & -x).bit_length() - 1 if x else 10**9
def tau(n):                                     # 2^C(n,2)/n! as an element of Z/2^K
    f = 1; e = n * (n - 1) // 2
    for i in range(2, n + 1):
        while i % 2 == 0: i //= 2; e -= 1
        f = f * i % MOD
    return (pow(f, -1, MOD) << e) % MOD if e < K else 0
d = 1
while tau(d + 1) != 0 or any(tau(j) for j in range(d + 2, d + 6)): d += 1
Q = [tau(j) for j in range(d + 1)]              # Q[0] = 1
# g = 1/Q:  g_n = -sum_{j=1}^d Q[j] g_{n-j},  g_0 = 1.  f(n) = -g_n for n >= 1.
def g_small(N):
    g = [1] + [0] * N
    for n in range(1, N + 1): g[n] = -sum(Q[j] * g[n - j] for j in range(1, min(d, n) + 1)) % MOD
    return g
def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(d)) % MOD for j in range(d)] for i in range(d)]
C = [[0] * d for _ in range(d)]
for j in range(d): C[0][j] = -Q[j + 1] % MOD
for i in range(1, d): C[i][i - 1] = 1
G0 = g_small(d)                                  # state vector (g_{d-1}, ..., g_0) at n = d-1
def g_at(n):
    if n < d: return G0[n]
    P = None; B = C; e = n - (d - 1)
    while e:
        if e & 1: P = B if P is None else matmul(B, P)
        B = matmul(B, B); e >>= 1
    st = [G0[d - 1 - i] for i in range(d)]
    return sum(P[0][j] * st[j] for j in range(d)) % MOD
f = lambda n: (-g_at(n)) % MOD
print(f"K = {K} bits; recurrence order d = {d}")
# sanity: compare with exact small values s(n)/n! from Wright's recurrence
s = [0, 1]; 
for n in range(2, 30):
    s.append(2 ** comb(n, 2) - sum(comb(n, k) * s[k] * 2 ** comb(n - k, 2) for k in range(1, n)))
from fractions import Fraction
ok = all((Fraction(s[n], 1) / __import__('math').factorial(n)).numerator * pow((Fraction(s[n], 1) / __import__('math').factorial(n)).denominator, -1, MOD) % MOD == f(n) for n in range(1, 30))
print("matches exact s(n)/n! for n < 30:", ok)
# zero of the class n = 2 mod 3, bit by bit, testing at large n
BIG = 10**6
r = 2; k = 0     # n* = r mod 2^k
while k < K - 8:
    hits = []
    for c in (r, r + (1 << k)):
        # n = c mod 2^(k+1), n = 2 mod 3, n large
        n = c + (1 << (k + 1)) * (((2 - c) * pow(1 << (k + 1), -1, 3)) % 3) + 3 * (1 << (k + 1)) * BIG
        if v2(f(n)) >= 3 + k: hits.append(c)
    if len(hits) != 1: print("ambiguous/failed at k =", k, hits); break
    r = hits[0]; k += 1
print(f"n* = {r} mod 2^{k}")
# rational reconstruction: find a/b = r mod 2^k with small |a|, b odd
def ratrec(r, m):
    a0, a1, b0, b1 = m, r, 0, 1
    import math
    bound = math.isqrt(m // 2)
    while a1 > bound:
        q = a0 // a1; a0, a1 = a1, a0 - q * a1; b0, b1 = b1, b0 - q * b1
    return a1, b1
a, b = ratrec(r, 1 << k)
print("rational reconstruction:", a, "/", b, " check:", (a * pow(b, -1, 1 << k) - r) % (1 << k) == 0 if b % 2 else "b even")
# quasi-period: v2(f(n + 3*2^j) - f(n)) for n in each class, large n
for j in range(1, 9):
    vals = [v2((f(n + 3 * 2 ** j) - f(n)) % MOD) for n in range(10**5, 10**5 + 6)]
    print(f"  shift 3*2^{j}: v2(f(n+3*2^{j}) - f(n)) for 6 consecutive n: {vals}")
# thresholds: least n0 with v2(f(n + 3*2^j) - f(n)) >= j + 2 for all n0 <= n <= 600 (exact small values)
gs = g_small(600 + 3 * 2**7)
fs = [(-x) % MOD for x in gs]
for j in range(1, 8):
    P = 3 * 2**j
    fails = [n for n in range(1, 601) if v2((fs[n + P] - fs[n]) % MOD) < j + 2]
    print(f"  quasi-period 3*2^{j} mod 2^{j+2}: failing n (1..600) = {fails[:12]}{' ...' if len(fails) > 12 else ''} (count {len(fails)})")
