#!/usr/bin/env python3
r"""S37 audit, part 7: the companion note special_primes_curling_genus_zero_langlands_20261015.md, with methods
independent of the session's script where possible.

 1. Non-generous primes in (2, 10^7): numba spf sieve + least primitive root, Python pow mod p^2.  6692367337.
 2. Orbit-count corollary: brute-force orbit counts of x -> g' x on Z_(p^2) \ {0} (g' = odd-order part of a primitive
    root g mod p) against E = (p-1)/d' (g primitive mod p^2) or (p+1)(p-1)/(2d') (g^(p-1) = 1 mod p^2), for EVERY
    primitive root class mod p^2 of every prime 5 <= p <= 61 (p = 40487 itself: a7b_orbit40487.py).
 3. Woodall / Cullen itineraries (n <= 80), including the exact length of the initial odd run of W_n.
 4. Leyland primes x^y + y^x <= 10^12 (all 2 <= x <= y), and 2^y + y^2 prime for y < 400.
 5. Supersingular: F_p-rational supersingular j counted by POINT COUNTING (a_p = 0), 5 <= p < 2000, against
    floor(p/12) + eps_p and against the Paley half-imbalance I(p) (p = 3 mod 4).
 6. Stern primes < 10^6 and odd composites < 10^6 not of the form q + 2a^2 (numpy sieve).
 7. Wagstaff primes (2^p + 1)/3, p < 1300 (BPSW).
"""
import time, math
import numpy as np
import numba as nb
import sympy as sp

t00 = time.time()


def log(m):
    print(f'[{time.time() - t00:6.1f} s] {m}', flush=True)


LIM = 10 ** 7


@nb.njit(cache=False)
def spf_sieve(n):
    spf = np.zeros(n + 1, np.int32)
    for i in range(2, n + 1):
        if spf[i] == 0:
            spf[i] = i
            if i * i <= n:
                for j in range(i * i, n + 1, i):
                    if spf[j] == 0:
                        spf[j] = i
    return spf


@nb.njit(cache=False)
def powmod(b, e, m):
    r = 1
    b %= m
    while e > 0:
        if e & 1:
            r = r * b % m
        b = b * b % m
        e >>= 1
    return r


@nb.njit(cache=False)
def least_roots(spf, n):
    out_p = np.zeros(700000, np.int64)
    out_g = np.zeros(700000, np.int64)
    k = 0
    fac = np.zeros(32, np.int64)
    for p in range(3, n):
        if spf[p] != p:
            continue
        m = p - 1
        nf = 0
        while m > 1:
            q = spf[m]
            fac[nf] = q
            nf += 1
            while m % q == 0:
                m //= q
        g = 2
        while True:
            ok = True
            for i in range(nf):
                if powmod(g, (p - 1) // fac[i], p) == 1:
                    ok = False
                    break
            if ok:
                break
            g += 1
        out_p[k] = p
        out_g[k] = g
        k += 1
    return out_p[:k], out_g[:k]


spf = spf_sieve(LIM)
ps, gs = least_roots(spf, LIM)
log(f'least primitive roots for {len(ps)} primes 3 <= p < 10^7')
nongen = [(int(p), int(g)) for p, g in zip(ps, gs) if pow(int(g), int(p) - 1, int(p) * int(p)) == 1]
log(f'1. non-generous primes 2 < p < 10^7: {nongen}')
big = 6692367337
fb = sp.factorint(big - 1)
gb = 2
while not all(pow(gb, (big - 1) // q, big) != 1 for q in fb):
    gb += 1
log(f'   p = {big}: p - 1 = {fb}, least primitive root {gb}, g^(p-1) mod p^2 == 1: {pow(gb, big - 1, big * big) == 1}, '
    f'prime: {sp.isprime(big)}')
log(f'   p = 40487: p - 1 = {sp.factorint(40486)}; 5^(p-1) mod p^2 = {pow(5, 40486, 40487 ** 2)}')


@nb.njit(cache=False)
def count_orbits(m, g):
    # orbits of x -> g x mod m on {1, .., m-1}; bitmap; also checks no orbit is closed under negation
    vis = np.zeros((m >> 3) + 1, np.uint8)
    norb = 0
    selfneg = 0
    for x in range(1, m):
        if (vis[x >> 3] >> (x & 7)) & 1:
            continue
        norb += 1
        y = x
        neg = m - x
        found = False
        while True:
            vis[y >> 3] |= np.uint8(1 << (y & 7))
            if y == neg:
                found = True
            y = (y * g) % m
            if y == x:
                break
        if found:
            selfneg += 1
    return norb, selfneg


def odd_part(n):
    while n % 2 == 0:
        n //= 2
    return n


# 2a: all primitive-root classes mod p^2, primes 5..61
nchk = 0
nfail = 0
nwief = 0
for p in sp.primerange(5, 62):
    m = p * p
    dprime = odd_part(p - 1)
    pr_modp = [g for g in range(1, p) if sp.n_order(g, p) == p - 1]
    for g0 in pr_modp:
        for k in range(p):
            g = g0 + k * p
            wief = pow(g, p - 1, m) == 1
            gp = pow(g, (p - 1) // dprime, m)
            norb, selfneg = count_orbits(m, gp)
            E = norb // 2
            Ef = (p + 1) * (p - 1) // (2 * dprime) if wief else (p - 1) // dprime
            nchk += 1
            nwief += wief
            nfail += (E != Ef) or selfneg != 0 or norb % 2 != 0
log(f'2. orbit law, all primitive-root lifts mod p^2 for 5 <= p <= 61: {nchk} bases ({nwief} with g^(p-1) = 1 mod p^2), '
    f'{nfail} failures')
p = 40487
q = 6692367337
dq = odd_part(q - 1)
log(f'   6692367337: (p-1)/d\' = {(q - 1) // dq}, non-generous E = (p+1)(p-1)/(2d\') = {(q + 1) * (q - 1) // (2 * dq)}')


# 3. Woodall / Cullen
def T(x):
    return x // 2 if x % 2 == 0 else (3 * x + 1) // 2


okW = okC = okC2 = True
runs = []
for n in range(1, 81):
    x = n * 2 ** n - 1
    for _ in range(n):
        okW &= x % 2 == 1
        x = T(x)
    okW &= x == n * 3 ** n - 1
    y = n * 2 ** n - 1
    run = 0
    while y % 2 == 1 and y > 1:
        y = T(y)
        run += 1
    v2 = (n & -n).bit_length() - 1
    runs.append((n, run, n + v2))
    c = n * 2 ** n + 1
    for j in range(0, n // 2 + 1):
        if 2 * j <= n:
            okC2 &= (c == 3 ** j * n * 2 ** (n - 2 * j) + 1)
        if 2 * j + 2 <= n:
            okC &= c % 2 == 1
            c = T(c)
            okC &= c % 2 == 0
            c = T(c)
log(f'3. Woodall: first n steps odd and T^n(W_n) = n 3^n - 1 (n <= 80): {okW}; initial odd run length = n + v2(n) '
    f'for all 2 <= n <= 80: {all(r == e for k, r, e in runs if k > 1)}; examples (n, run): {[(n, r) for n, r, e in runs[:8]]}')
log(f'   Cullen: (OE)^j reaches 3^j n 2^(n-2j) + 1 for all 2j <= n (n <= 80): {okC and okC2}')

# 4. Leyland
ley = []
x = 2
while 2 ** x + x * x <= 10 ** 12 or x <= 3:
    y = x
    while True:
        v = x ** y + y ** x
        if v > 10 ** 12:
            break
        if sp.isprime(v):
            ley.append((v, x, y))
        y += 1
    x += 1
    if x ** x * 2 > 10 ** 12:
        break
ley.sort()
log(f'4. Leyland primes x^y + y^x <= 10^12 (2 <= x <= y): {ley}')
ys = [y for y in range(1, 400) if sp.isprime(2 ** y + y * y)]
log(f'   2^y + y^2 prime, y < 400: {ys}; all y = 1 or 3 mod 6: {all(y == 1 or y % 6 == 3 for y in ys)}')


# 5. supersingular by point counting
def ss_rational_count(p):
    xs = np.arange(p, dtype=np.int64)
    chi = -np.ones(p, np.int64)
    chi[(xs[1:] ** 2) % p] = 1
    chi[0] = 0
    js = np.arange(p, dtype=np.int64)
    j1728 = 1728 % p
    A = (3 * js % p) * ((1728 - js) % p) % p
    B = (2 * js % p) * (((1728 - js) % p) ** 2 % p) % p
    A[0], B[0] = 0, 1
    A[j1728], B[j1728] = 1, 0
    x3 = xs ** 3 % p
    cnt = 0
    for lo in range(0, p, 256):
        hi = min(p, lo + 256)
        V = (x3[None, :] + (A[lo:hi, None] * xs[None, :]) % p + B[lo:hi, None]) % p
        ap = -chi[V].sum(axis=1)
        cnt += int(np.count_nonzero(ap == 0))
    return cnt


def nss(p):
    return p // 12 + {1: 0, 5: 1, 7: 1, 11: 2}[p % 12]


allrat = []
paley_ok = True
sat = []
for p in sp.primerange(5, 2000):
    c = ss_rational_count(p)
    tot = nss(p)
    if c == tot:
        allrat.append(int(p))
    if p % 4 == 3:
        I = sum(1 if pow(a, (p - 1) // 2, p) == 1 else -1 for a in range(1, (p - 1) // 2 + 1))
        pred = I if p % 8 == 7 else sp.Rational(2, 3) * I
        paley_ok &= (c == pred)
        if c == tot:
            sat.append(int(p))
log(f'5. primes 5 <= p < 2000 with all supersingular j F_p-rational (point counting): {allrat}')
log(f'   #F_p-rational ss j = I(p) (p = 7 mod 8) or (2/3) I(p) (p = 3 mod 8), all p = 3 mod 4 < 2000: {paley_ok}; '
    f'saturating p = 3 mod 4: {sat}')

# 6. Stern
M = 10 ** 6
isp = np.zeros(M + 1, bool)
isp[2:] = True
for i in range(2, int(M ** 0.5) + 1):
    if isp[i]:
        isp[i * i::i] = False
rep = np.zeros(M + 1, bool)
a = 1
while 2 * a * a <= M:
    s = 2 * a * a
    rep[s:] |= isp[:M + 1 - s]
    a += 1
stern = [int(n) for n in np.nonzero(isp & ~rep)[0]]
oddcomp = [int(n) for n in range(9, M + 1, 2) if not isp[n] and not rep[n]]
log(f'6. Stern primes < 10^6: {stern}; odd composites < 10^6 not q + 2a^2: {oddcomp}')

# 7. Wagstaff
wag = [int(p) for p in sp.primerange(3, 1300) if sp.isprime((2 ** p + 1) // 3)]
log(f'7. Wagstaff primes (2^p+1)/3, prime p < 1300 (BPSW): {wag}')
log(f'   first odd ancestors (2^m+1)/3 of 1 under 3x-1 with odd m that are NOT Wagstaff numbers (m composite), m < 30: '
    f'{[(m, (2 ** m + 1) // 3) for m in range(1, 30, 2) if not sp.isprime(m)]}')
