#!/usr/bin/env python3
"""Independent audit of THM-4669 (S43 referee, 2026-10-18).

Claim: for every odd prime p and r >= 1, n = (p^r - 1)/2,
    s7(n) == (-1)^n A(n)  (mod p^2),
A(n) = sum_k C(n,k)^2 C(n+k,k)^2, s7(n) = sum_k C(n,k)^2 C(n+k,k) C(2k,n).
Corollary: s7((p-1)/2) == (-1/p) a(p) (mod p^2), sum a(m) q^m = eta(2tau)^4 eta(4tau)^4.

Written from scratch (no session code imported). Parts:
  1. exact big-integer test at r = 1 for all odd primes p < PMAX1, plus mod p^3 failure census;
  2. exact test at Frobenius points r >= 2 (list below), plus mod p^3;
  3. corollary against eta(2tau)^4 eta(4tau)^4 (q-expansion by pentagonal series);
  4. the claimed failures at n = p - 1 and n = (3p-1)/2;
  5. proof steps p-adically at Frobenius points (valuations of the neglected terms);
  6. formal power series identities: Clausen 2F1(1/4,1/4;1;z)^2 = 3F2(1/2,1/2,1/2;1,1;z),
     quadratic transformation 2F1(1/4,1/4;1;4x(1-x)) = 2F1(1/2,1/2;1;x), and (E2).
Usage: python audit_thm4669.py [part ...]   (default: all)
"""
import sys, time
from math import comb
from fractions import Fraction as Fr

def primes_upto(N):
    s = bytearray([1]) * (N + 1)
    s[0:2] = b'\x00\x00'
    for i in range(2, int(N ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(N + 1) if s[i]]

def s7(n):
    return sum(comb(n, k) ** 2 * comb(n + k, k) * comb(2 * k, n) for k in range((n + 1) // 2, n + 1))

def apery(n):
    return sum((comb(n, k) * comb(n + k, k)) ** 2 for k in range(n + 1))

def vp_int(x, p):
    if x == 0:
        return 10 ** 9
    v = 0
    while x % p == 0:
        x //= p
        v += 1
    return v

def vp(x, p):
    x = Fr(x)
    if x == 0:
        return 10 ** 9
    return vp_int(x.numerator, p) - vp_int(x.denominator, p)

def part1(PMAX1=4000):
    t = time.time()
    ps = [p for p in primes_upto(PMAX1) if p > 2]
    fail2, hold3 = [], []
    for p in ps:
        n = (p - 1) // 2
        d = s7(n) - (-1) ** n * apery(n)
        if d % (p * p):
            fail2.append(p)
        if d % (p ** 3) == 0:
            hold3.append(p)
    print('[1] r = 1: odd primes 3 <= p < %d: %d tested; mod p^2 failures: %s; primes where it even holds mod p^3: %s  (%.1fs)'
          % (PMAX1, len(ps), fail2, hold3, time.time() - t), flush=True)

FROB = [(3, r) for r in range(2, 10)] + [(5, r) for r in range(2, 7)] + [(7, r) for r in range(2, 6)] + \
       [(11, 2), (11, 3), (11, 4), (13, 2), (13, 3), (13, 4), (17, 3), (19, 3), (23, 3), (29, 3), (31, 3)] + \
       [(p, 2) for p in primes_upto(150) if p > 13]

def part2():
    for p, r in FROB:
        t = time.time()
        n = (p ** r - 1) // 2
        d = s7(n) - (-1) ** n * apery(n)
        print('[2] p = %d, r = %d, n = %d: mod p^2 %s, mod p^3 %s, v_p(diff) = %d  (%.1fs)'
              % (p, r, n, d % (p * p) == 0, d % p ** 3 == 0, vp_int(d, p), time.time() - t), flush=True)

def eta_powers_series(N):
    """coefficients of prod_{n>=1} (1 - x^n) up to x^N via the pentagonal number theorem"""
    c = [0] * (N + 1)
    k = 0
    while True:
        done = True
        for kk in ((k, ) if k == 0 else (k, -k)):
            e = kk * (3 * kk - 1) // 2
            if e <= N:
                c[e] += (-1) ** (kk % 2)
                done = False
        if done and k > 0:
            break
        k += 1
    return c

def mul(a, b, N):
    import numpy as np
    A = np.array(a[:N + 1], dtype=object)
    B = np.array(b[:N + 1], dtype=object)
    out = [0] * (N + 1)
    nzA = [i for i in range(N + 1) if a[i]]
    for i in nzA:
        ai = a[i]
        for j in range(N + 1 - i):
            if b[j]:
                out[i + j] += ai * b[j]
    return out

def part3(PMAX3=3000):
    t = time.time()
    M = PMAX3 // 2 + 2
    P = eta_powers_series(M)               # prod (1 - x^n)
    P2 = [0] * (M + 1)                      # prod (1 - x^(2n))
    for i, c in enumerate(P):
        if 2 * i <= M:
            P2[2 * i] = c
    def powr(s, e):
        r = [1] + [0] * M
        for _ in range(e):
            r = mul(r, s, M)
        return r
    g = mul(powr(P, 4), powr(P2, 4), M)     # eta(2tau)^4 eta(4tau)^4 = q * g(q^2)
    def a(m):  # coefficient of q^m
        return g[(m - 1) // 2] if m % 2 == 1 else 0
    print('[3] eta(2tau)^4 eta(4tau)^4 = q + %s ...' % ' '.join('%+d q^%d' % (a(m), m) for m in range(3, 16, 2) if a(m)), flush=True)
    bad = []
    ps = [p for p in primes_upto(PMAX3) if p > 2]
    for p in ps:
        n = (p - 1) // 2
        leg = 1 if p % 4 == 1 else -1
        if (s7(n) - leg * a(p)) % (p * p):
            bad.append(p)
    print('[3] corollary s7((p-1)/2) == (-1/p) a(p) mod p^2 for odd primes p < %d: %d tested, failures %s  (%.1fs)'
          % (PMAX3, len(ps), bad, time.time() - t), flush=True)

def part4():
    ps = [p for p in primes_upto(300) if p > 3]
    for label, f in (('n = p-1', lambda p: p - 1), ('n = (3p-1)/2', lambda p: (3 * p - 1) // 2)):
        held2 = [p for p in ps if (s7(f(p)) - (-1) ** f(p) * apery(f(p))) % (p * p) == 0]
        held1 = [p for p in ps if (s7(f(p)) - (-1) ** f(p) * apery(f(p))) % p == 0]
        heldm = [p for p in ps if (s7(f(p)) + (-1) ** f(p) * apery(f(p))) % p == 0]
        print('[4] %s, primes 5..300: holds mod p^2 at %s; mod p at %s; with opposite sign mod p at %s'
              % (label, held2, held1, heldm), flush=True)
    ps = [p for p in primes_upto(300) if p > 3]
    print('[4] s7(p-1) mod p for p < 60:', [(p, s7(p - 1) % p) for p in ps if p < 60], flush=True)

def part5():
    """valuations of each neglected term in the proof at Frobenius points"""
    cases = [(3, 2), (3, 3), (3, 4), (3, 5), (5, 2), (5, 3), (7, 2), (7, 3), (11, 2), (13, 2), (17, 2), (19, 2), (23, 2)]
    for p, r in cases:
        t = time.time()
        n = (p ** r - 1) // 2
        eps = Fr(p ** r, 2)
        a = [Fr(comb(2 * k, k), 4 ** k) for k in range(n + 1)]
        C = [comb(n, k) for k in range(n + 1)]
        # (C1) ratio and difference
        minC1 = min(vp(C[k] * comb(n + k, k) - (-1) ** k * a[k] ** 2, p) for k in range(n + 1))
        # (C2): a_m - (-1)^m C(n,m)(1 - eps S_m); also the first-order term eps*S_m valuation
        S = Fr(0)
        minC2 = 10 ** 9
        minS = 10 ** 9
        Ss = [Fr(0)]
        for m in range(1, n + 1):
            S += Fr(1, n - (m - 1))
            Ss.append(S)
        for m in range(n + 1):
            minC2 = min(minC2, vp(a[m] - (-1) ** m * C[m] * (1 - eps * Ss[m]), p))
            if m:
                minS = min(minS, vp(eps * Ss[m], p))
        # (L) termwise: a_k^2 a_{n-k}^2 - a_k^4 - 2 eps C(n,k)^4 (S_k - S_{n-k}) should be in p^2 Z_p
        minLterm = min(vp(a[k] ** 2 * a[n - k] ** 2 - a[k] ** 4 - 2 * eps * C[k] ** 4 * (Ss[k] - Ss[n - k]), p)
                       for k in range(n + 1))
        anti = sum(C[k] ** 4 * (Ss[k] - Ss[n - k]) for k in range(n + 1))
        L = vp(sum(a[k] ** 2 * a[n - k] ** 2 for k in range(n + 1)) - sum(x ** 4 for x in a), p)
        st2 = vp(s7(n) - (-1) ** n * sum(a[k] ** 2 * a[n - k] ** 2 for k in range(n + 1)), p)
        print('[5] p=%d r=%d n=%d: min v_p (C1 diff) = %d; min v_p (C2 remainder) = %d; min v_p(eps S_m) = %d; '
              'min v_p (L termwise remainder) = %d; antisym sum = %s; v_p(L diff) = %d; v_p(statement 2 diff) = %d (%.1fs)'
              % (p, r, n, minC1, minC2, minS, minLterm, anti, L, st2, time.time() - t), flush=True)

def hyp_coeffs(params_num, params_den, N):
    """coefficients of pFq(num; den; z) as Fractions up to z^N"""
    c = [Fr(1)]
    for k in range(N):
        t = c[-1]
        for a in params_num:
            t *= (a + k)
        for b in params_den:
            t /= (b + k)
        t /= (k + 1)
        c.append(t)
    return c

def series_mul(a, b, N):
    out = [Fr(0)] * (N + 1)
    for i in range(N + 1):
        if a[i]:
            for j in range(N + 1 - i):
                out[i + j] += a[i] * b[j]
    return out

def compose_4x1mx(c, N):
    """sum_k c_k (4x(1-x))^k as a series in x up to x^N"""
    out = [Fr(0)] * (N + 1)
    for k in range(N + 1):
        if c[k] == 0:
            continue
        # (4x(1-x))^k = 4^k x^k sum_j C(k,j) (-x)^j
        for j in range(0, N - k + 1):
            out[k + j] += c[k] * 4 ** k * comb(k, j) * (-1) ** j
    return out

def part6(N=60):
    q = Fr(1, 4); h = Fr(1, 2)
    f14 = hyp_coeffs([q, q], [Fr(1)], N)
    f3 = hyp_coeffs([h, h, h], [Fr(1), Fr(1)], N)
    clausen = series_mul(f14, f14, N) == f3
    f12 = hyp_coeffs([h, h], [Fr(1)], N)
    quad = compose_4x1mx(f14, N) == f12
    lhs = compose_4x1mx(f3, N)
    rhs = series_mul(f12, f12, N)
    e2 = lhs == rhs
    a = [Fr(comb(2 * k, k), 4 ** k) for k in range(N + 1)]
    e2b = all(sum((-1) ** (n - k) * 4 ** k * a[k] ** 3 * comb(k, n - k) for k in range(n + 1)) == sum(a[k] ** 2 * a[n - k] ** 2 for k in range(n + 1)) for n in range(N + 1))
    e1 = all(comb(n, k) * comb(2 * k, n) == comb(2 * k, k) * comb(k, n - k) for n in range(200) for k in range(n + 1))
    print('[6] to order %d: Clausen 2F1(1/4,1/4;1;z)^2 = 3F2(1/2,1/2,1/2;1,1;z): %s; quadratic 2F1(1/4,1/4;1;4x(1-x)) = 2F1(1/2,1/2;1;x): %s; '
          '3F2(..;4x(1-x)) = 2F1(1/2,1/2;1;x)^2: %s; (E2) coefficientwise: %s; (E1) for n < 200: %s' % (N, clausen, quad, e2, e2b, e1), flush=True)

if __name__ == '__main__':
    parts = sys.argv[1:] or ['6', '4', '5', '3', '1', '2']
    for p in parts:
        {'1': part1, '2': part2, '3': part3, '4': part4, '5': part5, '6': part6}[p]()
