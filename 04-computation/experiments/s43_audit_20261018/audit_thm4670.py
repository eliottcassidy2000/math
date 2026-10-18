#!/usr/bin/env python3
"""Independent audit of THM-4670 (opus S43 version; S43 referee, 2026-10-18).

Claim: for every prime p >= 5, n = (p-1)/2:  sum_{k<=n} C(n+k,k)^3 == b(n) = sum_k C(n,k)^2 C(n+k,k)  (mod p^2).
Proof route: C(n+k,k) == a_k(1+pO_k), C(n,k) == (-1)^k a_k (1-pO_k) (mod p^2), reduction to E_n = S_n - H_n D_n == 0 (mod p),
derivative of Dixon's well-poised 3F2(1) at a = -n + d, b = c = -n.

Parts (all code written from scratch):
  1. theorem for primes 5 <= p < PMAX (exact integers), mod p^3 census, p = 3, and Frobenius points r = 2 (informational);
  2. Step 1 congruences and E_n == 0 (mod p) directly; sum a_k^3 O_k == E_n/2 (mod p);
  3. closed forms: odd n via sympy exact gamma values; even n formula; Dixon value;
  4. high-precision numerical check that Dixon holds identically in d near 0 and of the derivative bookkeeping (mpmath);
  5. the two classical guises of Beukers' zeta(2) supercongruence (Mortenson 2008 eqs (1.2), (1.3)):
       b((p-1)/2) == c(p) and sum_{k<=(p-1)/2} a_k^3 == c(p) (mod p^2), c(p) the coefficients of eta(4tau)^6,
     and the elementary identity  sum C(n+k,k)^3 - b(n) == 4 (sum a_k^3 - b(n))  (mod p^2), which makes THM-4670 an
     immediate corollary of those classical results.
Usage: python audit_thm4670.py [part ...]
"""
import sys, time
from math import comb, factorial
from fractions import Fraction as Fr

def primes_upto(N):
    s = bytearray([1]) * (N + 1)
    s[0:2] = b'\x00\x00'
    for i in range(2, int(N ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(N + 1) if s[i]]

def twin(n):
    return sum(comb(n + k, k) ** 3 for k in range(n + 1))

def bz2(n):
    return sum(comb(n, k) ** 2 * comb(n + k, k) for k in range(n + 1))

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

def red(x, m, p):
    x = Fr(x)
    assert x.denominator % p
    return x.numerator * pow(x.denominator, -1, m) % m

def H(m):
    return sum(Fr(1, j) for j in range(1, m + 1))

def part1(PMAX=5000):
    t = time.time()
    ps = [p for p in primes_upto(PMAX) if p >= 5]
    fail, hold3 = [], []
    for p in ps:
        n = (p - 1) // 2
        d = twin(n) - bz2(n)
        if d % (p * p):
            fail.append(p)
        if d % p ** 3 == 0:
            hold3.append(p)
    print('[1] primes 5 <= p < %d: %d tested; mod p^2 failures %s; holds even mod p^3 at %s (%.1fs)' % (PMAX, len(ps), fail, hold3, time.time() - t), flush=True)
    d3 = twin(1) - bz2(1)
    print('[1] p = 3 (n = 1): twin %d, b %d, difference %d, mod 9 = %d' % (twin(1), bz2(1), d3, d3 % 9), flush=True)
    for p in (3, 5, 7, 11, 13, 17, 19, 23):
        n = (p * p - 1) // 2
        d = twin(n) - bz2(n)
        print('[1] Frobenius point r = 2 (not claimed): p = %d, n = %d: v_p(twin - b) = %d' % (p, n, vp_int(d, p)), flush=True)

def part2(PMAX=1500):
    t = time.time()
    ok1 = okE = okO = True
    ps = [p for p in primes_upto(PMAX) if p >= 5]
    for p in ps:
        n = (p - 1) // 2
        m2 = p * p
        a = [Fr(comb(2 * k, k), 4 ** k) for k in range(n + 1)]
        O = [Fr(0)]
        for k in range(1, n + 1):
            O.append(O[-1] + Fr(1, 2 * k - 1))
        if p < 400:
            for k in range(n + 1):
                if (comb(n + k, k) - red(a[k] * (1 + p * O[k]), m2, p)) % m2 or \
                   (comb(n, k) - red((-1) ** k * a[k] * (1 - p * O[k]), m2, p)) % m2:
                    ok1 = False
        D = sum((-1) ** k * comb(n, k) ** 3 for k in range(n + 1))
        Hs = [Fr(0)]
        for j in range(1, n + 1):
            Hs.append(Hs[-1] + Fr(1, j))
        S = sum((-1) ** k * comb(n, k) ** 3 * Hs[n - k] for k in range(n + 1))
        E = S - Hs[n] * D
        if vp(E, p) < 1:
            okE = False
        if p < 400:
            lhs = red(sum(a[k] ** 3 * O[k] for k in range(n + 1)), p, p)
            if (lhs - red(E / 2, p, p)) % p:
                okO = False
    print('[2] Step 1 congruences (p < 400): %s; E_n == 0 mod p for all primes 5 <= p < %d: %s; sum a_k^3 O_k == E_n/2 mod p (p < 400): %s (%.1fs)'
          % (ok1, PMAX, okE, okO, time.time() - t), flush=True)

def part3(NMAX=61):
    import sympy as sp
    okodd = okeven = True
    for n in range(1, NMAX + 1):
        D = sum((-1) ** k * comb(n, k) ** 3 for k in range(n + 1))
        S = sum((-1) ** k * comb(n, k) ** 3 * H(n - k) for k in range(n + 1))
        if n % 2:
            g = sp.gamma(1 - sp.Rational(n, 2)) * sp.gamma(1 + sp.Rational(3 * n, 2)) / (3 * n * sp.gamma(1 + sp.Rational(n, 2)) ** 2)
            g = sp.nsimplify(sp.simplify(g))
            if D != 0 or sp.Rational(S.numerator, S.denominator) != g:
                okodd = False
                print('   odd mismatch at n =', n, S, g)
        else:
            m = n // 2
            if D != (-1) ** m * Fr(factorial(3 * m), factorial(m) ** 3):
                okeven = False
            if S - H(n) * D != Fr(D, 2) * (H(m) - H(3 * m)):
                okeven = False
                print('   even mismatch at n =', n)
    print('[3] odd n <= %d: D_n = 0 and S_n = Gamma(1-n/2)Gamma(1+3n/2)/(3n Gamma(1+n/2)^2) (sympy exact): %s; even n: Dixon value and E_n = (D_n/2)(H_m - H_3m): %s'
          % (NMAX, okodd, okeven), flush=True)
    # valuation claim for odd n: v_p(S_n) = 1 exactly when p = 2n+1 prime
    vals = []
    for n in range(1, 400, 2):
        p = 2 * n + 1
        if all(p % d for d in range(2, int(p ** 0.5) + 1)) and p > 3:
            S = sum((-1) ** k * comb(n, k) ** 3 * H(n - k) for k in range(n + 1)) if n < 60 else None
            if S is not None:
                vals.append((p, vp(S, p)))
    print('[3] v_p(S_n) for odd n, p = 2n+1 prime (n < 60):', vals, flush=True)

def part4():
    import mpmath as mp
    mp.mp.dps = 80
    def lhs(n, d):
        a = -n + d; b = c = -n
        s = mp.mpf(0); term = mp.mpf(1)
        for k in range(n + 1):
            s += term
            term = term * (a + k) * (b + k) * (c + k) / ((1 + a - b + k) * (1 + a - c + k) * (k + 1))
        return s
    def rhs(n, d):
        a = -n + d; b = c = -n
        G = mp.gamma; rg = mp.rgamma
        return G(1 + a / 2) * G(1 + a - b) * G(1 + a - c) * G(1 + a / 2 - b - c) * rg(1 + a) * rg(1 + a / 2 - b) * rg(1 + a / 2 - c) * rg(1 + a - b - c)
    worst = mp.mpf(0)
    for n in range(1, 25):
        for d in (mp.mpf('1e-3'), mp.mpf('-2.5e-4'), mp.mpf('0.37'), mp.mpf('1e-12')):
            r = rhs(n, d)
            diff = abs(lhs(n, d) - r) / max(1, abs(r))
            worst = max(worst, diff)
    print('[4] Dixon identity at a = -n + d, b = c = -n, n = 1..24, d in {1e-3, -2.5e-4, 0.37, 1e-12}: max relative error %s' % mp.nstr(worst, 5), flush=True)
    worstd = mp.mpf(0)
    for n in range(1, 25):
        der = mp.diff(lambda d: lhs(n, d), 0)
        D = sum((-1) ** k * comb(n, k) ** 3 for k in range(n + 1))
        Hn = [sum(mp.mpf(1) / j for j in range(1, m + 1)) for m in range(n + 1)]
        S = sum((-1) ** k * comb(n, k) ** 3 * Hn[n - k] for k in range(n + 1))
        T = sum((-1) ** k * comb(n, k) ** 3 * Hn[k] for k in range(n + 1))
        pred = -Hn[n] * D + S - 2 * T
        derR = mp.diff(lambda d: rhs(n, d), mp.mpf('1e-30')) if n % 2 == 0 else mp.diff(lambda d: rhs(n, d), 0)
        scale = max(1, abs(pred))
        worstd = max(worstd, abs(der - pred) / scale, abs(der - derR) / scale)
    print('[4] derivative bookkeeping LHS\'(0) = -H_n D_n + S_n - 2 T_n and = RHS\'(0), n <= 24: max relative error %s' % mp.nstr(worstd, 5), flush=True)

def eta_prod_coeffs(N):
    """prod (1 - x^n) to x^N (pentagonal)"""
    c = [0] * (N + 1)
    k = 0
    while True:
        hit = False
        for kk in ((0,) if k == 0 else (k, -k)):
            e = kk * (3 * kk - 1) // 2
            if e <= N:
                c[e] += -1 if kk % 2 else 1
                hit = True
        if not hit:
            break
        k += 1
    return c

def part5(PMAX=3000):
    t = time.time()
    M = PMAX // 4 + 2
    P = eta_prod_coeffs(M)
    P6 = [1] + [0] * M
    for _ in range(6):
        new = [0] * (M + 1)
        for i, x in enumerate(P6):
            if x:
                for j in range(M + 1 - i):
                    if P[j]:
                        new[i + j] += x * P[j]
        P6 = new
    def c(m):  # eta(4tau)^6 = q prod (1-q^(4n))^6
        return P6[(m - 1) // 4] if m % 4 == 1 else 0
    ps = [p for p in primes_upto(PMAX) if p >= 5]
    okcm = all((c(p) - (4 * min(x for x in range(1, int(p ** 0.5) + 1, 2) if any(x * x + y * y == p for y in range(0, int(p ** 0.5) + 1))) ** 2 - 2 * p)) == 0 for p in ps if p % 4 == 1 and p < 500)
    bad12, bad13, badtw, badid = [], [], [], []
    for p in ps:
        n = (p - 1) // 2
        m2 = p * p
        b = bz2(n) % m2
        h3 = sum(comb(2 * k, k) ** 3 * pow(64, -k, m2) for k in range(n + 1)) % m2
        tw = twin(n) % m2
        if (b - c(p)) % m2:
            bad12.append(p)
        if (h3 - c(p)) % m2:
            bad13.append(p)
        if (tw - c(p)) % m2:
            badtw.append(p)
        if (tw - b - 4 * (h3 - b)) % m2:
            badid.append(p)
    print('[5] eta(4tau)^6 = q prod(1-q^4n)^6: a(p) = 4x^2 - 2p (p = x^2 + y^2, x odd) for p = 1 mod 4 < 500: %s' % okcm, flush=True)
    print('[5] primes 5 <= p < %d: (1.2) b((p-1)/2) == c(p) mod p^2 failures %s; (1.3) sum_{k<=(p-1)/2} a_k^3 == c(p) failures %s; '
          'twin == c(p) failures %s; identity twin - b == 4 (sum a_k^3 - b) mod p^2 failures %s (%.1fs)'
          % (PMAX, bad12, bad13, badtw, badid, time.time() - t), flush=True)

if __name__ == '__main__':
    parts = sys.argv[1:] or ['3', '4', '2', '5', '1']
    for q in parts:
        {'1': part1, '2': part2, '3': part3, '4': part4, '5': part5}[q]()
