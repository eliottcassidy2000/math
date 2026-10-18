#!/usr/bin/env python3
"""S43b referee: the CM harmonic law (L_p) and the class-6 congruence (K), independent implementation.

(L_p)  6 sum_{k<=n} C(2k,k)^3 m^-k (H_2k - H_k) == q_p(m) sum_{k<=n} C(2k,k)^3 m^-k   (mod p),  n = (p-1)/2.
(K)    sum_{j<=n/2} a_j^3 a_2j (H_2j - 2H_j - 2H_4j + 2H_(n-j) + H_(n-2j)) == 0 (mod p),  a_k = C(2k,k)/4^k.

Method (different from the session's): discrete logarithms to a primitive root g mod p; all residues, harmonic numbers
and central binomials mod p are produced by vectorized numpy operations (exact, int64, no floating point).
Also: the class-3, class-6 and class-7 twins checked directly mod p^2 (exact integers), including the identity
    u6a - (-1)^n u6b == p * sum_j a_j^3 a_2j K_j   (mod p^2)
(the referee's own derivation of the class-6 reduction).
Usage: python audit_cm_law.py [PMAX_ALL] [PMAX_GEN] [PMAX_TWIN]
"""
import sys
import time
from math import comb
import numpy as np

CM = [1, -8, 16, -64, 64, 256, -512, 4096]


def primes_in(a, b):
    s = np.ones(b + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(b ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return [int(x) for x in np.nonzero(s)[0] if x >= a]


def prim_root(p):
    m = p - 1
    fs = []
    d = 2
    while d * d <= m:
        if m % d == 0:
            fs.append(d)
            while m % d == 0:
                m //= d
        d += 1
    if m > 1:
        fs.append(m)
    for g in range(2, p):
        if all(pow(g, (p - 1) // f, p) != 1 for f in fs):
            return g


def dlog_tables(p):
    g = prim_root(p)
    B = int((p - 1) ** 0.5) + 1
    small = np.empty(B, dtype=np.int64)
    x = 1
    for i in range(B):
        small[i] = x
        x = x * g % p
    gB = x
    big = np.empty(B, dtype=np.int64)
    y = 1
    for i in range(B):
        big[i] = y
        y = y * gB % p
    pw = ((big[:, None] * small[None, :]) % p).ravel()[:p - 1]      # pw[e] = g^e
    dl = np.zeros(p, dtype=np.int64)
    dl[pw] = np.arange(p - 1, dtype=np.int64)
    return pw, dl


def fermat_q(m, p):
    pp = p * p
    return ((pow(m % pp, p - 1, pp) - 1) // p) % p


def prime_data(p):
    pw, dl = dlog_tables(p)
    P1 = p - 1
    n = (p - 1) // 2
    idx = np.arange(1, p, dtype=np.int64)
    inv = np.zeros(p, dtype=np.int64)
    inv[1:] = pw[(-dl[idx]) % P1]
    H = np.zeros(p, dtype=np.int64)
    H[1:] = np.cumsum(inv[1:]) % p
    lf = np.zeros(p, dtype=np.int64)
    lf[1:] = np.cumsum(dl[1:]) % P1                        # log of i! for i <= p-1
    k = np.arange(n + 1, dtype=np.int64)
    Lc = (lf[2 * k] - 2 * lf[k]) % P1                       # log C(2k,k)
    dH = (H[2 * k] - H[k]) % p
    return dict(pw=pw, dl=dl, H=H, Lc=Lc, dH=dH, n=n, k=k)


def law_defect(p, D, m):
    if m % p == 0:
        return None
    P1 = p - 1
    e = (3 * D['Lc'] - D['k'] * int(D['dl'][m % p])) % P1
    t = D['pw'][e]
    F = int(t.sum() % p)
    h = int(((t * D['dH']) % p).sum() % p)
    return (6 * h - fermat_q(m, p) * F) % p, F


def K_value(p, D):
    n = D['n']
    P1 = p - 1
    j = np.arange(n // 2 + 1, dtype=np.int64)
    l4 = int(D['dl'][4])
    la = (D['Lc'] - D['k'] * l4) % P1                         # log a_k
    e = (3 * la[j] + la[2 * j]) % P1
    H = D['H']
    br = (H[2 * j] - 2 * H[j] - 2 * H[4 * j] + 2 * H[n - j] + H[n - 2 * j]) % p
    return int(((D['pw'][e] * br) % p).sum() % p)


def twins_modp2(p):
    """class 3, 6, 7 twins directly mod p^2, and the class-6 reduction identity."""
    pp = p * p
    n = (p - 1) // 2
    Cc = [comb(2 * j, j) % pp for j in range(n + 1)]
    C2nj = [comb(2 * n - j, j) % pp for j in range(n + 1)]
    F1 = sum(pow(c, 3, pp) for c in Cc) % pp
    T1 = sum(Cc[j] * C2nj[j] ** 2 for j in range(n + 1)) % pp
    c7a = sum(comb(2 * n - 2 * k, n - k) ** 2 * comb(n + k, 2 * k) * comb(2 * n, n + k) for k in range(n + 1)) % pp
    c7b = sum(comb(n + k, 2 * k) ** 3 * comb(2 * n, n + k) for k in range(n + 1)) % pp
    u6a = sum(comb(n + k, k) ** 3 * comb(n, 2 * k) for k in range(n // 2 + 1)) % pp
    u6b = sum(comb(n + k, k) * comb(2 * k, k) * comb(2 * k, n) * comb(2 * n - 2 * k, n - k)
              for k in range((n + 1) // 2, n + 1)) % pp
    sg = -1 if n % 2 else 1
    # K sum mod p with exact rationals reduced mod p
    inv = [0] + [pow(i, -1, p) for i in range(1, p)]
    H = [0] * p
    for i in range(1, p):
        H[i] = (H[i - 1] + inv[i]) % p
    a = [Cc[k] % p * pow(4, -k, p) % p for k in range(n + 1)]
    K = sum(pow(a[j], 3, p) * a[2 * j] * (H[2 * j] - 2 * H[j] - 2 * H[4 * j] + 2 * H[n - j] + H[n - 2 * j])
            for j in range(n // 2 + 1)) % p
    d6 = (u6a - sg * u6b) % pp
    return dict(c3=(T1 - F1) % pp == 0, c7=(c7a - c7b) % pp == 0, c6=d6 == 0,
                c6_identity=(d6 - p * K) % pp == 0, d6_over_p=d6 // p if d6 % p == 0 else None, K=K)


def main():
    PA = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
    PG = int(sys.argv[2]) if len(sys.argv) > 2 else 3000
    PT = int(sys.argv[3]) if len(sys.argv) > 3 else 1200
    t0 = time.time()
    windows = [(5, PA), (199500, 200000)]
    for lo, hi in windows:
        ps = primes_in(lo, hi)
        fails = {m: [] for m in CM}
        ss = {m: 0 for m in CM}
        kf = []
        cnt = {m: 0 for m in CM}
        for p in ps:
            if p < 5:
                continue
            D = prime_data(p)
            for m in CM:
                r = law_defect(p, D, m)
                if r is None:
                    continue
                cnt[m] += 1
                if r[0]:
                    fails[m].append(p)
                if r[1] == 0:
                    ss[m] += 1
            if K_value(p, D):
                kf.append(p)
        print('[D1] primes in [%d, %d]: %d primes' % (lo, hi, len(ps)))
        for m in CM:
            print('     m = %5d: tested %5d, (L_p) failures %s; F == 0 mod p at %d (%.3f)'
                  % (m, cnt[m], fails[m][:8] if fails[m] else 'none', ss[m], ss[m] / max(cnt[m], 1)))
        print('     (K) failures: %s   [%.0f s]' % (kf[:8] if kf else 'none', time.time() - t0), flush=True)
    # generic m
    gen = [m for m in range(-100, 101) if m not in CM and m != 0]
    ps = primes_in(5, PG)
    res = {}
    for p in ps:
        D = prime_data(p)
        for m in gen:
            r = law_defect(p, D, m)
            if r is None:
                continue
            a, b = res.get(m, (0, 0))
            res[m] = (a + 1, b + (r[0] == 0))
    worst = max(res.items(), key=lambda kv: kv[1][1] / kv[1][0])
    frac = sorted(b / a for a, b in res.values())
    print('[D2] generic m (all 0 < |m| <= 100 except the five CM values in range), primes 5..%d:' % PG)
    print('     max pass fraction %.3f (m = %d: %d of %d); median pass fraction %.3f; min %.3f'
          % (worst[1][1] / worst[1][0], worst[0], worst[1][1], worst[1][0], frac[len(frac) // 2], frac[0]))
    for m in (3, 5, 7, 10, 100, -3, 2, -1, 32, -100):
        a, b = res[m]
        print('     m = %4d: passes %d of %d (failures %d)' % (m, b, a, a - b))
    # direct mod p^2 twins
    ok3 = ok6 = ok7 = okid = 0
    bad = []
    pts = primes_in(5, PT)
    for p in pts:
        r = twins_modp2(p)
        ok3 += r['c3']
        ok6 += r['c6']
        ok7 += r['c7']
        okid += r['c6_identity']
        if not (r['c3'] and r['c6'] and r['c7'] and r['c6_identity']):
            bad.append(p)
    print('[D3] direct mod p^2 twins, %d primes 5..%d: class 3 %d, class 7 %d, class 6 %d hold; class-6 reduction identity'
          ' u6a - (-1)^n u6b == p*K (mod p^2) holds at %d; exceptions %s' % (len(pts), PT, ok3, ok7, ok6, okid, bad[:10]))
    print('total %.0f s' % (time.time() - t0))


if __name__ == '__main__':
    main()
