#!/usr/bin/env python3
"""S43b referee (written from scratch, imports no session module): THM-4685 Lemma 1, Lemma 2, part (i), Lemma 3,
part (ii) and the twin form, checked in exact modular arithmetic.

All p-adic quantities are p-adic integers here (denominators are units: powers of 2, odd numbers < p, m with p !| m),
so every congruence is tested on integer representatives mod p^2 (or mod p).

Usage: python audit_thm4685.py [P_MODP2] [P_LAW]
  P_MODP2: bound for the mod p^2 statements (Lemma 1, Lemma 2(d), part (i), twins)       default 3000
  P_LAW:   bound for the mod p statements at m = +-64 and the reflection identities       default 6000
"""
import sys
import random
from math import comb


def primes_in(a, b):
    s = bytearray([1]) * (b + 1)
    s[0:2] = b'\x00\x00'
    for i in range(2, int(b ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(range(i * i, b + 1, i)))
    return [i for i in range(max(a, 2), b + 1) if s[i]]


def fermat_q(m, p):
    pp = p * p
    return ((pow(m % pp, p - 1, pp) - 1) // p) % p


def tables(p, M):
    """inverses mod p of 1..M (M < p), harmonic numbers H_0..H_M mod p, odd sums O_0..O_{(M+1)//2} mod p."""
    inv = [0] * (M + 1)
    if M >= 1:
        inv[1] = 1
    for i in range(2, M + 1):
        inv[i] = (-(p // i) * inv[p % i]) % p
    H = [0] * (M + 1)
    for i in range(1, M + 1):
        H[i] = (H[i - 1] + inv[i]) % p
    return inv, H


def check_prime_modp2(p, report):
    """Lemma 1 (S1-S6, M1-M4, consequences), Lemma 2 (a)-(d), part (i), twin form identity, all at one prime."""
    pp = p * p
    n = (p - 1) // 2
    inv_p, H = tables(p, p - 1)               # H_j mod p, j <= p-1
    O = [0] * (n + 1)
    for k in range(1, n + 1):
        O[k] = (O[k - 1] + inv_p[2 * k - 1]) % p
    inv2 = pow(2, -1, pp)
    # a_k = prod_{j<=k} (2j-1)/(2j) mod p^2
    a = [1] * (n + 1)
    for k in range(1, n + 1):
        a[k] = a[k - 1] * (2 * k - 1) % pp * pow(2 * k, -1, pp) % pp
    inv4 = pow(4, -1, pp)
    # exact binomials mod p^2 (math.comb is exact)
    Cn = [comb(n, k) % pp for k in range(n + 1)]
    Cnk = [comb(n + k, k) % pp for k in range(n + 1)]
    Cnk2 = [comb(n + k, 2 * k) % pp for k in range(n + 1)]
    C2n = [comb(2 * n, j) % pp for j in range(2 * n + 1)]
    C2nj = [comb(2 * n - j, j) % pp for j in range(n + 1)]
    Ccen = [comb(2 * j, j) % pp for j in range(n + 1)]
    bad = []
    q2 = fermat_q(2, p)
    for k in range(n + 1):
        sg = -1 if k % 2 else 1
        if (Cn[k] - sg * a[k] * (1 - p * O[k])) % pp:
            bad.append(('S1', k))
        if (Cnk[k] - a[k] * (1 + p * O[k])) % pp:
            bad.append(('S2', k))
        if (Cnk2[k] - sg * a[k] * pow(inv4, k, pp)) % pp:
            bad.append(('S3', k))
        if (C2nj[k] - sg * Ccen[k] * (1 - p * (H[2 * k] - H[k]))) % pp:
            bad.append(('S5', k))
        sgn_n = -1 if n % 2 else 1
        if (a[n - k] - sgn_n * a[k] * (1 + p * (O[n - k] - O[k]))) % pp:
            bad.append(('S6', k))
        # mod p statements
        if (O[k] - (H[n - k] - H[n]) * inv_p[2]) % p:
            bad.append(('M1', k))
        if (a[k] - sg * Cn[k]) % p:
            bad.append(('M4', k))
        if ((O[n - k] - O[k]) - (H[k] - H[n - k]) * inv_p[2]) % p:
            bad.append(('Delta', k))
        if ((H[2 * k] - H[k]) - (H[n - k] - H[k] - H[n]) * inv_p[2]) % p:
            bad.append(('H2k-Hk', k))
        if (H[n + k] - H[n - k]) % p:
            bad.append(('Hn+k', k))
    for j in range(2 * n + 1):
        sg = -1 if j % 2 else 1
        if (C2n[j] - sg * (1 - p * H[j])) % pp:
            bad.append(('S4', j))
    for i in range(p):
        if (H[p - 1 - i] - H[i]) % p:
            bad.append(('M2', i))
    if (H[n] + 2 * q2) % p:
        bad.append(('M3', 0))
    # Lemma 2: D_n, S_n, T_n mod p
    Dn = sum((-1) ** k * Cn[k] ** 3 for k in range(n + 1)) % p
    Sn = sum((-1) ** k * Cn[k] ** 3 * H[n - k] for k in range(n + 1)) % p
    Tn = sum((-1) ** k * Cn[k] ** 3 * H[k] for k in range(n + 1)) % p
    F2 = sum(pow(a[k], 3, pp) for k in range(n + 1)) % pp      # F mod p^2
    F1 = F2 % p
    if (Tn - (-1) ** n * Sn) % p:
        bad.append(('T=(-1)^n S', 0))
    if n % 2 == 1 and (Dn % p or Sn % p):
        bad.append(('THM-4670 odd case', 0))
    if n % 2 == 0 and (Sn - H[n] * Dn) % p:
        bad.append(('THM-4670 even case', 0))
    La = sum(pow(a[k], 3, p) * O[k] for k in range(n + 1)) % p
    Lb = sum(pow(a[k], 3, p) * O[n - k] for k in range(n + 1)) % p
    Lc = sum(pow(a[k], 3, p) * H[n - k] for k in range(n + 1)) % p
    if La:
        bad.append(('L2a', La))
    if Lb:
        bad.append(('L2b', Lb))
    if (Lc - Sn) % p or (Lc - H[n] * F1) % p:
        bad.append(('L2c', 0))
    if (F1 - Dn) % p:
        bad.append(('F==D_n mod p', 0))
    # (d): sum a_k^3 Delta_k == (T_n - S_n)/2 and, for odd n, F == 0 mod p^2
    sD = sum(pow(a[k], 3, p) * (O[n - k] - O[k]) for k in range(n + 1)) % p
    if (sD - (Tn - Sn) * inv_p[2]) % p:
        bad.append(('L2d sumDelta', 0))
    if n % 2 == 1 and F2 % pp:
        bad.append(('L2d F==0 mod p^2', F2))
    # 2F == -3p sum a^3 Delta (odd n) as an identity mod p^2
    if n % 2 == 1 and (2 * F2 + 3 * p * sD) % pp:
        bad.append(('L2d 2F identity', 0))
    # part (i)
    four_n = pow(4, n, pp)
    sixteen_n = pow(16, n, pp)
    sg_n = -1 if n % 2 else 1
    ia = sum(Cnk[k] * Ccen[k] * (comb(2 * n - 2 * k, n - k) % pp) for k in range(n + 1)) % pp
    ia_mid = sg_n * four_n * sum(pow(a[k], 3, pp) * (1 + p * O[n - k]) for k in range(n + 1)) % pp
    ib = sum(pow(Cn[k], 3, pp) * C2n[n + k] for k in range(n + 1)) % pp
    ib_mid = sg_n * sum(pow(a[k], 3, pp) * (1 - 3 * p * O[k] - p * H[n + k]) for k in range(n + 1)) % pp
    ic = sum(Cn[k] * Cnk[k] ** 2 * C2n[n + k] for k in range(n + 1)) % pp
    ic_mid = sg_n * sum(pow(a[k], 3, pp) * (1 + p * O[k] - p * H[n + k]) for k in range(n + 1)) % pp
    bn = sum(Cn[k] ** 2 * Cnk[k] for k in range(n + 1)) % pp
    tw = sum(pow(Cnk[k], 3, pp) for k in range(n + 1)) % pp
    if (ia - ia_mid) % pp:
        bad.append(('i-a intermediate', 0))
    if (ia - sg_n * four_n * F2) % pp:
        bad.append(('i-a', 0))
    if (ib - ib_mid) % pp:
        bad.append(('i-b intermediate', 0))
    if (ib - sg_n * sixteen_n * F2) % pp:
        bad.append(('i-b', 0))
    if (1 - p * H[n] - sixteen_n) % pp:
        bad.append(('1-pH_n == 16^n', 0))
    if (ic - ic_mid) % pp:
        bad.append(('i-c intermediate', 0))
    if (ic - sg_n * sixteen_n * F2) % pp:
        bad.append(('i-c', 0))
    if (bn - F2) % pp or (tw - F2) % pp:
        bad.append(('i-d', 0))
    # class-2 partner C(n, n/2)^2 [n even] vs (i-a)
    partner = comb(n, n // 2) ** 2 % pp if n % 2 == 0 else 0
    if (partner - ia) % pp:
        bad.append(('class2 partner', 0))
    # twin form identity T_m == F_m - 2p h(64/m) (mod p^2) at m = +-64 and at two random m
    twin = {}
    for m in (64, -64, random.randrange(2, 10 ** 6), -random.randrange(2, 10 ** 6)):
        if m % p == 0:
            continue
        im = pow(m % pp, -1, pp)
        Fm = 0
        Tm = 0
        hm = 0
        w = 1
        for j in range(n + 1):
            c3 = pow(Ccen[j], 3, pp)
            Fm = (Fm + c3 * w) % pp
            Tm = (Tm + Ccen[j] * C2nj[j] ** 2 % pp * w) % pp
            hm = (hm + c3 * w % p * (H[2 * j] - H[j])) % p
            w = w * im % pp
        if (Tm - Fm + 2 * p * hm) % pp:
            bad.append(('twin identity', m))
        if m in (64, -64):
            twin[m] = (pow(4, p - 1, pp) * Tm - Fm) % pp == 0
            if not twin[m]:
                bad.append(('4^(p-1) T == F at m', m))
    report.append((p, bad[:6]))
    return not bad, F2


def law_modp(p, report, nrand=3):
    """part (ii) and Lemma 3 at one prime: L(+-64); Lemma 3 at random z; defect identity delta(m) vs delta(4096/m)."""
    n = (p - 1) // 2
    inv_p, H = tables(p, p - 1)
    # C(2k,k) mod p, k <= n
    c = [1] * (n + 1)
    for k in range(1, n + 1):
        c[k] = c[k - 1] * (2 * k) % p * (2 * k - 1) % p * inv_p[k] % p * inv_p[k] % p
    c3 = [x * x % p * x % p for x in c]
    inv64 = pow(64, -1, p)
    a3 = [c3[k] * pow(inv64, k, p) % p for k in range(n + 1)]       # a_k^3 mod p
    dH = [(H[2 * k] - H[k]) % p for k in range(n + 1)]

    def Fh(z):
        F = 0
        h = 0
        w = 1
        for k in range(n + 1):
            t = a3[k] * w % p
            F += t
            h += t * dH[k]
            w = w * z % p
        return F % p, h % p

    bad = []
    q2 = fermat_q(2, p)
    sgn = -1 if n % 2 else 1
    for z, m in ((1, 64), (p - 1, -64)):
        F, h = Fh(z)
        if (6 * h - fermat_q(m, p) * F) % p:
            bad.append(('L(%d)' % m, n % 2))
    # Lemma 3 at random z, and the defect identity
    for _ in range(nrand):
        # an integer m with p !| m; z = 64/m (mod p); the law's q_p(m) depends on m mod p^2, so keep m as an integer
        m = random.randrange(2, 10 ** 9)
        while m % p == 0:
            m = random.randrange(2, 10 ** 9)
        if random.random() < 0.5:
            m = -m
        z = 64 * pow(m % p, -1, p) % p
        zi = pow(z, -1, p)
        F, h = Fh(z)
        Fi, hi = Fh(zi)
        zn = pow(z, n, p)
        if (F - sgn * zn * Fi) % p:
            bad.append(('L3 F', z))
        if (h + sgn * zn * (hi + H[n] * Fi)) % p:
            bad.append(('L3 h', z))
        mm = 4096 * pow(m % (p * p), -1, p * p) % (p * p)       # the p-adic unit 4096/m, mod p^2
        dm = (6 * h - fermat_q(m, p) * F) % p
        dmm = (6 * hi - fermat_q(mm, p) * Fi) % p
        if (dm + sgn * zn * dmm) % p:
            bad.append(('defect identity', m))
    # q_p(4096/m) == 12 q_p(2) - q_p(m)
    for _ in range(2):
        m = random.randrange(2, 10 ** 6)
        if m % p == 0:
            continue
        mm = 4096 * pow(m, -1, p * p) % (p * p)
        if (fermat_q(mm, p) - 12 * q2 + fermat_q(m, p)) % p:
            bad.append(('q relation', m))
    if bad:
        report.append((p, bad[:6]))
    return not bad


def main():
    P2 = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
    PL = int(sys.argv[2]) if len(sys.argv) > 2 else 6000
    random.seed(20261018)
    rep = []
    ok_all = True
    zero_odd = 0
    nz_even = 0
    ps = primes_in(5, P2)
    for p in ps:
        ok, F2 = check_prime_modp2(p, rep)
        ok_all &= ok
        n = (p - 1) // 2
        if n % 2 == 1 and F2 == 0:
            zero_odd += 1
        if n % 2 == 0 and F2 % (p * p):
            nz_even += 1
    fails = [r for r in rep if r[1]]
    print('[A1] Lemma 1 (S1-S6 mod p^2, M1-M4 + consequences mod p), Lemma 2 (a)-(d), part (i) (i-a..i-d, with the',
          'intermediate congruences and exact signs (-1)^n, 4^n, 16^n), class-2 partner, twin identity:')
    print('     %d primes 5 <= p <= %d: %s' % (len(ps), P2, 'ALL HOLD' if not fails else 'FAILURES %s' % fails[:5]))
    print('     F == 0 mod p^2 at all %d primes p == 3 mod 4; F != 0 mod p^2 at %d of %d primes p == 1 mod 4'
          % (zero_odd, nz_even, sum(1 for p in ps if p % 4 == 1)))
    rep2 = []
    psl = primes_in(5, PL)
    good = sum(law_modp(p, rep2) for p in psl)
    print('[A2] part (ii): L(64), L(-64) (both parities of n), Lemma 3 (F and h) at 3 random z per prime, the defect',
          'identity delta(m) == -(-1)^n z^n delta(4096/m), q_p(4096/m) == 12 q_p(2) - q_p(m):')
    print('     %d primes 5 <= p <= %d: %d pass; failures: %s' % (len(psl), PL, good, rep2[:5] if rep2 else 'none'))
    # p = 3 sanity (excluded by hypothesis)
    p = 3
    pp = 9
    n = 1
    F = (1 + 8 * pow(64, -1, pp)) % pp
    ia = sum(comb(n + k, k) * comb(2 * k, k) * comb(2 * n - 2 * k, n - k) for k in range(2)) % pp
    print('[A3] p = 3 (excluded): F mod 9 = %d, (i-a) sum mod 9 = %d, (-1)^n 4^n F mod 9 = %d' % (F, ia, (-4 * F) % 9))


if __name__ == '__main__':
    main()
