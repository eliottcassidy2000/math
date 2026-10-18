#!/usr/bin/env python3
"""S43b referee: p-shift (Lucas-block) lemma, (B1) <=> (L_p), F_{p^2} vs F_p^2, KLMSY Lemma 10, the Conjecture-9 (s = 1)
corollary (D) of THM-4686, and THM-4686's inputs (reflection (R), KLMSY Theorem 2 on their model E_lambda, and the
2-isogeny between E_lambda and the -2 twist of E_{1/lambda}).  Written from scratch.

Notation: t_k = C(2k,k)^3 m^-k, F = sum_{k<=n} t_k, h = sum_{k<=n} t_k (H_2k - H_k), q = q_p(m), n = (p-1)/2,
F_N = sum_{k<=N} t_k ("a_N" in KLMSY), lambda = 64/m.
Usage: python audit_lucas_klmsy.py
"""
import sys
from math import comb
from fractions import Fraction

CM = [1, -8, 16, -64, 64, 256, -512, 4096]
GEN = [3, 5, 7, 10, -3, 2, -1, 100]


def primes_in(a, b):
    return [p for p in range(max(a, 2), b + 1) if all(p % d for d in range(2, int(p ** 0.5) + 1))]


def fq(m, p):
    pp = p * p
    return ((pow(m % pp, p - 1, pp) - 1) // p) % p


def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def central_modp2(p, K):
    """t-free data: list of (v, u) with C(2k,k) = p^v * unit, unit == u (mod p^2), for k = 0..K."""
    pp = p * p
    out = [(0, 1)]
    v, u = 0, 1
    for k in range(1, K + 1):
        for x, sgn in ((2 * k, 1), (2 * k - 1, 1), (k, -1), (k, -1)):
            e = 0
            while x % p == 0:
                x //= p
                e += 1
            v += sgn * e
            u = u * (x if sgn > 0 else pow(x, -1, pp)) % pp
        out.append((v, u))
    return out


def t_modp2(cen, m, p, K):
    """t_k mod p^2 for k <= K (t_k == 0 mod p^2 whenever p | C(2k,k), since then p^3 | C^3)."""
    pp = p * p
    w = pow(m % pp, -1, pp)
    t = []
    pw = 1
    for k in range(K + 1):
        v, u = cen[k]
        t.append(0 if v >= 1 else pow(u, 3, pp) * pw % pp)
        pw = pw * w % pp
    return t


def Fh(p, m):
    n = (p - 1) // 2
    pp = p * p
    H = [0] * p
    for i in range(1, p):
        H[i] = (H[i - 1] + pow(i, -1, p)) % p
    w = pow(m % pp, -1, pp)
    F = 0
    h = 0
    pw = 1
    for k in range(n + 1):
        t = comb(2 * k, k) ** 3 % pp * pw % pp
        F = (F + t) % pp
        h = (h + t * (H[2 * k] - H[k])) % p
        pw = pw * w % pp
    return F, h, H


def main():
    out = []
    # [1] p-shift lemma for all digits a in [0, 2p] (includes 2a >= p and a >= p), j <= n
    bad = 0
    tot = 0
    for p in primes_in(5, 31):
        pp = p * p
        n = (p - 1) // 2
        H = [0] * p
        for i in range(1, p):
            H[i] = (H[i - 1] + pow(i, -1, p)) % p
        for a in range(0, 2 * p + 1):
            for j in range(n + 1):
                N = a * p + j
                lhs = comb(2 * N, N) % pp
                rhs = comb(2 * a, a) * comb(2 * j, j) * (1 + 2 * a * p * (H[2 * j] - H[j])) % pp
                tot += 1
                bad += (lhs - rhs) % pp != 0
    print('[1] p-shift lemma C(2(ap+j),ap+j) == C(2a,a)C(2j,j)(1+2ap(H_2j-H_j)) mod p^2, 5<=p<=31, 0<=a<=2p, j<=n:'
          ' %d checks, %d failures' % (tot, bad))
    # [2] block identity block_a == t_a (F + a p (6h - qF)) and (B1) <=> (L_p)
    bad2 = 0
    disagree = 0
    rows = 0
    for p in primes_in(5, 61):
        pp = p * p
        n = (p - 1) // 2
        K = 2 * p * p + p
        cen = central_modp2(p, (2 * p + 1) * p + n)
        for m in CM + GEN:
            if m % p == 0:
                continue
            F, h, H = Fh(p, m)
            q = fq(m, p)
            t = t_modp2(cen, m, p, (2 * p + 1) * p + n)
            for a in range(0, 2 * p + 2):
                blk = sum(t[a * p + j] for j in range(n + 1)) % pp
                ta = t[a] if a < len(t) else None
                pred = ta * (F + a * p * ((6 * h - q * F) % p)) % pp
                if (blk - pred) % pp:
                    bad2 += 1
            B1 = (sum(t[k] for k in range(p, p + n + 1)) - 8 * pow(m % pp, -1, pp) * F) % pp == 0
            L = (6 * h - q * F) % p == 0
            rows += 1
            disagree += (B1 != L)
    print('[2] block_a == t_a (F + a p (6h - qF)) mod p^2 for 0<=a<=2p+1, p<=61, 8 CM + 8 generic m: %d failures;'
          ' (B1) <=> (L_p) truth table: %d (m,p) rows, %d disagreements' % (bad2, rows, disagree))
    # [3] KLMSY Lemma 10: a_{k+pn} == a_k a_n + p b_k sum_{i<=n} i t_i (mod p^2), (p-1)/2 <= k < p, all n <= 2p
    bad3 = 0
    tot3 = 0
    for p in primes_in(5, 23):
        pp = p * p
        K = p - 1 + p * 2 * p
        cen = central_modp2(p, K)
        H = [Fraction(0)]
        for i in range(1, 2 * p):
            H.append(H[-1] + Fraction(1, i))
        for m in (3, 5, -7, 1, 4096):
            if m % p == 0:
                continue
            t = t_modp2(cen, m, p, K)
            a = [0] * (K + 1)
            s = 0
            for k in range(K + 1):
                s = (s + t[k]) % pp
                a[k] = s
            q = fq(m, p)
            # b_k mod p for k < p: the terms with (p-1)/2 < i < p are 0 mod p (p^3 | C^3 against one p in H_2i)
            def bmod(k):
                tot_ = 0
                w = pow(m % p, -1, p)
                for i in range(min(k, (p - 1) // 2) + 1):
                    hh = H[2 * i] - H[i]
                    tot_ += comb(2 * i, i) ** 3 * pow(w, i, p) * (6 * hh.numerator * pow(hh.denominator, -1, p) - q)
                return tot_ % p
            S = [0]
            for i in range(1, 2 * p + 1):
                S.append((S[-1] + i * t[i]) % p)
            for k in range((p - 1) // 2, p):
                bk = bmod(k)
                for nn in range(0, 2 * p + 1):
                    tot3 += 1
                    if (a[k + p * nn] - a[k] * a[nn] - p * bk * S[nn]) % pp:
                        bad3 += 1
    print('[3] KLMSY Lemma 10 a_{k+pn} == a_k a_n + p b_k sum_{i<=n} i t_i (mod p^2), (p-1)/2<=k<p, n<=2p, 5<=p<=23:'
          ' %d checks, %d failures' % (tot3, bad3))
    # [4] F_{p^2} vs F_p^2 and the identity F_{p^2} - F_p^2 == p (6h - qF) sum_{a<=n} a t_a
    ok4 = {m: 0 for m in CM}
    cnt4 = {m: 0 for m in CM}
    weak = {m: [] for m in CM}
    badid = 0
    gen_eq_law = 0
    gen_rows = 0
    for p in primes_in(5, 151):
        pp = p * p
        n = (p - 1) // 2
        cen = central_modp2(p, p * p - 1)
        for m in CM + GEN:
            if m % p == 0:
                continue
            t = t_modp2(cen, m, p, p * p - 1)
            Fp = sum(t[:p]) % pp
            Fp2 = sum(t) % pp
            F, h, H = Fh(p, m)
            q = fq(m, p)
            S = sum(a * t[a] for a in range(n + 1)) % p
            if (Fp2 - Fp * Fp - p * ((6 * h - q * F) * S % p)) % pp:
                badid += 1
            if m in CM:
                cnt4[m] += 1
                ok4[m] += (Fp2 - Fp * Fp) % pp == 0
                if S == 0:
                    weak[m].append(p)
            else:
                gen_rows += 1
                gen_eq_law += ((Fp2 - Fp * Fp) % pp == 0) == ((6 * h - q * F) % p == 0)
    print('[4] identity F_{p^2} - F_p^2 == p(6h - qF) sum_{a<=n} a t_a (mod p^2), p<=151, 16 values of m: %d failures'
          % badid)
    for m in CM:
        print('     m = %5d: F_{p^2} == F_p^2 mod p^2 at %d of %d primes 5..151; primes with sum a t_a == 0 mod p'
              ' (converse unavailable): %s' % (m, ok4[m], cnt4[m], weak[m]))
    print('     generic m: [F_{p^2} == F_p^2] agrees with [(L_p)] in %d of %d (m,p) rows' % (gen_eq_law, gen_rows))
    # [5] THM-4686 corollary (D): F_{(Mp-1)/2} == F_{(p-1)/2} F_{(M-1)/2} mod p^2, odd M <= 2p+3 (digits a >= p), p <= 61
    bad5 = 0
    tot5 = 0
    genfail = 0
    gentot = 0
    for p in primes_in(5, 61):
        pp = p * p
        Ms = list(range(1, 2 * p + 4, 2))
        K = (max(Ms) * p - 1) // 2
        cen = central_modp2(p, K)
        for m in CM + [3, 5, 7, 10]:
            if m % p == 0:
                continue
            t = t_modp2(cen, m, p, K)
            pref = []
            s = 0
            for k in range(K + 1):
                s = (s + t[k]) % pp
                pref.append(s)
            for M in Ms:
                lhs = pref[(M * p - 1) // 2]
                rhs = pref[(p - 1) // 2] * pref[(M - 1) // 2] % pp
                if m in CM:
                    tot5 += 1
                    bad5 += (lhs - rhs) % pp != 0
                elif M > 1:
                    gentot += 1
                    genfail += (lhs - rhs) % pp != 0
    print('[5] (D) F_{(Mp-1)/2} == F_{(p-1)/2} F_{(M-1)/2} mod p^2, all odd M <= 2p+3, 5<=p<=61: CM %d checks, %d failures;'
          ' generic m = 3,5,7,10: %d of %d fail' % (tot5, bad5, genfail, gentot))
    # [6] THM-4686 inputs: (R), KLMSY Theorem 2 on E_lambda: y^2 = (x-1)(x^2 - 1/(1-lambda)), and the twist relation
    lams = [Fraction(64), Fraction(1, 64), Fraction(-8), Fraction(-1, 8), Fraction(4), Fraction(1, 4), Fraction(-1)]
    badR = 0
    badT2 = []
    badT2c = []
    nT2 = 0
    badtw = []
    for p in primes_in(5, 1000):
        pp = p * p
        n = (p - 1) // 2
        a3 = []
        for k in range(n + 1):
            a3.append(comb(2 * k, k) ** 3 * pow(64, -k, pp) % pp)       # a_k^3 = C(2k,k)^3/64^k
        O = [0] * (n + 1)
        for k in range(1, n + 1):
            O[k] = (O[k - 1] + pow(2 * k - 1, -1, p)) % p
        ap = {}
        for lam in lams:
            if lam.numerator % p == 0 or lam.denominator % p == 0 or (lam.denominator - lam.numerator) % p == 0:
                continue
            z = lam.numerator * pow(lam.denominator, -1, pp) % pp
            zi = pow(z, -1, pp)
            F = sum(a3[k] * pow(z, k, pp) for k in range(n + 1)) % pp
            Fi = sum(a3[k] * pow(zi, k, pp) for k in range(n + 1)) % pp
            Gi = sum(a3[k] * pow(zi, k, pp) % pp * (O[n - k] - O[k]) for k in range(n + 1)) % pp
            if (F - (-1) ** n * pow(z, n, pp) * (Fi + 3 * p * Gi)) % pp:
                badR += 1
            # point count on E_lambda (KLMSY model), good reduction away from p | 2*3*7*num*den
            c = pow((1 - z) % p, -1, p)
            if p > 7:
                s = 0
                for x in range(p):
                    s += legendre((x - 1) * (x * x - c), p)
                a_p = -s
                ap[lam] = a_p
                rhs = 0 if a_p % p == 0 else legendre((z - 1) % p, p) * (a_p * a_p - 2 * p)
                if (F - rhs) % pp:
                    badT2.append((str(lam), p))
                rhs2 = 0 if a_p % p == 0 else legendre((1 - z) % p, p) * (a_p * a_p - 2 * p)
                if (F - rhs2) % pp:
                    badT2c.append((str(lam), p))
                nT2 += 1
        for lam in lams:
            li = 1 / lam
            if lam in ap and li in ap and lam != li:
                if ap[lam] != legendre(-2, p) * ap[li]:
                    badtw.append((str(lam), p))
    print('[6] THM-4686: (R) mod p^2 at z in {64,1/64,-8,-1/8,4,1/4,-1}, 5<=p<1000: %d failures' % badR)
    print('    KLMSY Theorem 2 on E_lambda: y^2=(x-1)(x^2-1/(1-lambda)), %d (lambda,p) pairs, 7<p<1000, alpha^2 == a_p^2 - 2p:' % nT2)
    print('      as printed, F(lambda) == ((lambda-1)/p) alpha^2 (0 if supersingular): %d failures, first %s' % (len(badT2), badT2[:6]))
    print('      with ((1-lambda)/p) instead: %d failures %s' % (len(badT2c), badT2c[:6]))
    print('    a_p(E_lambda) == (-2/p) a_p(E_{1/lambda}) (E_lambda 2-isogenous to the -2 twist of E_{1/lambda}),'
          ' 7<p<1000: failures %s' % (badtw[:6] if badtw else 'none'))


if __name__ == '__main__':
    main()
