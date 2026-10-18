#!/usr/bin/env python3
"""Mechanisms behind the half-period twin classes (S43, after halfperiod_twin_class_forms_20261018.py).

With n = (p-1)/2, a_k = C(2k,k)/4^k, O_k = sum_{j<=k} 1/(2j-1), every library binomial has a p-adic shadow mod p^2:
    C(n,k) == (-1)^k a_k (1 - p O_k),  C(n+k,k) == a_k (1 + p O_k),  C(n+k,2k) == (-1)^k a_k / 4^k,
    C(2n,j) == (-1)^j (1 - p H_j),  C(n,2k) == a_2k (1 - p O_2k),  C(2n-j, j) == (-1)^j C(2j,j)(1 - p(H_2j - H_j)).
The twin relations then reduce to harmonic congruences mod p:
    classes 1, 2, 4 (4C):  sum_{k<=n} a_k^3 O_k == 0                                   (THM-4670 lemma)
    classes 3, 7  (7AB):   (H)  sum_{j<=n} C(2j,j)^3 (H_2j - H_j) == 0
    class 6       (4A):    (K)  sum_{j<=n/2} a_j^3 a_2j (H_2j - 2H_j - 2H_4j + 2H_(n-j) + H_(n-2j)) == 0
This script (1) checks the shadow formulas exactly mod p^2 for p < 200, (2) checks (H) and (K) for p < PMAX,
(3) prints the single harmonic sums sum C(2j,j)^3 H_j and sum C(2j,j)^3 H_2j against a_7(p) and q_p(2), q_p(3),
(4) prints the first terms of the class-3/6/7 sequences and guesses P-recurrences.
Usage: python halfperiod_twin_mechanisms_20261018.py [PMAX]
Session opus-2026-10-17-S43.
"""
import os, sys, json
from math import comb
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from halfperiod_twin_class_forms_20261018 import eta_quotient


def C(a, b):
    if b < 0 or a < 0 or b > a:
        return 0
    return comb(a, b)


def primes_upto(N):
    s = bytearray([1]) * (N + 1)
    s[0] = s[1] = 0
    for i in range(2, int(N ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(N + 1) if s[i]]


def harmonic_table(p, M):
    H = [0] * (M + 1)
    for i in range(1, M + 1):
        H[i] = (H[i - 1] + pow(i, p - 2, p)) % p
    return H


def central_mod(p, M):
    """C(2j,j) mod p for j <= M (M <= (p-1)/2, no p in numerators)."""
    c = [1] * (M + 1)
    for j in range(1, M + 1):
        c[j] = c[j - 1] * (2 * j) * (2 * j - 1) % p * pow(j * j, p - 2, p) % p
    return c


def main():
    PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
    out = {}
    # (1) shadow formulas mod p^2
    bad = []
    for p in primes_upto(200)[2:]:
        n = (p - 1) // 2
        pp = p * p
        def inv(x): return pow(x % pp, -1, pp)
        a = [Fraction(comb(2 * k, k), 4 ** k) for k in range(2 * n + 2)]
        O = [Fraction(0)]
        for j in range(1, 2 * n + 2):
            O.append(O[-1] + Fraction(1, 2 * j - 1))
        H = [Fraction(0)]
        for j in range(1, 2 * n + 2):
            H.append(H[-1] + Fraction(1, j))
        def md(fr): return fr.numerator * inv(fr.denominator) % pp
        for k in range(n + 1):
            checks = [(C(n, k), (-1) ** k * a[k] * (1 - p * O[k])), (C(n + k, k), a[k] * (1 + p * O[k])),
                      (C(n + k, 2 * k), (-1) ** k * a[k] / 4 ** k)]
            if 2 * k <= n:
                checks.append((C(n, 2 * k), a[2 * k] * (1 - p * O[2 * k])))
            if 2 * k < p:
                checks.append((C(2 * n - k, k), (-1) ** k * comb(2 * k, k) * (1 - p * (H[2 * k] - H[k]))))
            for j in (k, n + k):
                if j <= 2 * n:
                    checks.append((C(2 * n, j), (-1) ** j * (1 - p * H[j])))
            for lhs, rhs in checks:
                if (lhs - md(rhs)) % pp:
                    bad.append((p, k))
    print('(1) shadow formulas mod p^2 for 5 <= p < 200: %s' % ('all hold' if not bad else 'FAIL %s' % bad[:5]))
    out['shadows_ok'] = not bad
    # (2) (H) and (K), (3) single harmonic sums
    a7 = eta_quotient({1: 3, 7: 3}, 2000)
    Hfail, Kfail, Afail = [], [], []
    rows = []
    for p in primes_upto(PMAX)[2:]:
        n = (p - 1) // 2
        Ht = harmonic_table(p, 2 * n)
        cb = central_mod(p, n)
        sH = sum(pow(cb[j], 3, p) * (Ht[2 * j] - Ht[j]) for j in range(n + 1)) % p
        if sH:
            Hfail.append(p)
        # a_k = C(2k,k)/4^k mod p, a_2k needs C(4k,2k) with 4k <= 2n: fine for 2k <= n
        inv4 = pow(4, p - 2, p)
        Hn4 = harmonic_table(p, 2 * n)
        def a_mod(k):
            return cb_full[k] * pow(inv4, k, p) % p
        cb_full = central_mod(p, n)
        sK = 0
        for j in range(n // 2 + 1):
            term = pow(a_mod(j), 3, p) * a_mod(2 * j) % p
            br = (Hn4[2 * j] - 2 * Hn4[j] - 2 * Hn4[4 * j] + 2 * Hn4[n - j] + Hn4[n - 2 * j]) % p
            sK = (sK + term * br) % p
        if sK:
            Kfail.append(p)
        if p < 2000 and p != 7:
            s1 = sum(pow(cb[j], 3, p) * Ht[j] for j in range(n + 1)) % p
            s0 = sum(pow(cb[j], 3, p) for j in range(n + 1)) % p
            ap = a7[p - 1] % p
            q2 = (pow(2, p - 1, p * p) - 1) // p % p
            q3 = (pow(3, p - 1, p * p) - 1) // p % p
            if (s0 - ap) % p:
                Afail.append(p)
            rows.append((p, ap, s1, q2, q3))
    print('(2) (H) sum_{j<=n} C(2j,j)^3 (H_2j - H_j) == 0 mod p for 5 <= p < %d: %s' % (PMAX, 'holds' if not Hfail else 'FAILS at %s' % Hfail[:10]))
    print('    (K) class-6 harmonic congruence for 5 <= p < %d: %s' % (PMAX, 'holds' if not Kfail else 'FAILS at %s' % Kfail[:10]))
    print('    mod-p check sum C(2j,j)^3 == a_7(p): %s' % ('holds' if not Afail else 'fails %s' % Afail[:5]))
    out['H_fail'] = Hfail
    out['K_fail'] = Kfail
    # (3) fit sum C(2j,j)^3 H_j == (alpha q2 + beta q3 + gamma) a_7(p) mod p with small rationals
    fits = []
    cands = [Fraction(i, d) for d in (1, 2, 3, 4, 6, 8) for i in range(-24, 25)]
    cands = sorted(set(cands))
    for al in cands:
        for ga in cands:
            ok = True
            for p, ap, s1, q2, q3 in rows:
                rhs = (al.numerator * pow(al.denominator, -1, p) * q2 + ga.numerator * pow(ga.denominator, -1, p)) * ap
                if (s1 - rhs) % p:
                    ok = False
                    break
            if ok:
                fits.append((str(al), str(ga)))
    print('(3) sum_{j<=n} C(2j,j)^3 H_j == (alpha q_p(2) + gamma) a_7(p) mod p, fits (alpha, gamma):', fits)
    out['H_j_fits'] = fits
    # (4) first terms of the class 3/6/7 sequences
    seqs = {
        'g(2n)=sum C(2j,j) C(2n-j,j)^2': lambda n: sum(C(2 * j, j) * C(2 * n - j, j) ** 2 for j in range(n + 1)),
        'sum C(2k,k)^3 (k<=n)': lambda n: sum(C(2 * k, k) ** 3 for k in range(n + 1)),
        'c6a=sum C(n+k,k)^3 C(n,2k)': lambda n: sum(C(n + k, k) ** 3 * C(n, 2 * k) for k in range(n + 1)),
        'c6b=sum C(n+k,k)C(2k,k)C(2k,n)C(2n-2k,n-k)': lambda n: sum(C(n + k, k) * C(2 * k, k) * C(2 * k, n) * C(2 * n - 2 * k, n - k) for k in range(n + 1)),
        'c7a=sum C(2n-2k,n-k)^2 C(n+k,2k) C(2n,n+k)': lambda n: sum(C(2 * n - 2 * k, n - k) ** 2 * C(n + k, 2 * k) * C(2 * n, n + k) for k in range(n + 1)),
        'c7b=sum C(n+k,2k)^3 C(2n,n+k)': lambda n: sum(C(n + k, 2 * k) ** 3 * C(2 * n, n + k) for k in range(n + 1)),
    }
    out['terms'] = {}
    for nm, f in seqs.items():
        t = [f(n) for n in range(16)]
        out['terms'][nm] = t
        print('(4) %s: %s' % (nm, t[:12]))
        # guess a P-recurrence sum_{i=0}^{r} P_i(n) u(n+i) = 0, deg d, via exact linear algebra over Q
        for r in (1, 2, 3):
            found = None
            for d in (1, 2, 3, 4):
                nunk = (r + 1) * (d + 1)
                N = 40
                vals = [f(n) for n in range(N + r + 1)]
                rowsM = []
                for n0 in range(N):
                    rowsM.append([Fraction(vals[n0 + i] * n0 ** e) for i in range(r + 1) for e in range(d + 1)])
                if len(rowsM) < nunk + 4:
                    continue
                # nullspace via Gaussian elimination
                M = [row[:] for row in rowsM]
                ncol = nunk
                piv = []
                rI = 0
                for c in range(ncol):
                    pr = next((i for i in range(rI, len(M)) if M[i][c] != 0), None)
                    if pr is None:
                        continue
                    M[rI], M[pr] = M[pr], M[rI]
                    pv = M[rI][c]
                    M[rI] = [x / pv for x in M[rI]]
                    for i in range(len(M)):
                        if i != rI and M[i][c] != 0:
                            fac = M[i][c]
                            M[i] = [x - fac * y for x, y in zip(M[i], M[rI])]
                    piv.append(c)
                    rI += 1
                if len(piv) < ncol:
                    free = [c for c in range(ncol) if c not in piv][0]
                    vec = [Fraction(0)] * ncol
                    vec[free] = Fraction(1)
                    for i, c in enumerate(piv):
                        vec[c] = -M[i][free]
                    found = (r, d, [str(x) for x in vec])
                    break
            if found:
                print('      P-recurrence order %d degree %d: coefficients (u(n+i) n^e, i-major) %s' % found)
                out.setdefault('recurrences', {})[nm] = found
                break
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(out, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'halfperiod_twin_mechanisms_20261018.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
