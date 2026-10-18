#!/usr/bin/env python3
"""S43b referee: HYP-9311 (four Ramanujan families), written from scratch.

(1) the law sum_{k<p} A(k) m^-k w_k == q_p(m) sum_{k<p} A(k) m^-k (mod p) at the 33 CM parameters, 7 <= p <= PMAX,
    p !| 6m, with a per-term p-integrality check (H_N = R_N + H_{floor(N/p)}/p, R_N == H_{N mod p} mod p);
(2) the p-shift relation A(p+j) == A(1) (A(j) + p [A(j) w_j]) (mod p^2) for ALL 0 <= j < p (carries included);
(3) (B1') <=> (L_p^A) truth table at generic m;
(4) part (O): (B1') modulo p^3 at ordinary primes (F_p unit) / supersingular primes.
Usage: python audit_hyp9311.py [PMAX] [PMAX_O]
"""
import sys
from math import comb, isqrt
from fractions import Fraction

FAM = {
    'A1': ((lambda k: comb(2 * k, k) ** 3), [(6, 2), (-6, 1)]),
    'A2': ((lambda k: comb(2 * k, k) ** 2 * comb(3 * k, k)), [(3, 3), (2, 2), (-5, 1)]),
    'A3': ((lambda k: comb(2 * k, k) ** 2 * comb(4 * k, 2 * k)), [(4, 4), (-4, 1)]),
    'A4': ((lambda k: comb(2 * k, k) * comb(3 * k, k) * comb(6 * k, 3 * k)), [(6, 6), (-3, 3), (-3, 1)]),
}
CMP = {
    'A1': [-512, -64, -8, 1, 16, 64, 256, 4096],
    'A2': [-8640, -1728, -192, -27, 8, 64, 108, 216, 1458, 3375],
    'A3': [-82944, -12288, -3969, -1024, -144, 81, 256, 648, 2304, 20736],
    'A4': [-32768, -3375, 1728, 8000, 54000],
}
GEN = [2, 3, 5, 7, 10, 100, -1, -5]


def primes_in(a, b):
    return [p for p in range(max(a, 2), b + 1) if all(p % d for d in range(2, isqrt(p) + 1))]


def vp(x, p):
    if x == 0:
        return 99
    v = 0
    while x % p == 0:
        x //= p
        v += 1
    return v


def fq(m, p):
    pp = p * p
    return ((pow(m % pp, p - 1, pp) - 1) // p) % p


HS = [Fraction(0)]
for r in range(1, 20):
    HS.append(HS[-1] + Fraction(1, r))


def weighted_terms(Avals, wts, p):
    """[A(k) w_k] mod p for k < p, plus a flag list of non-integral k."""
    H = [0] * p
    for i in range(1, p):
        H[i] = (H[i - 1] + pow(i, -1, p)) % p
    out = []
    nonint = []
    for k in range(p):
        A = Avals[k]
        v = vp(A, p)
        WR = sum(c * H[(s * k) % p] for c, s in wts) % p
        WS = sum(c * HS[(s * k) // p] for c, s in wts)            # exact, tiny denominators
        if v == 0:
            if WS != 0:
                nonint.append(k)
            out.append(A % p * WR % p)
        elif v == 1:
            out.append((A // p) % p * (WS.numerator * pow(WS.denominator, -1, p)) % p)
        else:
            out.append(0)
    return out, nonint


def main():
    PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    PO = int(sys.argv[2]) if len(sys.argv) > 2 else 300
    primes = primes_in(7, PMAX)
    Kmax = 2 * max(PMAX, PO) + 2
    for fam, (Af, wts) in FAM.items():
        Avals = [Af(k) for k in range(Kmax)]
        # (1) the law at the CM parameters
        fails = {m: [] for m in CMP[fam]}
        tested = {m: 0 for m in CMP[fam]}
        nonint_any = []
        shift_bad = 0
        shift_tot = 0
        tt_rows = 0
        tt_dis = 0
        for p in primes:
            pp = p * p
            AW, nonint = weighted_terms(Avals, wts, p)
            if nonint:
                nonint_any.append((p, nonint[:4]))
            # (2) p-shift with carries, a = 1
            if p <= 61:
                for j in range(p):
                    shift_tot += 1
                    lhs = Avals[p + j] % pp
                    rhs = Avals[1] * (Avals[j] + p * AW[j]) % pp
                    shift_bad += (lhs - rhs) % pp != 0
            for m in CMP[fam] + (GEN if p <= 200 else []):
                if (6 * m) % p == 0:
                    continue
                w = pow(m % p, -1, p)
                pw = 1
                F = 0
                W = 0
                for k in range(p):
                    F = (F + Avals[k] % p * pw) % p
                    W = (W + AW[k] * pw) % p
                    pw = pw * w % p
                law = (W - fq(m, p) * F) % p == 0
                if m in CMP[fam]:
                    tested[m] += 1
                    if not law:
                        fails[m].append(p)
                if p <= 200:
                    # (B1') mod p^2
                    w2 = pow(m % pp, -1, pp)
                    pw2 = 1
                    S0 = S1 = 0
                    for k in range(2 * p):
                        t = Avals[k] % pp * pw2 % pp
                        if k < p:
                            S0 += t
                        else:
                            S1 += t
                        pw2 = pw2 * w2 % pp
                    B1 = (S1 - Avals[1] * w2 * S0) % pp == 0
                    tt_rows += 1
                    tt_dis += (B1 != law)
        nf = sum(len(v) for v in fails.values())
        print('%s: (1) law at %d CM parameters, primes 7..%d with p !| 6m: %d prime-instances, %d failures %s'
              % (fam, len(CMP[fam]), PMAX, sum(tested.values()), nf, {m: v for m, v in fails.items() if v}))
        print('    non-p-integral terms A(k) w_k (k < p): %s' % (nonint_any[:5] if nonint_any else 'none'))
        print('    (2) A(p+j) == A(1)(A(j) + p[A(j)w_j]) mod p^2 for all 0 <= j < p (carries included), p <= 61:'
              ' %d checks, %d failures' % (shift_tot, shift_bad))
        print('    (3) (B1\') <=> (L_p^A) truth table, CM + 8 generic m, p <= 200: %d rows, %d disagreements'
              % (tt_rows, tt_dis), flush=True)
        # (4) part (O)
        ordi = ordi_ok = ss = ss_ok = 0
        sporadic = []
        for p in primes_in(7, PO):
            P3 = p ** 3
            for m in CMP[fam]:
                if (6 * m) % p == 0:
                    continue
                w = pow(m % P3, -1, P3)
                pw = 1
                S0 = S1 = 0
                for k in range(2 * p):
                    t = Avals[k] % P3 * pw % P3
                    if k < p:
                        S0 += t
                    else:
                        S1 += t
                    pw = pw * w % P3
                ok3 = (S1 - Avals[1] * w * S0) % P3 == 0
                if S0 % p:
                    ordi += 1
                    ordi_ok += ok3
                else:
                    ss += 1
                    ss_ok += ok3
                    if ok3:
                        sporadic.append((m, p))
        print('    (4) (B1\') mod p^3, primes 7..%d: ordinary %d of %d hold; supersingular %d of %d hold %s'
              % (PO, ordi_ok, ordi, ss_ok, ss, sporadic[:10]), flush=True)


if __name__ == '__main__':
    main()
