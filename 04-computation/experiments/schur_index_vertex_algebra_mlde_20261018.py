#!/usr/bin/env python3
"""4d N=2 SCFTs -> vertex algebras (Beem-Lemos-Liendo-Peelaers-Rastelli-van Rees): Schur indices and modular linear
differential equations, computed from scratch (S43).

1. SU(2) with N_f = 4 hypermultiplets (flavour so(8)). The Schur index is the SU(2) Haar integral
       I(q) = (1/2) oint (1-z^2)(1-z^-2) prod_n (1-q^n)^2 (1-z^2 q^n)^2 (1-z^-2 q^n)^2
                         / prod_(n>=0) [(1-z q^(n+1/2))(1-z^-1 q^(n+1/2))]^8   dz/(2 pi i z),
   expanded exactly in q^(1/2). Prediction (BLLPRR; Beem-Rastelli 2017; Arakawa-Kawasetsu 2018): I is the vacuum
   character of so(8) at level -2, c_2d = -14, and f = q^(7/12) I satisfies the second-order MLDE
       D^(2) f = alpha E_4 f,   alpha = h (h - 1/6),  h = -c/24 = 7/12,  alpha = 35/144,
   with D^(2) = (q d/dq - E_2/6)(q d/dq).
2. The (A1, A2) Argyres-Douglas theory: Schur index = q^(11/60) H(q), H = prod 1/((1-q^(5n+2))(1-q^(5n+3))), the
   vacuum character of the Lee-Yang model M(2,5), c = -22/5 (corrected after S43 audit n1: G(q) =
   prod 1/((1-q^(5n+1))(1-q^(5n+4))) is the h = -1/5 character); check the MLDE with alpha = 11/3600.
3. Level-5 link to the repository's Apery zeta(2) numbers: Beukers' Hauptmodul t(tau) = q prod (1-q^n)^(5 (n/5))
   (the fifth power of the Rogers-Ramanujan continued fraction); we expand sum b_n t^n with b_n = sum C(n,k)^2 C(n+k,k)
   and test whether it is a weight-1 form on Gamma_1(5): an eigenform of the diamond/Hecke structure is checked through
   multiplicativity of its q-coefficients.
Usage: python schur_index_vertex_algebra_mlde_20261018.py [NQ]
Session opus-2026-10-17-S43.
"""
import sys, os, json
from fractions import Fraction as Fr
from collections import defaultdict


def schur_su2_nf4(NH):
    """Return I(q) coefficients in powers of q^(1/2) up to q^(NH/2) (NH = number of half-steps)."""
    # series as dict: (half-power h, z-power) -> int; truncated at h <= NH, |z| <= ZM
    ZM = NH + 4
    def mul(A, B):
        C = defaultdict(int)
        for (h1, z1), c1 in A.items():
            for (h2, z2), c2 in B.items():
                h = h1 + h2
                if h <= NH and abs(z1 + z2) <= ZM:
                    C[(h, z1 + z2)] += c1 * c2
        return {k: v for k, v in C.items() if v}
    S = {(0, 0): 1}
    # vector multiplet: prod_{n>=1} (1-q^n)^2 (1-z^2 q^n)^2 (1-z^-2 q^n)^2
    for n in range(1, NH // 2 + 1):
        for zz in (0, 2, -2):
            f = {(0, 0): 1, (2 * n, zz): -1}
            S = mul(S, f)
            S = mul(S, f)
    # hypers: prod_{n>=0} 1/((1 - z q^(n+1/2))(1 - z^-1 q^(n+1/2)))^8 ; 1/(1-x) = sum x^m
    for n in range(0, NH):
        h = 2 * n + 1
        if h > NH:
            break
        for zz in (1, -1):
            geo = {(m * h, m * zz): 1 for m in range(0, NH // h + 1)}
            for _ in range(8):
                S = mul(S, geo)
    # Haar measure (1/2)(1-z^2)(1-z^-2) = (1/2)(2 - z^2 - z^-2)
    meas = {(0, 0): 2, (0, 2): -1, (0, -2): -1}
    S = mul(S, meas)
    return [Fr(S.get((h, 0), 0), 2) for h in range(NH + 1)]


def sigma(n, k=1):
    return sum(d ** k for d in range(1, n + 1) if n % d == 0)


def E2(N):
    return [Fr(1)] + [Fr(-24 * sigma(n)) for n in range(1, N + 1)]


def E4(N):
    return [Fr(1)] + [Fr(240 * sigma(n, 3)) for n in range(1, N + 1)]


def mlde_residual(coeffs, h, alpha, N):
    """f = q^h sum c_n q^n. Residual coefficients of D^(2) f - alpha E_4 f (divided by q^h)."""
    e2, e4 = E2(N), E4(N)
    c = coeffs
    qf = [(h + n) * c[n] for n in range(N + 1)]                 # q d/dq f
    qqf = [(h + n) ** 2 * c[n] for n in range(N + 1)]           # (q d/dq)^2 f
    res = []
    for n in range(N + 1):
        t = qqf[n] - Fr(1, 6) * sum(e2[k] * qf[n - k] for k in range(n + 1)) - alpha * sum(e4[k] * c[n - k] for k in range(n + 1))
        res.append(t)
    return res


def main():
    NQ = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    out = {}
    I = schur_su2_nf4(2 * NQ)
    half = [I[h] for h in range(1, 2 * NQ + 1, 2)]
    integral = [I[h] for h in range(0, 2 * NQ + 1, 2)]
    print('1. SU(2) N_f = 4 Schur index: half-integer powers vanish: %s' % all(x == 0 for x in half))
    print('   I(q) = %s' % ' + '.join('%s q^%d' % (integral[n], n) for n in range(min(8, NQ + 1))))
    h = Fr(7, 12)
    alpha = h * (h - Fr(1, 6))
    res = mlde_residual(integral, h, alpha, NQ)
    print('   MLDE D^(2) f = (35/144) E_4 f with f = q^(7/12) I: residuals n = 0..%d all zero: %s; alpha = %s' % (NQ, all(r == 0 for r in res), alpha))
    out['su2_nf4'] = dict(I=[str(x) for x in integral], mlde_ok=all(r == 0 for r in res), alpha=str(alpha))
    # 2. Rogers-Ramanujan / Lee-Yang M(2,5): vacuum character q^(11/60) H(q), H = prod 1/((1-q^(5n+2))(1-q^(5n+3)))
    #    (no q^1 term: L_(-1)|0> = 0), and the h = -1/5 character q^(-1/60) G(q); both solve one MLDE, alpha = 11/3600.
    N = 60
    def rr(rs):
        F = [Fr(0)] * (N + 1)
        F[0] = Fr(1)
        for n in range(N + 1):
            for r in rs:
                k = 5 * n + r
                if k > N:
                    continue
                for m in range(k, N + 1):
                    F[m] += F[m - k]
        return F
    H = rr((2, 3)); G = rr((1, 4))
    a2 = Fr(11, 60) * (Fr(11, 60) - Fr(1, 6))
    okH = all(r == 0 for r in mlde_residual(H, Fr(11, 60), a2, N))
    okG = all(r == 0 for r in mlde_residual(G, Fr(-1, 60), a2, N))
    print('2. Lee-Yang: vacuum H(q) = %s ...; MLDE (alpha = %s) residuals zero to q^%d: vacuum q^(11/60)H %s, primary q^(-1/60)G %s'
          % ([str(x) for x in H[:12]], a2, N, okH, okG))
    out['lee_yang'] = dict(H=[str(x) for x in H[:30]], G=[str(x) for x in G[:30]], alpha=str(a2), mlde_vacuum=okH, mlde_primary=okG)
    # 3. Beukers level 5: t = q prod (1-q^n)^(5 (n/5)), sum b_n t^n
    def leg5(n):
        r = n % 5
        return 0 if r == 0 else (1 if r in (1, 4) else -1)
    M = 50
    t = [Fr(0)] * (M + 1)
    base = [Fr(0)] * (M + 1)
    base[0] = Fr(1)
    for n in range(1, M + 1):
        e = 5 * leg5(n)
        for _ in range(abs(e)):
            if e > 0:
                for m in range(M, n - 1, -1):
                    base[m] -= base[m - n]
            else:
                for m in range(n, M + 1):
                    base[m] += base[m - n]
    for m in range(M):
        t[m + 1] = base[m]
    from math import comb
    b = [sum(comb(n, k) ** 2 * comb(n + k, k) for k in range(n + 1)) for n in range(M + 1)]
    # sum b_n t^n
    tot = [Fr(0)] * (M + 1)
    tp = [Fr(0)] * (M + 1)
    tp[0] = Fr(1)
    for n in range(M + 1):
        for i in range(M + 1):
            tot[i] += b[n] * tp[i]
        # tp *= t
        new = [Fr(0)] * (M + 1)
        for i, x in enumerate(tp):
            if x:
                for j, yv in enumerate(t[:M + 1 - i]):
                    if yv:
                        new[i + j] += x * yv
        tp = new
    print('3. Beukers level 5: sum b_n t^n = %s' % [str(x) for x in tot[:16]])
    # Beukers: sum b_n t^n = 1 + sum_n Re[(3 - i) sum_{d|n} chi(d)] q^n, chi the character mod 5 with chi(2) = i
    chi = {0: 0, 1: 1, 2: 1j, 3: -1j, 4: -1}
    ok = True
    for n in range(1, M + 1):
        S = sum(chi[d % 5] for d in range(1, n + 1) if n % d == 0)
        val = ((3 - 1j) * S).real
        if abs(val - float(tot[n])) > 1e-9:
            ok = False
    print('   matches the weight-1 Eisenstein form 1 + sum Re[(3-i) sum_{d|n} chi(d)] q^n on Gamma_1(5) to q^%d: %s' % (M, ok))
    out['beukers_level5'] = dict(series=[str(x) for x in tot[:30]], eisenstein_match=ok)
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(out, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'schur_index_vertex_algebra_mlde_20261018.json'), 'w'),
              indent=1)


if __name__ == '__main__':
    main()
