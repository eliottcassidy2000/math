#!/usr/bin/env python3
"""S39 audit, item 4 (finite-field part) and the coordinator's addendum on THM-4663 (10).

All exact (integers, Legendre symbols); no floating point.
  A. R_k(c) (k <= 3) by an exact radial walk DP (model th^2 = th + r) for every prime 5 <= p <= 199, against
     #E'_c(F_p) and #L_c(F_p) (character sums, smooth models with two points at infinity), every c != 0, 1, 9;
     R_3(1) = 3p, R_3(9) = p + 1 - (p/3); 6 | R_3(c); the theta identities of statement 8.
  B. the four-bar moments for every prime 5 <= p < 700: sum a_c, sum a_c^2, sum eta_c, sum eta_c a_c
     (eta_c = eta(c(c - 4))), against THM-4663 (10) with a_p(f15), b_p from a2 (independent).
  C. the K3 reading: the degree-2 base change X along s^2 = c(c - 4) (c = 4/(1 - m^2)); the trace of Frobenius on
     H^1(P^1_m, j_* R^1) = -sum_(m smooth) a_(c(m)) - sum_(m multiplicative) eps_m, with eps from the Weierstrass model
     Y^2 = X^3 + (c+3)^2 X^2 + 8(c-1)^2(c+3) X + 16(c-1)^4 (eps_0 = eps_9 = (p/3), eps_1 = eps_oo = 1), compared with
     (p/3) p + a_p(f15) and with (p/3) p - a_p(f15): this decides the sign of the transcendental trace.
  D. the twisted second moment and the third moment: Lefschetz predictions
       sum_c eta_c a_c^2 = -(p + 1)(1 + (p/3) + (5/p)) - T_p,   T_p = trace on a 6-dim weight-3 H^1,
       sum_c a_c^3      = -(1 + (p/3))(4p + 2) - t_p,          t_p = trace on a 4-dim weight-4 H^1,
     checked through the Weil bounds |T_p| <= 6 p^(3/2), |t_p| <= 4 p^2; t_p mod p recorded.
  E. smoothness of the plane cubic of statement 8 modulo every prime 7 <= p < 200 (Groebner basis over F_p).
Usage: python -B a4b_fourbar_finite_fields.py
"""
import json
import numpy as np
import sympy
from sympy import primerange

F = json.load(open('a2_forms_to_3000.json'))
A15, B46 = F['f15'], F['b46']


def qr_table(p):
    t = -np.ones(p, dtype=np.int64)
    t[0] = 0
    t[(np.arange(1, p) ** 2) % p] = 1
    return t


def kr(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def counts_both(p, chi):
    """#E'_c(F_p) and #L_c(F_p) for every c in F_p (arrays)."""
    tt = np.arange(p, dtype=np.int64)
    E = np.zeros(p, dtype=np.int64); L = np.zeros(p, dtype=np.int64)
    for c in range(p):
        QE = ((tt * tt - 4 * c) % p) * ((tt - c - 1) % p) % p * ((tt - c + 3) % p) % p
        A = (tt * tt + (1 + c) * tt + 1) % p
        QL = (A * A - 4 * c * tt % p * ((tt + 1) * (tt + 1) % p)) % p
        E[c] = p + 2 + chi[QE].sum()
        L[c] = p + 2 + chi[QL].sum()
    return E, L


def radial(p):
    r = next(r for r in range(1, p) if pow((1 + 4 * r) % p, (p - 1) // 2, p) == p - 1)
    Nf = lambda a, b: (a * a + a * b - r * b * b) % p
    mu = [(a, b) for a in range(p) for b in range(p) if Nf(a, b) == 1]
    rep = {}
    for a in range(p):
        for b in range(p):
            rep.setdefault(Nf(a, b), (a, b))
    T = [[Nf(rep[c][0] - x, rep[c][1] - y) for (x, y) in mu] for c in range(p)]
    R = {0: [1] + [0] * (p - 1)}
    for k in range(1, 4):
        R[k] = [sum(R[k - 1][n] for n in T[c]) for c in range(p)]
    return R


# ------------------------------------------------------------------ A
badA, bad16, badth, bad6 = [], [], [], []
for p in primerange(5, 200):
    chi = qr_table(p)
    E, L = counts_both(p, chi)
    R = radial(p)
    e3 = kr(p, 3)
    for c in range(1, p):
        if c in (1, 9 % p):
            continue
        if not (R[3][c] == E[c] == L[c]):
            badA.append((p, c))
        if R[3][c] % 6:
            bad6.append((p, c))
    if R[3][1] != 3 * p or R[3][9 % p] != p + 1 - e3:
        bad16.append(p)
    th = lambda x, y, z: R[x][0] * R[y][0] * R[z][0] + (p + 1) * sum(R[x][c] * R[y][c] * R[z][c] for c in range(1, p))
    W3, W5 = th(1, 1, 1) if False else R[3][0], None
    W5 = R[2][0] * R[3][0] + (p + 1) * sum(R[2][c] * R[3][c] for c in range(1, p))
    chk = [th(2, 2, 2) == p ** 3 + 7 * p * p + 4 * p - 2, th(1, 2, 3) == 3 * (1 - e3) * p * (p + 1),
           th(1, 3, 3) == 9 * p * p * (p + 1), th(2, 2, 3) == 2 * W5 + (p * p - 1) * W3 - (p + 1) * R[3][4 % p]]
    if not all(chk):
        badth.append((p, chk))
print('A. exact radial DP, every prime 5 <= p < 200:')
print('   R_3(c) = #E\'_c(F_p) = #L_c(F_p) for every c not in {0,1,9}: failures', badA[:6])
print('   6 | R_3(c): failures', bad6[:6], '; R_3(1) = 3p and R_3(9) = p + 1 - (p/3): failures', bad16)
print('   theta(2,2,2), theta(1,2,3), theta(1,3,3) closed forms and the theta(2,2,3) identity: failures', badth)
assert not badA and not bad6 and not bad16 and not badth

# ------------------------------------------------------------------ B, C, D
badB, badK3plus, badK3minus = [], [], []
rowsD = []
for p in primerange(5, 700):
    chi = qr_table(p)
    E, L = counts_both(p, chi)
    a = p + 1 - E
    e3, e5 = kr(p, 3), kr(5, p)
    cs = [c for c in range(1, p) if c not in (1, 9 % p)]
    eta = {c: kr(c * (c - 4), p) for c in range(p)}
    S1 = sum(a[c] for c in cs); S2 = sum(a[c] ** 2 for c in cs)
    T0 = sum(eta[c] for c in cs); T1 = sum(eta[c] * a[c] for c in cs)
    ok = (S1 == -2 - 2 * e3, S2 == (p + 1) * (p - 4) - B46[p], T0 == -1 - e3 - e5,
          T1 == -1 - e3 - kr(-15, p) - e3 * p - A15[p])
    if not all(ok):
        badB.append((p, ok))
    if p >= 11:
        # C: base change along s^2 = c(c - 4), c = 4/(1 - m^2); m = oo -> c = 0; m = +-1 -> c = oo
        tr = 0
        eps = {0: e3, 1: 1, 9 % p: e3}
        for m in range(p):
            den = (1 - m * m) % p
            if den == 0:
                tr -= 1                     # c = oo: split I_3
                continue
            cm = 4 * pow(den, -1, p) % p
            if cm in eps:
                tr -= eps[cm]
            else:
                tr -= int(a[cm])
        tr -= eps[0]                        # m = oo: c = 0, I_2 after base change, eps = (p/3)
        if tr != e3 * p + A15[p]:
            badK3plus.append(p)
        if tr != e3 * p - A15[p]:
            badK3minus.append(p)
        # D
        T2 = sum(eta[c] * a[c] ** 2 for c in cs)
        S3 = sum(a[c] ** 3 for c in cs)
        Tp = -T2 - (p + 1) * (1 + e3 + e5)
        tp = -S3 - (1 + e3) * (4 * p + 2)
        rowsD.append((p, Tp, Tp / p ** 1.5, tp, tp / p ** 2, tp % p))
print('\nB. four-bar moments, every prime 5 <= p < 700 (sum over c not in {0,1,9}):')
print('   sum a_c = -2-2(p/3), sum a_c^2 = (p+1)(p-4) - b_p, sum eta_c = -1-(p/3)-(5/p), sum eta_c a_c = -1-(p/3)-(-15/p)-(p/3)p-a_p(f15):')
print('   failures', badB)
print('\nC. K3 base change: trace on H^1(P^1_m, j_*R^1) = (p/3) p + a_p(f15) fails at', badK3plus[:10],
      '; = (p/3) p - a_p(f15) fails at', len(badK3minus), 'primes (first', badK3minus[:6], ')')
mxT = max(abs(r[2]) for r in rowsD); mxt = max(abs(r[4]) for r in rowsD)
print('\nD. Weil-bound check of the Lefschetz predictions, primes 11 <= p < 700:')
print('   max |T_p|/p^(3/2) = %.3f (bound 6: twisted second moment = 6-dim weight-3 motive)' % mxT)
print('   max |t_p|/p^2     = %.3f (bound 4: third moment = 4-dim weight-4 motive)' % mxt)
print('   t_p mod p = 0 at %d of %d primes' % (sum(1 for r in rowsD if r[5] == 0), len(rowsD)))
print('   samples (p, T_p, T_p/p^1.5, t_p, t_p/p^2):', [(r[0], r[1], round(r[2], 3), r[3], round(r[4], 3)) for r in rowsD[:8]])
assert not badB and not badK3plus and mxT <= 6 and mxt <= 4

# ------------------------------------------------------------------ E
u, v = sympy.symbols('u v')
Fc = 2 * u ** 2 * v + 2 * u * v ** 2 - u ** 2 - 5 * u * v - v ** 2 + 2 * u + 2 * v
sing = []
for p in primerange(2, 200):
    G = sympy.groebner([Fc, sympy.diff(Fc, u), sympy.diff(Fc, v)], u, v, modulus=p, order='lex')
    if list(G) != [1]:
        sing.append(p)
print('\nE. the plane cubic has an affine singular point over Fbar_p exactly for p in', sing, '(p < 200)')
print('ALL CHECKS PASSED')
