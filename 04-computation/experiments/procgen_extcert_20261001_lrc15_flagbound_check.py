#!/usr/bin/env python3
"""procgen_extcert_20261001_lrc15_flagbound_check.py -- the orchestrator's independent check of the analytic half of
J. Allikvere, "Fifteen lonely runners" (Zenodo record 22667683; part of arXiv:2609.02604 "Fourteen and fifteen lonely
runners"): the lattice flag bound (Section 3), which lowers the finite-checking threshold for 14 speeds from 810.07 to
414.78, the forced divisor, and the prime mass of the 71 gates.

Written from the statements and proofs in paper_v2.tex (read by the orchestrator); none of the package's code is read
or run.  Only data is taken from the package: the 71 gate primes (SHA256SUMS.txt file names) and, for comparison, the
numbers printed in the paper / LRC15_AUDIT_SUMMARY.json.

Checks:
  A. Constants: B_r = (2r/((n+1)(n+1-r)))^2, c_i = B_i - B_{i-1}, K_14 = prod c_i equals the paper's exact rational;
     min c_{i+1}/c_i = 847/513 > 4/3 (also n = 13: 1175/688, n = 15: 3751/2349); A_14 = K_14^(-1/2) = 206083383792625.44..;
     14 log(A_14/28) = 414.77936455669981..; target 14 log(A_14/28) - log 360360 = 401.98450574593444..;
     lcm(2..15) = 360360; the MSS threshold 14(13 log 105 - log 14) = 810.07..; the n = 13 value 341.03...
  B. Lemma 3.6 (product minimization) for n = 14, d = 13, rho = 3/4, by brute force: every one of the C(25,13) = 5,200,300
     choices of 13 active constraints among the 12 ratio constraints t_{i+1} >= rho t_i and the 13 prefix constraints
     t_1 + .. + t_r >= B_r is solved; over all feasible vertices, sum log t_i >= sum log c_i, with equality only at t = c.
     (A concave nondecreasing objective on a pointed polyhedron whose recession cone lies in the nonnegative orthant is
     minimized at a vertex.)
  C. Lemma 3.2 (covolume in the ellipsoidal norm) and Lemma 3.1 (E inside D) on random primitive speed vectors
     (n = 3..7): the E-Gram determinant of an explicit basis of P Z^n equals 1/(S F_n(x))^2; the support functions satisfy
     h_E <= h_D at every vertex of the projected cross-polytope and in random directions.
  D. Lemma 3.7 (shape bound, n = 14): H(u) = (sum w(u_i)) prod a(u_i) > 784 on (0,1]^14 with max u_i = 1, by multistart
     L-BFGS-B (global minimum 798.31 at (1, t0, .., t0), t0 = 0.078207); the paper's one-variable bound H(t0) > 796.4; and
     R_14(x) < 1/2 directly on random points of the simplex (checks the reduction R^2 = n^2/H).
  E. The prime mass: the 71 gate primes are distinct primes > 15 (so coprime to 360360); sum log p = 408.8233173785861..
     > target; margin 6.8388 > log 569; the 52 'J(14,p) empty' gates alone give less than the target, so the 19 gates
     closed with the few-exception shift lemma (p = 89..241 except 239) are load-bearing.
"""
import itertools
import math
import os
import random
import sys
import time
from fractions import Fraction

import mpmath
import numpy as np
from scipy.optimize import minimize

mpmath.mp.dps = 60
HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.join(HERE, '..', '..', 'scratch', 'lrc_certs', 'lrc15')
OKS = []


def ok(c, msg):
    OKS.append(bool(c))
    print(('[OK] ' if c else '[FAIL] ') + msg, flush=True)


def flag_constants(n):
    d = n - 1
    B = [Fraction(0)] + [Fraction(2 * r, (n + 1) * (n + 1 - r)) ** 2 for r in range(1, d + 1)]
    c = [B[i] - B[i - 1] for i in range(1, d + 1)]
    K = Fraction(1)
    for x in c:
        K *= x
    ratios = [c[i + 1] / c[i] for i in range(d - 1)]
    return B, c, K, min(ratios)


def check_A():
    B, c, K, rmin = flag_constants(14)
    Kpaper = Fraction(2289670297, 97243124257250844122250000000000000000)
    A = mpmath.mpf(K.denominator) ** 0.5 / mpmath.mpf(K.numerator) ** 0.5
    thr = 14 * mpmath.log(A / 28)
    target = thr - mpmath.log(360360)
    ok(K == Kpaper, f'A: K_14 = prod c_i = {K} (= the paper\'s exact rational)')
    ok(rmin == Fraction(847, 513) and rmin > Fraction(4, 3), f'A: min c_(i+1)/c_i = {rmin} > 4/3 (n = 14)')
    ok(abs(A - mpmath.mpf('206083383792625.44')) < 0.01, f'A: A_14 = K_14^(-1/2) = {mpmath.nstr(A, 20)}')
    ok(abs(thr - mpmath.mpf('414.77936455669981')) < 1e-13, f'A: 14 log(A_14/28) = {mpmath.nstr(thr, 20)} (paper: 414.77936455669981...)')
    ok(abs(target - mpmath.mpf('401.9845057459344408819653635544')) < 1e-25,
       f'A: target = 14 log(A_14/28) - log 360360 = {mpmath.nstr(target, 25)} (LRC15_AUDIT_SUMMARY: 401.98450574593444088...)')
    ok(math.lcm(*range(2, 16)) == 360360 == 2 ** 3 * 3 ** 2 * 5 * 7 * 11 * 13, 'A: lcm(2..15) = 360360 = 2^3 3^2 5 7 11 13 (forced divisor, Lemma 3.9)')
    mss = 14 * (13 * mpmath.log(105) - mpmath.log(14))
    ok(abs(mss - mpmath.mpf('810.07')) < 0.01, f'A: the MSS threshold 14(13 log 105 - log 14) = {mpmath.nstr(mss, 8)}')
    _, _, K13, r13 = flag_constants(13)
    A13 = mpmath.mpf(K13.denominator) ** 0.5 / mpmath.mpf(K13.numerator) ** 0.5
    t13 = 13 * (mpmath.log(A13) - mpmath.log(13) + mpmath.log(mpmath.mpf(101) / 200))
    _, _, K15, r15 = flag_constants(15)
    A15 = mpmath.mpf(K15.denominator) ** 0.5 / mpmath.mpf(K15.numerator) ** 0.5
    t15 = 15 * mpmath.log(A15 / 30)
    ok(r13 == Fraction(1175, 688) and r15 == Fraction(3751, 2349) and abs(t13 - mpmath.mpf('341.03')) < 0.01
       and abs(t15 - mpmath.mpf('497.03')) < 0.01,
       f'A: n = 13: ratio min {r13}, threshold {mpmath.nstr(t13, 8)}; n = 15: ratio min {r15}, threshold {mpmath.nstr(t15, 8)} (Table 1)')
    return B, c


def check_B(B, c):
    d = 13
    rho = 0.75
    rows, rhs = [], []
    for i in range(d - 1):          # t_{i+1} - rho t_i >= 0
        a = np.zeros(d); a[i] = -rho; a[i + 1] = 1.0
        rows.append(a); rhs.append(0.0)
    for r in range(1, d + 1):       # t_1 + .. + t_r >= B_r
        a = np.zeros(d); a[:r] = 1.0
        rows.append(a); rhs.append(float(B[r]))
    Afull = np.array(rows); bfull = np.array(rhs)
    target = sum(math.log(float(x)) for x in c)
    cvec = np.array([float(x) for x in c])
    combos = itertools.combinations(range(25), 13)
    best = math.inf; best_t = None; nfeas = 0; nsing = 0; total = 0; near = 0
    CH = 40000
    while True:
        chunk = np.array(list(itertools.islice(combos, CH)), dtype=np.int64)
        if len(chunk) == 0:
            break
        total += len(chunk)
        M = Afull[chunk]; b = bfull[chunk]
        det = np.linalg.det(M)
        good = np.abs(det) > 1e-14
        nsing += int((~good).sum())
        M = M[good]; b = b[good]
        t = np.linalg.solve(M, b[..., None])[..., 0]
        slack = t @ Afull.T - bfull
        feas = (slack >= -1e-12 * (1 + np.abs(bfull))).all(axis=1) & (t > 0).all(axis=1)
        nfeas += int(feas.sum())
        if feas.any():
            tf = t[feas]
            val = np.log(tf).sum(axis=1)
            j = int(np.argmin(val))
            if val[j] < best:
                best, best_t = float(val[j]), tf[j]
            near += int((val < target + 1e-9).sum())
    at_c = best_t is not None and np.allclose(best_t, cvec, rtol=1e-9)
    ok(total == math.comb(25, 13) and best >= target - 1e-9 and at_c,
       f'B: Lemma 3.6 (n = 14): {total} active sets ({nsing} singular, {nfeas} feasible vertices); min sum log t = {best:.12f} >= sum log c_i = {target:.12f}, attained at t = c ({near} vertex solutions within 1e-9 of it)')


def unimodular_with_first_column(v):
    """integer matrix W with det +-1 and first column v (v primitive), by recording Euclid row operations"""
    n = len(v)
    w = list(v)
    M = [[int(i == j) for j in range(n)] for i in range(n)]      # M v = w throughout
    Minv = [[int(i == j) for j in range(n)] for i in range(n)]   # Minv w = v throughout (Minv = M^-1)
    def rowop(i, j, q):  # row_i -= q row_j  (w_i -= q w_j)
        w[i] -= q * w[j]
        M[i] = [a - q * b for a, b in zip(M[i], M[j])]
        for row in Minv:   # column_j += q column_i
            row[j] += q * row[i]
    def swap(i, j):
        w[i], w[j] = w[j], w[i]
        M[i], M[j] = M[j], M[i]
        for row in Minv:
            row[i], row[j] = row[j], row[i]
    while any(w[i] != 0 for i in range(1, n)):
        nz = [i for i in range(n) if w[i] != 0]
        i0 = min(nz, key=lambda i: abs(w[i]))
        if i0 != 0:
            swap(0, i0)
        for i in range(1, n):
            if w[i] != 0:
                rowop(i, 0, w[i] // w[0])
    if w[0] < 0:
        w[0] = -w[0]; M[0] = [-a for a in M[0]]
        for row in Minv:
            row[0] = -row[0]
    assert w[0] == 1, 'v not primitive'
    W = np.array(Minv, dtype=float)
    assert all(int(round(x)) == y for x, y in zip(W[:, 0], v))
    return W


def check_C():
    rng = random.Random(15)
    worst_cov = 0.0; worst_inc = -1.0; ntest = 0
    for trial in range(300):
        n = rng.randint(3, 7)
        while True:
            v = [rng.randint(1, 40) for _ in range(n)]
            if math.gcd(*v) == 1:
                break
        vv = np.array(v, dtype=float)
        S = vv.sum(); x = vv / S; u = vv / vv.max(); q = 1 + 2 * u - u * u
        F = math.sqrt(np.prod(q) * np.sum(x * x / q))
        P = np.eye(n) - np.outer(vv, vv) / vv.dot(vv)
        U = np.linalg.svd(P)[0][:, :n - 1]                       # orthonormal basis of H = v^perp
        W = unimodular_with_first_column(v)
        basis = P @ W[:, 1:]                                     # P w_2 .. P w_n: a basis of Lambda = P Z^n
        Y = U.T @ basis
        GE = np.linalg.inv(U.T @ np.diag(q) @ U)                 # E = P Q^(1/2) B: Gram matrix (U^T Q U)^(-1)
        cov = math.sqrt(np.linalg.det(Y.T @ GE @ Y))
        worst_cov = max(worst_cov, abs(cov * S * F - 1))
        # inclusion: h_E(y) = |Q^(1/2) y| <= h_D(y) = sum |y_i| for y in H
        ys = [(vv[j] * np.eye(n)[i] - vv[i] * np.eye(n)[j]) / (vv[i] + vv[j]) for i in range(n) for j in range(n) if i != j]
        ys += [P @ np.array([rng.gauss(0, 1) for _ in range(n)]) for _ in range(50)]
        for y in ys:
            worst_inc = max(worst_inc, math.sqrt(np.sum(q * y * y)) - np.sum(np.abs(y)))
        ntest += 1
    ok(worst_cov < 1e-9, f'C: Lemma 3.2: covol_E(P Z^n) = 1/(S F_n(x)) on {ntest} random primitive v (n = 3..7), max relative error {worst_cov:.1e}')
    ok(worst_inc <= 1e-12, f'C: Lemma 3.1: h_E <= h_D at every vertex of the projected cross-polytope and in 50 random directions per v (max h_E - h_D = {worst_inc:.2e})')


def H_of_u(u, n=14):
    q = 1 + 2 * u - u * u
    return np.sum(u * u / q) * np.prod(q / u ** (2.0 / n))


def check_D():
    n = 14
    rng = np.random.default_rng(714)
    def f(z):
        u = np.concatenate([np.clip(z, 1e-12, 1.0), [1.0]])
        return math.log(H_of_u(u))
    best = math.inf; bestu = None
    for s in range(400):
        z0 = rng.uniform(0.01, 1.0, n - 1) if s % 2 else rng.uniform(0.01, 0.3, n - 1)
        res = minimize(f, z0, method='L-BFGS-B', bounds=[(1e-9, 1.0)] * (n - 1))
        if res.fun < best:
            best, bestu = res.fun, res.x
    Hmin = math.exp(best)
    t = np.sort(bestu)
    # the paper's one-variable reduction and its bound H(t0) > 796.4
    Ht = lambda t: (1 + 2 * t + 25 * t * t) * (1 + 2 * t - t * t) ** 12 / t ** (13 / 7)
    ts = np.linspace(0.05, 0.1, 200001)
    vals = Ht(ts)
    j = int(np.argmin(vals))
    root = [r.real for r in np.roots([325, 23, 9, -1]) if abs(r.imag) < 1e-12 and r.real > 0]
    # exact rational comparison of the paper's displayed lower bound, raised to the 7th power
    num = (1 + 2 * Fraction(782, 10000) + 25 * Fraction(782, 10000) ** 2) * (1 + 2 * Fraction(782, 10000) - Fraction(782, 10000) ** 2) ** 12
    # num / 0.0783^(13/7) > 796.4  <=>  num^7 > 796.4^7 * 0.0783^13
    exact = num ** 7 > Fraction(7964, 10) ** 7 * Fraction(783, 10000) ** 13
    ok(Hmin > 784 and abs(Hmin - vals[j]) < 1e-3 and np.allclose(t, root[0], atol=1e-4),
       f'D: Lemma 3.7: multistart minimum of H on (0,1]^14 (max u = 1) = {Hmin:.4f} > 784 = 4 n^2, at u = (1, {t.min():.5f} .. {t.max():.5f}) (= the one-variable minimum below)')
    ok(0.0782 < root[0] < 0.0783 and abs(ts[j] - root[0]) < 1e-5 and vals[j] > 796.4 and exact,
       f'D: one coordinate 1, thirteen equal to t: H minimal at t0 = {root[0]:.6f} (root of 325t^3+23t^2+9t-1), H(t0) = {vals[j]:.4f}; the paper\'s displayed bound > 796.4 holds exactly')
    worst = 0.0
    for _ in range(20000):
        x = rng.dirichlet(np.full(n, rng.choice([0.05, 0.3, 1.0, 5.0])))
        x = np.maximum(x, 1e-300); x = x / x.sum()
        u = x / x.max(); qq = 1 + 2 * u - u * u
        F = math.sqrt(np.prod(qq) * np.sum(x * x / qq))
        R = n * math.exp(np.mean(np.log(x))) / F
        worst = max(worst, R)
    ok(worst < 0.5, f'D: R_14(x) < 1/2 on 20000 random simplex points (max {worst:.4f}; the bound from the minimum is {math.sqrt(196 / Hmin):.4f})')


def isprime(m):
    return m > 1 and all(m % k for k in range(2, int(m ** 0.5) + 1))


def check_E():
    names = [line.split()[1] for line in open(os.path.join(PKG, 'SHA256SUMS.txt')) if line.strip()]
    P = sorted(int(nm.split('_')[1]) for nm in names)
    target = 14 * mpmath.log((mpmath.mpf(K14.denominator) ** 0.5 / mpmath.mpf(K14.numerator) ** 0.5) / 28) - mpmath.log(360360)
    mass = mpmath.fsum(mpmath.log(p) for p in P)
    # claim type per gate: from the per-gate records of LRC15_AUDIT_SUMMARY.json (data; cross-checked against every
    # archive's own SUMMARY.json by the certificate audit)
    import json
    res = json.load(open(os.path.join(PKG, 'LRC15_AUDIT_SUMMARY.json')))['results']
    kind = {r['p']: r['claim'] for r in res}
    near = sorted(p for p in P if kind[p].startswith('divisibility'))
    strict_p = sorted(p for p in P if kind[p].startswith('J(K,p) empty'))
    strict = mpmath.fsum(mpmath.log(p) for p in strict_p)
    assert len(near) + len(strict_p) == 71
    ok(len(P) == len(set(P)) == 71 and all(isprime(p) and p > 15 for p in P),
       f'E: 71 distinct gate primes, all prime and > 15 (coprime to 360360), from {P[0]} to {P[-1]}')
    ok(abs(mass - mpmath.mpf('408.823317378586107568088088794557')) < 1e-25 and mass > target and mass - target > mpmath.log(P[-1]),
       f'E: sum log p = {mpmath.nstr(mass, 22)} > target {mpmath.nstr(target, 15)}; margin {mpmath.nstr(mass - target, 8)} > log {P[-1]} = {mpmath.nstr(mpmath.log(P[-1]), 6)}')
    ok(strict < target and len(near) == 19 and len(strict_p) == 52,
       f'E: the 52 strict J(14,p) = {{}} gates ({strict_p[0]}, {strict_p[1]}, ..) alone give {mpmath.nstr(strict, 10)} < target: the 19 gates {near} closed through the few-exception shift lemma are load-bearing')


K14 = None


def main():
    global K14
    t0 = time.time()
    print('==== A ====', flush=True)
    B, c = check_A()
    K14 = flag_constants(14)[2]
    print('==== B ====', flush=True); check_B(B, c)
    print('==== C ====', flush=True); check_C()
    print('==== D ====', flush=True); check_D()
    print('==== E ====', flush=True); check_E()
    print(f'elapsed {time.time() - t0:.0f} s')
    print('ALL CHECKS PASSED' if all(OKS) else 'SOME CHECK FAILED')


if __name__ == '__main__':
    main()
