"""Smith-form-free certificates for the HYP-9445 (1)-(2) torsion verdicts (S_2 and M_2), every strong class n <= 8.

torsion      <=  a prime p with rank_Fp(X) < rank_Q(X)   (p taken from the smallest nonzero invariant factor > 1
                 reported by FLINT, then CHECKED by an independent mod-p elimination and FLINT's rational rank)
torsion-free <=  |det X| = 1 (nonsingular), or gcd of several nonzero maximal (r x r) minors = 1 (singular),
                 since the r-th determinantal divisor divides every r x r minor.
"""
import random
import math
import numpy as np
import flint
from alib import (gentourng, adj_from_ascii, path_table, ham_cycles, m2_matrix, snf_diag, rank_mod_p, det_exact)

random.seed(20261029)
BIGP = 2147483647


def small_prime_factor(d):
    for q in range(2, 100000):
        if d % q == 0:
            return q
    return None


def independent_rows_cols(X, p=BIGP, order_r=None, order_c=None):
    """Greedy rows/cols giving a nonsingular r x r minor (mod a big prime, hence over Z)."""
    X = np.array(X, dtype=object)
    rows, cols = X.shape
    order_r = list(range(rows)) if order_r is None else order_r
    order_c = list(range(cols)) if order_c is None else order_c
    Y = np.array([[int(X[i, j]) % p for j in order_c] for i in order_r], dtype=object)
    # row-reduce Y mod p choosing pivots in the given orders
    piv_r, piv_c = [], []
    Y = Y.copy()
    used = [False] * rows
    for c in range(cols):
        pr = next((i for i in range(rows) if not used[i] and Y[i, c] % p), None)
        if pr is None:
            continue
        used[pr] = True
        piv_r.append(order_r[pr])
        piv_c.append(order_c[c])
        inv = pow(int(Y[pr, c]), p - 2, p)
        for i in range(rows):
            if i != pr and Y[i, c] % p:
                f = (int(Y[i, c]) * inv) % p
                Y[i] = (Y[i] - f * Y[pr]) % p
    return piv_r, piv_c


def certify_tf(X):
    fm = flint.fmpz_mat([[int(x) for x in row] for row in X])
    r = int(fm.rank())
    n = X.shape[0]
    if r == n:
        return abs(int(fm.det())) == 1, 'det'
    g = 0
    for t in range(40):
        orr = list(range(n)); occ = list(range(n))
        if t:
            random.shuffle(orr); random.shuffle(occ)
        R, C = independent_rows_cols(X, order_r=orr, order_c=occ)
        assert len(R) == r
        sub = flint.fmpz_mat([[int(X[i, j]) for j in C] for i in R])
        g = math.gcd(g, int(sub.det()))
        if g == 1:
            return True, 'minors(%d)' % (t + 1)
    return False, 'gcd=%d' % g


def certify_torsion(X):
    fm = flint.fmpz_mat([[int(x) for x in row] for row in X])
    rq = int(fm.rank())
    if rq == X.shape[0]:
        dt = abs(int(fm.det()))
        return dt > 1, 'nonsingular |det|=%d' % dt
    d = snf_diag(X)
    big = sorted(x for x in d if x > 1)
    if not big:
        return False, 'no factor > 1'
    q = small_prime_factor(big[0])
    if q is None:
        q = int(min(f for f, e in flint.fmpz(big[0]).factor()))
    rp = rank_mod_p(X, q) if q < 3 * 10 ** 9 else None
    if rp is None:
        return False, 'prime too large for int64 elimination'
    return rp < rq, 'p=%d rank_Q=%d rank_Fp=%d' % (q, rq, rp)


for n in range(3, 9):
    res = dict(S_tors_ok=0, M_tors_ok=0, n_hc2=0, S_tf_ok=0, S_tf_claims=0, M_tf_ok=0, M_tf_claims=0,
               S_tors_unique_ok=0, S_tors_unique=0, M_tors_unique_ok=0, M_tors_unique=0)
    fails = []
    for s in gentourng(n):
        A = adj_from_ascii(n, s)
        H = path_table(A)
        hc = ham_cycles(A, H)
        if hc == 0:
            continue
        M, pr = m2_matrix(A, H)
        m = M.shape[0] // 2
        S = np.array([[M[2 * a, 2 * b] - M[2 * a, 2 * b + 1] for b in range(m)] for a in range(m)], dtype=np.int64)
        if hc >= 2:
            res['n_hc2'] += 1
            okS, why = certify_torsion(S)
            okM, whyM = certify_torsion(M)
            res['S_tors_ok'] += okS
            res['M_tors_ok'] += okM
            if not (okS and okM):
                fails.append((s, hc, why, whyM))
        else:
            dS, dM = snf_diag(S), snf_diag(M)
            if any(x > 1 for x in dS):
                res['S_tors_unique'] += 1
                ok, why = certify_torsion(S)
                res['S_tors_unique_ok'] += ok
            else:
                res['S_tf_claims'] += 1
                ok, why = certify_tf(S)
                res['S_tf_ok'] += ok
                if not ok:
                    fails.append((s, hc, 'S tf not certified', why))
            if any(x > 1 for x in dM):
                res['M_tors_unique'] += 1
                ok, why = certify_torsion(M)
                res['M_tors_unique_ok'] += ok
            else:
                res['M_tf_claims'] += 1
                ok, why = certify_tf(M)
                res['M_tf_ok'] += ok
                if not ok:
                    fails.append((s, hc, 'M tf not certified', why))
    print('n=%d' % n, res, 'failures:', fails, flush=True)
