"""Audit item 2 at large K (independent class-based code; my own Phi, my own sparse class-lemma check).

(a) Class lemma, directly: for every level m and every class C = c mod 2^m, (sum_{a in C} C_K(a, .)) equals
    2^(K-1-m+m') 1_{Phi(C)} where m' is the precision of Phi(C) (scipy.sparse aggregation, exact int64).
(b) Trajectories a -> Phi^n(a): R_n = all classes of precision 1..K-n; multiplicity of each class of precision m >= 2 in
    R_n equals binom(K-m-1, n-1) (closed form derived by the auditor: choose the n-1 intermediate precisions in (m, K));
    hitting times, mean K/2, worst TV, average TV profile, Dobrushin, l^2 identity:
    sum_{m>=2} w_m(n) = P(Bin(K-2, 1/2) >= n), and the full singular spectrum of P (sqrt of tail sums).
(c) Universality: the same for p in {1, -1, 5, 7, 9, 11, 13, 15, 2^k+1 ...}.
"""
import math
import sys
import time
from fractions import Fraction

import numpy as np
import scipy.sparse as sps


def v2_vec(t, cap):
    """2-adic valuation by repeated halving (independent of the author's frexp trick); t > 0 int64 array; capped."""
    v = np.zeros_like(t)
    cur = t.copy()
    for _ in range(cap):
        even = (cur & 1) == 0
        if not even.any():
            break
        v[even] += 1
        cur[even] >>= 1
    return v


def build_C_sparse(K, p):
    n = 1 << (K - 1)
    mod = 1 << K
    a = 2 * np.arange(n, dtype=np.int64) + 1
    t = (p * a + 1) % mod                                   # in [0, 2^K); 0 only for a*
    special = t == 0
    tt = np.where(special, mod, t)
    v = v2_vec(tt, K + 1)
    v = np.where(special, K, v)
    rows, cols, vals = [], [], []
    for vv in range(1, K + 1):
        sel = np.nonzero(v == vv)[0]
        if len(sel) == 0:
            continue
        if vv == K:
            for i in sel:
                rows.append(np.full(n, i, dtype=np.int64))
                cols.append(np.arange(n, dtype=np.int64))
                vals.append(np.ones(n, dtype=np.int64))
            continue
        r = (tt[sel] >> vv) % (1 << (K - vv))              # odd residue mod 2^(K-v)
        k = np.arange(1 << vv, dtype=np.int64)
        b = (r[:, None] + k[None, :] * (1 << (K - vv)))     # all lifts mod 2^K
        rows.append(np.repeat(sel, 1 << vv))
        cols.append(((b - 1) // 2).ravel())
        vals.append(np.full(len(sel) << vv, 1 << (K - 1 - vv), dtype=np.int64))
    R = np.concatenate(rows); Cc = np.concatenate(cols); V = np.concatenate(vals)
    M = sps.csr_matrix((V, (R, Cc)), shape=(n, n), dtype=np.int64)
    return M, int(a[special][0])


def phi_cls(c, m, p):
    """class map: returns (c', m') with c' reduced mod 2^m' (root = (1, 1))."""
    mod = np.left_shift(np.int64(1), m)
    t = (p * c + 1) % mod
    zero = t == 0
    tt = np.where(zero, mod, t)
    v = v2_vec(tt, 64)
    det = (~zero) & (v < m)
    m2 = np.where(det, m - v, 1)
    c2 = np.where(det, (tt >> np.where(det, v, 0)) % np.left_shift(np.int64(1), m2), 1)
    c2 = np.where(m2 == 1, 1, c2)
    return c2.astype(np.int64), m2.astype(np.int64)


def class_lemma_sparse(K, p):
    C, astar = build_C_sparse(K, p)
    n = 1 << (K - 1)
    assert C.nnz == K * n
    rs = np.asarray(C.sum(axis=1)).ravel(); cs = np.asarray(C.sum(axis=0)).ravel()
    assert np.all(rs == n) and np.all(cs == n)
    a = 2 * np.arange(n, dtype=np.int64) + 1
    for m in range(1, K + 1):
        ncls = 1 << (m - 1)
        cls_of_a = (a % (1 << m) - 1) // 2                    # class index of a at level m
        A = sps.csr_matrix((np.ones(n, dtype=np.int64), (cls_of_a, np.arange(n))), shape=(ncls, n))
        S = (A @ C).tocsr()
        S.sort_indices()
        c = 2 * np.arange(ncls, dtype=np.int64) + 1
        c2, m2 = phi_cls(c, np.full(ncls, m, dtype=np.int64), p)
        # expected: row C has value 2^(K-1-m+m') on exactly the residues b = c2 mod 2^m2
        for i in range(ncls):
            lo, hi = S.indptr[i], S.indptr[i + 1]
            cols = S.indices[lo:hi]; vals = S.data[lo:hi]
            b = 2 * cols + 1
            mm = int(m2[i]); cc = int(c2[i])
            assert len(cols) == (1 << (K - mm)), (K, p, m, i)
            assert np.all((b - cc) % (1 << mm) == 0)
            assert np.all(vals == (1 << (K - 1 - m + mm)))
    return astar


def trajectories(K, p, profile=False):
    n = 1 << (K - 1)
    c = 2 * np.arange(n, dtype=np.int64) + 1
    m = np.full(n, K, dtype=np.int64)
    hit = np.full(n, -1, dtype=np.int64)
    worst_tv, avg_tv, dob, l2 = [], [], [], []
    tails = []
    for s in range(K):
        key = c * 64 + m
        uniq, counts = np.unique(key, return_counts=True)
        um = uniq % 64
        if s >= 1:
            # R_s = all classes of precision 1..K-s, each class (c mod 2^m) exactly once as a key
            assert len(uniq) == (1 << (K - s)) - 1 and set(np.unique(um).tolist()) == set(range(1, K - s + 1))
            for lev in range(2, K - s + 1):
                cl = counts[um == lev]
                assert len(cl) == (1 << (lev - 1))
                assert np.all(cl == math.comb(K - lev - 1, s - 1)), (K, p, s, lev, np.unique(cl))
            assert int(counts[um == 1][0]) == n - sum(math.comb(K - lev - 1, s - 1) * (1 << (lev - 1))
                                                      for lev in range(2, K - s + 1))
        mx = int(m.max())
        worst_tv.append(1 - Fraction(1, 1 << (mx - 1)))
        lev, cnt = np.unique(m, return_counts=True)
        avg_tv.append(sum(Fraction(int(k), n) * (1 - Fraction(1, 1 << (int(L) - 1))) for L, k in zip(lev, cnt)))
        dob.append(1 if mx >= 2 else 0)
        w = {int(L): Fraction(int(k), n) for L, k in zip(lev, cnt)}
        l2.append(sum(x for L, x in w.items() if L >= 2))
        if s == 1:
            # tail sums of w (singular values^2 of P on level-m wavelets): 1 - 2^(m-K) for m = 2..K-1, 0 for m = K
            tails = [sum(x for L, x in w.items() if L >= mm) for mm in range(2, K + 1)]
        hit[(m == 1) & (hit < 0)] = s
        c, m = phi_cls(c, m, p)
    assert np.all(m == 1) and np.all(hit >= 1)
    dist = np.bincount(hit, minlength=K)
    assert all(int(dist[h]) == 2 * math.comb(K - 2, h - 1) for h in range(1, K))
    assert Fraction(int(hit.sum()), n) == Fraction(K, 2)
    assert worst_tv == [1 - Fraction(1, 1 << (K - 1 - s)) for s in range(K)]
    assert dob == [1] * (K - 1) + [0]
    binom = [sum(Fraction(math.comb(K - 2, j), 1 << (K - 2)) for j in range(s, K - 1)) for s in range(K)]
    assert l2 == binom
    assert tails == [1 - Fraction(1 << mm, 1 << K) for mm in range(2, K)] + [0]
    return dict(dist=[int(x) for x in dist[1:]], avg_tv=avg_tv, l2=l2)


def main():
    t0 = time.time()
    print('A3: class-based independent checks')
    for K in range(2, 21):
        t1 = time.time()
        if K <= 18:
            astar = class_lemma_sparse(K, 3)
            cl = f'class lemma checked directly at all {K} levels (sparse aggregation), a* = {astar}; '
        else:
            cl = ''
        r = trajectories(K, 3)
        print(f'p = 3, K = {K:2d}: {cl}R_n = all classes of precision <= K-n with multiplicity binom(K-m-1, n-1); '
              f'hitting counts 2 binom(K-2, h-1), mean K/2, worst TV, Dobrushin, l^2 = binomial tails  [{time.time()-t1:.1f}s]')
        if K == 12:
            print('   K = 12 hitting counts', r['dist'])
            print('   K = 12 average TV profile', [round(float(x), 4) for x in r['avg_tv']])
        sys.stdout.flush()
    for K in (21, 22):
        t1 = time.time()
        trajectories(K, 3)
        print(f'p = 3, K = {K}: trajectory checks pass  [{time.time()-t1:.1f}s]')
    for p in (1, -1, 5, 7, 9, 11, 13, 15, 17, 33, 1023, -3, -5):
        Ks = range(2, 17)
        for K in Ks:
            if K <= 12:
                class_lemma_sparse(K, p)
            trajectories(K, p)
        print(f'p = {p}: class lemma (K <= 12) and all trajectory checks pass for K = 2..16')
        sys.stdout.flush()
    print(f'A3 ALL CHECKS PASSED  [{time.time()-t0:.0f}s]')


if __name__ == '__main__':
    main()
