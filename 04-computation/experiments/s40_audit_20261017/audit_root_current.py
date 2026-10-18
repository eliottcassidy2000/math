"""S40 audit: independent recomputation of the ROOT current a_N = H(4N) - H(N) on octaves X = 2^k (rows X <= N < 2X).

T(n) = n/2 (n even), (3n+1)/2 (n odd); d(n) = number of T-steps from n to 1 (full value stored, int16);
w(n) = 3[n = 2 mod 3] - 1; f = (-1)^d; H(M) = sum_(n <= M) w(n) f(n).
Exact integer outputs per octave: energy numerator sum a_N^2, tail numerator max_t t^2 #{N : |a_N| >= t}
(maximised over t in the set of values |a_N|), positive fraction, sign changes, min/max |a_N|; fold
I(N) = -3 sum_(m <= 2N odd) chi(m) f(m), B(N) = sum_(2N < m <= 4N odd) w(m) f(m): I + B = a_N at every row, and the
numerators of ||I||^2, ||B||^2, <I, B>. Spot check of d(n) by direct iteration in Python integers.
Compared with the session's JSONs (read-only).
Usage: python audit_root_current.py KMAX
"""
import sys, os, json, time, random
import numpy as np
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', '..', '..', '05-knowledge', 'results')


@njit(cache=True)
def stopping_depths(L):
    d = np.zeros(L, dtype=np.int16)
    for n in range(2, L):
        m = n
        s = 0
        while m >= n:
            if m & 1:
                m = (3 * m + 1) >> 1
            else:
                m >>= 1
            s += 1
        d[n] = s + d[m]
    return d


@njit(cache=True)
def prefixes(d):
    L = d.shape[0]
    H = np.zeros(L, dtype=np.int32)
    Pc = np.zeros(L, dtype=np.int32)
    Pw = np.zeros(L, dtype=np.int32)
    h = 0; pc = 0; pw_ = 0
    for n in range(1, L):
        f = 1 if (d[n] & 1) == 0 else -1
        r = n % 3
        w = 2 if r == 2 else -1
        h += w * f
        if n & 1:
            pw_ += w * f
            if r == 1:
                pc += f
            elif r == 2:
                pc -= f
        H[n] = h; Pc[n] = pc; Pw[n] = pw_
    return H, Pc, Pw


@njit(cache=True)
def octave(H, Pc, Pw, X):
    a = np.empty(X, dtype=np.int64)
    I = np.empty(X, dtype=np.int64)
    B = np.empty(X, dtype=np.int64)
    bad = 0
    for i in range(X):
        N = X + i
        a[i] = np.int64(H[4 * N]) - np.int64(H[N])
        I[i] = -3 * np.int64(Pc[2 * N])
        B[i] = np.int64(Pw[4 * N]) - np.int64(Pw[2 * N])
        if I[i] + B[i] != a[i]:
            bad += 1
    return a, I, B, bad


def d_direct(n):
    s = 0
    while n != 1:
        n = n // 2 if n % 2 == 0 else (3 * n + 1) // 2
        s += 1
    return s


def main():
    KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 22
    KMIN = 8
    t0 = time.time()
    L = 4 * 2 ** (KMAX + 1) + 1
    d = stopping_depths(L)
    random.seed(1)
    for n in [1, 2, 3, 7, 27, 97, 703, 871] + [random.randrange(2, L) for _ in range(2000)]:
        assert d[n] == d_direct(n), n
    H, Pc, Pw = prefixes(d)
    del d
    sess = {r['k']: r for r in json.load(open(os.path.join(RES, 'collatz_root_current_exact_20261017.json')))['rows']}
    sessf = {r['k']: r for r in json.load(open(os.path.join(RES, 'collatz_root_current_fold_full_20261017.json')))}
    rows = []
    mism = []
    for k in range(KMIN, KMAX + 1):
        X = 2 ** k
        a, I, B, bad = octave(H, Pc, Pw, X)
        e_num = sum(int(x) * int(x) for x in a) if X <= 2 ** 16 else int(np.sum(a * a))   # int64 safe: X*max^2 < 2^62 checked below
        mx = int(np.max(np.abs(a)))
        assert X * mx * mx < 2 ** 62
        absa = np.sort(np.abs(a))[::-1]
        K_num = int(np.max(absa * absa * np.arange(1, X + 1, dtype=np.int64)))
        ii = int(np.sum(I * I)); bb = int(np.sum(B * B)); ib = int(np.sum(I * B))
        assert ii + bb + 2 * ib == e_num
        pos = int(np.sum(a > 0))
        sg = np.sign(a)
        chg = int(np.sum(sg[1:] != sg[:-1]))
        row = dict(k=k, energy_num=e_num, K_num=K_num, den=X * X, frac_pos=pos / X, sign_changes=chg, min_abs=int(absa[-1]), max_abs=mx,
                   fold_failures=bad, II=ii, BB=bb, IB=ib, energy=e_num / X ** 2, K=K_num / X ** 2)
        rows.append(row)
        s = sess.get(k); sf = sessf.get(k)
        ok = (s is not None and s['energy_numerator'] == e_num and s['K_numerator'] == K_num and s['sign_changes'] == chg
              and abs(s['frac_pos'] - pos / X) < 1e-12 and s['min_abs'] == int(absa[-1]) and s['max_abs'] == mx)
        okf = (sf is not None and sf['II_num'] == ii and sf['BB_num'] == bb and sf['IB_num'] == ib and sf['failures'] == 0) if k <= 26 else None
        if not ok or okf is False:
            mism.append(k)
        print('k=%2d energy %.4f K %.4f pos %.3f chg %d min %d max %d | fold fails %d | II %.4f BB %.4f IB %.4f | session match %s fold match %s'
              % (k, e_num / X ** 2, K_num / X ** 2, pos / X, chg, int(absa[-1]), mx, bad, ii / X ** 2, bb / X ** 2, ib / X ** 2, ok, okf),
              flush=True)
    print('mismatching octaves:', mism, '; %.1f s' % (time.time() - t0))
    json.dump(rows, open(os.path.join(HERE, 'audit_root_current.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
