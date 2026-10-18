#!/usr/bin/env python3
"""S41 audit A14 (bonus; THM-4667 statement 7 'Matching' bullet, added after the audit request): compare the actual
transported coefficients c_d(t) (own direct transport) with the ideal K_d(t mod 3^(d+1), b(t)), band
M_d/2^(b+1) < t/N <= M_d/2^b, M_d = 4 (3/2)^d, for d <= 6 and N = 2^16, 2^18, at every seed t (ROOT included), using
own K_d tables (forward push, exact integers, bands < 64).
Usage: python a14_l2_matching_bonus.py
"""
import os, sys, json
sys.dont_write_bytecode = True
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
K0 = np.array([[-1, -1, 2], [0, -3, 3], [-1, -1, 2], [-3, 0, 3], [-1, -1, 2], [0, 0, 0]], dtype=np.int64)


def Fmap(Q, A):
    a = 1 if A % 2 == 1 else 2
    cA = (3 + 2 ** a - 2 ** A) // 3
    inv = pow(2, -A, Q)
    return lambda s: ((3 * s + cA) * inv) % Q


def tables(D, BR):
    K = np.array([K0[b % 6] for b in range(BR)]).T.copy()
    out = [K]
    for d in range(1, D + 1):
        Q = 3 ** (d + 1)
        q = Q // 3
        new = np.zeros((Q, BR), dtype=np.int64)
        for o in range(BR):
            f = Fmap(Q, o + 1)
            sg = 1 if (o + 1) % 2 == 0 else -1
            tgt = np.array([f(s) for s in range(q)])
            new[tgt, o:] += sg * K[:, :BR - o]
        K = new
        out.append(K)
    return out


def seed_vec(u):
    u = u.copy()
    m = (u % 8) == 5
    while m.any():
        u[m] = (u[m] - 1) // 4
        m = (u % 8) == 5
    return u


def aggregate(keys, vals):
    order = np.argsort(keys, kind='stable')
    k = keys[order]
    v = vals[order]
    starts = np.flatnonzero(np.r_[True, k[1:] != k[:-1]])
    s = np.add.reduceat(v, starts)
    return k[starts], s


def bands_vec(t, N, d):
    """b with 2^b 2^d t <= 4 3^d N < 2^(b+1) 2^d t, exact int64 correction of a float guess (vectorised)."""
    num = np.int64(4 * 3 ** d * N)
    den = (2 ** d) * t.astype(np.int64)
    b = np.floor(np.log2(num / den.astype(np.float64))).astype(np.int64)
    b = np.maximum(b, 0)
    for _ in range(3):
        hi = np.left_shift(den, b) > num
        b = np.where(hi & (b > 0), b - 1, b)
        lo = np.left_shift(den, b + 1) <= num
        b = np.where(lo, b + 1, b)
    assert np.all(np.left_shift(den, b) <= num) and np.all(np.left_shift(den, b + 1) > num)
    return b


def main():
    D = 6
    Ks = tables(D, 64)
    res = {}
    for N in (2 ** 16, 2 ** 18):
        m = np.arange(1, 4 * N + 1, 2, dtype=np.int64)
        r = m % 3
        k = np.where(m < 2 * N, -3 * np.where(r == 0, 0, np.where(r == 1, 1, -1)), np.where(r == 2, 2, -1)).astype(np.int64)
        S, C = aggregate(seed_vec(m), k)
        for d in range(0, D + 1):
            if d > 0:
                y = 3 * S + 1
                a = np.where((y & (-y)) == 4, 2, 1)
                S, C = aggregate(seed_vec(y >> a), np.where(a == 2, C, -C))
            Q = 3 ** (d + 1)
            # ideal value at every seed of the support, plus every seed where the ideal value is nonzero is needed too:
            # restrict to seeds t <= 4 N (3/2)^d (bands b >= 0); compare on the union of supports
            tmax = (4 * 3 ** d * N) // (2 ** d)
            allt = np.arange(1, tmax + 2 * int(1.5 ** d) + 8, 2, dtype=np.int64)
            allt = allt[allt % 8 != 5]
            inband = allt * 2 ** d <= 4 * 3 ** d * N          # bands b >= 0; beyond, the ideal coefficient is 0
            ideal = np.zeros(len(allt), dtype=np.int64)
            bands = bands_vec(allt[inband], N, d)
            assert bands.max() < 64
            ideal[inband] = Ks[d][allt[inband] % Q, bands]
            actual = np.zeros(len(allt), dtype=np.int64)
            pos = np.searchsorted(allt, S)
            inside = (pos < len(allt))
            inside[inside] &= allt[pos[inside]] == S[inside]
            assert np.all(inside | (C == 0)), 'support beyond the ideal band range'
            actual[pos[inside]] = C[inside]
            diff = actual - ideal
            bad = np.flatnonzero(diff)
            err_l2 = float(np.sqrt(np.sum(diff.astype(np.float64) ** 2)))
            off_root = [int(t) for t in allt[bad] if t != 1]
            res['%d_%d' % (N, d)] = dict(n_diff=len(bad), seeds=[int(t) for t in allt[bad]][:12], l2_err_over_N=err_l2 / N,
                                         off_root_l2_over_N=float(np.sqrt(np.sum(diff[allt != 1].astype(np.float64) ** 2))) / N)
            print('N = 2^%d d = %d: seeds where actual != ideal: %d %s; l2 error/N %.5f (off ROOT %.5f)'
                  % (N.bit_length() - 1, d, len(bad), [int(t) for t in allt[bad]][:8], err_l2 / N,
                     res['%d_%d' % (N, d)]['off_root_l2_over_N']), flush=True)
    json.dump(res, open(os.path.join(HERE, 'a14_l2_matching_bonus.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
