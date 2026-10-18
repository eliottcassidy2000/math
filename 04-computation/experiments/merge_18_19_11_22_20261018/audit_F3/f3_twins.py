#!/usr/bin/env python3
"""Audit F3, item (e): finite-N measurements with independent code (written from THM-4690's definitions).

Full transport V = {1}.  Sources: every odd m <= 4N carries K_N(m) = w(m) + thr(m)[m < 2N] (w = 2 if m = 2 mod 3 else
-1; thr = 1, -2, 1 for m = 0, 1, 2 mod 3) and sits on its seed (strip B(x) = 4x + 1 while m = 5 mod 8); c_0 is the
seed sum.  Transport P(s) = seed((3s + 1)/2^a), a = v_2(3s + 1) in {1, 2}, literal clock a + 2 * strips, sign (-1)^clock.
Synchronous steps: every source moves one step per round (1 -> 1 with clock 2).

(1) clock twins at depth k: pairs of sources at the same node after k steps with equal cumulative clock and different
    nodes after k - 1 steps; share with 3 | Delta, Delta = x_small - x_large, and sums of c0 c0' by Delta mod 3.
    Reported for sources with c_0 != 0 (the session's convention) and for all seeds.
(2) equal-depth meet sums at depths 1 and 2 by clock lag (pairs first meeting at that depth), per N, and the energy
    identity E_k - E_(k-1) = 2 * (sum over depth-k meet pairs), E_k = sum over nodes of (signed mass)^2.
Usage: python3 f3_twins.py LOG2N KMAX
"""
import sys, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def seed_strip(m):
    m = m.copy()
    strips = np.zeros_like(m)
    while True:
        sel = (m % 8) == 5
        if not sel.any():
            break
        m[sel] = (m[sel] - 1) // 4
        strips[sel] += 1
    return m, strips


def sources(N):
    m = np.arange(1, 4 * N + 1, 2, dtype=np.int64)
    r = m % 3
    w = np.where(r == 2, 2, -1)
    thr = np.choose(r, [1, -2, 1])
    K = w + np.where(m < 2 * N, thr, 0)
    s, _ = seed_strip(m)
    xs, inv = np.unique(s, return_inverse=True)
    c0 = np.zeros(len(xs), dtype=np.int64)
    np.add.at(c0, inv, K)
    return xs, c0


def step(y):
    v = 3 * y + 1
    a = np.where(v % 4 == 0, 2, 1)
    assert not ((v % 8 == 0) & (a == 2)).any(), 'a >= 3 at a seed'
    z = v >> a
    z, strips = seed_strip(z)
    return z, a + 2 * strips


def twins(xs, c0, ys, ts, k, N, out, label):
    key, prev, tt = ys[k], ys[k - 1], ts[k]
    o = np.lexsort((xs, tt, key))
    K_, P_, T_, X_, C_ = key[o], prev[o], tt[o], xs[o], c0[o]
    brk = np.flatnonzero((np.diff(K_) != 0) | (np.diff(T_) != 0)) + 1
    st = np.concatenate(([0], brk)); en = np.concatenate((brk, [len(K_)]))
    big = np.flatnonzero(en - st >= 2)
    npairs = [0, 0, 0]; sums = [0, 0, 0]
    for g in big:
        s0, e0 = st[g], en[g]
        for i in range(s0, e0):
            for j in range(i + 1, e0):
                if P_[i] == P_[j]:
                    continue
                dm = int((X_[i] - X_[j]) % 3)
                npairs[dm] += 1
                sums[dm] += int(C_[i]) * int(C_[j])
    n = sum(npairs)
    out.append('%s depth %d: twin pairs %d, share 3|Delta = %.4f (%d); sums/N: 3|D %+.4f, D=1 %+.4f, D=2 %+.4f, total %+.5f'
               % (label, k, n, npairs[0] / max(n, 1), npairs[0], sums[0] / N, sums[1] / N, sums[2] / N, sum(sums) / N))
    print(out[-1], flush=True)


def lag_sums(ys, ts, c0, k, N):
    """equal-depth meet sums at depth k by clock lag (pairs at the same node after k steps, different nodes after
    k - 1 steps), via aggregated (node, clock) and (node, prev, clock) masses."""
    sgn = np.where(ts[k] % 2 == 0, 1, -1)
    val = sgn * c0
    TB = 512
    assert ts[k].max() < TB // 2
    def agg(keys):
        u, inv = np.unique(keys, return_inverse=True)
        s = np.zeros(len(u), dtype=np.int64)
        np.add.at(s, inv, val)
        return u, s
    # node index compression
    nodes, nidx = np.unique(ys[k], return_inverse=True)
    pr_, pidx = np.unique(ys[k - 1], return_inverse=True)
    kT = nidx.astype(np.int64) * TB + ts[k]
    kS = (nidx.astype(np.int64) * len(pr_) + pidx) * TB + ts[k]
    uT, sT = agg(kT)
    uS, sS = agg(kS)
    res = {}
    for D in range(0, 40):
        def corr(u, s):
            pos = np.searchsorted(u, u + D)
            pos = np.minimum(pos, len(u) - 1)
            ok = u[pos] == u + D
            return int((s[ok] * s[pos[ok]]).sum())
        tot = corr(uT, sT) - corr(uS, sS)
        if D == 0:
            tot //= 2
        if tot:
            res[D] = tot / N
    # energies
    def energy(keys):
        u, inv = np.unique(keys, return_inverse=True)
        s = np.zeros(len(u), dtype=np.int64)
        np.add.at(s, inv, np.where(ts[k] % 2 == 0, 1, -1) * c0)
        return int((s * s).sum())
    return res


def energy_at(ys, ts, c0, k):
    sgn = np.where(ts[k] % 2 == 0, 1, -1)
    u, inv = np.unique(ys[k], return_inverse=True)
    s = np.zeros(len(u), dtype=np.int64)
    np.add.at(s, inv, sgn * c0)
    return int((s * s).sum())


def main():
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    KMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    N = 2 ** L
    xs, c0 = sources(N)
    out = ['N = 2^%d: %d seeds <= 4N, %d with c_0 != 0; E_0/N = %.6f (34/7 = %.6f); c_0 values %s'
           % (L, len(xs), int((c0 != 0).sum()), float((c0 * c0).sum()) / N, 34 / 7, sorted(set(c0.tolist())))]
    print(out[-1], flush=True)
    for label, sel in (('c0!=0', c0 != 0), ('all seeds', np.ones_like(c0, dtype=bool))):
        X, C = xs[sel], c0[sel]
        ys = [X]; ts = [np.zeros_like(X)]
        for k in range(KMAX):
            z, tau = step(ys[-1])
            ys.append(z); ts.append(ts[-1] + tau)
        for k in range(1, KMAX + 1):
            twins(X, C, ys, ts, k, N, out, label)
        if label == 'c0!=0':
            E = [energy_at(ys, ts, C, k) / N for k in range(0, 4)]
            out.append('energies E_k/N, k = 0..3: %s' % ['%.6f' % e for e in E])
            print(out[-1], flush=True)
            for k in (1, 2):
                res = lag_sums(ys, ts, C, k, N)
                tot = sum(res.values())
                out.append('depth-%d meet sums by lag (per N): total %.6f, (E_%d - E_%d)/2 = %.6f' % (k, tot, k, k - 1, (E[k] - E[k - 1]) / 2))
                print(out[-1], flush=True)
                out.append('   ' + ', '.join('lag %d: %+.5f' % (D, v) for D, v in sorted(res.items()) if abs(v) > 5e-6))
                print(out[-1], flush=True)
    open(os.path.join(HERE, 'f3_twins_N2p%d.out' % L), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
