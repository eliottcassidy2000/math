#!/usr/bin/env python3
"""Matched null for HYP-9282: random signs that keep the exact halving relation t(2m) = -t(m).

Collatz parity satisfies t(2m) = -t(m) exactly (tau(2m) = tau(m) + 1), which creates the octave-to-octave inheritance
C(A; N) = -C(A/2; N/2) + (odd part) of every class sum, and hence octave memory in the channels.  This null draws
t(n) = +-1 independently for odd n and sets t(2^v m) = (-1)^v t(m); it has the same inheritance and no odd-step
structure.  For each replicate it computes the channels g_(q,r)(N) = H_(q,r)(N)/sqrt(6N/3^q) (q <= 5) on the grid
N_j = round(2^(10+j/64)) and their octave means over the octaves 14..LOGMAX-1, exactly as for the Collatz data, and
reports the naive t-values (mean over octaves / (sd/sqrt(#octaves))), the counts |t| > 5 and |t| > 8, and the pooled
octave-mean autocorrelations by level.  The same statistics are computed for the Collatz census on the same octaves.
Session opus-2026-10-15-S37.  Usage: python collatz_bias_halving_null_20261015.py LOGMAX NREP collatz_ext.npz out.json
"""
import sys, json, math, time
import numpy as np
import numba as nb

LOGMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 30
NREP = int(sys.argv[2]) if len(sys.argv) > 2 else 4
CNPZ = sys.argv[3] if len(sys.argv) > 3 else None
OUT = sys.argv[4] if len(sys.argv) > 4 else 'halving_null.json'
M = 243


@nb.njit
def sample_null(seed, logmax, grid, Mod):
    np.random.seed(seed)
    NMAX = 1 << logmax
    t = np.zeros(NMAX + 1, dtype=np.int8)
    G = grid.shape[0]
    out = np.zeros((G, Mod), dtype=np.int64)
    C = np.zeros(Mod, dtype=np.int64)
    r = 1
    gi = 0
    for n in range(1, NMAX + 1):
        if n & 1:
            t[n] = 1 if np.random.random() < 0.5 else -1
        else:
            t[n] = -t[n >> 1]
        C[r] += t[n]
        while gi < G and grid[gi] == n:
            for a in range(Mod):
                out[gi, a] = C[a]
            gi += 1
        r += 1
        if r == Mod:
            r = 0
    return out


def channel_stats(grid, C, ks):
    N = grid.astype(np.float64)
    u = np.log2(N)
    lev = {0: C.sum(axis=1)[:, None].astype(np.float64)}
    for q in range(1, 6):
        m = 3 ** q
        lev[q] = np.stack([C[:, a::m].sum(axis=1) for a in range(m)], axis=1).astype(np.float64)
    cols, levs = [], []
    for q in range(1, 6):
        m, mp = 3 ** q, 3 ** (q - 1)
        for r in range(m):
            cols.append((3 * lev[q][:, r] - lev[q - 1][:, r % mp]) / np.sqrt(6 * N / m))
            levs.append(q)
    Gm = np.stack(cols, axis=1)
    OM = np.stack([Gm[(u >= k) & (u < k + 1)].mean(axis=0) for k in ks])
    mean = OM.mean(axis=0)
    sd = OM.std(axis=0, ddof=1)
    t = mean / (sd / math.sqrt(len(ks)))
    levs = np.array(levs)
    X = OM - mean
    v = (X ** 2).mean(axis=0)
    acf = {}
    for q in range(1, 6):
        s = levs == q
        acf[q] = [float(np.mean((X[L:, s] * X[:-L, s]).mean(axis=0) / v[s])) for L in (1, 2, 3)]
    return dict(n_t5=int((np.abs(t) > 5).sum()), n_t8=int((np.abs(t) > 8).sum()), n_t3=int((np.abs(t) > 3).sum()),
                max_abs_t=float(np.abs(t).max()), rms_mean_by_level={q: float(np.sqrt((mean[levs == q] ** 2).mean())) for q in range(1, 6)},
                sd_by_level={q: float(np.sqrt((sd[levs == q] ** 2).mean())) for q in range(1, 6)}, acf_by_level=acf,
                t_quantiles=[float(np.quantile(np.abs(t), p)) for p in (0.5, 0.9, 0.99)])


if __name__ == '__main__':
    raw = sorted(set(int(round(2 ** (10 + j / 64))) for j in range(0, (LOGMAX - 10) * 64 + 1)))
    grid = np.array([g for g in raw if g < (1 << LOGMAX)], dtype=np.int64)
    ks = list(range(14, LOGMAX))
    res = dict(logmax=LOGMAX, octaves=[ks[0], ks[-1]], null=[])
    for rep in range(NREP):
        t0 = time.time()
        C = sample_null(1000 + rep, LOGMAX, grid, M)
        st = channel_stats(grid, C, ks)
        res['null'].append(st)
        print(f'null replicate {rep}: |t|>5: {st["n_t5"]}, |t|>8: {st["n_t8"]}, |t|>3: {st["n_t3"]}, max |t| {st["max_abs_t"]:.1f}, '
              f'acf lag1 by level {[round(st["acf_by_level"][q][0], 2) for q in range(1, 6)]} ({time.time() - t0:.1f} s)', flush=True)
    if CNPZ:
        Z = np.load(CNPZ)
        g2 = Z['grid']
        sel = g2 < (1 << LOGMAX)
        st = channel_stats(g2[sel], Z['classsums'][sel], ks)
        res['collatz_same_octaves'] = st
        print(f'Collatz, same octaves: |t|>5: {st["n_t5"]}, |t|>8: {st["n_t8"]}, |t|>3: {st["n_t3"]}, max |t| {st["max_abs_t"]:.1f}, '
              f'acf lag1 by level {[round(st["acf_by_level"][q][0], 2) for q in range(1, 6)]}', flush=True)
        print('Collatz rms mean by level', st['rms_mean_by_level'], 'sd by level', st['sd_by_level'])
        print('null rms mean by level (rep 0)', res['null'][0]['rms_mean_by_level'], 'sd', res['null'][0]['sd_by_level'])
    json.dump(res, open(OUT, 'w'), indent=1)
