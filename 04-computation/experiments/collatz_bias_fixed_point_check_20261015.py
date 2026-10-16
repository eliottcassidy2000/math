#!/usr/bin/env python3
"""Self-consistency test of HYP-9282: do the measured square-root biases satisfy the fixed-point equation of THM-4658 (4b)?

For every cylinder A = a + 3^q Z_3 (q <= 5) let b(A) be the mean over octaves 2^14..2^(LOGMAX-1) of the octave-averaged
normalized excess [C(A; N) - Haar(A) S(N)] / sqrt(N), from log-uniform samples (64 per octave).  If C(A; N) - Haar(A)S(N)
= b(A) sqrt(N) + (mean-zero fluctuation) with b independent of the octave, THM-4658 (4b) gives
    b(A) = -2^(-1/2) b(A/2) - (3/2)^(1/2) [ b((3A+1)/2) - Haar(A) b(2 + 3Z_3) ].
The script reports, level by level (q = 1..4, using level q+1 on the right), the residual of this equation against the
size of b, and the correlation between b at level q and its prediction from levels q and q+1.  Octave-averaging and the
shifted scale ranges make the relation approximate; the residuals are compared with the octave-to-octave scatter.
CAUTION (result of the run): because THM-4658 statement 1 is an exact identity at every N, octave-averaged means satisfy
this fixed-point equation automatically up to end effects of relative size about 1/(number of octaves), whatever the
sequence.  The check therefore validates the code and the identity; it is not evidence for stationary main terms.
Session opus-2026-10-15-S37.  Usage: python collatz_bias_fixed_point_check_20261015.py LOGMAX out.json
"""
import sys, os, json, time, math
import importlib.util
import numpy as np

import numba as nb

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('cen', os.path.join(HERE, 'collatz_octave_memory_census_20261015.py'))
cen = importlib.util.module_from_spec(spec)
_argv = sys.argv
sys.argv = ['x', '24', '20', 'unused.json', 'unused.npz']
spec.loader.exec_module(cen)
sys.argv = _argv


@nb.njit
def getbit(P, n):
    return (P[n >> 6] >> np.uint64(n & 63)) & np.uint64(1)


@nb.njit
def sample(P, logmax, grid, M):
    NMAX = 1 << logmax
    G = grid.shape[0]
    out = np.zeros((G, M), dtype=np.int64)
    C = np.zeros(M, dtype=np.int64)
    r = 1
    gi = 0
    for N in range(1, NMAX + 1):
        t = 1 - 2 * int(getbit(P, N))
        C[r] += t
        while gi < G and grid[gi] == N:
            for a in range(M):
                out[gi, a] = C[a]
            gi += 1
        r += 1
        if r == M:
            r = 0
    return out
LOGMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 32
OUT = sys.argv[2] if len(sys.argv) > 2 else 'bias_fixed_point.json'
QMAX = 5

if __name__ == '__main__':
    t0 = time.time()
    P = cen.parity_bits(LOGMAX)
    raw = sorted(set(int(round(2 ** (10 + j / 64))) for j in range(0, (LOGMAX - 10) * 64 + 1)))
    grid = np.array([g for g in raw if g <= (1 << LOGMAX)], dtype=np.int64)
    C243 = sample(P, LOGMAX, grid, 3 ** QMAX).astype(np.float64)
    N = grid.astype(np.float64)
    u = np.log2(N)
    S = C243.sum(axis=1)
    ks = list(range(14, LOGMAX))
    b = {}
    sd = {}
    for q in range(0, QMAX + 1):
        m = 3 ** q
        Cq = np.stack([C243[:, a::m].sum(axis=1) for a in range(m)], axis=1)
        ex = (Cq - S[:, None] / m) / np.sqrt(N)[:, None]
        om = np.stack([ex[(u >= k) & (u < k + 1)].mean(axis=0) for k in ks])
        b[q] = om.mean(axis=0)
        sd[q] = om.std(axis=0, ddof=1) / math.sqrt(len(ks))
    print(f'biases computed to 2^{LOGMAX}: {time.time() - t0:.1f} s', flush=True)
    res = {'logmax': LOGMAX, 'levels': []}
    b2 = b[1][2]
    for q in range(1, QMAX):
        m = 3 ** q
        inv2 = pow(2, -1, m)
        inv2b = pow(2, -1, 3 * m)
        pred = np.zeros(m)
        for a in range(m):
            half = (a * inv2) % m
            child = (inv2b * (3 * a + 1)) % (3 * m)
            pred[a] = -2 ** -0.5 * b[q][half] - math.sqrt(1.5) * (b[q + 1][child] - b2 / m)
        resid = b[q] - pred
        corr = float(np.corrcoef(b[q], pred)[0, 1])
        rel = float(np.sqrt((resid ** 2).mean() / (b[q] ** 2).mean()))
        noise = float(np.sqrt((sd[q] ** 2).mean()))
        res['levels'].append(dict(q=q, corr_b_pred=corr, rms_b=float(np.sqrt((b[q] ** 2).mean())),
                                  rms_resid=float(np.sqrt((resid ** 2).mean())), rel_resid=rel, rms_standard_error=noise))
        print(f'level {q}: rms b = {np.sqrt((b[q] ** 2).mean()):.4f}, rms residual = {np.sqrt((resid ** 2).mean()):.4f} '
              f'(relative {rel:.3f}), rms standard error of b = {noise:.4f}, corr(b, prediction) = {corr:.4f}', flush=True)
    res['b_level1'] = [float(x) for x in b[1]]
    json.dump(res, open(OUT, 'w'), indent=1)
