#!/usr/bin/env python3
"""Residue-class biases of the Collatz parity census t(n) = (-1)^tau(n) at levels q = 1..5 (moduli 3..243).

For each q and residue r mod 3^q, the wavelet channel H_(q,r)(N) = 3 C_(3^q)[r](N) - C_(3^(q-1))[r mod 3^(q-1)](N),
where C_m[a](N) = sum_(n <= N, n = a mod m) t(n), is normalised by sqrt(N w_q), w_q = 6/3^q (its variance scale for
independent signs).  The script samples all class sums at N_j = round(2^(10 + j/64)) up to 2^LOGMAX, and for every
channel reports the mean of the normalised value over each octave 2^14..2^LOGMAX, the mean over octaves, a t-statistic
over octaves, the number of octaves with positive mean, and a linear drift per octave.  It also reports the
multiples-of-3^j splits D_j(N) = sum_(3^j || n) chi(n/3^j) t(n) (chi the nontrivial character mod 3) normalised by
sqrt(N), j = 0..3, and checks the exact identity C_m sums against the total S.
Parity bits come from collatz_octave_memory_census_20261015.parity_bits.  Session opus-2026-10-15-S37.
Usage: python collatz_residue_bias_deep_levels_20261015.py LOGMAX out.json
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

LOGMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 32
OUT = sys.argv[2] if len(sys.argv) > 2 else 'bias_deep.json'
QMAX = 5
M = 3 ** QMAX


@nb.njit(cache=True)
def sample(P, logmax, grid, M):
    NMAX = 1 << logmax
    G = grid.shape[0]
    out = np.zeros((G, M), dtype=np.int64)
    C = np.zeros(M, dtype=np.int64)
    r = 1
    gi = 0
    for N in range(1, NMAX + 1):
        t = 1 - 2 * int(cen.getbit(P, N))
        C[r] += t
        while gi < G and grid[gi] == N:
            for a in range(M):
                out[gi, a] = C[a]
            gi += 1
        r += 1
        if r == M:
            r = 0
    return out


if __name__ == '__main__':
    t0 = time.time()
    P = cen.parity_bits(LOGMAX)
    raw = sorted(set(int(round(2 ** (10 + j / 64))) for j in range(0, (LOGMAX - 10) * 64 + 1)))
    grid = np.array([g for g in raw if g <= (1 << LOGMAX)], dtype=np.int64)
    C243 = sample(P, LOGMAX, grid, M).astype(np.float64)
    print(f'class sums mod 243 at {len(grid)} grid points to 2^{LOGMAX}: {time.time() - t0:.1f} s', flush=True)
    N = grid.astype(np.float64)
    u = np.log2(N)
    S = C243.sum(axis=1)

    def classes(q):
        m = 3 ** q
        return np.stack([C243[:, a::m].sum(axis=1) for a in range(m)], axis=1)
    Cq = {q: classes(q) for q in range(0, QMAX + 1)}
    ks = list(range(14, LOGMAX))
    res = {'logmax': LOGMAX, 'octaves': ks, 'channels': [], 'splits': []}
    for q in range(1, QMAX + 1):
        wq = 6.0 / 3 ** q
        for r in range(3 ** q):
            H = 3 * Cq[q][:, r] - Cq[q - 1][:, r % 3 ** (q - 1)]
            g = H / np.sqrt(N * wq)
            om = np.array([g[(u >= k) & (u < k + 1)].mean() for k in ks])
            t = om.mean() / (om.std(ddof=1) / math.sqrt(len(om)))
            res['channels'].append(dict(q=q, r=r, mean=float(om.mean()), sd=float(om.std(ddof=1)), t=float(t),
                                        positive_octaves=int((om > 0).sum()), drift=float(np.polyfit(ks, om, 1)[0])))
    for j in range(0, 4):
        m1 = 3 ** (j + 1)
        a1, a2 = 3 ** j, 2 * 3 ** j
        D = Cq[j + 1][:, a1 % m1] - Cq[j + 1][:, a2 % m1]
        g = D / np.sqrt(N)
        om = np.array([g[(u >= k) & (u < k + 1)].mean() for k in ks])
        res['splits'].append(dict(j=j, mean_over_sqrtN=float(om.mean()), sd=float(om.std(ddof=1)),
                                  t=float(om.mean() / (om.std(ddof=1) / math.sqrt(len(om)))),
                                  positive_octaves=int((om > 0).sum()), drift=float(np.polyfit(ks, om, 1)[0]),
                                  octave_means=[float(x) for x in om]))
    res['top_channels'] = sorted(res['channels'], key=lambda c: -abs(c['t']))[:24]
    json.dump(res, open(OUT, 'w'), indent=1)
    print('multiples-of-3^j splits D_j(N)/sqrt(N): sum_(3^j || n) chi(n/3^j) t(n):')
    for s in res['splits']:
        print(f"  j = {s['j']}: mean {s['mean_over_sqrtN']:+.4f}, sd {s['sd']:.4f}, t {s['t']:+.1f}, "
              f"positive octaves {s['positive_octaves']}/{len(ks)}, drift {s['drift']:+.4f}/octave")
    print('strongest channels:')
    for c in res['top_channels'][:16]:
        print(f"  H_({c['q']},{c['r']:3d}): mean {c['mean']:+.3f} sd {c['sd']:.3f} t {c['t']:+6.1f} pos {c['positive_octaves']:2d} drift {c['drift']:+.4f}")
