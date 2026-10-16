#!/usr/bin/env python3
"""Exact positive memory of the chi-twisted multiples-of-3^j splits and the 2 + sqrt(2) amplification.

D_j(N) = sum_(n <= N, 3^j || n) chi(n/3^j) t(n), chi the nontrivial character mod 3, t(n) = (-1)^tau(n).
Since chi(2) = -1 and t(2m) = -t(m), every even summand equals the summand at m = n/2, so exactly
    D_j(N) = D_j(floor(N/2)) + O_j(N),   O_j(N) = sum over ODD n <= N with 3^j || n of chi(n/3^j) t(n),
and, since t(n) = -t((3n+1)/2) for odd n,
    O_j(N) = -[C(c_+ mod 3^(j+2); M) - C(c_- mod 3^(j+2); M)],  M = floor((3N+1)/2),
    c_+ = 2^(-1)(3^(j+1) + 1),  c_- = 2^(-1)(2*3^(j+1) + 1)  (mod 3^(j+2)).
If O_j(y) ~ o_j sqrt(y) then D_j(N) ~ o_j sqrt(N) / (1 - 2^(-1/2)) = (2 + sqrt 2) o_j sqrt(N).
This script checks both exact identities on N <= 2730 and, on log-uniform samples to 2^LOGMAX, reports the octave
means of D_j(N)/sqrt(N) and of O_j(N)/sqrt(N) and their ratio.  CAUTION: the ratio of the octave-averaged means is
2 + sqrt 2 up to end effects for ANY sequence obeying the exact identity, so it is a code check, not evidence for a
sqrt(N) main term (S37 audit; see the S37 note).
Session opus-2026-10-15-S37.  Usage: python collatz_twisted_split_memory_20261015.py LOGMAX out.json
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
OUT = sys.argv[2] if len(sys.argv) > 2 else 'twisted_split_memory.json'
J = 3


def chi(m):
    r = m % 3
    return 0 if r == 0 else (1 if r == 1 else -1)


@nb.njit
def sample(P, logmax, pts, M):
    # class sums mod M (= 3^(J+2)) at sorted sample points, plus D_j and O_j for j < J
    NMAX = 1 << logmax
    G = pts.shape[0]
    C = np.zeros(M, dtype=np.int64)
    D = np.zeros(4, dtype=np.int64)
    O = np.zeros(4, dtype=np.int64)
    outC = np.zeros((G, M), dtype=np.int64)
    outD = np.zeros((G, 4), dtype=np.int64)
    outO = np.zeros((G, 4), dtype=np.int64)
    r = 1
    gi = 0
    for N in range(1, NMAX + 1):
        t = 1 - 2 * int(cen.getbit(P, N))
        C[r] += t
        # 3-adic valuation and chi of the cofactor
        m = N
        v = 0
        while m % 3 == 0 and v < 4:
            m //= 3
            v += 1
        if v < 4:
            c = 1 if m % 3 == 1 else -1
            D[v] += c * t
            if N & 1:
                O[v] += c * t
        while gi < G and pts[gi] == N:
            for a in range(M):
                outC[gi, a] = C[a]
            for j in range(4):
                outD[gi, j] = D[j]
                outO[gi, j] = O[j]
            gi += 1
        r += 1
        if r == M:
            r = 0
    return outC, outD, outO


if __name__ == '__main__':
    t0 = time.time()
    P = cen.parity_bits(LOGMAX)
    M = 3 ** (J + 2)
    # exact identities on a full small range
    small = 1 << 20
    pts_small = np.arange(1, small + 1, dtype=np.int64)
    base = sorted(set(int(round(2 ** (10 + j / 32))) for j in range(0, (LOGMAX - 11) * 32 + 1)))
    base = [b for b in base if (3 * b + 1) // 2 <= (1 << LOGMAX)]
    extra = [(3 * b + 1) // 2 for b in base]
    pts = np.array(sorted(set(base) | set(extra) | set(range(1, 4097))), dtype=np.int64)
    Cg, Dg, Og = sample(P, LOGMAX, pts, M)
    idx = {int(p): i for i, p in enumerate(pts)}
    print(f'sampled {len(pts)} points to 2^{LOGMAX}: {time.time() - t0:.1f} s', flush=True)
    # exact checks on N <= 4096 using the dense prefix of the sample
    bad1 = bad2 = 0
    inv2 = {}
    for j in range(J):
        m = 3 ** (j + 2)
        i2 = pow(2, -1, m)
        cp = (i2 * (3 ** (j + 1) + 1)) % m
        cm = (i2 * (2 * 3 ** (j + 1) + 1)) % m
        inv2[j] = (cp, cm)
        for N in range(2, 2731):
            if N // 2 in idx and N in idx:
                if Dg[idx[N], j] != Dg[idx[N // 2], j] + Og[idx[N], j]:
                    bad1 += 1
    # direct exact check of O_j identity using class sums mod 3^(j+2) from mod 3^(J+2) sums
    for j in range(J):
        m = 3 ** (j + 2)
        cp, cm = inv2[j]
        for N in range(1, 2731):
            Mq = (3 * N + 1) // 2
            if Mq in idx and N in idx:
                cs = Cg[idx[Mq]]
                cpl = cs[cp::m].sum()
                cmi = cs[cm::m].sum()
                if Og[idx[N], j] != -(cpl - cmi):
                    bad2 += 1
    print(f'exact identity D_j(N) = D_j(floor(N/2)) + O_j(N) on N <= 2730: {bad1} failures; '
          f'O_j(N) = -[C(c+) - C(c-)](floor((3N+1)/2)): {bad2} failures', flush=True)
    res = dict(logmax=LOGMAX, classes={j: inv2[j] for j in range(J)}, exact_failures=[bad1, bad2], rows=[])
    ks = list(range(14, LOGMAX - 1))
    for j in range(J):
        dm = []
        om = []
        for k in ks:
            sel = [b for b in base if (1 << k) <= b < (1 << (k + 1))]
            dm.append(np.mean([Dg[idx[b], j] / math.sqrt(b) for b in sel]))
            om.append(np.mean([Og[idx[b], j] / math.sqrt(b) for b in sel]))
        dm = np.array(dm)
        om = np.array(om)
        ratio = dm.mean() / om.mean()
        res['rows'].append(dict(j=j, D_mean=float(dm.mean()), O_mean=float(om.mean()), ratio=float(ratio),
                                prediction=2 + math.sqrt(2), D_octaves=[float(x) for x in dm], O_octaves=[float(x) for x in om]))
        print(f'j = {j}: mean D_j/sqrt(N) {dm.mean():+.4f}, mean O_j/sqrt(N) {om.mean():+.4f}, ratio {ratio:.3f} '
              f'(2 + sqrt 2 = {2 + math.sqrt(2):.4f}); O_j positive octaves {(om > 0).sum()}/{len(ks)}', flush=True)
    json.dump(res, open(OUT, 'w'), indent=1)
