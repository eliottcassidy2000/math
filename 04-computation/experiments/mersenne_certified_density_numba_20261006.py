"""Certified 2-adic Mersenne switching density to larger template totals K (numba, parallel).
Session opus-2026-10-06-S18; note 05-knowledge/results/mersenne_switch_parity_f21_compression_20261006.md, section 4.
Extends THM-4556 (v) (K <= 20, share 0.1199) to K = 27.  Reproduce: python3 <this> 9 27
Odd a <-> X = 3^(a-1) mod 2^(K+1), X = 1 mod 8.  n's post-run start p = 2X - 1; partner (2^(a-D)-1) post-run start
q = 2X 3^(-D) - 1 (D odd).  Collision: prefixes with equal total, length difference D, equal numerator."""
import sys, time
import numpy as np
import numba as nb

KMIN = int(sys.argv[1]) if len(sys.argv) > 1 else 9
KMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 24


@nb.njit(cache=True)
def word_vals(x, M, totals, vals):
    n = 0
    cur = 0
    bits = M
    x &= (1 << M) - 1
    X = -1
    e = 0
    while bits > 0:
        if x & 1:
            if cur > 0:
                X = 3 * X + (1 << e)
                e += cur
                totals[n] = e
                vals[n] = X
                n += 1
            cur = 1
            x = (3 * x + 1) >> 1
        else:
            cur += 1
            x >>= 1
        bits -= 1
        if bits > 0:
            x &= (1 << bits) - 1
        else:
            x = 0
    return n


@nb.njit(parallel=True, cache=True)
def count_hits(K, inv3pows):
    M = K + 1
    mod = 1 << M
    nclass = mod // 8
    hits = np.zeros(nclass, dtype=np.uint8)
    for c in nb.prange(nclass):
        Xr = 1 + 8 * c
        tu = np.empty(64, dtype=np.int64)
        vu = np.empty(64, dtype=np.int64)
        tq = np.empty(64, dtype=np.int64)
        vq = np.empty(64, dtype=np.int64)
        p = (2 * Xr - 1) & (mod - 1)
        nu = word_vals(p, M, tu, vu)
        hit = 0
        D = 1
        while D < K and hit == 0:
            q = ((2 * Xr * inv3pows[D]) - 1) & (mod - 1)
            nq = word_vals(q, M, tq, vq)
            for i in range(nu):
                j = i + D
                if j < nq and tq[j] == tu[i] and vq[j] == vu[i]:
                    hit = 1
                    break
            D += 2
        hits[c] = hit
    return hits.sum(), nclass


for K in range(KMIN, KMAX + 1):
    M = K + 1
    mod = 1 << M
    inv3 = pow(3, -1, mod)
    inv3pows = np.array([pow(inv3, D, mod) for D in range(K + 2)], dtype=np.int64)
    t0 = time.time()
    h, ncl = count_hits(K, inv3pows)
    print(f'K = {K:2d}: certified share {h}/{ncl} = {h / ncl:.4f}   ({time.time() - t0:.1f}s)', flush=True)
