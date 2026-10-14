#!/usr/bin/env python3
"""Cross-check of the multi-run W_4 coefficients: the symbolic results P(q) + C Sigma(W_4) for runs (2,6), (4,4),
(2,2,4), (2,2,2,2) are evaluated at q = 7 and 11 and compared with brute-force DFS over injective assignments
(the sums must vanish whenever the configuration needs more distinct points than q)."""
import sys, time
import numpy as np
import numba as nb
sys.path.insert(0, '.')
from audit_A10_cores import injective_sum_paths, core_sum, chi_table
from audit_run_polys_symbolic import peval


@nb.njit(cache=True)
def brute_runs(q, chi, starts, npts):
    """sum over injective v: [npts] -> F_q of prod over path edges (i -> i+1 inside each run)"""
    used = np.zeros(q, dtype=np.int64)
    seq = np.zeros(npts, dtype=np.int64)
    nxt = np.zeros(npts + 1, dtype=np.int64)
    isstart = np.zeros(npts, dtype=np.int64)
    for s in starts:
        isstart[s] = 1
    sgn = np.zeros(npts + 1, dtype=np.int64)
    sgn[0] = 1
    tot = 0
    depth = 0
    nxt[0] = 0
    while depth >= 0:
        if depth == npts:
            tot += sgn[npts]
            depth -= 1
            used[seq[depth]] = 0
            nxt[depth] = seq[depth] + 1
            continue
        x = nxt[depth]
        while x < q and used[x] == 1:
            x += 1
        if x >= q:
            depth -= 1
            if depth >= 0:
                used[seq[depth]] = 0
                nxt[depth] = seq[depth] + 1
            continue
        seq[depth] = x
        used[x] = 1
        if isstart[depth] == 1:
            sgn[depth + 1] = sgn[depth]
        else:
            sgn[depth + 1] = sgn[depth] * chi[(x - seq[depth - 1]) % q]
        depth += 1
        nxt[depth] = 0
    return tot


def main():
    for lens in ([2, 6], [4, 4], [2, 2, 4], [2, 2, 2, 2]):
        res, _ = injective_sum_paths(lens)
        starts, off = [], 0
        for L in lens:
            starts.append(off)
            off += L + 1
        for q in (7, 11):
            chi = chi_table(q)
            val = peval(res.get('P', {}), q)
            for k, v in res.items():
                if k == 'P':
                    continue
                m = 1 + max(max(e) for e in k)
                s = int(core_sum(q, chi, m, np.array([e[0] for e in k], dtype=np.int64), np.array([e[1] for e in k], dtype=np.int64)))
                val += peval(v, q) * s
            t0 = time.time()
            b = int(brute_runs(q, chi, np.array(starts, dtype=np.int64), off)) if off <= q else 0
            print(f'runs {lens} ({off} points) at q = {q}: symbolic P + C Sigma(W_4) = {val};  brute force = {b};  agree {val == b} ({time.time() - t0:.1f} s)', flush=True)


if __name__ == '__main__':
    main()
