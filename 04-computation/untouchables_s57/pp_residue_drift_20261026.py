#!/usr/bin/env python3
"""Residue-class test of the Pollack-Pomerance class model (opus-2026-10-26-S57).

PP's model (Trans. AMS Ser. B 3 (2016), sec. 3.1) gives an even target k an untouchability probability that depends on
k only through a = gcd(k, A_y), A_y = lcm(1..y).  For a prime l the residues r = 1..l-1 of k mod l therefore receive
the SAME predicted share (by CRT the other coordinates of a are independent of which nonzero residue k has mod l).
We measure, for even k in dyadic windows (2^(e-1), 2^e], the untouchable share in each residue class mod l, l = 3..13,
and report the spread among the nonzero residues (max/min) and its drift with e.
Usage: python pp_residue_drift_20261026.py N bigdir [emin]
"""
import sys, os, json
import numpy as np
import numba as nb

HERE = os.path.dirname(os.path.abspath(__file__))
LS = np.array([3, 5, 7, 11, 13], np.int64)


@nb.njit(cache=False)
def tally(mult, lo_k, hi_k, LS, tot, unt):
    # even k in (lo_k, hi_k]
    j0 = lo_k // 2 + 1
    j1 = hi_k // 2
    for j in range(j0, j1 + 1):
        k = 2 * j
        z = 1 if mult[j] == 0 else 0
        for t in range(LS.shape[0]):
            r = k % LS[t]
            tot[t, r] += 1
            unt[t, r] += z


if __name__ == '__main__':
    N = int(float(sys.argv[1]))
    mult = np.load(os.path.join(sys.argv[2], f'mult_even_{N}.npy'), mmap_mode='r')
    emin = int(sys.argv[3]) if len(sys.argv) > 3 else 14
    e = emin
    rows = []
    while 2 ** e <= N:
        tot = np.zeros((LS.shape[0], 13), np.int64)
        unt = np.zeros((LS.shape[0], 13), np.int64)
        tally(mult, 2 ** (e - 1), 2 ** e, LS, tot, unt)
        row = {'e': e}
        line = f'2^{e:>2}:'
        for t, l in enumerate(LS):
            sh = [unt[t, r] / tot[t, r] for r in range(l)]
            nz = sh[1:]
            row[int(l)] = {'shares': [round(x, 5) for x in sh], 'nonzero_max_over_min': max(nz) / min(nz),
                           'argmax': 1 + int(np.argmax(nz)), 'argmin': 1 + int(np.argmin(nz))}
            line += f'  mod {l}: ' + ' '.join(f'{x:.3f}' for x in sh) + f' (nz max/min {max(nz) / min(nz):.3f})'
        rows.append(row)
        print(line, flush=True)
        e += 1
    json.dump({'N': N, 'rows': rows}, open(os.path.join(HERE, f'pp_residue_drift_{N}.json'), 'w'), indent=1)
