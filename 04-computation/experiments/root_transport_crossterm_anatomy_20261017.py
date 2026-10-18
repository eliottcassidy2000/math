#!/usr/bin/env python3
"""Anatomy of the ideal energy cross terms (S41): which merging pairs dissipate the transport energy?

In the ideal recursion K_d(r, b) = sum_j v_j with v_j = sigma_j K_(d-1)(s_j(r), b - o_j), the band-weighted energy is
E_d = E_d^diag + E_d^cross with E_d^diag = E_(d-1) exactly (diagonal isometry) and
E_d^cross = sum_b 2^-b sum_r sum_(j != j') v_j v_j' (normalised by 4 2^d (1 - 2^-P_d)). This script splits the cross
terms by the relation between the two fibre sources: 3-adic siblings (s_j = s_j' mod 3^(d-1), different last digit:
the zero-sum lifts), equal sources (s_j = s_j' mod 3^d) and unrelated; and by the lag j' - j.
Usage: python root_transport_crossterm_anatomy_20261017.py DMAX
Session opus-2026-10-17-S41.
"""
import sys, os, json
import numpy as np
from numba import njit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from root_transport_depth_ladder_20261017 import K0, step


@njit(cache=True)
def anatomy(prev, d):
    mod = 3 ** (d + 1); period = 2 * mod
    oldmod = mod // 3; oldperiod = 2 * oldmod; sib = oldmod // 3 if d >= 1 else 1
    diag = 0.0; c_sib = 0.0; c_eq = 0.0; c_other = 0.0
    lagsum = np.zeros(64)
    js = np.zeros(period, dtype=np.int64); ss = np.zeros(period, dtype=np.int64); vs = np.zeros(period, dtype=np.int64)
    for r in range(mod):
        for b in range(period):
            n = 0
            rho = r; j = 0
            while 2 * j <= b:
                if rho % 3 != 0:
                    a = 2 if rho % 3 == 1 else 1
                    off = (a - 1) + 2 * j
                    if off <= b:
                        s = ((2 ** a * rho - 1) // 3) % oldmod
                        v = (1 if a == 2 else -1) * prev[(b - off) % oldperiod, s]
                        if v != 0:
                            js[n] = j; ss[n] = s; vs[n] = v; n += 1
                rho = (4 * rho + 1) % mod; j += 1
            w = 0.5 ** b
            for i in range(n):
                diag += w * vs[i] * vs[i]
                for k in range(i + 1, n):
                    t = 2.0 * w * vs[i] * vs[k]
                    if ss[i] == ss[k]:
                        c_eq += t
                    elif sib >= 1 and ss[i] % sib == ss[k] % sib:
                        c_sib += t
                    else:
                        c_other += t
                    lag = js[k] - js[i]
                    if lag < 64:
                        lagsum[lag] += t
    norm = 4.0 * 2.0 ** d * (1.0 - 2.0 ** (-period))
    return diag / norm, c_sib / norm, c_eq / norm, c_other / norm, lagsum / norm


def main():
    DMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    K = K0.copy(); rows = []
    for d in range(1, DMAX + 1):
        diag, csib, ceq, coth, lag = anatomy(K, d)
        K = step(K, d)
        print('d = %d: diagonal %.6f, cross: siblings %.6f, equal %.6f, other %.6f (total cross %.6f); lag profile (1..8): %s'
              % (d, diag, csib, ceq, coth, csib + ceq + coth, ' '.join('%.4f' % x for x in lag[1:9])), flush=True)
        rows.append(dict(d=d, diag=diag, cross_siblings=csib, cross_equal=ceq, cross_other=coth, lag=lag[:16].tolist()))
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(rows, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'root_transport_crossterm_anatomy_20261017.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
