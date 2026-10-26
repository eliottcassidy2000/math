#!/usr/bin/env python3
"""Drift of the untouchable share per abundancy bin x = k/s(k) across scales (opus-2026-10-26-S57).
PP's model, read per number (class of k -> sigma(k)/k), predicts share exp(-x) in each bin, for every large scale.
Windows (2^(e-1), 2^e], e = 16, 18, ..., 30 (last capped at 1e9).
Usage: python pp_xbin_drift_20261026.py bigdir
"""
import sys, os, json, math
import numpy as np
import numba as nb
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from untouchable_sieve_20261026 import primes_upto
from untouchable_profile_20261026 import s_even_segment

@nb.njit(cache=False)
def binit(lo, cnt, sv, mult, nb_, tot, unt, sumexp):
    for i in range(cnt):
        k = lo + 2 * i
        x = k / sv[i]
        b = int(x * nb_ / 2.0)
        if b >= nb_:
            b = nb_ - 1
        tot[b] += 1
        sumexp[b] += math.exp(-x)
        if mult[k // 2] == 0:
            unt[b] += 1

if __name__ == '__main__':
    mult = np.load(os.path.join(sys.argv[1], 'mult_even_1000000000.npy'))
    primes = primes_upto(int(math.isqrt(2 * 10 ** 9)) + 2)
    nb_ = 10
    res = {}
    for e in range(16, 31, 2):
        lo, hi = 2 ** (e - 1) + 2, min(2 ** e, 10 ** 9)
        tot = np.zeros(nb_, np.int64); unt = np.zeros(nb_, np.int64); se = np.zeros(nb_)
        L = lo
        while L <= hi:
            cnt = min(1 << 22, (hi - L) // 2 + 1)
            sv = np.empty(cnt, np.int64)
            s_even_segment(L, cnt, primes, sv)
            binit(L, cnt, sv, mult, nb_, tot, unt, se)
            L += 2 * cnt
        res[e] = [(int(tot[b]), unt[b] / tot[b] if tot[b] else None, se[b] / tot[b] if tot[b] else None) for b in range(nb_)]
        print(f'2^{e}: ' + '  '.join(f'[{2*b/nb_:.1f},{2*(b+1)/nb_:.1f}) {r[1]:.3f}/{r[2]:.3f}' for b, r in enumerate(res[e]) if r[0] > 500), flush=True)
    print('(each entry: empirical untouchable share / PP mean exp(-x) in the bin)')
    json.dump({str(k): v for k, v in res.items()}, open(os.path.join(HERE, 'pp_xbin_drift.json'), 'w'), indent=1)
