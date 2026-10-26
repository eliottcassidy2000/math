#!/usr/bin/env python3
"""A6 per-depth hazard: from every entry offset just past the top of depth d-1, run the exact accelerated
dynamics until the orbit leaves depth d; count entries that HALT or are CAPTURED inside the window (S56).

Entry points W = T_{d-1} + o for o = 1..O (T_{d-1} = 24*2^(d-1) - 3(d-1) - 16 is the halting top of depth d-1).
The observed orbit enters depth d with overshoot o mostly in [1, ~300].  h_T(d), h_N(d) = fractions halting /
captured; their decay in d decides whether the expected number of target hits diverges (sum ~ c/d) or not.
Usage: python a6_window_hazard_20261025.py d_min d_max O
"""
import sys, os, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import a6_accel_20261025 as acc
from a6_fast_20261025 import depth_of


def traverse(W, d, max_sweeps=10 ** 15):
    """Run from W until depth > d. Returns ('HALT'|'CAPTURE'|'exit', W_exit, sweeps)."""
    st = acc.orbit(W, max_sweeps=max_sweeps, stop=lambda x, n: depth_of(0, x) > d)
    if st[0] == 'stopped':
        return 'exit', st[1], st[2]
    return st[0], st[1], st[2]


if __name__ == '__main__':
    dmin, dmax, O = (int(x) for x in sys.argv[1:4])
    res = {}
    for d in range(dmin, dmax + 1):
        T_prev = 24 * (1 << (d - 1)) - 3 * (d - 1) - 16
        t0 = time.time()
        cnt = {'HALT': [], 'CAPTURE': [], 'exit': 0, 'other': []}
        exits = {}
        for o in range(1, O + 1):
            W = T_prev + o
            assert depth_of(0, W) == d
            kind, Wx, n = traverse(W, d)
            if kind == 'exit':
                cnt['exit'] += 1
                exits[o] = Wx - (24 * (1 << d) - 3 * d - 16)
            elif kind in ('HALT', 'CAPTURE'):
                cnt[kind].append(o)
            else:
                cnt['other'].append((o, kind))
        hT, hN = len(cnt['HALT']) / O, len(cnt['CAPTURE']) / O
        res[d] = {'O': O, 'halt_offsets': cnt['HALT'], 'capture_offsets': cnt['CAPTURE'], 'other': cnt['other'],
                  'exit_overshoots_sample': sorted(set(exits.values()))[:40]}
        print(f'd={d}: h_T={hT:.4f} ({len(cnt["HALT"])}) h_N={hN:.4f} ({len(cnt["CAPTURE"])}) other={cnt["other"][:3]} '
              f'halt offsets {cnt["HALT"][:12]} capture offsets {cnt["CAPTURE"][:12]} ({time.time() - t0:.0f}s)',
              flush=True)
    json.dump(res, open(os.path.join(HERE, f'a6_window_hazard_{dmin}_{dmax}_{O}.json'), 'w'), indent=1, default=str)
