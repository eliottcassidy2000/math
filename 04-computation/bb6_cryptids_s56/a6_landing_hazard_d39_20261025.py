#!/usr/bin/env python3
"""A6 hazard at the big-jump landings (S56): at depths d = 3 (mod 4) the orbit enters depth d by a long jump from
depth d-1 to offset e_d past T_{d-1} (e_d ~ (4/9) 2^d).  Fates (HALT / CAPTURE / exit) of the entries
T_{d-1} + e_d + delta, |delta| <= R, measure the basin densities around the actual landing point.
Usage: python a6_landing_hazard_20261025.py R [workers]
"""
import sys, os, json, time
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from a6_transfer_20261025 import T

LANDING = {39: 244335917530}


def job(args):
    d, W = args
    import a6_accel_20261025 as acc
    from a6_fast_20261025 import depth_of
    st = acc.orbit(W, max_sweeps=10 ** 18, stop=lambda x, n: depth_of(0, x) > d)
    if st[0] == 'stopped':
        return d, W, 'exit', st[1] - T(d)
    return d, W, st[0], st[1]


if __name__ == '__main__':
    R = int(sys.argv[1])
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    out = {}
    with Pool(workers) as pool:
        for d, e in sorted(LANDING.items()):
            base = T(d - 1) + e
            t0 = time.time()
            rows = pool.map(job, [(d, base + k) for k in range(-R, R + 1)], chunksize=8)
            nH = [r[1] - base for r in rows if r[2] == 'HALT']
            nC = [r[1] - base for r in rows if r[2] == 'CAPTURE']
            ex = {}
            for r in rows:
                if r[2] == 'exit':
                    ex[r[3]] = ex.get(r[3], 0) + 1
            actual = [r for r in rows if r[1] == base][0]
            out[d] = {'landing_offset': e, 'R': R, 'halts': nH, 'captures': nC, 'exits': sorted(ex.items(), key=lambda kv: -kv[1])[:10],
                      'actual': actual[2:]}
            print(f'd={d} landing e_d={e} ({time.time() - t0:.0f}s): halt {len(nH)}/{2 * R + 1}, capture {len(nC)}/{2 * R + 1} '
                  f'(deltas {nC[:6]}), distinct exits {len(ex)} top {sorted(ex.items(), key=lambda kv: -kv[1])[:4]}; '
                  f'actual entry fate {actual[2:]}', flush=True)
            json.dump(out, open(os.path.join(HERE, f'a6_landing_hazard_d39_R{R}.json'), 'w'), indent=1, default=str)
