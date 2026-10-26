#!/usr/bin/env python3
"""A6 depth-transfer maps E_d: entry overshoot o (W = T_{d-1} + o) -> exit overshoot o' (first W past T_d) (S56).

Exact (certified accelerator).  Outcomes HALT / CAPTURE are recorded.  Parallel over entries.
Usage: python a6_transfer_20261025.py d_min d_max o_max [workers]
"""
import sys, os, json, time
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def T(d):
    return 24 * (1 << d) - 3 * d - 16


def job(args):
    d, o = args
    import a6_accel_20261025 as acc
    from a6_fast_20261025 import depth_of
    W = T(d - 1) + o
    st = acc.orbit(W, max_sweeps=10 ** 18, stop=lambda x, n: depth_of(0, x) > d)
    if st[0] == 'stopped':
        return d, o, 'exit', st[1] - T(d), st[2]
    return d, o, st[0], st[1], st[2]


if __name__ == '__main__':
    dmin, dmax, omax = (int(x) for x in sys.argv[1:4])
    workers = int(sys.argv[4]) if len(sys.argv) > 4 else 10
    out = {}
    with Pool(workers) as pool:
        for d in range(dmin, dmax + 1):
            t0 = time.time()
            rows = pool.map(job, [(d, o) for o in range(1, omax + 1)], chunksize=4)
            ex = {}
            halts, caps = [], []
            for _, o, kind, val, n in rows:
                if kind == 'exit':
                    ex.setdefault(val, []).append(o)
                elif kind == 'HALT':
                    halts.append(o)
                elif kind == 'CAPTURE':
                    caps.append((o, val))
            img = sorted(ex.items(), key=lambda kv: -len(kv[1]))
            out[d] = {'image': [(k, len(v), v[:6]) for k, v in img], 'halts': halts, 'captures': caps}
            print(f'd={d} ({time.time() - t0:.0f}s): {len(img)} distinct exits; top: '
                  + ', '.join(f'{k}x{len(v)}' for k, v in img[:6])
                  + f' | halts {len(halts)} {halts[:8]} | captures {len(caps)} {sorted(set(w for _, w in caps))[:4]}',
                  flush=True)
            json.dump(out, open(os.path.join(HERE, f'a6_transfer_{dmin}_{dmax}_{omax}.json'), 'w'), indent=1,
                      default=str)
