#!/usr/bin/env python3
"""A6 basin densities by sampling (S56): uniform random sweep starts W in the depth-d window, plus a neighbourhood
of a given landing point; exact fate before leaving the window (HALT / CAPTURE / exit), certified accelerator.

Usage: python a6_basin_sample_20261025.py d_min d_max n_samples [workers] [landing_file.json]
Output: per depth the halting and capture fractions with binomial standard errors, by window decile.
"""
import sys, os, json, time, random, math
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def T(d):
    return 24 * (1 << d) - 3 * d - 16


def fate(args):
    d, W = args
    import a6_accel_20261025 as acc
    from a6_fast_20261025 import depth_of
    st = acc.orbit(W, max_sweeps=10 ** 18, stop=lambda x, n: depth_of(0, x) > d)
    if st[0] == 'stopped':
        return d, W, 'exit', st[1] - T(d)
    return d, W, st[0], st[1]


if __name__ == '__main__':
    dmin, dmax, ns = (int(x) for x in sys.argv[1:4])
    workers = int(sys.argv[4]) if len(sys.argv) > 4 else 10
    out = {}
    with Pool(workers) as pool:
        for d in range(dmin, dmax + 1):
            rng = random.Random(1000 + d)
            lo, hi = T(d - 1) + 1, T(d)
            Ws = [rng.randrange(lo, hi + 1) for _ in range(ns)]
            t0 = time.time()
            rows = pool.map(fate, [(d, W) for W in Ws], chunksize=8)
            nH = sum(1 for r in rows if r[2] == 'HALT')
            nC = sum(1 for r in rows if r[2] == 'CAPTURE')
            exits = {}
            for r in rows:
                if r[2] == 'exit':
                    exits[r[3]] = exits.get(r[3], 0) + 1
            dec = {}
            for r in rows:
                k = int(10 * (r[1] - lo) / (hi - lo + 1))
                dec.setdefault(k, [0, 0, 0])
                dec[k][0] += 1
                dec[k][1] += r[2] == 'HALT'
                dec[k][2] += r[2] == 'CAPTURE'
            se = lambda k: math.sqrt(max(k, 1) * (ns - k) / ns) / ns
            out[d] = {'n': ns, 'halt': nH, 'capture': nC, 'distinct_exits': len(exits),
                      'top_exits': sorted(exits.items(), key=lambda kv: -kv[1])[:8],
                      'deciles': dec}
            print(f'd={d} ({time.time() - t0:.0f}s): halt {nH}/{ns} = {nH / ns:.4f}+-{se(nH):.4f}, capture {nC}/{ns} = '
                  f'{nC / ns:.4f}+-{se(nC):.4f}, distinct exits {len(exits)}, top '
                  f'{sorted(exits.items(), key=lambda kv: -kv[1])[:4]}; deciles (n,H,C): '
                  + ' '.join(f'{k}:{v[1]}/{v[2]}' for k, v in sorted(dec.items())), flush=True)
            json.dump(out, open(os.path.join(HERE, f'a6_basin_sample_{dmin}_{dmax}_{ns}.json'), 'w'), indent=1,
                      default=str)
