#!/usr/bin/env python3
"""Referee: plain A6 orbit with the referee's phase map; record (sweep index, W) of the first orbit point at each
depth, and near-miss gaps to T_d, N_d (exact, every sweep), for comparison with the session's accelerated run."""
import sys, os, time, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from a6_my_fast import my_fast_sweep as _mfs
def my_sweep(W): return (_mfs(W),)


def depth(W):
    j = 0
    while not (W + 15 + 3 * j < 24 << j):
        j += 1
    return j


if __name__ == '__main__':
    N = int(float(sys.argv[1]))
    W = 0
    cur = 0
    t0 = time.time()
    entries = {0: (0, 0)}
    gaps = {}
    T = lambda d: 24 * 2 ** d - 3 * d - 16
    Nd = lambda d: 18 * 2 ** d - 3 * d - 4
    i = 0
    dT = T(cur)
    while i < N:
        W2 = my_sweep(W)[0]
        i += 1
        if W2 == 'HALT':
            print('HALT', i, W)
            break
        # near misses (straddling steps) for T and N of all depths crossed
        for d in range(max(1, depth(W) - 1), depth(W2) + 2):
            for name, X in (('T', T(d)), ('N', Nd(d))):
                if W <= X < W2 or W < X <= W2:
                    gaps[(name, d)] = (X - W, W2 - X)
        W = W2
        if W > dT:
            d = depth(W)
            for dd in range(cur + 1, d + 1):
                entries[dd] = (i, W)
            cur = d
            dT = T(cur)
            print(f'depth {cur} entry sweep {i} W {W} e_d {W - T(cur - 1)} ({time.time() - t0:.0f}s)', flush=True)
    json.dump({'entries': entries, 'gaps': {f'{k[0]}{k[1]}': v for k, v in gaps.items()}, 'sweeps': i, 'W': W},
              open(os.path.join(HERE, 'a6_entries.json'), 'w'), indent=1)
    print('done', i, W)
