#!/usr/bin/env python3
"""A6 exact basins inside one depth window: for EVERY sweep start W of depth d, its fate before leaving the window
(HALT, CAPTURE, or the exit overshoot past T_d), by memoised backward processing (S56).

Window of depth d: T_{d-1} < W <= T_d.  fate(W) = fate(S(W)) when S(W) stays in the window and S(W) > W
(checked; decreasing steps are handled by a second pass).  Prints basin fractions and the number of distinct exits.
Usage: python a6_basins_20261025.py d_min d_max
"""
import sys, os, json, time
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import a6_accel_20261025 as acc


def T(d):
    return 24 * (1 << d) - 3 * d - 16


def basins(d):
    lo, hi = T(d - 1) + 1, T(d)
    n = hi - lo + 1
    fate = [None] * n          # int >= 0: exit overshoot; -1 HALT; -2 CAPTURE; -3 pending (S(W) < W)
    succ = {}
    back = 0
    for W in range(hi, lo - 1, -1):
        i = W - lo
        try:
            W2, _ = acc.sweep(W)
        except acc.Halt:
            fate[i] = -1
            continue
        except acc.Capture:
            fate[i] = -2
            continue
        if W2 > hi:
            fate[i] = W2 - hi
        elif W2 > W:
            fate[i] = fate[W2 - lo]
        else:
            back += 1
            succ[i] = W2 - lo
            fate[i] = -3
    # resolve backward steps (W2 <= W) by iteration
    changed = True
    while changed:
        changed = False
        for i, j in list(succ.items()):
            if j < 0:
                continue
            if fate[j] != -3:
                fate[i] = fate[j]
                del succ[i]
                changed = True
    return fate, back, len(succ)


if __name__ == '__main__':
    dmin, dmax = int(sys.argv[1]), int(sys.argv[2])
    out = {}
    for d in range(dmin, dmax + 1):
        t0 = time.time()
        fate, back, unresolved = basins(d)
        n = len(fate)
        c = Counter(fate)
        nH, nC = c.get(-1, 0), c.get(-2, 0)
        exits = Counter(f for f in fate if f is not None and f >= 0)
        top = exits.most_common(6)
        # basin of halting among the first 400 entries (the natural entry range)
        e400 = Counter(fate[:400])
        out[d] = {'n': n, 'halt_basin': nH, 'capture_basin': nC, 'backsteps': back, 'unresolved': unresolved,
                  'distinct_exits': len(exits), 'top_exits': top, 'entries400': sorted(e400.items())[:10]}
        print(f'd={d} window {n}: halt basin {nH} ({nH / n:.2e}), capture basin {nC} ({nC / n:.2e}), '
              f'back-steps {back}, unresolved {unresolved}, distinct exits {len(exits)}, top {top[:4]} '
              f'({time.time() - t0:.0f}s)', flush=True)
    json.dump(out, open(os.path.join(HERE, f'a6_basins_{dmin}_{dmax}.json'), 'w'), indent=1, default=str)
