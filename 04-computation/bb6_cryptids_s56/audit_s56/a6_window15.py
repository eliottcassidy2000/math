#!/usr/bin/env python3
"""Referee: exhaustive fates of all sweep starts in the depth-15 window with the referee's own sweep map
(memoised forward-to-exit), locating halting and captured basins relative to T_14 and the orbit landing e_15."""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from a6_my_fast import my_fast_sweep
T = lambda d: 24 * 2 ** d - 3 * d - 16
d = int(sys.argv[1]) if len(sys.argv) > 1 else 15
lo, hi = T(d - 1) + 1, T(d)
fate = {}
nonret = []
for W in range(hi, lo - 1, -1):
    try:
        W2 = my_fast_sweep(W, cap=200000)
    except RuntimeError:
        fate[W] = 'NR'        # non-returning within cap (candidate capture)
        nonret.append(W)
        continue
    if W2 == 'HALT':
        fate[W] = 'H'
    elif W2 > hi:
        fate[W] = ('X', W2 - hi)
    elif W2 > W:
        fate[W] = fate[W2]
    else:
        fate[W] = ('BACK', W2)
H = sorted(W - T(d - 1) for W, f in fate.items() if f == 'H')
C = sorted(W - T(d - 1) for W, f in fate.items() if f == 'NR')
print(f'd={d}: window size {hi - lo + 1}; halting basin {len(H)}; non-returning basin {len(C)}; back-steps',
      sum(1 for f in fate.values() if isinstance(f, tuple) and f[0] == 'BACK'))
def ranges(xs):
    out = []
    for x in xs:
        if out and x - out[-1][1] <= 60:
            out[-1][1] = x; out[-1][2] += 1
        else:
            out.append([x, x, 1])
    return out
print('halting basin offsets (clusters [first,last,count]):', ranges(H)[:20])
print('non-returning basin offsets (clusters):', ranges(C)[:20])
