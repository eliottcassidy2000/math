#!/usr/bin/env python3
"""A6: exhaustive halting set H and capture set C (certified non-returning sweeps) for W < 2^E (S56).

For every sweep start W < 2^E: HALT, CAPTURE (certified translated cycler, three-point lemma) or return W'.
We print per-depth counts and fit the families {alpha*2^d + beta*d + gamma} that contain the elements.
Usage: python a6_target_sets_20261025.py [E]
"""
import sys, os, json
from collections import defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import a6_accel_20261025 as acc
from a6_fast_20261025 import depth_of

if __name__ == '__main__':
    E = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    H, C = [], []
    certs = {}
    for W in range(1 << E):
        try:
            acc.sweep(W)
        except acc.Halt:
            H.append(W)
        except acc.Capture as c:
            C.append(W)
            b, D, cyc, e = c.args[0]
            certs[W] = (b, D, len(cyc), [l for l in cyc])
        except RuntimeError as r:
            print('UNRESOLVED (cap) at', W, r.args, flush=True)
    byd = defaultdict(lambda: [[], []])
    for W in H:
        byd[depth_of(0, W)][0].append(W)
    for W in C:
        byd[depth_of(0, W)][1].append(W)
    print('depth: #H #C | H elements as (W - 24*2^d) | C elements as (W - 18*2^d)')
    for d in sorted(byd):
        h, c = byd[d]
        print(d, len(h), len(c), '| H:', [W - 24 * (1 << d) for W in h], '| C:', [W - 18 * (1 << d) for W in c])
    # affine-exponential families: group by (W - a*2^d) linear in d for a in a small set
    print('\ncapture certificates (first 30): W -> (base state, translation, period length)')
    for W in C[:30]:
        b, D, p, cyc = certs[W]
        print(W, b, D, p, cyc)
    json.dump({'E': E, 'H': H, 'C': C, 'certs': {str(k): [list(v[0]), list(v[1]), v[2], [list(x) for x in v[3]]] for k, v in certs.items()}},
              open(os.path.join(HERE, f'a6_target_sets_E{E}.json'), 'w'), indent=1, default=str)
