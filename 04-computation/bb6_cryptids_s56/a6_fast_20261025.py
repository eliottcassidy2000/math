#!/usr/bin/env python3
"""A6 exact accelerated simulator: phase map + jump over P(8, .) e0ueven runs (S56).

Phase P(v0, w0) = J(0, v0, w0); see a6_phase_map_20261025.py for the derivation of the end cases.
Run acceleration (proved in the session note): for v0 = 8 at depth j = 0 mod 4 the end case e0ueven maps
P(8, w0) -> P(8, w0 + 3j + 20) and keeps eps = 0 and u even, so the whole run at depth j is jumped in O(1).
sweep(W) returns (W', signature) or ('HALT', signature) or raises Infinite after `cap` phases.
"""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from a6_phase_map_20261025 import phase


class TooManyPhases(Exception):
    pass


def depth_of(v0, w0):
    A = v0 + 6
    j = 0
    if w0 + 2 * A + 3 >= 4 * A:
        j = max(0, (w0 + 2 * A + 3) // (4 * A)).bit_length() - 1
        while j > 0 and w0 + 2 * A + 3 * (j - 1) + 3 < 4 * A * (1 << (j - 1)):
            j -= 1
        while not (w0 + 2 * A + 3 * j + 3 < 4 * A * (1 << j)):
            j += 1
    return j


def sweep(W, cap=10 ** 6, record=False):
    v0, w0 = 0, W
    sig = [] if record else None
    n = 0
    while True:
        n += 1
        if n > cap:
            raise TooManyPhases((v0, w0))
        if v0 == 8:
            j = depth_of(8, w0)
            if j % 4 == 0 and j > 0:
                # e0ueven iff w0 + 3j odd and u even; check by one phase evaluation
                res, jj, lab = phase(8, w0)
                if lab == 'e0ueven':
                    delta = 3 * j + 20
                    kstar = (56 * (1 << j) - 32 - 3 * j - w0) // delta
                    if kstar >= 1:
                        if record:
                            sig.append(('run', j, kstar + 1))
                        w0 = w0 + (kstar + 1) * delta
                        continue
        res, j, lab = phase(v0, w0)
        if record:
            sig.append((v0 if v0 in (0, 8) else 'x', j, lab))
        if res[0] == 'HALT':
            return 'HALT', sig
        if res[0] == 'ret':
            return res[1], sig
        v0, w0 = res[1], res[2]


if __name__ == '__main__':
    # self-test against the plain phase map on all W < 20000 and on the orbit's first 20000 sweeps
    from a6_phase_map_20261025 import sweep_phase
    def fast_or_inf(W):
        try:
            return sweep(W, cap=20000)[0]
        except TooManyPhases:
            return 'INF'

    def plain_or_inf(W):
        try:
            return sweep_phase(W)[0]
        except RuntimeError:
            return 'INF'
    bad = []
    infs = []
    for W in range(20000):
        a = fast_or_inf(W)
        if a == 'INF':
            infs.append(W)
            continue
        if a != plain_or_inf(W):
            bad.append(W)
    print('fast sweep == phase map for W < 20000 (returning sweeps):', not bad, bad[:5])
    print('non-returning (>20000 phases) W < 20000:', len(infs), infs[:30])
    W = 0
    ok = True
    for i in range(20000):
        a = sweep(W)[0]
        b = sweep_phase(W)[0]
        if a != b:
            ok = False
            print('mismatch at', W)
            break
        W = a
    print('orbit first 20000 sweeps agree:', ok, 'W =', W)
