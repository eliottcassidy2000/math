#!/usr/bin/env python3
"""Referee: plain A6 sweep with one closed-form jump derived by the referee:
P(8, w0) at depth j = 0 (mod 4), j > 0, end case e0ueven  ->  P(8, w0 + 3j + 20), same depth, eps and u parity kept;
the run continues while w0 + k(3j+20) + 31 + 3j < 56 * 2^j.  Verified against my_phase stepping below."""
import sys, os, random
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from a6_my_orbit import my_phase


def depth_v(v0, w0):
    A = v0 + 6
    j = 0
    while not (w0 + 2 * A + 3 * j + 3 < (4 * A) << j):
        j += 1
    return j


def is_e0ueven(v0, w0, j):
    A = v0 + 6
    v = (A << j) - 6
    w = w0 - 2 * A * ((1 << j) - 1) + 3 * j
    if w < 3:
        return False
    t, eps = divmod(w - 3, 2)
    u = v + 3 - t
    return eps == 0 and u % 2 == 0


def my_fast_sweep(W, cap=10 ** 8):
    v0, w0 = 0, W
    for n in range(cap):
        if v0 == 8:
            r, j = my_phase(8, w0)
            if j > 0 and j % 4 == 0 and r[0] == 'P' and r[1] == 8 and is_e0ueven(8, w0, j):
                delta = 3 * j + 20
                K = (56 * (1 << j) - 32 - 3 * j - w0) // delta
                if K >= 1:
                    w0 += (K + 1) * delta
                    continue
        r, j = my_phase(v0, w0)
        if r[0] == 'H':
            return 'HALT'
        if r[0] == 'R':
            return r[1]
        v0, w0 = r[1], r[2]
    raise RuntimeError(('cap', v0, w0))


if __name__ == '__main__':
    rng = random.Random(5)
    bad = 0
    for _ in range(3000):
        j = rng.choice([4, 8, 12])
        w0 = rng.randrange(0, 56 << j)
        a = my_fast_sweep_part = None
        # compare the run jump with explicit stepping from P(8, w0) until v0 != 8 or depth change
        v0, w = 8, w0
        r, jj = my_phase(8, w)
        if not (jj == j and jj % 4 == 0 and is_e0ueven(8, w, jj)):
            continue
        # explicit
        x = w
        steps = 0
        while True:
            r, jj = my_phase(8, x)
            if r[0] == 'P' and r[1] == 8 and jj == j and is_e0ueven(8, x, jj):
                x = r[2]
                steps += 1
            else:
                break
        delta = 3 * j + 20
        K = (56 * (1 << j) - 32 - 3 * j - w) // delta
        jumped = w + (K + 1) * delta
        if K >= 1 and jumped != x:
            bad += 1
    print('run jump vs explicit stepping:', 'OK' if bad == 0 else f'{bad} BAD')
