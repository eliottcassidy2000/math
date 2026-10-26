#!/usr/bin/env python3
"""Referee: own closed-form A6 phase map (derived by the referee; cross-checked against raw rules), plain sweep
orbit from W = 0, and Q-monotonicity check (2 v0 + w0 nondecreasing within a sweep, bounded by the next W)."""
import sys, os, time, random
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def my_phase(v0, w0):
    """Returns ('H',) | ('R', W') | ('P', v0', w0'), depth j."""
    A = v0 + 6
    # depth: least j with w0 + 2A + 3j + 3 < 4A 2^j
    j = 0
    x = w0 + 2 * A + 3
    # start from an estimate
    if x >= 4 * A:
        j = max(0, (x // (4 * A)).bit_length() - 2)
        while not (x + 3 * j < (4 * A) << j):
            j += 1
    v = (A << j) - 6
    w = w0 - 2 * A * ((1 << j) - 1) + 3 * j
    if w <= 2:
        if v % 2 == 0:
            return (('R', 2 * v + 12), ('P', 8, 2 * v + 5), ('R', 2 * v + 18))[w], j
        return (('P', 8, 2 * v + 3), ('R', 2 * v + 14), ('R', 2 * v + 22))[w], j
    t, eps = divmod(w - 3, 2)
    u = v + 3 - t
    if eps:
        if u == 1:
            return ('H',), j
        if u == 2:
            return ('R', 4 * t + 10), j
        return (('P', u + 4, 4 * t + 1) if u % 2 == 0 else ('P', u, 4 * t + 5)), j
    if u == 1:
        return ('R', 4 * t + 12), j
    return (('P', 8, 2 * u + 4 * t + 1) if u % 2 == 0 else ('R', 2 * u + 4 * t + 10)), j


def my_sweep(W, cap=10 ** 7, track=False):
    v0, w0 = 0, W
    qmax = 0
    for n in range(cap):
        r, j = my_phase(v0, w0)
        if track:
            qmax = max(qmax, 2 * v0 + w0 + 3 * j)
        if r[0] == 'H':
            return 'HALT', n + 1, qmax
        if r[0] == 'R':
            return r[1], n + 1, qmax
        v0, w0 = r[1], r[2]
    raise RuntimeError('cap')


if __name__ == '__main__':
    # (a) cross-check my_phase against the raw rules (via a6_phase_check helpers) on random inputs
    from a6_phase_check import raw_until, my_depth
    from tm_rules_check import HALT
    rng = random.Random(11)
    bad = 0
    for _ in range(5000):
        sz = rng.choice([10, 1000, 10 ** 9, 2 ** 70])
        v0, w0 = rng.randrange(sz), rng.randrange(sz * 8)
        r, j = my_phase(v0, w0)
        tgt = HALT if r[0] == 'H' else ((0, 0, r[1]) if r[0] == 'R' else (0, r[1], r[2]))
        ok, passed, n = raw_until((0, v0, w0), tgt)
        if not ok or my_depth(v0, w0)[0] != j:
            bad += 1
    print('my_phase vs raw rules, 5000 random:', 'OK' if bad == 0 else f'{bad} BAD', flush=True)
    # (b) plain orbit
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 3 * 10 ** 6
    checkpoints = {3 * 10 ** 6, 2 * 10 ** 7}
    W = 0
    t0 = time.time()
    qviol = 0
    for i in range(1, N + 1):
        W2, nph, qmax = my_sweep(W, track=(i % 1000 == 0))
        if W2 == 'HALT':
            print('HALT at sweep', i, 'W', W)
            break
        if i % 1000 == 0 and qmax > W2 + 3 * 64:
            qviol += 1
        W = W2
        if i in checkpoints:
            print(f'plain orbit (referee phase map): after {i} sweeps W = {W} ({time.time() - t0:.0f}s); Q-bound violations {qviol}', flush=True)
    print('done', N, W)
