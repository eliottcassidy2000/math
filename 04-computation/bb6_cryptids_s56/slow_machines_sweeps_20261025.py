#!/usr/bin/env python3
"""Sweep-level statistics of the slow-growing DeepMind BB(6) Cryptids A1, A2, A5, A6, A10 (S56).

For each machine we log every halting opportunity (the configuration at which an exact small value decides halting,
see the channel derivations in the session note) together with a scale (the largest parameter), and estimate the
heuristic halting hazard  sum over opportunities of P(exact hit) ~ sum 1/scale.  We also log the scale at each sweep
start to measure the growth law (linear vs geometric in the number of sweeps).
Channels:
  A6 : loop start J(u, v0, w), u >= 1, w >= 2: halts iff w = 2u - 1.        margin = w - (2u - 1)
  A2 : D(0, b odd, w): halts iff w = 2.                                      margin = w - 2
  A5 : D1(0, z odd, w): halts iff w = 2.                                     margin = w - 2
  A10: loop E(a, b, c), b >= 1, c >= 3 ending with b' = 0: halts iff c = 2b + 3 (margin = c - 2b - 3); also E(0,0,1).
  A1 : loop A(a, b, 0, d), a >= 3: ends at a' = 0 iff 3 | a and a/3 <= d; halts iff d = a + b + 2
       (margin = d - a - b - 2, recorded when 3 | a); and A(2, b, 0, d): halts iff d = b + 4 (margin = d - b - 4).
Usage: python slow_machines_sweeps_20261025.py [seconds_per_machine]
"""
import sys, time, json, os, math
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__))
from deepmind_cryptids_20261025 import step, START, HALT


def channel(nm, p):
    if nm == 'A6':
        u, v, w = p
        if u >= 1 and w >= 2:
            return w - (2 * u - 1), max(u, v, w)
    elif nm == 'A2':
        a, b, w = p
        if a == 0 and b % 2 == 1:
            return w - 2, max(b, w)
    elif nm == 'A5':
        a, z, w = p
        if a == 0 and z % 2 == 1:
            return w - 2, max(z, w)
    elif nm == 'A10':
        a, b, c = p
        if b >= 1 and c >= 3 and b <= (c - 1) // 2:
            return c - 2 * b - 3, max(a, b, c)
        if (a, b, c) == (0, 0, 1):
            return 0, 1
    elif nm == 'A1':
        a, b, c, d = p
        if c == 0 and a >= 3 and a % 3 == 0 and a // 3 <= d:
            return d - a - b - 2, max(a, b, d)
        if c == 0 and a == 2 and d > 0:
            return d - b - 4, max(b, d)
    return None


def sweep_start(nm, p):
    if nm == 'A6':
        return p[0] == 0 and p[1] == 0
    if nm == 'A2':
        return p[0] == 0 and p[1] == 1
    if nm == 'A5':
        return p[0] == 0 and p[1] == 2
    if nm == 'A10':
        return p[0] == 1 and p[1] == 0
    if nm == 'A1':
        return p[2] == 0 and p[3] == 0
    return False


if __name__ == '__main__':
    secs = float(sys.argv[1]) if len(sys.argv) > 1 else 60
    out = {}
    for nm in ('A6', 'A2', 'A5', 'A10', 'A1'):
        p = START[nm]
        n = 0
        n_opp = 0
        hazard = 0.0
        small_margin = Counter()
        sweeps = []
        rel = []
        t0 = time.time()
        halted = False
        while time.time() - t0 < secs:
            for _ in range(5000):
                if p == HALT:
                    halted = True
                    break
                ch = channel(nm, p)
                if ch is not None:
                    mg, sc = ch
                    n_opp += 1
                    hazard += 1.0 / max(1, sc)
                    if abs(mg) <= 6:
                        small_margin[mg] += 1
                    if n_opp % 997 == 0:
                        rel.append(mg / max(1, sc))
                if sweep_start(nm, p):
                    sweeps.append((n, max(p)))
                p = step(nm, p)
                n += 1
            if halted:
                break
        # growth law: scale at sweep starts vs sweep index
        S = len(sweeps)
        pts = [(i, sweeps[i][1]) for i in range(S) if sweeps[i][1] > 0]
        fit = None
        if len(pts) > 100:
            # log-log slope over the second half
            half = pts[len(pts) // 2:]
            x0, y0 = math.log(half[0][0] + 1), math.log(half[0][1])
            x1, y1 = math.log(half[-1][0] + 1), math.log(half[-1][1])
            fit = (y1 - y0) / (x1 - x0) if x1 > x0 else None
        rel_sorted = sorted(rel)
        res = {'halted': halted, 'macro_steps': n, 'opportunities': n_opp, 'hazard_sum_1_over_scale': hazard,
               'small_margins': sorted(small_margin.items()), 'sweeps': S,
               'final_scale_bits': max(p).bit_length() if p != HALT else None,
               'loglog_slope_scale_vs_sweep': fit,
               'sweep_samples': sweeps[::max(1, S // 12)][:14],
               'margin_over_scale_quantiles': [rel_sorted[int(q * (len(rel_sorted) - 1))] for q in (0, .1, .25, .5, .75, .9, 1)] if rel_sorted else None}
        out[nm] = res
        print(nm, json.dumps(res, default=str)[:900], flush=True)
    json.dump(out, open(os.path.join(HERE, f'slow_machines_sweeps_{int(secs)}s.json'), 'w'), indent=1, default=str)
