#!/usr/bin/env python3
"""A6 = 1RB0LF_1RC1RA_0RD0RB_1LE---_0LA1LE_1RD1LA: the sweep map J(0,0,W) -> J(0,0,W') and the halting set (S56).

S(W) := next visit to J(0,0,.) from J(0,0,W) under the paper's simplified rules (A.6), or HALT.  We compute:
  - the exact halting set H = {W : the sweep from J(0,0,W) halts before returning} on windows [X, 2X);
  - its density h(X) = |H cap [X,2X)| / X, to test h(X) ~ c/X;
  - the distribution of S(W) - W and S(W)/W on windows (growth law of the sweep map);
  - the actual orbit W_0 = 0, W_{n+1} = S(W_n): growth vs n, and the heuristic expected number of halts
    sum_n h(W_n) (with h interpolated from the measured densities).
Usage: python a6_sweep_map_20261025.py
"""
import sys, json, os, math, random
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__))
from deepmind_cryptids_20261025 import step, HALT


def sweep(W, cap=100000):
    p = step('A6', (0, 0, W))
    n = 1
    while p != HALT and not (p[0] == 0 and p[1] == 0) and n < cap:
        p = step('A6', p)
        n += 1
    if p == HALT:
        return HALT, n
    if n >= cap:
        return None, n
    return p[2], n


if __name__ == '__main__':
    out = {}
    # (1) halting-set density on dyadic windows (exhaustive for small X, sampled for large)
    dens = []
    for e in range(4, 25):
        X = 1 << e
        if e <= 16:
            Ws = range(X, 2 * X)
            nW = X
        else:
            rng = random.Random(e)
            Ws = [rng.randrange(X, 2 * X) for _ in range(60000)]
            nW = len(Ws)
        h = 0
        ratios = []
        for W in Ws:
            r, _ = sweep(W)
            if r == HALT:
                h += 1
            elif r is not None:
                ratios.append(r / W)
        frac = h / nW
        ratios.sort()
        med = ratios[len(ratios) // 2] if ratios else None
        meanlog = sum(math.log(x) for x in ratios) / len(ratios) if ratios else None
        dens.append({'X': X, 'halting_fraction': frac, 'X_times_fraction': frac * X, 'median_ratio': med,
                     'mean_log_ratio': meanlog, 'sampled': e > 16})
        print(f'X=2^{e}: halting fraction {frac:.3e} (x X = {frac * X:.3f}); S(W)/W median {med:.4f}, '
              f'E log ratio {meanlog:.4f}', flush=True)
    out['density'] = dens
    # (2) the actual orbit
    W = 0
    orbit = [W]
    t_sweeps = 0
    while len(orbit) < 400000:
        r, _ = sweep(W)
        if r == HALT:
            orbit.append('HALT')
            break
        W = r
        orbit.append(W)
    L = len(orbit)
    print(f'orbit: {L} sweeps, halted: {orbit[-1] == "HALT"}, last W = {orbit[-1]}')
    # growth law: W_n vs n
    samples = [(n, orbit[n]) for n in (10, 100, 1000, 10000, 50000, 100000, 200000, 300000, L - 1) if n < L]
    print('orbit samples (n, W_n):', samples)
    # heuristic expected halts: sum c/W_n with c = median of X*fraction over large windows
    cs = [d['X_times_fraction'] for d in dens if d['X'] >= 1 << 12 and d['halting_fraction'] > 0]
    c = sorted(cs)[len(cs) // 2] if cs else None
    if c:
        exp_halts = sum(c / max(1, w) for w in orbit[1:] if isinstance(w, int) and w > 0)
        tail = sum(c / max(1, w) for w in orbit[L // 2:] if isinstance(w, int) and w > 0)
        print(f'c (median X*density, X >= 4096) = {c:.3f}; heuristic expected halts over the orbit = {exp_halts:.3f}; '
              f'over its second half = {tail:.4f}')
        out['c'] = c
        out['expected_halts'] = exp_halts
        out['expected_halts_second_half'] = tail
    out['orbit_samples'] = samples
    json.dump(out, open(os.path.join(HERE, 'a6_sweep_map.json'), 'w'), indent=1, default=str)
