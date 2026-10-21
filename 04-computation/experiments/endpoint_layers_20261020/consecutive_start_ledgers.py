#!/usr/bin/env python3
"""Consecutive-start-port pattern of the closed endpoint-k layer: I = {2, ..., k} (x and the next k-2 middle positions), so
J = r_P(I) = {n-1-k, ..., n-3} (the k-1 positions just before c = n-2; I and J overlap at the centre for odd n). Two ledgers, split by epsilon = [b -> y] (= [x -> c]):
  epsilon = 0: the inherited closed-no-chord ledger Q1 = (c, W, b, z, a), Q2 = (z, c, W, b, a), Q3 = (a, z, b, W, c)  (W = middle minus b, c);
  epsilon = 1: the port-block insertion  Q1 = (a, b, x, U, z, J, c), Q2 = (a, x, U, z, J, b, c), Q3 = (a, x, U, J, b, c, z)
               with U = (3, ..., n-k-2) and J = (n-1-k, ..., n-3) in path order.
Checks every closed standard-blue presentation of the pattern at the given (n, k): directedness, blackness of each move and the
literal target K_P(T), via the packet's evaluate (all guards re-tested)."""
import sys, time
from collections import Counter
from pathlib import Path
PK = Path('/Users/e/Documents/GitHub/math-wt-mersenne-561/05-knowledge/results/endpoint_three_repairs_20261020h')
sys.path.insert(0, str(PK)); sys.path.insert(0, '.')
import probe as P
from endpoint_four_bank import endpoint_k
def ledgers(n, k):
    a, b, x, c, z = 0, 1, 2, n - 2, n - 1
    W = tuple(range(2, n - 2)); U = tuple(range(3, n - k - 1)); J = tuple(range(n - k - 1, n - 2))
    e0 = ((c,) + W + (b, z, a), (z, c) + W + (b, a), (a, z, b) + W + (c,))
    e1 = ((a, b, x) + U + (z,) + J + (c,), (a, x) + U + (z,) + J + (b, c), (a, x) + U + J + (b, c, z))
    return e0, e1
for n, k in [(int(sys.argv[i]), int(sys.argv[i + 1])) for i in range(1, len(sys.argv), 2)] or [(9, 4), (10, 4), (11, 4), (11, 5)]:
    I = tuple(range(2, k + 1)); e0, e1 = ledgers(n, k); stats = Counter(); t0 = time.time(); bad = []
    for t, u, v in endpoint_k(n, k):
        if u != I: continue
        if P.arc(t, 1, n - 2): stats['open_middle'] += 1; continue   # pattern is the closed middle
        eps = P.arc(t, 1, n - 3)
        assert eps == P.arc(t, 2, n - 2), 'blue reflection eps'
        w = e1 if eps else e0
        r = P.evaluate(t, w)
        key = ('eps1' if eps else 'eps0', 'ok' if r['success'] else ('illegal' if 'invalid_step' in r else 'blue'))
        stats[key] += 1
        if not r['success'] and len(bad) < 3: bad.append((eps, r.get('invalid_step'), r.get('blue_steps'), r.get('missing_arcs')))
    print(f"n = {n}, k = {k}, I = {I}: {dict(stats)}  ({time.time()-t0:.0f}s)  failures: {bad}", flush=True)
