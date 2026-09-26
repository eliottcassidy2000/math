#!/usr/bin/env python3
"""
procgen_drift_20260926_local.py -- local search for small class-(i) flip sets of q n + 1 (upper bounds on
delta_k(q)), drift lane 2026-09-26 (session collatz-procgen-20260922).

Moves (all class-(i) tests are exact: the C engine's potential at the threshold F = best lower approximation of
log_q 2 with denominator <= 2^k, which is equivalent to "no expanding simple cycle"; every set that is reported is
re-verified by an edge-checked integer potential):
  * lift:   a level-k set R -> level k+1 set {r, r + 2^k : r in R} (same map T, same class);
  * prune:  try to restore sigma = + at each flipped residue (random order), incrementally;
  * kick:   remove m random flips, repair (flip odd nodes on harvested expanding cycles until class (i)), prune;
            accept the result if it is not larger.
"""
import os
import sys
import time
import random

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import procgen_drift_20260926_lib as L  # noqa: E402


def lift(S, k):
    """level-k set -> level-(k+1) set with the same map"""
    return set(S) | set(r + (1 << k) for r in S)


def prune(k, q, S, rng, F=None, rounds=2):
    F = F if F is not None else L.best_lower_approx(1 << k, q)
    fl = L.flip_array(k, S)
    E = L.Engine(k, q, fl, F)
    L.check(E.solve(), "prune: start set is not class (i)")
    for _ in range(rounds):
        order = [int(v) for v in np.nonzero(E.flip())[0]]
        rng.shuffle(order)
        changed = False
        for v in order:
            if E.toggle(v):
                changed = True
        if not changed:
            break
    return set(int(v) for v in np.nonzero(E.flip())[0])


def repair(k, q, S, rng, F=None, max_steps=200, M=40):
    """flip/unflip odd nodes on harvested expanding cycles until class (i); returns the set or None"""
    F = F if F is not None else L.best_lower_approx(1 << k, q)
    cur = set(S)
    for _ in range(max_steps):
        fl = L.flip_array(k, cur)
        cycs = L.disjoint_cycles(k, q, fl, M, F)
        if not cycs:
            return cur
        for cyc in cycs:
            cand = [v for v in cyc if (v & 1) and v not in cur]
            if cand and rng.random() < 0.97:
                cur.add(rng.choice(cand))
            else:
                odd = [v for v in cyc if v & 1]
                v = rng.choice(odd)
                if v in cur:
                    cur.discard(v)
                else:
                    cur.add(v)
    return None


def kick_search(k, q, S, rng, time_limit=60.0, mmax=3, F=None, log=None):
    F = F if F is not None else L.best_lower_approx(1 << k, q)
    best = set(S)
    t0 = time.time()
    tries = 0
    while time.time() - t0 < time_limit:
        tries += 1
        cur = set(best)
        m = rng.randint(1, mmax)
        for v in rng.sample(sorted(cur), min(m, len(cur))):
            cur.discard(v)
        rep = repair(k, q, cur, rng, F)
        if rep is None:
            continue
        rep = prune(k, q, rep, rng, F)
        if len(rep) <= len(best):
            if len(rep) < len(best) and log:
                log(f"      kick: {len(best)} -> {len(rep)} (try {tries}, {time.time() - t0:.0f}s)")
            best = rep
    return best


def verify_set(k, q, S):
    fl = L.flip_array(k, S)
    st = L.classify(k, q, fl)
    L.check(st[0] == 'I', f"set is not class (i) (q={q}, k={k})")
    return st[1], st[2]


if __name__ == '__main__':
    q = int(sys.argv[1])
    k0 = int(sys.argv[2])
    k1 = int(sys.argv[3])
    tl = float(sys.argv[4]) if len(sys.argv) > 4 else 30.0
    import json
    start = json.loads(sys.argv[5]) if len(sys.argv) > 5 else None
    rng = random.Random(12345)
    S = set(start) if start else None
    for k in range(k0, k1 + 1):
        if S is None:
            raise SystemExit("need a start set")
        if max(S) >= (1 << k):
            pass
        else:
            while max(S) < (1 << (k - 1)) and k > 1:
                S = lift(S, k - 1)
                break
        S = prune(k, q, S, rng)
        S = kick_search(k, q, S, rng, time_limit=tl, log=lambda *a: print(*a, flush=True))
        verify_set(k, q, S)
        print(k, len(S), round(len(S) / 2 ** (k - 1), 4), round(k * len(S) / 2 ** (k - 1), 3), flush=True)
        S = lift(S, k)
