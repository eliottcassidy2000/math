#!/usr/bin/env python3
"""
lrc14_Bk_floor_decisive_kps-S5-wf.py  (kps-2026-06-18-S5-wf)

THE DECISIVE QUESTION for the bounded-spread reduction:
  For fixed k (esp. k=13), does  min_{|E|=k, spread=s} mu(E)  FLOOR at a positive c0 as s->inf,
  or keep DROPPING toward 0?

  If it floors  -> the reduction is sound (inf mu(E) = a positive minimum at bounded spread).
  If it -> 0    -> the reduction FAILS (no uniform floor; LRC(14) route via pure mu is broken,
                   though the G_P intersection rho* may still floor).

Method: fast vectorized numeric mu (numpy), simulated-annealing search over E at each spread,
many restarts. We push spread up to ~8k and watch the min-mu trajectory.
"""
import numpy as np, random, sys
from math import gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

X = (np.arange(200000) + 0.5) / 200000   # sample grid on [0,1)

def mu_num(E):
    # E: tuple of ints. points = outer(E, X) mod 1, per-column maxgap.
    P = np.mod(np.outer(np.array(E, dtype=np.float64), X), 1.0)  # shape (k, Nx)
    P.sort(axis=0)
    gaps = np.empty_like(P)
    gaps[:-1] = P[1:] - P[:-1]
    gaps[-1] = P[0] + 1.0 - P[-1]
    mg = gaps.max(axis=0)
    return float(np.mean(mg > 2.0/7.0))

def anneal_min(k, s, restarts=40, iters=400, seed=0):
    rng = random.Random(seed)
    best = 1.0; bestE = None
    if s - 1 < k - 2:
        return None, None
    for r in range(restarts):
        cur = sorted(rng.sample(range(1, s), k-2))
        E = tuple([0] + cur + [s])
        cm = mu_num(E)
        if cm < best and reduce(gcd, E) == 1:
            best = cm; bestE = E
        T = 0.05
        for it in range(iters):
            idx = rng.randint(0, k-3)
            nv = rng.randint(1, s-1)
            if nv in cur: continue
            cand = sorted(cur[:idx] + [nv] + cur[idx+1:])
            Ec = tuple([0] + cand + [s])
            mc = mu_num(Ec)
            if mc < cm or rng.random() < np.exp(-(mc-cm)/max(T,1e-4)):
                cur = cand; cm = mc; E = Ec
                if cm < best and reduce(gcd, E) == 1:
                    best = cm; bestE = E
            T *= 0.995
    return best, bestE

if __name__ == "__main__":
    k = int(sys.argv[1]) if len(sys.argv) > 1 else 13
    print(f"=== DECISIVE FLOOR TEST k={k} (vectorized numeric mu, annealing) ===")
    spreads = [k-1, k+3, 2*k, 3*k, 4*k, 6*k, 8*k]
    for s in spreads:
        b, E = anneal_min(k, s, restarts=30, iters=300, seed=s*13+k)
        print(f"  spread={s:3d} (~{s/k:.1f}k): min mu_num ~ {b:.4f}   E={E}")
    print("DONE.")
