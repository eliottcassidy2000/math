#!/usr/bin/env python3
"""
kind-pasteur-2026-07-06-S20b: THE DECISIVE EXPERIMENT for the l >= 7 residual.

The torus-split rung (HYP-4247) kills l <= 6 lifted at rho = 2/25 (density
2*rho*l < 1).  At l = 7: 2*rho*7 = 28/25 > 1 -- density does NOT forbid
covering.  THE QUESTION that decides the endgame strategy:

  Can 7 (or more) DISTINCT-frequency combs {theta : ||r_i theta|| < rho},
  rho = 2/25, actually COVER the circle?

  * YES => the l >= 7 residual needs the HEIGHT/finiteness argument
    (J-K finite census + equidistribution threshold); density lane ends at 6.
  * NO (positive uncovered floor) => a Newman-shaped covering-impossibility
    to formalize; the gap is empty at l = 7 too, structurally.
"""
import numpy as np

RHO = 2.0 / 25.0
rng = np.random.default_rng(2026070620)

def uncovered_fraction(freqs, shifts, grid=100000):
    th = (np.arange(grid) + 0.5) / grid
    covered = np.zeros(grid, dtype=bool)
    for r, s in zip(freqs, shifts):
        x = r * th + s
        f = x - np.round(x)
        covered |= (np.abs(f) < RHO)
    return 1.0 - covered.mean()

# ---------- (1) exhaustive small distinct integer frequency sets, zero shift ----------
print("=== (1) distinct integer freqs 1..12, zero shift, l = 7 ===", flush=True)
from itertools import combinations
best = (1.0, None)
for combo in combinations(range(1, 13), 7):
    uf = uncovered_fraction(list(combo), [0.0]*7, grid=40000)
    if uf < best[0]:
        best = (uf, combo)
print(f"  scanned {len(list(combinations(range(1,13),7)))} sets; min uncovered = {best[0]:.5f} at freqs {best[1]}", flush=True)

# ---------- (2) shift-optimized covering attempts ----------
print("=== (2) shift-optimized covering attempts, l = 7 ===", flush=True)
def optimize_shifts(freqs, iters=3000, restarts=6, grid=30000):
    bestuf = 1.0
    for _ in range(restarts):
        shifts = rng.random(len(freqs))
        cur = uncovered_fraction(freqs, shifts, grid)
        for _ in range(iters):
            k = rng.integers(len(freqs))
            old = shifts[k]
            shifts[k] = rng.random()
            new = uncovered_fraction(freqs, shifts, grid)
            if new <= cur:
                cur = new
            else:
                shifts[k] = old
        bestuf = min(bestuf, cur)
    return bestuf

for freqs in [[1,2,3,4,5,6,7], [1,2,4,8,16,32,64],
              list(best[1]), [1,2,3,5,7,11,13], [6,7,8,9,10,11,12]]:
    uf = optimize_shifts(freqs)
    print(f"  freqs {freqs}: min uncovered over shifts = {uf:.5f}", flush=True)

# ---------- (3) all-equal frequency sanity ----------
print("=== (3) all r_i = 1, evenly spaced, l = 7 (expect 0) ===", flush=True)
shifts7 = np.array([-(i/7) for i in range(7)])
print(f"  uncovered = {uncovered_fraction([1]*7, shifts7):.5f}", flush=True)

# ---------- (4) the REAL question: DISTINCT freqs, aggressive cover search ----------
print("=== (4) DISTINCT frequencies, aggressive cover search, l = 7,8,9 ===", flush=True)
def search_distinct_cover(l, fmax=25, trials=2000, grid=25000):
    bestuf = 1.0; bestcfg = None
    for _ in range(trials):
        freqs = rng.choice(np.arange(1, fmax+1), size=l, replace=False)
        shifts = rng.random(l)
        uf = uncovered_fraction(freqs, shifts, grid)
        if uf < bestuf:
            bestuf = uf; bestcfg = (freqs.copy(), shifts.copy())
    if bestcfg:
        freqs, shifts = bestcfg
        cur = bestuf
        for _ in range(5000):
            k = rng.integers(l)
            old = shifts[k]; shifts[k] = rng.random()
            new = uncovered_fraction(freqs, shifts, grid)
            if new <= cur: cur = new
            else: shifts[k] = old
        bestuf = min(bestuf, cur); bestcfg = (freqs, shifts)
    return bestuf, bestcfg

for l in [7, 8, 9]:
    uf, cfg = search_distinct_cover(l)
    tag = "COVER FOUND" if uf < 1e-4 else "positive uncovered floor"
    print(f"  l={l}: min uncovered = {uf:.6f}  [{tag}]  freqs={sorted(cfg[0].tolist()) if cfg else None}", flush=True)

print(flush=True)
print("VERDICT:", flush=True)
print("  l=7 distinct uncovered ~ 0 => residual is HEIGHT-based (J-K finite), density lane ends at 6.", flush=True)
print("  positive floor persists => Newman covering-impossibility lives (formalize phi>0 at 2/25).", flush=True)
