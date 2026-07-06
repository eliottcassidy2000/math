#!/usr/bin/env python3
"""
kind-pasteur-2026-07-06-S20c: THE DISTINCT-FREQUENCY COVERING THRESHOLD at rho = 2/25.

S20b established: distinct-freq combs leave a POSITIVE uncovered floor at l = 7,8,9.
Only REPEATED frequencies (all r=1) cover.  Since LRC runners are DISTINCT speeds,
the lifted combs have DISTINCT frequencies -- so the covering case is excluded.

THIS determines the endgame: find L* = min l for which SOME distinct-frequency
set covers the circle at rho = 2/25.  In the torus split the base has >= 1 runner
(the citation needs a nonempty cited side to give a clear t0), so lifted <= 11.
  * If L* >= 12: the (A) residual is EMPTY -- distinct lifted combs can NEVER
    cover, the gap (1/13, 2/25) has NO coupled 2-torus value at all.
  * If L* <= 11: the residual survives at l in [L*, 11]; those need the height
    argument.

We optimize HARD (multi-restart annealing + basin polish, fine grid) and also
test large/coprime frequencies (the independence-favorable regime).

We also fit the coprime free-fraction model (21/25)^l = P(independent uncovered)
as the theoretical backstop and locate where it crosses grid resolution.
"""
import numpy as np

RHO = 2.0 / 25.0
rng = np.random.default_rng(2026070621)

def uncovered_fraction(freqs, shifts, grid=60000):
    th = (np.arange(grid) + 0.5) / grid
    covered = np.zeros(grid, dtype=bool)
    for r, s in zip(freqs, shifts):
        x = r * th + s
        covered |= (np.abs(x - np.round(x)) < RHO)
    return 1.0 - covered.mean()

def hard_cover_search(l, fmax, trials=1500, anneal=8000, grid=40000):
    """Best (min uncovered) over random distinct-freq sets + shift annealing."""
    bestuf, bestcfg = 1.0, None
    for _ in range(trials):
        freqs = rng.choice(np.arange(1, fmax+1), size=l, replace=False)
        shifts = rng.random(l)
        uf = uncovered_fraction(freqs, shifts, grid)
        if uf < bestuf:
            bestuf, bestcfg = uf, (freqs.copy(), shifts.copy())
    # anneal shifts of the best config
    if bestcfg is not None:
        freqs, shifts = bestcfg
        cur = uncovered_fraction(freqs, shifts, grid)
        T = 0.05
        for it in range(anneal):
            k = rng.integers(l)
            old = shifts[k]; shifts[k] = rng.random()
            new = uncovered_fraction(freqs, shifts, grid)
            if new <= cur or rng.random() < np.exp((cur - new)/max(T,1e-6)):
                cur = new
            else:
                shifts[k] = old
            T *= 0.9997
        bestuf = min(bestuf, cur); bestcfg = (freqs, shifts)
    return bestuf, bestcfg

print("=== distinct-frequency covering threshold at rho = 2/25 ===", flush=True)
print(f"  independence model (21/25)^l:", flush=True)
for l in range(6, 16):
    print(f"    l={l:2d}: (21/25)^l = {(21/25)**l:.5f}", flush=True)
print(flush=True)

print("  HARD search (min uncovered found), small freqs (fmax=l+3):", flush=True)
Lstar_small = None
for l in range(6, 13):
    uf, cfg = hard_cover_search(l, fmax=l+3, trials=1200, anneal=6000)
    covered = uf < 3e-4
    if covered and Lstar_small is None:
        Lstar_small = l
    print(f"    l={l:2d}: min uncovered = {uf:.6f}  {'<-- COVER' if covered else ''}  freqs={sorted(cfg[0].tolist())}", flush=True)

print(flush=True)
print("  HARD search, LARGE distinct freqs (fmax=60, independence-favorable):", flush=True)
Lstar_large = None
for l in range(8, 15):
    uf, cfg = hard_cover_search(l, fmax=60, trials=1500, anneal=6000)
    covered = uf < 3e-4
    if covered and Lstar_large is None:
        Lstar_large = l
    print(f"    l={l:2d}: min uncovered = {uf:.6f}  {'<-- COVER' if covered else ''}", flush=True)

print(flush=True)
print(f"VERDICT: L*(small freqs) = {Lstar_small}, L*(large freqs) = {Lstar_large}", flush=True)
print("  If both >= 12: the (A) residual is EMPTY (lifted <= 11 can never cover distinct).", flush=True)
print("  The relevant bound for LRC: total 12 runners, base >= 1 => lifted <= 11.", flush=True)
