#!/usr/bin/env python3
"""
kind-pasteur-2026-07-06-S20f: DOES CONSECUTIVE tile at 2/25?  (validating my
S20b/c 'distinct-freq dead through l=14' against mac-mini S5's self-correction
'>=10 consecutive DO tile at 2/25').

My S20b/c used RANDOM frequency samples + annealing; they may have MISSED the
consecutive-tiling configuration.  mac-mini found 7-9 empty (no cover) but
>=10 CONSECUTIVE frequencies DO cover.  If true, my 'dead through l=14' is an
OVERCLAIM (holds only for the sampled sets, and only l<=9 unconditionally).

Here we test CONSECUTIVE frequency blocks {1,...,l} directly, with heavy shift
optimization AND the natural 'staggered' shift guess (shift_i = phase making the
combs interleave).  We report min uncovered per l, and locate the tiling
threshold L_tile for consecutive.
"""
import numpy as np

RHO = 2.0 / 25.0
rng = np.random.default_rng(2026070622)

def uncovered_fraction(freqs, shifts, grid=200000):
    th = (np.arange(grid) + 0.5) / grid
    covered = np.zeros(grid, dtype=bool)
    for r, s in zip(freqs, shifts):
        x = r * th + s
        covered |= (np.abs(x - np.round(x)) < RHO)
    return 1.0 - covered.mean()

def anneal(freqs, restarts=12, iters=9000, grid=60000):
    best = 1.0; bestsh = None
    for rs in range(restarts):
        # restart 0 = zero shifts; 1 = staggered; else random
        if rs == 0:
            shifts = np.zeros(len(freqs))
        elif rs == 1:
            shifts = np.array([ (i * 0.5) % 1.0 for i in range(len(freqs)) ])
        else:
            shifts = rng.random(len(freqs))
        cur = uncovered_fraction(freqs, shifts, grid)
        T = 0.04
        for it in range(iters):
            k = rng.integers(len(freqs))
            old = shifts[k]; shifts[k] = rng.random()
            new = uncovered_fraction(freqs, shifts, grid)
            if new <= cur or rng.random() < np.exp((cur-new)/max(T,1e-6)):
                cur = new
            else:
                shifts[k] = old
            T *= 0.9997
        if cur < best:
            best, bestsh = cur, shifts.copy()
    return best, bestsh

print("=== CONSECUTIVE {1,...,l} at rho = 2/25: min uncovered (heavy anneal) ===", flush=True)
L_tile = None
for l in range(7, 14):
    freqs = list(range(1, l+1))
    uf, sh = anneal(freqs)
    tiles = uf < 1e-4
    if tiles and L_tile is None:
        L_tile = l
    # verify the best at a FINER grid to rule out grid-resolution false-zero
    uf_fine = uncovered_fraction(freqs, sh, grid=1000000) if sh is not None else uf
    print(f"  l={l:2d} freqs 1..{l}: min uncovered = {uf:.6f} (fine grid {uf_fine:.6f}) "
          f"{'<-- TILES' if tiles else ''}", flush=True)

print(flush=True)
print(f"CONSECUTIVE tiling threshold L_tile = {L_tile}", flush=True)
print("  If L_tile <= 11: my S20b/c 'dead through l=14' was an OVERCLAIM (random search", flush=True)
print("  missed the consecutive tiling); the covering-impossibility holds only l <= L_tile-1,", flush=True)
print("  and l >= L_tile consecutive needs opus phase-orbit / the finite census (NOT covering).", flush=True)
print("  mac-mini S5 self-correction (>=10 consecutive tile) would be CONFIRMED.", flush=True)
