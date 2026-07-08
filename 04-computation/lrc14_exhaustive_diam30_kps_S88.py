#!/usr/bin/env python3
"""kps-2026-07-08-S88 -- EXTEND the k=11 exhaustive covering floor from prim-diam 24 (klein-S184)
to prim-diam 30, closing the small-scale tail where the AP+interior extremals live (opus-S155/kps-S87).

Method = klein-S184's batched numpy grid D3 (N grid, margin >> grid error), collecting every shape
with grid-D3 below a generous threshold for EXACT (Farey) re-verification.  For each PRIMITIVE 11-set
with primitive diameter D in [25,30]: grid-D3; report per-D min + all low shapes.  If min D3 >= bar,
the tail floor is exhaustively confirmed through prim-diam 30 (the d=3 AP+interior sits at 27)."""
import numpy as np
from math import gcd
from functools import reduce
from itertools import combinations
import sys
sys.path.insert(0, "."); sys.path.insert(0, "04-computation")

TH = 1/7; M = 6/7; BAR = 83549/252252
N = 200                                   # grid; error ~1/N ~ 5e-3 << margin ~0.12 (scan only; low shapes exact-verified)
x = (np.arange(N) + 0.5) / N
THRESH = 0.55                             # collect grid-D3 < THRESH for exact re-verify

def batch_floor(shapes):
    ph = (shapes[:, :, None] * x[None, None, :]) % 1.0     # B x 11 x N
    ph.sort(axis=1)
    g = np.diff(ph, axis=1); wrap = ph[:, :1, :] + 1.0 - ph[:, -1:, :]
    gaps = np.concatenate([g, wrap], axis=1)
    W = np.maximum(gaps - TH, 0).sum(axis=1)               # B x N
    m1 = W.mean(1); m2 = (W*W).mean(1); m3 = (W**3).mean(1)
    den = m2 - m3/M
    D3 = np.where(den > 0, m1/M + (m1 - m2/M)**2/np.where(den > 0, den, 1), m1/M)
    return D3

def main():
    print(f"EXTEND exhaustive k=11 covering floor to prim-diam 30 (kps-S88); bar={BAR:.5f}; grid N={N}", flush=True)
    gmin = (9.9, None); low_shapes = []
    for D in range(25, 31):
        mn = (9.9, None); cnt = 0; buf = []
        def flush():
            nonlocal mn, cnt
            if not buf: return
            arr = np.array(buf, dtype=np.float64)
            d3 = batch_floor(arr); cnt += len(buf)
            i = int(d3.argmin())
            if d3[i] < mn[0]: mn = (float(d3[i]), buf[i])
            for k in np.where(d3 < THRESH)[0]:
                low_shapes.append((float(d3[k]), buf[int(k)]))
            buf.clear()
        for mids in combinations(range(1, D), 9):
            E = (0,) + mids + (D,)
            if reduce(gcd, E) != 1: continue
            # reflection canonicalization (E and D-E have identical D3): keep the lex-smaller
            rev = (0,) + tuple(D - m for m in mids[::-1]) + (D,)
            if rev < E: continue
            buf.append(E)
            if len(buf) >= 8000: flush()
        flush()
        if mn[0] < gmin[0]: gmin = mn
        ok = "OK" if mn[0] >= BAR else "**BELOW BAR**"
        print(f"  prim-diam {D}: {cnt:>9} primitive shapes; min grid-D3={mn[0]:.5f} (margin {mn[0]-BAR:+.5f}) [{ok}]; minimizer={mn[1]}", flush=True)
    print(f"\nGLOBAL grid-min D3 over prim-diam 25..30 = {gmin[0]:.5f} at {gmin[1]} (margin {gmin[0]-BAR:+.5f})", flush=True)

    # EXACT re-verify the low shapes
    print(f"\nEXACT verification of the {len(low_shapes)} shapes with grid-D3 < {THRESH}:", flush=True)
    from lrc14_d3_exact_verify_klein_S184 import D3 as D3e
    low_shapes.sort()
    exact_min = (9.9, None); seen = set()
    for gd, E in low_shapes[:400]:                     # exact-verify the lowest (dedup)
        if E in seen: continue
        seen.add(E)
        de = float(D3e(E))
        if de < exact_min[0]: exact_min = (de, E)
    print(f"  EXACT min over the low shapes = {exact_min[0]:.6f} at {exact_min[1]} (margin {exact_min[0]-BAR:+.6f})", flush=True)
    print(f"  # distinct low shapes exact-checked: {len(seen)}", flush=True)
    verdict = "ALL prim-diam 25..30 CLEAR (exact min >= bar)" if exact_min[0] >= BAR else "SOME BELOW BAR -- investigate"
    print(f"\nVERDICT: {verdict}", flush=True)
    # show the 8 lowest exact
    print("  lowest exact D3 shapes:")
    exact_list = sorted({E: float(D3e(E)) for gd, E in low_shapes[:200]}.items(), key=lambda kv: kv[1])[:8]
    for E, de in exact_list:
        print(f"    D3={de:.6f}  prim-diam {max(E)-min(E)}  {E}", flush=True)

if __name__ == "__main__":
    main()
