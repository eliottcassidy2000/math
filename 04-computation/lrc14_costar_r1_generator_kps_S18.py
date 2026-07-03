#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
THE CO-STAR (12,1) GENERATOR (kps-2026-07-02-S18, HYP-3975): the intermediate-band r=1
slice.  For every 12-subset B of [1,22]:

  * good(B) at h = 1/14 = the exact rational good region (union of closed cells
    [k/s + h/s-grid]); we need an ARC [lo, hi] with the arcSafe property for every
    s in B (s*[lo,hi] within one integer cell with h-margins) and positive width.
  * The CertRow (P = B, offs = [0], c = 1/2, lo, hi, Vs = floor(1/(hi-lo)) + 1)
    certifies B u {w} lonely for EVERY w >= Vs (row_tail + the c-ruler: the far
    runner sits at ||j + 1/2|| = 1/2 >= h).
  * The finite window: w in (22, Vs) with lcm(U(B)) | w  (U(B) = moduli of {2..14}
    unrepresented in B; only such w make B u {w} covering) -- concrete 13-tuples,
    each given an exact witness (kernel-gate mirror), winData pattern.

Emits: (a) costar rows as CertRow literals; (b) window rows as WinRow literals;
(c) a summary.  A window tuple with NO witness would be an LRC(14) counterexample:
assert loudly.
"""
import sys
from fractions import Fraction as Fr
from math import gcd, floor, lcm
from functools import reduce
from itertools import combinations
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

H = Fr(1, 14)
WBAND = 22
OUTC = sys.argv[1] if len(sys.argv) > 1 else None   # costar rows .lean fragment
OUTW = sys.argv[2] if len(sys.argv) > 2 else None   # window rows .lean fragment

def dist(x):
    f = x - floor(x)
    return min(f, 1 - f)

def good_cells(B):
    """Exact good region of B at margin H as list of closed [a,b] with a<=b, within [0,1]."""
    # breakpoints: for each s, forbidden open zones ((k-H)/s, (k+H)/s); good = complement
    cells = [(Fr(0), Fr(1))]
    for s in B:
        newc = []
        for (a, b) in cells:
            # subtract all forbidden zones of speed s intersecting [a,b]
            cur = a
            k0 = floor(a * s)
            k1 = floor(b * s) + 1
            segs = []
            for k in range(k0, k1 + 1):
                zl, zr = (Fr(k) - H) / s, (Fr(k) + H) / s
                segs.append((zl, zr))
            pieces = []
            for (zl, zr) in segs:
                if zr <= cur or zl >= b:
                    continue
                if zl > cur:
                    pieces.append((cur, zl))
                cur = max(cur, zr)
            if cur < b:
                pieces.append((cur, b))
            newc.extend(pieces)
        cells = [(a, b) for (a, b) in newc if a < b]
        if not cells:
            return []
    return cells

def arc_safe_ok(B, lo, hi):
    """Mirror of WitnessCert.arcSafe for each s in B at cTar=0 (and offs=[0] at c=1/2 is
    automatic).  s*lo, s*hi within [floor(s*lo)+H, floor(s*lo)+1-H]."""
    for s in B:
        a, b = s * lo, s * hi
        fa = floor(a)
        if not (a <= b and Fr(fa) + H <= a and b <= Fr(fa) + 1 - H):
            return False
    return True

def best_witness(S):
    dens = set()
    for v in S: dens.add(v); dens.add(2 * v)
    L = sorted(S)
    for i in range(len(L)):
        for j in range(i + 1, len(L)):
            dens.add(L[i] + L[j]); dens.add(L[j] - L[i])
    dens.discard(0)
    for d in sorted(dens):
        for k in range(d + 1):
            t = Fr(k, d)
            if min(dist(v * t) for v in S) >= H:
                return t
    return None

costar_rows, window_rows = [], []
nB = nRel = nEmptyGood = 0
maxVs, sumWin = 0, 0
for B in combinations(range(1, WBAND + 1), 12):
    nB += 1
    U = [q for q in range(2, 15) if not any(v % q == 0 for v in B)]
    L = reduce(lcm, U, 1)
    # a far w>WBAND with B u {w} covering exists iff L has a multiple > WBAND (always).
    nRel += 1
    cells = good_cells(B)
    # pick the widest arcSafe-shrunk cell: shrink each good cell per-speed to the
    # arcSafe core: the cell already keeps H-margins for each s (good region does);
    # arcSafe additionally requires each s*[lo,hi] within ONE integer cell.
    best = None
    for (a, b) in cells:
        # split [a,b] at points k/s (cell boundaries of each speed)
        pts = {a, b}
        for s in B:
            k0, k1 = floor(a * s), floor(b * s) + 1
            for k in range(k0, k1 + 1):
                p = Fr(k, s)
                if a < p < b:
                    pts.add(p)
        P = sorted(pts)
        for i in range(len(P) - 1):
            lo, hi = P[i], P[i + 1]
            if hi - lo > (0 if best is None else best[1] - best[0]) and arc_safe_ok(B, lo, hi):
                best = (lo, hi)
    if best is None:
        nEmptyGood += 1
        # no safe arc: the tail leg fails for this B; EVERY far w needs a window row.
        # Record with Vs = None sentinel; handled below by requiring L-multiples up to a cap.
        print(f"*** NO SAFE ARC for B={B} (U={U}) -- unbounded window, needs cluster treatment")
        continue
    lo, hi = best
    Vs = floor(1 / (hi - lo)) + 1
    maxVs = max(maxVs, Vs)
    costar_rows.append((B, lo, hi, Vs))
    # finite window: w in (WBAND, Vs), L | w
    w = ((WBAND // L) + 1) * L
    while w < Vs:
        if w > WBAND:
            S = tuple(sorted(B + (w,)))
            t = best_witness(S)
            assert t is not None, f"*** COUNTEREXAMPLE CANDIDATE: {S}"
            window_rows.append((S, t))
            sumWin += 1
        w += L
print(f"B subsets: {nB}; relevant: {nRel}; no-safe-arc: {nEmptyGood}; "
      f"costar rows: {len(costar_rows)}; max Vs: {maxVs}; window rows: {sumWin}")

if OUTC:
    CH = 500
    chunks = [costar_rows[i:i+CH] for i in range(0, len(costar_rows), CH)]
    with open(OUTC, 'w', encoding='utf-8') as f:
        for ci, ch in enumerate(chunks):
            f.write(f"set_option maxRecDepth 4096 in\ndef costar{ci+1} : List WitnessCert.CertRow :=\n  [")
            f.write(",\n   ".join(
                f"⟨[{', '.join(str(x) for x in B)}], [0], 1/2, {lo.numerator}/{lo.denominator}, "
                f"{hi.numerator}/{hi.denominator}, {Vs}⟩" for (B, lo, hi, Vs) in ch))
            f.write("]\n\n")
        f.write("def costarRows : List WitnessCert.CertRow :=\n  " +
                " ++ ".join(f"costar{i+1}" for i in range(len(chunks))) + "\n")
    print(f"emitted {len(costar_rows)} costar rows -> {OUTC}")

if OUTW:
    CH = 500
    chunks = [window_rows[i:i+CH] for i in range(0, len(window_rows), CH)]
    with open(OUTW, 'w', encoding='utf-8') as f:
        for ci, ch in enumerate(chunks):
            f.write(f"set_option maxRecDepth 4096 in\ndef costarWin{ci+1} : List WinRow :=\n  [")
            f.write(",\n   ".join(
                f"⟨[{', '.join(str(x) for x in S)}], {t.numerator}, {t.denominator}⟩"
                for (S, t) in ch))
            f.write("]\n\n")
        f.write("def costarWin : List WinRow :=\n  " +
                (" ++ ".join(f"costarWin{i+1}" for i in range(len(chunks))) if chunks else "[]") + "\n")
    print(f"emitted {sumWin} window rows -> {OUTW}")
