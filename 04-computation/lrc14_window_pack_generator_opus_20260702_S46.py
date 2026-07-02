#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
THE BOUNDED-WINDOW PACK GENERATOR (hwindow slices): enumerate every PRIMITIVE COVERING
13-tuple with max <= W0, compute its exact max-min rational witness (complete breakpoint
grid), and emit Lean rows through kps's rational-witness gate (LRCRatWitness.lean:
lonely_of_ratWitness, kernel-pure soundness; per-row native_decide).

Slice 1 (W0 = 18): 966 rows, 0 failures -> LRCWindowPack1.lean.
Scale-up: rerun with larger W0; split output per band (W0, W1] to keep per-file
native_decide counts manageable; the dispatch's HNF quotient dedupes dilations upstream
(we emit primitive tuples only).

A FAILURE ROW (max-min < 1/14) would be an LRC(14) COUNTEREXAMPLE -- the generator asserts
loudly rather than skipping.  opus-2026-07-02-S46.
"""
import sys, io
from fractions import Fraction as Fr
from math import gcd, floor
from functools import reduce
from itertools import combinations
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

W0 = int(sys.argv[1]) if len(sys.argv) > 1 else 18
WLO = int(sys.argv[2]) if len(sys.argv) > 2 else 0   # band lower cut: require max > WLO
OUT = sys.argv[3] if len(sys.argv) > 3 else None

def dist(x):
    f = x - floor(x)
    return min(f, 1 - f)

def best_witness(S):
    dens = set()
    for v in S: dens.add(v); dens.add(2 * v)
    L = sorted(S)
    for i in range(len(L)):
        for j in range(i + 1, len(L)):
            dens.add(L[i] + L[j]); dens.add(L[j] - L[i])
    dens.discard(0)
    best, argt = Fr(-1), None
    for d in sorted(dens):
        for k in range(d + 1):
            t = Fr(k, d)
            m = min(dist(v * t) for v in S)
            if m > best: best, argt = m, t
    return best, argt

rows, fails = [], 0
for S in combinations(range(1, W0 + 1), 13):
    if max(S) <= WLO: continue
    if reduce(gcd, S) != 1: continue
    if not all(any(v % q == 0 for v in S) for q in range(2, 15)): continue
    mm, t = best_witness(S)
    if mm >= Fr(1, 14): rows.append((S, t))
    else:
        fails += 1
        print(f"*** COUNTEREXAMPLE CANDIDATE (report immediately): {S} max-min {mm}")
print(f"band ({WLO}, {W0}]: covering primitive 13-tuples: {len(rows)} rows, {fails} failures")

if OUT and fails == 0:
    body = "\n\n".join(
        f"theorem winRow_{i} : Lonely 14 (![{','.join(str(x) for x in S)}] : Fin 13 → ℤ) "
        f"((({t.numerator}/{t.denominator} : ℚ)) : ℝ) :=\n  lonely_of_ratWitness (by native_decide)"
        for i, (S, t) in enumerate(rows))
    hdr = ("/-\nCopyright (c) 2026 The TournamentH7 project contributors. All rights reserved.\n"
           "Released under Apache 2.0 license as described in the file LICENSE.\n"
           "Authors: opus (LRC multi-agent project, generated pack)\n-/\n"
           "import TournamentH7.LRCKernelGate\n\n"
           f"/-! # Bounded-window census band ({WLO}, {W0}] -- generated pack -/\n\n"
           "namespace LonelyRunner\nnamespace WindowPack\n\nopen LonelyRunner\n\n")
    io.open(OUT, "w", encoding="utf-8", newline="\n").write(hdr + body + "\n\nend WindowPack\nend LonelyRunner\n")
    print("emitted", OUT)
