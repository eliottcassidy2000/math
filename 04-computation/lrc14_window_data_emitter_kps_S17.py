#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
THE WINDOW DATA EMITTER (kps-2026-07-02-S17, HYP-3974): the full (0, 20] bounded band as
PLAIN-LIST Lean data rows for the pack-list dispatcher (LRCWindowData.lean).

Same enumeration + exact-witness search as opus's S46 generator, but emits
`⟨[s1,...,s13], num, den⟩ : WinRow` DATA (no `![...]` tuple literals — the S47 crash was
in per-row theorem elaboration; data lists + one gate evaluation avoid it), chunked for
maxRecDepth.  A failure row (max-min < 1/14) would be an LRC(14) counterexample: assert
loudly.
"""
import sys
from fractions import Fraction as Fr
from math import gcd, floor
from functools import reduce
from itertools import combinations
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

W0 = int(sys.argv[1]) if len(sys.argv) > 1 else 20
OUT = sys.argv[2] if len(sys.argv) > 2 else None

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
    for d in sorted(dens):
        for k in range(d + 1):
            t = Fr(k, d)
            if min(dist(v * t) for v in S) >= Fr(1, 14):
                return t
    return None

rows, fails = [], 0
for S in combinations(range(1, W0 + 1), 13):
    if reduce(gcd, S) != 1: continue
    if not all(any(v % q == 0 for v in S) for q in range(2, 15)): continue
    t = best_witness(S)
    if t is not None:
        rows.append((S, t))
    else:
        fails += 1
        print(f"*** COUNTEREXAMPLE CANDIDATE (report immediately): {S}")

print(f"window (0, {W0}]: covering primitive 13-tuples: {len(rows)} rows, {fails} failures")
assert fails == 0

# sanity: verify the integer kernel-gate check mirrors: den <= 14*min(r, den-r) for each speed
for S, t in rows:
    num, den = t.numerator, t.denominator
    for s in S:
        r = (s * num) % den
        assert den <= 14 * min(r, den - r), (S, t, s)
print("kernel-gate integer mirror: all rows pass")

if OUT:
    CH = 500
    chunks = [rows[i:i+CH] for i in range(0, len(rows), CH)]
    with open(OUT, 'w', encoding='utf-8') as f:
        for ci, ch in enumerate(chunks):
            f.write(f"set_option maxRecDepth 4096 in\ndef winData{ci+1} : List WinRow :=\n  [")
            f.write(",\n   ".join(
                f"⟨[{', '.join(str(x) for x in S)}], {t.numerator}, {t.denominator}⟩"
                for S, t in ch))
            f.write("]\n\n")
        f.write("/-- The full (0, 20] window data: all chunks. -/\n")
        f.write("def winData : List WinRow :=\n  " +
                " ++ ".join(f"winData{i+1}" for i in range(len(chunks))) + "\n")
    print(f"emitted {len(rows)} rows in {len(chunks)} chunks -> {OUT}")
