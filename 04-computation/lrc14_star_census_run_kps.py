#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
HYP-3960: THE (star)-CENSUS RUN — DAG critical-path item 3 (kps-assigned).
kind-pasteur-2026-07-02-S3.

Universe (the admissible single-level slice): shapes (P, U) with P ⊆ {1..13}, |P| = 13 − k,
U ∋ 0, U ⊆ {0..14}, |U| = k, for k = 13, 12, 11, 10.  Per shape: search for a WITNESS-ARC
CERTIFICATE (c, lo, hi) — an arc on which every p ∈ P is safe at target 0 and every o ∈ U is safe
at target c (band h = 1/14) — via a c-grid sweep with the float interval engine (P-arcs
precomputed per P), refined 4x on failure.  PARTITION:
   [CERTIFIED single-level (width > 0; V* = ceil(1/width))] | [NEEDS NESTING / investigate].
The k=13 page is then EXACTLY verified (Fractions) and emitted as Lean CertRow literals for the
mass-table module (LRCCertTable.lean).
"""
import sys, itertools, time, math
from fractions import Fraction as Fr
from math import gcd, floor
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

H = 1.0/14.0

def clipf(arcs, v, c):
    out = []
    for (lo, hi) in arcs:
        k0 = int(floor(lo*v - c)) - 1; k1 = int(floor(hi*v - c)) + 1
        for k in range(k0, k1+1):
            a = (k + c + H)/v; b = (k + 1 + c - H)/v
            l = lo if lo > a else a
            r = hi if hi < b else b
            if l < r - 1e-14: out.append((l, r))
    return out

def P_arcs(P):
    arcs = [(0.0, 1.0)]
    for p in sorted(P): arcs = clipf(arcs, p, 0.0)
    return arcs

def widest_cert(Parcs, U, NC):
    best = None
    for ic in range(NC):
        c = (ic + 0.5)/NC
        if not (H <= c <= 1 - H): continue      # the o = 0 constraint
        arcs = Parcs
        ok = True
        for o in sorted(U):
            if o == 0: continue
            arcs = clipf(arcs, o, c)
            if not arcs: ok = False; break
        if not ok: continue
        for (lo, hi) in arcs:
            if best is None or hi - lo > best[0]: best = (hi - lo, c, lo, hi)
    return best

def census_k(k, sample_limit=None, NC=96):
    smalls = list(itertools.combinations(range(1, 14), 13 - k))
    offsets = [ (0,) + u for u in itertools.combinations(range(1, 15), k - 1) ]
    total = len(smalls) * len(offsets)
    import random
    rng = random.Random(3960 + k)
    if sample_limit and total > sample_limit:
        pairs = [ (rng.choice(smalls), rng.choice(offsets)) for _ in range(sample_limit) ]
        mode = f"SAMPLED {sample_limit}/{total}"
    else:
        pairs = [ (P, U) for P in smalls for U in offsets ]
        mode = f"FULL {total}"
    t0 = time.time()
    fails = []
    minw = None; maxVs = 0
    Pcache = {}
    n = 0
    for P, U in pairs:
        n += 1
        if P not in Pcache: Pcache[P] = P_arcs(P)
        b = widest_cert(Pcache[P], U, NC)
        if b is None:
            b = widest_cert(Pcache[P], U, NC*4)
        if b is None:
            fails.append((P, U)); continue
        w = b[0]
        if minw is None or w < minw[0]: minw = (w, P, U, b)
        Vs = int(1.0/w) + 1
        if Vs > maxVs: maxVs = Vs
    dt = time.time() - t0
    print(f"  k={k:2d} [{mode}]: certified {n - len(fails)}/{n}   min width={minw[0]:.6f} "
          f"(V* up to {maxVs}) at P={minw[1]} U={minw[2]}   fails={len(fails)}   ({dt:.0f}s)", flush=True)
    for P, U in fails[:6]:
        print(f"      NEEDS-NESTING/INVESTIGATE: P={P} U={U}")
    return fails, minw

print("="*104)
print(" THE (star)-CENSUS: admissible single-level shapes, window-14 offsets, k = 13..10")
print("="*104)
allfails = {}
for k in (13, 12, 11):
    fails, minw = census_k(k)
    allfails[k] = fails
fails10, minw10 = census_k(10, sample_limit=120000)
allfails[10] = fails10

# ---------------- exact verification + Lean emission for the k=13 page ----------------
print("\n" + "="*104)
print(" EXACT VERIFICATION (Fractions) + Lean rows: the k=13 page (P = empty, 91 offset patterns)")
print("="*104)
hq = Fr(1, 14)
def clipq(arcs, v, c):
    out = []
    for (lo, hi) in arcs:
        k0 = floor(lo*v - c) - 1; k1 = floor(hi*v - c) + 1
        for k in range(k0, k1+1):
            a = (k + c + hq)/v; b = (k + 1 + c - hq)/v
            l = max(lo, a); r = min(hi, b)
            if l < r: out.append((l, r))
    return out
def exact_cert(U, cf):
    # snap c to a nearby simple rational and recompute exactly
    c = Fr(round(cf*800), 800)
    if not (hq <= c <= 1 - hq): return None
    arcs = [(Fr(0), Fr(1))]
    for o in sorted(U):
        if o == 0: continue
        arcs = clipq(arcs, o, c)
        if not arcs: return None
    best = max(arcs, key=lambda ab: ab[1]-ab[0])
    return c, best[0], best[1]
rows = []
for U in [ (0,) + u for u in itertools.combinations(range(1, 15), 12) ]:
    b = widest_cert([(0.0,1.0)], U, 192)
    ex = exact_cert(U, b[1]) if b else None
    if ex is None:
        print(f"  U={U}: exact snap FAILED (needs finer c) -- flagged")
        continue
    c, lo, hi = ex
    w = hi - lo
    Vs = int(1/w) + 1
    rows.append((U, c, lo, hi, Vs))
print(f"  exact-certified k=13 rows: {len(rows)}/91   max V* = {max(r[4] for r in rows)}")
with open("../05-knowledge/results/lrc14_star_census_k13_leanrows_kps.txt", "w", encoding="utf-8") as f:
    for U, c, lo, hi, Vs in rows:
        f.write(f"  ⟨[{', '.join(str(u) for u in U)}], {c.numerator}/{c.denominator}, "
                f"{lo.numerator}/{lo.denominator}, {hi.numerator}/{hi.denominator}, {Vs}⟩,\n")
print("  Lean row literals written to lrc14_star_census_k13_leanrows_kps.txt")
print("DONE.")
