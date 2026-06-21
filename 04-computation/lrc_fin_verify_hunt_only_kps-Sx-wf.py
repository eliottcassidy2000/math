#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Standalone adversarial HUNT (V2/V2b) for the glue-and-global-sup verification.
Writes its own output file so there is no tee buffering loss. EXACT rationals."""
import itertools, random
from fractions import Fraction as F
from math import gcd

CAPS = {8: F(2243, 5880), 9: F(1979, 4004), 10: F(55, 91), 11: F(66, 91), 12: F(6, 7)}
INNER = set(range(1, 7))

def sector_of(p):
    return int((p % 1) * 7)

def p0_inner(E):
    E = sorted(set(int(e) for e in E))
    bps = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        for m in range(0, 7 * e + 1):
            bps.add(F(m, 7 * e))
    bps = sorted(bps)
    tot = F(0)
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0:
            continue
        xm = (x0 + x1) / 2
        if INNER <= {sector_of(e * xm) for e in E}:
            tot += x1 - x0
    return tot

def primitive(E):
    g = 0
    for e in E:
        g = gcd(g, int(e))
    return g == 1

OUT = open("05-knowledge/results/lrc_fin_verify_hunt_only_kps-Sx-wf.out", "w", encoding="utf-8")
def w(*a):
    print(*a); print(*a, file=OUT); OUT.flush()

champ = {k: None for k in range(8, 13)}
viols = []

def consider(E, tag):
    E = tuple(sorted(set(int(e) for e in E)))
    if 0 not in E:
        E = tuple(sorted(set((0,) + E)))
    k = len(E)
    if k not in CAPS or not primitive(E):
        return
    p = p0_inner(E)
    r = p / CAPS[k]
    if p > CAPS[k]:
        viols.append((tag, E, float(p), float(CAPS[k])))
    if champ[k] is None or r > champ[k][0]:
        champ[k] = (r, E, p, tag)

rng = random.Random(20260620)

w("== V2b SPLIT vs SINGLE BLOCK (does splitting raise the cover?) ==")
# single block consec_m proxy vs two widely separated consec clusters
for k in range(8, 11):
    sb = p0_inner(range(k))  # the pinned single block, the upper reference
    worst_split = F(0); arg = None
    for sep in (10, 14, 17, 20, 50, 100, 1000, 10000):
        for sz in range(1, k):
            E = tuple(list(range(sz)) + list(range(sep, sep + (k - sz))))
            if len(set(E)) != k:
                continue
            p = p0_inner(E)
            if p > worst_split:
                worst_split = p; arg = (sep, sz, E)
    w(f"  k={k}: pinned single-block p0={float(sb):.5f}  best 2-cluster split p0={float(worst_split):.5f}  "
      f"split<pinned? {worst_split < sb}  arg={arg}")

w("\n== V2 HUNT ==")
# (a) exhaustive span bands 15..20 k=8,9 ; 15..17 k=10
w("  (a) exhaustive span bands ...")
for s in range(15, 21):
    for k in (8, 9):
        interior = list(range(1, s)); need = k - 2
        if need > len(interior):
            continue
        for combo in itertools.combinations(interior, need):
            consider((0,) + combo + (s,), f"band-s{s}-k{k}")
for s in range(15, 18):
    interior = list(range(1, s)); need = 8
    if need <= len(interior):
        for combo in itertools.combinations(interior, need):
            consider((0,) + combo + (s,), f"band-s{s}-k10")
# (b) dilated APs d!=1 (scale-invariance stress)
w("  (b) dilated APs d!=1 ...")
for k in range(8, 13):
    for d in (1, 2, 3, 4, 5, 6, 7, 9, 11, 13, 17, 23):
        consider(tuple(d * i for i in range(k)), f"dilAP-d{d}-k{k}")
# (c) geometric / two-cluster resonance
w("  (c) geometric / two-cluster resonance ...")
for k in range(8, 13):
    for sep in (7, 10, 14, 20, 49, 50, 100):
        for sz in range(1, k):
            consider(tuple(list(range(sz)) + list(range(sep, sep + (k - sz)))), f"2clu-sep{sep}-sz{sz}-k{k}")
    for ratio in (2, 3):
        pts = []; v = 1
        while len(pts) < k:
            pts.append(v); pts.append(v + 1); v *= ratio
        consider(tuple([0] + pts[:k - 1]), f"geom-r{ratio}-k{k}")
# (d) multi-far collars random
w("  (d) multi-far collars (random) ...")
for k in range(8, 13):
    for _ in range(8000):
        nfar = rng.randint(2, 4); ncore = k - nfar
        if ncore < 2:
            continue
        core = sorted(set(rng.sample(range(0, 15), ncore)) | {0})
        far = sorted(rng.sample(range(15, 400), k - len(core)))
        consider(tuple(core) + tuple(far), f"collar-{nfar}far-k{k}")
# (e) large random wide rows
w("  (e) large random wide rows ...")
for k in range(8, 13):
    for _ in range(8000):
        E = sorted(rng.sample(range(1, 500), k - 1))
        consider((0,) + tuple(E), f"rand-k{k}")

w("\n  k :  worst p0/cap     p0        cap       via                  champion E")
wo = F(0)
for k in range(8, 13):
    r, E, p, tag = champ[k]
    wo = max(wo, r)
    fl = "  *** EXCEEDS CAP ***" if r > 1 else ""
    w(f"  {k:2d}:  {float(r):.5f}        {float(p):.5f}   {float(CAPS[k]):.5f}   {tag:18s}  {E}{fl}")
w(f"\n  GLOBAL worst p0/cap: {float(wo):.5f}  ({'NO VIOLATION' if wo <= 1 else 'VIOLATION'})")
w(f"  cap violations: {len(viols)}")
for v in viols[:10]:
    w("   ", v)
w("DONE.")
OUT.close()
