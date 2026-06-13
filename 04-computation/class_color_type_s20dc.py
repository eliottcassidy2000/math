#!/usr/bin/env python3
"""
class_color_type_s20dc.py — PURE BLUE / PURE BLACK / MIXED class classification
kind-pasteur-2026-03-23-S20dc

Each tiling is either BLUE (grid-symmetric) or BLACK (not grid-symmetric).
Each iso class contains some number of each type.
This gives a CLASS COLOR TYPE:
  PURE BLUE = all tilings grid-symmetric (class has reversal anti-aut for ALL reps)
  PURE BLACK = no tilings grid-symmetric
  MIXED = both types present

Also: the BLUE FRACTION of a class = gs_count / fiber.
This is a continuous measure of how "blue" a class is.

Investigate: how does class color type relate to SC/NS, H, degree, etc.?
"""

import sys
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
import numpy as np
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  CLASS COLOR TYPE: PURE BLUE / PURE BLACK / MIXED")
print("  kind-pasteur-2026-03-23-S20dc")
print("=" * 80)

for n in [4, 5, 6, 7]:
    print(f"\n{'#'*70}")
    print(f"  n = {n}")
    print(f"{'#'*70}")
    t0 = time.time()

    tiles = [(r,c) for r in range(1,n-1) for c in range(1,n-r)]
    m = len(tiles)
    tile_arcs = [(r+c, c-1) for r,c in tiles]

    explorer_tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            explorer_tiles.append((x, y))
    tk = {(x,y): i for i, (x,y) in enumerate(explorer_tiles)}
    tm = [tk.get((n-y+1, n-x+1), -1) for x,y in explorer_tiles]

    def is_gs(bits):
        for i in range(m):
            ti = tm[i]
            if ti != i and ti >= 0 and ((bits >> i) & 1) != ((bits >> ti) & 1):
                return False
        return True

    def build_t(tb):
        a = [[0]*n for _ in range(n)]
        for j in range(1, n): a[j][j-1] = 1
        for k in range(m):
            i, j = tile_arcs[k]
            if tb & (1 << k): a[i][j] = 1
            else: a[j][i] = 1
        return a

    def canon(a):
        sc = [sum(a[i][j] for j in range(n)) for i in range(n)]
        sg = defaultdict(list)
        for v in range(n): sg[sc[v]].append(v)
        gs = [sg[s] for s in sorted(set(sc))]
        if all(len(g)==1 for g in gs):
            p = [g[0] for g in gs]
            return tuple(a[p[i]][p[j]] for i in range(n) for j in range(n))
        best = None
        def gp(gs2):
            if not gs2: yield []; return
            for p in permutations(gs2[0]):
                for r in gp(gs2[1:]): yield list(p)+r
        for p in gp(gs):
            f = tuple(a[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or f < best: best = f
        return best

    def Hdp(a):
        dp = {}
        for v in range(n): dp[(1<<v, v)] = 1
        for S in range(1, 1<<n):
            for v in range(n):
                if not (S & (1<<v)): continue
                val = dp.get((S,v), 0)
                if val == 0: continue
                for u in range(n):
                    if S & (1<<u): continue
                    if a[v][u]:
                        dp[(S|(1<<u), u)] = dp.get((S|(1<<u), u), 0) + val
        return sum(dp.get(((1<<n)-1, v), 0) for v in range(n))

    # Build all tilings
    class_data = defaultdict(lambda: {'tilings': [], 'gs': 0, 'H': None, 'sc': None})
    for tb in range(1 << m):
        T = build_t(tb)
        cn = canon(T)
        gs = is_gs(tb)
        cd = class_data[cn]
        cd['tilings'].append(tb)
        if gs: cd['gs'] += 1
        if cd['H'] is None:
            cd['H'] = Hdp(T)
            Tc = [[1-T[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]
            cd['sc'] = (canon(Tc) == cn)
            cd['c3'] = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
                         if (T[i][j] and T[j][k] and T[k][i]) or (T[i][k] and T[k][j] and T[j][i]))

        if (tb+1) % 5000 == 0:
            print(f"    {tb+1}/{1<<m} ({time.time()-t0:.0f}s)")

    classes = sorted(class_data.values(), key=lambda x: x['H'])
    print(f"  {len(classes)} classes, built in {time.time()-t0:.0f}s")

    # Classify each class
    pure_blue = [c for c in classes if c['gs'] == len(c['tilings'])]
    pure_black = [c for c in classes if c['gs'] == 0]
    mixed = [c for c in classes if 0 < c['gs'] < len(c['tilings'])]

    print(f"\n  CLASS COLOR TYPE:")
    print(f"    PURE BLUE (all tilings grid-sym): {len(pure_blue)}")
    print(f"    PURE BLACK (no grid-sym tilings): {len(pure_black)}")
    print(f"    MIXED (some of each): {len(mixed)}")

    # Cross-tabulate with SC/NS
    print(f"\n  COLOR TYPE vs SC/NS:")
    ct = Counter()
    for c in classes:
        if c['gs'] == len(c['tilings']): color = 'PURE_BLUE'
        elif c['gs'] == 0: color = 'PURE_BLACK'
        else: color = 'MIXED'
        ct[(color, 'SC' if c['sc'] else 'NS')] += 1

    print(f"    {'':>12} {'SC':>6} {'NS':>6} {'total':>6}")
    for color in ['PURE_BLUE', 'MIXED', 'PURE_BLACK']:
        sc = ct.get((color, 'SC'), 0)
        ns = ct.get((color, 'NS'), 0)
        print(f"    {color:>12} {sc:6d} {ns:6d} {sc+ns:6d}")

    # Blue fraction distribution
    print(f"\n  BLUE FRACTION (gs/fiber) distribution:")
    fracs = [c['gs']/len(c['tilings']) for c in classes]
    frac_dist = Counter(round(f, 3) for f in fracs)
    for fval in sorted(frac_dist.keys()):
        print(f"    bf={fval:.3f}: {frac_dist[fval]} classes")

    # H distribution by color type
    print(f"\n  H RANGES by color type:")
    for label, subset in [('PURE BLUE', pure_blue), ('MIXED', mixed), ('PURE BLACK', pure_black)]:
        if subset:
            Hs = [c['H'] for c in subset]
            print(f"    {label:>12}: H in [{min(Hs)}, {max(Hs)}], avg={np.mean(Hs):.1f}")

    # CONJECTURE: PURE BLUE classes are exactly the classes where the
    # reversal permutation is an anti-automorphism for ALL labeled tournaments
    # in the class that contain P_0.
    #
    # A class is PURE BLUE iff EVERY tiling in its fiber is grid-symmetric.
    # This means: for every labeled tournament T in the class with P_0 as a Ham path,
    # the reversal sigma(v)=n+1-v is an anti-automorphism of T.
    #
    # QUESTION: Is PURE BLUE = SC? Or is PURE BLUE stricter?
    print(f"\n  IS PURE BLUE = SC?")
    pb_sc = sum(1 for c in pure_blue if c['sc'])
    pb_nsc = sum(1 for c in pure_blue if not c['sc'])
    print(f"    PURE BLUE that are SC: {pb_sc}")
    print(f"    PURE BLUE that are NOT SC: {pb_nsc}")
    if pb_nsc > 0:
        print(f"    ANSWER: NO — some pure blue classes are not SC!")
    elif pb_sc < sum(1 for c in classes if c['sc']):
        sc_pb = sum(1 for c in classes if c['sc'] and c['gs'] == len(c['tilings']))
        sc_total = sum(1 for c in classes if c['sc'])
        print(f"    SC that are PURE BLUE: {sc_pb}/{sc_total}")
        print(f"    ANSWER: PURE BLUE is STRICTER than SC")
    else:
        print(f"    ANSWER: PURE BLUE = SC exactly!")

print(f"\n{'='*70}")
print(f"  CROSS-n SUMMARY")
print(f"{'='*70}")

print("""
  At each n, every tiling is DEFINITIVELY blue (grid-sym) or black (not).
  This is a perfect binary partition of the tiling space.

  At the class level, the split into PURE BLUE / PURE BLACK / MIXED
  depends on whether the class contains any grid-symmetric tilings.

  PURE BLUE classes are the 'most symmetric' — they have the reversal
  as an anti-automorphism for EVERY base-path representative.

  The relationship between PURE BLUE and SC:
  PURE BLUE => SC (grid-sym implies reversal anti-aut implies SC)
  SC does NOT => PURE BLUE (SC tournaments may use non-reversal anti-auts)
""")

print(f"  DONE.")
print("=" * 80)
