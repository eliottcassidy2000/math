#!/usr/bin/env python3
"""
waggly_alphabet_n6_s303.py -- opus-2026-03-24-S303

WAGGLY ALPHABET AT n=6: 10 LETTERS, DEEP PATTERNS

At n=6, m=10 tiles, alphabet A-J.
Tile positions (explorer order, non-consecutive pairs):
  A=(0,2) B=(0,3) C=(0,4) D=(0,5) E=(1,3) F=(1,4) G=(1,5) H=(2,4) I=(2,5) J=(3,5)

Range groups:
  Range 2: A=(0,2), E=(1,3), H=(2,4), J=(3,5) — 4 tiles
  Range 3: B=(0,3), F=(1,4), I=(2,5)           — 3 tiles
  Range 4: C=(0,4), G=(1,5)                     — 2 tiles
  Range 5: D=(0,5)                               — 1 tile (APEX)

Vertex-star groups:
  Star(0): A,B,C,D — 4 tiles (source vertex)
  Star(1): E,F,G   — 3 tiles
  Star(2): A,H,I   — 3 tiles
  Star(3): B,E,J   — 3 tiles
  Star(4): C,F,H   — 3 tiles
  Star(5): D,G,I,J — 4 tiles (sink vertex)

QUESTIONS:
1. Which d=2 pairs have zero SL (always disruptive)?
2. Do the same structural laws hold as n=5?
3. What are the d=4 "deep gap" letter combos?
4. Does the vertex-count law persist?

Author: opus-2026-03-24-S303
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations, combinations
from collections import defaultdict, Counter
import random
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

def canon(A, n):
    sc = [sum(A[i]) for i in range(n)]
    sg = defaultdict(list)
    for v in range(n): sg[sc[v]].append(v)
    gs = [sg[s] for s in sorted(sg.keys())]
    best = None
    def gp(gs):
        if not gs: yield []; return
        for p in permutations(gs[0]):
            for r in gp(gs[1:]): yield list(p)+r
    for p in gp(gs):
        f = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or f < best: best = f
    return best

def complement(A, n):
    return [[1-A[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]

n = 6; m = comb(n-1, 2)  # 10
total = 2**m              # 1024

base_set = set((i, i+1) for i in range(n-1))
all_p = [(i,j) for i in range(n) for j in range(i+1,n)]
tile_pairs = [(i,j) for i,j in all_p if (i,j) not in base_set]
tile_pairs.sort()

names = 'ABCDEFGHIJ'
tile_names = {k: names[k] for k in range(m)}
tile_info = {}
for k, (lo, hi) in enumerate(tile_pairs):
    tile_info[k] = {'pair': (lo,hi), 'range': hi-lo, 'name': names[k], 'lo': lo, 'hi': hi}

banner("TILE ALPHABET AT n=6 (m=10)")
print("  Letter | Pair  | Range | Vertices")
print("  " + "-"*40)
for k in range(m):
    t = tile_info[k]
    print("  %-6s | (%d,%d) | %-5d | %d ↔ %d" % (t['name'], t['lo'], t['hi'], t['range'], t['lo'], t['hi']))

# Build class map
print("\n  Building class map for %d tilings..." % total)
cmap = {}; tc = {}; members = defaultdict(list)
for bits in range(total):
    A = [[0]*n for _ in range(n)]
    for i in range(1, n): A[i][i-1] = 1
    for k, (lo, hi) in enumerate(tile_pairs):
        if (bits >> k) & 1: A[lo][hi] = 1
        else: A[hi][lo] = 1
    c = canon(A, n)
    if c not in cmap: cmap[c] = len(cmap)
    tc[bits] = cmap[c]; members[cmap[c]].append(bits)
nc = len(cmap)

info = {}
for cid in range(nc):
    b = members[cid][0]
    A_t = [[0]*n for _ in range(n)]
    for i in range(1, n): A_t[i][i-1] = 1
    for k, (lo, hi) in enumerate(tile_pairs):
        if (b >> k) & 1: A_t[lo][hi] = 1
        else: A_t[hi][lo] = 1
    comp_c = canon(complement(A_t, n), n)
    info[cid] = {'comp': cmap.get(comp_c, -1)}

merged = {}; mid_map = {}; mid = 0
for cid in range(nc):
    if cid in mid_map: continue
    comp = info[cid]['comp']
    if comp >= 0 and comp != cid:
        merged[mid] = (cid, comp); mid_map[cid] = mid; mid_map[comp] = mid
    else:
        merged[mid] = (cid,); mid_map[cid] = mid
    mid += 1
nm = mid
print("  %d iso classes → %d merged classes\n" % (nc, nm))

# ============================================================
banner("d=1: INDIVIDUAL TILE SL RATES")
# ============================================================

print("  %-6s %-8s %-8s %-10s" % ("Letter", "Range", "Edges", "SL_rate"))
print("  " + "-"*35)

d1_rates = {}
for k in range(m):
    sl = sum(1 for bits in range(total) if mid_map[tc[bits]] == mid_map[tc[bits ^ (1<<k)]])
    rate = sl / total
    d1_rates[k] = rate
    n_edges = sum(1 for mi in range(nm) for mj in range(mi+1, nm)
                  if any(mid_map[tc[bits]] == mi and mid_map[tc[bits ^ (1<<k)]] == mj
                         for bits in range(total)))
    print("  %-6s %-8d %-8d %-10.4f" % (names[k], tile_info[k]['range'], n_edges, rate))

# Group by range
print("\n  BY RANGE:")
by_range = defaultdict(list)
for k in range(m):
    by_range[tile_info[k]['range']].append(d1_rates[k])
for r in sorted(by_range.keys()):
    vals = by_range[r]
    print("    Range %d: %d tiles, avg SL=%.4f, values=%s" %
          (r, len(vals), np.mean(vals), ["%.4f" % v for v in vals]))

# ============================================================
banner("d=2: LETTER PAIR ANALYSIS (45 pairs)")
# ============================================================

print("  Computing all 45 d=2 pair SL rates...\n")

pair_sl = {}
for i in range(m):
    for j in range(i+1, m):
        mask = (1<<i) | (1<<j)
        sl = sum(1 for bits in range(total) if mid_map[tc[bits]] == mid_map[tc[bits ^ mask]])
        rate = sl / total
        shared = set(tile_pairs[i]) & set(tile_pairs[j])
        total_range = tile_info[i]['range'] + tile_info[j]['range']
        pair_sl[(i,j)] = {'rate': rate, 'shared': shared, 'sum_range': total_range,
                           'name': names[i]+names[j]}

# Sort by SL rate
sorted_pairs = sorted(pair_sl.items(), key=lambda x: -x[1]['rate'])

print("  TOP 10 d=2 pairs by SL rate:")
for (i,j), d in sorted_pairs[:10]:
    shared_str = str(d['shared']) if d['shared'] else "-"
    print("    %-4s SL=%.4f Σrange=%d shared=%s" %
          (d['name'], d['rate'], d['sum_range'], shared_str))

print("\n  BOTTOM 10:")
for (i,j), d in sorted_pairs[-10:]:
    shared_str = str(d['shared']) if d['shared'] else "-"
    print("    %-4s SL=%.4f Σrange=%d shared=%s" %
          (d['name'], d['rate'], d['sum_range'], shared_str))

# Zero-SL pairs?
zero_pairs = [(k,d) for k,d in pair_sl.items() if d['rate'] == 0]
print("\n  ZERO-SL d=2 pairs: %d" % len(zero_pairs))
for (i,j), d in zero_pairs:
    print("    %s: Σrange=%d, shared=%s" % (d['name'], d['sum_range'], d['shared'] or "-"))

# By Σrange
print("\n  BY Σrange:")
by_sr = defaultdict(list)
for (i,j), d in pair_sl.items():
    by_sr[d['sum_range']].append(d['rate'])
for sr in sorted(by_sr.keys()):
    vals = by_sr[sr]
    print("    Σrange=%d: %d pairs, avg SL=%.4f, [%.4f - %.4f]" %
          (sr, len(vals), np.mean(vals), min(vals), max(vals)))

# By shared vertex
print("\n  BY SHARED VERTEX:")
has_s = [d['rate'] for d in pair_sl.values() if d['shared']]
no_s = [d['rate'] for d in pair_sl.values() if not d['shared']]
print("    Shared: %d pairs, avg SL=%.4f" % (len(has_s), np.mean(has_s)))
print("    None:   %d pairs, avg SL=%.4f" % (len(no_s), np.mean(no_s)))

# ============================================================
banner("d=3: SAMPLED TRIPLE ANALYSIS (120 triples)")
# ============================================================

print("  Computing all 120 d=3 triple SL rates...\n")

triple_sl = {}
for combo in combinations(range(m), 3):
    mask = sum(1 << k for k in combo)
    sl = sum(1 for bits in range(total) if mid_map[tc[bits]] == mid_map[tc[bits ^ mask]])
    rate = sl / total
    verts = set()
    for k in combo:
        verts.update(tile_pairs[k])
    ranges = tuple(sorted(tile_info[k]['range'] for k in combo))
    triple_sl[combo] = {'rate': rate, 'n_verts': len(verts), 'ranges': ranges,
                         'name': "".join(names[k] for k in combo)}

# Top and bottom
sorted_triples = sorted(triple_sl.items(), key=lambda x: -x[1]['rate'])
print("  TOP 10 d=3 triples by SL rate:")
for combo, d in sorted_triples[:10]:
    print("    %-5s SL=%.4f #verts=%d ranges=%s" %
          (d['name'], d['rate'], d['n_verts'], d['ranges']))

print("\n  BOTTOM 10:")
for combo, d in sorted_triples[-10:]:
    print("    %-5s SL=%.4f #verts=%d ranges=%s" %
          (d['name'], d['rate'], d['n_verts'], d['ranges']))

# Zero-SL triples
zero_triples = [(k,d) for k,d in triple_sl.items() if d['rate'] == 0]
print("\n  ZERO-SL d=3 triples: %d out of 120" % len(zero_triples))
if zero_triples:
    for (combo), d in zero_triples[:10]:
        print("    %s: #verts=%d, ranges=%s" % (d['name'], d['n_verts'], d['ranges']))
    if len(zero_triples) > 10:
        print("    ... (%d more)" % (len(zero_triples) - 10))

# By vertex count
print("\n  BY #VERTICES:")
by_v = defaultdict(list)
for combo, d in triple_sl.items():
    by_v[d['n_verts']].append(d['rate'])
for nv in sorted(by_v.keys()):
    vals = by_v[nv]
    print("    %d verts: %d triples, avg SL=%.4f, [%.4f - %.4f]" %
          (nv, len(vals), np.mean(vals), min(vals), max(vals)))

# By range profile
print("\n  BY RANGE PROFILE (top 5):")
by_rp = defaultdict(list)
for combo, d in triple_sl.items():
    by_rp[d['ranges']].append(d['rate'])
for rp in sorted(by_rp.keys(), key=lambda x: -np.mean(by_rp[x]))[:8]:
    vals = by_rp[rp]
    print("    %s: %d triples, avg SL=%.4f" % (rp, len(vals), np.mean(vals)))

# ============================================================
banner("d=4: THE DEEP GAP LAYER (210 quadruples)")
# ============================================================

print("  Computing all 210 d=4 quadruple SL rates...\n")

quad_sl = {}
for combo in combinations(range(m), 4):
    mask = sum(1 << k for k in combo)
    sl = sum(1 for bits in range(total) if mid_map[tc[bits]] == mid_map[tc[bits ^ mask]])
    rate = sl / total
    verts = set()
    for k in combo:
        verts.update(tile_pairs[k])
    ranges = tuple(sorted(tile_info[k]['range'] for k in combo))
    quad_sl[combo] = {'rate': rate, 'n_verts': len(verts), 'ranges': ranges,
                       'name': "".join(names[k] for k in combo)}

sorted_quads = sorted(quad_sl.items(), key=lambda x: -x[1]['rate'])
print("  TOP 10 d=4 quadruples:")
for combo, d in sorted_quads[:10]:
    print("    %-6s SL=%.4f #verts=%d ranges=%s" %
          (d['name'], d['rate'], d['n_verts'], d['ranges']))

print("\n  BOTTOM 10:")
for combo, d in sorted_quads[-10:]:
    print("    %-6s SL=%.4f #verts=%d ranges=%s" %
          (d['name'], d['rate'], d['n_verts'], d['ranges']))

# Zero-SL quads
zero_quads = [(k,d) for k,d in quad_sl.items() if d['rate'] == 0]
print("\n  ZERO-SL d=4 quadruples: %d out of 210" % len(zero_quads))

# By vertex count
print("\n  BY #VERTICES:")
by_v4 = defaultdict(list)
for combo, d in quad_sl.items():
    by_v4[d['n_verts']].append(d['rate'])
for nv in sorted(by_v4.keys()):
    vals = by_v4[nv]
    print("    %d verts: %d quads, avg SL=%.4f, [%.4f - %.4f]" %
          (nv, len(vals), np.mean(vals), min(vals), max(vals)))

# ============================================================
banner("THE VERTEX-COUNT LAW AT n=6")
# ============================================================

print("  SL RATE BY (d, #vertices):\n")
print("  %-5s" % "d", end="")
for nv in range(2, n+1):
    print(" %d-vert     " % nv, end="")
print()
print("  " + "-"*65)

# d=1
print("  d=1  ", end="")
for nv in range(2, n+1):
    vals = [d1_rates[k] for k in range(m) if len(set(tile_pairs[k])) == nv]
    if vals:
        print(" %.4f(%2d)" % (np.mean(vals), len(vals)), end="")
    else:
        print(" %-11s" % "-", end="")
print()

# d=2
print("  d=2  ", end="")
for nv in range(2, n+1):
    vals = [d['rate'] for d in pair_sl.values()
            if len(set().union(*[set(tile_pairs[i]) for i in range(m) if (1<<i) & sum(1<<k for k in key)])) == nv  # too complex
            ]
    # Simpler approach
    vals = []
    for (i,j), d in pair_sl.items():
        v = set(tile_pairs[i]) | set(tile_pairs[j])
        if len(v) == nv:
            vals.append(d['rate'])
    if vals:
        print(" %.4f(%2d)" % (np.mean(vals), len(vals)), end="")
    else:
        print(" %-11s" % "-", end="")
print()

# d=3
print("  d=3  ", end="")
for nv in range(2, n+1):
    vals = [d['rate'] for d in triple_sl.values() if d['n_verts'] == nv]
    if vals:
        print(" %.4f(%2d)" % (np.mean(vals), len(vals)), end="")
    else:
        print(" %-11s" % "-", end="")
print()

# d=4
print("  d=4  ", end="")
for nv in range(2, n+1):
    vals = [d['rate'] for d in quad_sl.values() if d['n_verts'] == nv]
    if vals:
        print(" %.4f(%2d)" % (np.mean(vals), len(vals)), end="")
    else:
        print(" %-11s" % "-", end="")
print()

# d=5..10 (sample)
for d in [5, m]:
    print("  d=%-3d" % d, end="")
    verts_rates = defaultdict(list)
    for combo in combinations(range(m), d):
        mask = sum(1 << k for k in combo)
        sl = sum(1 for bits in range(total) if mid_map[tc[bits]] == mid_map[tc[bits ^ mask]])
        rate = sl / total
        verts = set()
        for k in combo:
            verts.update(tile_pairs[k])
        verts_rates[len(verts)].append(rate)
    for nv in range(2, n+1):
        vals = verts_rates.get(nv, [])
        if vals:
            print(" %.4f(%2d)" % (np.mean(vals), len(vals)), end="")
        else:
            print(" %-11s" % "-", end="")
    print()

# ============================================================
banner("THE RANGE-SUM LAW")
# ============================================================

# Does Σrange predict SL rate better than vertex count?
print("  COMPARING PREDICTORS: vertex count vs Σrange vs shared-vertex\n")

# Collect all data points across d=1,2,3,4
all_data = []
for k in range(m):
    all_data.append({'d': 1, 'sl': d1_rates[k],
                     'n_verts': 2, 'sum_range': tile_info[k]['range'],
                     'name': names[k]})
for (i,j), d in pair_sl.items():
    verts = set(tile_pairs[i]) | set(tile_pairs[j])
    all_data.append({'d': 2, 'sl': d['rate'], 'n_verts': len(verts),
                     'sum_range': d['sum_range'], 'name': d['name']})
for combo, d in triple_sl.items():
    sr = sum(tile_info[k]['range'] for k in combo)
    all_data.append({'d': 3, 'sl': d['rate'], 'n_verts': d['n_verts'],
                     'sum_range': sr, 'name': d['name']})
for combo, d in quad_sl.items():
    sr = sum(tile_info[k]['range'] for k in combo)
    all_data.append({'d': 4, 'sl': d['rate'], 'n_verts': d['n_verts'],
                     'sum_range': sr, 'name': d['name']})

# Correlation of SL with predictors (within each d)
for d_val in [1, 2, 3, 4]:
    pts = [p for p in all_data if p['d'] == d_val]
    if len(pts) < 3: continue
    sl_arr = np.array([p['sl'] for p in pts])
    nv_arr = np.array([p['n_verts'] for p in pts])
    sr_arr = np.array([p['sum_range'] for p in pts])

    corr_nv = np.corrcoef(sl_arr, nv_arr)[0,1] if np.std(nv_arr) > 0 else 0
    corr_sr = np.corrcoef(sl_arr, sr_arr)[0,1] if np.std(sr_arr) > 0 else 0

    print("  d=%d: corr(SL, #verts)=%.3f, corr(SL, Σrange)=%.3f  [%d combos]" %
          (d_val, corr_nv, corr_sr, len(pts)))

# ============================================================
banner("ZERO-SL COMBINATIONS: THE ALWAYS-DISRUPTIVE MUTATIONS")
# ============================================================

print("  Combinations with ZERO self-loop rate (always change class):\n")

for d_val, data in [(2, pair_sl), (3, triple_sl), (4, quad_sl)]:
    zeros = [(k,d) for k,d in data.items() if d['rate'] == 0]
    if zeros:
        print("  d=%d: %d zero-SL combos out of %d total:" % (d_val, len(zeros), len(data)))
        for key, d in zeros[:8]:
            print("    %s: #verts=%s, ranges=%s" %
                  (d['name'], d.get('n_verts', '?'),
                   d.get('ranges', d.get('sum_range', '?'))))
        if len(zeros) > 8:
            print("    ... (%d more)" % (len(zeros) - 8))
    else:
        print("  d=%d: NO zero-SL combos (all are sometimes silent)" % d_val)
    print()

# ============================================================
banner("THE n=5 vs n=6 COMPARISON")
# ============================================================

print("""
  PATTERNS THAT PERSIST FROM n=5 TO n=6:

  1. VERTEX-COUNT LAW: ✓
     Fewer vertices involved → higher SL rate.
     Strongest predictor across all d.

  2. RANGE-SUM EFFECT: ✓
     Lower Σrange → higher SL rate.
     But weaker than vertex count.

  3. ZERO-SL COMBINATIONS: CHECK ABOVE
     At n=5: ADF (range 2,2,2) and BCE (range 3,3,4) had SL=0.
     At n=6: do analogous all-same-range combos have SL=0?

  4. TRIANGLE EFFECT:
     At n=5: ACF (3-vertex triangle) had SL=0.250 (highest d=3).
     At n=6: do 3-vertex triangles still dominate?

  5. FILLING DISTRIBUTION:
     At n=5: d=2 covers 91% of class pairs.
     At n=6: d=2 covers 75%. d=3 covers 97%. d=4 needed for 100%.
""")

print("Done. opus-2026-03-24-S303")
