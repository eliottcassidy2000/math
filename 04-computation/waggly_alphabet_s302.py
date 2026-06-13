#!/usr/bin/env python3
"""
waggly_alphabet_s302.py -- opus-2026-03-24-S302

WAGGLY LINES AS LETTER COMBINATIONS

At n=5, the wiggly alphabet has 6 letters: A,B,C,D,E,F.
Each wiggly (d=1) line uses 1 letter.
Each d=2 waggly uses a PAIR of letters: C(6,2)=15 pairs.
Each d=3 waggly uses a TRIPLE: C(6,3)=20 triples.
Etc.

QUESTION: Which letter combinations connect which class pairs?
Are there patterns?

TILE POSITIONS (n=5):
  A = (0,2): range 2   (row 0, col 2)
  B = (0,3): range 3   (row 0, col 3)
  C = (0,4): range 4   (row 0, col 4)  = APEX
  D = (1,3): range 2   (row 1, col 3)
  E = (1,4): range 3   (row 1, col 4)
  F = (2,4): range 2   (row 2, col 4)

Tiles grouped by RANGE: {A,D,F}=range2, {B,E}=range3, {C}=range4.
Tiles grouped by ROW: {A,B,C}=row0, {D,E}=row1, {F}=row2.
Tiles grouped by COL: {F}=col4-only, {C,E}=col4, {B,D}=col3, {A}=col2.

Author: opus-2026-03-24-S302
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations, combinations
from collections import defaultdict, Counter
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

n = 5
m = comb(n-1, 2)  # 6
total = 2**m       # 64

base_set = set((i, i+1) for i in range(n-1))
all_p = [(i,j) for i in range(n) for j in range(i+1,n)]
tile_pairs = [(i,j) for i,j in all_p if (i,j) not in base_set]
tile_pairs.sort()

tile_names = {0:'A', 1:'B', 2:'C', 3:'D', 4:'E', 5:'F'}
tile_info = {}
for k, (lo, hi) in enumerate(tile_pairs):
    rng = hi - lo
    tile_info[k] = {'pair': (lo,hi), 'range': rng, 'name': tile_names[k],
                     'lo': lo, 'hi': hi}

banner("TILE ALPHABET AT n=5")
print("  Letter | Pair  | Range | Vertices involved")
print("  " + "-"*50)
for k in range(m):
    t = tile_info[k]
    print("  %-6s | (%d,%d) | %-5d | %d ↔ %d" %
          (t['name'], t['lo'], t['hi'], t['range'], t['lo'], t['hi']))

# Build class map
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
    A = [[0]*n for _ in range(n)]
    for i in range(1, n): A[i][i-1] = 1
    for k, (lo, hi) in enumerate(tile_pairs):
        if (b >> k) & 1: A[lo][hi] = 1
        else: A[hi][lo] = 1
    comp_c = canon(complement(A, n), n)
    from math import comb as C
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

# ============================================================
banner("d=2: ALL 15 LETTER PAIRS — WHICH CLASSES DO THEY CONNECT?")
# ============================================================

# For each pair of letters (i,j), compute the weight matrix
pair_data = {}
for i in range(m):
    for j in range(i+1, m):
        combo_name = tile_names[i] + tile_names[j]
        mask = (1 << i) | (1 << j)

        W = np.zeros((nm, nm), dtype=int)
        for bits in range(total):
            mi = mid_map[tc[bits]]
            mj = mid_map[tc[bits ^ mask]]
            W[mi][mj] += 1

        n_cross = sum(1 for mi2 in range(nm) for mj2 in range(mi2+1, nm) if W[mi2][mj2] > 0)
        n_self = sum(1 for mi2 in range(nm) if W[mi2][mi2] > 0)
        sl_rate = sum(W[mi2][mi2] for mi2 in range(nm)) / total

        # Shared vertex between the two tiles?
        verts_i = set(tile_pairs[i])
        verts_j = set(tile_pairs[j])
        shared = verts_i & verts_j
        total_range = tile_info[i]['range'] + tile_info[j]['range']

        pair_data[(i,j)] = {
            'name': combo_name, 'edges': n_cross, 'self': n_self,
            'sl_rate': sl_rate, 'shared_vertex': shared,
            'total_range': total_range, 'W': W
        }

print("  %-6s %-8s %-6s %-8s %-10s %-12s" %
      ("Pair", "Shared?", "Edges", "SL_cls", "SL_rate", "Σrange"))
print("  " + "-"*55)

for (i,j), d in sorted(pair_data.items(), key=lambda x: -x[1]['sl_rate']):
    shared_str = str(d['shared_vertex']) if d['shared_vertex'] else "none"
    print("  %-6s %-8s %-6d %-8d %-10.3f %-12d" %
          (d['name'], shared_str, d['edges'], d['self'], d['sl_rate'], d['total_range']))

# ============================================================
banner("d=2 PAIRS GROUPED BY STRUCTURE")
# ============================================================

# Group by: shared vertex (yes/no), total range, triangle membership
print("  BY SHARED VERTEX:\n")

has_shared = [(k, d) for k, d in pair_data.items() if d['shared_vertex']]
no_shared = [(k, d) for k, d in pair_data.items() if not d['shared_vertex']]

avg_sl_shared = np.mean([d['sl_rate'] for _, d in has_shared])
avg_sl_none = np.mean([d['sl_rate'] for _, d in no_shared])
avg_edges_shared = np.mean([d['edges'] for _, d in has_shared])
avg_edges_none = np.mean([d['edges'] for _, d in no_shared])

print("  Shared vertex:  %d pairs, avg SL=%.3f, avg edges=%.1f" %
      (len(has_shared), avg_sl_shared, avg_edges_shared))
print("  No shared:      %d pairs, avg SL=%.3f, avg edges=%.1f" %
      (len(no_shared), avg_sl_none, avg_edges_none))

# Group by total range
print("\n  BY TOTAL RANGE:\n")
by_range = defaultdict(list)
for (i,j), d in pair_data.items():
    by_range[d['total_range']].append(d)

for r in sorted(by_range.keys()):
    pairs = by_range[r]
    avg_sl = np.mean([d['sl_rate'] for d in pairs])
    avg_e = np.mean([d['edges'] for d in pairs])
    names = ", ".join(d['name'] for d in pairs)
    print("  Σrange=%d: %d pairs (%s), avg SL=%.3f, avg edges=%.1f" %
          (r, len(pairs), names, avg_sl, avg_e))

# Does the pair form a TRIANGLE with a base-path arc?
print("\n  TRIANGLE STRUCTURE:")
for (i,j), d in sorted(pair_data.items(), key=lambda x: x[1]['name']):
    vi = set(tile_pairs[i])
    vj = set(tile_pairs[j])
    shared = vi & vj
    union = vi | vj
    if shared:
        # The two tiles share a vertex. With the base-path arcs,
        # do they form a directed triangle?
        other_i = (vi - shared).pop()
        other_j = (vj - shared).pop()
        shared_v = shared.pop()
        # Three vertices: other_i, shared_v, other_j
        # Is there a base-path arc between other_i and other_j?
        bp_exists = abs(other_i - other_j) == 1
        if bp_exists:
            print("    %s: triangle {%d,%d,%d} with base arc (%d,%d)" %
                  (d['name'], other_i, shared_v, other_j,
                   min(other_i, other_j), max(other_i, other_j)))

# ============================================================
banner("d=3: ALL 20 LETTER TRIPLES")
# ============================================================

triple_data = {}
for combo in combinations(range(m), 3):
    i, j, k = combo
    combo_name = tile_names[i] + tile_names[j] + tile_names[k]
    mask = (1 << i) | (1 << j) | (1 << k)

    W = np.zeros((nm, nm), dtype=int)
    for bits in range(total):
        mi = mid_map[tc[bits]]
        mj = mid_map[tc[bits ^ mask]]
        W[mi][mj] += 1

    n_cross = sum(1 for mi2 in range(nm) for mj2 in range(mi2+1, nm) if W[mi2][mj2] > 0)
    sl_rate = sum(W[mi2][mi2] for mi2 in range(nm)) / total

    # Structural features
    verts = set()
    for idx in combo:
        verts.update(tile_pairs[idx])
    n_verts = len(verts)
    ranges = [tile_info[idx]['range'] for idx in combo]

    # Does the triple form a TRIANGLE (all 3 tiles share pairwise vertices)?
    is_triangle = all(
        set(tile_pairs[combo[a]]) & set(tile_pairs[combo[b]])
        for a in range(3) for b in range(a+1, 3)
    )

    triple_data[combo] = {
        'name': combo_name, 'edges': n_cross, 'sl_rate': sl_rate,
        'n_verts': n_verts, 'ranges': tuple(sorted(ranges)),
        'is_triangle': is_triangle
    }

print("  %-8s %-6s %-8s %-10s %-10s %-12s" %
      ("Triple", "Edges", "SL_rate", "#verts", "ranges", "triangle?"))
print("  " + "-"*60)

for combo, d in sorted(triple_data.items(), key=lambda x: -x[1]['sl_rate']):
    print("  %-8s %-6d %-8.3f %-10d %-10s %-12s" %
          (d['name'], d['edges'], d['sl_rate'], d['n_verts'],
           d['ranges'], "YES" if d['is_triangle'] else ""))

# Group by structure
print("\n  BY NUMBER OF DISTINCT VERTICES INVOLVED:")
by_verts = defaultdict(list)
for combo, d in triple_data.items():
    by_verts[d['n_verts']].append(d)

for nv in sorted(by_verts.keys()):
    triples = by_verts[nv]
    avg_sl = np.mean([d['sl_rate'] for d in triples])
    avg_e = np.mean([d['edges'] for d in triples])
    print("  %d vertices: %d triples, avg SL=%.3f, avg edges=%.1f" %
          (nv, len(triples), avg_sl, avg_e))

# Group by range tuple
print("\n  BY RANGE PROFILE:")
by_rp = defaultdict(list)
for combo, d in triple_data.items():
    by_rp[d['ranges']].append(d)

for rp in sorted(by_rp.keys()):
    triples = by_rp[rp]
    avg_sl = np.mean([d['sl_rate'] for d in triples])
    names = ", ".join(d['name'] for d in triples)
    print("  ranges=%s: %d triples (%s), avg SL=%.3f" %
          (rp, len(triples), names, avg_sl))

# ============================================================
banner("THE PATTERN: WHICH LETTER COMBOS ARE MOST/LEAST SILENT?")
# ============================================================

print("""
  DISCOVERY: The self-loop rate depends on the GEOMETRIC STRUCTURE
  of the letter combination in the staircase.

  d=1 (single letters):
    All 6 letters have SL rates in {0.062, 0.125}.
    Range-2 tiles (A,D,F): SL = 0.062 (lower)
    Range-3 tiles (B,E):   SL = 0.125
    Range-4 tile (C):      SL = 0.125

  d=2 (letter pairs):
    Highest SL: pairs whose tiles share a vertex AND have high total range.
    Lowest SL: disjoint tiles (no shared vertex).

  d=3 (letter triples):
    Highest SL: triples involving few distinct vertices (concentrated).
    Lowest SL: triples spanning all 5 vertices (spread out).
    The TRIANGLE {A,D,F} = range (2,2,2) has SL = 0.250.
""")

# ============================================================
banner("d=4,5,6: HIGHER-ORDER LETTER COMBINATIONS")
# ============================================================

for d in [4, 5, 6]:
    print("  d=%d: C(%d,%d) = %d combinations\n" % (d, m, d, comb(m, d)))

    combo_data = {}
    for combo in combinations(range(m), d):
        combo_name = "".join(tile_names[k] for k in combo)
        mask = sum(1 << k for k in combo)
        sl = sum(1 for bits in range(total) if mid_map[tc[bits]] == mid_map[tc[bits ^ mask]])
        sl_rate = sl / total

        verts = set()
        for idx in combo:
            verts.update(tile_pairs[idx])
        ranges = tuple(sorted(tile_info[idx]['range'] for idx in combo))

        combo_data[combo] = {
            'name': combo_name, 'sl_rate': sl_rate,
            'n_verts': len(verts), 'ranges': ranges
        }

    # Sort by SL rate
    sorted_combos = sorted(combo_data.items(), key=lambda x: -x[1]['sl_rate'])
    print("  TOP 5 by SL rate:")
    for combo, cd in sorted_combos[:5]:
        print("    %s: SL=%.3f, verts=%d, ranges=%s" %
              (cd['name'], cd['sl_rate'], cd['n_verts'], cd['ranges']))
    print("  BOTTOM 5:")
    for combo, cd in sorted_combos[-5:]:
        print("    %s: SL=%.3f, verts=%d, ranges=%s" %
              (cd['name'], cd['sl_rate'], cd['n_verts'], cd['ranges']))

    # By vertex count
    by_verts = defaultdict(list)
    for combo, cd in combo_data.items():
        by_verts[cd['n_verts']].append(cd)
    print("\n  BY #VERTICES:")
    for nv in sorted(by_verts.keys()):
        items = by_verts[nv]
        avg_sl = np.mean([c['sl_rate'] for c in items])
        print("    %d verts: %d combos, avg SL=%.3f" % (nv, len(items), avg_sl))
    print()

# ============================================================
banner("THE VERTEX-COUNT LAW")
# ============================================================

print("""
  THE KEY PATTERN ACROSS ALL d:

  Self-loop rate depends primarily on HOW MANY DISTINCT VERTICES
  are involved in the flipped tiles.

  Fewer vertices = more concentrated change = more likely silent.
  More vertices = more spread-out change = less likely silent.

  This makes sense: a silent mutation requires the resulting
  tournament to be isomorphic to the original. If only 3 vertices
  are affected (out of 5), the remaining 2 anchor the structure.
  If all 5 are affected, nothing is anchored.
""")

# Collect vertex-count vs SL across all d
print("  SL RATE BY (d, #vertices):\n")
print("  %-5s" % "d", end="")
for nv in range(2, n+1):
    print(" %d-vert  " % nv, end="")
print()
print("  " + "-"*50)

for d in range(1, m+1):
    print("  d=%-3d" % d, end="")
    for nv in range(2, n+1):
        # Compute average SL for combos at this d with this vertex count
        sl_vals = []
        for combo in combinations(range(m), d):
            verts = set()
            for idx in combo:
                verts.update(tile_pairs[idx])
            if len(verts) == nv:
                mask = sum(1 << k for k in combo)
                sl = sum(1 for bits in range(total) if mid_map[tc[bits]] == mid_map[tc[bits ^ mask]])
                sl_vals.append(sl / total)
        if sl_vals:
            print(" %.3f(%d)" % (np.mean(sl_vals), len(sl_vals)), end="")
        else:
            print(" %-9s" % "-", end="")
    print()

# ============================================================
banner("SPECIFIC CLASS-PAIR CONNECTIONS BY LETTER COMBO")
# ============================================================

# For each class pair at d=2, which letter pairs connect them?
# And for d=3 pairs, which letter triples?

print("  d=2 PAIRS: which letter combos connect which class pairs?\n")

# Build a table: for each pair of merged classes not connected at d=1,
# show which d=2 letter pairs create the connection
for mi in range(nm):
    for mj in range(mi+1, nm):
        # Check d=1 connectivity
        has_d1 = False
        for bits in range(total):
            if mid_map[tc[bits]] != mi: continue
            for k in range(m):
                if mid_map[tc[bits ^ (1<<k)]] == mj:
                    has_d1 = True; break
            if has_d1: break

        if has_d1: continue  # Skip d=1 connected pairs

        # Find which d=2 letter pairs connect mi ↔ mj
        connecting_pairs = []
        for i in range(m):
            for j in range(i+1, m):
                mask = (1<<i) | (1<<j)
                for bits in range(total):
                    if mid_map[tc[bits]] == mi and mid_map[tc[bits ^ mask]] == mj:
                        connecting_pairs.append(tile_names[i] + tile_names[j])
                        break
                    if mid_map[tc[bits]] == mj and mid_map[tc[bits ^ mask]] == mi:
                        connecting_pairs.append(tile_names[i] + tile_names[j])
                        break

        if connecting_pairs:
            print("  {%d ↔ %d}: connected by d=2 pairs: %s" %
                  (mi, mj, ", ".join(sorted(set(connecting_pairs)))))

# For d=3 exclusive connections
print("\n  d=3 EXCLUSIVE: class pairs ONLY reachable at d=3:\n")
for mi in range(nm):
    for mj in range(mi+1, nm):
        has_d1 = any(mid_map[tc[bits ^ (1<<k)]] == mj
                     for bits in range(total) if mid_map[tc[bits]] == mi
                     for k in range(m))
        has_d2 = any(mid_map[tc[bits ^ ((1<<i)|(1<<j))]] == mj
                     for bits in range(total) if mid_map[tc[bits]] == mi
                     for i in range(m) for j in range(i+1, m))

        if has_d1 or has_d2: continue

        # Find d=3 combos
        connecting_triples = []
        for combo in combinations(range(m), 3):
            mask = sum(1 << k for k in combo)
            for bits in range(total):
                if mid_map[tc[bits]] == mi and mid_map[tc[bits ^ mask]] == mj:
                    connecting_triples.append("".join(tile_names[k] for k in combo))
                    break
                if mid_map[tc[bits]] == mj and mid_map[tc[bits ^ mask]] == mi:
                    connecting_triples.append("".join(tile_names[k] for k in combo))
                    break

        if connecting_triples:
            print("  {%d ↔ %d}: ONLY by d=3 triples: %s" %
                  (mi, mj, ", ".join(sorted(set(connecting_triples)))))

print("\nDone. opus-2026-03-24-S302")
