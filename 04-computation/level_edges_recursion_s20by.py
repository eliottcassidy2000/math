#!/usr/bin/env python3
"""
level_edges_recursion_s20by.py -- kind-pasteur-2026-03-22-S20by

CLOSING TWO MORE GAPS:
GAP 3: Level edge formula (0, 0, 1, 15)
GAP 5: The n -> n-2 recursion

LEVEL EDGES: edges between iso classes with the SAME H value.
These are the "non-DAG" edges -- the degenerate cases where
an arc flip changes the iso class but NOT H.

For this to happen: flip(T, e) must have:
  - Different iso class from T (not a self-loop)
  - Same H value (the flip preserves H)

From the OCF: H = 1 + 2*alpha_1 + 4*alpha_2 + ...
So level edges require: the flip changes the cycle structure
but keeps the TOTAL OCF value unchanged.

THE n -> n-2 RECURSION: The PoS (most ambiguous) score class at n
contains tournaments that, when you remove the source and sink,
give tournaments on n-2 vertices. Does this give a map G_n -> G_{n-2}?

Author: kind-pasteur-2026-03-22-S20by
"""
import sys
import numpy as np
from math import comb, factorial
from collections import defaultdict
from itertools import permutations
sys.stdout.reconfigure(line_buffering=True)

def count_hp(A, n):
    dp = defaultdict(int)
    for v in range(n): dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[(mask, v)] == 0: continue
            for w in range(n):
                if mask & (1 << w): continue
                if A[v][w]: dp[(mask | (1 << w), w)] += dp[(mask, v)]
    return sum(dp[((1 << n) - 1, v)] for v in range(n))

def canonical(A, n):
    best = None
    for perm in permutations(range(n)):
        form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best: best = form
    return best

print("=" * 70)
print("  CLOSING TWO GAPS: LEVEL EDGES + n->n-2 RECURSION")
print("=" * 70)

# ================================================================
# GAP 3: LEVEL EDGES
# ================================================================
print(f"\n{'='*70}")
print(f"  GAP 3: LEVEL EDGE ANALYSIS")
print(f"{'='*70}\n")

for n in [5, 6]:
    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)

    # Build iso classes
    canon_map = defaultdict(list)
    H_map = {}
    for bits in range(2**m):
        A = np.zeros((n,n), dtype=np.int8)
        for k, (i,j) in enumerate(pairs):
            if (bits >> k) & 1: A[i][j] = 1
            else: A[j][i] = 1
        H = count_hp(A, n)
        H_map[bits] = H
        cf = canonical(A, n)
        canon_map[cf].append(bits)

    classes = []
    cf_to_id = {}
    for cf, members in sorted(canon_map.items(), key=lambda x: H_map[x[1][0]]):
        cid = len(classes)
        cf_to_id[cf] = cid
        classes.append({'id': cid, 'H': H_map[members[0]], 'size': len(members),
                       'members': set(members)})

    bits_to_class = {}
    for c in classes:
        for b in c['members']:
            bits_to_class[b] = c['id']

    # Find all level edges
    level_edges = set()
    for c in classes:
        T = list(c['members'])[0]
        for k in range(m):
            nb = T ^ (1 << k)
            nb_class = bits_to_class[nb]
            if nb_class != c['id'] and classes[nb_class]['H'] == c['H']:
                edge = (min(c['id'], nb_class), max(c['id'], nb_class))
                level_edges.add(edge)

    print(f"  n={n}: {len(level_edges)} level edges")

    # Analyze each level edge
    for i, j in sorted(level_edges):
        ci, cj = classes[i], classes[j]
        print(f"    ({i},{j}): H={ci['H']}, sizes=({ci['size']},{cj['size']})")

    # Which H values have level edges?
    level_H = set()
    for i, j in level_edges:
        level_H.add(classes[i]['H'])
    print(f"  Level edges occur at H values: {sorted(level_H)}")

    # How many iso classes at each level-H?
    for h in sorted(level_H):
        classes_at_h = [c for c in classes if c['H'] == h]
        n_classes = len(classes_at_h)
        # How many level edges at this H?
        edges_at_h = [(i,j) for i,j in level_edges if classes[i]['H'] == h]
        print(f"    H={h}: {n_classes} classes, {len(edges_at_h)} level edges")

    print()

# ================================================================
# GAP 5: THE n -> n-2 RECURSION
# ================================================================
print(f"{'='*70}")
print(f"  GAP 5: THE n -> n-2 RECURSION")
print(f"{'='*70}\n")

# At n=5: the PoS class (1,2,2,2,3) has 3 iso classes with H=11,13,15.
# Remove the source (score 1) and sink (score 3) vertices.
# The remaining 3 vertices form an (n-2)=3 tournament.
# Does this give a map from PoS iso classes at n=5 to iso classes at n=3?

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)

# Rebuild classes
canon_map = defaultdict(list)
H_map = {}
for bits in range(2**m):
    A = np.zeros((n,n), dtype=np.int8)
    for k, (i,j) in enumerate(pairs):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    H = count_hp(A, n)
    H_map[bits] = H
    s = tuple(sorted(A.sum(axis=1).astype(int)))
    cf = canonical(A, n)
    canon_map[cf].append((bits, s))

classes = []
for cf, members in sorted(canon_map.items(), key=lambda x: H_map[x[0][0]]):
    bits0, score0 = members[0]
    cid = len(classes)
    classes.append({'id': cid, 'H': H_map[bits0], 'score': score0,
                   'size': len(members), 'members': [(b,s) for b,s in members]})

# PoS class: score (1,2,2,2,3)
pos_classes = [c for c in classes if c['score'] == (1,2,2,2,3)]
print(f"  n={n}: PoS class (1,2,2,2,3) has {len(pos_classes)} iso classes")

for c in pos_classes:
    print(f"    Class {c['id']}: H={c['H']}, size={c['size']}")

    # For a representative, find source (score 1) and sink (score 3)
    bits, score = c['members'][0]
    A = np.zeros((n,n), dtype=np.int8)
    for k, (i,j) in enumerate(pairs):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1

    scores = A.sum(axis=1).astype(int)
    # Score 1 = almost-source (beats 1, loses to 3)
    # Score 3 = almost-sink (beats 3, loses to 1)
    source = list(scores).index(1)  # first vertex with score 1
    sink = [v for v in range(n) if scores[v] == 3][0]  # first with score 3

    # Extract the middle tournament (n-2 = 3 vertices)
    middle_verts = [v for v in range(n) if v != source and v != sink]
    n_mid = len(middle_verts)
    A_mid = np.zeros((n_mid, n_mid), dtype=np.int8)
    for i_m, v in enumerate(middle_verts):
        for j_m, w in enumerate(middle_verts):
            A_mid[i_m][j_m] = A[v][w]

    H_mid = count_hp(A_mid, n_mid)
    cf_mid = canonical(A_mid, n_mid)

    print(f"      Source={source}(score {scores[source]}), Sink={sink}(score {scores[sink]})")
    print(f"      Middle vertices: {middle_verts}, scores: {[scores[v] for v in middle_verts]}")
    print(f"      Middle tournament H = {H_mid}")
    print(f"      Middle canonical: {cf_mid[:9]}...")

# The n=3 iso classes:
n3 = 3
pairs3 = [(i,j) for i in range(n3) for j in range(i+1, n3)]
canon3_map = defaultdict(list)
for bits in range(2**len(pairs3)):
    A = np.zeros((n3,n3), dtype=np.int8)
    for k, (i,j) in enumerate(pairs3):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    H = count_hp(A, n3)
    cf = canonical(A, n3)
    canon3_map[cf].append(H)

print(f"\n  n=3 iso classes:")
for cf, Hs in canon3_map.items():
    print(f"    H={Hs[0]}, size={len(Hs)}")

print(f"""
  THE RECURSION:
  PoS iso classes at n=5 map to iso classes at n=3 via "remove source and sink."
  - PoS class with H=11 -> middle H=?
  - PoS class with H=13 -> middle H=?
  - PoS class with H=15 -> middle H=?

  If H_mid determines the PoS H: we have a FORMULA connecting levels.
  The source-sink embedding formula: H_n = (2^k + 1) * H_{n-2}
  where k = (n-1)/2 = 2 at n=5. So H_5 = 5 * H_3.
  H_3 = 1: H_5 = 5. H_3 = 3: H_5 = 15.
  But PoS has H = 11, 13, 15 -- not just 5 and 15.
  So the recursion is NOT just multiplication by 2^k+1.

  The CORRECTION: the source-sink formula works for the T->S extreme
  (sink beats source), giving H = (2^k+1)*H_mid = 5*H_mid.
  For S->T (source beats sink): H = 2*H_mid + 2n-1 (different formula).
  The intermediate cases (PoS with source-sink in various configurations)
  give H values between these extremes.
""")
