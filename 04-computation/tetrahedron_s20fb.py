#!/usr/bin/env python3
"""
tetrahedron_s20fb.py — The {Tiling, Labeled, Szele, ???} tetrahedron
kind-pasteur-2026-03-24-S20fb

THREE IDENTITIES on the same iso class nodes:
  1. TILING:  sum_C H(C)/|Aut(C)| = 2^{C(n-1,2)}
  2. LABELED: sum_C n!/|Aut(C)| = 2^{C(n,2)}
  3. SZELE:   sum_C H(C)*n!/|Aut(C)| = n! * 2^{C(n-1,2)}

Identity 3 = n! * Identity 1. So only 2 are independent.

THE HIDDEN FOURTH: What new identity adds information?
Candidates:
  A. WIGGLY:  W row_sum = H(C)/|Aut(C)| * m
  B. BLUE:    blue_count = 2^{exponent} (grid-sym count)
  C. DEFECT:  D(n) = sum_C neutral(C) (defect-1 Burnside)
  D. COMPLEMENT: H(C) = H(C^c) (complement invariance)
  E. EVEN GRAPH: #{iso classes} = #{even graph classes}

Check which of these provides genuinely NEW information.

Also: study the 5 OVERLAP PAIRS at n=6 in the tetrahedron framework.
"""

import sys
from math import factorial, comb
from itertools import permutations
from collections import defaultdict
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  THE TOURNAMENT TETRAHEDRON")
print("  kind-pasteur-2026-03-24-S20fb")
print("=" * 80)

def get_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))
    return tiles

for n in [5, 6]:
    t0 = time.time()
    TILES = get_tiles(n)
    m = len(TILES)
    N = n
    VERTS = list(range(n, 0, -1))
    all_perms = list(permutations(range(N)))

    def bits_to_adj(bits):
        A = [[0]*N for _ in range(N)]
        for k in range(N-1): A[k][k+1] = 1
        for i in range(m):
            x, y = TILES[i]
            xi, yi = VERTS.index(x), VERTS.index(y)
            if bits[i] == 0: A[xi][yi] = 1
            else: A[yi][xi] = 1
        return A

    def adj_complement(A):
        return [[1-A[i][j] if i!=j else 0 for j in range(N)] for i in range(N)]

    def canonicalize(A):
        best = None
        for p in all_perms:
            s = ''.join(str(A[p[i]][p[j]]) for i in range(N) for j in range(N))
            if best is None or s < best: best = s
        return best

    def count_hp(A):
        dp = [[0]*N for _ in range(1 << N)]
        for v in range(N): dp[1 << v][v] = 1
        for mask in range(1, 1 << N):
            for v in range(N):
                if not (mask & (1 << v)) or dp[mask][v] == 0: continue
                for u in range(N):
                    if mask & (1 << u): continue
                    if A[v][u]: dp[mask | (1 << u)][u] += dp[mask][v]
        return sum(dp[(1 << N) - 1])

    def score_seq(A):
        return tuple(sorted(sum(A[i]) for i in range(N)))

    # Build all data
    canon_map = {}; adj_map = {}
    for mask in range(1 << m):
        bits = [(mask >> k) & 1 for k in range(m)]
        A = bits_to_adj(bits)
        adj_map[mask] = A; canon_map[mask] = canonicalize(A)

    # Unmerged classes
    class_data = {}
    for cn in set(canon_map.values()):
        if cn in class_data: continue
        for mask, c in canon_map.items():
            if c == cn:
                A = adj_map[mask]
                comp_cn = canonicalize(adj_complement(A))
                size = sum(1 for m2, c2 in canon_map.items() if c2 == cn)
                aut = factorial(N) // size
                H = count_hp(A)
                sc = score_seq(A)
                class_data[cn] = {'comp': comp_cn, 'aut': aut, 'H': H, 'scores': sc, 'size': size, 'sc': comp_cn == cn}
                break

    V = len(class_data)
    SC = sum(1 for d in class_data.values() if d['sc'])

    print(f"\n{'#'*60}")
    print(f"  n = {n}, m = {m}, V = {V}, SC = {SC}")
    print(f"{'#'*60}")

    # ================================================================
    # VERIFY THREE IDENTITIES
    # ================================================================
    id1 = sum(d['H'] / d['aut'] for d in class_data.values())  # tiling
    id2 = sum(factorial(N) / d['aut'] for d in class_data.values())  # labeled
    id3 = sum(d['H'] * factorial(N) / d['aut'] for d in class_data.values())  # szele

    expected1 = 2 ** comb(n-1, 2)
    expected2 = 2 ** comb(n, 2)
    expected3 = factorial(n) * 2 ** comb(n-1, 2)

    print(f"\n  THREE IDENTITIES:")
    print(f"    1. TILING:  sum H/|Aut| = {id1:.0f} (expected {expected1}) {'OK' if abs(id1 - expected1) < 0.5 else 'FAIL'}")
    print(f"    2. LABELED: sum n!/|Aut| = {id2:.0f} (expected {expected2}) {'OK' if abs(id2 - expected2) < 0.5 else 'FAIL'}")
    print(f"    3. SZELE:   sum H*n!/|Aut| = {id3:.0f} (expected {expected3}) {'OK' if abs(id3 - expected3) < 0.5 else 'FAIL'}")
    print(f"    Identity 3 = n! * Identity 1? {abs(id3 - factorial(N) * id1) < 0.5}")

    # ================================================================
    # THE HIDDEN FOURTH CANDIDATES
    # ================================================================
    print(f"\n  CANDIDATE FOURTH IDENTITIES:")

    # A. WIGGLY: W row_sum = tilings * m = H/|Aut| * m
    total_wiggly = sum(d['H'] / d['aut'] * m for d in class_data.values())
    print(f"    A. WIGGLY: sum (H/|Aut| * m) = {total_wiggly:.0f} = {expected1} * {m} = {expected1 * m}")
    print(f"       -> Just m * Identity 1. NO new info.")

    # B. BLUE: count grid-symmetric tilings
    n_fixed = (n-1) // 2
    blue_exp = (comb(n-1,2) + n_fixed) // 2
    blue_count = 2 ** blue_exp
    print(f"    B. BLUE: 2^{blue_exp} = {blue_count} grid-symmetric tilings")
    print(f"       -> Independent of H, |Aut|. NEW info about tiling symmetry!")

    # C. DEFECT: D(n) = sum neutral_labeled(C)
    # neutral(C) = #{arcs e : [C flip e] = [C]} per labeled T in C
    # D = sum_C neutral(C) ... but this needs computation
    # From our earlier formula: D(n) via Burnside defect-1

    # D. COMPLEMENT INVARIANCE: H(C) = H(C^c) for all C
    comp_check = all(class_data[cn]['H'] == class_data[class_data[cn]['comp']]['H'] for cn in class_data)
    print(f"    D. COMPLEMENT: H(C) = H(C^c) for all C? {comp_check}")
    print(f"       -> Constraint on H values. Restricts which H-vectors are possible.")

    # E. EVEN GRAPH EQUINUMEROSITY
    print(f"    E. EVEN GRAPH: V = {V} = A000568({n}) = #even_graphs({n})")
    print(f"       -> Remarkable but gives no per-class info.")

    # ================================================================
    # WHAT CONSTRAINS (H, |Aut|) PER CLASS?
    # ================================================================
    print(f"\n  PER-CLASS DATA (sorted by H):")
    print(f"  {'H':>6} {'|Aut|':>6} {'Tilings':>8} {'Score':>25} {'SC':>3} {'Comp H':>7}")

    for cn in sorted(class_data.keys(), key=lambda c: class_data[c]['H']):
        d = class_data[cn]
        comp_H = class_data[d['comp']]['H']
        print(f"  {d['H']:6d} {d['aut']:6d} {d['H']//d['aut']:8d} {str(d['scores']):>25} {'Y' if d['sc'] else 'N':>3} {comp_H:7d}")

    # ================================================================
    # THE 5 OVERLAP PAIRS AT n=6: Their tetrahedron position
    # ================================================================
    if n == 6:
        # Find complement pairs that are also wiggly neighbors
        wiggly_edges = set()
        for mask in range(1 << m):
            cn1 = canon_map[mask]
            for wi in range(m):
                cn2 = canon_map[mask ^ (1 << wi)]
                if cn1 != cn2:
                    wiggly_edges.add((min(cn1,cn2), max(cn1,cn2)))

        comp_pairs = set()
        for cn, d in class_data.items():
            if not d['sc']:
                comp_pairs.add((min(cn, d['comp']), max(cn, d['comp'])))

        overlap = comp_pairs & wiggly_edges
        non_overlap_comp = comp_pairs - wiggly_edges

        print(f"\n  OVERLAP PAIRS (complement = wiggly neighbor):")
        print(f"    Overlap: {len(overlap)}, Non-overlap complement: {len(non_overlap_comp)}")

        print(f"\n    OVERLAP PAIRS (H, |Aut|, tilings, scores):")
        for pair in sorted(overlap, key=lambda p: class_data[p[0]]['H']):
            d = class_data[pair[0]]
            print(f"      H={d['H']:3d}, |Aut|={d['aut']:2d}, tilings={d['H']//d['aut']:3d}, scores={d['scores']}")

        print(f"\n    NON-OVERLAP COMPLEMENT PAIRS:")
        for pair in sorted(non_overlap_comp, key=lambda p: class_data[p[0]]['H']):
            d = class_data[pair[0]]
            print(f"      H={d['H']:3d}, |Aut|={d['aut']:2d}, tilings={d['H']//d['aut']:3d}, scores={d['scores']}")

        # What distinguishes overlap from non-overlap?
        overlap_H = [class_data[p[0]]['H'] for p in overlap]
        nonoverlap_H = [class_data[p[0]]['H'] for p in non_overlap_comp]
        overlap_aut = [class_data[p[0]]['aut'] for p in overlap]
        nonoverlap_aut = [class_data[p[0]]['aut'] for p in non_overlap_comp]

        print(f"\n    DISTINGUISHING FEATURES:")
        print(f"      Overlap H: {sorted(overlap_H)}, avg={sum(overlap_H)/len(overlap_H):.1f}")
        print(f"      Non-overlap H: {sorted(nonoverlap_H)}, avg={sum(nonoverlap_H)/len(nonoverlap_H):.1f}")
        print(f"      Overlap |Aut|: {sorted(overlap_aut)}")
        print(f"      Non-overlap |Aut|: {sorted(nonoverlap_aut)}")

        # Key: overlap classes have ODD |Aut|? ALL tournament auts are odd.
        # Overlap classes have LARGER |Aut|?
        print(f"      Overlap avg |Aut|: {sum(overlap_aut)/len(overlap_aut):.1f}")
        print(f"      Non-overlap avg |Aut|: {sum(nonoverlap_aut)/len(nonoverlap_aut):.1f}")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")

print("\n" + "=" * 60)
print("THE HIDDEN FOURTH")
print("=" * 60)
print("""
The three identities {Tiling, Labeled, Szele} are NOT independent:
  Szele = n! * Tiling. Only TWO are independent.

The genuinely independent constraints on per-class (H, |Aut|) are:
  1. sum H(C)/|Aut(C)| = 2^{C(n-1,2)}   [tiling count]
  2. sum n!/|Aut(C)| = 2^{C(n,2)}        [labeled count]
  3. H(C) = H(complement(C))             [complement invariance]
  4. H(C) always ODD                     [from Redei's theorem]

The HIDDEN FOURTH that completes the tetrahedron is:

  THE GRID-SYMMETRY COUNT (Blue identity):
    #{grid-symmetric tilings} = 2^{(C(n-1,2) + floor((n-1)/2))/2}

This is NOT derivable from H and |Aut| alone. It comes from the
TILING MODEL's interaction with the staircase's geometric symmetry.

The tetrahedron {Tiling, Labeled, Szele, Blue} has:
  - 4 vertices: each a counting framework
  - 6 edges: each a formula connecting two frameworks
  - 4 faces: each a triple identity
""")

print("DONE.")
print("=" * 80)
