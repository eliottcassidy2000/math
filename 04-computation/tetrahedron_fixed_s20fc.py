#!/usr/bin/env python3
"""
tetrahedron_fixed_s20fc.py — FIXED tetrahedron (correct |Aut| computation)
kind-pasteur-2026-03-24-S20fc

BUG FIX: |Aut| must be computed by counting automorphisms directly,
not from n!/size (which gives n!*|Aut|/H, not |Aut|).

The fiber identity: tilings(C) = H(C)/|Aut(C)|
So: |Aut(C)| = H(C) / tilings(C) = H(C) / size(C)
where size = #{tilings mapping to class C}
"""

import sys
from math import factorial, comb
from itertools import permutations
from collections import defaultdict
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  FIXED TOURNAMENT TETRAHEDRON")
print("  kind-pasteur-2026-03-24-S20fc")
print("=" * 80)

def get_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))
    return tiles

for n in [4, 5, 6]:
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

    def count_aut(A):
        """Count actual automorphisms of tournament A."""
        count = 0
        for p in all_perms:
            is_aut = True
            for i in range(N):
                for j in range(N):
                    if A[p[i]][p[j]] != A[i][j]:
                        is_aut = False
                        break
                if not is_aut: break
            if is_aut: count += 1
        return count

    def score_seq(A):
        return tuple(sorted(sum(A[i]) for i in range(N)))

    # Build all data
    canon_map = {}; adj_map = {}
    for mask in range(1 << m):
        bits = [(mask >> k) & 1 for k in range(m)]
        A = bits_to_adj(bits)
        adj_map[mask] = A
        canon_map[mask] = canonicalize(A)

    # Unmerged classes with CORRECT properties
    class_data = {}
    for cn in set(canon_map.values()):
        if cn in class_data: continue
        for mask, c in canon_map.items():
            if c == cn:
                A = adj_map[mask]
                comp_cn = canonicalize(adj_complement(A))
                size = sum(1 for m2, c2 in canon_map.items() if c2 == cn)
                H = count_hp(A)
                aut = count_aut(A)  # CORRECT |Aut|
                sc = score_seq(A)
                class_data[cn] = {
                    'comp': comp_cn, 'aut': aut, 'H': H,
                    'scores': sc, 'size': size, 'sc': comp_cn == cn
                }
                # Verify fiber identity
                assert size == H // aut or abs(size - H/aut) < 0.001, f"Fiber identity fails: size={size}, H/aut={H/aut}"
                break

    V = len(class_data)
    SC = sum(1 for d in class_data.values() if d['sc'])

    print(f"\n{'#'*60}")
    print(f"  n = {n}, m = {m}, V = {V}, SC = {SC}")
    print(f"{'#'*60}")

    # VERIFY THREE IDENTITIES with correct |Aut|
    id1 = sum(d['H'] / d['aut'] for d in class_data.values())
    id2 = sum(factorial(N) / d['aut'] for d in class_data.values())
    id3 = sum(d['H'] * factorial(N) / d['aut'] for d in class_data.values())

    exp1 = 2 ** comb(n-1, 2)
    exp2 = 2 ** comb(n, 2)
    exp3 = factorial(n) * exp1

    print(f"\n  THREE IDENTITIES (corrected |Aut|):")
    print(f"    1. TILING:  sum H/|Aut| = {id1:.0f} (expected {exp1}) {'OK' if abs(id1 - exp1) < 0.5 else 'FAIL'}")
    print(f"    2. LABELED: sum n!/|Aut| = {id2:.0f} (expected {exp2}) {'OK' if abs(id2 - exp2) < 0.5 else 'FAIL'}")
    print(f"    3. SZELE:   sum H*n!/|Aut| = {id3:.0f} (expected {exp3}) {'OK' if abs(id3 - exp3) < 0.5 else 'FAIL'}")

    # Per-class table
    print(f"\n  {'H':>6} {'|Aut|':>6} {'Tilings':>8} {'n!/|Aut|':>10} {'Score':>25} {'SC':>3}")
    for cn in sorted(class_data.keys(), key=lambda c: class_data[c]['H']):
        d = class_data[cn]
        print(f"  {d['H']:6d} {d['aut']:6d} {d['size']:8d} {factorial(N)//d['aut']:10d} {str(d['scores']):>25} {'Y' if d['sc'] else 'N':>3}")

    # ================================================================
    # OVERLAP ANALYSIS AT n=6
    # ================================================================
    if n == 6:
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

        print(f"\n  OVERLAP AT n=6: {len(overlap)} pairs")
        print(f"\n    OVERLAP (complement + wiggly neighbor):")
        for pair in sorted(overlap, key=lambda p: class_data[p[0]]['H']):
            d = class_data[pair[0]]
            print(f"      H={d['H']:3d}, |Aut|={d['aut']}, tilings={d['size']}, scores={d['scores']}")

        print(f"\n    NON-OVERLAP complement pairs:")
        for pair in sorted(comp_pairs - overlap, key=lambda p: class_data[p[0]]['H']):
            d = class_data[pair[0]]
            print(f"      H={d['H']:3d}, |Aut|={d['aut']}, tilings={d['size']}, scores={d['scores']}")

        # KEY DISTINGUISHER: overlap pairs have PALINDROMIC scores
        overlap_palindromic = all(
            d['scores'][i] + d['scores'][N-1-i] == N-1
            for pair in overlap
            for d in [class_data[pair[0]]]
            for i in range(N//2)
        )
        nonoverlap_palindromic = [
            all(class_data[pair[0]]['scores'][i] + class_data[pair[0]]['scores'][N-1-i] == N-1 for i in range(N//2))
            for pair in comp_pairs - overlap
        ]
        print(f"\n    All overlap palindromic? {overlap_palindromic}")
        print(f"    Non-overlap palindromic count: {sum(nonoverlap_palindromic)}/{len(nonoverlap_palindromic)}")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")

print("\nDONE.")
print("=" * 80)
