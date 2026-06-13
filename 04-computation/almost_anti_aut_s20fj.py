#!/usr/bin/env python3
"""
almost_anti_aut_s20fj.py — Which SC tournaments have almost-anti-automorphisms?
kind-pasteur-2026-03-24-S20fj

An "almost-anti-automorphism" of T is sigma in S_n with defect(T, sigma) = C(n,2)-1.
This means sigma(T) agrees with T at exactly 1 arc.
Equivalently: sigma(T) and T^c differ at exactly 1 arc.

At odd n: ONLY SC tournaments have this property (SC filtering theorem).
Question: WHICH SC tournaments? What characterizes them?

Also: at n=6 (even), which NON-SC tournaments have this? How do they
relate to the 5 overlap pairs?

Compute per-class: does this class have an almost-anti-automorphism?
Correlate with H, |Aut|, scores, position in metagraph.
"""

import sys
from math import factorial, comb
from itertools import permutations
from collections import defaultdict
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  ALMOST-ANTI-AUTOMORPHISM CLASSIFICATION")
print("  kind-pasteur-2026-03-24-S20fj")
print("=" * 80)

def get_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))
    return tiles

for n in [5, 6, 7]:
    t0 = time.time()
    N = n
    M = comb(n, 2)
    TILES = get_tiles(n)
    m = len(TILES)
    VERTS = list(range(n, 0, -1))
    all_perms = list(permutations(range(N)))

    tile_verts = []
    for i, (x, y) in enumerate(TILES):
        tile_verts.append((VERTS.index(x), VERTS.index(y)))

    def bits_to_adj(bits):
        A = [[0]*N for _ in range(N)]
        for k in range(N-1): A[k][k+1] = 1
        for i in range(m):
            xi, yi = tile_verts[i]
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

    def count_aut(A):
        return sum(1 for p in all_perms if all(A[p[i]][p[j]] == A[i][j] for i in range(N) for j in range(N)))

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

    def compute_defect(A, sigma):
        sigma_inv = [0]*N
        for i in range(N): sigma_inv[sigma[i]] = i
        d = 0
        for i in range(N):
            for j in range(i+1, N):
                if A[sigma_inv[i]][sigma_inv[j]] != A[i][j]:
                    d += 1
        return d

    def has_almost_anti_aut(A):
        """Check if T has sigma with defect = M-1."""
        for sigma in all_perms:
            if compute_defect(A, sigma) == M - 1:
                return True
        return False

    print(f"\n{'#'*60}")
    print(f"  n = {n}, M = {M}, target = {M-1}")
    print(f"{'#'*60}")

    # Build iso classes
    canon_map = {}
    for mask in range(1 << m):
        bits = [(mask >> k) & 1 for k in range(m)]
        A = bits_to_adj(bits)
        canon_map[mask] = canonicalize(A)

    class_reps = {}
    for mask in range(1 << m):
        cn = canon_map[mask]
        if cn not in class_reps:
            class_reps[cn] = mask

    # For each class: compute properties and check almost-anti-aut
    V = len(class_reps)
    results = []

    for cn, mask in class_reps.items():
        bits = [(mask >> k) & 1 for k in range(m)]
        A = bits_to_adj(bits)
        comp_cn = canonicalize(adj_complement(A))
        is_sc = (comp_cn == cn)
        sc = tuple(sorted(sum(A[i]) for i in range(N)))
        palindromic = all(sc[i] + sc[N-1-i] == N-1 for i in range(N//2))
        H = count_hp(A)
        aut = count_aut(A)

        has_aaa = has_almost_anti_aut(A)

        results.append({
            'cn': cn, 'sc': is_sc, 'H': H, 'aut': aut,
            'scores': sc, 'palindromic': palindromic, 'aaa': has_aaa
        })

    # Classify
    sc_with_aaa = [r for r in results if r['sc'] and r['aaa']]
    sc_without_aaa = [r for r in results if r['sc'] and not r['aaa']]
    nsc_with_aaa = [r for r in results if not r['sc'] and r['aaa']]
    nsc_without_aaa = [r for r in results if not r['sc'] and not r['aaa']]

    print(f"\n  Classification:")
    print(f"    SC with almost-anti-aut:  {len(sc_with_aaa)}")
    print(f"    SC without:               {len(sc_without_aaa)}")
    print(f"    NSC with almost-anti-aut: {len(nsc_with_aaa)}")
    print(f"    NSC without:              {len(nsc_without_aaa)}")
    print(f"    Total: {V}")

    # Detail for classes with almost-anti-aut
    print(f"\n  Classes WITH almost-anti-automorphism:")
    print(f"  {'H':>6} {'|Aut|':>6} {'SC':>3} {'Scores':>25}")
    for r in sorted(sc_with_aaa + nsc_with_aaa, key=lambda x: x['H']):
        print(f"  {r['H']:6d} {r['aut']:6d} {'Y' if r['sc'] else 'N':>3} {str(r['scores'])}")

    # What distinguishes SC with vs without aaa?
    if sc_with_aaa and sc_without_aaa:
        H_with = [r['H'] for r in sc_with_aaa]
        H_without = [r['H'] for r in sc_without_aaa]
        aut_with = [r['aut'] for r in sc_with_aaa]
        aut_without = [r['aut'] for r in sc_without_aaa]

        print(f"\n  SC CLASSES: WITH vs WITHOUT almost-anti-aut:")
        print(f"    H range WITH:    {min(H_with)}-{max(H_with)} (avg {sum(H_with)/len(H_with):.1f})")
        print(f"    H range WITHOUT: {min(H_without)}-{max(H_without)} (avg {sum(H_without)/len(H_without):.1f})")
        print(f"    |Aut| WITH:    {sorted(aut_with)}")
        print(f"    |Aut| WITHOUT: {sorted(aut_without)}")

    # At even n: which NSC have aaa?
    if nsc_with_aaa:
        print(f"\n  NON-SC classes with almost-anti-aut (= overlap pairs):")
        for r in sorted(nsc_with_aaa, key=lambda x: x['H']):
            print(f"    H={r['H']}, |Aut|={r['aut']}, scores={r['scores']}")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")

print("\nDONE.")
print("=" * 80)
