#!/usr/bin/env python3
"""
Verify: does t3 parity determine the bipartition at n=7?

We already have the n=7 skeleton edges. We need:
1. The bipartition coloring (from BFS)
2. The t3 value for each SC class

For t3, we need the class representatives. Use refined invariant to
build classes faster, but only need SC classes that appear in skeleton.

kind-pasteur-2026-03-06-S25h
"""
import re
from itertools import permutations, combinations
from collections import defaultdict, deque
import time

def tournament_from_tiling(n, tiling_bits):
    A = [[0]*n for _ in range(n)]
    for i in range(n-1):
        A[i][i+1] = 1
    idx = 0
    for i in range(n):
        for j in range(i+2, n):
            if (tiling_bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def num_tiling_bits(n):
    return n*(n-1)//2 - (n-1)

def count_3cycles(A, n):
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]:
            t3 += 1
        if A[i][k] and A[k][j] and A[j][i]:
            t3 += 1
    return t3

def canonical(A, n):
    best = None
    for perm in permutations(range(n)):
        form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i != j)
        if best is None or form < best:
            best = form
    return best

def score_seq(A, n):
    return tuple(sorted(sum(A[i]) for i in range(n)))

def converse(A, n):
    return [[A[j][i] for j in range(n)] for i in range(n)]

def is_SC(A, n):
    return canonical(A, n) == canonical(converse(A, n), n)

def tiling_transpose_pairs(n):
    edges = []
    for i in range(n):
        for j in range(i+2, n):
            edges.append((i, j))
    edge_to_idx = {e: idx for idx, e in enumerate(edges)}
    pairs = []
    fixed = []
    seen = set()
    for idx, (i, j) in enumerate(edges):
        if idx in seen:
            continue
        ti, tj = n-1-j, n-1-i
        if ti > tj:
            ti, tj = tj, ti
        if (ti, tj) == (i, j):
            fixed.append(idx)
            seen.add(idx)
        elif (ti, tj) in edge_to_idx:
            tidx = edge_to_idx[(ti, tj)]
            pairs.append((idx, tidx))
            seen.add(idx)
            seen.add(tidx)
    return pairs, fixed

def gen_gs_tilings(n, pairs, fixed):
    gs_dof = len(pairs) + len(fixed)
    result = []
    for free_val in range(2**gs_dof):
        bits = 0
        for k, (idx1, idx2) in enumerate(pairs):
            if (free_val >> k) & 1:
                bits |= (1 << idx1) | (1 << idx2)
        for k, fidx in enumerate(fixed):
            if (free_val >> (len(pairs) + k)) & 1:
                bits |= (1 << fidx)
        result.append(bits)
    return result

def flip_tiling(tiling_bits, m):
    return tiling_bits ^ ((1 << m) - 1)

n = 7
m = num_tiling_bits(n)
print(f"BIPARTITION INVARIANT VERIFICATION at n={n}")
print(f"  Tiling bits: {m}")

pairs, fixed = tiling_transpose_pairs(n)
gs_tilings = gen_gs_tilings(n, pairs, fixed)
print(f"  GS tilings: {len(gs_tilings)}")

# Build class database from ALL tilings
t0 = time.time()
print(f"  Building class database...", flush=True)

canon_db = {}
class_list = []
bits_to_class = {}

for bits in range(2**m):
    A = tournament_from_tiling(n, bits)
    c = canonical(A, n)
    if c not in canon_db:
        canon_db[c] = len(class_list)
        class_list.append({
            'canon': c, 'rep': A, 'sc': is_SC(A, n),
            'scores': score_seq(A, n),
            't3': count_3cycles(A, n),
            'gs_tilings': set()
        })
    bits_to_class[bits] = canon_db[c]

t1 = time.time()
print(f"  Built in {t1-t0:.1f}s, {len(class_list)} classes")

# Mark GS tilings
for bits in gs_tilings:
    idx = bits_to_class[bits]
    class_list[idx]['gs_tilings'].add(bits)

# Build GS flip graph
adj = defaultdict(set)
for bits in gs_tilings:
    c_from = bits_to_class[bits]
    flipped = flip_tiling(bits, m)
    c_to = bits_to_class[flipped]
    if c_from != c_to:
        adj[c_from].add(c_to)
        adj[c_to].add(c_from)

sc_in_skel = set()
for v in adj:
    sc_in_skel.add(v)
    for u in adj[v]:
        sc_in_skel.add(u)

# BFS bipartition
color = {}
is_bip = True
for start in sc_in_skel:
    if start in color:
        continue
    color[start] = 0
    queue = deque([start])
    while queue:
        u = queue.popleft()
        for v in adj[u]:
            if v not in color:
                color[v] = 1 - color[u]
                queue.append(v)
            elif color[v] == color[u]:
                is_bip = False

print(f"\n  BIPARTITE: {is_bip}")

if is_bip:
    side_A = sorted(v for v in sc_in_skel if color[v] == 0)
    side_B = sorted(v for v in sc_in_skel if color[v] == 1)

    t3_parity_A = set(class_list[v]['t3'] % 2 for v in side_A)
    t3_parity_B = set(class_list[v]['t3'] % 2 for v in side_B)

    print(f"  Side A ({len(side_A)} classes): t3 parities = {t3_parity_A}")
    print(f"  Side B ({len(side_B)} classes): t3 parities = {t3_parity_B}")

    if t3_parity_A == {0} and t3_parity_B == {1}:
        print(f"  ** CONFIRMED: t3 parity determines bipartition (A=even, B=odd) **")
    elif t3_parity_A == {1} and t3_parity_B == {0}:
        print(f"  ** CONFIRMED: t3 parity determines bipartition (A=odd, B=even) **")
    elif len(t3_parity_A) == 1 and len(t3_parity_B) == 1 and t3_parity_A != t3_parity_B:
        print(f"  ** CONFIRMED: t3 parity determines bipartition **")
    else:
        print(f"  ** t3 parity does NOT determine bipartition at n={n} **")
        # Check individual mismatches
        mismatches = 0
        for v in side_A:
            if class_list[v]['t3'] % 2 != list(t3_parity_A)[0] if len(t3_parity_A)==1 else -1:
                mismatches += 1
        print(f"  Side A t3 values: {sorted(class_list[v]['t3'] for v in side_A)}")
        print(f"  Side B t3 values: {sorted(class_list[v]['t3'] for v in side_B)}")

    # Also check: NSC class t3 parities
    nsc_indices = [i for i, c in enumerate(class_list) if not c['sc']]
    t3_nsc = [class_list[i]['t3'] for i in nsc_indices]
    t3_nsc_parities = set(t % 2 for t in t3_nsc)
    print(f"\n  NSC classes: {len(nsc_indices)}, t3 parities = {t3_nsc_parities}")

    # For each NSC pair: do they have the same or different t3 parity?
    seen = set()
    same_parity = 0
    diff_parity = 0
    for i in nsc_indices:
        A_op = converse(class_list[i]['rep'], n)
        partner = canon_db[canonical(A_op, n)]
        pair = (min(i, partner), max(i, partner))
        if pair in seen:
            continue
        seen.add(pair)
        if class_list[i]['t3'] % 2 == class_list[partner]['t3'] % 2:
            same_parity += 1
        else:
            diff_parity += 1

    print(f"\n  NSC pair t3 parity:")
    print(f"    Same parity (both on same side): {same_parity}")
    print(f"    Different parity (on opposite sides): {diff_parity}")

t2 = time.time()
print(f"\n  Total: {t2-t0:.1f}s")
print("DONE")
