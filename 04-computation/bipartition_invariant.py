#!/usr/bin/env python3
"""
What invariant determines the bipartition of the blue line skeleton?

The skeleton is bipartite at n=5 (4+4) and n=7 (44+44).
Check: is it the parity of t3? Or some other cycle invariant?

Also check n=3 (trivially bipartite, 1 edge).

kind-pasteur-2026-03-06-S25h
"""
from itertools import permutations, combinations
from collections import defaultdict, deque

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

def count_3cycles(A, n):
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]:
            t3 += 1
        if A[i][k] and A[k][j] and A[j][i]:
            t3 += 1
    return t3

def count_5cycles(A, n):
    """Count directed 5-cycles."""
    t5 = 0
    for verts in combinations(range(n), 5):
        for perm in permutations(verts):
            if all(A[perm[i]][perm[(i+1) % 5]] for i in range(5)):
                t5 += 1
    return t5 // 5  # Each 5-cycle counted 5 times (cyclic rotations)

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

def analyze(n):
    m = num_tiling_bits(n)
    print(f"\n{'='*60}")
    print(f"BIPARTITION INVARIANT at n={n}")
    print(f"{'='*60}")

    pairs, fixed = tiling_transpose_pairs(n)
    gs_tilings = gen_gs_tilings(n, pairs, fixed)

    # Build class database
    canon_db = {}
    class_list = []
    bits_to_class = {}

    for bits in range(2**m):
        A = tournament_from_tiling(n, bits)
        c = canonical(A, n)
        if c not in canon_db:
            canon_db[c] = len(class_list)
            class_list.append({
                'canon': c, 'rep': A, 'tilings': set(), 'sc': is_SC(A, n),
                'scores': score_seq(A, n), 'gs_tilings': set(),
                't3': count_3cycles(A, n)
            })
        idx = canon_db[c]
        class_list[idx]['tilings'].add(bits)
        bits_to_class[bits] = idx

    for bits in gs_tilings:
        idx = bits_to_class[bits]
        class_list[idx]['gs_tilings'].add(bits)

    sc_indices = [i for i, c in enumerate(class_list) if c['sc']]

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

    # Check bipartiteness
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

    print(f"  BIPARTITE: {is_bip}")

    if not is_bip:
        print("  Not bipartite, skipping invariant check")
        return

    side_A = sorted(v for v in sc_in_skel if color[v] == 0)
    side_B = sorted(v for v in sc_in_skel if color[v] == 1)

    # Check t3 parity
    t3_A = [class_list[v]['t3'] for v in side_A]
    t3_B = [class_list[v]['t3'] for v in side_B]
    t3_parity_A = set(t % 2 for t in t3_A)
    t3_parity_B = set(t % 2 for t in t3_B)

    print(f"\n  Side A ({len(side_A)} classes):")
    print(f"    t3 values: {sorted(t3_A)}")
    print(f"    t3 parities: {t3_parity_A}")

    print(f"\n  Side B ({len(side_B)} classes):")
    print(f"    t3 values: {sorted(t3_B)}")
    print(f"    t3 parities: {t3_parity_B}")

    if t3_parity_A == {0} and t3_parity_B == {1}:
        print(f"\n  ** t3 PARITY DETERMINES BIPARTITION: A=even, B=odd **")
    elif t3_parity_A == {1} and t3_parity_B == {0}:
        print(f"\n  ** t3 PARITY DETERMINES BIPARTITION: A=odd, B=even **")
    else:
        print(f"\n  t3 parity does NOT determine bipartition")

    # Check other potential invariants
    # Total score sum modulo something?
    # Score sum is always n*(n-1)/2 for tournaments, so that's constant.

    # Check: Hamming weight of representative GS tiling?
    print(f"\n  Representative GS tiling HW:")
    for side_name, side in [("A", side_A), ("B", side_B)]:
        hws = []
        for v in side:
            for b in class_list[v]['gs_tilings']:
                hws.append(bin(b).count('1'))
        hw_parities = set(h % 2 for h in hws)
        print(f"    Side {side_name}: HW parities = {hw_parities}")

    # Check: sum of squared scores mod 2?
    print(f"\n  Sum of squared scores:")
    for side_name, side in [("A", side_A), ("B", side_B)]:
        ssq = [sum(s**2 for s in class_list[v]['scores']) for v in side]
        ssq_mod2 = set(s % 2 for s in ssq)
        print(f"    Side {side_name}: sum(scores^2) mod 2 = {ssq_mod2}")

    # Check: D = C(n,3) + 2*t3 = total dominance count
    # D mod 2 = C(n,3) mod 2 + 0 = C(n,3) mod 2 (constant)
    # So D doesn't distinguish.

    # Check t3 mod 4?
    t3_mod4_A = set(t % 4 for t in t3_A)
    t3_mod4_B = set(t % 4 for t in t3_B)
    print(f"\n  t3 mod 4:")
    print(f"    Side A: {sorted(t3_mod4_A)}")
    print(f"    Side B: {sorted(t3_mod4_B)}")

    # Print detailed class info
    print(f"\n  Detailed class info:")
    for side_name, side in [("A", side_A), ("B", side_B)]:
        for v in side:
            c = class_list[v]
            gs_hws = sorted(bin(b).count('1') for b in c['gs_tilings'])
            print(f"    Side {side_name}: class {v:2d}, scores={c['scores']}, t3={c['t3']}, GS_HWs={gs_hws}")

for n in [3, 5]:
    analyze(n)

print("\nDONE")
