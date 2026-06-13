#!/usr/bin/env python3
"""
NSC class t3 parity analysis.

Questions:
1. Do NSC pairs (class i, partner j) have the same or opposite t3 parity?
2. If opposite: each partner connects to a different "side" of the skeleton
3. What about at EVEN n?

Also: the skeleton at even n — is it still bipartite?

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
        if A[i][j] and A[j][k] and A[k][i]: t3 += 1
        if A[i][k] and A[k][j] and A[j][i]: t3 += 1
    return t3

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

for n in [4, 5, 6]:
    m = num_tiling_bits(n)
    print(f"\n{'='*60}")
    print(f"n={n}")
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
                'rep': A, 'sc': is_SC(A, n),
                'scores': score_seq(A, n),
                't3': count_3cycles(A, n),
                'gs_tilings': set(), 'tilings': set()
            })
        idx = canon_db[c]
        class_list[idx]['tilings'].add(bits)
        bits_to_class[bits] = idx

    for bits in gs_tilings:
        idx = bits_to_class[bits]
        class_list[idx]['gs_tilings'].add(bits)

    sc_indices = [i for i, c in enumerate(class_list) if c['sc']]
    nsc_indices = [i for i, c in enumerate(class_list) if not c['sc']]
    print(f"  {len(class_list)} classes, {len(sc_indices)} SC, {len(nsc_indices)} NSC")

    # Build GS flip graph
    adj = defaultdict(set)
    same_class = 0
    for bits in gs_tilings:
        c_from = bits_to_class[bits]
        flipped = flip_tiling(bits, m)
        c_to = bits_to_class[flipped]
        if c_from != c_to:
            adj[c_from].add(c_to)
            adj[c_to].add(c_from)
        else:
            same_class += 1

    print(f"  Same-class GS flips: {same_class}")

    sc_in_skel = set()
    for v in adj:
        sc_in_skel.add(v)
        for u in adj[v]:
            sc_in_skel.add(u)

    # Bipartiteness
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

    print(f"  Skeleton bipartite: {is_bip}")

    if is_bip and sc_in_skel:
        side_A = [v for v in sc_in_skel if color[v] == 0]
        side_B = [v for v in sc_in_skel if color[v] == 1]
        t3_parity_A = set(class_list[v]['t3'] % 2 for v in side_A)
        t3_parity_B = set(class_list[v]['t3'] % 2 for v in side_B)
        print(f"  Side A ({len(side_A)}): t3 parities = {t3_parity_A}")
        print(f"  Side B ({len(side_B)}): t3 parities = {t3_parity_B}")

    # NSC pair t3 parity
    print(f"\n  NSC PAIR t3 PARITY:")
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

        t3_i = class_list[i]['t3']
        t3_p = class_list[partner]['t3']

        if t3_i % 2 == t3_p % 2:
            same_parity += 1
        else:
            diff_parity += 1

    print(f"    Same t3 parity: {same_parity}")
    print(f"    Different t3 parity: {diff_parity}")

    # For NSC classes with t3 values:
    # The converse T^op has the SAME t3 as T (since reversing all edges
    # maps CW 3-cycles to CCW and vice versa, preserving count)
    # So NSC pair members ALWAYS have the same t3!
    print(f"\n  NOTE: T and T^op always have same t3 (reversing all edges")
    print(f"  maps CW<->CCW), so NSC pairs MUST have same t3 parity.")

    if is_bip and sc_in_skel:
        # Therefore both members of an NSC pair connect to the SAME side!
        print(f"  => NSC pair members are on the SAME side of the bipartition")
        print(f"  => They connect to the skeleton from the same side")

        # Which side has more NSC classes?
        nsc_on_A = sum(1 for i in nsc_indices if class_list[i]['t3'] % 2 == list(t3_parity_A)[0])
        nsc_on_B = len(nsc_indices) - nsc_on_A
        print(f"  NSC on side A (t3 parity {t3_parity_A}): {nsc_on_A}")
        print(f"  NSC on side B (t3 parity {t3_parity_B}): {nsc_on_B}")
