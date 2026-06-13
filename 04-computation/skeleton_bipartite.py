#!/usr/bin/env python3
"""
Is the blue line skeleton BIPARTITE?

If so, it has two natural "sides" and NSC classes may connect
preferentially to one side — this would explain the user's observation
that "most paired iso classes connect via black line to only one side."

Also check: what invariant determines the bipartition?

kind-pasteur-2026-03-06-S25h
"""
from itertools import permutations
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

def is_bipartite(adj, vertices):
    """Check bipartiteness and return partition. BFS."""
    color = {}
    for start in vertices:
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
                    return False, None
    return True, color

def analyze(n):
    m = num_tiling_bits(n)
    print(f"\n{'='*60}")
    print(f"BIPARTITE ANALYSIS at n={n}")
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
                'scores': score_seq(A, n), 'gs_tilings': set()
            })
        idx = canon_db[c]
        class_list[idx]['tilings'].add(bits)
        bits_to_class[bits] = idx

    for bits in gs_tilings:
        idx = bits_to_class[bits]
        class_list[idx]['gs_tilings'].add(bits)

    num_classes = len(class_list)
    sc_indices = [i for i, c in enumerate(class_list) if c['sc']]
    nsc_indices = [i for i, c in enumerate(class_list) if not c['sc']]

    # Build GS flip graph
    adj = defaultdict(set)
    gs_edges = defaultdict(int)
    for bits in gs_tilings:
        c_from = bits_to_class[bits]
        flipped = flip_tiling(bits, m)
        c_to = bits_to_class[flipped]
        if c_from != c_to:
            adj[c_from].add(c_to)
            adj[c_to].add(c_from)
            edge = (min(c_from, c_to), max(c_from, c_to))
            gs_edges[edge] += 1

    sc_in_skel = set()
    for (i, j) in gs_edges:
        sc_in_skel.add(i)
        sc_in_skel.add(j)

    bip, coloring = is_bipartite(adj, sc_in_skel)
    print(f"  Skeleton: {len(sc_in_skel)} vertices, {len(gs_edges)} edges")
    print(f"  BIPARTITE: {bip}")

    if bip and coloring:
        side_A = [v for v in sc_in_skel if coloring[v] == 0]
        side_B = [v for v in sc_in_skel if coloring[v] == 1]
        print(f"  Side A: {len(side_A)} classes, Side B: {len(side_B)} classes")

        # What distinguishes side A from side B?
        # Check: HW parity of GS tilings
        print(f"\n  Side A classes:")
        for v in sorted(side_A)[:10]:
            gs_hws = [bin(b).count('1') for b in class_list[v]['gs_tilings']]
            print(f"    Class {v}: scores={class_list[v]['scores']}, GS HWs={sorted(gs_hws)}")

        print(f"\n  Side B classes:")
        for v in sorted(side_B)[:10]:
            gs_hws = [bin(b).count('1') for b in class_list[v]['gs_tilings']]
            print(f"    Class {v}: scores={class_list[v]['scores']}, GS HWs={sorted(gs_hws)}")

        # Check HW parity pattern
        side_A_parities = set()
        side_B_parities = set()
        for v in side_A:
            for b in class_list[v]['gs_tilings']:
                side_A_parities.add(bin(b).count('1') % 2)
        for v in side_B:
            for b in class_list[v]['gs_tilings']:
                side_B_parities.add(bin(b).count('1') % 2)
        print(f"\n  Side A HW parities: {side_A_parities}")
        print(f"  Side B HW parities: {side_B_parities}")

        # NSC sidedness: for each NSC pair, which side do they connect to?
        print(f"\n  NSC SIDE PREFERENCE:")
        seen_pairs = set()
        one_side_count = 0
        both_sides_count = 0
        no_sc_count = 0

        for i in nsc_indices:
            A_op = converse(class_list[i]['rep'], n)
            partner = canon_db[canonical(A_op, n)]

            pair_key = (min(i, partner), max(i, partner))
            if pair_key in seen_pairs:
                continue
            seen_pairs.add(pair_key)

            # For this NSC class, find nearest SC class by min Hamming distance
            min_dist_A = float('inf')
            min_dist_B = float('inf')
            nearest_A = None
            nearest_B = None

            for bits_i in class_list[i]['tilings']:
                for v in side_A:
                    for bits_v in class_list[v]['tilings']:
                        d = bin(bits_i ^ bits_v).count('1')
                        if d < min_dist_A:
                            min_dist_A = d
                            nearest_A = v
                for v in side_B:
                    for bits_v in class_list[v]['tilings']:
                        d = bin(bits_i ^ bits_v).count('1')
                        if d < min_dist_B:
                            min_dist_B = d
                            nearest_B = v

            # Same for partner
            min_dist_A_p = float('inf')
            min_dist_B_p = float('inf')

            for bits_p in class_list[partner]['tilings']:
                for v in side_A:
                    for bits_v in class_list[v]['tilings']:
                        d = bin(bits_p ^ bits_v).count('1')
                        if d < min_dist_A_p:
                            min_dist_A_p = d
                for v in side_B:
                    for bits_v in class_list[v]['tilings']:
                        d = bin(bits_p ^ bits_v).count('1')
                        if d < min_dist_B_p:
                            min_dist_B_p = d

            # Which side is each member of the pair closer to?
            i_side = "A" if min_dist_A < min_dist_B else ("B" if min_dist_B < min_dist_A else "=")
            p_side = "A" if min_dist_A_p < min_dist_B_p else ("B" if min_dist_B_p < min_dist_A_p else "=")

            if i_side != "=" and p_side != "=" and i_side != p_side:
                one_side_count += 1
                tag = "OPPOSITE-SIDED"
            elif i_side == p_side and i_side != "=":
                both_sides_count += 1
                tag = "same-sided"
            else:
                no_sc_count += 1
                tag = "equidistant"

            print(f"    Pair ({i},{partner}) scores={class_list[i]['scores']}: "
                  f"class {i}->side {i_side} (dA={min_dist_A},dB={min_dist_B}), "
                  f"class {partner}->side {p_side} (dA={min_dist_A_p},dB={min_dist_B_p}) "
                  f"-> {tag}")

        print(f"\n  SIDE PREFERENCE SUMMARY:")
        print(f"    Opposite-sided (pair members on different sides): {one_side_count}")
        print(f"    Same-sided (pair members on same side): {both_sides_count}")
        print(f"    Equidistant: {no_sc_count}")

analyze(5)
print("\nDONE")
