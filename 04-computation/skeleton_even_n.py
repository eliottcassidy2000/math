#!/usr/bin/env python3
"""
Even-n skeleton: NOT bipartite. What's the topology?

At even n, t3 parity is preserved by GS flip (since n-2 is even).
So same-class flips can occur, and the skeleton may have odd cycles.

Explore:
1. The odd cycle structure at n=4,6
2. What invariant DOES change under flip at even n?
3. Connection to Mobius strip / non-orientable topology

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

for n in [4, 6]:
    m = num_tiling_bits(n)
    print(f"\n{'='*60}")
    print(f"EVEN-n SKELETON at n={n}")
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
                'gs_tilings': set()
            })
        bits_to_class[bits] = canon_db[c]

    for bits in gs_tilings:
        idx = bits_to_class[bits]
        class_list[idx]['gs_tilings'].add(bits)

    sc_indices = [i for i, c in enumerate(class_list) if c['sc']]
    nsc_indices = [i for i, c in enumerate(class_list) if not c['sc']]
    print(f"  {len(class_list)} classes, {len(sc_indices)} SC, {len(nsc_indices)} NSC")
    print(f"  {len(gs_tilings)} GS tilings")

    # Build GS flip graph
    adj = defaultdict(set)
    gs_edges = defaultdict(int)
    same_class_count = 0

    for bits in gs_tilings:
        c_from = bits_to_class[bits]
        flipped = flip_tiling(bits, m)
        c_to = bits_to_class[flipped]
        if c_from == c_to:
            same_class_count += 1
        else:
            adj[c_from].add(c_to)
            adj[c_to].add(c_from)
            edge = (min(c_from, c_to), max(c_from, c_to))
            gs_edges[edge] += 1

    sc_in_skel = set()
    for v in adj:
        sc_in_skel.add(v)
        for u in adj[v]:
            sc_in_skel.add(u)

    print(f"  Same-class GS flips: {same_class_count}")
    print(f"  Skeleton: {len(sc_in_skel)} vertices, {len(gs_edges)} edges")

    # Edges
    print(f"  Edges:")
    for (i, j), w in sorted(gs_edges.items()):
        print(f"    {i}(t3={class_list[i]['t3']}) <-[{w}]-> {j}(t3={class_list[j]['t3']})")

    # Self-loops (same-class flips)
    for i in sc_indices:
        self_count = 0
        for bits in class_list[i]['gs_tilings']:
            flipped = flip_tiling(bits, m)
            if bits_to_class[flipped] == i:
                self_count += 1
        if self_count > 0:
            print(f"    {i}(t3={class_list[i]['t3']}) self-flip: {self_count} tilings")

    # Check bipartiteness
    color = {}
    is_bip = True
    odd_edge = None
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
                    odd_edge = (u, v)

    print(f"\n  Bipartite: {is_bip}")
    if not is_bip and odd_edge:
        u, v = odd_edge
        print(f"  Odd edge: {u}(t3={class_list[u]['t3']}) -- {v}(t3={class_list[v]['t3']})")
        print(f"  Same t3 parity? {class_list[u]['t3'] % 2 == class_list[v]['t3'] % 2}")

    # Find shortest odd cycle via BFS
    if not is_bip:
        min_odd_cycle = float('inf')
        for start in sc_in_skel:
            dist = {start: 0}
            parent = {start: None}
            queue = deque([start])
            while queue:
                u = queue.popleft()
                for v in adj[u]:
                    if v not in dist:
                        dist[v] = dist[u] + 1
                        parent[v] = u
                        queue.append(v)
                    elif dist[v] >= dist[u]:
                        cycle_len = dist[u] + dist[v] + 1
                        if cycle_len % 2 == 1 and cycle_len < min_odd_cycle:
                            min_odd_cycle = cycle_len
        if min_odd_cycle < float('inf'):
            print(f"  Shortest odd cycle: length {min_odd_cycle}")

    # t3 analysis: at even n, what DO the edges look like?
    t3_diffs = []
    for (i, j) in gs_edges:
        t3_diffs.append(class_list[i]['t3'] - class_list[j]['t3'])
    print(f"\n  t3 diff across edges: {sorted(set(t3_diffs))}")

    # Check: at even n, does flip preserve t3 or just its parity?
    for bits in gs_tilings[:5]:
        A = tournament_from_tiling(n, bits)
        Af = tournament_from_tiling(n, flip_tiling(bits, m))
        t3_b = count_3cycles(A, n)
        t3_a = count_3cycles(Af, n)
        print(f"  GS tiling bits={bits}: t3={t3_b} -> flip t3={t3_a}, sum={t3_b+t3_a}, diff={t3_a-t3_b}")
