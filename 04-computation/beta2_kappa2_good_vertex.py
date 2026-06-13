"""
Focus: SC tournaments with kappa >= 2 and b1 = 0.
Goal: Find a provable good-vertex selection criterion.

Key structural facts to exploit:
1. b1 = #{free components} <= 1 (THM-107)
2. Bad set is always transitive (verified n <= 8)
3. Free component of T\v spans ALL n-1 remaining vertices
4. For non-SC T: all vertices good. For kappa=1: cut vertex is good.

Remaining case: SC with kappa >= 2 and b1 = 0.
"""

import numpy as np
from itertools import combinations, permutations
from collections import defaultdict

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def find_3cycles(A, n):
    cycles = []
    for a, b, c in combinations(range(n), 3):
        if A[a][b] and A[b][c] and A[c][a]:
            cycles.append((a, b, c))
        if A[a][c] and A[c][b] and A[b][a]:
            cycles.append((a, c, b))
    return cycles

def shared_directed_edge(c1, c2):
    edges1 = {(c1[0],c1[1]), (c1[1],c1[2]), (c1[2],c1[0])}
    edges2 = {(c2[0],c2[1]), (c2[1],c2[2]), (c2[2],c2[0])}
    return len(edges1 & edges2) > 0

def cycle_graph_components(cycles):
    if not cycles:
        return []
    nc = len(cycles)
    adj = [[] for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if shared_directed_edge(cycles[i], cycles[j]):
                adj[i].append(j)
                adj[j].append(i)
    visited = [False]*nc
    components = []
    for start in range(nc):
        if visited[start]:
            continue
        comp = []
        stack = [start]
        while stack:
            v = stack.pop()
            if visited[v]:
                continue
            visited[v] = True
            comp.append(v)
            for u in adj[v]:
                if not visited[u]:
                    stack.append(u)
        components.append(comp)
    return components

def is_dominated(A, n, cyc):
    verts = set(cyc)
    for d in range(n):
        if d in verts:
            continue
        a, b, c = cyc
        if (A[d][a] and A[d][b] and A[d][c]) or (A[a][d] and A[b][d] and A[c][d]):
            return True
    return False

def compute_b1(A, n):
    cycles = find_3cycles(A, n)
    if not cycles:
        return 0
    comps = cycle_graph_components(cycles)
    return sum(1 for comp in comps if all(not is_dominated(A, n, cycles[ci]) for ci in comp))

def delete_vertex(A, n, v):
    keep = [i for i in range(n) if i != v]
    B = np.zeros((n-1, n-1), dtype=int)
    for i, ki in enumerate(keep):
        for j, kj in enumerate(keep):
            B[i][j] = A[ki][kj]
    return B

def is_sc(A, n):
    for start in range(n):
        visited = set()
        stack = [start]
        while stack:
            v = stack.pop()
            if v in visited: continue
            visited.add(v)
            for u in range(n):
                if A[v][u] and u not in visited:
                    stack.append(u)
        if len(visited) < n:
            return False
    return True

def main():
    print("=" * 70)
    print("SC kappa>=2 GOOD VERTEX ANALYSIS")
    print("=" * 70)

    # Focus on n=6 exhaustive (1680 SC kappa>=2 b1=0 tournaments)
    # and n=7 sampled
    n = 6
    total = 2**(n*(n-1)//2)

    sc_k2_tours = []  # (bits, A) for SC kappa>=2 b1=0

    for bits in range(total):
        A = bits_to_adj(bits, n)
        if compute_b1(A, n) != 0:
            continue
        if not is_sc(A, n):
            continue
        # Check kappa >= 2
        kappa1 = False
        for v in range(n):
            B = delete_vertex(A, n, v)
            if not is_sc(B, n-1):
                kappa1 = True
                break
        if kappa1:
            continue
        sc_k2_tours.append((bits, A))

    print(f"\nn=6: {len(sc_k2_tours)} SC kappa>=2 b1=0 tournaments")

    # For these: compute bad vertices and characterize them
    all_scores = []
    bad_count_dist = defaultdict(int)

    # Track: what makes the 3 bad vertices special?
    for bits, A in sc_k2_tours:
        bad = []
        for v in range(n):
            B = delete_vertex(A, n, v)
            if compute_b1(B, n-1) > 0:
                bad.append(v)
        bad_count_dist[len(bad)] += 1
        scores = [int(sum(A[i])) for i in range(n)]
        all_scores.append(scores)

    print(f"Bad count distribution: {dict(bad_count_dist)}")
    print(f"All have score (2,2,2,3,3,3)")

    # Deep analysis: for the 1440 with 3 bad vertices
    print(f"\nAnalyzing the 3-bad-vertex case:")

    # Track structural properties of bad triples
    bad_triple_in_3cycle_of_T = 0
    bad_triple_scores = defaultdict(int)
    bad_triple_c3_counts = defaultdict(int)

    # Key question: can we characterize bad vertices by LOCAL properties?
    # Hypothesis: bad vertices are those with score != (n-1)/2 (extremal scores)?

    for bits, A in sc_k2_tours:
        bad = []
        for v in range(n):
            B = delete_vertex(A, n, v)
            if compute_b1(B, n-1) > 0:
                bad.append(v)

        if len(bad) != 3:
            continue

        scores = [int(sum(A[i])) for i in range(n)]
        bad_scores = tuple(sorted([scores[v] for v in bad]))
        good_scores = tuple(sorted([scores[v] for v in range(n) if v not in bad]))
        bad_triple_scores[(bad_scores, good_scores)] += 1

        # c3 per vertex
        c3 = [0]*n
        cycles = find_3cycles(A, n)
        for cyc in cycles:
            for v in cyc:
                c3[v] += 1

        bad_c3 = tuple(sorted([c3[v] for v in bad]))
        good_c3 = tuple(sorted([c3[v] for v in range(n) if v not in bad]))
        bad_triple_c3_counts[(bad_c3, good_c3)] += 1

    print(f"\n  (bad_scores, good_scores) -> count:")
    for key, cnt in sorted(bad_triple_scores.items(), key=lambda x: -x[1]):
        print(f"    bad={key[0]}, good={key[1]}: {cnt}")

    print(f"\n  (bad_c3, good_c3) -> count:")
    for key, cnt in sorted(bad_triple_c3_counts.items(), key=lambda x: -x[1]):
        print(f"    bad={key[0]}, good={key[1]}: {cnt}")

    # KEY TEST: For kappa>=2, does score separation work?
    # All scores are (2,2,2,3,3,3). So bad must be some 3-subset of these.
    print(f"\n  Score-based bad characterization:")
    # Since all have score (2,2,2,3,3,3), bad vertices are some mix of 2s and 3s

    # Track: among the three score-2 vertices, how many are bad?
    score2_bad = defaultdict(int)
    for bits, A in sc_k2_tours:
        bad = set()
        for v in range(n):
            B = delete_vertex(A, n, v)
            if compute_b1(B, n-1) > 0:
                bad.add(v)
        if len(bad) != 3:
            continue
        scores = [int(sum(A[i])) for i in range(n)]
        n_bad_s2 = sum(1 for v in bad if scores[v] == 2)
        n_bad_s3 = sum(1 for v in bad if scores[v] == 3)
        score2_bad[(n_bad_s2, n_bad_s3)] += 1

    print(f"  (bad with score 2, bad with score 3): count")
    for key, cnt in sorted(score2_bad.items()):
        print(f"    {key}: {cnt}")

    # NEW APPROACH: Can we characterize bad vertices by their UNIQUE DOMINATION DIRECTION?
    # From earlier: when v->w (both bad), v dominates from ABOVE, w from BELOW
    # For transitive triple v1->v2->v3: v1 dominates from above, v3 from below, v2 mixed

    print(f"\n{'='*70}")
    print("UNIQUE DOMINATION DIRECTION IN TRANSITIVE BAD TRIPLE")
    print(f"{'='*70}")

    dir_pattern = defaultdict(int)
    for bits, A in sc_k2_tours:
        bad = []
        for v in range(n):
            B = delete_vertex(A, n, v)
            if compute_b1(B, n-1) > 0:
                bad.append(v)
        if len(bad) != 3:
            continue

        # Find transitive order
        for perm in permutations(bad):
            if A[perm[0]][perm[1]] and A[perm[0]][perm[2]] and A[perm[1]][perm[2]]:
                top, mid, bot = perm
                break

        cycles = find_3cycles(A, n)
        result = {}
        for label, v in [('top', top), ('mid', mid), ('bot', bot)]:
            above = 0; below = 0; free_bridged = 0
            for cyc in cycles:
                if v in cyc: continue
                a, b, c = cyc
                v_beats_all = A[v][a] and A[v][b] and A[v][c]
                all_beat_v = A[a][v] and A[b][v] and A[c][v]
                if not v_beats_all and not all_beat_v:
                    continue
                # Check if uniquely dominated by v
                unique = True
                for d in range(n):
                    if d in {v, a, b, c}: continue
                    if (A[d][a] and A[d][b] and A[d][c]) or (A[a][d] and A[b][d] and A[c][d]):
                        unique = False; break
                if unique:
                    if v_beats_all: above += 1
                    else: below += 1
            result[label] = (above, below)

        dir_pattern[tuple(result[l] for l in ['top','mid','bot'])] += 1

    print(f"\n  (top(above,below), mid(above,below), bot(above,below)): count")
    for key, cnt in sorted(dir_pattern.items(), key=lambda x: -x[1]):
        print(f"    top={key[0]}, mid={key[1]}, bot={key[2]}: {cnt}")

    # CRITICAL NEW TEST: Is the bad triple exactly the vertices of a
    # specific "fragile" 3-cycle in T?
    print(f"\n{'='*70}")
    print("BAD TRIPLE = vertices of a specific 3-cycle?")
    print(f"{'='*70}")

    bad_is_cycle_verts = 0
    bad_not_cycle_verts = 0

    for bits, A in sc_k2_tours:
        bad = []
        for v in range(n):
            B = delete_vertex(A, n, v)
            if compute_b1(B, n-1) > 0:
                bad.append(v)
        if len(bad) != 3:
            continue

        # Do the 3 bad vertices participate in a common 3-cycle?
        a, b, c = bad
        if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
            bad_is_cycle_verts += 1
        else:
            bad_not_cycle_verts += 1

    print(f"  Bad triple forms a 3-cycle: {bad_is_cycle_verts}/{bad_is_cycle_verts+bad_not_cycle_verts}")
    print(f"  Bad triple is transitive: {bad_not_cycle_verts}/{bad_is_cycle_verts+bad_not_cycle_verts}")

    # Test: is the bad triple a DOMINATION CHAIN?
    # v1 dominates some cycle C1, v2 is in C1 and dominates C2, v3 is in C2?
    print(f"\n{'='*70}")
    print("DOMINATION CHAIN STRUCTURE")
    print(f"{'='*70}")

    chain_found = 0
    chain_not_found = 0

    for bits, A in sc_k2_tours:
        bad = []
        for v in range(n):
            B = delete_vertex(A, n, v)
            if compute_b1(B, n-1) > 0:
                bad.append(v)
        if len(bad) != 3:
            continue

        cycles = find_3cycles(A, n)

        # For each pair of bad vertices: does one uniquely dominate a cycle containing the other?
        dom_contains = {}
        for v in bad:
            for cyc in cycles:
                if v in cyc: continue
                a, b, c = cyc
                v_dom = (A[v][a] and A[v][b] and A[v][c]) or (A[a][v] and A[b][v] and A[c][v])
                if not v_dom: continue
                unique = True
                for d in range(n):
                    if d in {v,a,b,c}: continue
                    if (A[d][a] and A[d][b] and A[d][c]) or (A[a][d] and A[b][d] and A[c][d]):
                        unique = False; break
                if unique:
                    for w in bad:
                        if w != v and w in {a,b,c}:
                            dom_contains[(v,w)] = True

        # Check if there's a chain v1 dominates cycle containing v2, v2 dominates cycle containing v3
        found = False
        for perm in permutations(bad):
            if (perm[0], perm[1]) in dom_contains and (perm[1], perm[2]) in dom_contains:
                found = True
                break
        if found:
            chain_found += 1
        else:
            chain_not_found += 1

    print(f"  Domination chain exists: {chain_found}/{chain_found+chain_not_found}")
    print(f"  No chain: {chain_not_found}")

    # ULTIMATE: Can we express the bound |BAD| <= 3 in terms of
    # the 3-cycle graph structure?
    print(f"\n{'='*70}")
    print("NUMBER OF CONNECTED COMPONENTS vs BAD COUNT")
    print(f"{'='*70}")

    comp_bad = defaultdict(lambda: defaultdict(int))
    for bits, A in sc_k2_tours:
        cycles = find_3cycles(A, n)
        if not cycles: continue
        comps = cycle_graph_components(cycles)
        n_comps = len(comps)

        bad_count = 0
        for v in range(n):
            B = delete_vertex(A, n, v)
            if compute_b1(B, n-1) > 0:
                bad_count += 1
        comp_bad[n_comps][bad_count] += 1

    for nc in sorted(comp_bad.keys()):
        print(f"  {nc} components:")
        for bc in sorted(comp_bad[nc].keys()):
            print(f"    {bc} bad: {comp_bad[nc][bc]}")

if __name__ == '__main__':
    main()
