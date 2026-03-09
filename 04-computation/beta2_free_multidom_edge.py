"""
Key lemma for proving good-vertex existence:

CLAIM: If there exist a free cycle F and a cycle D with |dom(D)| >= 2
such that F and D share a directed edge, then every vertex NOT in F is good.

Proof sketch:
- For v not in F: F is in freed(v).
- D remains dominated in T\v (|dom(D)| >= 2, so removing v leaves at least 1 dominator).
- F and D share a directed edge, so iso_edges(v) > 0.
- Since bad => iso_edges = 0, v is good.

This gives n-3 good vertices from a single (F, D) pair!

REMAINING CASE: Every cycle adjacent to a free cycle is either free or uniquely dominated.

Test: How often does the remaining case actually occur?
If it NEVER occurs for n >= some threshold, the lemma suffices.
"""
import numpy as np
from itertools import combinations
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

def random_tournament(n, rng):
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def find_3cycles(A, n):
    cycles = []
    for a, b, c in combinations(range(n), 3):
        if A[a][b] and A[b][c] and A[c][a]:
            cycles.append((a, b, c))
        if A[a][c] and A[c][b] and A[b][a]:
            cycles.append((a, c, b))
    return cycles

def get_dom(A, n, cyc):
    doms = []
    a, b, c = cyc
    for d in range(n):
        if d in {a,b,c}: continue
        if (A[d][a] and A[d][b] and A[d][c]) or (A[a][d] and A[b][d] and A[c][d]):
            doms.append(d)
    return doms

def shared_directed_edge(c1, c2):
    edges1 = {(c1[0],c1[1]), (c1[1],c1[2]), (c1[2],c1[0])}
    edges2 = {(c2[0],c2[1]), (c2[1],c2[2]), (c2[2],c2[0])}
    return len(edges1 & edges2) > 0

def compute_b1(A, n):
    cycles = find_3cycles(A, n)
    if not cycles:
        return 0
    nc = len(cycles)
    adj = [[] for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if shared_directed_edge(cycles[i], cycles[j]):
                adj[i].append(j)
                adj[j].append(i)
    visited = [False]*nc
    b1 = 0
    for start in range(nc):
        if visited[start]: continue
        comp = []
        stack = [start]
        while stack:
            v = stack.pop()
            if visited[v]: continue
            visited[v] = True
            comp.append(v)
            for u in adj[v]:
                if not visited[u]: stack.append(u)
        if all(not get_dom(A, n, cycles[ci]) for ci in comp):
            b1 += 1
    return b1

def delete_vertex(A, n, v):
    keep = [i for i in range(n) if i != v]
    B = np.zeros((n-1, n-1), dtype=int)
    for i, ki in enumerate(keep):
        for j, kj in enumerate(keep):
            B[i][j] = A[ki][kj]
    return B

def main():
    rng = np.random.RandomState(42)

    print("=" * 70)
    print("FREE-MULTIDOM EDGE LEMMA: Does the condition always hold?")
    print("=" * 70)

    for n in [5, 6, 7, 8, 9, 10]:
        if n <= 6:
            total = 2**(n*(n-1)//2)
            sample = range(total)
            method = "exhaustive"
        else:
            nsamp = {7: 10000, 8: 5000, 9: 2000, 10: 1000}[n]
            sample = range(nsamp)
            method = f"sampled({nsamp})"

        has_free_multidom_edge = 0
        no_free_multidom_edge = 0
        total_b1_0 = 0
        no_edge_but_good_exists = 0
        no_edge_and_all_bad = 0  # would be counterexample

        for idx in sample:
            if n <= 6:
                A = bits_to_adj(idx, n)
            else:
                A = random_tournament(n, rng)

            if compute_b1(A, n) != 0:
                continue
            total_b1_0 += 1

            cycles = find_3cycles(A, n)
            if not cycles:
                continue

            # Classify cycles
            free_cycs = []
            unique_dom_cycs = []  # dom size = 1
            multi_dom_cycs = []  # dom size >= 2

            for cyc in cycles:
                d = get_dom(A, n, cyc)
                if len(d) == 0:
                    free_cycs.append(cyc)
                elif len(d) == 1:
                    unique_dom_cycs.append(cyc)
                else:
                    multi_dom_cycs.append(cyc)

            # Check: does any free cycle share a directed edge with any multi-dom cycle?
            found = False
            for fc in free_cycs:
                for mc in multi_dom_cycs:
                    if shared_directed_edge(fc, mc):
                        found = True
                        break
                if found:
                    break

            if found:
                has_free_multidom_edge += 1
            else:
                no_free_multidom_edge += 1

                # In this case, does a good vertex still exist?
                has_good = False
                for v in range(n):
                    B = delete_vertex(A, n, v)
                    if compute_b1(B, n-1) == 0:
                        has_good = True
                        break
                if has_good:
                    no_edge_but_good_exists += 1
                else:
                    no_edge_and_all_bad += 1

        pct_edge = 100*has_free_multidom_edge/max(total_b1_0, 1)
        print(f"\nn={n} ({method}): {total_b1_0} with b1=0")
        print(f"  free-multidom edge EXISTS: {has_free_multidom_edge} ({pct_edge:.1f}%)")
        print(f"  free-multidom edge ABSENT: {no_free_multidom_edge}")
        if no_free_multidom_edge > 0:
            print(f"    of those, good vertex exists: {no_edge_but_good_exists}")
            print(f"    of those, ALL bad (counterex): {no_edge_and_all_bad}")

    # PART 2: In the "no free-multidom edge" case, what's the structure?
    # Free cycles only neighbor free or uniquely-dominated cycles.
    print("\n" + "=" * 70)
    print("PART 2: Structure when no free-multidom edge")
    print("-" * 50)

    for n in [6, 7]:
        if n <= 6:
            total = 2**(n*(n-1)//2)
            sample = range(total)
            method = "exhaustive"
        else:
            nsamp = 5000
            sample = range(nsamp)
            method = f"sampled({nsamp})"

        structures = defaultdict(int)

        for idx in sample:
            if n <= 6:
                A = bits_to_adj(idx, n)
            else:
                A = random_tournament(n, rng)

            if compute_b1(A, n) != 0:
                continue

            cycles = find_3cycles(A, n)
            if not cycles:
                continue

            free_cycs = [cyc for cyc in cycles if not get_dom(A, n, cyc)]
            multi_dom_cycs = [cyc for cyc in cycles if len(get_dom(A, n, cyc)) >= 2]

            # Check for free-multidom edge
            has_edge = any(shared_directed_edge(fc, mc)
                         for fc in free_cycs for mc in multi_dom_cycs)
            if has_edge:
                continue

            # Analyze structure
            n_free = len(free_cycs)
            n_udom = len([cyc for cyc in cycles if len(get_dom(A, n, cyc)) == 1])
            n_mdom = len(multi_dom_cycs)

            # How many uniquely-dominated cycles are adjacent to free cycles?
            udom_adj_free = 0
            for cyc in cycles:
                d = get_dom(A, n, cyc)
                if len(d) != 1:
                    continue
                if any(shared_directed_edge(cyc, fc) for fc in free_cycs):
                    udom_adj_free += 1

            structures[(n_free, n_udom, n_mdom, udom_adj_free)] += 1

            # Also: count bad vertices
            bad_count = sum(1 for v in range(n)
                          if compute_b1(delete_vertex(A, n, v), n-1) > 0)
            structures[f"bad={bad_count}"] += 1

        print(f"\nn={n} ({method}): No free-multidom edge cases")
        print(f"  (free, unique_dom, multi_dom, udom_adj_free) distribution:")
        for k, v in sorted(structures.items()):
            if isinstance(k, tuple):
                print(f"    free={k[0]}, udom={k[1]}, mdom={k[2]}, udom_adj_free={k[3]}: {v}")
        print(f"  Bad vertex count:")
        for k, v in sorted(structures.items()):
            if isinstance(k, str) and k.startswith("bad"):
                print(f"    {k}: {v}")

    # PART 3: In the no-free-multidom-edge case, the free cycles only neighbor
    # uniquely-dominated cycles. The unique dominator of each such cycle is a
    # single vertex d. When d is removed, the cycle becomes free.
    # For v != d and v not in the free cycle: the free cycle is in freed(v),
    # and the uniquely-dominated cycle is STILL dominated (by d, since d != v).
    # So there IS an isolation edge! Wait, no... dom(D) = {d}, and v != d,
    # so D is still dominated by d in T\v. But D is adjacent to a free cycle F.
    # If v is not in F: F is freed for v. D is dom-remaining for v.
    # F and D share an edge. So iso_edges(v) > 0.
    # Therefore v is good!
    #
    # The vertices NOT guaranteed good: v in F (3 vertices) or v = d (1 vertex).
    # Total at-risk: at most 3+1 = 4 per (F, D) pair. Since n >= 5, n-4 >= 1.
    #
    # But WAIT: we need d not in F. If d is in F, then d is both in the free
    # cycle F and the unique dominator of D. But d can't be IN F (F is a cycle
    # that d dominates from outside). Actually no: d dominates D, not F.
    # D is uniquely dominated by d. F is free (no dominator).
    # F and D are adjacent (share a directed edge).
    # d is outside D (by definition of external dominator).
    # d could be inside F or outside F.
    # If d is in F: then d is a vertex of the free cycle F.
    #   For v = d: F is removed from freed(v) since d is in F.
    #   D becomes free since dom(D) = {d} and d is removed.
    #   But D is adjacent to... other free cycles perhaps.
    # If d is outside F: then for v = d:
    #   F is in freed(d) (F is free and d not in F).
    #   D is in freed(d) (dom(D)={d}, so D becomes free).
    #   F and D are in the same freed set, connected by their shared edge.
    #   But we're looking at whether d is bad, not whether F-D is an isolation edge.
    #
    # The key point: for any v that is NOT in F and is NOT the unique dominator d:
    #   F is in freed(v), D is in dom_remaining(v), they share an edge.
    #   So iso_edges(v) > 0, hence v is good.
    #
    # Vertices at risk: v in F (3 vertices) or v = d (1 vertex, possibly in F).
    # If d is in F: at-risk = 3 (the 3 vertices of F, one of which is d).
    # If d not in F: at-risk = 4 (3 from F + 1 for d).
    #
    # In both cases: n - at_risk >= n - 4 >= 1 good vertices (for n >= 5).
    #
    # BUT: this argument works for a SINGLE (F, D) pair. Different (F, D) pairs
    # give different at-risk sets. We need SOME v not at-risk for ANY such pair.

    # Actually, we just need ANY free cycle F adjacent to ANY dominated cycle D
    # (whether uniquely dom or multi-dom). The argument above works for multi-dom
    # case directly. For unique-dom case, the argument gives n-4 guaranteed good.
    # For multi-dom case, the argument gives n-3 guaranteed good.

    # So: if ANY free cycle is adjacent to ANY dominated cycle, we get good vertices.
    # When can a free cycle NOT be adjacent to any dominated cycle?
    # Only when the free cycle is in an all-free connected component!
    # But b1(T) = 0 means no all-free component exists.
    # Therefore every free cycle IS adjacent to some dominated cycle
    # (possibly via chain of free cycles, but eventually reaches dominated).

    # Wait, the free cycle doesn't need to be DIRECTLY adjacent to a dominated cycle.
    # It just needs to be in the same connected component, which it is (b1=0).
    # But for the isolation-edge argument, we need DIRECT adjacency.

    # So the question is: does every free cycle have a dominated neighbor?
    # If the free cycles form a connected subgraph, and this subgraph connects
    # to at least one dominated cycle (since b1=0), then SOME free cycle is
    # directly adjacent to a dominated cycle.

    # Therefore: if b1(T) = 0 and free cycles exist, there always exists
    # a free cycle F directly adjacent to a dominated cycle D.

    print("\n" + "=" * 70)
    print("PART 3: KEY LEMMA - free cycle always adj to dominated cycle (given b1=0)")
    print("-" * 50)

    for n in [5, 6, 7, 8, 9]:
        if n <= 6:
            total = 2**(n*(n-1)//2)
            sample = range(total)
            method = "exhaustive"
        else:
            nsamp = {7: 10000, 8: 5000, 9: 2000}[n]
            sample = range(nsamp)
            method = f"sampled({nsamp})"

        always_adj = 0
        not_adj = 0
        no_free = 0
        total_b1_0 = 0

        for idx in sample:
            if n <= 6:
                A = bits_to_adj(idx, n)
            else:
                A = random_tournament(n, rng)

            if compute_b1(A, n) != 0:
                continue
            total_b1_0 += 1

            cycles = find_3cycles(A, n)
            if not cycles:
                no_free += 1
                continue

            free_cycs = [cyc for cyc in cycles if not get_dom(A, n, cyc)]
            if not free_cycs:
                no_free += 1
                continue

            dom_cycs = [cyc for cyc in cycles if get_dom(A, n, cyc)]

            # Does any free cycle share a directed edge with any dominated cycle?
            found = any(shared_directed_edge(fc, dc)
                       for fc in free_cycs for dc in dom_cycs)

            if found:
                always_adj += 1
            else:
                not_adj += 1
                print(f"  COUNTEREXAMPLE at n={n}! No free-dom adj edge")

        print(f"n={n} ({method}): {total_b1_0} b1=0 tournaments")
        print(f"  free adj dom: {always_adj}")
        print(f"  no free: {no_free}")
        print(f"  free NOT adj dom: {not_adj}")

if __name__ == '__main__':
    main()
