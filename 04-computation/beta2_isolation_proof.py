"""
Correct characterization of bad vertices:

v is bad iff in the restricted cycle graph G' (cycles not through v),
there exists a connected component C' such that:
1. All cycles in C' are free in T\v (i.e., free in T or dom={v})
2. C' spans n-1 vertices
3. C' is NOT connected to any dominated cycle in T\v

Condition 3 is the ISOLATION requirement: the freed component must be
disconnected from all remaining dominated cycles in G'.

Key question: how many edges connect freed cycles to remaining dominated
cycles in G'? Each such edge PREVENTS the freed set from being isolated.

PROOF STRATEGY: Show that for SOME vertex v, the freed set has many edges
to remaining dominated cycles, preventing isolation.
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
    print("ISOLATION ANALYSIS: freed set must be isolated from remaining dom")
    print("=" * 70)

    # PART 1: For each bad vertex, verify isolation.
    # For each good vertex, show non-isolation.
    print("\nPART 1: Verify bad = freed set isolated")
    print("-" * 50)

    for n in [5, 6, 7]:
        if n <= 6:
            total = 2**(n*(n-1)//2)
            sample = range(total)
            method = "exhaustive"
        else:
            nsamp = 5000
            sample = range(nsamp)
            method = f"sampled({nsamp})"

        match = 0
        mismatch = 0

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

            for v in range(n):
                # Compute freed and remaining-dom in restricted graph
                freed_idx = []
                dom_remaining_idx = []
                restricted = []

                for i, cyc in enumerate(cycles):
                    if v in set(cyc):
                        continue
                    restricted.append(cyc)
                    doms = get_dom(A, n, cyc)
                    # Free in T\v: free in T (doms empty) or dom={v}
                    doms_minus_v = [d for d in doms if d != v]
                    if not doms_minus_v:
                        freed_idx.append(len(restricted) - 1)
                    else:
                        dom_remaining_idx.append(len(restricted) - 1)

                # Build adjacency in restricted graph
                nr = len(restricted)
                r_adj = [set() for _ in range(nr)]
                for i in range(nr):
                    for j in range(i+1, nr):
                        if shared_directed_edge(restricted[i], restricted[j]):
                            r_adj[i].add(j)
                            r_adj[j].add(i)

                # Check: is there a connected component of freed cycles
                # that is isolated from dom_remaining?
                freed_set = set(freed_idx)
                dom_set = set(dom_remaining_idx)

                # Find connected components among freed cycles
                freed_visited = set()
                has_isolated_spanning_comp = False

                for start in freed_idx:
                    if start in freed_visited:
                        continue
                    # BFS only through freed cycles
                    comp = []
                    stack = [start]
                    touches_dom = False
                    while stack:
                        x = stack.pop()
                        if x in freed_visited:
                            continue
                        freed_visited.add(x)
                        comp.append(x)
                        for y in r_adj[x]:
                            if y in freed_set and y not in freed_visited:
                                stack.append(y)
                            if y in dom_set:
                                touches_dom = True

                    if not touches_dom:
                        # This component is isolated from dom cycles
                        span = set()
                        for ci in comp:
                            span.update(restricted[ci])
                        if len(span) >= n - 1:
                            has_isolated_spanning_comp = True
                            break

                # Actual test
                B = delete_vertex(A, n, v)
                actually_bad = compute_b1(B, n-1) > 0

                if has_isolated_spanning_comp == actually_bad:
                    match += 1
                else:
                    mismatch += 1
                    if mismatch <= 5:
                        print(f"  MISMATCH n={n}: v={v}, predicted={has_isolated_spanning_comp}, "
                              f"actual={actually_bad}")
                        print(f"    freed={len(freed_idx)}, dom_rem={len(dom_remaining_idx)}")

        print(f"n={n} ({method}): {match} match, {mismatch} mismatch")

    # PART 2: Count "isolation edges" per vertex
    # isolation_edges(v) = #{edges between freed(v) and remaining dom in restricted graph}
    # If isolation_edges(v) > 0, then v is likely good (freed set not isolated)
    print("\n" + "=" * 70)
    print("PART 2: Isolation edge count (freed-to-dom edges in restricted graph)")
    print("-" * 50)

    for n in [6, 7, 8]:
        if n <= 6:
            total = 2**(n*(n-1)//2)
            sample = range(total)
            method = "exhaustive"
        else:
            nsamp = {7: 5000, 8: 2000}[n]
            sample = range(nsamp)
            method = f"sampled({nsamp})"

        iso_edges_bad = defaultdict(int)
        iso_edges_good = defaultdict(int)

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

            for v in range(n):
                freed = []
                dom_rem = []
                restricted = []

                for cyc in cycles:
                    if v in set(cyc):
                        continue
                    restricted.append(cyc)
                    doms = get_dom(A, n, cyc)
                    doms_minus_v = [d for d in doms if d != v]
                    if not doms_minus_v:
                        freed.append(len(restricted) - 1)
                    else:
                        dom_rem.append(len(restricted) - 1)

                # Count edges between freed and dom_rem
                freed_s = set(freed)
                dom_s = set(dom_rem)
                iso_edges = 0
                for fi in freed:
                    for di in dom_rem:
                        if shared_directed_edge(restricted[fi], restricted[di]):
                            iso_edges += 1

                B = delete_vertex(A, n, v)
                is_bad = compute_b1(B, n-1) > 0

                if is_bad:
                    iso_edges_bad[iso_edges] += 1
                else:
                    iso_edges_good[iso_edges] += 1

        print(f"\nn={n} ({method}):")
        print(f"  BAD vertices - isolation edges: {dict(sorted(iso_edges_bad.items()))}")
        print(f"  GOOD vertices - isolation edges: {dict(sorted(iso_edges_good.items()))}")

        if iso_edges_bad:
            max_bad = max(iso_edges_bad.keys())
            print(f"  Max isolation edges for BAD: {max_bad}")
        if iso_edges_good:
            min_good_with_edges = min(k for k in iso_edges_good.keys() if k > 0) if any(k > 0 for k in iso_edges_good.keys()) else "N/A"
            print(f"  Min isolation edges > 0 for GOOD: {min_good_with_edges}")

    # PART 3: Key test - does max(isolation_edges) always pick a good vertex?
    print("\n" + "=" * 70)
    print("PART 3: Selection rule: max isolation edges")
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

        success = 0
        total_tested = 0

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

            # Check if any bad vertex exists
            has_bad = False
            for v in range(n):
                B = delete_vertex(A, n, v)
                if compute_b1(B, n-1) > 0:
                    has_bad = True
                    break
            if not has_bad:
                continue

            total_tested += 1

            # Compute isolation edges for each vertex
            iso_per_v = []
            for v in range(n):
                restricted = []
                freed = []
                dom_rem = []
                for cyc in cycles:
                    if v in set(cyc):
                        continue
                    restricted.append(cyc)
                    doms = get_dom(A, n, cyc)
                    doms_minus_v = [d for d in doms if d != v]
                    if not doms_minus_v:
                        freed.append(len(restricted) - 1)
                    else:
                        dom_rem.append(len(restricted) - 1)

                iso_edges = sum(1 for fi in freed for di in dom_rem
                               if shared_directed_edge(restricted[fi], restricted[di]))
                iso_per_v.append(iso_edges)

            # Select vertex with MAX isolation edges
            selected = max(range(n), key=lambda v: iso_per_v[v])
            B = delete_vertex(A, n, selected)
            if compute_b1(B, n-1) == 0:
                success += 1

        pct = 100*success/max(total_tested, 1)
        marker = " ***" if pct == 100 else ""
        print(f"  n={n} ({method}): {success}/{total_tested} ({pct:.1f}%){marker}")

    # PART 4: Test combined rule: max(isolation_edges), break ties by min(freed_count)
    print("\n" + "=" * 70)
    print("PART 4: Combined selection rule")
    print("-" * 50)

    for n in [7, 8, 9]:
        nsamp = {7: 10000, 8: 5000, 9: 2000}[n]
        sample = range(nsamp)

        rules = {
            'max_iso': lambda ie, fc: (ie[0], -fc[0]),  # wrong, use tuple version
        }

        success_iso = 0
        success_combo = 0
        total_tested = 0

        for _ in sample:
            A = random_tournament(n, rng)
            if compute_b1(A, n) != 0:
                continue

            cycles = find_3cycles(A, n)
            if not cycles:
                continue

            has_bad = False
            for v in range(n):
                B = delete_vertex(A, n, v)
                if compute_b1(B, n-1) > 0:
                    has_bad = True
                    break
            if not has_bad:
                continue
            total_tested += 1

            iso_per_v = []
            freed_per_v = []
            for v in range(n):
                restricted = []
                freed = []
                dom_rem = []
                for cyc in cycles:
                    if v in set(cyc):
                        continue
                    restricted.append(cyc)
                    doms = get_dom(A, n, cyc)
                    doms_minus_v = [d for d in doms if d != v]
                    if not doms_minus_v:
                        freed.append(len(restricted) - 1)
                    else:
                        dom_rem.append(len(restricted) - 1)

                iso_edges = sum(1 for fi in freed for di in dom_rem
                               if shared_directed_edge(restricted[fi], restricted[di]))
                iso_per_v.append(iso_edges)
                freed_per_v.append(len(freed))

            # Rule 1: max isolation edges
            sel1 = max(range(n), key=lambda v: iso_per_v[v])
            B = delete_vertex(A, n, sel1)
            if compute_b1(B, n-1) == 0:
                success_iso += 1

            # Rule 2: max isolation edges, break ties by min freed
            sel2 = max(range(n), key=lambda v: (iso_per_v[v], -freed_per_v[v]))
            B = delete_vertex(A, n, sel2)
            if compute_b1(B, n-1) == 0:
                success_combo += 1

        print(f"n={n}: max_iso={success_iso}/{total_tested} "
              f"({100*success_iso/max(total_tested,1):.1f}%), "
              f"combo={success_combo}/{total_tested} "
              f"({100*success_combo/max(total_tested,1):.1f}%)")

if __name__ == '__main__':
    main()
