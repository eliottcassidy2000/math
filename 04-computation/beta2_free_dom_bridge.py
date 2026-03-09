"""
Precise characterization of bad vertices via free-dom bridge.

v is bad iff: removing v's cycles from the cycle graph disconnects some set of
free cycles from ALL dominated cycles, forming a spanning free component.

Formally: In the cycle graph G, let F = free cycles, D = dominated cycles.
Consider the bipartite graph B(F, D) where F_i ~ D_j if they share a directed edge
AND at least one of them goes through v (i.e., the connection goes via v).

Wait, the correct formulation: v is bad iff in the restricted cycle graph
G' = G[cycles not through v], there's a connected component C' that is all-free
AND spans n-1 vertices.

The all-free condition: C' consists only of free cycles (cycles that remain
free after removing v). A cycle C not through v is free in T\v iff it was free
in T (since unique_dom(v)=0, if v doesn't uniquely dominate it, its domination
status doesn't change... wait, it DOES change if v was one of multiple dominators).

CORRECTION: Domination in T\v is different from T!
- A cycle C not through v is dominated in T\v iff some vertex d != v, d not in C,
  dominates C in T\v (which means d dominates C in T as well, since edges don't change).
- But: d might be in C in T\v (d was relabeled)... no, C's vertex set doesn't change.
- Actually: C is dominated in T iff some d outside C has d->all of C or all of C->d.
  C is dominated in T\v iff some d != v outside C has d->all of C or all of C->d.
  The ONLY difference: v can no longer serve as dominator.

So: C is dominated in T\v iff C is dominated in T by some vertex other than v.
Equivalently: C is free in T\v iff C was free in T OR dom_T(C) = {v} (uniquely v-dominated).

For v to be bad: the free-in-T\v cycles not through v must include a spanning connected component.

Free-in-T\v not through v = {free-in-T not through v} UNION {dom={v} not through v}

So the mechanisms are unified: the freed set is free-in-T + newly-freed-by-v-deletion.

For v to be bad, this freed set must:
1. Form a connected subgraph in the restricted cycle graph
2. Span all n-1 vertices except v

This is the CORRECT and COMPLETE characterization.

Now: for a PROOF, we want to show that not all n vertices can satisfy this.

Key quantity: Let's track for each vertex v:
- freed(v) = {free cycles not through v} + {cycles with dom={v}}
- span(v) = vertices covered by freed(v)
- conn(v) = whether freed(v) is connected in restricted graph

v is bad iff |span(v)| = n-1 AND conn(v) = True.

Now consider two vertices u, v. The freed cycles for u and v share all free cycles
that go through neither u nor v. The uniquely-dominated cycles are disjoint.

Let F = set of all free cycles in T.
|F not through v| = |F| - |F thru v|

If all n vertices are bad, then for each v:
freed(v) spans V\{v} and is connected.

This means: for each pair (u, v) with u != v, there exists a cycle in freed(v)
containing u. This cycle is either:
(a) A free cycle containing u but not v
(b) A cycle with dom={v}, containing u but not v

This is the coverage requirement.

Key constraint from connectivity: the freed set must be CONNECTED.
For a graph on ~|F| + few unique-dom cycles, connectivity requires enough
edge-sharing between these cycles.

Let me compute: what fraction of vertex pairs (u, v) are connected by WHICH
type of freed cycle? And what does connectivity actually require?
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

def is_connected_cycles(cyc_list):
    """Check if a list of cycles is connected via shared directed edges."""
    if not cyc_list:
        return False
    nc = len(cyc_list)
    if nc == 1:
        return True
    adj = [[] for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if shared_directed_edge(cyc_list[i], cyc_list[j]):
                adj[i].append(j)
                adj[j].append(i)
    visited = [False]*nc
    stack = [0]
    while stack:
        v = stack.pop()
        if visited[v]: continue
        visited[v] = True
        for u in adj[v]:
            if not visited[u]: stack.append(u)
    return all(visited)

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
    print("FREE-DOM BRIDGE: Precise characterization via freed set")
    print("=" * 70)

    # PART 1: Verify: v bad iff freed(v) is connected and spans n-1 vertices
    print("\nPART 1: Verification of freed-set characterization")
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
                # Compute freed(v) = free not thru v + dom={v} not thru v
                freed = []
                for cyc in cycles:
                    if v in set(cyc):
                        continue
                    doms = get_dom(A, n, cyc)
                    if len(doms) == 0 or doms == [v]:
                        freed.append(cyc)

                # Check: spans n-1 vertices?
                spanned = set()
                for cyc in freed:
                    spanned.update(cyc)
                spans = (len(spanned) >= n - 1)

                # Check: connected?
                conn = is_connected_cycles(freed) if freed else False

                # Actual test: is v bad?
                B = delete_vertex(A, n, v)
                actually_bad = compute_b1(B, n-1) > 0

                # v is bad iff freed is connected AND spans n-1
                predicted_bad = spans and conn

                if predicted_bad == actually_bad:
                    match += 1
                else:
                    mismatch += 1
                    if mismatch <= 3:
                        print(f"  MISMATCH at n={n}: v={v}, predicted={predicted_bad}, actual={actually_bad}")
                        print(f"    freed={len(freed)} cycles, spans={len(spanned)}, conn={conn}")

        print(f"n={n} ({method}): {match} match, {mismatch} mismatch")

    # PART 2: Given the freed-set characterization, analyze what limits |BAD|.
    # For v bad: freed(v) connected + spans n-1.
    # freed(v) = F\v + D_v where F\v = free not thru v, D_v = uniquely-v-dom not thru v
    #
    # Key: F\v is the SAME for all v (minus different exclusions).
    # F\v = F - {cycles in F that go through v}
    #
    # If all n vertices are bad:
    # - Each freed(v) is connected and spans n-1.
    # - freed(v) = (F minus F_thru_v) union D_v
    # - For EACH pair (u, w) with u,w != v: u and w are connected in freed(v)
    #   (connected component means any two vertices in the span are reachable)
    #
    # Consider the FREE CYCLES graph: G_F = subgraph of cycle graph restricted to F.
    # For vertex v: the restricted free graph G_F\v = G_F minus cycles through v.
    # If this is connected and spans n-1 vertices, AND b1=0 requires no free component
    # in T, then G_F is NOT all of these things simultaneously for ALL v.
    #
    # Actually G_F might not be connected in T (b1=0 means it's mixed with dom cycles).
    # So freed(v) = G_F\v + D_v.

    print("\n" + "=" * 70)
    print("PART 2: Free cycle graph properties")
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

        stats = defaultdict(int)

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

            # Free cycles
            free_cycles = [cyc for cyc in cycles if not get_dom(A, n, cyc)]

            # Free cycle graph: connected?
            if free_cycles:
                f_connected = is_connected_cycles(free_cycles)
                f_span = set()
                for cyc in free_cycles:
                    f_span.update(cyc)
                f_span_size = len(f_span)
            else:
                f_connected = False
                f_span_size = 0

            stats[f"f_connected={f_connected}"] += 1
            stats[f"f_span={f_span_size}"] += 1

            # Compute bad count
            bad_count = 0
            for v in range(n):
                B = delete_vertex(A, n, v)
                if compute_b1(B, n-1) > 0:
                    bad_count += 1
            stats[f"bad={bad_count}"] += 1

            # Key question: if F spans all n vertices and is connected,
            # can ALL vertices be bad?
            if f_span_size == n and f_connected and bad_count == n:
                stats["F_full_and_all_bad"] += 1

            # For each v: does F\v span n-1 and is F\v connected?
            all_v_F_ok = True
            for v in range(n):
                f_not_v = [cyc for cyc in free_cycles if v not in set(cyc)]
                if not f_not_v:
                    all_v_F_ok = False
                    break
                s = set()
                for cyc in f_not_v:
                    s.update(cyc)
                if len(s) < n - 1:
                    all_v_F_ok = False
                    break
                if not is_connected_cycles(f_not_v):
                    all_v_F_ok = False
                    break

            stats[f"all_v_F_ok={all_v_F_ok}"] += 1

        print(f"\nn={n} ({method}):")
        for k in sorted(stats.keys()):
            print(f"  {k}: {stats[k]}")

    # PART 3: Connection between free cycles spanning and badness
    # Hypothesis: if free cycles alone span n-1 and are connected for v,
    # then v is bad (since freed(v) >= F\v which already spans and connects).
    # Is this true?
    print("\n" + "=" * 70)
    print("PART 3: Does F\\v spanning+connected imply v is bad?")
    print("-" * 50)

    for n in [5, 6, 7]:
        if n <= 6:
            total = 2**(n*(n-1)//2)
            sample = range(total)
        else:
            nsamp = 5000
            sample = range(nsamp)

        implies = 0
        not_implies = 0

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

            free_cycles = [cyc for cyc in cycles if not get_dom(A, n, cyc)]

            for v in range(n):
                f_not_v = [cyc for cyc in free_cycles if v not in set(cyc)]
                if not f_not_v:
                    continue
                s = set()
                for cyc in f_not_v:
                    s.update(cyc)
                if len(s) < n - 1:
                    continue
                if not is_connected_cycles(f_not_v):
                    continue

                # F\v spans n-1 and is connected. Is v bad?
                B = delete_vertex(A, n, v)
                if compute_b1(B, n-1) > 0:
                    implies += 1
                else:
                    not_implies += 1

        print(f"  n={n}: F\\v spans+connected => v bad: {implies} yes, {not_implies} no")

    # PART 4: What about the converse? v bad implies F\v spans+connected?
    print("\n" + "=" * 70)
    print("PART 4: v bad implies F\\v spans n-1?")
    print("-" * 50)

    for n in [5, 6, 7]:
        if n <= 6:
            total = 2**(n*(n-1)//2)
            sample = range(total)
        else:
            nsamp = 5000
            sample = range(nsamp)

        bad_f_spans = 0
        bad_f_not_spans = 0

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

            free_cycles = [cyc for cyc in cycles if not get_dom(A, n, cyc)]

            for v in range(n):
                B = delete_vertex(A, n, v)
                if compute_b1(B, n-1) == 0:
                    continue  # v is good

                # v is bad. Does F\v span n-1?
                f_not_v = [cyc for cyc in free_cycles if v not in set(cyc)]
                s = set()
                for cyc in f_not_v:
                    s.update(cyc)
                if len(s) >= n - 1:
                    bad_f_spans += 1
                else:
                    bad_f_not_spans += 1

        print(f"  n={n}: v bad => F\\v spans n-1: {bad_f_spans} yes, {bad_f_not_spans} no")
        if bad_f_not_spans > 0:
            print(f"    NOT TRUE: bad vertex can exist without F\\v spanning")

    # PART 5: Summary of structural constraints
    print("\n" + "=" * 70)
    print("PART 5: Freed-set counting at n=6")
    print("-" * 50)

    n = 6
    total = 2**(n*(n-1)//2)

    freed_size_bad = defaultdict(int)
    freed_size_good = defaultdict(int)
    free_count_dist = defaultdict(int)
    udom_count_dist = defaultdict(int)

    for bits in range(total):
        A = bits_to_adj(bits, n)
        if compute_b1(A, n) != 0:
            continue

        cycles = find_3cycles(A, n)
        if not cycles:
            continue

        free_cycles = [cyc for cyc in cycles if not get_dom(A, n, cyc)]
        free_count_dist[len(free_cycles)] += 1

        # Unique dom cycles per vertex
        for v in range(n):
            udom = sum(1 for cyc in cycles if v not in set(cyc) and get_dom(A, n, cyc) == [v])
            udom_count_dist[udom] += 1

        for v in range(n):
            freed = []
            for cyc in cycles:
                if v in set(cyc):
                    continue
                doms = get_dom(A, n, cyc)
                if len(doms) == 0 or doms == [v]:
                    freed.append(cyc)

            B = delete_vertex(A, n, v)
            is_bad = compute_b1(B, n-1) > 0

            if is_bad:
                freed_size_bad[len(freed)] += 1
            else:
                freed_size_good[len(freed)] += 1

    print(f"  Free cycle count distribution: {dict(free_count_dist)}")
    print(f"  Unique-dom per vertex distribution: {dict(udom_count_dist)}")
    print(f"  |freed(v)| for BAD v: {dict(freed_size_bad)}")
    print(f"  |freed(v)| for GOOD v: {dict(freed_size_good)}")

if __name__ == '__main__':
    main()
