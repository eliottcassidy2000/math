"""
iterated_dc_expansion.py
kind-pasteur-2026-03-07-S35

Iterate deletion-contraction to decompose H(T) as a sum over edge subsets.

For edge set E = {e_1, ..., e_m} of tournament T (m = n(n-1)/2):
  H(T) = sum_{S in 2^E} H_ham(T//S)

where T//S is the minor obtained by contracting all edges in S and deleting all
edges in E\\S. H_ham counts Hamiltonian paths in the minor.

The minor T//S has at most n - |S| vertices (less if contractions merge vertices
already merged).

Most terms are 0 (minors without Ham paths). Which subsets S contribute?
- S must include a set of edges forming a Hamiltonian path in the contraction structure.

We explore:
1. What fraction of 2^m terms are nonzero?
2. Which subsets S contribute? Are they related to spanning paths?
3. Is there a nice closed-form or interpretation of the expansion?
"""

import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
from itertools import permutations, combinations
from collections import defaultdict


def tournament_from_bits(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj


def ham_paths_dp(adj, n):
    if n == 0:
        return 1
    if n == 1:
        return 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if dp[mask][v] == 0 or not (mask & (1 << v)):
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))


def get_directed_edges(adj, n):
    """Return list of directed edges as (u, v, undirected_id)."""
    edges = []
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if adj[i][j]:
                edges.append((i, j, idx))
            else:
                edges.append((j, i, idx))
            idx += 1
    return edges


def compute_minor(adj, n, contract_set, edges):
    """
    Compute the minor T//S where S is the set of edges to contract.

    Process contractions one by one. Each contraction merges two vertices.
    Use Union-Find to track which original vertices are merged.

    After all contractions, build the minor:
    - Vertices = equivalence classes of Union-Find
    - For each edge NOT in S: it's deleted (no edge in minor)
    - For edges in S: they're contracted (vertices merged)
    - Between two classes A, B: edge (A,B) iff the contracted edge goes A->B

    Wait, this needs more thought. In deletion-contraction:
    - Contracted edges: the merged vertex inherits IN from tail, OUT from head
    - But with multiple contractions, order matters if edges share vertices

    Actually for the Tutte-like expansion, the minor T//S should be:
    - Vertices: connected components of the graph induced by S edges
      (using undirected connectivity)
    - Edges: for each contracted edge (u,v) in S, merge u and v
    - After merging, remaining edges in S that connect already-merged vertices
      become self-loops (ignored for Ham paths)

    Let me use a simpler approach: process edges in order.
    """
    # Union-Find
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    # Track edge directions for the contracted edges
    # For each contracted edge (u,v), w = merged = inherits IN from u, OUT from v
    # With multiple contractions, we need to track the "representative" vertex

    # Process contractions
    contraction_order = []
    for eidx in contract_set:
        u, v, _ = edges[eidx]
        pu, pv = find(u), find(v)
        if pu != pv:
            # Merge: the representative inherits IN from u's rep, OUT from v's rep
            union(u, v)
            contraction_order.append((u, v))

    # Find unique classes
    classes = {}
    for i in range(n):
        root = find(i)
        if root not in classes:
            classes[root] = len(classes)

    new_n = len(classes)
    if new_n == 0:
        return [[]], 0

    # Build adjacency for the minor
    # For the minor, we need to determine edge directions between classes
    # This is tricky with multiple contractions. Let's use a voting/priority scheme.

    # Approach: for each pair of classes (A, B), the edge direction is determined
    # by the FIRST contraction that bridges them, using the convention:
    # IN from tail's class, OUT from head's class.

    # Actually, simpler: for each original edge (u,v) NOT in S (deleted),
    # we have NO edge in the minor. For edges IN S (contracted), the vertices
    # are merged. Between two classes, there's an edge iff... hmm.

    # Actually in deletion-contraction, ALL non-contracted edges are DELETED.
    # The minor has ONLY the contracted edges (which become internal to classes)
    # plus... wait, no.

    # In the standard Tutte expansion: minor T//S = contract S, delete E\S.
    # After contracting all edges in S, we get a graph on k = new_n vertices.
    # Then we delete all edges in E\S from this contracted graph.
    # But the contracted graph inherits edges from T that connect different classes.
    # These non-S edges are then deleted.

    # So the minor T//S is a graph on k vertices with NO edges (all deleted).
    # H(no edges, k vertices) = 0 for k >= 2, 1 for k = 1.

    # Wait, that can't be right. The Tutte expansion sums T//S over all S,
    # and most terms would be 0. Let me reconsider.

    # Actually, in the standard Tutte polynomial, the expansion is:
    # T(G;x,y) = sum_A (x-1)^{r(E)-r(A)} (y-1)^{|A|-r(A)}
    # where r is the rank function. This is NOT the same as deletion-contraction
    # iterated naively.

    # For our H(T) deletion-contraction, iterating gives:
    # H(T) = H(T\e1) + H(T/e1)
    # For each term, we apply DC to the NEXT edge.
    # The result depends on order because contraction changes the graph.

    # Let me just compute H(T) by iterating DC in a fixed order and tracking
    # what the final terms look like.

    # Simple recursive approach for small n:
    return None, new_n  # placeholder


def iterate_dc_recursive(adj, n, edges, depth=0, label=""):
    """
    Recursively apply deletion-contraction.
    Returns list of (terminal_adj, terminal_n, sign, label) tuples.
    """
    if not edges:
        # No more edges: terminal graph
        return [(adj, n, label)]

    # Pick next edge
    u, v, eidx = edges[0]
    remaining = edges[1:]

    # DELETION: remove edge u->v
    adj_del = [row[:] for row in adj]
    adj_del[u][v] = 0

    # Map remaining edges (no vertex change)
    del_remaining = []
    for (eu, ev, ei) in remaining:
        del_remaining.append((eu, ev, ei))

    del_terms = iterate_dc_recursive(adj_del, n, del_remaining, depth+1, label + "D")

    # CONTRACTION: merge u,v into w
    # w inherits IN from u, OUT from v
    others = [x for x in range(n) if x != u and x != v]
    new_n = n - 1

    # Map old vertices to new
    w_new = 0
    omap = {x: i+1 for i, x in enumerate(others)}

    con_adj = [[0]*new_n for _ in range(new_n)]
    for x in others:
        for y in others:
            con_adj[omap[x]][omap[y]] = adj[x][y]
    for x in others:
        if adj[x][u]:
            con_adj[omap[x]][w_new] = 1
        if adj[v][x]:
            con_adj[w_new][omap[x]] = 1

    # Map remaining edges to new vertex indices
    con_remaining = []
    for (eu, ev, ei) in remaining:
        if eu == u or eu == v:
            new_eu = w_new
        else:
            new_eu = omap.get(eu, -1)
        if ev == u or ev == v:
            new_ev = w_new
        else:
            new_ev = omap.get(ev, -1)
        if new_eu == -1 or new_ev == -1:
            continue
        if new_eu == new_ev:
            continue  # self-loop from contraction
        con_remaining.append((new_eu, new_ev, ei))

    con_terms = iterate_dc_recursive(con_adj, new_n, con_remaining, depth+1, label + "C")

    return del_terms + con_terms


def test_iterated_dc(n):
    """Test iterated DC expansion at small n."""
    num_bits = n * (n - 1) // 2

    print(f"\n=== n={n}: Iterated DC Expansion ===")

    # Just test a few tournaments
    for bits in range(min(8, 1 << num_bits)):
        adj = tournament_from_bits(n, bits)
        H_T = ham_paths_dp(adj, n)
        edges = get_directed_edges(adj, n)

        terminals = iterate_dc_recursive(adj, n, edges)

        # Compute H for each terminal
        total = 0
        nonzero = 0
        for (t_adj, t_n, label) in terminals:
            h = ham_paths_dp(t_adj, t_n)
            total += h
            if h > 0:
                nonzero += 1

        num_terms = len(terminals)
        print(f"  bits={bits}: H(T)={H_T}, DC expansion sum={total}, "
              f"terms={num_terms}, nonzero={nonzero}/{num_terms}")

        if total != H_T:
            print(f"  MISMATCH!")

        # Analyze terminal structure
        if bits < 2:
            sizes = defaultdict(int)
            for (t_adj, t_n, label) in terminals:
                h = ham_paths_dp(t_adj, t_n)
                sizes[(t_n, h)] += 1
            print(f"    Terminal (size, H) counts: {dict(sizes)}")


def analyze_terminal_structure(n):
    """Analyze what the terminal graphs look like after full DC expansion."""
    num_bits = n * (n - 1) // 2

    print(f"\n=== n={n}: Terminal Structure Analysis ===")

    adj = tournament_from_bits(n, 0)
    edges = get_directed_edges(adj, n)
    terminals = iterate_dc_recursive(adj, n, edges)

    # Statistics
    h_dist = defaultdict(int)
    size_dist = defaultdict(int)

    for (t_adj, t_n, label) in terminals:
        h = ham_paths_dp(t_adj, t_n)
        h_dist[h] += 1
        size_dist[t_n] += 1

    print(f"  Total terminal terms: {len(terminals)}")
    print(f"  Size distribution: {dict(sorted(size_dist.items()))}")
    print(f"  H distribution: {dict(sorted(h_dist.items()))}")

    # How many contractions in each terminal?
    c_dist = defaultdict(int)
    for (t_adj, t_n, label) in terminals:
        num_c = label.count('C')
        c_dist[num_c] += 1
    print(f"  # contractions distribution: {dict(sorted(c_dist.items()))}")

    # Verify: all terms with H > 0 have terminal size = 1
    all_size1 = all(t_n == 1 for (t_adj, t_n, label) in terminals if ham_paths_dp(t_adj, t_n) > 0)
    print(f"  All nonzero terminals have size 1: {all_size1}")


# Run tests
for n in [3, 4]:
    test_iterated_dc(n)

for n in [3, 4]:
    analyze_terminal_structure(n)

print("\nDONE")
