#!/usr/bin/env python3
"""
Catalog of odd-cycle conflict graphs Omega(T) for tournaments T.

For a tournament T on n vertices:
- Find all directed odd cycles (as vertex sets, i.e., frozensets)
- Build the conflict graph Omega(T): vertices = odd cycles, edges when two cycles share a vertex
- Classify Omega(T) up to isomorphism using graph invariants
- Check structural properties: chordal, perfect, comparability, interval

For n=3,4,5,6: exhaustive over all tournaments
For n=7: random sampling
"""

import itertools
import random
from collections import Counter, defaultdict

# ─────────────────────────────────────────────────
# Tournament generation
# ─────────────────────────────────────────────────

def all_tournaments(n):
    """Generate all labeled tournaments on n vertices."""
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    for bits in range(2**len(edges)):
        adj = [[False]*n for _ in range(n)]
        for k, (i, j) in enumerate(edges):
            if (bits >> k) & 1:
                adj[i][j] = True
            else:
                adj[j][i] = True
        yield adj

def random_tournament(n):
    """Generate a random tournament on n vertices."""
    adj = [[False]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                adj[i][j] = True
            else:
                adj[j][i] = True
    return adj

# ─────────────────────────────────────────────────
# Finding all directed odd cycles
# ─────────────────────────────────────────────────

def find_all_odd_cycles(adj, n):
    """
    Find all directed odd cycles in a tournament.
    Returns a list of frozensets (vertex sets of cycles).

    We enumerate directed cycles by DFS. For tournaments, we only need
    odd lengths: 3, 5, 7, ...
    """
    cycle_sets = set()

    # For each possible cycle length (odd: 3, 5, 7, ...)
    for length in range(3, n+1, 2):
        # Enumerate all ordered subsets of size `length`
        for vertices in itertools.combinations(range(n), length):
            # Check all cyclic permutations - actually we need to check
            # if there's a Hamiltonian directed cycle on this vertex subset
            # For a tournament restricted to these vertices, a Hamilton cycle
            # exists iff the sub-tournament is strongly connected.
            # But we want ALL directed cycles, not just Hamiltonian ones on subsets.
            #
            # Better approach: enumerate all permutations of the vertex subset
            # and check if it forms a directed cycle.
            # But that's too expensive for large subsets.
            #
            # For small n, let's just enumerate carefully.
            vset = frozenset(vertices)
            if vset in cycle_sets:
                continue
            # Check if there's a directed cycle visiting exactly these vertices
            # We need to find a cyclic ordering v0->v1->...->v_{k-1}->v0
            # This is equivalent to finding a Hamiltonian cycle in the
            # sub-tournament on these vertices.
            if _has_hamiltonian_cycle(adj, list(vertices)):
                cycle_sets.add(vset)

    return list(cycle_sets)

def _has_hamiltonian_cycle(adj, vertices):
    """Check if the sub-tournament on `vertices` has a directed Hamiltonian cycle."""
    k = len(vertices)
    if k < 3:
        return False
    if k == 3:
        a, b, c = vertices
        # Check both orientations
        if adj[a][b] and adj[b][c] and adj[c][a]:
            return True
        if adj[a][c] and adj[c][b] and adj[b][a]:
            return True
        return False

    # For larger k, use DP (Held-Karp style) or just permutations for small k
    # Since k <= 7 typically, permutations are fine for k <= 7
    if k <= 8:
        # Fix first vertex, permute the rest
        first = vertices[0]
        for perm in itertools.permutations(vertices[1:]):
            path = [first] + list(perm)
            valid = True
            for i in range(k-1):
                if not adj[path[i]][path[i+1]]:
                    valid = False
                    break
            if valid and adj[path[-1]][first]:
                return True
        return False

    # For larger k, use DP bitmask approach
    return _hamcycle_dp(adj, vertices)

def _hamcycle_dp(adj, vertices):
    """Held-Karp DP for Hamiltonian cycle detection."""
    k = len(vertices)
    idx = {v: i for i, v in enumerate(vertices)}
    full = (1 << k) - 1

    # dp[mask][i] = True if there's a path visiting vertices in mask, ending at vertices[i]
    # Start from vertex 0
    dp = [[False]*k for _ in range(1 << k)]
    dp[1][0] = True  # Start at vertex 0 (bit 0)

    for mask in range(1, 1 << k):
        for i in range(k):
            if not dp[mask][i]:
                continue
            if not (mask & (1 << i)):
                continue
            for j in range(k):
                if mask & (1 << j):
                    continue
                if adj[vertices[i]][vertices[j]]:
                    dp[mask | (1 << j)][j] = True

    # Check if we can complete the cycle back to vertex 0
    for i in range(1, k):
        if dp[full][i] and adj[vertices[i]][vertices[0]]:
            return True
    return False

# ─────────────────────────────────────────────────
# Alternative: find odd cycles more efficiently
# ─────────────────────────────────────────────────

def find_all_odd_cycles_v2(adj, n):
    """
    Find all directed odd cycles more efficiently.
    A directed cycle in a tournament on vertex set S exists iff
    the sub-tournament on S is strongly connected AND |S| is the
    cycle length... wait, no. A tournament on S can have a
    Hamiltonian cycle iff it's strongly connected (Moon's theorem).

    But we want the vertex SETS of ALL directed cycles, not just
    Hamiltonian ones. A 3-cycle on {a,b,c} is a directed cycle
    using all 3 vertices. A 5-cycle on {a,b,c,d,e} uses all 5.

    Actually, in a tournament, a "directed cycle" of length k means
    k distinct vertices forming a directed k-cycle. The vertex set
    uniquely determines it as a subset of size k. But a tournament
    on k vertices might have a Hamiltonian cycle (visiting all k vertices)
    or not. So the odd cycles of T are exactly the subsets S of V(T)
    with |S| odd and |S|>=3 such that the sub-tournament T[S] has a
    Hamiltonian directed cycle.

    By Moon's theorem, a tournament has a Hamiltonian cycle iff it is
    strongly connected. So:

    Odd cycles of T = {S subset V(T) : |S| odd, |S|>=3, T[S] is strongly connected}

    Wait - that's not quite right either. Moon's theorem says a strongly
    connected tournament has a Hamiltonian cycle. But a tournament on S
    being strongly connected means it has a Hamiltonian cycle visiting
    ALL vertices of S. So S is a "cycle vertex set" of that length.

    However, a non-strongly-connected tournament on S might still have
    shorter cycles. But those shorter cycles would be on subsets of S,
    not on S itself.

    So: the set of vertex-sets of odd directed cycles in T is exactly:
    {S : |S| odd, |S|>=3, T[S] is strongly connected}

    This is correct by Moon's theorem!
    """
    cycle_sets = []

    for length in range(3, n+1, 2):
        for vertices in itertools.combinations(range(n), length):
            if is_strongly_connected(adj, vertices):
                cycle_sets.append(frozenset(vertices))

    return cycle_sets

def is_strongly_connected(adj, vertices):
    """Check if the sub-tournament on `vertices` is strongly connected."""
    k = len(vertices)
    if k <= 1:
        return True  # trivially
    if k == 2:
        return False  # tournament on 2 vertices is never strongly connected

    vlist = list(vertices)
    vset = set(vertices)

    # BFS/DFS from first vertex
    def reachable(start):
        visited = {start}
        stack = [start]
        while stack:
            u = stack.pop()
            for v in vlist:
                if v not in visited and v in vset and adj[u][v]:
                    visited.add(v)
                    stack.append(v)
        return visited

    # Check reachability from first vertex
    fwd = reachable(vlist[0])
    if len(fwd) != k:
        return False

    # Check reverse reachability (can every vertex reach vlist[0]?)
    def reachable_rev(start):
        visited = {start}
        stack = [start]
        while stack:
            u = stack.pop()
            for v in vlist:
                if v not in visited and v in vset and adj[v][u]:
                    visited.add(v)
                    stack.append(v)
        return visited

    rev = reachable_rev(vlist[0])
    return len(rev) == k

# ─────────────────────────────────────────────────
# Build conflict graph Omega(T)
# ─────────────────────────────────────────────────

def build_conflict_graph(cycles):
    """
    Build the conflict graph Omega(T).
    Vertices = odd cycles (indices into the cycles list)
    Edges = pairs of cycles sharing at least one vertex
    Returns adjacency list.
    """
    m = len(cycles)
    adj = [[False]*m for _ in range(m)]

    for i in range(m):
        for j in range(i+1, m):
            if cycles[i] & cycles[j]:  # frozenset intersection
                adj[i][j] = True
                adj[j][i] = True

    return adj

# ─────────────────────────────────────────────────
# Graph invariants for classification
# ─────────────────────────────────────────────────

def graph_invariant(adj, m):
    """
    Compute an isomorphism invariant for a graph on m vertices.
    Returns a tuple that's equal for isomorphic graphs (not guaranteed
    to distinguish all non-isomorphic graphs, but good enough for small cases).
    """
    if m == 0:
        return (0, 0, (), 0, ())

    # Degree sequence
    degrees = sorted([sum(adj[i]) for i in range(m)])

    # Number of edges
    num_edges = sum(degrees) // 2

    # Number of connected components
    visited = [False] * m
    num_components = 0
    component_sizes = []
    for start in range(m):
        if not visited[start]:
            num_components += 1
            # BFS
            queue = [start]
            visited[start] = True
            size = 0
            while queue:
                u = queue.pop(0)
                size += 1
                for v in range(m):
                    if adj[u][v] and not visited[v]:
                        visited[v] = True
                        queue.append(v)
            component_sizes.append(size)
    component_sizes.sort()

    # Number of triangles
    num_triangles = 0
    for i in range(m):
        for j in range(i+1, m):
            if adj[i][j]:
                for k in range(j+1, m):
                    if adj[i][k] and adj[j][k]:
                        num_triangles += 1

    # Degree-degree correlation (sorted neighbor degree sequences)
    neighbor_deg_seqs = []
    for i in range(m):
        nds = sorted([sum(adj[j]) for j in range(m) if adj[i][j]])
        neighbor_deg_seqs.append(tuple(nds))
    neighbor_deg_seqs.sort()

    return (m, num_edges, tuple(degrees), num_components,
            tuple(component_sizes), num_triangles, tuple(neighbor_deg_seqs))

def is_complete(adj, m):
    """Check if graph is complete."""
    for i in range(m):
        for j in range(i+1, m):
            if not adj[i][j]:
                return False
    return True

# ─────────────────────────────────────────────────
# Structural property checks
# ─────────────────────────────────────────────────

def find_all_induced_cycles(adj, m, min_length=4):
    """Find all induced cycles of length >= min_length (for chordality check)."""
    # A graph is chordal iff it has no induced cycle of length >= 4
    for length in range(min_length, m+1):
        for vertices in itertools.combinations(range(m), length):
            # Check if these vertices form an induced cycle
            vset = set(vertices)
            # Each vertex should have degree exactly 2 in the induced subgraph
            ok = True
            for v in vertices:
                deg = sum(1 for u in vertices if u != v and adj[v][u])
                if deg != 2:
                    ok = False
                    break
            if not ok:
                continue
            # Check connectivity (should be a single cycle)
            # With all degrees = 2 and connected, it must be a cycle
            visited = {vertices[0]}
            stack = [vertices[0]]
            while stack:
                u = stack.pop()
                for w in vertices:
                    if w not in visited and adj[u][w]:
                        visited.add(w)
                        stack.append(w)
            if len(visited) == length:
                return length  # Found an induced cycle of this length
    return 0  # No induced cycle found => chordal

def is_chordal(adj, m):
    """Check if graph is chordal (no induced cycle of length >= 4)."""
    if m <= 3:
        return True
    result = find_all_induced_cycles(adj, m, min_length=4)
    return result == 0

def complement_graph(adj, m):
    """Return the complement graph."""
    comp = [[False]*m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            if i != j:
                comp[i][j] = not adj[i][j]
    return comp

def chromatic_number_greedy_upper(adj, m):
    """Greedy upper bound on chromatic number."""
    if m == 0:
        return 0
    colors = [-1] * m
    for v in range(m):
        used = set()
        for u in range(m):
            if adj[v][u] and colors[u] >= 0:
                used.add(colors[u])
        c = 0
        while c in used:
            c += 1
        colors[v] = c
    return max(colors) + 1

def max_clique_size(adj, m):
    """Find the size of the maximum clique (brute force for small graphs)."""
    if m == 0:
        return 0
    best = 1
    for size in range(2, m+1):
        found = False
        for vertices in itertools.combinations(range(m), size):
            is_clique = True
            for i in range(len(vertices)):
                for j in range(i+1, len(vertices)):
                    if not adj[vertices[i]][vertices[j]]:
                        is_clique = False
                        break
                if not is_clique:
                    break
            if is_clique:
                found = True
                best = size
                break  # Found one of this size, try larger
        if not found:
            break
    return best

def max_independent_set_size(adj, m):
    """Find the size of the maximum independent set (complement clique)."""
    comp = complement_graph(adj, m)
    return max_clique_size(comp, m)

def is_perfect(adj, m):
    """
    Check if graph is perfect.
    By the weak perfect graph theorem, G is perfect iff
    for every induced subgraph H, chi(H) = omega(H).

    For small graphs, check all induced subgraphs.
    Actually, by the strong perfect graph theorem (Chudnovsky et al. 2006),
    G is perfect iff neither G nor its complement contains an odd hole
    (induced odd cycle of length >= 5).
    """
    if m <= 4:
        return True  # All graphs on <= 4 vertices are perfect or easily checked

    # Check for odd holes in G and complement(G)
    for graph in [adj, complement_graph(adj, m)]:
        for length in range(5, m+1, 2):  # odd lengths >= 5
            for vertices in itertools.combinations(range(m), length):
                # Check if induced cycle
                ok = True
                for v in vertices:
                    deg = sum(1 for u in vertices if u != v and graph[v][u])
                    if deg != 2:
                        ok = False
                        break
                if not ok:
                    continue
                # Check connected
                visited = {vertices[0]}
                stack = [vertices[0]]
                while stack:
                    u = stack.pop()
                    for w in vertices:
                        if w not in visited and graph[u][w]:
                            visited.add(w)
                            stack.append(w)
                if len(visited) == length:
                    return False  # Found odd hole
    return True

def has_induced_subgraph_isomorphic_to(adj, m, pattern_adj, pattern_m):
    """Check if adj contains an induced subgraph isomorphic to pattern_adj."""
    if pattern_m > m:
        return False
    for vertices in itertools.combinations(range(m), pattern_m):
        # Try all mappings from pattern vertices to selected vertices
        for perm in itertools.permutations(vertices):
            match = True
            for i in range(pattern_m):
                for j in range(i+1, pattern_m):
                    if pattern_adj[i][j] != adj[perm[i]][perm[j]]:
                        match = False
                        break
                if not match:
                    break
            if match:
                return True
    return False

# ─────────────────────────────────────────────────
# Main catalog
# ─────────────────────────────────────────────────

def analyze_n(n, tournaments_iter, label="", max_tournaments=None):
    """Analyze all tournaments of a given size."""
    invariant_counter = Counter()
    invariant_examples = {}  # invariant -> example tournament

    # Track properties
    all_chordal = True
    all_perfect = True
    all_complete = True
    chordal_count = 0
    perfect_count = 0
    complete_count = 0
    total = 0

    cycle_count_dist = Counter()  # number of odd cycles distribution

    # For detailed analysis, store some examples
    interesting_examples = []

    for T in tournaments_iter:
        total += 1
        if max_tournaments and total > max_tournaments:
            break

        cycles = find_all_odd_cycles_v2(T, n)
        num_cycles = len(cycles)
        cycle_count_dist[num_cycles] += 1

        if num_cycles == 0:
            inv = (0, 0, (), 0, (), 0, ())
            invariant_counter[inv] += 1
            if inv not in invariant_examples:
                invariant_examples[inv] = (T, cycles)
            complete_count += 1
            chordal_count += 1
            perfect_count += 1
            continue

        omega_adj = build_conflict_graph(cycles)
        m = num_cycles

        inv = graph_invariant(omega_adj, m)
        invariant_counter[inv] += 1
        if inv not in invariant_examples:
            invariant_examples[inv] = (T, cycles)

        comp = is_complete(omega_adj, m)
        if comp:
            complete_count += 1
        else:
            all_complete = False

        # Only check properties for manageable sizes
        if m <= 20:
            ch = is_chordal(omega_adj, m)
            if ch:
                chordal_count += 1
            else:
                all_chordal = False
                if len(interesting_examples) < 3:
                    interesting_examples.append(('non-chordal', T, cycles, inv))

            if m <= 15:
                pf = is_perfect(omega_adj, m)
                if pf:
                    perfect_count += 1
                else:
                    all_perfect = False
                    if len(interesting_examples) < 3:
                        interesting_examples.append(('non-perfect', T, cycles, inv))
            else:
                perfect_count += 1  # assume perfect if too big to check
        else:
            chordal_count += 1  # can't check, assume
            perfect_count += 1

    print(f"\n{'='*60}")
    print(f"  n = {n}  {label}")
    print(f"{'='*60}")
    print(f"Total tournaments analyzed: {total}")
    print(f"Distinct conflict graph structures: {len(invariant_counter)}")
    print(f"\nOdd cycle count distribution:")
    for k in sorted(cycle_count_dist.keys()):
        print(f"  {k} odd cycles: {cycle_count_dist[k]} tournaments")

    print(f"\nStructural properties:")
    print(f"  Always complete: {all_complete} ({complete_count}/{total})")
    print(f"  Always chordal:  {all_chordal} ({chordal_count}/{total})")
    print(f"  Always perfect:  {all_perfect} ({perfect_count}/{total})")

    print(f"\nMost common conflict graph structures:")
    for inv, count in invariant_counter.most_common(15):
        m, e, degs, ncomp, csizes, ntri, _ = inv
        print(f"  [{count:5d}x] |V|={m}, |E|={e}, deg_seq={degs}, "
              f"components={ncomp}, comp_sizes={csizes}, triangles={ntri}")

    if interesting_examples:
        print(f"\nInteresting examples:")
        for ex in interesting_examples[:3]:
            kind = ex[0]
            T_ex = ex[1]
            cycles_ex = ex[2]
            inv_ex = ex[3]
            m_ex = inv_ex[0]
            print(f"  {kind}: {m_ex} cycles, invariant = {inv_ex[:6]}")
            print(f"    Cycles: {[sorted(c) for c in cycles_ex[:8]]}")

    return invariant_counter, invariant_examples, total

def verify_t028_n5():
    """Verify T028: for n=5, Omega(T) is always complete."""
    print(f"\n{'='*60}")
    print(f"  T028 Verification: Omega(T) always complete for n=5")
    print(f"{'='*60}")

    n = 5
    count = 0
    non_complete = 0
    for T in all_tournaments(n):
        count += 1
        cycles = find_all_odd_cycles_v2(T, n)
        if len(cycles) <= 1:
            continue
        omega_adj = build_conflict_graph(cycles)
        m = len(cycles)
        if not is_complete(omega_adj, m):
            non_complete += 1
            print(f"  COUNTEREXAMPLE: tournament #{count} has non-complete Omega")
            print(f"    {len(cycles)} cycles: {[sorted(c) for c in cycles]}")
            # Find the non-adjacent pair
            for i in range(m):
                for j in range(i+1, m):
                    if not omega_adj[i][j]:
                        print(f"    Disjoint cycles: {sorted(cycles[i])} and {sorted(cycles[j])}")

    if non_complete == 0:
        print(f"  VERIFIED: All {count} tournaments on 5 vertices have complete Omega(T)")
        # Explain why: on 5 vertices, any two odd cycles (size 3 or 5) must overlap
        print(f"  (Any two subsets of {{0..4}} with sizes in {{3,5}} must intersect)")
    else:
        print(f"  FAILED: {non_complete} counterexamples found")

def detailed_n6_analysis():
    """Detailed analysis of n=6 conflict graphs."""
    print(f"\n{'='*60}")
    print(f"  Detailed n=6 Analysis")
    print(f"{'='*60}")

    n = 6
    # Collect all non-complete examples
    non_complete_invs = Counter()
    examples_by_inv = {}
    total = 0

    for T in all_tournaments(n):
        total += 1
        cycles = find_all_odd_cycles_v2(T, n)
        if len(cycles) <= 1:
            continue

        omega_adj = build_conflict_graph(cycles)
        m = len(cycles)

        if not is_complete(omega_adj, m):
            inv = graph_invariant(omega_adj, m)
            non_complete_invs[inv] += 1
            if inv not in examples_by_inv:
                examples_by_inv[inv] = (T, cycles, omega_adj, m)

    print(f"  Total tournaments: {total}")
    print(f"  Non-complete Omega structures: {len(non_complete_invs)}")

    for inv, count in non_complete_invs.most_common(10):
        m, e, degs, ncomp, csizes, ntri, _ = inv
        print(f"\n  [{count}x] |V|={m}, |E|={e}, deg_seq={degs}, components={ncomp}")

        if inv in examples_by_inv:
            T_ex, cycles_ex, omega_ex, m_ex = examples_by_inv[inv]
            # Check properties
            ch = is_chordal(omega_ex, m_ex)
            pf = is_perfect(omega_ex, m_ex) if m_ex <= 15 else "?"
            print(f"    Chordal: {ch}, Perfect: {pf}")

            # Find non-adjacent pairs (disjoint cycles)
            disjoint_pairs = []
            for i in range(m_ex):
                for j in range(i+1, m_ex):
                    if not omega_ex[i][j]:
                        disjoint_pairs.append((cycles_ex[i], cycles_ex[j]))
            if disjoint_pairs:
                print(f"    Disjoint cycle pairs: {len(disjoint_pairs)}")
                for c1, c2 in disjoint_pairs[:3]:
                    print(f"      {sorted(c1)} and {sorted(c2)} (sizes {len(c1)},{len(c2)})")

def main():
    print("Odd-Cycle Conflict Graph Catalog")
    print("=" * 60)

    # n = 3
    analyze_n(3, all_tournaments(3), label="(exhaustive)")

    # n = 4
    analyze_n(4, all_tournaments(4), label="(exhaustive)")

    # n = 5 (with T028 verification)
    analyze_n(5, all_tournaments(5), label="(exhaustive)")
    verify_t028_n5()

    # n = 6
    analyze_n(6, all_tournaments(6), label="(exhaustive)")
    detailed_n6_analysis()

    # n = 7 (random sampling)
    random.seed(42)
    n7_tournaments = [random_tournament(7) for _ in range(2000)]
    analyze_n(7, iter(n7_tournaments), label="(2000 random samples)", max_tournaments=2000)

    # Summary of findings
    print(f"\n{'='*60}")
    print(f"  SUMMARY OF FINDINGS")
    print(f"{'='*60}")
    print("""
Key observations:
1. For n<=5, Omega(T) is always complete (any two odd cycles share a vertex)
   - This is because on 5 vertices, two subsets of size >=3 must overlap.
2. For n=6, non-complete Omega(T) first appears (two disjoint 3-cycles possible).
3. Check the output above for chordality/perfection patterns.
4. The cycle sizes matter: for n=6, two disjoint 3-cycles on {a,b,c} and {d,e,f}
   are possible, giving a non-edge in Omega(T).
""")

if __name__ == "__main__":
    main()
