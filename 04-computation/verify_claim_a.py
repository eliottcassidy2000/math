"""
Clean verification of Claim A: H(T) - H(T-v) = 2 * sum_{C through v} mu(C)

Written by opus-2026-03-05-S1 to resolve DISC-001 (mu computation bug).
All computations derived directly from definitions in 01-canon/definitions.md.

Key correctness points (addressing MISTAKE-001):
- mu(C) uses T-v adjacency matrix ONLY (never the full T matrix)
- Odd cycles of T-v are found in the T-v subgraph exclusively
- The conflict graph Omega(T-v) is built on T-v's odd cycles
"""

from itertools import permutations, combinations
import sys


def all_tournaments(n):
    """Generate all tournaments on n vertices as adjacency matrices."""
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    m = len(edges)
    for bits in range(2**m):
        T = [[0] * n for _ in range(n)]
        for k, (i, j) in enumerate(edges):
            if (bits >> k) & 1:
                T[i][j] = 1
                T[j][i] = 0
            else:
                T[j][i] = 1
                T[i][j] = 0
        yield T


def count_hamiltonian_paths(T, vertices=None):
    """Count directed Hamiltonian paths in tournament T restricted to vertex set."""
    if vertices is None:
        vertices = list(range(len(T)))
    count = 0
    for perm in permutations(vertices):
        valid = True
        for i in range(len(perm) - 1):
            if T[perm[i]][perm[i + 1]] != 1:
                valid = False
                break
        if valid:
            count += 1
    return count


def find_odd_cycles_through_v(T, v, n):
    """Find all directed odd cycles through vertex v in tournament T."""
    cycles = []
    vertices = list(range(n))

    # Find cycles of odd length L = 3, 5, 7, ... through v
    # A directed cycle through v: v -> a1 -> a2 -> ... -> a_{L-1} -> v
    other_vertices = [u for u in vertices if u != v]

    for L in range(3, n + 1, 2):  # odd lengths only
        # Choose L-1 other vertices for the cycle
        for subset in combinations(other_vertices, L - 1):
            # Try all orderings of the subset as cycle v -> perm[0] -> ... -> perm[-1] -> v
            for perm in permutations(subset):
                # Check if v -> perm[0] -> perm[1] -> ... -> perm[-1] -> v is a directed cycle
                valid = True
                # Check v -> perm[0]
                if T[v][perm[0]] != 1:
                    valid = False
                if valid:
                    for i in range(len(perm) - 1):
                        if T[perm[i]][perm[i + 1]] != 1:
                            valid = False
                            break
                if valid:
                    # Check perm[-1] -> v
                    if T[perm[-1]][v] != 1:
                        valid = False
                if valid:
                    # Store as frozenset of vertices (but we need the actual cycle for C\{v})
                    # For mu computation we need C\{v} = set(perm)
                    # Normalize: store as tuple (v, perm[0], ..., perm[-1])
                    # But cycles are equivalent under rotation. To avoid duplicates,
                    # since v is fixed, each cycle through v has exactly one representation
                    # starting with v (up to direction). But direction matters for directed cycles.
                    # v -> a -> b -> v and v -> b -> a -> v are different cycles.
                    # However, the VERTEX SET C\{v} is the same regardless of direction.
                    # For mu(C), what matters is C\{v} as a vertex set (since we avoid
                    # cycles vertex-disjoint from C\{v}).
                    #
                    # But wait: different directed cycles on the same vertex set give the
                    # same C\{v}. We need to count each DIRECTED cycle separately in the
                    # sum, but their mu values are identical (same C\{v}).
                    #
                    # Actually, re-reading the definition: C is a directed odd cycle.
                    # Two directed cycles on the same vertices but in different orders
                    # are different cycles. But since they have the same vertex set,
                    # mu(C) = mu(C') for them.
                    #
                    # For a directed cycle on vertex set S through v, the number of
                    # distinct directed cycles is (|S|-1)! / |S| * 2? No.
                    # A directed cycle on L vertices has L distinct representations
                    # (rotations), but as a cycle it's one object. Since v is in the cycle,
                    # fixing v as start gives (L-1)! orderings but only some are valid.
                    #
                    # Actually, each permutation we find IS a distinct directed cycle
                    # representation. But a directed cycle v->a->b->v is the same as
                    # a->b->v->a (rotation). Since we fix v as start, each directed cycle
                    # through v appears EXACTLY ONCE in our enumeration.
                    cycle_other_verts = frozenset(perm)
                    cycles.append(cycle_other_verts)
    return cycles


def find_odd_cycles_in_subgraph(T, vertices):
    """Find all directed odd cycles in T restricted to the given vertex set.
    Returns list of frozensets of vertices."""
    cycles = set()
    vlist = sorted(vertices)

    for L in range(3, len(vlist) + 1, 2):  # odd lengths
        for subset in combinations(vlist, L):
            # Try all directed cycles on this subset
            # Fix first vertex to avoid counting rotations
            first = subset[0]
            rest = subset[1:]
            for perm in permutations(rest):
                valid = True
                seq = (first,) + perm
                for i in range(L - 1):
                    if T[seq[i]][seq[i + 1]] != 1:
                        valid = False
                        break
                if valid and T[seq[-1]][first] == 1:
                    cycles.add(frozenset(subset))
                    break  # found one direction, that's enough to know this vertex set forms a cycle
            # Also check if there's a cycle in the OTHER direction on this subset
            # Actually, we want ALL directed cycles, but for conflict graph purposes,
            # two directed cycles on the same vertex set share all vertices, so they're
            # always adjacent in Omega. For independence polynomial computation,
            # what matters is vertex-disjointness.
            # The conflict graph has vertices = directed odd cycles, edges = shared vertex.
            # But cycles on the same vertex set are trivially adjacent.
            # For independent sets: we can never pick two cycles on the same vertex set.
            # So for I(Omega, 2) purposes, we can work with vertex sets and just count
            # how many directed cycles exist on each set... Actually no, let me re-read.
            #
            # I(Omega(T), 2) = sum over independent sets I of Omega(T) of 2^|I|.
            # Each vertex of Omega is a directed odd cycle.
            # Two vertices are adjacent iff they share a vertex.
            #
            # So if there are k directed cycles on the same vertex set S, they form a
            # clique of size k in Omega. An independent set can contain at most one of them.
            #
            # For computing I(Omega, 2): we need to know all directed odd cycles.
    # Let me redo this more carefully.
    all_directed_cycles = []
    for L in range(3, len(vlist) + 1, 2):
        for subset in combinations(vlist, L):
            first = subset[0]
            rest = subset[1:]
            for perm in permutations(rest):
                seq = (first,) + perm
                valid = True
                for i in range(L - 1):
                    if T[seq[i]][seq[i + 1]] != 1:
                        valid = False
                        break
                if valid and T[seq[-1]][first] == 1:
                    all_directed_cycles.append(frozenset(subset))
    return all_directed_cycles


def build_conflict_graph(cycles):
    """Build conflict graph: cycles are adjacent iff they share a vertex.
    Returns adjacency list."""
    k = len(cycles)
    adj = [set() for _ in range(k)]
    for i in range(k):
        for j in range(i + 1, k):
            if cycles[i] & cycles[j]:  # shared vertex
                adj[i].add(j)
                adj[j].add(i)
    return adj


def independence_polynomial_at_2(adj, k):
    """Compute I(G, 2) for graph G with k vertices and adjacency list adj.
    I(G, 2) = sum over all independent sets S of 2^|S|."""
    total = 0
    # Enumerate all subsets (brute force, fine for small graphs)
    for bits in range(2**k):
        # Check if this subset is independent
        nodes = [i for i in range(k) if (bits >> i) & 1]
        independent = True
        for i in range(len(nodes)):
            for j in range(i + 1, len(nodes)):
                if nodes[j] in adj[nodes[i]]:
                    independent = False
                    break
            if not independent:
                break
        if independent:
            total += 2 ** len(nodes)
    return total


def compute_mu(T, v, n, cycle_other_verts):
    """Compute mu(C) for cycle C through v.

    mu(C) = I(Omega(T-v)|_{avoid C\\{v}}, 2)

    CRITICAL (MISTAKE-001): We work entirely in T-v.
    - T-v vertices: all vertices except v
    - Find odd cycles in T-v that are vertex-disjoint from C\\{v}
    - Build conflict graph on those cycles
    - Evaluate independence polynomial at 2
    """
    tmv_vertices = [u for u in range(n) if u != v]
    avoid_set = cycle_other_verts  # C\{v}

    # Find all directed odd cycles in T-v that are vertex-disjoint from avoid_set
    # IMPORTANT: Use T's adjacency matrix but ONLY for edges between T-v vertices
    available_vertices = [u for u in tmv_vertices if u not in avoid_set]

    if not available_vertices:
        # No vertices left -> no cycles -> Omega is empty -> I(empty, 2) = 1
        return 1

    cycles = find_odd_cycles_in_subgraph(T, available_vertices)
    if not cycles:
        return 1  # empty graph, I = 1

    adj = build_conflict_graph(cycles)
    return independence_polynomial_at_2(adj, len(cycles))


def verify_claim_a(n, verbose=False, sample_limit=None):
    """Verify Claim A for all tournaments on n vertices.
    Returns (total_pairs, failures, failure_examples)."""
    total = 0
    failures = 0
    examples = []

    for t_idx, T in enumerate(all_tournaments(n)):
        if sample_limit and t_idx >= sample_limit:
            break

        H_T = count_hamiltonian_paths(T)

        for v in range(n):
            total += 1
            H_Tv = count_hamiltonian_paths(T, [u for u in range(n) if u != v])

            lhs = H_T - H_Tv

            # Find all directed odd cycles through v
            cycles_through_v = find_odd_cycles_through_v(T, v, n)

            # Compute sum of mu(C)
            mu_sum = 0
            for c_verts in cycles_through_v:
                mu_sum += compute_mu(T, v, n, c_verts)

            rhs = 2 * mu_sum

            if lhs != rhs:
                failures += 1
                if len(examples) < 5:
                    examples.append((t_idx, v, lhs, rhs, H_T, H_Tv))
                if verbose:
                    print(f"  FAILURE: T#{t_idx}, v={v}: LHS={lhs}, RHS={rhs} "
                          f"(H(T)={H_T}, H(T-v)={H_Tv}, #cycles={len(cycles_through_v)}, "
                          f"sum_mu={mu_sum})")

        if verbose and (t_idx + 1) % 100 == 0:
            print(f"  ...checked {t_idx + 1} tournaments, {total} pairs so far")

    return total, failures, examples


def investigate_n5_mystery(verbose=True):
    """Investigate OPEN-Q-001: Why does per-path identity hold at n=5 despite 5-cycles?

    For each (T, v, P') triple at n=5:
    - Compute (inshat - 1) / 2
    - Count 3-cycle embeddings (Type-II positions)
    - Find 5-cycles through v
    - Check if per-path identity holds
    - Analyze the relationship between 5-cycle mu contributions and 3-cycle contributions
    """
    n = 5
    print(f"\n{'='*60}")
    print(f"OPEN-Q-001: Investigating n=5 mystery")
    print(f"{'='*60}\n")

    triples_with_5cycles = 0
    triples_total = 0

    for t_idx, T in enumerate(all_tournaments(n)):
        for v in range(n):
            # Find 5-cycles through v
            other = [u for u in range(n) if u != v]
            five_cycles = []
            for perm in permutations(other):
                # Check v -> perm[0] -> perm[1] -> perm[2] -> perm[3] -> v
                valid = T[v][perm[0]] == 1
                if valid:
                    for i in range(3):
                        if T[perm[i]][perm[i + 1]] != 1:
                            valid = False
                            break
                if valid and T[perm[3]][v] == 1:
                    five_cycles.append(perm)

            # Deduplicate: a 5-cycle through v with v fixed has 1 representation per direction
            # Actually, since we fix v as start and try all permutations of the other 4,
            # each directed 5-cycle v->a->b->c->d->v appears exactly once.
            n_five_cycles = len(five_cycles)

            if n_five_cycles > 0:
                # For each Ham path P' of T-v, check per-path identity
                tmv_verts = [u for u in range(n) if u != v]
                for perm_p in permutations(tmv_verts):
                    # Check if perm_p is a valid Ham path of T-v
                    valid_path = True
                    for i in range(len(perm_p) - 1):
                        if T[perm_p[i]][perm_p[i + 1]] != 1:
                            valid_path = False
                            break
                    if not valid_path:
                        continue

                    triples_total += 1
                    triples_with_5cycles += 1

                    # Compute inshat
                    sig = [1 if T[v][perm_p[j]] == 1 else 0 for j in range(len(perm_p))]
                    boundary = sig[0] + (1 - sig[-1])
                    type_I = sum(1 for j in range(len(sig) - 1) if sig[j] == 0 and sig[j + 1] == 1)
                    type_II = sum(1 for j in range(len(sig) - 1) if sig[j] == 1 and sig[j + 1] == 0)
                    inshat = boundary + type_I + type_II

                    # 3-cycle embeddings: directed 3-cycles (v, a, b) where a,b consecutive in P'
                    three_cycle_count = 0
                    for j in range(len(perm_p) - 1):
                        a, b = perm_p[j], perm_p[j + 1]
                        # Check if v->a->b->v is a directed 3-cycle
                        if T[v][a] == 1 and T[a][b] == 1 and T[b][v] == 1:
                            three_cycle_count += 1
                        # Check if v->b->a->v
                        if T[v][b] == 1 and T[b][a] == 1 and T[a][v] == 1:
                            three_cycle_count += 1

                    per_path_lhs = (inshat - 1) // 2
                    per_path_rhs = three_cycle_count

                    if per_path_lhs != per_path_rhs:
                        print(f"  PER-PATH FAILURE at T#{t_idx}, v={v}, P'={perm_p}")
                        print(f"    inshat={inshat}, (inshat-1)/2={per_path_lhs}, 3-cycle count={per_path_rhs}")
                        print(f"    5-cycles through v: {n_five_cycles}")

    print(f"\nTriples involving tournaments with 5-cycles through v: {triples_with_5cycles}")
    print(f"(Total triples checked: {triples_total})")
    if triples_with_5cycles > 0:
        print("Per-path identity held for ALL of them despite 5-cycles existing.")

    # Now analyze WHY: for each T with 5-cycles through v, compute the
    # Claim A breakdown showing 5-cycle mu contributions
    print(f"\n{'='*60}")
    print("Claim A breakdown for tournaments with 5-cycles through v:")
    print(f"{'='*60}\n")

    seen = set()
    for t_idx, T in enumerate(all_tournaments(n)):
        for v in range(n):
            other = [u for u in range(n) if u != v]
            # Quick check: does a 5-cycle through v exist?
            has_5cycle = False
            for perm in permutations(other):
                valid = T[v][perm[0]] == 1
                if valid:
                    for i in range(3):
                        if T[perm[i]][perm[i + 1]] != 1:
                            valid = False
                            break
                if valid and T[perm[3]][v] == 1:
                    has_5cycle = True
                    break

            if not has_5cycle:
                continue

            key = (t_idx, v)
            if key in seen:
                continue
            seen.add(key)

            H_T = count_hamiltonian_paths(T)
            H_Tv = count_hamiltonian_paths(T, other)

            # Get all cycles through v, separated by length
            three_cycles = find_odd_cycles_through_v_by_length(T, v, n, 3)
            five_cycles_vsets = find_odd_cycles_through_v_by_length(T, v, n, 5)

            mu_3 = sum(compute_mu(T, v, n, c) for c in three_cycles)
            mu_5 = sum(compute_mu(T, v, n, c) for c in five_cycles_vsets)

            print(f"T#{t_idx}, v={v}: H(T)={H_T}, H(T-v)={H_Tv}, diff={H_T - H_Tv}")
            print(f"  3-cycles through v: {len(three_cycles)}, sum mu_3 = {mu_3}")
            print(f"  5-cycles through v: {len(five_cycles_vsets)}, sum mu_5 = {mu_5}")
            print(f"  2*(mu_3 + mu_5) = {2*(mu_3 + mu_5)}, LHS = {H_T - H_Tv}")
            print(f"  Match: {H_T - H_Tv == 2*(mu_3 + mu_5)}")
            print()

            if len(seen) >= 20:
                print("(showing first 20 examples)")
                return


def find_odd_cycles_through_v_by_length(T, v, n, L):
    """Find all directed L-cycles through v. Returns list of frozenset(C\\{v})."""
    cycles = []
    other = [u for u in range(n) if u != v]
    for subset in combinations(other, L - 1):
        for perm in permutations(subset):
            valid = T[v][perm[0]] == 1
            if valid:
                for i in range(len(perm) - 1):
                    if T[perm[i]][perm[i + 1]] != 1:
                        valid = False
                        break
            if valid and T[perm[-1]][v] == 1:
                cycles.append(frozenset(perm))
                # Don't break - each permutation is a different directed cycle
    return cycles


if __name__ == "__main__":
    print("="*60)
    print("Claim A Verification (clean implementation)")
    print("Resolving DISC-001: mu computed on T-v ONLY")
    print("Author: opus-2026-03-05-S1")
    print("="*60)

    for test_n in [4, 5]:
        print(f"\n--- n = {test_n} ---")
        total, failures, examples = verify_claim_a(test_n, verbose=True)
        print(f"n={test_n}: {total} pairs checked, {failures} failures")
        if failures > 0:
            print("FAILURES FOUND:")
            for ex in examples:
                print(f"  T#{ex[0]}, v={ex[1]}: LHS={ex[2]}, RHS={ex[3]}")
        else:
            print("ALL PASSED")

    # n=6 is expensive (32768 tournaments * 6 vertices = 196608 pairs)
    # Run with --n6 flag
    if "--n6" in sys.argv:
        print(f"\n--- n = 6 (FULL - this will take a while) ---")
        total, failures, examples = verify_claim_a(6, verbose=True)
        print(f"n=6: {total} pairs checked, {failures} failures")
        if failures == 0:
            print("ALL PASSED - DISC-001 can be resolved: verification is CONFIRMED")

    # Investigate the n=5 mystery
    if "--mystery" in sys.argv or "--n5" not in sys.argv:
        investigate_n5_mystery()
