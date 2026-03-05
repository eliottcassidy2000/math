#!/usr/bin/env python3
"""
Tournament Computation Library
================================
Correct implementations of all core mathematical objects from definitions.md.

MISTAKE-001 compliance:
    mu(C) computation ALWAYS operates on T-v's adjacency matrix, never the full T.
    This is enforced by the delete_vertex -> find_odd_cycles pipeline.

Usage:
    from tournament_lib import *

    T = random_tournament(6)
    print(hamiltonian_path_count(T))          # H(T)
    ok, lhs, rhs = verify_claim_a(T, 0)      # check Claim A for vertex 0
    self_test()                                # run sanity checks

Representation:
    A tournament on n vertices {0, ..., n-1} is a list of n lists (adjacency matrix).
    T[i][j] = 1 means arc i -> j.  T[i][j] + T[j][i] = 1 for i != j.
"""

from itertools import permutations, combinations


# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------

def tournament_from_bits(n, bits):
    """Construct tournament from integer encoding of upper triangle.
    Bit k (k-th least significant) corresponds to pair (i,j) with i<j
    in lexicographic order.  bit=1 means i->j, bit=0 means j->i."""
    T = [[0]*n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << k):
                T[i][j] = 1
            else:
                T[j][i] = 1
            k += 1
    return T


def all_tournaments(n):
    """Generate all tournaments on n vertices."""
    m = n * (n - 1) // 2
    for bits in range(1 << m):
        yield tournament_from_bits(n, bits)


def random_tournament(n, rng=None):
    """Generate a uniformly random tournament on n vertices."""
    import random
    if rng is None:
        rng = random
    T = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                T[i][j] = 1
            else:
                T[j][i] = 1
    return T


def opposite_tournament(T):
    """Return T^op (all arcs reversed)."""
    n = len(T)
    return [[T[j][i] for j in range(n)] for i in range(n)]


# ---------------------------------------------------------------------------
# Vertex deletion
# ---------------------------------------------------------------------------

def delete_vertex(T, v):
    """Remove vertex v from tournament T.
    Returns (T_minus_v, old_labels) where old_labels[new_idx] = old_idx.

    IMPORTANT (MISTAKE-001): Always use this to get T-v.
    Never index into the original T when computing on T-v."""
    n = len(T)
    old_labels = [i for i in range(n) if i != v]
    new_n = n - 1
    Tv = [[0]*new_n for _ in range(new_n)]
    for ni in range(new_n):
        for nj in range(new_n):
            Tv[ni][nj] = T[old_labels[ni]][old_labels[nj]]
    return Tv, old_labels


# ---------------------------------------------------------------------------
# Hamiltonian paths
# ---------------------------------------------------------------------------

def hamiltonian_path_count(T):
    """Count directed Hamiltonian paths in tournament T.  (= H(T))
    Bitmask DP, O(2^n * n^2).  Practical for n <= ~20."""
    n = len(T)
    if n == 0:
        return 0
    if n == 1:
        return 1

    full = (1 << n) - 1
    # dp[mask][v] = paths through exactly the vertices in mask, ending at v
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1

    for mask in range(1, 1 << n):
        for v in range(n):
            cnt = dp[mask][v]
            if not (mask & (1 << v)) or cnt == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if T[v][u]:
                    dp[mask | (1 << u)][u] += cnt

    return sum(dp[full][v] for v in range(n))


# ---------------------------------------------------------------------------
# Odd cycles
# ---------------------------------------------------------------------------

def find_odd_cycles(T):
    """Find all directed odd cycles in tournament T.
    Returns list of tuples of vertices.  Each cycle is canonical:
    minimum vertex is placed first (eliminates rotational duplicates)."""
    n = len(T)
    cycles = []
    for length in range(3, n + 1, 2):
        for verts in combinations(range(n), length):
            first = verts[0]  # minimum, fixed as start
            for perm in permutations(verts[1:]):
                path = (first,) + perm
                valid = True
                for i in range(length):
                    if not T[path[i]][path[(i + 1) % length]]:
                        valid = False
                        break
                if valid:
                    cycles.append(path)
    return cycles


# ---------------------------------------------------------------------------
# Conflict graph and independence polynomial
# ---------------------------------------------------------------------------

def conflict_graph(cycles):
    """Build conflict graph on a list of directed cycles.
    Two cycles are adjacent iff they share at least one vertex.
    Returns adjacency matrix (list of lists)."""
    m = len(cycles)
    adj = [[0]*m for _ in range(m)]
    vsets = [set(c) for c in cycles]
    for i in range(m):
        for j in range(i + 1, m):
            if vsets[i] & vsets[j]:
                adj[i][j] = 1
                adj[j][i] = 1
    return adj


def independence_poly_at(adj, x):
    """Evaluate independence polynomial I(G, x) by enumeration.
    adj is an adjacency matrix (list of lists).
    O(2^m * m) where m = number of vertices."""
    m = len(adj)
    if m == 0:
        return 1  # empty graph: I = 1 (just the empty set)

    # Precompute neighbor bitmasks for faster independence checks
    nbr = [0] * m
    for i in range(m):
        for j in range(m):
            if adj[i][j]:
                nbr[i] |= 1 << j

    total = 0
    for mask in range(1 << m):
        # Check independence: no two selected vertices are adjacent
        ok = True
        seen = 0
        temp = mask
        while temp:
            v = (temp & -temp).bit_length() - 1  # lowest set bit
            if nbr[v] & seen:
                ok = False
                break
            seen |= 1 << v
            temp &= temp - 1
        if ok:
            total += x ** bin(mask).count('1')
    return total


# ---------------------------------------------------------------------------
# mu(C) -- the cycle weight
# ---------------------------------------------------------------------------

def mu(T, v, cycle, _tv_cache=None):
    """Compute mu(C) for directed odd cycle C through vertex v in tournament T.

    mu(C) = I(Omega(T-v)|_{avoid C\\{v}}, 2)

    Per definitions.md and MISTAKE-001:
      1. Form T-v
      2. Find all odd cycles of T-v
      3. Keep only cycles vertex-disjoint from C\\{v}
      4. Build conflict graph on those
      5. Return I(conflict_graph, 2)

    _tv_cache: optional (Tv, old_labels, tv_cycles) to avoid recomputation
               when calling mu multiple times for the same (T, v).
    """
    if _tv_cache is not None:
        Tv, old_labels, tv_cycles = _tv_cache
    else:
        Tv, old_labels = delete_vertex(T, v)
        tv_cycles = find_odd_cycles(Tv)

    # Map C\{v} to T-v labels
    old_to_new = {old: new for new, old in enumerate(old_labels)}
    avoid = {old_to_new[u] for u in cycle if u != v}

    # Filter to cycles vertex-disjoint from avoid set
    filtered = [c for c in tv_cycles if not (set(c) & avoid)]

    if not filtered:
        return 1  # I(empty graph, 2) = 1

    cg = conflict_graph(filtered)
    return independence_poly_at(cg, 2)


# ---------------------------------------------------------------------------
# Claim verification
# ---------------------------------------------------------------------------

def verify_claim_a(T, v):
    """Verify Claim A: H(T) - H(T-v) = 2 * sum_{C through v} mu(C).
    Returns (passed, lhs, rhs)."""
    ht = hamiltonian_path_count(T)
    Tv, old_labels = delete_vertex(T, v)
    htv = hamiltonian_path_count(Tv)
    lhs = ht - htv

    all_cycles = find_odd_cycles(T)
    cycles_v = [c for c in all_cycles if v in set(c)]

    if not cycles_v:
        return lhs == 0, lhs, 0

    # Cache T-v data for mu calls
    tv_cycles = find_odd_cycles(Tv)
    cache = (Tv, old_labels, tv_cycles)

    total_mu = sum(mu(T, v, c, _tv_cache=cache) for c in cycles_v)
    rhs = 2 * total_mu
    return lhs == rhs, lhs, rhs


def verify_claim_b(T, v):
    """Verify Claim B: I(Omega(T),2) - I(Omega(T-v),2) = 2 * sum_{C through v} mu(C).
    Returns (passed, lhs, rhs)."""
    all_cycles_T = find_odd_cycles(T)
    cg_T = conflict_graph(all_cycles_T)
    i_T = independence_poly_at(cg_T, 2)

    Tv, old_labels = delete_vertex(T, v)
    tv_cycles = find_odd_cycles(Tv)
    cg_Tv = conflict_graph(tv_cycles)
    i_Tv = independence_poly_at(cg_Tv, 2)

    lhs = i_T - i_Tv

    cycles_v = [c for c in all_cycles_T if v in set(c)]
    cache = (Tv, old_labels, tv_cycles)
    total_mu = sum(mu(T, v, c, _tv_cache=cache) for c in cycles_v)
    rhs = 2 * total_mu
    return lhs == rhs, lhs, rhs


def verify_redei(T):
    """Verify Redei's theorem: H(T) is odd.  Returns (passed, H(T))."""
    h = hamiltonian_path_count(T)
    return h % 2 == 1, h


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

def self_test():
    """Quick sanity checks against known values."""
    # --- n=3 transitive: 0->1, 0->2, 1->2 ---
    T3t = [[0, 1, 1], [0, 0, 1], [0, 0, 0]]
    assert hamiltonian_path_count(T3t) == 1, "H(transitive T3) should be 1"
    assert find_odd_cycles(T3t) == [], "transitive T3 has no cycles"

    # --- n=3 cyclic: 0->1, 1->2, 2->0 ---
    T3c = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
    assert hamiltonian_path_count(T3c) == 3, "H(cyclic T3) should be 3"
    cyc3 = find_odd_cycles(T3c)
    assert len(cyc3) == 1, "cyclic T3 has exactly one 3-cycle"
    assert set(cyc3[0]) == {0, 1, 2}

    # --- Claim A for n=3 cyclic, each vertex ---
    for v in range(3):
        ok, lhs, rhs = verify_claim_a(T3c, v)
        assert ok, f"Claim A failed for cyclic T3, v={v}: {lhs} != {rhs}"
        assert lhs == 2 and rhs == 2

    # --- Claim A for n=3 transitive, each vertex ---
    for v in range(3):
        ok, lhs, rhs = verify_claim_a(T3t, v)
        assert ok, f"Claim A failed for transitive T3, v={v}: {lhs} != {rhs}"

    # --- Redei's theorem for all n=4 tournaments ---
    count4 = 0
    for T in all_tournaments(4):
        ok, h = verify_redei(T)
        assert ok, f"Redei failed: H(T)={h} is even"
        count4 += 1
    assert count4 == 64, f"Expected 64 tournaments on 4 vertices, got {count4}"

    # --- Claim A for all n=4 tournaments ---
    failures = 0
    pairs = 0
    for T in all_tournaments(4):
        for v in range(4):
            ok, lhs, rhs = verify_claim_a(T, v)
            pairs += 1
            if not ok:
                failures += 1
    assert failures == 0, f"Claim A: {failures}/{pairs} failures at n=4"
    assert pairs == 256

    # --- Claim B for all n=4 tournaments ---
    failures_b = 0
    for T in all_tournaments(4):
        for v in range(4):
            ok, _, _ = verify_claim_b(T, v)
            if not ok:
                failures_b += 1
    assert failures_b == 0, f"Claim B: {failures_b}/256 failures at n=4"

    # --- Independence polynomial sanity ---
    # Empty graph on 3 vertices: I = 1 + 3x + 3x^2 + x^3 = (1+x)^3
    adj_empty = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    assert independence_poly_at(adj_empty, 2) == 27  # 3^3
    assert independence_poly_at(adj_empty, 1) == 8   # 2^3

    # Complete graph K3: only singletons and empty set: I = 1 + 3x
    adj_k3 = [[0, 1, 1], [1, 0, 1], [1, 1, 0]]
    assert independence_poly_at(adj_k3, 2) == 7  # 1 + 3*2
    assert independence_poly_at(adj_k3, 1) == 4  # 1 + 3

    # I(empty graph, x) = 1
    assert independence_poly_at([], 2) == 1

    print("All self-tests passed.")


if __name__ == "__main__":
    self_test()
