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


def hamiltonian_path_count_ocf(T):
    """Count directed Hamiltonian paths via OCF: H(T) = I(Omega(T), 2).
    Finds all odd directed cycles, builds conflict graph, evaluates I at x=2.
    Faster than bitmask DP when cycle count is manageable (typical for n <= 14).
    For large n with many cycles, falls back to bitmask DP."""
    n = len(T)
    if n <= 1:
        return 1
    cycles = find_odd_cycles(T)
    if n <= 14:
        return independence_poly_at_fast(cycles, 2)
    else:
        cg = conflict_graph(cycles)
        return independence_poly_at(cg, 2)


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


def independence_poly_at_fast(cycles, x, max_indep_size=None):
    """Fast I(Omega, x) for conflict graphs of cycle collections.
    Exploits the fact that independent sets correspond to vertex-disjoint
    cycle collections. For n <= 11, max independent set size <= 3.
    For n <= 14, max size <= 4.  Much faster than brute-force 2^m."""
    m = len(cycles)
    if m == 0:
        return 1
    if max_indep_size is None:
        # Estimate: max VD odd cycles in n vertices is floor(n/3)
        total_verts = len(set(v for c in cycles for v in c))
        max_indep_size = total_verts // 3

    vsets = [frozenset(c) for c in cycles]
    # Build adjacency bitmasks
    adj_bits = [0] * m
    for a in range(m):
        for b in range(a + 1, m):
            if vsets[a] & vsets[b]:
                adj_bits[a] |= 1 << b
                adj_bits[b] |= 1 << a

    total = 1 + x * m  # size 0 and 1

    if max_indep_size >= 2:
        pairs = []
        for a in range(m):
            for b in range(a + 1, m):
                if not (adj_bits[a] & (1 << b)):
                    pairs.append((a, b))
                    total += x ** 2

        if max_indep_size >= 3:
            triples = []
            for a, b in pairs:
                for c in range(b + 1, m):
                    if not (adj_bits[a] & (1 << c)) and not (adj_bits[b] & (1 << c)):
                        triples.append((a, b, c))
                        total += x ** 3

            if max_indep_size >= 4:
                for a, b, c in triples:
                    for d in range(c + 1, m):
                        if (not (adj_bits[a] & (1 << d)) and
                            not (adj_bits[b] & (1 << d)) and
                            not (adj_bits[c] & (1 << d))):
                            total += x ** 4

    return total


def independence_poly_at(adj, x):
    """Evaluate independence polynomial I(G, x) by enumeration.
    adj is an adjacency matrix (list of lists).
    O(2^m * m) where m = number of vertices.  Use independence_poly_at_fast
    for conflict graphs of cycles when m is large."""
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
# Arc-flip operations (for THM-013 / OPEN-Q-009)
# ---------------------------------------------------------------------------

def flip_arc(T, i, j):
    """Return T' obtained by flipping arc i->j to j->i.
    Requires T[i][j] == 1.  Returns a new tournament."""
    assert T[i][j] == 1, f"No arc {i}->{j} to flip"
    n = len(T)
    Tp = [row[:] for row in T]
    Tp[i][j] = 0
    Tp[j][i] = 1
    return Tp


def adj_count(T, i, j):
    """Count Hamiltonian paths of T with i immediately before j.
    adj_T(i,j) = #{P in Ham(T) : P[k]=i, P[k+1]=j for some k}.

    Method: adj(i,j) = sum_S dp_fwd[S][i] * dp_fwd_rev[V\\S][j]
    where dp_fwd[S][i] = paths through S ending at i (forward DP on T)
    and dp_fwd_rev[R][j] = paths through R ending at j (forward DP on T^op)
                         = paths through R starting at j in T."""
    n = len(T)
    if not T[i][j]:
        return 0
    full = (1 << n) - 1

    # Forward DP on T: dp[mask][v] = paths through mask ending at v
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

    # Forward DP on T^op: dp_r[mask][v] = paths through mask ending at v in T^op
    #                                    = paths through mask starting at v in T
    dp_r = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp_r[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            cnt = dp_r[mask][v]
            if not (mask & (1 << v)) or cnt == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if T[u][v]:  # reversed: arc u->v in T means v->u in T^op
                    dp_r[mask | (1 << u)][u] += cnt

    # adj(i,j) = sum over S containing i (not j) of dp[S][i] * dp_r[V\S][j]
    total = 0
    for mask in range(1, 1 << n):
        if not (mask & (1 << i)):
            continue
        if mask & (1 << j):
            continue
        prefix_count = dp[mask][i]
        if prefix_count == 0:
            continue
        rest = full ^ mask
        if not (rest & (1 << j)):
            continue
        suffix_count = dp_r[rest][j]
        total += prefix_count * suffix_count
    return total


def compute_s_x(T, i, j, x):
    """s_x = 1 - T[x][i] - T[j][x] for x in V\\{i,j}."""
    return 1 - T[x][i] - T[j][x]


def count_directed_5_cycles_through_arc(T, i, j):
    """Count directed 5-cycles in T that use the arc i->j.
    A 5-cycle through arc i->j: i->j->a->b->c->i for some {a,b,c} subset of V\\{i,j}."""
    n = len(T)
    if not T[i][j]:
        return 0
    others = [x for x in range(n) if x != i and x != j]
    count = 0
    for a, b, c in permutations(others, 3):
        if T[j][a] and T[a][b] and T[b][c] and T[c][i]:
            count += 1
    return count


def verify_thm013(T, i, j):
    """Verify THM-013 adjacency identity for flipping arc i->j.

    Identity: adj_T(i,j) - adj_{T'}(j,i) = -2*sum_x s_x*H(B_x) + 2*(D5-C5)

    Returns (passed, lhs, rhs)."""
    n = len(T)
    assert T[i][j] == 1, f"No arc {i}->{j}"

    Tp = flip_arc(T, i, j)

    lhs = adj_count(T, i, j) - adj_count(Tp, j, i)

    others = [x for x in range(n) if x != i and x != j]

    # Compute sum_x s_x * H(B_x)
    sum_sH = 0
    for x in others:
        sx = compute_s_x(T, i, j, x)
        # B_x = V\{i,j,x}
        bx_verts = [w for w in range(n) if w != i and w != j and w != x]
        if not bx_verts:
            h_bx = 1 if len(bx_verts) == 0 else 0
        else:
            # Build sub-tournament on bx_verts
            bx_n = len(bx_verts)
            Bx = [[0]*bx_n for _ in range(bx_n)]
            for ai in range(bx_n):
                for bi in range(bx_n):
                    Bx[ai][bi] = T[bx_verts[ai]][bx_verts[bi]]
            h_bx = hamiltonian_path_count(Bx)
        sum_sH += sx * h_bx

    D5 = count_directed_5_cycles_through_arc(T, i, j)
    C5 = count_directed_5_cycles_through_arc(Tp, j, i)

    rhs = -2 * sum_sH + 2 * (D5 - C5)

    return lhs == rhs, lhs, rhs


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


def verify_ocf(T):
    """Verify OCF: H(T) = I(Omega(T), 2).  Returns (passed, H, I)."""
    h = hamiltonian_path_count(T)
    cycles = find_odd_cycles(T)
    n = len(T)
    if n <= 14:
        i_val = independence_poly_at_fast(cycles, 2)
    else:
        cg = conflict_graph(cycles)
        i_val = independence_poly_at(cg, 2)
    return h == i_val, h, i_val


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
