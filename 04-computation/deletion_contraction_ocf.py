#!/usr/bin/env python3
"""
Investigate whether Mitrovic's deletion-contraction for the noncommuting
Rédei-Berge function can give an inductive proof of OCF.

The key identity (arXiv:2504.20968, Theorem 3.7):
    W_X = W_{X\e} - W_{X/e}↑
where e is an edge, X\e is edge deletion, X/e is edge contraction.

For tournaments T on n vertices with edge e=(u,v):
- T\e: digraph on n vertices with edge e removed (NOT a tournament)
- T/e: digraph on n-1 vertices with u,v merged

Question: How does OCF (H(T) = I(Omega(T), 2)) behave under deletion-contraction?

Specifically: if H(T) = I(Omega(T), 2) for all tournaments on < n vertices,
does the deletion-contraction decomposition force it for n vertices?

Note: T\e and T/e are NOT tournaments, so we need OCF for general digraphs
or a tournament-specific reduction.

kind-pasteur-2026-03-06-S19
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count
from itertools import combinations, permutations

def build_all_tournaments(n):
    """Generate all tournaments on n vertices."""
    m = n*(n-1)//2
    tournaments = []
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        tournaments.append(T)
    return tournaments

def count_ham_paths_digraph(adj, n):
    """Count Hamiltonian paths in a general digraph using DP."""
    # dp[mask][v] = number of Ham paths through vertices in mask, ending at v
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            c = dp.get((mask, v), 0)
            if c == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + c
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def delete_edge(T, u, v):
    """Delete edge (u,v) from tournament T. Result is NOT a tournament."""
    n = len(T)
    D = [row[:] for row in T]
    D[u][v] = 0
    return D

def contract_edge(T, u, v):
    """Contract edge (u,v): merge u and v into one vertex.
    The merged vertex w has:
    - outgoing edges to z if v->z (follow the head)
    - incoming edges from z if z->u (follow the tail)
    For z where u->z and v->z, w->z
    For z where z->u and z->v, z->w
    For z where u->z and z->v (or v->z and z->u): ambiguous in general

    For Mitrovic: contraction merges u,v into one vertex labeled u.
    (w,u) in E' iff (w,u) in E
    (u,w) in E' iff (v,w) in E
    """
    n = len(T)
    # New vertices: 0,...,n-1 minus v, relabeled
    # Vertex u becomes the merged vertex
    # We keep all vertices except v
    verts = [i for i in range(n) if i != v]
    m = len(verts)
    D = [[0]*m for _ in range(m)]
    for i_idx, i in enumerate(verts):
        for j_idx, j in enumerate(verts):
            if i_idx == j_idx:
                continue
            # If i is the merged vertex (was u)
            if i == u:
                # (merged, j): edge exists iff (v, j) in original
                D[i_idx][j_idx] = T[v][j]
            elif j == u:
                # (i, merged): edge exists iff (i, u) in original
                D[i_idx][j_idx] = T[i][u]
            else:
                D[i_idx][j_idx] = T[i][j]
    return D

def find_odd_cycles_digraph(adj, n):
    """Find all directed odd cycles in a digraph."""
    cycles = []
    for k in range(3, n+1, 2):
        for combo in combinations(range(n), k):
            verts = list(combo)
            v0 = verts[0]
            # DP: count directed Ham cycles on this subset
            dp = {}
            dp[(1, 0)] = 1
            for mask in range(1, 1 << k):
                if not (mask & 1):
                    continue
                for vi in range(k):
                    if not (mask & (1 << vi)):
                        continue
                    c = dp.get((mask, vi), 0)
                    if c == 0:
                        continue
                    for ui in range(k):
                        if mask & (1 << ui):
                            continue
                        if adj[verts[vi]][verts[ui]]:
                            key = (mask | (1 << ui), ui)
                            dp[key] = dp.get(key, 0) + c
            full = (1 << k) - 1
            num_cycles = 0
            for vi in range(1, k):
                c = dp.get((full, vi), 0)
                if c > 0 and adj[verts[vi]][verts[0]]:
                    num_cycles += c
            for _ in range(num_cycles):
                cycles.append(frozenset(combo))
    return cycles

def compute_ip_at_2(adj, n):
    """Compute I(Omega(digraph), 2) where Omega is the odd-cycle conflict graph."""
    # Find all directed odd cycles
    cycle_list = []
    for k in range(3, n+1, 2):
        for combo in combinations(range(n), k):
            verts = list(combo)
            v0 = verts[0]
            dp = {}
            dp[(1, 0)] = 1
            for mask in range(1, 1 << k):
                if not (mask & 1):
                    continue
                for vi in range(k):
                    if not (mask & (1 << vi)):
                        continue
                    c = dp.get((mask, vi), 0)
                    if c == 0:
                        continue
                    for ui in range(k):
                        if mask & (1 << ui):
                            continue
                        if adj[verts[vi]][verts[ui]]:
                            key = (mask | (1 << ui), ui)
                            dp[key] = dp.get(key, 0) + c
            full = (1 << k) - 1
            num_dir_cycles = 0
            for vi in range(1, k):
                c = dp.get((full, vi), 0)
                if c > 0 and adj[verts[vi]][verts[0]]:
                    num_dir_cycles += c
            # Each directed cycle is a vertex of Omega
            for _ in range(num_dir_cycles):
                cycle_list.append(frozenset(combo))

    if not cycle_list:
        return 1  # I(empty, 2) = 1

    # Build Omega: two cycles adjacent iff they share a vertex
    m = len(cycle_list)
    adj_omega = [[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            if cycle_list[i] & cycle_list[j]:
                adj_omega[i][j] = 1
                adj_omega[j][i] = 1

    # Compute I(Omega, 2) = sum over independent sets S of 2^|S|
    total = 0
    for mask in range(1 << m):
        # Check if mask is an independent set
        is_indep = True
        verts_in = [i for i in range(m) if mask & (1 << i)]
        for a in range(len(verts_in)):
            for b in range(a+1, len(verts_in)):
                if adj_omega[verts_in[a]][verts_in[b]]:
                    is_indep = False
                    break
            if not is_indep:
                break
        if is_indep:
            total += 2 ** len(verts_in)
    return total

# ============================================================
# Test: Does deletion-contraction preserve OCF?
# ============================================================
print("=" * 70)
print("DELETION-CONTRACTION AND OCF")
print("=" * 70)

# For each tournament T on n=4 vertices and each edge e:
# Check: H(T) = H(T\e) - H(T/e) (for Hamiltonian paths)
# Check: I(Omega(T), 2) = ??? for T\e and T/e

n = 4
m = n*(n-1)//2
print(f"\nn={n}: Testing deletion-contraction for Ham paths")

tested = 0
dc_works_h = 0
dc_works_ocf = 0

for bits in range(min(1 << m, 64)):  # Test first 64 tournaments
    T = tournament_from_bits(n, bits)
    H_T = hamiltonian_path_count(T)

    # Test each edge
    for u in range(n):
        for v in range(n):
            if u >= v or not T[u][v]:
                continue

            T_del = delete_edge(T, u, v)
            T_con = contract_edge(T, u, v)

            H_del = count_ham_paths_digraph(T_del, n)
            H_con = count_ham_paths_digraph(T_con, n-1)

            tested += 1

            # Check: H(T) = H(T\e) - H(T/e) ?
            if H_T == H_del - H_con:
                dc_works_h += 1

print(f"  Tested {tested} (tournament, edge) pairs")
print(f"  H(T) = H(T\\e) - H(T/e): {dc_works_h}/{tested} ({100*dc_works_h/tested:.1f}%)")

# The deletion-contraction is for W_X (noncommuting), not directly for H.
# Let's check if there's a simpler relationship.

print(f"\n  Checking H(T), H(T\\e), H(T/e) values:")
for bits in range(min(1 << m, 8)):
    T = tournament_from_bits(n, bits)
    H_T = hamiltonian_path_count(T)
    for u in range(n):
        for v in range(u+1, n):
            if not T[u][v]:
                continue
            T_del = delete_edge(T, u, v)
            T_con = contract_edge(T, u, v)
            H_del = count_ham_paths_digraph(T_del, n)
            H_con = count_ham_paths_digraph(T_con, n-1)
            print(f"    bits={bits}, e=({u},{v}): H(T)={H_T}, H(T\\e)={H_del}, H(T/e)={H_con}, "
                  f"del-con={H_del-H_con}, diff={H_T-(H_del-H_con)}")
            break  # Just first edge per tournament
        break

# Now test OCF for the non-tournament digraphs
print(f"\n\nChecking OCF for T\\e (non-tournament digraphs) at n={n}:")
ocf_holds_del = 0
ocf_fails_del = 0
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    for u in range(n):
        for v in range(u+1, n):
            if not T[u][v]:
                continue
            T_del = delete_edge(T, u, v)
            H_del = count_ham_paths_digraph(T_del, n)
            I_del = compute_ip_at_2(T_del, n)
            if H_del == I_del:
                ocf_holds_del += 1
            else:
                ocf_fails_del += 1
                if ocf_fails_del <= 3:
                    print(f"  FAIL: bits={bits}, e=({u},{v}): H(T\\e)={H_del}, I(Omega(T\\e),2)={I_del}")
            break
        break

total_del = ocf_holds_del + ocf_fails_del
print(f"  OCF for T\\e: {ocf_holds_del}/{total_del} hold ({100*ocf_holds_del/total_del:.1f}%)")

# Test OCF for T/e (contracted digraph, n-1 vertices)
print(f"\nChecking OCF for T/e (contracted digraphs) at n={n}->n-1={n-1}:")
ocf_holds_con = 0
ocf_fails_con = 0
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    for u in range(n):
        for v in range(u+1, n):
            if not T[u][v]:
                continue
            T_con = contract_edge(T, u, v)
            H_con = count_ham_paths_digraph(T_con, n-1)
            I_con = compute_ip_at_2(T_con, n-1)
            if H_con == I_con:
                ocf_holds_con += 1
            else:
                ocf_fails_con += 1
                if ocf_fails_con <= 3:
                    print(f"  FAIL: bits={bits}, e=({u},{v}): H(T/e)={H_con}, I(Omega(T/e),2)={I_con}")
            break
        break

total_con = ocf_holds_con + ocf_fails_con
print(f"  OCF for T/e: {ocf_holds_con}/{total_con} hold ({100*ocf_holds_con/total_con:.1f}%)")

print("\nDone.")
