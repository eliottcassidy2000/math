#!/usr/bin/env python3
"""
opus-2026-04-05-S25: Explore deletion-contraction for tournaments.

The Mitrovic NC Rédei-Berge function satisfies W_X = W_{X\e} - W_{X/e}^up
for any arc e. In commutative variables, this becomes:
  U_D = U_{D\e} - U_{D/e}  (with some caveats about "up" operation)

For tournament counting (H = number of Hamiltonian paths):
  Deletion: D\e = remove arc e (may create a non-tournament!)
  Contraction: D/e = contract arc e (identify head and tail)

But for TOURNAMENTS, deleting an arc breaks the tournament property.
So the standard graph-theoretic DC doesn't directly apply.

HOWEVER: we can use the TILING MODEL instead.
In the tiling model, flipping a tile (arc not on the base path) gives
another valid tournament. This is NOT deletion-contraction, but it IS
a local operation with algebraic structure.

ALTERNATIVE: Use the OCF directly.
H(T) = I(Ω(T), 2) = Σ_{S indep in Ω(T)} 2^|S|
This satisfies: I(G, x) = I(G\v, x) + x·I(G-N[v], x)
(deletion-contraction for independence polynomials)

So we have a DC recursion on the CONFLICT GRAPH Ω(T)!
Pick any cycle C in Ω(T):
  I(Ω, x) = I(Ω\C, x) + x·I(Ω - N[C], x)

where Ω\C = delete vertex C from Ω, Ω-N[C] = delete C and all its neighbors.

This gives: H(T) = I(Ω\C, 2) + 2·I(Ω-N[C], 2)
            = (H(T) minus 2× contribution of cycles through C's vertices)

Let me implement and explore this recursion.
"""

from itertools import combinations, permutations
from collections import defaultdict
import time

def all_tournaments(n):
    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)
    for bits in range(2**m):
        T = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(pairs):
            if (bits >> k) & 1: T[i][j] = 1
            else: T[j][i] = 1
        yield T

def hamiltonian_path_count(T):
    n = len(T)
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                if T[v][u]: dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[(1 << n) - 1])

def find_odd_cycles(T):
    n = len(T)
    cycles = set()
    for L in range(3, n+1, 2):
        for verts in combinations(range(n), L):
            v0 = verts[0]
            for perm in permutations(verts[1:]):
                path = (v0,) + perm
                ok = True
                for i in range(L):
                    if T[path[i]][path[(i+1) % L]] != 1:
                        ok = False; break
                if ok:
                    mi = path.index(min(path))
                    canon = tuple(path[(mi+j) % L] for j in range(L))
                    cycles.add(canon)
                    break
    return list(cycles)

def build_conflict_graph(cycles):
    k = len(cycles)
    adj = [[False]*k for _ in range(k)]
    vsets = [set(c) for c in cycles]
    for i in range(k):
        for j in range(i+1, k):
            if vsets[i] & vsets[j]:
                adj[i][j] = adj[j][i] = True
    return adj, vsets

def independence_poly_eval(adj, x):
    """Evaluate I(G, x) by enumeration."""
    k = len(adj)
    total = 0
    for mask in range(2**k):
        verts = [i for i in range(k) if (mask >> i) & 1]
        ok = True
        for a in range(len(verts)):
            for b in range(a+1, len(verts)):
                if adj[verts[a]][verts[b]]:
                    ok = False; break
            if not ok: break
        if ok:
            total += x ** len(verts)
    return total

def delete_vertex(adj, v):
    """Delete vertex v from the graph."""
    k = len(adj)
    new_adj = []
    for i in range(k):
        if i == v: continue
        row = []
        for j in range(k):
            if j == v: continue
            row.append(adj[i][j])
        new_adj.append(row)
    return new_adj

def delete_closed_neighborhood(adj, v):
    """Delete vertex v and all its neighbors (closed neighborhood)."""
    k = len(adj)
    to_remove = {v} | {j for j in range(k) if adj[v][j]}
    new_adj = []
    for i in range(k):
        if i in to_remove: continue
        row = []
        for j in range(k):
            if j in to_remove: continue
            row.append(adj[i][j])
        new_adj.append(row)
    return new_adj

# ═══════════════════════════════════════════════════════
# EXPERIMENT 1: Independence polynomial DC on Ω(T)
# ═══════════════════════════════════════════════════════
print("="*70)
print("EXPERIMENT 1: DC recursion I(Ω,x) = I(Ω\\C,x) + x·I(Ω-N[C],x)")
print("="*70)

for n in [5, 6]:
    print(f"\nn={n}:")
    count = 0
    for T in all_tournaments(n):
        count += 1
        if count > 5: break  # just a few examples

        H = hamiltonian_path_count(T)
        cycles = find_odd_cycles(T)
        if len(cycles) == 0:
            print(f"  H={H}: transitive (no cycles)")
            continue

        adj, vsets = build_conflict_graph(cycles)

        # Verify I(Ω, 2) = H
        I_2 = independence_poly_eval(adj, 2)
        assert I_2 == H, f"OCF mismatch: I={I_2}, H={H}"

        # Pick the first cycle as the DC pivot
        v = 0  # vertex index in conflict graph
        cycle_v = cycles[v]

        # DC recursion
        adj_del = delete_vertex(adj, v)
        adj_nb = delete_closed_neighborhood(adj, v)

        I_del = independence_poly_eval(adj_del, 2)
        I_nb = independence_poly_eval(adj_nb, 2)

        # Verify: I(Ω, 2) = I(Ω\v, 2) + 2·I(Ω-N[v], 2)
        reconstructed = I_del + 2 * I_nb
        ok = "✓" if reconstructed == H else "✗"

        print(f"  H={H:3d}: |Ω|={len(cycles)}, pivot={cycle_v}, "
              f"I(Ω\\v)={I_del}, I(Ω-N[v])={I_nb}, "
              f"reconstructed={reconstructed} {ok}")

# ═══════════════════════════════════════════════════════
# EXPERIMENT 2: DC tree depth — how deep until base case?
# ═══════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("EXPERIMENT 2: DC tree depth (recursive until empty graph)")
print("="*70)

def dc_depth(adj, memo=None):
    """Compute I(G, 2) using DC recursion, tracking depth."""
    if memo is None:
        memo = {}

    k = len(adj)
    if k == 0:
        return 1, 0

    # Convert to hashable key
    key = tuple(tuple(row) for row in adj)
    if key in memo:
        return memo[key]

    # Pick vertex with highest degree (to reduce graph fastest)
    degrees = [sum(1 for j in range(k) if adj[i][j]) for i in range(k)]
    v = degrees.index(max(degrees))

    adj_del = delete_vertex(adj, v)
    adj_nb = delete_closed_neighborhood(adj, v)

    I_del, d1 = dc_depth(adj_del, memo)
    I_nb, d2 = dc_depth(adj_nb, memo)

    result = I_del + 2 * I_nb
    depth = max(d1, d2) + 1

    memo[key] = (result, depth)
    return result, depth

for n in [5, 6]:
    print(f"\nn={n}:")
    depths = []
    for T in all_tournaments(n):
        cycles = find_odd_cycles(T)
        if len(cycles) == 0:
            depths.append(0)
            continue
        adj, _ = build_conflict_graph(cycles)
        _, depth = dc_depth(adj)
        depths.append(depth)

    from collections import Counter
    depth_dist = Counter(depths)
    print(f"  DC tree depth distribution: {dict(sorted(depth_dist.items()))}")
    print(f"  Max depth: {max(depths)}, Mean: {sum(depths)/len(depths):.2f}")

# ═══════════════════════════════════════════════════════
# EXPERIMENT 3: DC as a tournament factorization
# ═══════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("EXPERIMENT 3: What does DC reveal about tournament structure?")
print("="*70)

# At each DC step, we split H(T) into:
#   H = (cycles not using C's vertices) + 2×(disjoint from C's neighborhood)
#
# The "2×" comes from the independence polynomial evaluation at x=2.
# In terms of tournament paths: adding cycle C as an independent member
# DOUBLES the contribution of the disjoint part.

# Interesting: the DC tree encodes a HIERARCHICAL decomposition of H.
# Let me trace through a specific example.

n = 6
print(f"\nDetailed DC trace for a specific tournament at n={n}:")

# Use the regular tournament (max H)
T = [[0]*n for _ in range(n)]
# Regular tournament on 6 vertices: score (2,2,2,3,3,3)
# Use a specific construction: cyclic tournament C_6 mod 3
for i in range(n):
    for j in range(n):
        if i != j:
            if (j - i) % n in [1, 2, 3]:  # beat the next 3
                T[i][j] = 1

H = hamiltonian_path_count(T)
scores = tuple(sorted(sum(T[i]) for i in range(n)))
print(f"  Tournament: scores={scores}, H={H}")

cycles = find_odd_cycles(T)
print(f"  Odd cycles: {len(cycles)}")
for c in sorted(cycles, key=lambda x: (len(x), x)):
    print(f"    {c} (len {len(c)})")

adj, vsets = build_conflict_graph(cycles)
print(f"  Conflict graph: {len(adj)} vertices, "
      f"{sum(1 for i in range(len(adj)) for j in range(i+1,len(adj)) if adj[i][j])} edges")

# DC trace
def dc_trace(adj, cycles, indent=0):
    k = len(adj)
    prefix = "  " * indent

    if k == 0:
        print(f"{prefix}→ base case: I = 1")
        return 1

    degrees = [sum(1 for j in range(k) if adj[i][j]) for i in range(k)]
    v = degrees.index(max(degrees))
    c = cycles[v] if v < len(cycles) else f"vertex_{v}"

    print(f"{prefix}DC on cycle {c} (degree={degrees[v]}/{k-1}):")

    # Delete vertex
    remaining_del = [cycles[i] for i in range(k) if i != v]
    adj_del = delete_vertex(adj, v)
    print(f"{prefix}  Ω\\v ({k-1} vertices):")
    I_del = dc_trace(adj_del, remaining_del, indent + 2)

    # Delete closed neighborhood
    removed = {v} | {j for j in range(k) if adj[v][j]}
    remaining_nb = [cycles[i] for i in range(k) if i not in removed]
    adj_nb = delete_closed_neighborhood(adj, v)
    print(f"{prefix}  Ω-N[v] ({len(adj_nb)} vertices, removed {len(removed)}):")
    I_nb = dc_trace(adj_nb, remaining_nb, indent + 2)

    result = I_del + 2 * I_nb
    print(f"{prefix}  = {I_del} + 2×{I_nb} = {result}")
    return result

if len(adj) <= 20:
    print(f"\n  DC trace:")
    result = dc_trace(adj, cycles, indent=2)
    print(f"\n  Final: I(Ω, 2) = {result}, H = {H}")
    assert result == H

print("\n\nDONE.")
