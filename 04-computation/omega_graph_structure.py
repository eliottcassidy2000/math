#!/usr/bin/env python3
"""
omega_graph_structure.py — Structure of the odd-cycle intersection graph Ω(T)

THE KEY GAP: We need to show that spectral concentration → cycle clustering
→ higher independence polynomial Z(Ω, 2) = H.

Approach: Compare structural properties of Ω(Interval) vs Ω(Paley)
to identify which graph property drives the higher Z.

Candidates:
1. Maximum degree Δ(Ω): Scott-Sokal bound Z ≤ Π(1+λ)^{1/(d_v+1)}
2. Block structure: Does Ω(Int) decompose into clique-like blocks?
3. Fractional independence: Hoffman/Lovász bounds
4. Edge density: sparser Ω → higher Z
5. Clique cover number: smaller cover → more efficient packing

Author: opus-2026-03-12-S64
"""

import numpy as np
from itertools import combinations
from collections import defaultdict
import time

def make_tournament(p, S):
    A = np.zeros((p, p), dtype=np.int8)
    for i in range(p):
        for s in S:
            A[i][(i + s) % p] = 1
    return A

def get_QR(p):
    return sorted(set(pow(a, 2, p) for a in range(1, p)) - {0})

def find_all_directed_odd_cycles(A, p, max_len=None):
    """Find all directed odd cycles as (length, vertex_frozenset) pairs.
    For circulant tournaments, enumerate cycles through vertex 0
    and generate all rotations."""
    if max_len is None:
        max_len = p
    n = p
    cycles_through_0 = []

    for L in range(3, max_len + 1, 2):
        # Find directed L-cycles through vertex 0 (L edges, L vertices)
        # Path: 0 → v₁ → v₂ → ... → v_{L-1} → 0
        # This requires L steps total.
        dp = {(1, 0): 1}  # bit 0 set, at vertex 0

        for step in range(1, L + 1):
            new_dp = {}
            for (mask, v), count in dp.items():
                for w in range(n):
                    if step == L:
                        # Last step: must return to 0
                        if w == 0 and A[v][0]:
                            key = (mask, 0)
                            new_dp[key] = new_dp.get(key, 0) + count
                    else:
                        # Intermediate step: visit unvisited vertex (not 0)
                        if w != 0 and not (mask & (1 << w)) and A[v][w]:
                            new_mask = mask | (1 << w)
                            key = (new_mask, w)
                            new_dp[key] = new_dp.get(key, 0) + count
            dp = new_dp

        for (mask, v), count in dp.items():
            if v == 0 and count > 0:
                vset = frozenset(i for i in range(n) if i == 0 or (mask & (1 << i)))
                # Each directed cycle through 0 is counted once per direction
                # For tournaments (exactly one direction per edge), each cycle
                # corresponds to a unique directed cycle.
                # Actually, each k-cycle through 0 is found once because we
                # fixed the starting vertex as 0.
                cycles_through_0.append((L, vset, count))

    # For circulant tournaments on Z_p, cycles are invariant under rotation.
    # A cycle with vertex set V through 0 generates p cycles with vertex sets
    # {v+r mod p : v ∈ V} for r=0,...,p-1, but many are the same cycle.
    # Each cycle of length k visits k vertices, so the orbit has size p/gcd(p,orbit_stuff).
    # Since p is prime, each cycle has orbit size p (unless the cycle has
    # some special symmetry, which is rare).

    # For our purpose (vertex-disjointness), we need the vertex sets.
    # Generate all rotations of cycles through 0.
    all_cycle_vsets = set()
    for k, vset, count in cycles_through_0:
        for r in range(p):
            rotated = frozenset((v + r) % p for v in vset)
            all_cycle_vsets.add(rotated)

    return list(all_cycle_vsets)

def build_omega_graph(cycle_vsets):
    """Build the odd-cycle intersection graph: vertices = cycles,
    edges = pairs sharing a vertex."""
    n = len(cycle_vsets)
    adj = [[False]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if cycle_vsets[i] & cycle_vsets[j]:
                adj[i][j] = adj[j][i] = True
    return adj

def count_independent_sets(adj, n, max_check=None):
    """Count independent sets by size. Brute force for small n."""
    if max_check is None:
        max_check = n
    if n > 25:
        return None
    alpha = defaultdict(int)
    alpha[0] = 1
    for mask in range(1, 1 << n):
        k = bin(mask).count('1')
        if k > max_check:
            continue
        bits = [i for i in range(n) if mask & (1 << i)]
        indep = True
        for a in range(len(bits)):
            for b in range(a+1, len(bits)):
                if adj[bits[a]][bits[b]]:
                    indep = False
                    break
            if not indep:
                break
        if indep:
            alpha[k] += 1
    return dict(alpha)

# ========================================================================
print("=" * 72)
print("PART I: Ω GRAPH STRUCTURE FOR p=7")
print("=" * 72)

p = 7
m = (p - 1) // 2
S_int = list(range(1, m + 1))
QR = get_QR(p)

for name, S in [("Interval", S_int), ("Paley", QR)]:
    A = make_tournament(p, S)
    t0 = time.time()
    cycles = find_all_directed_odd_cycles(A, p)
    t1 = time.time()

    print(f"\n  {name} (S={S}):")
    print(f"    Odd cycles: {len(cycles)} [{t1-t0:.2f}s]")

    # Cycle sizes
    size_dist = defaultdict(int)
    for c in cycles:
        size_dist[len(c)] += 1
    for sz in sorted(size_dist):
        print(f"      Length {sz}: {size_dist[sz]} cycles")

    # Build Ω
    adj = build_omega_graph(cycles)
    nc = len(cycles)

    # Degree sequence
    degrees = [sum(adj[i]) for i in range(nc)]
    print(f"    |V(Ω)| = {nc}")
    print(f"    |E(Ω)| = {sum(degrees)//2}")
    print(f"    Degree: min={min(degrees)}, max={max(degrees)}, avg={np.mean(degrees):.2f}")
    print(f"    Edge density = {sum(degrees) / (nc*(nc-1)) if nc > 1 else 0:.4f}")

    # Degree distribution
    deg_dist = defaultdict(int)
    for d in degrees:
        deg_dist[d] += 1
    print(f"    Degree distribution: {dict(sorted(deg_dist.items()))}")

    # Independence polynomial
    alpha = count_independent_sets(adj, nc)
    if alpha:
        H = sum(alpha.get(k, 0) * 2**k for k in range(nc+1))
        print(f"    Independence polynomial α_k and 2^k·α_k:")
        for k in sorted(alpha):
            if alpha[k] > 0:
                print(f"      α_{k} = {alpha[k]:>6}, 2^{k}·α_{k} = {alpha[k]*2**k:>8}")
        print(f"    H = Z(Ω, 2) = {H}")
        print(f"    α(Ω) = max independent set size = {max(k for k in alpha if alpha[k] > 0)}")

    # Scott-Sokal bound: Z ≤ Π_v (1+λ)^{1/(d_v+1)}
    log_ss_bound = sum(np.log(3) / (d + 1) for d in degrees)  # λ=2
    print(f"    Scott-Sokal upper bound Z(Ω,2) ≤ {np.exp(log_ss_bound):.2f}")

    # Hoffman bound: α(Ω) ≥ n * (-λ_min) / (λ_max - λ_min)
    if nc > 0:
        adj_mat = np.array([[1 if adj[i][j] else 0 for j in range(nc)] for i in range(nc)], dtype=float)
        if nc > 1:
            eigs = np.linalg.eigvalsh(adj_mat)
            lambda_max = eigs[-1]
            lambda_min = eigs[0]
            hoffman_bound = nc * (-lambda_min) / (lambda_max - lambda_min) if lambda_max > lambda_min else nc
            print(f"    Hoffman bound: α(Ω) ≥ {hoffman_bound:.2f}")
            print(f"    Spectral gap of Ω: λ_max={lambda_max:.2f}, λ_min={lambda_min:.2f}")

    # Block structure: connected components of complement graph?
    # Or: clustering coefficient
    if nc > 1:
        triangles = 0
        total_triples = 0
        for i in range(nc):
            neighbors = [j for j in range(nc) if adj[i][j]]
            for a in range(len(neighbors)):
                for b in range(a+1, len(neighbors)):
                    total_triples += 1
                    if adj[neighbors[a]][neighbors[b]]:
                        triangles += 1
        clustering = triangles / total_triples if total_triples > 0 else 0
        print(f"    Clustering coefficient: {clustering:.4f}")

# ========================================================================
print("\n" + "=" * 72)
print("PART II: Ω GRAPH STRUCTURE FOR p=11")
print("=" * 72)

p = 11
m = (p - 1) // 2
S_int = list(range(1, m + 1))
QR = get_QR(p)

for name, S in [("Interval", S_int), ("Paley", QR)]:
    A = make_tournament(p, S)
    t0 = time.time()
    # Only find cycles up to length 7 (length 9 and 11 too slow for full enum)
    cycles = find_all_directed_odd_cycles(A, p, max_len=7)
    t1 = time.time()

    print(f"\n  {name} (S={S[:3]}...):")
    print(f"    Odd cycles (len ≤ 7): {len(cycles)} [{t1-t0:.1f}s]")

    size_dist = defaultdict(int)
    for c in cycles:
        size_dist[len(c)] += 1
    for sz in sorted(size_dist):
        print(f"      Length {sz}: {size_dist[sz]} cycles")

    adj = build_omega_graph(cycles)
    nc = len(cycles)

    degrees = [sum(adj[i]) for i in range(nc)]
    print(f"    |V(Ω≤7)| = {nc}")
    print(f"    |E(Ω≤7)| = {sum(degrees)//2}")
    print(f"    Degree: min={min(degrees)}, max={max(degrees)}, avg={np.mean(degrees):.2f}")
    print(f"    Edge density = {sum(degrees) / (nc*(nc-1)) if nc > 1 else 0:.4f}")

    # Check if brute-force independence polynomial is feasible
    if nc <= 22:
        alpha = count_independent_sets(adj, nc)
        if alpha:
            H_partial = sum(alpha.get(k, 0) * 2**k for k in range(nc+1))
            print(f"    Partial Z(Ω≤7, 2) = {H_partial}")
            print(f"    α_k: ", end="")
            for k in sorted(alpha):
                if alpha[k] > 0:
                    print(f"α_{k}={alpha[k]} ", end="")
            print()
    else:
        print(f"    Too many cycles ({nc}) for brute force independence polynomial")

    # Graph properties
    if nc > 1:
        adj_mat = np.array([[1 if adj[i][j] else 0 for j in range(nc)] for i in range(nc)], dtype=float)
        eigs = np.linalg.eigvalsh(adj_mat)
        print(f"    Spectral: λ_max={eigs[-1]:.2f}, λ_min={eigs[0]:.2f}")
        hoffman_bound = nc * (-eigs[0]) / (eigs[-1] - eigs[0]) if eigs[-1] > eigs[0] else nc
        print(f"    Hoffman bound: α(Ω≤7) ≥ {hoffman_bound:.2f}")

# ========================================================================
print("\n" + "=" * 72)
print("PART III: WHY BLOCK STRUCTURE HELPS — THE TURÁN ANALOGY")
print("=" * 72)
print("""
OBSERVATION from p=7 data above:
  - Interval Ω: cycles cluster → Ω has BLOCK structure
  - Paley Ω: cycles spread → Ω is more "random-like"

The TURÁN GRAPH T(n,r) maximizes edges with no (r+1)-clique.
By complement, the Turán complement MINIMIZES the independence number
for given edge count.

But we want the MAXIMUM Z(Ω, 2), not just max α.
For the independence polynomial at fugacity λ=2:
  Z(G, 2) = Σ_k α_k 2^k

This is maximized when α_k is large for large k, which requires
large independent sets.

KEY INSIGHT: For a graph with vertex set partitioned into cliques
(block structure), the independence polynomial factors:
  Z(G, 2) ≈ Π_{blocks B} (1 + 2 · choice from B)

For COMPLETE MULTIPARTITE graph K_{n₁,...,n_r}:
  α_k = Σ_{i₁<...<i_k} 1 (choose one vertex from each of k blocks)
  Z = Π_i (1 + n_i · λ)

For a graph with partial block structure (not exactly multipartite):
  Z is between the multipartite value (max) and the complete graph value (min).

PREDICTION: Ω(Interval) has MORE multipartite-like structure than Ω(Paley),
leading to higher Z = H.

This can be quantified by the MULTIPARTITE DISTANCE:
  d(G) = min_r || adj(G) - adj(K_{n₁,...,n_r}) ||_F
""")

# For p=7, check multipartite structure
print("p=7: Checking multipartite structure of Ω")
p = 7
m = (p - 1) // 2

for name, S in [("Interval", list(range(1, m+1))), ("Paley", get_QR(p))]:
    A = make_tournament(p, S)
    cycles = find_all_directed_odd_cycles(A, p)
    nc = len(cycles)
    adj = build_omega_graph(cycles)

    # Check if complement is a disjoint union of cliques
    # (i.e., if Ω is complete multipartite)
    # Complement: two cycles NON-adjacent in Ω iff disjoint
    # So complement edges = pairs of DISJOINT cycles
    # Ω is complete multipartite iff the NON-edges form a union of cliques
    # (i.e., mutual disjointness is an equivalence relation)

    # Check transitivity of disjointness
    violations = 0
    for i in range(nc):
        for j in range(i+1, nc):
            if not adj[i][j]:  # i,j disjoint
                for k in range(j+1, nc):
                    if not adj[i][k] and adj[j][k]:  # i,k disjoint but j,k NOT
                        violations += 1
                    if adj[i][k] and not adj[j][k]:
                        violations += 1

    print(f"\n  {name}: {nc} cycles, transitivity violations = {violations}")

    # Find connected components of complement
    comp_adj = [[not adj[i][j] and i != j for j in range(nc)] for i in range(nc)]
    visited = [False] * nc
    components = []
    for start in range(nc):
        if visited[start]:
            continue
        comp = []
        stack = [start]
        while stack:
            v = stack.pop()
            if visited[v]:
                continue
            visited[v] = True
            comp.append(v)
            for w in range(nc):
                if comp_adj[v][w] and not visited[w]:
                    stack.append(w)
        components.append(comp)

    print(f"    Complement components (disjoint cycle groups): {[len(c) for c in components]}")

    # Check if each component is a clique in complement (= independent set in Ω)
    for i, comp in enumerate(components):
        is_clique = True
        for a in range(len(comp)):
            for b in range(a+1, len(comp)):
                if adj[comp[a]][comp[b]]:
                    is_clique = False
                    break
            if not is_clique:
                break
        status = "CLIQUE" if is_clique else "NOT clique"
        if len(comp) > 1:
            print(f"      Component {i}: size {len(comp)}, {status}")

    # Show cycle vertex sets per component
    for i, comp in enumerate(components):
        if len(comp) <= 4:
            cycle_sets = [cycles[c] for c in comp]
            print(f"      Component {i}: {[sorted(s) for s in cycle_sets]}")

# ========================================================================
print("\n" + "=" * 72)
print("PART IV: VERTEX USAGE CONCENTRATION")
print("=" * 72)

# Another way to see clustering: how evenly are vertices used by odd cycles?
# For Interval: cycles concentrate on nearby vertices → uneven usage
# For Paley: cycles spread → even usage

print("Vertex participation in odd cycles (p=7):")
p = 7
for name, S in [("Interval", list(range(1, m+1))), ("Paley", get_QR(p))]:
    A = make_tournament(p, S)
    cycles = find_all_directed_odd_cycles(A, p)

    vertex_count = defaultdict(int)
    for c in cycles:
        for v in c:
            vertex_count[v] += 1

    counts = [vertex_count[v] for v in range(p)]
    print(f"\n  {name}: vertex participation = {counts}")
    print(f"    Mean = {np.mean(counts):.2f}, Std = {np.std(counts):.2f}")
    print(f"    Concentration (std/mean) = {np.std(counts)/np.mean(counts):.4f}")

# ========================================================================
print("\n" + "=" * 72)
print("PART V: KAHN-GALVIN-TETALI DIRECTION")
print("=" * 72)
print("""
THEOREM (Kahn, 2001; Galvin-Tetali, 2004): For a d-regular graph G:
  Z(G, λ)^{1/|V|} ≤ Z(K_{d,d}, λ)^{1/(2d)}

The bound is tight for disjoint unions of K_{d,d}.

For OUR problem: Ω is NOT regular. But we can use the ENTROPY bound:

THEOREM (Galvin-Tetali): For any graph G,
  Z(G, λ) ≤ exp(Σ_v log Z(K_{1,d_v}, λ) / (d_v + 1))

where d_v = degree of vertex v.

This means: vertices with LOWER degree contribute MORE to Z.
So if Ω(Interval) has more LOW-degree vertices than Ω(Paley),
it has higher Z = H.

But does it? Let's check the degree distributions.
""")

# Compare degree distributions
print("Degree distribution comparison (p=7):")
p = 7
for name, S in [("Interval", list(range(1, m+1))), ("Paley", get_QR(p))]:
    A = make_tournament(p, S)
    cycles = find_all_directed_odd_cycles(A, p)
    adj = build_omega_graph(cycles)
    nc = len(cycles)
    degrees = sorted(sum(adj[i]) for i in range(nc))

    # Galvin-Tetali bound
    gt_bound = 0
    for d in degrees:
        # Z(K_{1,d}, λ) = 1 + (d+1)λ for star K_{1,d}
        # Actually K_{1,d} has 1+d vertices, Z = (1+λ)^d + λ(1+λ)^{d} ?
        # No. K_{1,d} (star): independent sets = {}, {center}, or any subset of leaves.
        # Z(K_{1,d}, λ) = (1+λ)^d + λ  (choose any subset of d leaves, or choose center)
        z_star = (1 + 2)**d + 2  # λ=2
        gt_bound += np.log(z_star) / (d + 1)

    print(f"  {name}: degrees = {degrees}")
    print(f"    Galvin-Tetali bound: Z ≤ {np.exp(gt_bound):.2f}")
    alpha_dict = count_independent_sets(adj, nc)
    if alpha_dict:
        print(f"    Actual Z = H = {sum(alpha_dict.get(k,0)*2**k for k in range(nc+1))}")
    else:
        print(f"    (too many vertices for brute force)")

print("\nDONE.")
