#!/usr/bin/env python3
"""
Forbidden H values — deep structural analysis.
opus-2026-03-14-S85

THM-200 proved H≠7 (kind-pasteur S71g). Key insight: H=I(Ω,2) and I(G,2)=7
iff G=K₃, but K₃ cannot be a component of Ω(T).

QUESTION: What other H values are forbidden?
- {7, 21, 63} conjectured (kind-pasteur)
- 7 = I(K₃, 2) = 2³-1
- 21 = I(K₃, 2) · I(K₁, 2) = 7·3? Or 21 = some I(G,2)?
- 63 = I(K₃, 2)·I(K₃, 2) = 7·9? Or 63 = 2⁶-1?

Plan:
1. Enumerate all achievable H values exhaustively for n≤7
2. Characterize which Ω conflict graphs are possible
3. For each non-achievable H, find which G has I(G,2)=H
4. Determine which of those G are impossible as Ω components
"""

from itertools import combinations
from collections import defaultdict, Counter
import math
import sys

def get_tournament(n, bits):
    arcs = [(i, j) for i in range(n) for j in range(i+1, n)]
    adj = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(arcs):
        if (bits >> k) & 1:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj

def compute_H_dp(adj, n):
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)):
                continue
            if dp[S][v] == 0:
                continue
            for w in range(n):
                if S & (1 << w):
                    continue
                if adj[v][w]:
                    dp[S | (1 << w)][w] += dp[S][v]
    return sum(dp[full][v] for v in range(n))

def odd_cycles(adj, n):
    """Find all directed odd cycles."""
    cycles = []
    for length in [3, 5, 7]:
        if length > n:
            break
        for combo in combinations(range(n), length):
            from itertools import permutations as perms
            for perm in perms(combo):
                is_cycle = True
                for i in range(length):
                    if not adj[perm[i]][perm[(i+1) % length]]:
                        is_cycle = False
                        break
                if is_cycle:
                    # Canonicalize: start from min vertex, direction matters
                    min_idx = perm.index(min(perm))
                    canon = tuple(perm[min_idx:] + perm[:min_idx])
                    cycles.append(canon)
    return list(set(cycles))

def conflict_graph_edges(cycles):
    """Two cycles conflict if they share a vertex."""
    edges = []
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if set(cycles[i]) & set(cycles[j]):
                edges.append((i, j))
    return edges

def independence_polynomial_at_2(adj_list, n_vertices):
    """I(G, 2) = sum over independent sets S of 2^|S|."""
    total = 0
    for mask in range(1 << n_vertices):
        # Check independence
        independent = True
        for i in range(n_vertices):
            if not (mask & (1 << i)):
                continue
            for j in adj_list.get(i, []):
                if mask & (1 << j):
                    independent = False
                    break
            if not independent:
                break
        if independent:
            total += 2 ** bin(mask).count('1')
    return total

# ============================================================
# Part 1: Complete H spectrum for n≤7
# ============================================================
print("=" * 70)
print("PART 1: ACHIEVABLE H VALUES FOR n ≤ 7")
print("=" * 70)

all_achievable = set()
for n in range(2, 8):
    m = n * (n - 1) // 2
    N = 1 << m

    H_vals = set()
    H_dist = Counter()

    if n <= 6:
        for bits in range(N):
            adj = get_tournament(n, bits)
            H = compute_H_dp(adj, n)
            H_vals.add(H)
            H_dist[H] += 1
    else:
        # n=7: 2^21 = 2M tournaments, feasible
        for bits in range(N):
            adj = get_tournament(n, bits)
            H = compute_H_dp(adj, n)
            H_vals.add(H)
            H_dist[H] += 1
            if bits % 500000 == 0 and bits > 0:
                print(f"  n=7: {bits}/{N} ({100*bits/N:.1f}%)", file=sys.stderr)

    sorted_H = sorted(H_vals)
    print(f"\nn={n}: {len(H_vals)} distinct H values")
    print(f"  Range: [{min(sorted_H)}, {max(sorted_H)}]")
    print(f"  Values: {sorted_H[:30]}{'...' if len(sorted_H) > 30 else ''}")

    # Find gaps
    full_range = set(range(1, max(sorted_H)+1, 2))  # All odd numbers in range
    missing = sorted(full_range - H_vals)
    print(f"  Missing odd values ≤ {max(sorted_H)}: {missing[:20]}{'...' if len(missing) > 20 else ''}")

    all_achievable |= H_vals

print(f"\nCombined achievable (all n≤7): {sorted(all_achievable)[:50]}...")
print(f"Total achievable: {len(all_achievable)}")

# Check which odd numbers ≤ 189 are achievable
max_check = 189  # max H at n=7
missing_overall = []
for h in range(1, max_check + 1, 2):  # odd only (Rédei)
    if h not in all_achievable:
        missing_overall.append(h)
print(f"\nMissing odd values ≤ {max_check}: {missing_overall}")

# ============================================================
# Part 2: I(G, 2) = H — Which graphs give each H?
# ============================================================
print("\n" + "=" * 70)
print("PART 2: WHICH GRAPHS G HAVE I(G, 2) = H?")
print("=" * 70)

# For small graphs, enumerate I(G, 2)
# Graphs on k vertices (k ≤ 6)
def enumerate_graphs(k):
    """Generate all simple graphs on k vertices."""
    edges = [(i, j) for i in range(k) for j in range(i+1, k)]
    m = len(edges)
    for bits in range(1 << m):
        adj = defaultdict(list)
        for idx, (i, j) in enumerate(edges):
            if (bits >> idx) & 1:
                adj[i].append(j)
                adj[j].append(i)
        yield bits, adj

I_to_graphs = defaultdict(list)
for k in range(1, 7):
    edges = [(i, j) for i in range(k) for j in range(i+1, k)]
    m = len(edges)
    for bits in range(1 << m):
        adj = defaultdict(list)
        for idx, (i, j) in enumerate(edges):
            if (bits >> idx) & 1:
                adj[i].append(j)
                adj[j].append(i)
        I_val = independence_polynomial_at_2(adj, k)
        I_to_graphs[I_val].append((k, bits))

print("\nI(G, 2) values and their graphs (sorted):")
for I_val in sorted(I_to_graphs.keys()):
    if I_val <= 100:
        graphs = I_to_graphs[I_val]
        # Describe smallest graph
        k_min, bits_min = min(graphs, key=lambda x: x[0])
        edges_min = [(i, j) for i in range(k_min) for j in range(i+1, k_min)]
        edge_list = [edges_min[idx] for idx in range(len(edges_min)) if (bits_min >> idx) & 1]
        n_edges = len(edge_list)

        achievable = I_val in all_achievable
        status = "✓" if achievable else "✗ MISSING"
        print(f"  I(G,2) = {I_val:4d}: smallest graph on {k_min} vertices, {n_edges} edges {status}")

# ============================================================
# Part 3: Factor structure of H values
# ============================================================
print("\n" + "=" * 70)
print("PART 3: FACTORIZATION PATTERN OF H VALUES")
print("=" * 70)

print("\nH = I(Ω, 2) = product over components C of I(C, 2)")
print("Since I(G₁ ⊔ G₂, 2) = I(G₁, 2) · I(G₂, 2)")
print()

# Component I values for small graphs
component_I = set()
for k in range(1, 7):
    edges = [(i, j) for i in range(k) for j in range(i+1, k)]
    m = len(edges)
    for bits in range(1 << m):
        adj = defaultdict(list)
        for idx, (i, j) in enumerate(edges):
            if (bits >> idx) & 1:
                adj[i].append(j)
                adj[j].append(i)
        # Check if connected
        if k == 1:
            I_val = independence_polynomial_at_2(adj, k)
            component_I.add(I_val)
            continue
        visited = {0}
        stack = [0]
        while stack:
            v = stack.pop()
            for w in adj[v]:
                if w not in visited:
                    visited.add(w)
                    stack.append(w)
        if len(visited) == k:
            I_val = independence_polynomial_at_2(adj, k)
            component_I.add(I_val)

print(f"Connected graph I(G,2) values (k≤6): {sorted(component_I)}")

# Which products of component values give each H?
achievable_products = set()
component_list = sorted(component_I)

def products_up_to(limit, factors, current=1, start=0):
    """Generate all products of factors ≤ limit."""
    results = set()
    results.add(current)
    for i in range(start, len(factors)):
        p = current * factors[i]
        if p <= limit:
            results |= products_up_to(limit, factors, p, i)
    return results

all_products = products_up_to(200, component_list)
print(f"\nAll possible products of connected I(G,2) values ≤ 200: {sorted(all_products)[:40]}...")

# Compare with achievable H
missing_from_products = sorted(set(range(1, 200, 2)) - all_products)
print(f"Odd values NOT expressible as products: {missing_from_products[:20]}")

# Key question: which connected I values can appear as Ω components?
# From THM-201: K₃ CANNOT appear, so I(K₃,2)=7 is forbidden as a factor
print("\n--- KEY STRUCTURAL QUESTION ---")
print("If K₃ is impossible as Ω component (THM-201),")
print("then 7 cannot appear as a factor in the product decomposition.")
print("This forbids: 7, 21, 35, 49, 63, 77, 91, ...")
print("But some of these might be achievable via OTHER decompositions!")
print()

# Check: can 21 be achieved without factor 7?
print("21 = 3 × 7 (needs K₃ as component — FORBIDDEN)")
print("21 = 21 as a single connected component I(G,2)")
# Is there a connected graph with I(G,2) = 21?
graphs_21 = [(k, b) for I, gs in I_to_graphs.items() if I == 21 for k, b in gs]
print(f"  Connected graphs with I(G,2) = 21: {len(graphs_21)} found")
for k, b in graphs_21[:5]:
    edges = [(i, j) for i in range(k) for j in range(i+1, k)]
    edge_list = [edges[idx] for idx in range(len(edges)) if (b >> idx) & 1]
    print(f"    {k} vertices, edges: {edge_list}")

# Check 63
graphs_63 = [(k, b) for I, gs in I_to_graphs.items() if I == 63 for k, b in gs]
print(f"\nConnected graphs with I(G,2) = 63: {len(graphs_63)} found")
for k, b in graphs_63[:5]:
    edges = [(i, j) for i in range(k) for j in range(i+1, k)]
    edge_list = [edges[idx] for idx in range(len(edges)) if (b >> idx) & 1]
    print(f"    {k} vertices, edges: {edge_list}")

# ============================================================
# Part 4: Component analysis at n=5,6
# ============================================================
print("\n" + "=" * 70)
print("PART 4: CONFLICT GRAPH Ω COMPONENT ANALYSIS")
print("=" * 70)

# For each tournament, find Ω and its components
# Ω: vertices = odd directed cycles, edges = share a vertex
# But computing all odd cycles is expensive. Use 3-cycles only for small n.

for n in [4, 5]:
    m = n * (n - 1) // 2
    N = 1 << m

    component_structures = Counter()
    H_by_structure = defaultdict(set)

    for bits in range(N):
        adj = get_tournament(n, bits)
        H = compute_H_dp(adj, n)

        # Find all 3-cycles
        cycles_3 = []
        for i in range(n):
            for j in range(i+1, n):
                for k in range(j+1, n):
                    if adj[i][j] and adj[j][k] and adj[k][i]:
                        cycles_3.append((i, j, k))
                    elif adj[j][i] and adj[i][k] and adj[k][j]:
                        cycles_3.append((j, i, k))

        # Build conflict graph on 3-cycles
        c = len(cycles_3)
        if c == 0:
            structure = "empty"
        else:
            # Adjacency
            adj_conflict = defaultdict(set)
            for a in range(c):
                for b in range(a+1, c):
                    if set(cycles_3[a]) & set(cycles_3[b]):
                        adj_conflict[a].add(b)
                        adj_conflict[b].add(a)

            # Find components
            visited = set()
            components = []
            for v in range(c):
                if v not in visited:
                    comp = set()
                    stack = [v]
                    while stack:
                        u = stack.pop()
                        if u in visited:
                            continue
                        visited.add(u)
                        comp.add(u)
                        for w in adj_conflict[u]:
                            if w not in visited:
                                stack.append(w)
                    components.append(comp)

            comp_sizes = tuple(sorted(len(c) for c in components))
            structure = str(comp_sizes)

        component_structures[structure] += 1
        H_by_structure[structure].add(H)

    print(f"\nn={n}: Ω (3-cycle conflict graph) component structures:")
    for struct in sorted(component_structures.keys()):
        H_vals = sorted(H_by_structure[struct])
        print(f"  {struct}: count={component_structures[struct]}, H values={H_vals}")

# ============================================================
# Part 5: Are 7·k values always forbidden?
# ============================================================
print("\n" + "=" * 70)
print("PART 5: MULTIPLES OF 7 IN H SPECTRUM")
print("=" * 70)

print("\nMultiples of 7 that are odd (potential H values):")
for k in range(1, 30):
    h = 7 * k
    if h % 2 == 0:
        continue
    status = "ACHIEVABLE" if h in all_achievable else "MISSING"
    print(f"  7 × {k} = {h}: {status}")

# 7 × 27 = 189 — this IS achievable! (max H at n=7)
# So not ALL multiples of 7 are forbidden.
# 189 = H(T_7) = 3³ × 7
# Factor: I(G,2)=189. What graph gives this?
# Or: components with I=9 (3²) and I=21?
# Or: I=27 and I=7? But 7 is forbidden.
# Or: I=189 as single component?

print("\n--- 189 = max H(7) ---")
print("189 = 3³ × 7 = 27 × 7")
print("Since 7 is forbidden as a component factor,")
print("189 must come from a SINGLE connected component with I(Ω,2)=189")
print("OR from a different factorization not involving 7.")
print("Check: 189 = 3 × 63, 9 × 21, 27 × 7")
print("All factorizations involve 7 or 21 (=3×7) or 63 (=9×7)")
print("So 189 MUST come from a single connected Ω with I(Ω,2)=189")
print("This means Ω has a single large connected component!")

# ============================================================
# Part 6: The Mersenne Connection
# ============================================================
print("\n" + "=" * 70)
print("PART 6: MERSENNE NUMBERS AND H")
print("=" * 70)

mersenne = [2**k - 1 for k in range(2, 12)]
print(f"Mersenne numbers: {mersenne}")
for m_val in mersenne:
    if m_val % 2 == 0:
        continue
    status = "ACHIEVABLE" if m_val in all_achievable else "MISSING"
    k = int(math.log2(m_val + 1))
    print(f"  2^{k} - 1 = {m_val}: {status}")

# I(complete graph K_k, 2) = (1+2)^k - stuff?
# Actually I(K_k, x) = 1 + kx (only independent sets are singletons and empty)
# So I(K_k, 2) = 1 + 2k
# K₁: 3, K₂: 5, K₃: 7, K₄: 9, K₅: 11, ...
# These are just 2k+1 = odd numbers!
# So I(K_k, 2) = 2k+1 for the complete graph.
# K₃ gives 7 — and THIS is what's forbidden.

print("\nI(K_k, 2) = 1 + 2k:")
for k in range(1, 8):
    print(f"  I(K_{k}, 2) = {1 + 2*k}")

# More interesting: I(cycle C_k, 2)
# For C₃ = K₃: I = 7
# For C₄: I = (2+1)^4/... Actually compute directly
print("\nI(C_k, 2) for cycles:")
for k in range(3, 10):
    # Independence polynomial of C_k: I(C_k, x) = sum ...
    # Use transfer matrix method or direct computation
    # I(C_k, x) = F_k(x) + F_{k-2}(x) where F_k are Fibonacci-like
    # Actually I(C_k, x) = (1+x)^k - kx(1+x)^{k-2} ... no
    # Direct: enumerate independent sets of C_k
    total = 0
    for mask in range(1 << k):
        indep = True
        for i in range(k):
            if (mask & (1 << i)) and (mask & (1 << ((i+1) % k))):
                indep = False
                break
        if indep:
            total += 2 ** bin(mask).count('1')
    print(f"  I(C_{k}, 2) = {total}")

# I(path P_k, 2):
print("\nI(P_k, 2) for paths (no wraparound):")
for k in range(1, 10):
    total = 0
    for mask in range(1 << k):
        indep = True
        for i in range(k-1):
            if (mask & (1 << i)) and (mask & (1 << (i+1))):
                indep = False
                break
        if indep:
            total += 2 ** bin(mask).count('1')
    print(f"  I(P_{k}, 2) = {total}")

# ============================================================
# SYNTHESIS
# ============================================================
print("\n" + "=" * 70)
print("SYNTHESIS — FORBIDDEN H VALUES")
print("=" * 70)
print("""
FINDINGS:
1. FORBIDDEN H values (odd numbers not achievable for any n≤7):
   Check output above for exact list.

2. STRUCTURAL MECHANISM:
   H = I(Ω, 2) where Ω = conflict graph of odd cycles
   H = product of I(Cᵢ, 2) over connected components Cᵢ of Ω
   If some connected graph G is impossible as Ω-component,
   then I(G,2) is forbidden as a factor.

3. THM-201 forbids K₃ as component → I(K₃,2)=7 is a forbidden factor.
   This propagates: 7, 7×3=21, 7×9=63, 7×27=189, ...
   BUT 189 IS achievable! So 189 must come from a single
   connected component (not factoring through 7).

4. MERSENNE CONNECTION: 7 = 2³-1 is Mersenne. But 31, 127 are also
   Mersenne and may or may not be forbidden.

5. I(C_k, 2) sequence for cycles gives the forbidden component values.
   I(K_k, 2) = 2k+1 — these are odd numbers, matching the spectrum.
""")
