#!/usr/bin/env python3
"""
pi_zeta_89c.py — Tournament zeta functions and Euler products

opus-2026-03-14-S89c

Thread 1: H(P_p) factorizations — are they always products of distinct primes?
Thread 2: Tournament zeta function — does it have an Euler product?
Thread 3: Hard-core model on G(T) — critical fugacity and π
Thread 4: More general divisibility: does n | H(T) for regular T on n vertices?
"""

import math
import itertools
import random
from collections import Counter
from fractions import Fraction

def count_H(adj, n):
    dp = [0] * ((1 << n) * n)
    for v in range(n):
        dp[(1 << v) * n + v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            val = dp[mask * n + v]
            if val == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[(mask | (1 << u)) * n + u] += val
    full = (1 << n) - 1
    return sum(dp[full * n + v] for v in range(n))

def is_qr(a, p):
    if a % p == 0:
        return False
    return pow(a, (p-1)//2, p) == 1

def paley_tournament(p):
    adj = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and is_qr(j - i, p):
                adj[i][j] = 1
    return adj

def cyclic_tournament(n):
    adj = [[0]*n for _ in range(n)]
    forward = set(range(1, n//2 + 1))
    for i in range(n):
        for j in range(n):
            if i != j and (j - i) % n in forward:
                adj[i][j] = 1
    return adj

def factor(n):
    """Return prime factorization as list of (p, e) pairs."""
    factors = []
    d = 2
    while d * d <= n:
        e = 0
        while n % d == 0:
            n //= d
            e += 1
        if e > 0:
            factors.append((d, e))
        d += 1
    if n > 1:
        factors.append((n, 1))
    return factors

# ═══════════════════════════════════════════════════════════════
print("=" * 70)
print("PART 1: H(P_p) Factorizations")
print("=" * 70)

paley_h = {}
for p in [3, 7, 11, 19]:
    if p % 4 != 3:
        continue
    adj = paley_tournament(p)
    H = count_H(adj, p)
    paley_h[p] = H
    facts = factor(H)
    fact_str = " × ".join(f"{p_}^{e}" if e > 1 else str(p_) for p_, e in facts)
    print(f"  H(P_{p}) = {H} = {fact_str}")
    print(f"    Distinct prime factors: {[p_ for p_, e in facts]}")
    print(f"    Number of prime factors (with mult): {sum(e for _, e in facts)}")
    print(f"    p divides H: {H % p == 0}")

# H(P_3) = 3
# H(P_7) = 3³ × 7
# H(P_11) = 5 × 7 × 11 × 13 × 19
# H(P_19) = ?

# Factor H(P_19)
H19 = paley_h[19]
facts19 = factor(H19)
fact_str = " × ".join(f"{p_}^{e}" if e > 1 else str(p_) for p_, e in facts19)
print(f"\n  H(P_19) = {H19}")
print(f"         = {fact_str}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 2: Divisibility by p — General Circulant Tournaments")
print("=" * 70)

# THM-212: p | H(T) for any circulant tournament on p vertices
# This follows because the cyclic shift σ: i↦i+1 acts freely on HPs.
# A Hamiltonian path π is a permutation of Z/p. σ maps π to σ∘π.
# σ∘π = π iff π(i)+1 = π(i) for all i, which is impossible.
# So orbits have size p, hence p | H.

# Test: does this extend to n | H for ANY regular tournament on n vertices (n odd)?
# Not all regular tournaments are circulant, but the divisibility might still hold.

print("\nn | H(T) for all regular tournaments on n vertices:")
for n in [3, 5, 7]:
    regular_h = []
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)
    target_score = (n-1)//2

    for bits in range(1 << m):
        adj = [[0]*n for _ in range(n)]
        for idx, (i, j) in enumerate(pairs):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        # Check if regular
        scores = [sum(adj[i][j] for j in range(n)) for i in range(n)]
        if all(s == target_score for s in scores):
            H = count_H(adj, n)
            regular_h.append(H)

    print(f"\n  n={n}: {len(regular_h)} regular tournaments")
    h_mod_n = Counter(h % n for h in regular_h)
    print(f"    H mod {n}: {dict(sorted(h_mod_n.items()))}")
    all_div = all(h % n == 0 for h in regular_h)
    print(f"    ALL divisible by {n}: {all_div}")
    if not all_div:
        counterexamples = [h for h in regular_h if h % n != 0]
        print(f"    Counterexamples: {counterexamples[:5]}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 3: H mod p for ALL tournaments on p vertices")
print("=" * 70)

# What is the distribution of H mod p over ALL tournaments?
for p in [3, 5, 7]:
    h_mod_dist = Counter()
    for bits in range(1 << (p*(p-1)//2)):
        adj = [[0]*p for _ in range(p)]
        pairs = [(i, j) for i in range(p) for j in range(i+1, p)]
        for idx, (i, j) in enumerate(pairs):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        H = count_H(adj, p)
        h_mod_dist[H % p] += 1

    N = 1 << (p*(p-1)//2)
    print(f"\n  n=p={p}: H mod {p} distribution:")
    for r in range(p):
        count = h_mod_dist.get(r, 0)
        print(f"    H ≡ {r} mod {p}: {count} ({count*100/N:.1f}%)")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 4: The Tournament Zeta Function and Euler Product")
print("=" * 70)

# Define for a tournament T:
# ζ_T(s) = Σ_{C: directed odd cycle} |C|^{-s}
# = Σ_k t_k · k^{-s}
# This is analogous to ζ(s) = Σ n^{-s}

# For an "Euler product" to exist, we'd need ζ_T(s) = Π over "prime cycles"
# The "prime" cycles are those not decomposable into smaller ones.
# In tournament context, this doesn't have a clean product structure
# because cycles can overlap (share vertices).

# However, for INDEPENDENT cycles (vertex-disjoint), there IS a product:
# IP(G, x) = Π over connected components of G of IP(component, x)
# If G has components G_1, ..., G_k, then IP(G,x) = Π IP(G_i, x)

# For the odd-cycle disjointness graph G(T), the connected components
# correspond to "clusters" of mutually non-disjoint odd cycles.

print("\nOdd-cycle disjointness graph components:")
for n in range(3, 7):
    pairs_list = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs_list)

    component_sizes = []
    for bits in range(0, min(1 << m, 100)):
        adj = [[0]*n for _ in range(n)]
        for idx, (i, j) in enumerate(pairs_list):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1

        # Find all odd cycles
        cycles = []
        for length in range(3, n+1, 2):
            for combo in itertools.combinations(range(n), length):
                for perm in itertools.permutations(combo):
                    if all(adj[perm[i]][perm[(i+1)%length]] for i in range(length)):
                        min_idx = min(range(length), key=lambda i: perm[i])
                        normalized = perm[min_idx:] + perm[:min_idx]
                        if normalized not in cycles:
                            cycles.append(normalized)

        if len(cycles) == 0:
            component_sizes.append([])
            continue

        # Build disjointness graph and find connected components
        # Adjacency: two cycles are adjacent if they share a vertex
        cycle_verts = [set(c) for c in cycles]
        nc = len(cycles)
        visited = [False]*nc
        components = []
        for start in range(nc):
            if visited[start]:
                continue
            comp = []
            stack = [start]
            while stack:
                node = stack.pop()
                if visited[node]:
                    continue
                visited[node] = True
                comp.append(node)
                for other in range(nc):
                    if not visited[other] and cycle_verts[node] & cycle_verts[other]:
                        stack.append(other)
            components.append(len(comp))

        component_sizes.append(sorted(components, reverse=True))

    if n <= 4:
        print(f"\n  n={n}: first few G(T) component structures:")
        for i, cs in enumerate(component_sizes[:10]):
            print(f"    T_{i}: components = {cs}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 5: Hard-Core Model — Fugacity and Phase Transition")
print("=" * 70)

# H(T) = Z(G(T), λ=2) = partition function at fugacity 2
# The hard-core model on a graph G has a phase transition at
# λ_c(G) = 1/(Δ-1)^{Δ-1}/Δ^Δ ≈ 1/(eΔ) for large Δ
# where Δ = max degree of G

# For the odd-cycle disjointness graph G(T), what is the max degree?
# Cycle C is adjacent to all cycles sharing a vertex with C.

print("\nMax degree of G(T) for random tournaments:")
for n in [5, 7, 9, 11]:
    max_degrees = []
    samples = 50 if n <= 7 else 20
    for _ in range(samples):
        adj = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    adj[i][j] = 1
                else:
                    adj[j][i] = 1

        # Find all directed odd cycles (just 3-cycles for speed)
        cycles = []
        for a, b, c in itertools.combinations(range(n), 3):
            if adj[a][b] and adj[b][c] and adj[c][a]:
                cycles.append(frozenset([a,b,c]))
            if adj[a][c] and adj[c][b] and adj[b][a]:
                cycles.append(frozenset([a,b,c]))

        if not cycles:
            max_degrees.append(0)
            continue

        # Degree of each cycle (sharing vertices)
        max_deg = 0
        for i, c1 in enumerate(cycles):
            deg = sum(1 for j, c2 in enumerate(cycles) if i != j and c1 & c2)
            max_deg = max(max_deg, deg)
        max_degrees.append(max_deg)

    avg_max = sum(max_degrees) / len(max_degrees) if max_degrees else 0
    print(f"  n={n}: avg max degree Δ(G) = {avg_max:.1f}")
    if avg_max > 1:
        lc = (avg_max - 1)**(avg_max-1) / avg_max**avg_max if avg_max > 1 else float('inf')
        print(f"    Estimated λ_c ≈ {lc:.4f}")
        print(f"    We evaluate at λ=2: {'supercritical' if 2 > lc else 'subcritical'}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 6: H(P_19) = ? — Full Factorization")
print("=" * 70)

# We already computed H(P_19) = 1172695746915
H = 1172695746915
facts = factor(H)
print(f"  H(P_19) = {H}")
print(f"         = {' × '.join(f'{p}^{e}' if e > 1 else str(p) for p, e in facts)}")
print(f"  Prime factors: {[p for p, e in facts]}")
print(f"  Number of distinct prime factors: {len(facts)}")
print(f"  Squarefree? {all(e == 1 for _, e in facts)}")

# H(P_3) = 3: 1 prime factor, squarefree
# H(P_7) = 3³ × 7: 2 distinct primes, NOT squarefree
# H(P_11) = 5 × 7 × 11 × 13 × 19: 5 primes, squarefree
# H(P_19) = ?

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 7: The Residue H mod p² for Paley")
print("=" * 70)

for p in [3, 7, 11, 19]:
    if p % 4 != 3:
        continue
    H = paley_h[p]
    print(f"  H(P_{p}) mod p = {H % p}, mod p² = {H % p**2}, mod 2p = {H % (2*p)}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 8: Maximum H Achieved by Doubly Regular Tournaments?")
print("=" * 70)

# A tournament is doubly regular if every pair of vertices has the same
# number of common outneighbors. Paley tournaments are doubly regular.
# For p ≡ 3 mod 4: common outneighbors = (p-3)/4

# Does doubly regular → max H?
# We saw Paley is max at p=3,7,11 but NOT at p=19!
# At p=19: H(C_19) > H(P_19)
# C_19 is also doubly regular? Let me check.

for p in [7, 11, 19]:
    print(f"\n  n={p}:")

    # Check doubly regular for cyclic tournament
    adj_c = cyclic_tournament(p)
    common_c = set()
    for i in range(p):
        for j in range(i+1, p):
            c = sum(1 for k in range(p) if k != i and k != j and adj_c[i][k] and adj_c[j][k])
            common_c.add(c)
    print(f"    Cyclic C_{p}: doubly regular? {len(common_c) == 1} (common values: {common_c})")

    # Check for Paley
    adj_p = paley_tournament(p)
    common_p = set()
    for i in range(p):
        for j in range(i+1, p):
            c = sum(1 for k in range(p) if k != i and k != j and adj_p[i][k] and adj_p[j][k])
            common_p.add(c)
    print(f"    Paley P_{p}: doubly regular? {len(common_p) == 1} (common values: {common_p})")

    H_c = count_H(adj_c, p)
    H_p = count_H(adj_p, p)
    print(f"    H(C_{p}) = {H_c}, H(P_{p}) = {H_p}")
    print(f"    Winner: {'Cyclic' if H_c > H_p else 'Paley' if H_p > H_c else 'Tie'}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 9: π in H(P_p) / p! — Stirling's Approximation")
print("=" * 70)

# H(P_p)/p! is the fraction of all permutations that are Hamiltonian paths
# By Stirling: p! ≈ √(2πp)(p/e)^p
# So H(P_p)/p! ≈ H(P_p) / (√(2πp)(p/e)^p)

print("\nH(P_p)/p! and the √(2πp) connection:")
for p in [3, 7, 11, 19]:
    if p % 4 != 3:
        continue
    H = paley_h[p]
    pf = math.factorial(p)
    ratio = H / pf
    stirling_correction = ratio * math.sqrt(2 * math.pi * p)
    log_ratio = math.log(ratio)
    print(f"  p={p:2d}: H/p! = {ratio:.8e}")
    print(f"         H/p! × √(2πp) = {stirling_correction:.8e}")
    print(f"         log(H/p!) = {log_ratio:.4f}")
    # Compare to -(p-1)log2 + log(const)
    pred = -(p-1)*math.log(2) + math.log(2.4)  # max/mean ≈ 2.4
    print(f"         -(p-1)ln2 + ln(2.4) = {pred:.4f}")

print("\n" + "=" * 70)
print("END — opus-2026-03-14-S89c zeta")
print("=" * 70)
