#!/usr/bin/env python3
"""
per_W_analysis.py — The permanent of W(x) and its relation to F(T,x).

W(x)[u][v] = x if adj[u][v], 1 if adj[v][u], 0 if u=v.

per(W(x)) = sum_{sigma in S_n, sigma derangement} prod_{i} W[i][sigma(i)]
           = sum_{sigma, no fixed point} x^{# i with adj[i][sigma(i)]=1}

This counts DERANGEMENTS weighted by forward edges.
A derangement sigma decomposes into disjoint cycles.
Each cycle of length L contributes x^{# forward edges in cycle}.

QUESTION: How does per(W(x)) relate to F(T,x)?

per(W(x)) counts CYCLE COVERS (no fixed points).
F(T,x) counts Hamiltonian PATHS.

These are connected by Berge's theorem and the OCF!
I(Omega(T), 2) = H(T) via OCF.
The independence polynomial of Omega counts vertex-disjoint odd cycle collections.
Per(W) at x=1 counts derangements = CYCLE COVERS (all cycles, not just odd).

Author: opus-2026-03-07-S45
"""
from itertools import permutations, combinations
import math

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

def compute_F(adj, n):
    F = [0]*n
    for P in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
        F[fwd] += 1
    return F

def compute_per_W(adj, n):
    """per(W(x)) as polynomial in x. Degree up to n (all arcs forward in cycle cover)."""
    per = [0]*(n+1)
    for sigma in permutations(range(n)):
        # Check derangement (no fixed points)
        if any(sigma[i] == i for i in range(n)):
            continue
        fwd = sum(1 for i in range(n) if adj[i][sigma[i]])
        per[fwd] += 1
    return per

def compute_det_W(adj, n):
    """det(W(x)) as polynomial in x."""
    det = [0]*(n+1)
    for sigma in permutations(range(n)):
        if any(sigma[i] == i for i in range(n)):
            continue
        fwd = sum(1 for i in range(n) if adj[i][sigma[i]])
        # Sign of permutation
        inv = sum(1 for i in range(n) for j in range(i+1, n) if sigma[i] > sigma[j])
        sign = (-1)**inv
        det[fwd] += sign
    return det

def cycle_decomposition(sigma):
    """Decompose permutation into cycles."""
    n = len(sigma)
    visited = [False]*n
    cycles = []
    for i in range(n):
        if visited[i]:
            continue
        cycle = []
        j = i
        while not visited[j]:
            visited[j] = True
            cycle.append(j)
            j = sigma[j]
        cycles.append(tuple(cycle))
    return cycles

def compute_odd_per_W(adj, n):
    """per(W(x)) restricted to derangements with ALL ODD cycles."""
    per = [0]*(n+1)
    for sigma in permutations(range(n)):
        if any(sigma[i] == i for i in range(n)):
            continue
        cycles = cycle_decomposition(sigma)
        if any(len(c) % 2 == 0 for c in cycles):
            continue
        fwd = sum(1 for i in range(n) if adj[i][sigma[i]])
        per[fwd] += 1
    return per

# ============================================================
# per(W(x)) vs F(T,x) comparison
# ============================================================
print("=" * 60)
print("per(W(x)) vs F(T,x)")
print("=" * 60)

for n in [4, 5]:
    m = n*(n-1)//2
    print(f"\nn={n}:")
    seen = set()
    for bits in range(1 << m):
        adj = tournament_from_bits(n, bits)
        F = compute_F(adj, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)

        per = compute_per_W(adj, n)
        det = compute_det_W(adj, n)
        odd_per = compute_odd_per_W(adj, n)

        # Trim trailing zeros
        while per and per[-1] == 0:
            per.pop()
        while det and det[-1] == 0:
            det.pop()
        while odd_per and odd_per[-1] == 0:
            odd_per.pop()

        # Is per palindromic?
        d = len(per) - 1
        pal = all(per[k] == per[d-k] for k in range((d+1)//2 + 1)) if d >= 0 else True

        # Is odd_per palindromic?
        d2 = len(odd_per) - 1
        odd_pal = all(odd_per[k] == odd_per[d2-k] for k in range((d2+1)//2 + 1)) if d2 >= 0 else True

        print(f"  F={F}")
        print(f"  per(W)={per} (palindromic={pal})")
        print(f"  odd_per(W)={odd_per} (palindromic={odd_pal})")
        print(f"  det(W)={det}")

        # per(W) at x=1 vs subfactorial
        per_at_1 = sum(per)
        subfact = sum((-1)**k * math.factorial(n) // math.factorial(k) for k in range(n+1))
        # Actually subfactorial D_n = n! * sum_{k=0}^{n} (-1)^k / k!
        D_n = round(math.factorial(n) * sum((-1)**k / math.factorial(k) for k in range(n+1)))
        print(f"  per(W)(1) = {per_at_1}, subfactorial D_{n} = {D_n}")

        # odd_per at x=2
        odd_per_at_2 = sum(odd_per[k] * 2**k for k in range(len(odd_per)))
        print(f"  odd_per(W)(2) = {odd_per_at_2}")

        # H = F at x=... well, F(1) = n!
        # I(Omega, 2) should equal H (OCF)
        # odd_per counts derangements with all odd cycles
        # 2^{#cycles} * (count per independent set) should give... hmm

        print()

# ============================================================
# RELATIONSHIP: odd_per(W, 2) and I(Omega, 2)
# ============================================================
print("=" * 60)
print("ODD CYCLE DERANGEMENTS vs OCF")
print("=" * 60)

# The OCF: H(T) = I(Omega(T), 2) = sum_{S independent in Omega} 2^|S|
# where S ranges over collections of vertex-disjoint odd cycles.

# odd_per(W, x) counts derangements with ALL odd cycles, weighted by x^{fwd}.
# At x=1: odd_per(W, 1) = # derangements with all odd cycles.
# At x=2: ???

# The connection: a derangement with all odd cycles IS a collection of
# vertex-disjoint odd cycles covering ALL vertices. This is a MAXIMUM
# independent set in Omega (covering all vertices).
# I(Omega, 2) counts ALL independent sets (including partial covers), not just maximal.

# So odd_per(W, 1) counts only the all-vertex-covering collections.
# I(Omega, 2) - 1 counts all nonempty independent sets.

# Actually: I(Omega, 2) = sum_{k=0}^{max} a_k * 2^k where a_k = # indep sets of size k.
# This includes the empty set (k=0, contributes 1).
# odd_per(W, 1) counts the derangements, which are the maximum-coverage collections.
# The number of such derangements IS a_max * (number of orientation choices).

# Actually no. An independent set in Omega is a SET of cycles (vertex sets).
# For each such set, there may be multiple derangements (multiple cycle orderings).
# Each 3-cycle vertex set {a,b,c} has TWO directed 3-cycles: a->b->c->a and a->c->b->a.
# But only ONE is present in the tournament.
# A 5-cycle vertex set may have several directed 5-cycles in the tournament.

# So the relationship is: odd_per(W, 1) counts ORIENTED odd-cycle covers,
# while I(Omega, 2) counts UNORIENTED cycle collections with weight 2^|S|.

# For a single 3-cycle in a tournament, there is exactly one directed 3-cycle.
# So odd_per counts the same objects as the independent-set sum, but with different weighting.

n = 5
m = n*(n-1)//2
print(f"\nn={n}:")

seen = set()
for bits in range(1 << m):
    adj = tournament_from_bits(n, bits)
    F = compute_F(adj, n)
    key = tuple(F)
    if key in seen:
        continue
    seen.add(key)

    odd_per = compute_odd_per_W(adj, n)
    odd_per_1 = sum(odd_per)

    # Count 3-cycles
    c3 = 0
    for triple in combinations(range(n), 3):
        a, b, c = triple
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (adj[a][c] and adj[c][b] and adj[b][a]):
            c3 += 1

    # Count 5-cycles (directed)
    c5_directed = 0
    for perm in permutations(range(n)):
        if all(adj[perm[i]][perm[(i+1)%5]] for i in range(5)):
            c5_directed += 1
    c5_directed //= 5  # Each 5-cycle counted 5 times (rotation)

    print(f"  F={F}, c3={c3}, c5_dir={c5_directed}")
    print(f"  odd_per(1) = {odd_per_1}")
    print(f"  odd_per = {odd_per[:6]}")

    # At n=5, an all-odd derangement must be either:
    # - one 5-cycle, or
    # - one 3-cycle + one isolated... wait, all vertices must be covered.
    # n=5: 5-cycle or 3-cycle + (impossible, 2 vertices left can't form odd cycle).
    # So odd_per(W, 1) = c5_directed at n=5.
    print(f"  c5_dir matches odd_per(1)? {c5_directed == odd_per_1}")
    print()
