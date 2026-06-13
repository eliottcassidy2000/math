#!/usr/bin/env python3
"""
H=21 impossibility analysis with CORRECT cycle counting.

BUG FIX: must count ALL directed cycles, not just one per vertex set.
Multiple directed Hamiltonian cycles can exist on the same vertex subset.

opus-2026-03-07-S38
"""
from itertools import combinations, permutations
from collections import Counter, defaultdict

def held_karp(n, adj):
    dp = [[0]*n for _ in range(1<<n)]
    for v in range(n):
        dp[1<<v][v] = 1
    for S in range(1, 1<<n):
        for v in range(n):
            if not (S & (1<<v)): continue
            if dp[S][v] == 0: continue
            for u in range(n):
                if S & (1<<u): continue
                if adj[v] & (1<<u):
                    dp[S|(1<<u)][u] += dp[S][v]
    full = (1<<n)-1
    return sum(dp[full][v] for v in range(n))


def count_directed_cycles(n, adj_bits):
    """Count ALL directed odd cycles and vertex-disjoint pairs.

    A directed cycle = cyclic sequence (v0,v1,...,v_{k-1}) with v0->v1->...->v_{k-1}->v0,
    considered up to cyclic rotation.
    """
    edge_set = set()
    for i in range(n):
        for j in range(n):
            if i != j and (adj_bits[i] >> j & 1):
                edge_set.add((i, j))

    all_cycles = []  # list of (frozenset(vertices), canonical_tuple)

    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            # Find ALL directed Hamiltonian cycles on this vertex set
            seen = set()
            for p in permutations(verts):
                if all((p[i], p[(i+1) % length]) in edge_set for i in range(length)):
                    # Canonicalize: start from min vertex
                    min_idx = list(p).index(min(p))
                    canonical = tuple(list(p)[min_idx:] + list(p)[:min_idx])
                    if canonical not in seen:
                        seen.add(canonical)
                        all_cycles.append((frozenset(verts), canonical))

    alpha1 = len(all_cycles)

    # alpha_2: count pairs of vertex-disjoint cycles
    alpha2 = 0
    for i in range(len(all_cycles)):
        for j in range(i+1, len(all_cycles)):
            if len(all_cycles[i][0] & all_cycles[j][0]) == 0:
                alpha2 += 1

    return alpha1, alpha2


def analyze_n(n):
    """Exhaustive analysis at vertex count n."""
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)

    H_to_alpha = defaultdict(set)
    alpha_to_H = defaultdict(set)
    H_counter = Counter()
    alpha_counter = Counter()
    mismatches = 0

    for bits in range(1 << m):
        adj = [0] * n
        for k, (i, j) in enumerate(edges):
            if bits & (1 << k):
                adj[j] |= (1 << i)
            else:
                adj[i] |= (1 << j)

        H = held_karp(n, adj)
        a1, a2 = count_directed_cycles(n, adj)
        ocf = 1 + 2*a1 + 4*a2

        H_counter[H] += 1
        H_to_alpha[H].add((a1, a2))
        alpha_to_H[(a1, a2)].add(H)
        alpha_counter[(a1, a2)] += 1

        if H != ocf:
            mismatches += 1
            if mismatches <= 3:
                print(f"  MISMATCH: bits={bits}, H={H}, OCF={ocf}, a1={a1}, a2={a2}")

    return H_to_alpha, alpha_counter, H_counter, mismatches


# Quick verification at n=5
print("=== n=5 verification ===")
H5, ac5, hc5, mm5 = analyze_n(5)
print(f"Mismatches: {mm5}")
print(f"Achievable H: {sorted(hc5.keys())}")
print(f"H=21 achievable: {21 in hc5}")
for (a1, a2) in sorted(ac5.keys()):
    ocf = 1 + 2*a1 + 4*a2
    print(f"  (a1={a1}, a2={a2}): OCF={ocf}, count={ac5[(a1, a2)]}")

# n=6 analysis
print("\n=== n=6 analysis ===")
H6, ac6, hc6, mm6 = analyze_n(6)
print(f"Mismatches: {mm6}")
print(f"Achievable H: {sorted(hc6.keys())}")
print(f"H=21 achievable: {21 in hc6}")

print(f"\nH values near 21:")
for h in range(17, 26):
    if h in hc6:
        print(f"  H={h}: count={hc6[h]}, alpha pairs={sorted(H6[h])}")
    else:
        print(f"  H={h}: NOT achievable")

print(f"\nAll achievable (a1, a2) at n=6:")
for (a1, a2) in sorted(ac6.keys()):
    ocf = 1 + 2*a1 + 4*a2
    print(f"  (a1={a1}, a2={a2}): OCF={ocf}, count={ac6[(a1, a2)]}")

# H=21 decomposition check
print(f"\n=== H=21 decomposition analysis ===")
print("Required: 2*a1 + 4*a2 = 20")
for a2 in range(6):
    remainder = 20 - 4*a2
    if remainder < 0:
        break
    if remainder % 2 != 0:
        continue
    a1 = remainder // 2
    achievable_6 = (a1, a2) in ac6
    print(f"  (a1={a1}, a2={a2}): n=6 {'YES' if achievable_6 else 'NO'}")
