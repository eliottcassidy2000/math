#!/usr/bin/env python3
"""
Theoretical analysis of H=21 impossibility.

H(T) = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...
For H=21: 2*a1 + 4*a2 + 8*a3 + ... = 20

Decompositions (a1, a2, a3):
  (10, 0, 0), (8, 1, 0), (6, 2, 0), (4, 3, 0), (2, 4, 0), (0, 5, 0)
  (6, 0, 1), (4, 1, 1), (2, 2, 1), (0, 3, 1)
  (2, 0, 2), (0, 1, 2)

For each, we need to check:
1. Is there a tournament T on SOME n with alpha_k(Omega(T)) matching?
2. The key constraint: alpha_2 = number of disjoint odd cycle PAIRS
   If alpha_2 >= 1, need two vertex-disjoint odd cycles (min 6 vertices for two 3-cycles)
   If alpha_2 >= 2, need at least one 3-cycle pair + another cycle
   alpha_3: three pairwise disjoint odd cycles (min 9 vertices)

For n-vertex constraints:
- Max alpha_1(n) = total number of directed odd cycles ~ grows fast
- alpha_2 at n=6: max 1 (exactly one pair of disjoint 3-cycles using all 6 vertices)
  Wait, can we have disjoint 3-cycles using only 6 vertices? Yes: {0,1,2} and {3,4,5}
  But also 3-cycle + 5-cycle? Needs 8 vertices.
  Two disjoint 5-cycles? Needs 10 vertices.

The challenge for H=21:
- At small n (n=3,4,5): max H is 3, 5, 15. Too small.
- At n=6: H ranges from 1 to 45. H=21 would need to be achievable.
  OCF: H = 1 + 2*(t3+t5) + 4*bc where bc = #{disjoint 3-cycle pairs}
  Max bc = 1 at n=6. Max t3+t5 = ?

  For H=21 at n=6: 2*(t3+t5) + 4*bc = 20
  If bc=0: t3+t5 = 10
  If bc=1: t3+t5 = 8

  What are the achievable (t3, t5, bc) at n=6?

- At n=7: H ranges from 1 to 189. H=21 would need to be achievable.
  OCF: H = 1 + 2*(t3+t5+t7) + 4*bc
  where bc = #{disjoint 3-cycle pairs}
  (No disjoint 3+5 or 5+5 at n=7 since 3+5=8>7, 5+5=10>7)
  (No disjoint 3+3+3 at n=7 since 9>7)
  So alpha_2 = bc, alpha_k = 0 for k >= 3.

  For H=21 at n=7: 2*(t3+t5+t7) + 4*bc = 20
  If bc=0: t3+t5+t7 = 10
  If bc=1: t3+t5+t7 = 8
  If bc=2: t3+t5+t7 = 6
  If bc=3: t3+t5+t7 = 4
  ...

Let me compute the achievable (alpha_1, alpha_2) pairs at n=6 and n=7
by exhaustive enumeration.

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

def count_odd_cycles_and_pairs(n, adj_bits):
    """Return (alpha_1, alpha_2, cycle_details) for tournament."""
    edge_set = set()
    for i in range(n):
        for j in range(n):
            if i != j and (adj_bits[i] >> j & 1):
                edge_set.add((i, j))

    # Find all directed odd cycles (by vertex set)
    cycle_vsets = []
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            for p in permutations(verts):
                if all((p[i], p[(i+1) % length]) in edge_set for i in range(length)):
                    cycle_vsets.append(frozenset(verts))
                    break

    alpha1 = len(cycle_vsets)
    alpha2 = 0
    for i in range(len(cycle_vsets)):
        for j in range(i+1, len(cycle_vsets)):
            if len(cycle_vsets[i] & cycle_vsets[j]) == 0:
                alpha2 += 1

    return alpha1, alpha2

def analyze_n(n):
    """Exhaustive analysis at vertex count n."""
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)

    H_to_alpha = defaultdict(set)  # H -> set of (a1, a2) pairs
    alpha_counter = Counter()  # (a1, a2) -> count of tournaments

    for bits in range(1 << m):
        adj = [0] * n
        for k, (i, j) in enumerate(edges):
            if bits & (1 << k):
                adj[j] |= (1 << i)
            else:
                adj[i] |= (1 << j)

        H = held_karp(n, adj)
        a1, a2 = count_odd_cycles_and_pairs(n, adj)

        H_to_alpha[H].add((a1, a2))
        alpha_counter[(a1, a2)] += 1

    return H_to_alpha, alpha_counter

print("=== n=5 ===")
H_to_alpha5, ac5 = analyze_n(5)
print(f"Achievable H values: {sorted(H_to_alpha5.keys())}")
print(f"H=21 achievable: {21 in H_to_alpha5}")
print(f"\nAll (alpha_1, alpha_2) -> H:")
for (a1, a2) in sorted(ac5.keys()):
    H = 1 + 2*a1 + 4*a2
    print(f"  (a1={a1}, a2={a2}): OCF={H}, count={ac5[(a1, a2)]}")

print("\n=== n=6 ===")
H_to_alpha6, ac6 = analyze_n(6)
print(f"Achievable H values: {sorted(H_to_alpha6.keys())}")
print(f"H=21 achievable: {21 in H_to_alpha6}")
print(f"\nH values near 21:")
for h in [17,19,21,23,25]:
    if h in H_to_alpha6:
        print(f"  H={h}: alpha pairs = {sorted(H_to_alpha6[h])}")
    else:
        print(f"  H={h}: NOT achievable")

print(f"\nAll achievable (alpha_1, alpha_2) at n=6:")
for (a1, a2) in sorted(ac6.keys()):
    H = 1 + 2*a1 + 4*a2
    print(f"  (a1={a1}, a2={a2}): OCF={H}, count={ac6[(a1, a2)]}")

# Key question: which (a1, a2) that give H=21 are achievable?
print(f"\n=== H=21 decomposition analysis ===")
print("Required: 2*a1 + 4*a2 = 20")
for a2 in range(6):
    a1 = (20 - 4*a2) // 2
    if a1 < 0:
        break
    if (20 - 4*a2) % 2 != 0:
        continue
    achievable_5 = (a1, a2) in ac5
    achievable_6 = (a1, a2) in ac6
    print(f"  (a1={a1}, a2={a2}): n=5 {'YES' if achievable_5 else 'NO'}, "
          f"n=6 {'YES' if achievable_6 else 'NO'}")
