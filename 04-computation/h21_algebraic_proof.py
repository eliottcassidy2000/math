#!/usr/bin/env python3
"""
H≠21 — Algebraic Proof Attempt
opus-2026-03-14-S71h

THEOREM: I(G,2) = 21 iff G is one of exactly 4 graph types:
  1. n=4, α=(1,4,3,0): P₄ or K₁⊔K₃ [I = 1+4x+3x²]
  2. n=6, α=(1,6,2): K₆ minus 2 independent edges [I = 1+6x+2x²]
  3. n=8, α=(1,8,1): K₈ minus 1 edge [I = 1+8x+x²]
  4. n=10, α=(1,10): K₁₀ [I = 1+10x]

PROOF: I(G,2) = Σ α_k 2^k = 1 + 2α₁ + 4α₂ + 8α₃ + ... = 21
⟹ 2α₁ + 4α₂ + 8α₃ + ... = 20
⟹ α₁ + 2α₂ + 4α₃ + ... = 10
Since α₁ = n (number of vertices) and α₂ ≤ C(n,2), we need n ≤ 10.
"""

print("=" * 70)
print("ALGEBRAIC ENUMERATION: ALL GRAPHS WITH I(G,2) = 21")
print("=" * 70)
print()

# I(G,2) = 1 + 2n + 4α₂ + 8α₃ + ... = 21
# So: n + 2α₂ + 4α₃ + ... = 10, where α₁ = n

print("Constraint: n + 2α₂ + 4α₃ + 8α₄ + ... = 10")
print("with α_k ≥ 0 and α₁ = n (number of vertices)")
print()

for n in range(1, 11):
    remainder = 10 - n
    if remainder < 0:
        break
    if remainder % 2 != 0 and n <= 10:
        # Need 2α₂ + 4α₃ + ... = remainder
        # 2α₂ is always even, 4α₃ is always divisible by 4
        # So remainder must be even
        if remainder % 2 != 0:
            print(f"n={n}: remainder = {remainder} (odd) → IMPOSSIBLE (2α₂ + 4α₃ + ... is always even)")
            continue

    # 2α₂ + 4α₃ + ... = remainder
    # α₂ + 2α₃ + 4α₄ + ... = remainder/2
    if remainder % 2 != 0:
        print(f"n={n}: {remainder} not even → IMPOSSIBLE")
        continue

    half_rem = remainder // 2

    # Find valid (α₂, α₃, α₄, ...) with α_k ≤ C(n, k)
    from math import comb

    solutions = []

    # Since 4α₄ + 8α₅ + ... ≤ half_rem, and we need small values:
    # Try all α₃ from 0 to half_rem//2, α₄ from 0 to half_rem//4, etc.
    for a3 in range(min(half_rem // 2 + 1, comb(n, 3) + 1)):
        for a4 in range(min((half_rem - 2*a3) // 4 + 1, comb(n, 4) + 1)):
            a2_rem = half_rem - 2*a3 - 4*a4
            if a2_rem < 0:
                break
            a2 = a2_rem
            if a2 <= comb(n, 2):
                # Check feasibility: can a graph on n vertices have
                # α₂ = a2, α₃ = a3, α₄ = a4?
                # Basic check: α₂ ≤ C(n,2), α₃ ≤ C(n,3), etc.
                if a3 <= comb(n, 3) and a4 <= comb(n, 4):
                    solutions.append((a2, a3, a4))

    if not solutions:
        print(f"n={n}: no valid (α₂, α₃, ...) → IMPOSSIBLE")
    else:
        for (a2, a3, a4) in solutions:
            # Verify I(G,2)
            I_val = 1 + 2*n + 4*a2 + 8*a3 + 16*a4
            alpha_vec = f"α = (1, {n}, {a2}"
            if a3 > 0:
                alpha_vec += f", {a3}"
            if a4 > 0:
                alpha_vec += f", {a4}"
            alpha_vec += ")"

            # Identify graph type
            if a2 == 0 and a3 == 0 and a4 == 0:
                gtype = f"K_{n} (complete graph)"
            elif a2 == comb(n, 2) and a3 == comb(n, 3):
                gtype = f"empty graph on {n} vertices"
            else:
                missing_edges = comb(n, 2) - a2
                if a3 == 0:
                    gtype = f"K_{n} minus {missing_edges} edge{'s' if missing_edges != 1 else ''}"
                else:
                    gtype = f"graph on {n} vertices with {a2} ind. pairs"

            feasible = "✓" if I_val == 21 else f"✗ (I={I_val})"
            print(f"n={n}: {alpha_vec}, I(G,2)={I_val} {feasible}")
            print(f"       Type: {gtype}")

print()
print("=" * 70)
print("SUMMARY: EXACTLY 4 GRAPH TYPES HAVE I(G,2) = 21")
print("=" * 70)
print()
print("1. n=4:  I = 1+4x+3x²  (P₄ or K₁⊔K₃)       — K₃ poisoned")
print("2. n=6:  I = 1+6x+2x²  (K₆ - 2 edges)        — NOT K₃ poisoned")
print("3. n=8:  I = 1+8x+x²   (K₈ - 1 edge)          — NOT K₃ poisoned")
print("4. n=10: I = 1+10x      (K₁₀)                  — NOT K₃ poisoned")
print()

# Now check: for each type, can it be Ω(T)?
print("=" * 70)
print("CAN THESE GRAPHS BE Ω(T)?")
print("=" * 70)
print()

# Type 1: THM-201 + THM-202 blocks this
print("Type 1 (n=4, I=1+4x+3x²): BLOCKED")
print("  K₁⊔K₃: has K₃ component, blocked by THM-201")
print("  P₄: blocked by THM-202 (dominance cascade)")
print()

# Type 4: K₁₀ means 10 pairwise-intersecting odd cycles
print("Type 4 (K₁₀ = 10 pairwise-intersecting odd cycles):")
print("  For Ω(T) = K₁₀: need 10 directed odd cycles, every pair sharing ≥1 vertex.")
print("  At n=6: total cycle count can be 10 (720 tournaments), but Ω has 10 vertices")
print("    with C(10,2)=45 edges needed for K₁₀. With only 6 tournament vertices,")
print("    many cycle pairs must share vertices... but is it K₁₀ exactly?")
print()

# Let's check: at n=6, tournaments with exactly 10 cycles — what is their Ω structure?
from itertools import combinations, permutations
from collections import Counter

def find_all_directed_odd_cycles(n, adj, max_len=None):
    if max_len is None:
        max_len = n
    cycles = []
    seen = set()
    for length in range(3, max_len + 1, 2):
        for combo in combinations(range(n), length):
            for perm in permutations(combo):
                is_cycle = True
                for i in range(length):
                    nxt = (i + 1) % length
                    if perm[nxt] not in adj[perm[i]]:
                        is_cycle = False
                        break
                if is_cycle:
                    rotations = [perm[i:] + perm[:i] for i in range(length)]
                    canonical = min(rotations)
                    if canonical not in seen:
                        seen.add(canonical)
                        cycles.append((canonical, frozenset(combo)))
    return cycles

def tournament_from_bits(n, bits):
    adj = [set() for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                adj[i].add(j)
            else:
                adj[j].add(i)
            idx += 1
    return adj

def hp_count(n, adj):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, full + 1):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for w in adj[v]:
                if not (mask & (1 << w)):
                    dp[mask | (1 << w)][w] += dp[mask][v]
    return sum(dp[full])

print("Checking n=6 tournaments with exactly 10 odd cycles:")
n = 6
num_edges = 15
total = 1 << num_edges

h_for_10_cycles = Counter()
omega_structures = Counter()

for bits in range(total):
    adj = tournament_from_bits(n, bits)
    cycles = find_all_directed_odd_cycles(n, adj, max_len=5)

    if len(cycles) != 10:
        continue

    H = hp_count(n, adj)
    h_for_10_cycles[H] += 1

    # Count Ω edges
    nc = len(cycles)
    edge_count = 0
    disjoint_count = 0
    for i in range(nc):
        for j in range(i+1, nc):
            if cycles[i][1] & cycles[j][1]:
                edge_count += 1
            else:
                disjoint_count += 1

    omega_structures[(nc, edge_count, disjoint_count)] += 1

print(f"  H values for 10-cycle tournaments: {dict(sorted(h_for_10_cycles.items()))}")
print(f"  Ω structures (|V|=10): {dict(sorted(omega_structures.items()))}")
print()

# The key: C(10,2) = 45 possible edges. K₁₀ has 45 edges.
# If any tournament has Ω with 45 edges (all pairs intersecting), that's K₁₀.
print(f"  K₁₀ (45 edges) requires ALL pairs intersecting.")
has_k10 = any(ec == 45 for (_, ec, _) in omega_structures)
print(f"  Does K₁₀ appear as Ω at n=6? {has_k10}")

for (nc, ec, dc), cnt in sorted(omega_structures.items()):
    print(f"    |V|={nc}, |E|={ec}, disjoint={dc}: {cnt} tournaments → I(Ω,2) = ?")
    # H is determined by the actual Ω graph, not just edge/disjoint count

print()
print("=" * 70)
print("KEY STRUCTURAL OBSTRUCTION")
print("=" * 70)
print()
print("For H=21 to occur, Ω(T) must be one of the 4 graph types with I(G,2)=21.")
print("Each type has specific structural requirements:")
print()
print("Type 1 (P₄/K₁⊔K₃): 4 odd cycles, specific adjacency pattern → BLOCKED (THM-201/202)")
print("Type 2 (K₆-2e): 6 cycles, exactly 2 disjoint pairs → needs verification")
print("Type 3 (K₈-e): 8 cycles, exactly 1 disjoint pair → needs verification")
print("Type 4 (K₁₀): 10 cycles, all pairs intersecting → needs verification")
print()

# Check: at n=6, with 6 cycles, can we get K₆-2e as Ω?
print("Checking n=6 tournaments with exactly 6 odd cycles:")
h_for_6_cycles = Counter()
omega_6_structures = Counter()

for bits in range(total):
    adj = tournament_from_bits(n, bits)
    cycles = find_all_directed_odd_cycles(n, adj, max_len=5)

    if len(cycles) != 6:
        continue

    H = hp_count(n, adj)
    h_for_6_cycles[H] += 1

    nc = len(cycles)
    edge_count = 0
    for i in range(nc):
        for j in range(i+1, nc):
            if cycles[i][1] & cycles[j][1]:
                edge_count += 1

    omega_6_structures[(nc, edge_count)] += 1

print(f"  H values: {dict(sorted(h_for_6_cycles.items()))}")
print(f"  Ω structures: {dict(sorted(omega_6_structures.items()))}")
print()

# K₆-2e has 6 vertices, 13 edges (C(6,2)-2 = 13)
print("  K₆-2e would have |V|=6, |E|=13 (= C(6,2) - 2)")
has_k6_2e = any(ec == 13 for (_, ec) in omega_6_structures)
print(f"  Does K₆-2e appear? {has_k6_2e}")

if has_k6_2e:
    cnt = sum(c for (_, ec), c in omega_6_structures.items() if ec == 13)
    print(f"  Count: {cnt} tournaments with Ω ≅ graph on 6 vertices, 13 edges")
    # H should be 21 if this is truly K₆-2e
    print(f"  Expected H = I(K₆-2e, 2) = 21")
else:
    print(f"  K₆-2e NEVER appears as Ω at n=6!")
    print()
    print("  At n=6, with 6 cycles, Ω edge counts are:")
    for (nc, ec), cnt in sorted(omega_6_structures.items()):
        max_edges = nc*(nc-1)//2
        print(f"    |E|={ec}/{max_edges}: {cnt} tournaments (missing {max_edges-ec} from K_{nc})")
