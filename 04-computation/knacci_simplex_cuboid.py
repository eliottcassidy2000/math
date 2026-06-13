#!/usr/bin/env python3
"""
knacci_simplex_cuboid.py — opus-2026-03-14-S71g

Exploring the user's insight: k-nacci → 2, weighted k-nacci → 3,
simplices as (x+1)^n, cuboids as (x+2)^n, and packing them inside each other.

CONNECTIONS TO OCF:
- H(T) = I(Ω(T), 2) — evaluated at x=2
- Digraph world: each arc independent → (1+x)^m simplex
- Tournament world: paired arcs (x + (1-x)) → (1+2x) per pair → (1+2x)^{m/2}
  at x=1: (1+2)^{m/2} = 3^{m/2} — the cuboid!
- The k-nacci sequence: 1,1,2,4,8,16,... (k=2: Fibonacci→φ, k→∞: ratio→2)
- Weighted k-nacci: ratio → 3 (the tournament evaluation point!)

The "packing simplices inside cuboids" metaphor:
- A tournament is a constrained digraph
- The independence polynomial of a digraph is a "simplex" polynomial
- The tournament constraint converts it to a "cuboid" polynomial
- Packing = the constraint reducing degrees of freedom
"""

import numpy as np
from math import comb, factorial

# ============================================================
# Part 1: k-nacci sequences and their limits
# ============================================================
print("=" * 60)
print("K-NACCI SEQUENCES AND LIMIT RATIOS")
print("=" * 60)

for k in [2, 3, 4, 5, 6, 7, 10, 20, 50]:
    # k-nacci: a[n] = a[n-1] + a[n-2] + ... + a[n-k]
    # Initial: a[0]=1, a[1]=1, a[2]=2, ..., a[k-1]=2^{k-2}
    # Actually, standard k-nacci starts with k-1 zeros and one 1.
    # Let me use the "tribonacci-style" generalization.
    a = [0] * (k-1) + [1]  # k-1 zeros followed by a 1
    for _ in range(100):
        a.append(sum(a[-k:]))
    ratio = a[-1] / a[-2] if a[-2] > 0 else float('inf')
    print(f"  k={k:2d}: limit ratio = {ratio:.10f}")

print(f"\n  As k→∞, the k-nacci ratio → 2.")
print(f"  This is because a[n] ≈ a[n-1] + a[n-2] + ... + a[n-k]")
print(f"  ≈ a[n-1] + a[n-1]/r + a[n-1]/r² + ... ≈ a[n-1] · r/(r-1)")
print(f"  So r = r/(r-1) → r²-r=r → r²=2r → r=2.")

# ============================================================
# Part 2: Weighted k-nacci and the limit 3
# ============================================================
print(f"\n{'='*60}")
print("WEIGHTED K-NACCI: limit ratio → 3")
print(f"{'='*60}")

# Weighted k-nacci: a[n] = 2*a[n-1] + 2*a[n-2] + ... + 2*a[n-k]
# (weight 2 per term)
for k in [2, 3, 5, 10, 20]:
    a = [0] * (k-1) + [1]
    for _ in range(100):
        a.append(2 * sum(a[-k:]))
    ratio = a[-1] / a[-2] if a[-2] > 0 else float('inf')
    print(f"  k={k:2d}, weight=2: limit ratio = {ratio:.10f}")

print(f"\n  With weight w=2: r = w·r/(r-1) → r(r-1) = 2r → r² = 3r → r = 3!")
print(f"  This is the TOURNAMENT evaluation point!")

# General: weight w gives r = w+1
for w in [1, 2, 3, 4]:
    for k in [20]:
        a = [0] * (k-1) + [1]
        for _ in range(100):
            a.append(w * sum(a[-k:]))
        ratio = a[-1] / a[-2] if a[-2] > 0 else float('inf')
    print(f"  weight={w}: limit ratio = {ratio:.6f} (expected {w+1})")

# ============================================================
# Part 3: Simplex (x+1)^n and Cuboid (x+2)^n
# ============================================================
print(f"\n{'='*60}")
print("SIMPLEX (x+1)^n vs CUBOID (x+2)^n")
print(f"{'='*60}")

print("""
DIGRAPH:
  Each pair (i,j) has INDEPENDENT arcs: xᵢⱼ, xⱼᵢ ∈ {0,1}
  Per pair: (1+xᵢⱼ)(1+xⱼᵢ) = 4 states (no arc, i→j, j→i, both)
  Independence polynomial per pair: 1 + x (one arc contributes x to weight)
  Total: ∏ (1+x) = (1+x)^{n(n-1)} = SIMPLEX-type

TOURNAMENT:
  Each pair has CONSTRAINED arcs: xᵢⱼ + xⱼᵢ = 1
  Per pair: 2 states (i→j or j→i)
  For odd cycles: each pair contributes (x for forward) or (x for backward)
  The constraint couples the two directions.

  Key: I(Ω(T), x) at x=2.
  The "2" comes from: each odd cycle in an independent set can be traversed
  in 1 direction (since tournament determines direction), and the factor 2
  accounts for... what exactly?
""")

# The (x+1)^n simplex at x=1: 2^n (each vertex present or absent)
# The (x+2)^n cuboid at x=1: 3^n (each vertex has 3 states)
# The OCF: I(Ω, 2) = sum of 2^|S| over independent sets S

# For a graph with no edges (all cycles independent): I(G, x) = (1+x)^|V|
# At x=2: I = 3^|V| — the cuboid!
# For K_n (all edges): I(K_n, x) = 1 + nx
# At x=2: I = 1+2n — linear

print("Independence polynomial of edgeless graph on m vertices:")
for m in range(6):
    val = 3**m
    print(f"  I(K̄_{m}, 2) = (1+2)^{m} = 3^{m} = {val}")

print(f"\nFor edgeless Ω (no conflicts): H = 3^(#cycles)")
print(f"  This is the CUBOID evaluation!")
print(f"  Each independent cycle contributes a factor of 3.")

# ============================================================
# Part 4: Packing simplices inside cuboids
# ============================================================
print(f"\n{'='*60}")
print("PACKING: SIMPLICES INSIDE CUBOIDS")
print(f"{'='*60}")

print("""
The "packing" metaphor:

1. DIGRAPH → TOURNAMENT (adding constraint xⱼᵢ = 1-xᵢⱼ):
   - Reduces the state space from (1+x)^{n(n-1)} to something smaller
   - The constraint "packs" the simplex into a cuboid

2. Numerically:
   - Simplex: (x+1)^m, at x=1 gives 2^m choices
   - Cuboid: (x+2)^{m/2}, at x=1 gives 3^{m/2} choices
   - Tournament has m/2 = C(n,2)/2 paired edges, but each has 2 (not 3) states
   - The "3" appears in the independence polynomial: each independent cycle
     contributes factor 3 = 1+2 to I(Ω, 2)

3. The k-nacci connection:
   - k-nacci ratio → 2 as k→∞ (digraph world)
   - Weighted k-nacci ratio → 3 as k→∞ (tournament world)
   - The weight=2 is the OCF evaluation point!

4. Why (x+2)^n = cuboid?
   (x+2)^n = sum_k C(n,k) x^k 2^{n-k}
   Each position has 3 states: included-in-x (contributes x), or one of 2 "base" states.
   For tournaments: each arc pair has arc i→j (contributes to cycles in one way)
   or j→i (contributes differently), and the "2" counts the 2 options.
""")

# ============================================================
# Part 5: Concrete example — n=5 tournament
# ============================================================
print(f"{'='*60}")
print("EXAMPLE: n=5 tournament, simplex vs cuboid")
print(f"{'='*60}")

# Regular tournament on 5 vertices (all scores = 2)
# H = 15 (maximum at n=5)
# Ω has 5 three-cycles and 5 five-cycles (all pairwise sharing vertices since n=5)
# I(Ω, 2) = 1 + 2*10 = 21? Let me check.

# Actually, let me compute for a specific tournament
def make_tournament_bits(bits, n):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A

def count_hp(A, n):
    full = (1 << n) - 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)) or dp[S][v] == 0: continue
            for u in range(n):
                if not (S & (1 << u)) and A[v][u]:
                    dp[S | (1 << u)][u] += dp[S][v]
    return sum(dp[full][v] for v in range(n))

n = 5
# Transitive: H=1, Ω empty
# Regular: H=15, Ω has many cycles
# Let's look at H = I(Ω,2) values
h_counts = {}
for bits in range(1024):
    A = make_tournament_bits(bits, n)
    H = count_hp(A, n)
    h_counts[H] = h_counts.get(H, 0) + 1

print(f"\n  H values at n=5: {sorted(h_counts.keys())}")
print(f"  H=1 (transitive, Ω=empty): I(Ø,2) = 1 = (1+2)^0 = 3^0")
print(f"  H=3 (one 3-cycle, Ω=K₁): I(K₁,2) = 1+2 = 3 = 3^1")
print(f"  H=9 (Ω=?): I=9 = 1+8=1+2·4 (4 cycles, all pairwise sharing)")
print(f"     or I=9 = 1+2·1+4·... or I=9 = 3^2 (2 indep cycles)")

# For H=9: 9 = 3^2 if Ω has 2 vertex-disjoint cycles (independent)
# 9 = 1 + 2*4 if Ω=K₄ (4 pairwise sharing cycles)
# Both give I=9! Let's check which structure appears.

for bits in range(1024):
    A = make_tournament_bits(bits, n)
    H = count_hp(A, n)
    if H != 9:
        continue
    t3 = 0
    for a in range(5):
        for b in range(a+1,5):
            for c in range(b+1,5):
                if A[a][b] and A[b][c] and A[c][a]: t3 += 1
                if A[a][c] and A[c][b] and A[b][a]: t3 += 1
    print(f"  bits={bits}: H=9, t₃={t3}")
    break

# ============================================================
# Part 6: The (2,3) duality
# ============================================================
print(f"\n{'='*60}")
print("THE (2,3) DUALITY: GRAPH vs TOURNAMENT")
print(f"{'='*60}")

print("""
  GRAPH/DIGRAPH (simplex = (x+1)^n):
    - n independent dimensions
    - Each dimension: present (x) or absent (1)
    - At x=1: 2^n states (binary cube)
    - k-nacci limit: 2

  TOURNAMENT (cuboid = (x+2)^n):
    - n paired dimensions (each pair constrained)
    - Each pair: state A (x) or state B (2-x)... hmm
    - Independent cycles contribute 3 = 1+2 each
    - k-nacci limit (weighted): 3

  The OCF evaluation point x=2 IS the cuboid constant:
    I(Ω, 2) = sum_k α_k · 2^k
    If Ω is edgeless: I = (1+2)^{|V|} = 3^{|V|}
    This is the number of independent sets weighted by 2^{size}.

  The Grinberg-Stanley theorem says this weighted count
  equals the number of Hamiltonian paths.

  PACKING METAPHOR:
    (x+1)^n is a simplex (each arc independent)
    (x+2)^{n/2} is a cuboid (each pair contributes factor x+2=4 when x=2)
    Wait: at x=2, (2+1)^n = 3^n (simplex) and (2+2)^{n/2} = 4^{n/2} = 2^n (cuboid)

    Actually, let me reconsider. The user said:
    "simplices as (x+1)^n and cuboids as (x+2)^n"

    At x=1: simplex = 2^n, cuboid = 3^n.
    At x=2: simplex = 3^n, cuboid = 4^n.

    The OCF evaluates at x=2.
    I(Ω, 2) = sum α_k 2^k.
    For independent Ω: I = 3^m (where m = #cycles).
    For complete Ω: I = 1 + 2m.

    The simplex (1+x)^m at x=2 = 3^m = I(edgeless, 2).
    The cuboid (2+x)^m at x=2 = 4^m.

    So the I polynomial for edgeless Ω IS the simplex at x=2!
    And the cuboid would be I with some additional structure.

  BETTER INTERPRETATION:
    H(T) = I(Ω, 2).
    For a tournament with m independent cycles: H = 3^m.
    For a digraph: H could be anything (arcs independent).
    Simplex volume = 1/n! × side^n. Cuboid volume = side₁ × ... × sideₙ.
    Simplex fits inside cuboid when each sideᵢ ≥ appropriate bound.

    The "packing" is: a tournament (with its constraints) has
    H value that's always a value of I(G, 2) for some graph G.
    The constraint PACKS the HP count into a structured polynomial form.
""")

# ============================================================
# Part 7: Powers of 3 as "cuboid" H values
# ============================================================
print(f"{'='*60}")
print("POWERS OF 3: THE CUBOID H VALUES")
print(f"{'='*60}")

print("""
  H = 3^k means Ω has k independent (vertex-disjoint) cycles with
  no other odd cycles. This is the "maximally cuboid" configuration.

  3^0 = 1: transitive tournament (no cycles) — TRIVIAL CUBOID
  3^1 = 3: one 3-cycle + transitive complement
  3^2 = 9: two disjoint 3-cycles + transitive rest (n≥6)
  3^3 = 27: three disjoint 3-cycles (n=9, verified!)
  3^k: k disjoint 3-cycles on 3k vertices + transitive on remaining

  These are the "block-transitive" tournaments!
  H(T(B₁,...,Bₖ)) = ∏ H(Bᵢ) = 3^k when each block is a 3-cycle.

  The H=7 impossibility shows: NOT all values ≡ 1 mod 2 are achievable.
  Specifically: 7 ∉ {I(G, 2) : G = Ω(T) for some tournament T}.
  This is because Ω(T) has structural constraints that graphs G don't.
""")

# What fraction of odd numbers ≤ 100 are achievable?
achievable = set()
for n in range(3, 7):  # n=3..6 exhaustive (fast)
    total_edges = n*(n-1)//2
    total_t = 2**total_edges
    for bits in range(total_t):
        A = make_tournament_bits(bits, n)
        achievable.add(count_hp(A, n))
# n=7: sample
import random
random.seed(42)
n = 7
for _ in range(50000):
    bits = random.randint(0, 2**21-1)
    A = make_tournament_bits(bits, n)
    achievable.add(count_hp(A, n))

# Add products
ach_list = sorted(achievable)
products = set()
for a in ach_list:
    for b in ach_list:
        if a*b <= 200:
            products.add(a*b)
achievable.update(products)
# One more round
for a in sorted(achievable):
    if a > 200: continue
    for b in [1,3,5,9,11,13,15]:
        if a*b <= 200:
            achievable.add(a*b)

odd_achievable = sorted(h for h in achievable if h % 2 == 1 and h <= 100)
odd_not = sorted(h for h in range(1, 101, 2) if h not in achievable)
print(f"\n  Odd H ≤ 100 achievable: {len(odd_achievable)}/50")
print(f"  Odd H ≤ 100 NOT achievable: {odd_not}")
print(f"  Permanently forbidden: {{7, 21}} (proved for 7, conjectured for 21)")

print(f"\n{'='*60}")
print("DONE")
print(f"{'='*60}")
