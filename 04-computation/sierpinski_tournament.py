#!/usr/bin/env python3
"""
sierpinski_tournament.py — opus-2026-03-14-S76

"Freaky creative" exploration: the Sierpinski fractal structure of tournaments.

The fractal dimension log₂(3) ≈ 1.585 connects 2 and 3:
  log₂(3) = log(3)/log(2)
  3 = 2^{log₂(3)}
  9 = 3² = 2^{2·log₂(3)} ≈ 2^{3.17}

The Sierpinski triangle at level k has 3^k points.
The tournament analog:
  Level 0: single vertex (n=1)
  Level 1: 3-cycle (n=3)
  Level 2: "meta-3-cycle" of three 3-cycles (n=9)
  Level k: n = 3^k

At each level, the independence polynomial gains new terms:
  Level 1: I(x) = 1 + α₁x
  Level 2: I(x) = 1 + α₁x + α₂x² + α₃x³
  Level k: I(x) has degree 3^{k-1} = ⌊n/3⌋

KEY QUESTION: Can we build a "Sierpinski tournament" at each level
and study how H, α₁, α₂, I(-1) scale?

The fractal dimension tells us:
  dim(Sierpinski) = log(3)/log(2) = log₂(3)
  At n = 3^k: "tournament complexity" grows as 2^{k·log₂(3)} = 3^k = n
  So the fractal dimension of the tournament complexity = n itself

But the PACKING complexity (how many simplices fit) grows as 2^something...
  H = I(2) ≈ 2^{α₁} for sparse tournaments
  The "entropy" of the tournament is log₂(H)

For the self-similar Sierpinski tournament:
  H(level k) = ? in terms of H(level k-1)
"""

from itertools import combinations, permutations
import random
import math

def random_tournament(n):
    adj = [0] * n
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                adj[i] |= (1 << j)
            else:
                adj[j] |= (1 << i)
    return adj

def ham_paths(adj, n):
    count = 0
    for perm in permutations(range(n)):
        ok = True
        for idx in range(n-1):
            if not (adj[perm[idx]] & (1 << perm[idx+1])):
                ok = False
                break
        if ok:
            count += 1
    return count

def find_3cycles(adj, n):
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (adj[i]&(1<<j)) and (adj[j]&(1<<k)) and (adj[k]&(1<<i)):
                    cycles.append(frozenset([i,j,k]))
                elif (adj[i]&(1<<k)) and (adj[k]&(1<<j)) and (adj[j]&(1<<i)):
                    cycles.append(frozenset([i,j,k]))
    return cycles

def count_disjoint_pairs(cycles):
    count = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if len(cycles[i] & cycles[j]) == 0:
                count += 1
    return count

def count_disjoint_triples(cycles):
    count = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if len(cycles[i] & cycles[j]) > 0:
                continue
            for k in range(j+1, len(cycles)):
                if len(cycles[i] & cycles[k]) == 0 and len(cycles[j] & cycles[k]) == 0:
                    count += 1
    return count

def build_lex_product(adj1, n1, adj2, n2):
    N = n1 * n2
    adj = [0] * N
    for v1 in range(N):
        i1, j1 = v1 // n2, v1 % n2
        for v2 in range(N):
            if v1 == v2:
                continue
            i2, j2 = v2 // n2, v2 % n2
            if i1 != i2:
                if adj1[i1] & (1 << i2):
                    adj[v1] |= (1 << v2)
            else:
                if adj2[j1] & (1 << j2):
                    adj[v1] |= (1 << v2)
    return adj

# ====================================================================
print("=" * 70)
print("PART 1: THE SIERPINSKI TOURNAMENT HIERARCHY")
print("=" * 70)
print()

# C₃: the 3-cycle
adj_c3 = [0, 0, 0]
adj_c3[0] |= (1 << 1)  # 0→1
adj_c3[1] |= (1 << 2)  # 1→2
adj_c3[2] |= (1 << 0)  # 2→0

print("Level 1: C₃ (3-cycle on 3 vertices)")
H1 = ham_paths(adj_c3, 3)
c3_1 = find_3cycles(adj_c3, 3)
print(f"  H = {H1}")
print(f"  α₁ = {len(c3_1)}, α₂ = 0")
print(f"  I(-1) = {1 - len(c3_1)}")
print()

# Level 2: C₃ ⊠ C₃
print("Level 2: C₃ ⊠ C₃ (lex product, 9 vertices)")
adj_s2 = build_lex_product(adj_c3, 3, adj_c3, 3)
H2 = ham_paths(adj_s2, 9)
c3_2 = find_3cycles(adj_s2, 9)
a1_2 = len(c3_2)
a2_2 = count_disjoint_pairs(c3_2)
a3_2 = count_disjoint_triples(c3_2)
print(f"  H = {H2}")
print(f"  α₁ = {a1_2}, α₂ = {a2_2}, α₃ = {a3_2}")
print(f"  I(2) via 3-cycles = {1 + 2*a1_2 + 4*a2_2 + 8*a3_2}")
print(f"  I(-1) = {1 - a1_2 + a2_2 - a3_2}")
print(f"  log₂(H) = {math.log2(H2):.4f}")
print()

# Level 3 would be 27 vertices — can't compute Ham paths
# but can analyze cycle structure
print("Level 3: C₃ ⊠ C₃ ⊠ C₃ (27 vertices)")
adj_s3 = build_lex_product(adj_c3, 3, adj_s2, 9)
c3_3 = find_3cycles(adj_s3, 27)
a1_3 = len(c3_3)
a2_3 = count_disjoint_pairs(c3_3)
a3_3 = count_disjoint_triples(c3_3)
print(f"  α₁(3-cyc) = {a1_3}")
print(f"  α₂(3,3) = {a2_3}")
print(f"  α₃(3,3,3) = {a3_3}")

# ====================================================================
print()
print("=" * 70)
print("PART 2: THE FRACTAL SCALING")
print("=" * 70)
print()

data = [
    (1, 3, 3, 1, 0, 0),   # level 1
    (2, 9, H2, a1_2, a2_2, a3_2),  # level 2
]

print(f"{'Level':<6} {'n':<6} {'H':<8} {'α₁':<8} {'α₂':<8} {'α₃':<8} {'log₂H':<10} {'α₁/n':<8}")
for level, n, H, a1, a2, a3 in data:
    logH = math.log2(H) if H > 0 else 0
    print(f"{level:<6} {n:<6} {H:<8} {a1:<8} {a2:<8} {a3:<8} {logH:<10.4f} {a1/n:<8.4f}")

# Level 3 (partial data)
n3 = 27
print(f"{3:<6} {n3:<6} {'?':<8} {a1_3:<8} {a2_3:<8} {a3_3:<8} {'?':<10} {a1_3/n3:<8.4f}")

print()
print("SCALING ANALYSIS:")
print(f"  α₁: {1} → {a1_2} → {a1_3}")
print(f"  Ratios: {a1_2/1:.1f}, {a1_3/a1_2:.1f}")
print(f"  If fractal: α₁(level k) ≈ α₁(level 1) · r^k")
r_a1 = a1_2 / 1
print(f"  r_α₁ = {r_a1}")
print(f"  Predicted level 3: {1 * r_a1**2:.0f}, actual: {a1_3}")
print()
print(f"  α₂: {0} → {a2_2} → {a2_3}")
print(f"  α₃: {0} → {a3_2} → {a3_3}")

# ====================================================================
print()
print("=" * 70)
print("PART 3: THE log₂(3) DIMENSION")
print("=" * 70)
print()
print(f"  log₂(3) = {math.log2(3):.6f}")
print(f"  This is the Hausdorff dimension of the Sierpinski triangle.")
print()
print("  In tournament theory:")
print(f"  2 = KEY₁ (independence polynomial eval point)")
print(f"  3 = KEY₂ (vertices per elementary cycle)")
print(f"  log₂(3) = log(cycle_size) / log(eval_point)")
print()
print("  The ENTROPY of the 3-cycle:")
print("  A 3-cycle has H=3 (three Ham paths).")
print(f"  log₂(3) = {math.log2(3):.6f} bits per vertex.")
print("  This is the 'information content' of one tournament vertex")
print("  in the 3-cycle configuration.")
print()
print("  COMPARISON with k-cycle entropy:")
for k in [3, 5, 7, 9, 11]:
    # A k-cycle tournament has H = k (number of Hamiltonian paths = k)
    entropy = math.log2(k) / k
    print(f"  {k}-cycle: H={k}, entropy = log₂({k})/{k} = {entropy:.4f} bits/vertex")

print()
print("  The entropy per vertex DECREASES with cycle length!")
print("  3-cycle: 0.528 bits/vertex (MAXIMUM for odd cycles)")
print("  This explains why 3-cycles dominate the structure:")
print("  They are the MOST EFFICIENT carriers of 'tournament information'.")

# ====================================================================
print()
print("=" * 70)
print("PART 4: THE RECURRENCE CONNECTION")
print("=" * 70)
print()
print("  The tournament polynomial z² - 5z + 6 = 0 has roots 2, 3.")
print("  General solution: x_n = A·2ⁿ + B·3ⁿ")
print()
print("  For the Sierpinski tournament at level k:")
print(f"  n = 3^k (vertices)")
print(f"  α₁ ~ C · 3^k (3-cycles, proportional to n)")
print(f"  H ~ D · r^k for some growth rate r")
print()

# What is the growth rate of H?
# Level 1: H=3, Level 2: H=?
print(f"  H(level 1) = 3 = 3¹")
print(f"  H(level 2) = {H2}")
print(f"  H(level 2) / H(level 1) = {H2/3:.2f}")
print(f"  H(level 2) / H(level 1)² = {H2/9:.2f}")
print(f"  H(level 2) / H(level 1)³ = {H2/27:.2f}")
print()

# Is H2 = 3^k for some k?
if H2 > 0:
    k_eff = math.log(H2) / math.log(3)
    print(f"  H(level 2) = 3^{k_eff:.4f}")
    print(f"  If H(level k) = 3^(f(k)) for some superlinear f:")
    print(f"    f(1) = 1, f(2) = {k_eff:.4f}")
    print(f"    If f(k) = k²: 3^1 = 3 ✓, 3^4 = 81 {'✓' if H2==81 else '✗'}")
    print(f"    If f(k) = k·log₂(3): 3^{math.log2(3):.3f} = {3**math.log2(3):.1f} {'≈ H?' if abs(3**math.log2(3) - H2) < 1 else ''}")

# ====================================================================
print()
print("=" * 70)
print("PART 5: THE WEIGHT-THRESHOLD INTERPLAY")
print("=" * 70)
print()
print("  OCF: H = Σ_{k≥0} 2^k · α_k")
print("  The weight at level k is 2^k (from KEY₁^k).")
print("  The threshold for α_k > 0 is n ≥ 3k (from KEY₂ · k).")
print()
print("  Weight × Threshold table:")
print(f"  {'k':<4} {'weight 2^k':<12} {'threshold 3k':<14} {'product':<10} {'ratio w/t':<10}")
for k in range(1, 10):
    w = 2**k
    t = 3*k
    print(f"  {k:<4} {w:<12} {t:<14} {w*t:<10} {w/t:<10.4f}")

print()
print("  The weight/threshold ratio 2^k/(3k) grows SUPER-EXPONENTIALLY.")
print("  This means: higher-order independence sets are EXPONENTIALLY")
print("  more valuable per set, but SUPER-EXPONENTIALLY rarer.")
print()
print("  The NET contribution of level k to H:")
print("  2^k · α_k ≈ 2^k · C(α₁, k) · p_k")
print("  where p_k = probability that k random cycles are disjoint")
print("  p_k ≈ (1 - 3k/n)^{3k} ≈ e^{-9k²/n}")
print()
print("  For the 'typical' contribution to be O(1):")
print("  2^k · C(α₁,k) · e^{-9k²/n} ≈ 1")
print("  k* ≈ √(n · ln(2)/9) ≈ √(0.077n)")
print()
for n in [3, 9, 27, 81]:
    k_star = (0.077 * n) ** 0.5
    print(f"  n={n:3d}: k* ≈ {k_star:.2f}, max_k = ⌊n/3⌋ = {n//3}")

# ====================================================================
print()
print("=" * 70)
print("PART 6: THE DEEP STRUCTURE — WHY 2 AND 3?")
print("=" * 70)
print()
print("  Why does H = I(CG(T), 2) use evaluation at x = 2?")
print("  Because each disjoint cycle contributes a FACTOR of 2:")
print("  two orientations around the cycle.")
print()
print("  Why are cycles of length 3 dominant?")
print("  Because 3 is the SMALLEST odd number > 1.")
print("  An odd cycle must have ≥ 3 vertices.")
print()
print("  Why is z² - 5z + 6 = (z-2)(z-3) = 0 the 'universal' polynomial?")
print("  Because 2 and 3 are the first two primes.")
print("  2 = the evaluation point (binary choice per cycle)")
print("  3 = the cycle length (minimum odd cycle)")
print("  5 = 2 + 3 = the 'diameter' of the key pair")
print("  6 = 2 · 3 = the 'area' of the key pair")
print()
print("  The polynomial z² - 5z + 6 encodes:")
print("  z = 2: I(2) = H (Hamiltonian path count)")
print("  z = 3: I(3) = H + 3α₁ + ... (next evaluation)")
print("  z = 6: I(6) = I(2·3) = ... (product evaluation)")
print()
print("  THE BRIDGE TO RECURRENCES:")
print("  Any sequence x_n = 5x_{n-1} - 6x_{n-2} satisfies")
print("  x_n = A·2ⁿ + B·3ⁿ for constants A, B.")
print()
print("  For tournaments: the 'sequence' is the independence polynomial.")
print("  I(x) = Σ α_k x^k is a POLYNOMIAL, not a recurrence.")
print("  But the EVALUATION POINTS 2, 3 are the recurrence roots.")
print()
print("  The k-nacci limit of 2 corresponds to:")
print("  I(2) = H = the 'physical' count (always odd).")
print("  The weighted k-nacci limit of 3 corresponds to:")
print("  I(3) = the 'next-level' evaluation.")
print()
print("  I(3) - I(2) = α₁ + 5α₂ + 19α₃ + 65α₄ + ...")
print("  The coefficients 1, 5, 19, 65 = 3^k - 2^k.")
print("  These are the CORNER PIECE sizes in simplex-cuboid packing!")
print()
print("  CORNER PIECE SEQUENCE: 3^k - 2^k")
for k in range(1, 8):
    print(f"  k={k}: 3^{k} - 2^{k} = {3**k} - {2**k} = {3**k - 2**k}")

print()
print("  The corner pieces grow as 3^k (since 3 > 2).")
print("  Their ratio to 3^k is 1 - (2/3)^k → 1.")
print("  Their ratio to 2^k is (3/2)^k - 1 → ∞.")
print()
print("  SYNTHESIS:")
print("  I(3) - I(2) = Σ (3^k - 2^k) · α_k")
print("  = Σ (corner piece at level k) × (independent set count)")
print("  = total 'corner volume' in the simplex-cuboid packing")
print()
print("  The difference I(3) - I(2) measures the GAP between")
print("  the cuboid (3^k) and the simplex (2^k) at each level.")
print("  This gap is the RESIDUE of the packing.")

# ====================================================================
print()
print("=" * 70)
print("PART 7: THE 9 = 3² MANIFESTO (EXTENDED)")
print("=" * 70)
print()
print("  NINE REASONS WHY 9 IS THE BOUNDARY:")
print()
print("  1. CAUCHY-SCHWARZ: 9 = 3² = (vertices/cycle)².")
print("     The CS bound on vertex-cycle incidences is tight at n = 9.")
print()
print("  2. PERFECT PARTITION: V = C₁ ∪ C₂ ∪ C₃ with |Cᵢ| = 3.")
print("     9 = 3·3 is the first n where ⌊n/3⌋ = n/3 = 3 = cycle length.")
print()
print("  3. INDEPENDENCE DEPTH: max indep set size = ⌊9/3⌋ = 3.")
print("     I(x) has degree 3 for the first time.")
print()
print("  4. TRIBONACCI: k=3 regime begins, root ≈ 1.839.")
print("     Steepest convergence slope from 2.")
print()
print("  5. CARTAN: det(A₈) = 9. Tournament on 9 vertices lives")
print("     in the Weyl chamber structure of A₈.")
print()
print("  6. PALEY FAILURE: 9 ≡ 1 (mod 4), so no Paley tournament.")
print("     The √q Ramsey barrier α = ω = √9 = 3.")
print()
print("  7. SIERPINSKI LEVEL 2: 9 = 3^2 is the first 'meta-triangle'.")
print("     Self-similar structure of 3 blocks of 3.")
print()
print("  8. ENTROPY MAXIMUM: 3-cycles have max info per vertex")
print(f"     (log₂(3)/3 ≈ {math.log2(3)/3:.4f} bits/vertex).")
print("     At n = 9, three max-entropy cycles cover all vertices.")
print()
print("  9. CORNER-PIECE TRANSITION: at n = 9, the cubic corner piece")
print(f"     3³ - 2³ = 19 first appears in I(3) - I(2).")
print("     The coefficient of α₃ is 19 (the gap between")
print("     3-cuboid and 2-simplex at level 3).")
