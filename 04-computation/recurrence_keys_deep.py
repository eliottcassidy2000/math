#!/usr/bin/env python3
"""
recurrence_keys_deep.py — opus-2026-03-14-S75

Deep exploration: "2 and 3 are the keys to the universe"

The user's prompt calls for seeing EVERYTHING as recurrences, with:
- k-nacci → 2, weighted k-nacci → 3
- How 5,6 relate to 7=2²+3 and 8=2³
- 10,11 as shifts of 1
- Alternating sum non-negativity and simplex packing

This script explores the recurrence structure underlying I(x), H(T),
and the tournament polynomial z²-5z+6=0.
"""

from itertools import permutations, combinations
from math import factorial, gcd
from functools import lru_cache

print("=" * 70)
print("PART 1: THE FUNDAMENTAL RECURRENCE z² = 5z - 6")
print("=" * 70)
print()
print("  z² - 5z + 6 = (z-2)(z-3) = 0")
print()
print("  The general solution: a(n) = A·2ⁿ + B·3ⁿ")
print()
print("  Special sequences:")
print("  a(n) = 2ⁿ:     1, 2, 4,  8, 16,  32,  64, 128, 256, 512, 1024 ...")
print("  a(n) = 3ⁿ:     1, 3, 9, 27, 81, 243, 729 ...")
print("  a(n) = 3ⁿ-2ⁿ:  0, 1, 5, 19, 65, 211, 665 ...")
print("  a(n) = 2·3ⁿ-2ⁿ: 1, 4, 14, 46, 146, 454 ... (sum 3^k+2^k?)")
print()

# Check the recurrence for various sequences
def check_rec(name, seq):
    """Check if a(n) = 5a(n-1) - 6a(n-2)"""
    ok = all(seq[i] == 5*seq[i-1] - 6*seq[i-2] for i in range(2, len(seq)))
    print(f"  {name}: {seq[:8]}... satisfies z²=5z-6? {ok}")

check_rec("2^n", [2**n for n in range(10)])
check_rec("3^n", [3**n for n in range(10)])
check_rec("3^n-2^n", [3**n-2**n for n in range(10)])
check_rec("5^n", [5**n for n in range(10)])
check_rec("6^n", [6**n for n in range(10)])

print()
print("  5^n does NOT satisfy z²=5z-6. It satisfies z=5z-0? No.")
print("  5^n satisfies a(n)=5·a(n-1), a 1st-order recurrence.")
print("  6^n satisfies a(n)=6·a(n-1).")
print()

# The key: I(x) as a function of α₁ (with α₂ etc fixed)
print("  I(x) = 1 + α₁·x + α₂·x² + α₃·x³ + ...")
print("  Fixing α₂,α₃,..., I(x) is LINEAR in α₁.")
print("  So I(2) and I(3) are linear in α₁:")
print("    I(2) = (1+4α₂+8α₃+...) + 2α₁")
print("    I(3) = (1+9α₂+27α₃+...) + 3α₁")
print()
print("  The RATIO I(3)/I(2) approaches 3/2 as α₁ → ∞")
print("  because the linear terms dominate: 3α₁/(2α₁) = 3/2")
print()

print("=" * 70)
print("PART 2: k-NACCI SEQUENCES AND THE KEY 2")
print("=" * 70)
print()
print("  FIBONACCI: a(n) = a(n-1) + a(n-2)")
print("    Ratio → φ = (1+√5)/2 ≈ 1.618")
print()
print("  TRIBONACCI: a(n) = a(n-1) + a(n-2) + a(n-3)")
print("    Ratio → ≈ 1.839")
print()
print("  k-NACCI: a(n) = a(n-1) + a(n-2) + ... + a(n-k)")
print("    Ratio → 2 as k → ∞")
print()

# Compute k-nacci ratios
for k in range(2, 15):
    seq = [0]*k + [1]
    for i in range(100):
        seq.append(sum(seq[-(k):]))
    ratio = seq[-1]/seq[-2] if seq[-2] > 0 else 0
    err = abs(ratio - 2)
    halfk = (1/2)**k
    print(f"  k={k:2d}: ratio = {ratio:.10f}, |ratio-2| = {err:.2e}, (1/2)^k = {halfk:.2e}, err/(1/2)^k = {err/halfk:.4f}" if halfk > 0 else "")

print()
print("  The error |ratio_k - 2| ≈ C · (1/2)^k")
print("  k-nacci approaches the KEY 2 at rate (1/2)^k")
print()

print("=" * 70)
print("PART 3: WEIGHTED k-NACCI AND THE KEY 3")
print("=" * 70)
print()
print("  WEIGHTED k-NACCI: a(n) = 2·a(n-1) + 2·a(n-2) + ... + 2·a(n-k)")
print("  Wait, what weight gives 3?")
print()
print("  Consider: a(n) = c·[a(n-1) + a(n-2) + ... + a(n-k)]")
print("  The ratio r satisfies r^k = c·(r^{k-1} + ... + r + 1) = c·(r^k-1)/(r-1)")
print("  So r^k(r-1) = c(r^k-1)")
print("  As k→∞: r→c+1 (since r^k dominates)")
print()

# For c=2: ratio → 3
print("  c=1: standard k-nacci, ratio → 2 = 1+1")
print("  c=2: weighted k-nacci, ratio → 3 = 2+1")
print("  c=m: ratio → m+1")
print()
print("  The KEY 2 is (1+1)-nacci, the KEY 3 is (2+1)-nacci!")
print()

# Verify c=2 weighted k-nacci → 3
for k in range(2, 15):
    seq = [0]*k + [1]
    for i in range(100):
        seq.append(2 * sum(seq[-(k):]))
    ratio = seq[-1]/seq[-2] if seq[-2] > 0 else 0
    err = abs(ratio - 3)
    twothirds_k = (2/3)**k
    print(f"  k={k:2d}: ratio = {ratio:.10f}, |ratio-3| = {err:.2e}, (2/3)^k = {twothirds_k:.2e}, err/(2/3)^k = {err/twothirds_k:.4f}" if twothirds_k > 0 else "")

print()
print("  The error |ratio_k - 3| ≈ C · (2/3)^k")
print("  Weighted k-nacci approaches the KEY 3 at rate (2/3)^k")
print()

print("=" * 70)
print("PART 4: THE CONVERGENCE RATES AND TOURNAMENT STRUCTURE")
print("=" * 70)
print()
print("  k-nacci → 2 at rate (1/2)^k")
print("  weighted k-nacci → 3 at rate (2/3)^k")
print()
print("  Note: 1/2 and 2/3 are the RECIPROCALS of the keys!")
print("  The convergence rate of each key is the reciprocal of the other:")
print("    key 2: convergence rate = 1/3... no wait")
print("    rate to 2: (1/2)^k = inverse of key 2")
print("    rate to 3: (2/3)^k = (key 2)/(key 3)... the ratio of keys!")
print()
print("  (1/2)^k: the binary shrinkage rate")
print("  (2/3)^k: the ternary-vs-binary residual rate")
print()
print("  TOURNAMENT CONNECTION:")
print("  α_k decays roughly as (1/2)^k of α_1 in the 'generic' case.")
print("  This means the higher-order terms in I(x) = 1 + α₁x + α₂x² + ...")
print("  decay geometrically, making I(x) CONVERGENT for |x| < 2.")
print()
print("  The radius of convergence of the 'normalized' independence")
print("  polynomial is related to the inverse of the max eigenvalue")
print("  of the conflict graph, which is bounded by α₁.")
print()

print("=" * 70)
print("PART 5: 5 AND 6 — SUM AND PRODUCT OF KEYS")
print("=" * 70)
print()
print("  5 = 2 + 3 (sum of keys)")
print("  6 = 2 · 3 (product of keys)")
print("  The characteristic polynomial: z² - 5z + 6 = 0")
print()
print("  5 in tournament theory:")
print("    - n=5 is the FIRST non-trivial tournament size")
print("    - At n=5: H ∈ {1,3,5,9,11,13,15} — exactly 7 values")
print("    - The max H at n=5 is 15 = 2⁴-1 = 2^{n-1}-1")
print("    - I(5) = 1 + 5α₁ + 25α₂ + ... (discriminant Δ=1!)")
print()
print("  6 in tournament theory:")
print("    - n=6 is where α₂ first appears (disjoint cycle pairs)")
print("    - At n=6: H does NOT determine (α₁,α₂)")
print("    - 6 = 3! = volume of 3-cube in simplex units")
print("    - The product 2·3=6 is the 'unit cell' of the tournament lattice")
print()

# Compute I(5) and I(6) for n=5 tournaments
print("  I(x) values at n=5 (α₂=0, I(x)=1+α₁·x):")
print(f"  {'H':>4s} {'α₁':>4s} {'I(5)':>6s} {'I(6)':>6s} {'I(5)mod6':>8s} {'I(6)mod5':>8s}")
for h in [1,3,5,9,11,13,15]:
    a1 = (h-1)//2
    i5 = 1 + 5*a1
    i6 = 1 + 6*a1
    print(f"  {h:4d} {a1:4d} {i5:6d} {i6:6d} {i5%6:8d} {i6%5:8d}")

print()
print("  I(5) mod 6: always 1. Because I(5)=1+5α₁, and 5≡-1(mod6),")
print("  so I(5)=1+(-1)^1·α₁ = 1-α₁ mod 6... no.")
print("  Actually: I(5)=1+5α₁. 5α₁ mod 6: depends on α₁ mod 6.")
print("  Hmm, let me check: 1+5·0=1, 1+5·1=6, 1+5·2=11, 1+5·3=16, 1+5·4=21, 1+5·5=26, 1+5·6=31, 1+5·7=36")
vals = [1+5*a for a in range(8)]
mods = [v%6 for v in vals]
print(f"  I(5) mod 6 for α₁=0..7: {mods}")
print(f"  Pattern: {mods[:6]} repeats with period 6")
print()

print("  I(6) mod 5: always 1! Because I(6)=1+6α₁, and 6≡1(mod5),")
print("  so I(6)≡1+α₁ (mod5). For α₁=0,1,...,6: {}", [((1+a)%5) for a in range(7)])
mods5 = [(1+6*a)%5 for a in range(8)]
print(f"  Actually I(6) mod 5 for α₁=0..7: {mods5}")
print(f"  Hmm: 1+6α₁ mod 5 = 1+α₁ mod 5: {[(1+a)%5 for a in range(8)]}")
print("  So I(6)≡1+α₁ (mod 5), NOT always 1.")
print()

print("  KEY INSIGHT: I(x) mod (x+1) and I(x) mod (x-1):")
print("    I(x) mod (x-1): substitute x≡1 → I(1) = 1+α₁+α₂+... = total indep sets")
print("    I(x) mod (x+1): substitute x≡-1 → I(-1) = 1-α₁+α₂-... = Euler char")
print()
print("  So: I(5) mod 6 = I(-1) mod 6 (since 5≡-1 mod 6)")
print("       I(5) mod 4 = I(1) mod 4 (since 5≡1 mod 4)")
print()

# Verify
print("  Verification:")
for h in [1,3,5,9,11,13,15]:
    a1 = (h-1)//2
    i5 = 1 + 5*a1
    im1 = 1 - a1  # I(-1)
    i1 = 1 + a1   # I(1)
    print(f"    H={h}: I(5)={i5}, I(5)mod6={i5%6}, I(-1)mod6={im1%6}  |  I(5)mod4={i5%4}, I(1)mod4={i1%4}")

print()

print("=" * 70)
print("PART 6: 7=2²+3 AND 8=2³ — THE NEXT LEVEL")
print("=" * 70)
print()
print("  7 = 2² + 3 = 4 + 3")
print("  8 = 2³ = 8")
print()
print("  7 in tournament theory:")
print("    - Φ₃(2) = 7 (3rd cyclotomic polynomial at 2)")
print("    - H=7 is FORBIDDEN (the 'Forbidden 3': α₁=3 impossible at n≤6)")
print("    - n=7: first n where dc3=3 with dc5=0 is possible")
print("    - 7 values of H at n=5: {1,3,5,9,11,13,15}")
print("    - The H spectrum at n=7 starts {1,3,5,9,11,13,15,...} = same as n=5!")
print()
print("  8 in tournament theory:")
print("    - 8 = 2³ = coefficient of α₃ in I(2)")
print("    - At n≥9: I(2) = 1 + 2α₁ + 4α₂ + 8α₃ + ...")
print("    - 8 is the weight of 'triple disjoint cycles'")
print("    - 8 = max(H) at n=4 (wait: max H at n=4 is 5, not 8)")
print("    - Hmm, 8 appears: at n=5, H∈{1,3,5,9,11,13,15} — 8 is NOT an H value!")
print("    - 8 = 2³ is a STRUCTURAL number, not a spectral number")
print()

# 7 = 2²+3 decomposition
print("  7 = 2² + 3:")
print("  In I(x): the coefficient of α₁ is x, of α₂ is x².")
print("  At x=2: contribution of (α₁=1, α₂=1) to H is 2+4 = 6")
print("  At x=3: contribution of (α₁=1, α₂=1) to I(3) is 3+9 = 12")
print("  7 = 2²+3 is the contribution at x=2 from α₂=1 AND α₁=1... wait")
print("  Actually 4+3 is just a number-theoretic decomposition.")
print()
print("  MORE PRECISELY:")
print("  7 = Φ₃(2) = 1+2+4 = 2³-1 = sum of powers of 2 up to 2²")
print("  This is 'all bits set in 3-bit binary': 111₂ = 7")
print("  In independence polynomial: 1 + 2·1 + 4·1 = 7 corresponds to")
print("  α₁=1, α₂=1 (one cycle plus one disjoint pair)")
print("  But wait: if you have a disjoint pair, you have at least 2 cycles,")
print("  so α₁ ≥ 2. So (α₁=1, α₂=1) is IMPOSSIBLE.")
print("  THIS is why 7 = 1+2+4 can't arise as H with α_k ∈ {0,1}!")
print()

# Explore: which H values have α₁=1+α₂ (minimal α₁ for given α₂)?
print("  For α₂ > 0, we need α₁ ≥ 2 (at least 2 cycles for a disjoint pair).")
print("  Actually: α₂ independent sets of size 2 need ≥ 2 cycles for EACH pair,")
print("  but different pairs can share cycles. So α₁ ≥ 2 when α₂ ≥ 1.")
print()
print("  If α₂=1: need ≥ 2 cycles. Minimum H = 1 + 2·2 + 4·1 = 9.")
print("  If α₂=2: need ≥ 3 cycles. Minimum H = 1 + 2·3 + 4·2 = 15.")
print("  If α₂=k: need ≥ (k+1) cycles. Min H = 1 + 2(k+1) + 4k = 6k+3.")
print()
for k in range(5):
    min_a1 = k + 1  # minimum cycles for k disjoint pairs
    min_H = 1 + 2*min_a1 + 4*k
    print(f"  α₂={k}: min α₁={min_a1}, min H = {min_H}")

print()

print("=" * 70)
print("PART 7: 10 AND 11 AS SHIFTED ONES")
print("=" * 70)
print()
print("  10 = 1·10¹ + 0·10⁰ = '10' in base 10")
print("  11 = 1·10¹ + 1·10⁰ = '11' in base 10")
print()
print("  In base 2: 10₂ = 2, 11₂ = 3 (THE KEYS!)")
print("  In base 3: 10₃ = 3, 11₃ = 4")
print()
print("  '10' in base b = b (the base itself)")
print("  '11' in base b = b+1")
print()
print("  So '10' and '11' are the UNIVERSAL representation of")
print("  (base, base+1) — the first 'non-trivial' pair.")
print()
print("  For the tournament recurrence z²-5z+6=0:")
print("  The roots are 2='10₂' and 3='11₂'")
print("  In base 10: 10 and 11 are similarly 'shifted ones'")
print()
print("  CYCLOTOMIC VALUES:")
print("    Φ₁(2) = 1   (trivial)")
print("    Φ₂(2) = 3   (key)")
print("    Φ₃(2) = 7   (forbidden)")
print("    Φ₄(2) = 5   (sum of keys)")
print("    Φ₅(2) = 31  (Mersenne prime)")
print("    Φ₆(2) = 3   (= Φ₂(2)!)")
print("   Φ₁₀(2) = 11  ('11' = shifted 1)")
print("   Φ₁₂(2) = 13")
print()
print("  Where does 10 appear?")
print("  Φ_d(2) = 10 for some d?")
# Check
from sympy import factorint
for d in range(1, 50):
    # Compute Φ_d(2) using the product formula
    # Φ_d(x) = prod_{k|d} (x^k - 1)^{μ(d/k)}
    pass

# Direct: 2^d - 1 = prod_{k|d} Φ_k(2)
# So Φ_d(2) = (2^d - 1) / prod_{k|d, k<d} Φ_k(2)
@lru_cache(maxsize=None)
def cyclotomic_at_2(d):
    """Compute Φ_d(2)"""
    val = 2**d - 1
    for k in range(1, d):
        if d % k == 0:
            val //= cyclotomic_at_2(k)
    return val

print("  Φ_d(2) for d=1..30:")
for d in range(1, 31):
    v = cyclotomic_at_2(d)
    mark = ""
    if v == 10: mark = " ← 10!"
    if v == 11: mark = " ← 11!"
    if v in [1,2,3,5,7]: mark = f" ← {['','1','KEY','KEY','','SUM','','FORB'][v]}" if v <= 7 else ""
    if v == 1: mark = " ← trivial"
    if v == 2: mark = " ← KEY"
    if v == 3: mark = " ← KEY"
    if v == 5: mark = " ← SUM of keys"
    if v == 7: mark = " ← FORBIDDEN"
    print(f"    Φ_{d:2d}(2) = {v}{mark}")

print()
print("  Φ_d(2) = 10 never appears! Because 10=2·5, and if Φ_d(2)=10,")
print("  then 2|Φ_d(2), but Φ_d(2) is always odd for d≥2 (since Φ_d(2)|(2^d-1) which is odd).")
print("  And Φ₁(2)=1. So 10 CANNOT be a cyclotomic value at 2.")
print()
print("  11 = Φ₁₀(2). Period 10 = 2·5 = product of (key 2) and (sum of keys 5).")
print()

print("=" * 70)
print("PART 8: THE ALTERNATING SUM AND SIMPLEX PACKING — PRECISE CONNECTION")
print("=" * 70)
print()
print("  I(-1) = 1 - α₁ + α₂ - α₃ + ...")
print()
print("  GEOMETRIC INTERPRETATION via Inclusion-Exclusion:")
print("  Each odd cycle C contributes a 'twist' to the simplex arrangement.")
print("  The SIGNED count of simplices is:")
print("    I(-1) = [empty set] - [single cycles] + [disjoint pairs] - [triples] + ...")
print("           = 1 - α₁ + α₂ - α₃ + ...")
print()
print("  This is EXACTLY the Euler characteristic of the independence complex!")
print()
print("  For SIMPLEX PACKING:")
print("  Given H simplices arranged in n-space, we can ask:")
print("  What is the 'effective dimension' of their arrangement?")
print()
print("  The Euler characteristic I(-1) measures the 'net algebraic count'")
print("  after canceling overlapping regions with alternating signs.")
print()
print("  HILBERT'S 3RD PROBLEM:")
print("  Dehn's invariant for the tournament polytope O(T) is determined by")
print("  the dihedral angles of the simplices.")
print("  In a SIMPLEX ARRANGEMENT, the dihedral angles are all π/2 or π/3")
print("  (since the simplices come from hyperplane arrangements x_σ(i)>x_σ(j)).")
print()
print("  The dihedral angle between adjacent simplices (sharing a facet) is")
print("  determined by a SINGLE arc reversal in the tournament.")
print("  An arc reversal changes one inequality x_i > x_j to x_j > x_i.")
print()
print("  DEHN INVARIANT OF O(T):")
print("  D(O(T)) = Σ_{facets F} vol(F) ⊗ θ(F)")
print("  where θ(F) is the dihedral angle at facet F.")
print()
print("  For the ORDER POLYTOPE, every internal dihedral angle is π")
print("  (the facets between adjacent simplices are flat).")
print("  So the Dehn invariant only has contributions from BOUNDARY facets.")
print()
print("  BOUNDARY FACETS of O(T):")
print("  These are the facets of the n-cube: x_i = 0 or x_i = 1.")
print("  The dihedral angle at each boundary facet is π/2.")
print()
print("  D(O(T)) = Σ_{boundary facets} vol(F) ⊗ π/2")
print("          = (total boundary facet volume) ⊗ π/2")
print()
print("  Since π/2 is irrational/π (well, π/2 is 1/2 in units of π,")
print("  so π/2 ⊗_Q R/πQ = 1/2 which is NOT zero).")
print()
print("  The total boundary facet volume depends on H and the arrangement.")
print("  For the transitive tournament: O(T) = standard simplex,")
print("  boundary facet volume = (n+1)/(n-1)! (sum of n+1 facets).")
print()

print("=" * 70)
print("PART 9: THE 5-6 PATTERN IN TOURNAMENT SPECTRA")
print("=" * 70)

# Generate all tournaments at n=3,4,5,6 and collect H, I(-1)
def gen_tournaments(n):
    """Generate all tournaments on n vertices as adjacency tuples."""
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    results = []
    for bits in range(2**len(edges)):
        adj = [[False]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if bits & (1 << idx):
                adj[i][j] = True
            else:
                adj[j][i] = True
        results.append(adj)
    return results

def count_ham_paths(adj, n):
    """Count Hamiltonian paths via DP."""
    dp = [[0]*(1<<n) for _ in range(n)]
    for v in range(n):
        dp[v][1<<v] = 1
    for mask in range(1, 1<<n):
        for v in range(n):
            if not (mask & (1<<v)): continue
            if dp[v][mask] == 0: continue
            for u in range(n):
                if mask & (1<<u): continue
                if adj[v][u]:
                    dp[u][mask|(1<<u)] += dp[v][mask]
    return sum(dp[v][(1<<n)-1] for v in range(n))

def find_odd_cycles(adj, n):
    """Find all chordless odd cycles (3-cycles for small n)."""
    cycles = []
    # 3-cycles
    for i in range(n):
        for j in range(i+1,n):
            for k in range(j+1,n):
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    cycles.append(frozenset([i,j,k]))
                elif adj[j][i] and adj[i][k] and adj[k][j]:
                    cycles.append(frozenset([i,j,k]))
    # 5-cycles for n≥5
    if n >= 5:
        for combo in combinations(range(n), 5):
            verts = list(combo)
            for perm in permutations(verts):
                is_cycle = True
                for idx in range(5):
                    if not adj[perm[idx]][perm[(idx+1)%5]]:
                        is_cycle = False
                        break
                if is_cycle:
                    # Check chordless
                    chordless = True
                    for idx in range(5):
                        v1 = perm[idx]
                        v2 = perm[(idx+2)%5]
                        if adj[v1][v2]:
                            chordless = False
                            break
                    if chordless:
                        cycles.append(frozenset(combo))
                        break  # one direction suffices for detection
    return list(set(cycles))

def compute_alpha(adj, n):
    """Compute α₁, α₂ of independence polynomial of conflict graph."""
    cycles = find_odd_cycles(adj, n)
    alpha1 = len(cycles)
    # Count independent pairs (vertex-disjoint)
    alpha2 = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if len(cycles[i] & cycles[j]) == 0:
                alpha2 += 1
    return alpha1, alpha2

print()
print("  H spectrum pattern:")
for n in range(3, 7):
    tours = gen_tournaments(n)
    h_vals = set()
    for adj in tours:
        h = count_ham_paths(adj, n)
        h_vals.add(h)
    h_sorted = sorted(h_vals)
    print(f"    n={n}: {h_sorted}")
    # Check 5-6 pattern: differences
    diffs = [h_sorted[i+1]-h_sorted[i] for i in range(len(h_sorted)-1)]
    print(f"         diffs: {diffs}")

print()
print("  n=5 diffs: [2,2,4,2,2,2] — gaps of 2 and 4")
print("  n=6 diffs: more complex, but...")
print()

# The 5-6 pattern in H spectrum
print("  At n=5: H mod 5 values:")
tours5 = gen_tournaments(5)
h_mod5 = {}
for adj in tours5:
    h = count_ham_paths(adj, 5)
    r = h % 5
    h_mod5.setdefault(r, set()).add(h)
for r in sorted(h_mod5):
    print(f"    H ≡ {r} (mod 5): {sorted(h_mod5[r])}")

print()
print("  At n=5: H mod 6 values:")
h_mod6 = {}
for adj in tours5:
    h = count_ham_paths(adj, 5)
    r = h % 6
    h_mod6.setdefault(r, set()).add(h)
for r in sorted(h_mod6):
    print(f"    H ≡ {r} (mod 6): {sorted(h_mod6[r])}")

print()
print("  PATTERN: H mod 6 only takes values {1, 3, 5} at n=5 — always ODD mod 6.")
print("  Because H is always odd (Rédei) and H mod 3 alternates: 1,0,2,0,2,1,0")
print()

print("=" * 70)
print("PART 10: ALTERNATING SUM NON-NEGATIVITY AND KRUSKAL-KATONA")
print("=" * 70)
print()
print("  The claim: α₁ - α₂ + α₃ - ... ≥ 0 for ALL tournaments.")
print("  Equivalently: I(-1) ≤ 1.")
print()
print("  This is a constraint on the f-vector of the independence complex")
print("  of the conflict graph CG(T).")
print()
print("  The Kruskal-Katona theorem constrains f-vectors of simplicial complexes.")
print("  But α₁ ≥ α₂ is STRONGER than what KK gives in general.")
print()
print("  WHY should α₁ ≥ α₂ hold for TOURNAMENT conflict graphs?")
print("  Not all graphs G have this property!")
print()

# Counterexample: complete bipartite graph K_{3,3}
# Independent sets: singles and pairs from same side
# Actually for K_{m,m}: α₁ = 2m vertices, α₂ = C(m,2)*2 + m*m pairs...
# Let me think about this differently.
# For the conflict graph CG(T), vertices = odd cycles, edges = sharing a vertex.

# Key property: CG(T) is a CLIQUE on all cycles through any vertex v.
# This constrains the structure heavily!

print("  KEY STRUCTURAL PROPERTY of CG(T):")
print("  All odd cycles through vertex v form a CLIQUE in CG(T).")
print("  (Because any two cycles sharing v share a vertex.)")
print()
print("  This means: the NEIGHBORHOOD of every vertex in CG(T)")
print("  contains a large clique (all cycles through any shared vertex).")
print()
print("  CONSEQUENCE: The independence number is bounded.")
print("  An independent set in CG(T) can contain at most ONE cycle")
print("  through any given vertex v.")
print()
print("  With n vertices in the tournament, an independent set of k")
print("  disjoint cycles uses at least 3k vertices (each cycle has ≥3 vertices).")
print("  So α_k = 0 for k > ⌊n/3⌋.")
print()
print("  More precisely: a set of k disjoint odd cycles uses")
print("  at least 3k vertices (if all are 3-cycles) and at most n vertices.")
print("  So k ≤ ⌊n/3⌋.")
print()

# The α₁ ≥ α₂ argument:
print("  WHY α₁ ≥ α₂:")
print("  Each disjoint PAIR {C₁, C₂} contributes to α₂.")
print("  But C₁ and C₂ individually contribute to α₁.")
print("  So α₂ ≤ C(α₁, 2) trivially.")
print("  But we need the stronger α₁ ≥ α₂.")
print()
print("  CLAIM: In any graph G where every vertex has a neighborhood")
print("  containing a large clique, α₁(G) ≥ α₂(G).")
print("  (Not obvious — this is specific to tournament CGs.)")
print()

# Direct verification at n=6
print("  VERIFICATION at n=6:")
n = 6
tours6 = gen_tournaments(n)
violations = 0
max_ratio = 0
for adj in tours6:
    a1, a2 = compute_alpha(adj, n)
    if a2 > a1:
        violations += 1
    if a1 > 0:
        max_ratio = max(max_ratio, a2/a1)

print(f"    Checked {len(tours6)} tournaments")
print(f"    α₂ > α₁ violations: {violations}")
print(f"    Max ratio α₂/α₁: {max_ratio:.4f}")
print()

print("=" * 70)
print("PART 11: THE RECURRENCE VIEW OF EVERYTHING")
print("=" * 70)
print()
print("  THESIS: Every key quantity in tournament theory satisfies")
print("  a recurrence with characteristic polynomial (z-2)(z-3) = z²-5z+6.")
print()
print("  1. H(T) = I(2) = 1 + 2α₁ + 4α₂ + 8α₃ + ...")
print("     View as: H = Σ 2^k α_k (evaluation at z=2)")
print()
print("  2. I(3) = 1 + 3α₁ + 9α₂ + 27α₃ + ...")
print("     View as: I(3) = Σ 3^k α_k (evaluation at z=3)")
print()
print("  3. For FIXED α-sequence, the function f(x) = Σ x^k α_k satisfies:")
print("     f(x) at x=2: H")
print("     f(x) at x=3: I(3)")
print("     The RATIO f(3)/f(2) → 3/2 generically")
print()
print("  4. The DIFFERENCE I(3) - I(2) = α₁ + 5α₂ + 19α₃ + ...")
print("     = Σ (3^k - 2^k) α_k")
print("     3^k - 2^k satisfies the SAME recurrence z²=5z-6!")
print()
print("  5. The SUM I(3) + I(2) = 2 + 5α₁ + 13α₂ + 35α₃ + ...")
print("     = 2 + Σ (3^k + 2^k) α_k")
print("     3^k + 2^k ALSO satisfies z²=5z-6!")
print()

# Verify
print("  Verification that 3^k+2^k satisfies a(n)=5a(n-1)-6a(n-2):")
seq = [3**k + 2**k for k in range(10)]
print(f"    Sequence: {seq}")
for i in range(2, 8):
    lhs = seq[i]
    rhs = 5*seq[i-1] - 6*seq[i-2]
    print(f"    a({i}) = {lhs}, 5a({i-1})-6a({i-2}) = {rhs}, match: {lhs==rhs}")

print()
print("  EVERY linear combination A·2^k + B·3^k satisfies z²=5z-6.")
print("  So EVERY evaluation I(x) is connected to I(2) and I(3) via this recurrence.")
print()

# The ultimate recurrence: I(x) satisfies a recurrence in x?
print("  Does I(x) satisfy a recurrence in x?")
print("  I(x) = 1 + α₁x + α₂x² is a POLYNOMIAL, not a recurrence.")
print("  But the COEFFICIENTS {1, α₁, α₂, ...} can be viewed as")
print("  the Hadamard product of {1, α₁, α₂, ...} with {1, x, x², ...}")
print()
print("  The generating function Σ α_k z^k has radius of convergence")
print("  related to the max eigenvalue of CG(T).")
print()

print("=" * 70)
print("PART 12: THE TRINITY — 2, 3, AND ∞")
print("=" * 70)
print()
print("  THREE REGIMES OF THE INDEPENDENCE POLYNOMIAL:")
print()
print("  I(-1) ← TOPOLOGY    (Euler characteristic, Dehn invariant)")
print("  I(1)  ← COUNTING    (total independent sets)")
print("  I(2)  ← GEOMETRY    (Hamiltonian paths, volume)")
print("  I(3)  ← TERNARY     (3-colorings of independent sets)")
print("  I(∞)  ← ASYMPTOTICS (independence number)")
print()
print("  The key RATIOS:")
print("    I(2)/I(1) = H/(1+α₁+α₂+...) ← 'binary amplification'")
print("    I(3)/I(2) = I(3)/H ← 'ternary/binary ratio → 3/2'")
print("    I(2)/I(-1) = H/(1-α₁+α₂-...) ← 'volume/topology ratio'")
print()

# Compute these ratios at n=5
print("  Ratios at n=5 (α₂=0):")
print(f"  {'H':>4s} {'α₁':>4s} {'I(-1)':>5s} {'I(1)':>4s} {'I(2)':>4s} {'I(3)':>4s} {'I2/I1':>7s} {'I3/I2':>7s}")
for h in [1,3,5,9,11,13,15]:
    a1 = (h-1)//2
    im1 = 1 - a1
    i1 = 1 + a1
    i2 = h
    i3 = 1 + 3*a1
    r21 = i2/i1 if i1 > 0 else float('inf')
    r32 = i3/i2 if i2 > 0 else float('inf')
    print(f"  {h:4d} {a1:4d} {im1:5d} {i1:4d} {i2:4d} {i3:4d} {r21:7.4f} {r32:7.4f}")

print()
print("  I(2)/I(1) → 2 as α₁ → ∞ (binary amplification → factor of 2)")
print("  I(3)/I(2) → 3/2 as α₁ → ∞ (ternary vs binary → 50% more)")
print("  These limits are the KEYS: 2 and 3/2 = 3/2.")
print()

print("=" * 70)
print("PART 13: SIMPLEX PACKING AND NON-NEGATIVITY — THE DEEP WHY")
print("=" * 70)
print()
print("  WHY is I(-1) ≤ 1 for tournament conflict graphs?")
print()
print("  SIMPLEX PACKING ARGUMENT:")
print("  O(T) = union of H simplices Δ_σ, one per Hamiltonian path σ.")
print("  These simplices form a TRIANGULATION of O(T).")
print("  The Euler characteristic of a triangulated polytope is 1")
print("  (convex polytope ≅ ball, χ(ball) = 1).")
print()
print("  BUT: χ(O(T)) = 1 is about the POLYTOPE, not the independence complex.")
print("  The independence complex Δ(CG(T)) has χ(Δ) = I(-1).")
print()
print("  Are they related? YES, via the NERVE THEOREM:")
print("  If the H simplices {Δ_σ} cover O(T) and all intersections are")
print("  contractible, then χ(O(T)) = χ(nerve of {Δ_σ}).")
print()
print("  The nerve of {Δ_σ} is the simplicial complex whose vertices")
print("  are the H Hamiltonian paths, and whose faces are sets of")
print("  paths with non-empty common intersection.")
print()
print("  But the independence complex is about the CYCLE STRUCTURE,")
print("  not the path structure. The connection is more subtle.")
print()
print("  TOPOLOGICAL ARGUMENT (sketch):")
print("  The conflict graph CG(T) has chromatic number ≤ n-2")
print("  (since the max independent set uses ≤ ⌊n/3⌋ cycles,")
print("  and the total number of cycles is bounded).")
print("  By Kruskal-Katona + chromatic number bounds,")
print("  the alternating sum α₁-α₂+... is controlled.")
print()
print("  DEEPER: The independence complex of CG(T) might be")
print("  SHELLABLE (admits a shelling order), which would imply")
print("  the h-vector is non-negative, and in particular χ ≤ 1.")
print()

# Check: is the independence complex shellable at n=5?
print("  The independence complex at n=5:")
print("  CG(T) has α₁ vertices and α₂ edges.")
print("  Since α₂=0 at n=5, the independence complex is")
print("  a disjoint union of α₁ points (0-dimensional).")
print("  Its Euler characteristic is α₁ (number of components).")
print("  Wait: χ(α₁ points) = α₁. But I(-1) = 1 - α₁.")
print("  Hmm, that's the REDUCED Euler char: χ̃ = χ - 1 = α₁ - 1.")
print("  And I(-1) = 1 + χ̃ = χ. No wait...")
print()
print("  Let me be careful:")
print("  Independence complex Δ = {S ⊆ V(CG): S is independent}")
print("  This includes the EMPTY SET as a face.")
print("  f_{-1} = 1 (empty face), f_0 = α₁ (vertices), f_1 = α₂ (edges), ...")
print("  χ(Δ) = f_{-1} - f_0 + f_1 - ... = 1 - α₁ + α₂ - ... = I(-1)")
print()
print("  So I(-1) = χ(Δ) (Euler characteristic of independence complex).")
print("  For a contractible complex: χ = 1.")
print("  I(-1) ≤ 1 means χ(Δ) ≤ 1.")
print()
print("  At n=5: Δ = {∅} ∪ {{c}: c ∈ cycles} (since α₂=0)")
print("  This is a discrete set of α₁ points plus the empty face.")
print("  χ = 1 - α₁. So I(-1) = 1 - α₁ ∈ {1, 0, -1, -3, -4, -5, -6}")
print()
print("  At n=6: Δ has edges too (disjoint pairs = α₂).")
print("  χ = 1 - α₁ + α₂.")

print()
print("=" * 70)
print("PART 14: THE 2-3 BRIDGE VIA BERNOULLI NUMBERS")
print("=" * 70)
print()
print("  The Bernoulli numbers B_n satisfy:")
print("    Σ B_n x^n/n! = x/(e^x - 1)")
print()
print("  At x=1: Σ B_n/n! = 1/(e-1)")
print("  At x=log2: involves powers of 2")
print("  At x=log3: involves powers of 3")
print()
print("  The Euler numbers E_n (not Bernoulli) satisfy:")
print("    Σ E_n x^n/n! = 1/cosh(x) = 2/(e^x + e^{-x})")
print()
print("  Alternating permutations are counted by Euler numbers.")
print("  TANGENT NUMBERS = E_{2n+1} (odd-index Euler numbers)")
print()
print("  CONNECTION: The deformed Eulerian numbers a_k(T) = A(n,k) + correction.")
print("  The correction depends on the independence polynomial!")
print("  So the 2-3 recurrence acts on the CORRECTION to the Eulerian distribution.")
print()

# Eulerian numbers at n=5
print("  Eulerian numbers A(5,k) for k=0..4:")
# A(5,0)=1, A(5,1)=26, A(5,2)=66, A(5,3)=26, A(5,4)=1
eulerian = [1, 26, 66, 26, 1]
print(f"    A(5,k) = {eulerian}")
print(f"    Sum = {sum(eulerian)} = 5! = {factorial(5)}")
print()
print("  For a tournament T with α₁ odd cycles:")
print("    a_k(T) = A(5,k) + correction_k")
print("    The corrections sum to 0 (total is always n!)")
print("    The corrections are controlled by α₁")
print()

print()
print("=" * 70)
print("DONE — Keys to the Universe Summary")
print("=" * 70)
print()
print("  2 and 3 are the roots of z² - 5z + 6 = 0")
print("  5 = 2+3 (sum, first non-trivial n)")
print("  6 = 2·3 (product, first α₂ appears)")
print("  7 = 2²+3 = Φ₃(2) (forbidden H value)")
print("  8 = 2³ (weight of triple cycles in I(2))")
print("  10 = NOT a cyclotomic value (even = impossible)")
print("  11 = Φ₁₀(2) (period 10 = 2·5)")
print()
print("  k-nacci → 2 at rate (1/2)^k")
print("  weighted(2) k-nacci → 3 at rate (2/3)^k")
print("  Convergence rates are key/key: 1/2 and 2/3")
print()
print("  I(-1) = Euler characteristic ≤ 1 (topology)")
print("  I(1) = total independent sets (counting)")
print("  I(2) = H = Hamiltonian paths (geometry)")
print("  I(3) = ternary evaluation (the other key)")
print()
print("  Alternating sum non-negativity ⟺ α₁ ≥ α₂")
print("  ⟺ simplex dominance in independence complex")
print("  ⟺ tournament polytope has bounded topological complexity")
print()
