#!/usr/bin/env python3
"""
i_mod_cyclotomic.py — opus-2026-03-14-S74
Investigate the stunning result: I(7) ≡ 1 (mod 7) ALWAYS at n=5.

Is I(p) ≡ 1 (mod p) when p = Φ_d(2)?
What about I(x) ≡ 1 (mod x) in general?
Does this extend to n > 5?

This could be a THEOREM about independence polynomials of conflict graphs.
"""

import numpy as np
from itertools import combinations, permutations
from math import gcd, comb
from collections import Counter
import random

def enumerate_tournaments(n):
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    for bits in range(2**len(edges)):
        adj = [[0]*n for _ in range(n)]
        for idx, (i, j) in enumerate(edges):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        yield adj

def random_tournament(n):
    adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
    return adj

def count_hamiltonian_paths(adj, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))

def get_directed_odd_cycles(adj, n):
    cycles = []
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            for perm in permutations(verts):
                is_cycle = True
                for i in range(length):
                    if not adj[perm[i]][perm[(i+1) % length]]:
                        is_cycle = False
                        break
                if is_cycle:
                    min_idx = perm.index(min(perm))
                    canon = tuple(perm[min_idx:] + perm[:min_idx])
                    cycles.append(canon)
    return list(set(cycles))

def build_conflict_graph(cycles):
    nc = len(cycles)
    vsets = [set(c) for c in cycles]
    adj_cg = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if vsets[i] & vsets[j]:
                adj_cg[i][j] = True
                adj_cg[j][i] = True
    return adj_cg

def independence_polynomial_at(adj_cg, nc, x):
    total = 0
    for mask in range(1 << nc):
        verts = [i for i in range(nc) if mask & (1 << i)]
        is_indep = True
        for i in range(len(verts)):
            for j in range(i+1, len(verts)):
                if adj_cg[verts[i]][verts[j]]:
                    is_indep = False
                    break
            if not is_indep:
                break
        if is_indep:
            total += x ** len(verts)
    return total

# ====================================================================
# PART 1: I(x) mod x at n=5 — EXHAUSTIVE
# ====================================================================
print("=" * 70)
print("PART 1: I(x) mod x AT n=5 — ALL 1024 TOURNAMENTS")
print("=" * 70)

n = 5
x_values = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]

print(f"\n  Computing I(x) for x in {x_values}...")

all_data = []
for adj in enumerate_tournaments(n):
    H = count_hamiltonian_paths(adj, n)
    cycles = get_directed_odd_cycles(adj, n)
    alpha1 = len(cycles)
    cg = build_conflict_graph(cycles)
    nc = len(cycles)

    I_vals = {}
    for x in x_values:
        I_vals[x] = independence_polynomial_at(cg, nc, x)

    all_data.append((H, alpha1, I_vals))

print(f"  Done. Checking I(x) mod x:")
print(f"\n  {'x':>3}  {'I(x)≡? (mod x)':>20}  {'all ≡ 1?':>10}  notes")

for x in x_values:
    residues = Counter()
    for H, a1, I_vals in all_data:
        residues[I_vals[x] % x] += 1
    all_one = (len(residues) == 1 and 1 in residues)
    res_str = ", ".join(f"{r}:{c}" for r, c in sorted(residues.items()))
    notes = ""
    if all_one:
        notes = "✓ THEOREM CANDIDATE"
    print(f"  {x:>3}  {res_str:>40}  {'YES' if all_one else 'no':>10}  {notes}")

# ====================================================================
# PART 2: PROVE I(x) ≡ 1 (mod x) ALGEBRAICALLY
# ====================================================================
print("\n" + "=" * 70)
print("PART 2: ALGEBRAIC PROOF THAT I(x) ≡ 1 (mod x)")
print("=" * 70)

print("""
  I(x) = 1 + x·α₁ + x²·α₂ + x³·α₃ + ...
       = 1 + x(α₁ + xα₂ + x²α₃ + ...)
       ≡ 1 (mod x)

  THIS IS TRIVIALLY TRUE FOR ALL x AND ALL TOURNAMENTS.

  The independence polynomial I(Ω, x) = Σ_{k≥0} α_k · x^k
  always has constant term α₀ = 1 (the empty set).
  So I(x) - 1 is divisible by x.

  THEREFORE: I(x) ≡ 1 (mod x) is ALWAYS true.

  This is NOT specific to cyclotomic values — it holds for ALL x.
  The forbidden residues we saw in Part 6 of the previous script
  are about I(x) mod PRIMES p ≠ x, not I(x) mod x.
""")

# ====================================================================
# PART 3: THE REAL QUESTION — I(x) mod p for p ≠ x
# ====================================================================
print("=" * 70)
print("PART 3: THE REAL FORBIDDEN RESIDUES — I(x) mod p, p ≠ x")
print("=" * 70)

print("\n  From Part 6 of cyclotomic_tournament_deep.py:")
print("  I(5) mod 5: only residue 1 — trivially true (I(x)≡1 mod x)")
print("  I(7) mod 7: only residue 1 — trivially true!")
print("  I(5) mod 7: forbidden residue 2")
print("  I(7) mod 5: forbidden residue 2")
print("  I(8) mod 5: forbidden residue 0")
print("  I(6) mod 5: forbidden residue 4")
print()
print("  The REAL question: which (x, p) pairs have forbidden residues?")

print(f"\n  {'x':>3} {'mod p':>6}  {'forbidden residues':>30}  notes")
for x in [2, 3, 5, 6, 7, 8, 11]:
    for p in [3, 5, 7, 11, 13]:
        if p == x:
            continue
        residues = Counter()
        for H, a1, I_vals in all_data:
            residues[I_vals[x] % p] += 1
        missing = [r for r in range(p) if residues[r] == 0]
        if missing:
            notes = ""
            if x == 2 and p == 7:
                notes = "← H mod 7"
            print(f"  {x:>3}   mod {p:>2}  forbidden: {missing}  {notes}")

# ====================================================================
# PART 4: I(x) ≡ 1 (mod x) IMPLIES WHAT FOR H?
# ====================================================================
print("\n" + "=" * 70)
print("PART 4: CONSEQUENCES OF I(x) ≡ 1 (mod x)")
print("=" * 70)

print("""
  Since I(x) ≡ 1 (mod x) for ALL x:

  At x = 2: I(2) = H ≡ 1 (mod 2)  ← RÉDEI'S THEOREM!
  At x = 3: I(3) ≡ 1 (mod 3)
  At x = 5: I(5) ≡ 1 (mod 5)
  At x = p: I(p) ≡ 1 (mod p) for every prime p

  Rédei's theorem (H is odd) is the x=2 SPECIAL CASE of this
  trivial algebraic identity!

  But there's more. I(x) = 1 + x·Q(x) where Q(x) = α₁ + α₂x + ...
  So I(x)/x = 1/x + Q(x), which means:
    (I(x) - 1)/x = Q(x) = α₁ + α₂x + α₃x² + ...

  This is the "reduced" independence polynomial.
  Q(2) = (H-1)/2 = α₁ + 2α₂ + 4α₃ + ...
  Q(3) = (I(3)-1)/3 = α₁ + 3α₂ + 9α₃ + ...
""")

# Now check Q(x) mod small primes
print("  Q(x) = (I(x)-1)/x mod small primes at n=5:")
print(f"\n  {'x':>3} {'mod p':>6}  {'Q(x) forbidden residues':>30}")
for x in [2, 3, 5, 7]:
    for p in [3, 5, 7, 11]:
        residues = Counter()
        for H, a1, I_vals in all_data:
            Qx = (I_vals[x] - 1) // x
            residues[Qx % p] += 1
        missing = [r for r in range(p) if residues[r] == 0]
        if missing:
            print(f"  {x:>3}   mod {p:>2}  forbidden: {missing}")

# ====================================================================
# PART 5: EXTEND TO n=6 — SAMPLING
# ====================================================================
print("\n" + "=" * 70)
print("PART 5: I(x) MOD PRIMES AT n=6 — EXHAUSTIVE")
print("=" * 70)

n = 6
num_edges = comb(n, 2)  # 15 edges, 2^15 = 32768 tournaments
print(f"  n={n}: {2**num_edges} tournaments, {num_edges} edges")
print(f"  Computing I(x) for x ∈ {{2,3,5,7}}...")

residue_data = {x: {p: Counter() for p in [3,5,7,11]} for x in [2,3,5,7]}
count6 = 0

edges = [(i, j) for i in range(n) for j in range(i+1, n)]
for bits in range(2**num_edges):
    adj = [[0]*n for _ in range(n)]
    for idx, (i, j) in enumerate(edges):
        if bits & (1 << idx):
            adj[i][j] = 1
        else:
            adj[j][i] = 1

    H = count_hamiltonian_paths(adj, n)
    cycles = get_directed_odd_cycles(adj, n)
    alpha1 = len(cycles)
    cg = build_conflict_graph(cycles)
    nc = len(cycles)

    for x in [2, 3, 5, 7]:
        Ix = independence_polynomial_at(cg, nc, x)
        for p in [3, 5, 7, 11]:
            if p != x:
                residue_data[x][p][Ix % p] += 1

    count6 += 1
    if count6 % 5000 == 0:
        print(f"    {count6}/{2**num_edges}...")

print(f"  Done ({count6} tournaments).")
print(f"\n  Forbidden residues at n=6:")
print(f"  {'x':>3} {'mod p':>6}  {'forbidden':>20}  {'n=5 had':>15}")

# Compare with n=5
n5_forbidden = {}
for x in [2, 3, 5, 7]:
    for p in [3, 5, 7, 11]:
        if p == x:
            continue
        res5 = Counter()
        for H, a1, I_vals in all_data:
            res5[I_vals[x] % p] += 1
        missing5 = [r for r in range(p) if res5[r] == 0]
        n5_forbidden[(x, p)] = missing5

for x in [2, 3, 5, 7]:
    for p in [3, 5, 7, 11]:
        if p == x:
            continue
        missing6 = [r for r in range(p) if residue_data[x][p][r] == 0]
        missing5 = n5_forbidden.get((x, p), [])
        if missing5 or missing6:
            lifted = set(missing5) - set(missing6)
            notes = ""
            if lifted:
                notes = f"LIFTED: {lifted}"
            print(f"  {x:>3}   mod {p:>2}  forbidden: {str(missing6):>15}  n=5: {str(missing5):>15}  {notes}")

# ====================================================================
# PART 6: DEEPER — WHAT CONTROLS I(x) mod p?
# ====================================================================
print("\n" + "=" * 70)
print("PART 6: STRUCTURE OF I(x) mod p — THE OCF MODULAR THEORY")
print("=" * 70)

print("""
  I(x) mod p depends on (α₁, α₂, ...) mod p and on x mod p.

  Key insight: I(x) mod p = I(x mod p) mod p
  (since x^k ≡ (x mod p)^k mod p)

  So I(7) mod 7 ≡ I(0) mod 7 = 1  ← trivially!
  And I(7) mod 5 ≡ I(2) mod 5 = H mod 5  ← same as H mod 5!
  And I(7) mod 3 ≡ I(1) mod 3 = (1+α₁+α₂+...) mod 3

  This is crucial: I(x) mod p only depends on x mod p.

  CONSEQUENCES:
  I(7) mod 5 = I(2) mod 5 = H mod 5  [since 7≡2 mod 5]
  I(7) mod 3 = I(1) mod 3            [since 7≡1 mod 3]
  I(8) mod 3 = I(2) mod 3 = H mod 3  [since 8≡2 mod 3]
  I(8) mod 5 = I(3) mod 5            [since 8≡3 mod 5]

  So forbidden residues of I(7) mod 5 ARE the forbidden residues of H mod 5!
  And forbidden residues of I(8) mod 5 ARE the forbidden residues of I(3) mod 5.
""")

# Verify
print("  Verification at n=5:")
for x, p, equiv_x in [(7, 5, 2), (7, 3, 1), (8, 3, 2), (8, 5, 3), (5, 3, 2), (5, 7, 5), (11, 3, 2), (11, 7, 4)]:
    residues_x = Counter()
    residues_eq = Counter()
    for H, a1, I_vals in all_data:
        residues_x[I_vals[x] % p] += 1
        residues_eq[I_vals[equiv_x] % p] += 1
    match = (residues_x == residues_eq)
    print(f"  I({x}) mod {p} = I({equiv_x}) mod {p} [{x}≡{equiv_x} mod {p}]: {'✓' if match else '✗'}")

# ====================================================================
# PART 7: I(1) = 1 + α₁ + α₂ + ... = TOTAL INDEPENDENT SETS
# ====================================================================
print("\n" + "=" * 70)
print("PART 7: I(1) = TOTAL NUMBER OF INDEPENDENT SETS")
print("=" * 70)

print("""
  I(1) = Σ_k α_k = total number of independent sets in CG(T)
  (including the empty set)

  At n=5: I(1) = 1 + α₁ + α₂

  This is the "x=1 world" — the Lucas/Fibonacci world!

  Since I(x) mod p for x ≡ 1 (mod p) gives I(1) mod p,
  the values I(7) mod 3, I(8) mod 7, etc. all reduce to I(1) mod p.
""")

i1_data = Counter()
for H, a1, I_vals in all_data:
    # I(1) = 1 + a1 + a2 where a2 = (H-1-2*a1)/4
    a2 = (H - 1 - 2*a1) // 4
    i1 = 1 + a1 + a2
    i1_data[i1] += 1

print(f"  I(1) spectrum at n=5:")
for val in sorted(i1_data.keys()):
    print(f"    I(1) = {val}: {i1_data[val]} tournaments")

print(f"\n  I(1) mod 3:")
mod3 = Counter()
for val, cnt in i1_data.items():
    mod3[val % 3] += cnt
print(f"    {dict(sorted(mod3.items()))}")
missing = [r for r in range(3) if mod3[r] == 0]
print(f"    Forbidden: {missing if missing else 'none'}")

print(f"\n  I(1) mod 5:")
mod5 = Counter()
for val, cnt in i1_data.items():
    mod5[val % 5] += cnt
print(f"    {dict(sorted(mod5.items()))}")

# ====================================================================
# PART 8: THE RECURRENCE WEB — I(x) SATISFIES (5,6)-RECURRENCE?
# ====================================================================
print("\n" + "=" * 70)
print("PART 8: DOES I(x) SATISFY THE (5,6)-RECURRENCE IN x?")
print("=" * 70)

print("""
  Question: Does I(x+1) = 5·I(x) - 6·I(x-1) hold for all x?

  I(x) = 1 + α₁x + α₂x² (at n=5)
  I(x+1) = 1 + α₁(x+1) + α₂(x+1)²
  5I(x) = 5 + 5α₁x + 5α₂x²
  6I(x-1) = 6 + 6α₁(x-1) + 6α₂(x-1)²

  5I(x) - 6I(x-1) = 5 + 5α₁x + 5α₂x² - 6 - 6α₁x + 6α₁ - 6α₂x² + 12α₂x - 6α₂
                   = -1 + 6α₁ - 6α₂ + x(-α₁ + 12α₂) + x²(-α₂)

  I(x+1) = 1 + α₁x + α₁ + α₂x² + 2α₂x + α₂
          = 1 + α₁ + α₂ + x(α₁ + 2α₂) + x²α₂

  These are NOT equal in general!
  5I(x)-6I(x-1) ≠ I(x+1) unless special conditions on α₁, α₂.

  Required: -1 + 6α₁ - 6α₂ = 1 + α₁ + α₂  → 5α₁ = 2 + 7α₂
            -α₁ + 12α₂ = α₁ + 2α₂  → 2α₁ = 10α₂ → α₁ = 5α₂
            -α₂ = α₂ → α₂ = 0

  α₂ = 0 AND α₁ = 0. Only the trivial tournament!

  So I(x) does NOT satisfy the (5,6)-recurrence in x.
  This makes sense: I is a polynomial of degree ≤ max_indep_size,
  not an exponential.
""")

# But what about a DIFFERENT recurrence?
print("  However: I(x) IS a polynomial in x of degree ≤ α_max.")
print("  At n=5: degree ≤ 2, so I(x) is determined by 3 points.")
print("  The 'recurrence' is: I(x) = quadratic interpolation!")
print()
print("  Lagrange form through (0,1), (2,H), (3,I₃):")
print("    I(x) = 1·(x-2)(x-3)/((0-2)(0-3))")
print("         + H·(x-0)(x-3)/((2-0)(2-3))")
print("         + I₃·(x-0)(x-2)/((3-0)(3-2))")
print("         = (x²-5x+6)/6 - H·x(x-3)/2 + I₃·x(x-2)/3")

# Verify Lagrange formula
print("\n  Verification:")
errors = 0
for H, a1, I_vals in all_data[:10]:
    I3 = I_vals[3]
    for x in [5, 7, 8, 11]:
        lagrange = ((x**2-5*x+6)/6 - H*x*(x-3)/2 + I3*x*(x-2)/3)
        if abs(lagrange - I_vals[x]) > 0.01:
            errors += 1
            print(f"    ERROR: x={x}, H={H}, I₃={I3}, Lagrange={lagrange}, actual={I_vals[x]}")

if errors == 0:
    print("  Lagrange interpolation: ✓ ALL CORRECT (checked first 10)")

# ====================================================================
# PART 9: THE CONNECTION — WHY 5 AND 6 APPEAR IN LAGRANGE
# ====================================================================
print("\n" + "=" * 70)
print("PART 9: 5 AND 6 IN THE LAGRANGE FORMULA")
print("=" * 70)

print("""
  The Lagrange formula through (0,1), (2,H), (3,I₃):

    I(x) = (x²-5x+6)/6 - Hx(x-3)/2 + I₃·x(x-2)/3

  The numerator of the first term is x² - 5x + 6 = (x-2)(x-3)!
  And the denominator is 6 = 2·3!

  So: I(x) = (x-2)(x-3)/(2·3) - Hx(x-3)/2 + I₃·x(x-2)/3

  The nodes are at x = 0, 2, 3 — zero and the two keys!
  The polynomial (x-2)(x-3) = x² - 5x + 6 IS the tournament polynomial!

  Rewrite:
    I(x) = (x-2)(x-3)/6 + Hx(3-x)/2 + I₃x(x-2)/3

  At x = 2: (0)/6 + H·2·1/2 + I₃·0 = H ✓
  At x = 3: 0 + 0 + I₃·3·1/3 = I₃ ✓
  At x = 0: (-2)(-3)/6 + 0 + 0 = 1 ✓

  The TOURNAMENT POLYNOMIAL z²-5z+6 appears as the Lagrange basis
  polynomial for the x=0 node, when interpolating through the keys!

  This is the deep structural reason why (5,6) appears:
  The independence polynomial I(x), evaluated at the keys 2 and 3,
  gives H and I₃. The Lagrange interpolation formula automatically
  generates z²-5z+6 as its first basis polynomial.
""")

# ====================================================================
# PART 10: AT n=7+ — DEGREE OF I(x) AND IMPLICATIONS
# ====================================================================
print("=" * 70)
print("PART 10: DEGREE OF I(x) AT HIGHER n")
print("=" * 70)

print("""
  At n=5: I(x) has degree ≤ 2 (max independent set size = 2)
    → 3 points determine I(x) completely
    → Lagrange through (0,1), (2,H), (3,I₃) gives everything

  At n=7: I(x) has degree ≤ 3 (α₃ can appear)
    → 4 points needed: (0,1), (2,H), (3,I₃), (x₄, I(x₄))
    → The natural fourth point is x₄ = 5 (next hierarchy number after 3)

  At n=9: I(x) has degree ≤ 4
    → 5 points: (0,1), (2,H), (3,I₃), (5,I₅), (7,I₇)?

  PATTERN: The "natural" interpolation nodes are:
    0 (always), 2, 3, 5, 7, 11, ...
    These are 0 followed by the hierarchy numbers!

  Each new level of the OCF (a new α_k becoming nonzero)
  requires one more evaluation to determine I(x) fully.

  The natural evaluation points are:
    x = 0: gives constant term (always 1)
    x = 2: gives H (Hamiltonian paths)
    x = 3: gives I₃ (the "3-world")
    x = 5: gives I₅ (the "5-world")
    x = 7: gives I₇ (the "7-world")

  THESE ARE THE PRIMES (and 0)!

  More precisely: 0 and the primes 2,3,5,7,... are the NATURAL
  evaluation points for the independence polynomial, because:
  1. They give integer values (I(x) is a polynomial with integer coefficients)
  2. They are coprime to each other (giving maximal information)
  3. The Lagrange basis through them involves the simplest denominators
""")

# What's the max independent set size at n=7?
print("  Checking max independent set sizes by sampling at n=7...")
random.seed(42)
max_indep_sizes = []
for _ in range(200):
    adj = random_tournament(7)
    cycles = get_directed_odd_cycles(adj, 7)
    nc = len(cycles)
    if nc == 0:
        max_indep_sizes.append(0)
        continue
    cg = build_conflict_graph(cycles)
    # Find max independent set size
    max_size = 0
    for mask in range(1 << min(nc, 20)):  # limit for large CGs
        verts = [i for i in range(min(nc, 20)) if mask & (1 << i)]
        is_indep = True
        for i in range(len(verts)):
            for j in range(i+1, len(verts)):
                if cg[verts[i]][verts[j]]:
                    is_indep = False
                    break
            if not is_indep:
                break
        if is_indep:
            max_size = max(max_size, len(verts))
    max_indep_sizes.append(max_size)

size_counts = Counter(max_indep_sizes)
print(f"  Max independent set size distribution at n=7 (200 samples):")
for s in sorted(size_counts.keys()):
    print(f"    size {s}: {size_counts[s]} tournaments ({100*size_counts[s]/200:.1f}%)")

print(f"\n  Max degree of I(x) at n=7: {max(max_indep_sizes)}")

# ====================================================================
# PART 11: SYNTHESIS
# ====================================================================
print("\n" + "=" * 70)
print("PART 11: SYNTHESIS — THE POLYNOMIAL-CYCLOTOMIC-TOURNAMENT TRINITY")
print("=" * 70)

print("""
  THREE LAYERS OF STRUCTURE:

  1. THE INDEPENDENCE POLYNOMIAL I(x):
     I(x) = 1 + α₁x + α₂x² + ... (finite polynomial)
     I(x) ≡ 1 (mod x) always → Rédei's theorem is x=2 case
     Degree = max independent set size in CG(T)

  2. THE TOURNAMENT POLYNOMIAL z²-5z+6:
     Roots are the keys 2,3
     Appears as Lagrange basis through (0,2,3)
     Its coefficients (5,6) encode the universe
     The recurrence x(n)=5x(n-1)-6x(n-2) generates 2^n, 3^n

  3. THE CYCLOTOMIC DICTIONARY Φ_d(2):
     Φ₁=1, Φ₂=3, Φ₃=7, Φ₄=5, Φ₅=31, Φ₁₀=11
     The hierarchy IS the cyclotomic evaluation at x=2
     2^n - 1 = ∏_{d|n} Φ_d(2) gives Mersenne factorization
     ord_p(2) = d ↔ period of OCF coefficients mod p = d

  THE BRIDGE between layers:
     I(Φ_d(2)) ≡ I(2^d mod Φ_d(2)) ≡ I(0) ≡ 1 (mod Φ_d(2))
     [since 2^d ≡ 0 mod (2^d-1) factors involving Φ_d]

  This is trivially true but CONNECTS the three worlds:
     The cyclotomic primes Φ_d(2) are exactly the primes p
     where 2 has SHORT multiplicative order (d | p-1).
     Short period → MORE constraints on H mod p →
     forbidden residues PERSIST to larger n.

  HIERARCHY OF CONSTRAINTS:
     d=1: Φ₁(2)=1 — no constraint
     d=2: Φ₂(2)=3 — period 2: H≡1 mod 2 (Rédei)
     d=3: Φ₃(2)=7 — period 3: H≡0 mod 7 forbidden at n≤6
     d=4: Φ₄(2)=5 — period 4: H≡2 mod 5 forbidden at n≤5
     d=5: Φ₅(2)=31 — period 5: ???
     d=10: Φ₁₀(2)=11 — period 10: weakest constraint
""")

print("\nDone.")
