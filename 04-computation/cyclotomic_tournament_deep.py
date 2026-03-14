#!/usr/bin/env python3
"""
cyclotomic_tournament_deep.py — opus-2026-03-14-S74
Deep exploration of the cyclotomic dictionary at x=2 and its tournament meaning.

Key discovery from five_six_to_seven_eight.py:
  7 = Φ₃(2), 3 = Φ₂(2), 5 = Φ₄(2), 11 = Φ₁₀(2), 31 = Φ₅(2)

The ENTIRE hierarchy of tournament numbers appears as cyclotomic evaluations at x=2.
This script explores WHY, and what cyclotomic structure means for tournaments.

The factorization 2^n - 1 = ∏_{d|n} Φ_d(2) gives Mersenne numbers as products
of hierarchy elements. This must connect to the OCF structure.
"""

import numpy as np
from itertools import combinations, permutations
from math import gcd, factorial, comb
from functools import reduce

# ====================================================================
# PART 1: COMPLETE CYCLOTOMIC DICTIONARY AT x=2
# ====================================================================
print("=" * 70)
print("PART 1: CYCLOTOMIC DICTIONARY — Φ_d(2) FOR d = 1..30")
print("=" * 70)

def cyclotomic_poly_at_x(d, x):
    """Evaluate Φ_d(x) using the product formula:
    Φ_d(x) = ∏_{k=1..d, gcd(k,d)=1} (x - e^{2πik/d})

    For integer x, we can use: Φ_d(x) = ∏_{k|d} (x^k - 1)^{μ(d/k)}
    where μ is Möbius function.
    """
    def mobius(n):
        """Möbius function."""
        if n == 1:
            return 1
        # Factor n
        factors = []
        temp = n
        d2 = 2
        while d2 * d2 <= temp:
            if temp % d2 == 0:
                factors.append(d2)
                temp //= d2
                if temp % d2 == 0:
                    return 0  # squared factor
            d2 += 1
        if temp > 1:
            factors.append(temp)
        return (-1) ** len(factors)

    # Φ_d(x) = ∏_{k|d} (x^k - 1)^{μ(d/k)}
    result_num = 1
    result_den = 1
    for k in range(1, d + 1):
        if d % k == 0:
            mu = mobius(d // k)
            val = x**k - 1
            if mu == 1:
                result_num *= val
            elif mu == -1:
                result_den *= val
    return result_num // result_den

print(f"\n  {'d':>3}  {'Φ_d(2)':>12}  {'φ(d)':>4}  {'prime?':>6}  {'notes'}")
print(f"  {'---':>3}  {'------':>12}  {'----':>4}  {'------':>6}  {'-----'}")

hierarchy_numbers = set()
for d in range(1, 31):
    val = cyclotomic_poly_at_x(d, 2)
    euler_phi = sum(1 for k in range(1, d + 1) if gcd(k, d) == 1)
    is_prime = val > 1 and all(val % p != 0 for p in range(2, int(val**0.5) + 1))

    notes = []
    if val in [1, 2, 3, 5, 7, 11, 15, 31]:
        notes.append("← HIERARCHY")
        hierarchy_numbers.add(val)
    if d <= 12 and val == 2**d - 1 and d > 0:
        # Check if Mersenne prime
        pass

    # Check for tournament significance
    if val == 7:
        notes.append("n=7 barrier")
    elif val == 15:
        notes.append("max H at n=5")
    elif val == 31:
        notes.append("Mersenne prime, 2⁵-1")
    elif val == 11:
        notes.append("L(5)=J(5)")
    elif val == 5:
        notes.append("F(5)=5")
    elif val == 3:
        notes.append("Rédei parity")

    p_str = "yes" if is_prime else ""
    n_str = ", ".join(notes)
    print(f"  {d:>3}  {val:>12}  {euler_phi:>4}  {p_str:>6}  {n_str}")

# ====================================================================
# PART 2: MERSENNE FACTORIZATION VIA CYCLOTOMICS
# ====================================================================
print("\n" + "=" * 70)
print("PART 2: MERSENNE FACTORIZATION 2^n - 1 = ∏_{d|n} Φ_d(2)")
print("=" * 70)

for n in range(1, 16):
    mersenne = 2**n - 1
    divisors = [d for d in range(1, n + 1) if n % d == 0]
    factors = [cyclotomic_poly_at_x(d, 2) for d in divisors]
    product = reduce(lambda a, b: a * b, factors)
    assert product == mersenne, f"Product mismatch at n={n}"

    factor_str = " · ".join(f"Φ_{d}={f}" for d, f in zip(divisors, factors))
    print(f"  n={n:>2}: 2^{n}-1 = {mersenne:>8} = {factor_str}")

# ====================================================================
# PART 3: THE (5,6)-RECURRENCE IN TOURNAMENT ARITHMETIC
# ====================================================================
print("\n" + "=" * 70)
print("PART 3: THE (5,6)-RECURRENCE x(n) = 5x(n-1) - 6x(n-2)")
print("=" * 70)

print("\n  This recurrence has characteristic polynomial z²-5z+6 = (z-2)(z-3).")
print("  General solution: x(n) = A·2^n + B·3^n")
print("\n  Different initial conditions give different sequences:")

sequences = {
    "2^n":     (1, 2),     # A=1, B=0
    "3^n":     (1, 3),     # A=0, B=1
    "3^n-2^n": (0, 1),     # A=-1, B=1 (Jacobsthal-like)
    "3^n+2^n": (2, 7),     # A=1, B=1 (sum)
    "(3^n-2^n)/1": (0, 1), # Same as 3^n-2^n
}

# More interesting: what tournament-meaningful sequences satisfy this?
print(f"\n  {'name':>15}  {'x(0)':>5}  {'x(1)':>5}  x(2..8)")
print(f"  {'-'*15}  {'-'*5}  {'-'*5}  {'-'*30}")

def recurrence_56(x0, x1, n_terms=9):
    seq = [x0, x1]
    for _ in range(n_terms - 2):
        seq.append(5 * seq[-1] - 6 * seq[-2])
    return seq

interesting = [
    ("2^n", 1, 2),
    ("3^n", 1, 3),
    ("3^n - 2^n", 0, 1),
    ("3^n + 2^n", 2, 7),
    ("2·3^n - 2^n", 1, 4),
    ("H=1 (init)", 1, 5),
    ("H=3 (init)", 3, 9),  # 3·2^n + 0·3^n? No: A·1+B·1=3, A·2+B·3=9 → A=0,B=3 → 3^{n+1}
    ("OCF: 1,1+2α₁", 1, None),
]

for name, x0, x1 in interesting:
    if x1 is None:
        continue
    seq = recurrence_56(x0, x1)
    seq_str = ", ".join(f"{s}" for s in seq[2:])
    print(f"  {name:>15}  {x0:>5}  {x1:>5}  {seq_str}")

# ====================================================================
# PART 4: EVALUATIONS I(x) AT x = 5, 6, 7, 8 — EXACT FORMULAS
# ====================================================================
print("\n" + "=" * 70)
print("PART 4: I(x) AT x = 5, 6, 7, 8 — FORMULAS FROM OCF")
print("=" * 70)

print("""
  The OCF gives I(Ω, x) = Σ_S x^|S| where S ranges over independent sets.
  With x=2: H = I(2) = 1 + 2α₁ + 4α₂ + 8α₃ + ...

  For general x:
    I(x) = 1 + x·α₁ + x²·α₂ + x³·α₃ + ...

  where α_k = number of independent sets of size k in CG(T).

  So the FULL information is encoded in the sequence (α₁, α₂, α₃, ...).

  At n=5 (where α₂ is the max level):
    I(x) = 1 + x·α₁ + x²·α₂

  Relations:
    I(2) = H = 1 + 2α₁ + 4α₂
    I(3) = 1 + 3α₁ + 9α₂
    I(5) = 1 + 5α₁ + 25α₂
    I(6) = 1 + 6α₁ + 36α₂
    I(7) = 1 + 7α₁ + 49α₂
    I(8) = 1 + 8α₁ + 64α₂
""")

# Express I(5), I(6), I(7), I(8) in terms of H and α₁
print("  FORMULAS (using α₂ = (H - 1 - 2α₁)/4):")
print()
for x_val in [3, 5, 6, 7, 8, 11]:
    # I(x) = 1 + x·α₁ + x²·(H-1-2α₁)/4
    # = 1 + x·α₁ + x²H/4 - x²/4 - x²α₁/2
    # = 1 - x²/4 + α₁(x - x²/2) + x²H/4
    # = (4 - x² + 4xα₁ - 2x²α₁ + x²H) / 4
    # = (x²H + (4x - 2x²)α₁ + 4 - x²) / 4
    a = x_val**2
    b = 4*x_val - 2*x_val**2
    c = 4 - x_val**2
    print(f"    I({x_val}) = ({a}H + ({b})α₁ + ({c})) / 4")
    print(f"         = ({a}H {b:+d}·α₁ {c:+d}) / 4")
    # Simplify
    g = gcd(gcd(abs(a), abs(b)), abs(c))
    if g > 1 and 4 % g == 0:
        print(f"         = ({a//g}H {b//g:+d}·α₁ {c//g:+d}) / {4//g}")
    print()

# ====================================================================
# PART 5: VERIFY I(x) FORMULAS AT n=5 (all 1024 tournaments)
# ====================================================================
print("=" * 70)
print("PART 5: VERIFY I(x) FORMULAS AT n=5")
print("=" * 70)

def enumerate_tournaments(n):
    """Generate all tournaments on n vertices."""
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    for bits in range(2**len(edges)):
        adj = [[0]*n for _ in range(n)]
        for idx, (i, j) in enumerate(edges):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        yield adj

def count_hamiltonian_paths(adj, n):
    """Count Hamiltonian paths using DP on bitmask."""
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
    """Get all directed odd cycles (as frozen sets of directed edges)."""
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
                    # Canonical form: start from min vertex, going in the direction that gives smaller second vertex
                    min_idx = perm.index(min(perm))
                    canon = tuple(perm[min_idx:] + perm[:min_idx])
                    cycles.append(canon)
    # Remove duplicates
    return list(set(cycles))

def build_conflict_graph(cycles):
    """Build conflict graph: edge between cycles that share a vertex."""
    nc = len(cycles)
    vsets = [set(c) for c in cycles]
    adj_cg = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if vsets[i] & vsets[j]:
                adj_cg[i][j] = True
                adj_cg[j][i] = True
    return adj_cg

def independence_polynomial(adj_cg, nc, x):
    """Compute I(CG, x) = sum over independent sets S of x^|S|."""
    total = 0
    for mask in range(1 << nc):
        # Check independence
        is_indep = True
        verts = [i for i in range(nc) if mask & (1 << i)]
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

n = 5
print(f"\n  Enumerating all {2**comb(n,2)} tournaments at n={n}...")

errors = {x: 0 for x in [2, 3, 5, 6, 7, 8, 11]}
count = 0
data = []

for adj in enumerate_tournaments(n):
    H = count_hamiltonian_paths(adj, n)
    cycles = get_directed_odd_cycles(adj, n)
    alpha1 = len(cycles)

    # Build CG and compute I(x) directly
    cg = build_conflict_graph(cycles)
    nc = len(cycles)

    I_vals = {}
    for x in [2, 3, 5, 6, 7, 8, 11]:
        I_vals[x] = independence_polynomial(cg, nc, x)

    # Check H = I(2)
    if I_vals[2] != H:
        errors[2] += 1

    # Check formulas: I(x) = (x²H + (4x-2x²)α₁ + 4-x²) / 4
    for x in [3, 5, 6, 7, 8, 11]:
        formula = (x**2 * H + (4*x - 2*x**2) * alpha1 + 4 - x**2)
        if formula % 4 != 0:
            errors[x] += 1
        elif formula // 4 != I_vals[x]:
            errors[x] += 1

    data.append((H, alpha1, I_vals))
    count += 1

print(f"  Checked {count} tournaments.")
print(f"\n  Verification results:")
for x in [2, 3, 5, 6, 7, 8, 11]:
    status = "✓ ALL CORRECT" if errors[x] == 0 else f"✗ {errors[x]} errors"
    print(f"    I({x:>2}) formula: {status}")

# ====================================================================
# PART 6: THE LATTICE OF EVALUATIONS — I(x) mod small primes
# ====================================================================
print("\n" + "=" * 70)
print("PART 6: I(x) MOD SMALL PRIMES AT n=5")
print("=" * 70)

from collections import Counter

for x in [2, 5, 6, 7, 8]:
    for p in [3, 5, 7, 11]:
        residues = Counter()
        for H, a1, I_vals in data:
            residues[I_vals[x] % p] += 1
        missing = [r for r in range(p) if residues[r] == 0]
        if missing:
            print(f"  I({x}) mod {p}: FORBIDDEN residues {missing}")
            dist = ", ".join(f"{r}:{residues[r]}" for r in range(p) if residues[r] > 0)
            print(f"    Distribution: {dist}")

# ====================================================================
# PART 7: 7 = Φ₃(2) — CUBE ROOTS OF UNITY MOD 7 IN TOURNAMENTS
# ====================================================================
print("\n" + "=" * 70)
print("PART 7: 7 = Φ₃(2) — MEANING OF CYCLOTOMIC STRUCTURE")
print("=" * 70)

print("""
  WHY 7 = Φ₃(2):

  Φ₃(x) = x² + x + 1. At x=2: 4+2+1 = 7.

  The roots of Φ₃(x) are the primitive 3rd roots of unity: ω, ω².
  So 7 = (2-ω)(2-ω²) where ω = e^{2πi/3}.

  In F₇: since Φ₃(2) = 0 mod 7, we have 2²+2+1 ≡ 0 (mod 7).
  This means 2 is a primitive cube root of unity in F₇!
  Indeed: 2³ = 8 ≡ 1 (mod 7), and ord₇(2) = 3.

  For tournaments: H ≡ 1 (mod 2) always (Rédei).
  Since ord₇(2) = 3, the powers of 2 cycle through {1, 2, 4} mod 7.

  In the OCF: H = 1 + 2α₁ + 4α₂ + 8α₃ + ...
  Mod 7:      H ≡ 1 + 2α₁ + 4α₂ + 1·α₃ + 2α₄ + 4α₅ + ...  (mod 7)

  The coefficients cycle with period 3: {2, 4, 1, 2, 4, 1, ...}
  This is because 2³ ≡ 1 (mod 7), which IS the cyclotomic condition!
""")

# Verify the period-3 cycling
print("  Powers of 2 mod 7:")
for k in range(10):
    print(f"    2^{k} = {2**k} ≡ {2**k % 7} (mod 7)")

print(f"\n  Period = 3 = order of the cyclotomic polynomial Φ₃.")
print(f"  The period-3 cycling means: H mod 7 depends only on")
print(f"  (α₁ mod 7, α₂ mod 7, (α₁+α₂+α₃) mod 7, ...)")

# Now do the same for 5 = Φ₄(2) and 3 = Φ₂(2)
print(f"\n  COMPARISON WITH OTHER HIERARCHY PRIMES:")
print(f"  3 = Φ₂(2): ord₃(2) = 2, period-2 cycling: {{2,1,2,1,...}}")
print(f"  5 = Φ₄(2): ord₅(2) = 4, period-4 cycling: {{2,4,3,1,...}}")
print(f"  7 = Φ₃(2): ord₇(2) = 3, period-3 cycling: {{2,4,1,...}}")
print(f"  11= Φ₁₀(2): ord₁₁(2)= 10, period-10 cycling")
print(f"  31= Φ₅(2): ord₃₁(2) = 5, period-5 cycling")

print(f"\n  DEEP PATTERN: Φ_d(2) = p ↔ ord_p(2) = d")
print(f"  The cyclotomic index d IS the multiplicative order of 2 mod p!")
print(f"  This is NOT a coincidence — it's the definition of cyclotomic polynomials.")
print(f"\n  TOURNAMENT MEANING:")
print(f"  H mod p has OCF coefficients cycling with period d = ord_p(2).")
print(f"  Shorter period = MORE constraints on H mod p.")
print(f"  ord₃(2)=2: very constrained → only H≡1(mod2) survives")
print(f"  ord₇(2)=3: 3-periodic → still quite constrained at small n")
print(f"  ord₅(2)=4: 4-periodic → less constrained, but still interesting")

# ====================================================================
# PART 8: THE (5,6) RECURRENCE IN I-POLYNOMIAL SPACE
# ====================================================================
print("\n" + "=" * 70)
print("PART 8: (5,6)-RECURRENCE CONNECTING I(2) AND I(3)")
print("=" * 70)

print("""
  Since I(x) = 1 + x·α₁ + x²·α₂ (at n=5), we have:

    I(2) = H = 1 + 2α₁ + 4α₂
    I(3) = 1 + 3α₁ + 9α₂

  So: I(3) - I(2) = α₁ + 5α₂
  And: 3·I(2) - 2·I(3) = 1 - 6α₂

  More generally, for the (5,6)-recurrence:
    x(n) = 5x(n-1) - 6x(n-2)  has roots 2, 3.

  Apply this to I(x) viewed as evaluating at x=2,3,4,...:
    I(x+1) = 5·I(x) - 6·I(x-1)?

  NO — I is a polynomial in x, not a linear recurrence in x.
  But we CAN write:
    I(x) = 1 + α₁·x + α₂·x² = quadratic in x

  Three evaluations determine it completely:
    I(0) = 1 (always)
    I(2) = H
    I(3) = (9H - 5 - 6α₁)/4  (from the formula)

  So I(x) is determined by (H, α₁), which is determined by
  any two evaluations at x ≠ 0.
""")

# The key insight: I(2) and I(3) together determine EVERYTHING at n=5
print("  I(2)=H and I(3) determine all I(x) at n=5:")
print(f"\n  {'H':>5} {'α₁':>5} {'I(2)':>5} {'I(3)':>5} {'I(5)':>6} {'I(6)':>6} {'I(7)':>7} {'I(8)':>7}")

seen = set()
for H, a1, I_vals in data:
    key = (H, a1)
    if key not in seen:
        seen.add(key)
        print(f"  {H:>5} {a1:>5} {I_vals[2]:>5} {I_vals[3]:>5} {I_vals[5]:>6} {I_vals[6]:>6} {I_vals[7]:>7} {I_vals[8]:>7}")

# ====================================================================
# PART 9: THE 6-SPECTRUM AND ITS TOURNAMENT MEANING
# ====================================================================
print("\n" + "=" * 70)
print("PART 9: I(6) = I(2·3) — THE PRODUCT-KEY EVALUATION")
print("=" * 70)

print("""
  I(6) is special because 6 = 2·3 (product of keys).

  At n=5: I(6) = 1 + 6α₁ + 36α₂
         = (36H + (-48)α₁ + (-32)) / 4
         = 9H - 12α₁ - 8

  When α₂=0: I(6) = 1 + 6α₁ = 1 + 6·(H-1)/2 = 3H - 2 ✓

  The I(6) spectrum at n=5:
""")

i6_counts = Counter()
for H, a1, I_vals in data:
    i6_counts[I_vals[6]] += 1

print(f"  {'I(6)':>6}  {'count':>5}")
for val in sorted(i6_counts.keys()):
    print(f"  {val:>6}  {i6_counts[val]:>5}")

print(f"\n  I(6) values: {sorted(i6_counts.keys())}")
print(f"  All ≡ 1 (mod 6)? {all(v % 6 == 1 for v in i6_counts.keys())}")
print(f"  All ≡ 1 (mod 2)? {all(v % 2 == 1 for v in i6_counts.keys())}")
print(f"  All ≡ 1 (mod 3)? {all(v % 3 == 1 for v in i6_counts.keys())}")

# ====================================================================
# PART 10: INTERPOLATION — RECOVERING (α₁, α₂) FROM TWO EVALUATIONS
# ====================================================================
print("\n" + "=" * 70)
print("PART 10: VANDERMONDE RECOVERY — (α₁,α₂) FROM (I(a), I(b))")
print("=" * 70)

print("""
  I(a) = 1 + a·α₁ + a²·α₂
  I(b) = 1 + b·α₁ + b²·α₂

  Solve:  [a  a²] [α₁]   [I(a)-1]
          [b  b²] [α₂] = [I(b)-1]

  Determinant = ab(b-a) = ab·V(a,b) where V is Vandermonde.

  Minimal |det| over positive integers a<b:
    (a,b)=(1,2): det = 1·2·1 = 2
    (a,b)=(1,3): det = 1·3·2 = 6
    (a,b)=(2,3): det = 2·3·1 = 6

  The KEYS (2,3) give det = 6 = 2·3 = product of keys!
  This is V(2,3) = 6, the Vandermonde determinant.

  Inversion formulas from (H, I(3)):
    α₁ = (3(H-1) - (I(3)-1)) / 6
    α₂ = ((I(3)-1) - 2(H-1)) / 6

  Or equivalently:
    α₁ = (3H - I(3) - 2) / 6
    α₂ = (I(3) - 2H + 1) / 6
""")

# Verify
print("  Verification:")
all_ok = True
for H, a1, I_vals in data:
    I3 = I_vals[3]
    a1_recovered = (3*H - I3 - 2) / 6
    a2_recovered = (I3 - 2*H + 1) / 6
    if a1_recovered != a1:
        all_ok = False
        break
    # Also check α₂
    a2_direct = (H - 1 - 2*a1) // 4
    if a2_recovered != a2_direct:
        all_ok = False
        break

print(f"  α₁ recovery from (H, I(3)): {'✓ ALL CORRECT' if all_ok else '✗ ERRORS'}")

# ====================================================================
# PART 11: THE FIVE DECOMPOSITIONS OF 7 AND 8
# ====================================================================
print("\n" + "=" * 70)
print("PART 11: DECOMPOSITIONS OF 7 AND 8 IN THE (2,3) UNIVERSE")
print("=" * 70)

print("""
  SEVEN = 2² + 3 = 2·2 + 3     (quadratic + linear)
        = Φ₃(2)                  (cyclotomic)
        = 2³ - 1                  (Mersenne)
        = L(4) = Lucas(4)         (in the x=1 world)
        = C(7,1) = 7             (trivial but: 7 vertices!)

  EIGHT = 2³                     (pure exponential)
        = 2·3 + 2 = 6 + 2       (product of keys + first key)
        = 3² - 1                  (near-square of second key)
        = F(6) = Fibonacci(6)     (in the Fib world)
        = C(8,3)/7 = 56/7 = 8   (3-cycles per vertex at n=8)

  KEY STRUCTURAL DIFFERENCE:
  7 = Φ₃(2) is IRREDUCIBLE over Z — it's a prime, a cyclotomic unit.
  8 = 2³ is COMPOSITE (as a prime factorization) — pure power of first key.

  In tournament theory:
  7: α₁+2α₂ ≡ 3 (mod 7) impossible at n≤6 → H never ≡ 0 (mod 7)
  8: H values at n=5 are {1,3,5,9,11,13,15} — all odd, none divisible by 8
     At n=5, max H = 15, so 8 divides no H value because H ≤ 15 and H odd.
     At n=6, H can reach ~80+, so 8 doesn't divide odd numbers → NEVER.
     Wait: H is always odd! So H mod 8 ∈ {1,3,5,7} always.
""")

print("  H mod 8 at n=5:")
mod8 = Counter()
for H, a1, I_vals in data:
    mod8[H % 8] += 1
print(f"  {dict(sorted(mod8.items()))}")

print("\n  H mod 7 at n=5:")
mod7 = Counter()
for H, a1, I_vals in data:
    mod7[H % 7] += 1
print(f"  {dict(sorted(mod7.items()))}")
forbidden_7 = [r for r in range(7) if mod7[r] == 0]
print(f"  Forbidden: {forbidden_7}")

# ====================================================================
# PART 12: THE RATIO I(3)/I(2) = I(3)/H AND ITS CONVERGENCE
# ====================================================================
print("\n" + "=" * 70)
print("PART 12: THE RATIO I(3)/H AND ITS ASYMPTOTICS")
print("=" * 70)

print("""
  At n=5: I(3)/H = (9H - 5 - 6α₁) / (4H)

  If α₁ ≈ c·H for large H:
    I(3)/H ≈ (9 - 6c)/4

  For c = (H-1)/(2H) ≈ 1/2:
    I(3)/H ≈ (9-3)/4 = 3/2

  This is the "3/2 ratio" observed in cycle_multiplicity.out!

  More precisely: I(3)/I(2) → (3/2)^? as n → ∞

  If α₂ << α₁ << H, then:
    I(3)/I(2) ≈ (1 + 3α₁)/(1 + 2α₁) → 3/2 as α₁ → ∞

  The 3/2 ratio IS the ratio of the keys!
""")

ratios = []
for H, a1, I_vals in data:
    if H > 0:
        ratios.append(I_vals[3] / H)

print(f"  At n=5:")
print(f"    I(3)/H range: [{min(ratios):.4f}, {max(ratios):.4f}]")
print(f"    I(3)/H mean:  {np.mean(ratios):.4f}")
print(f"    3/2 = 1.5000")
print(f"    Convergence to 3/2 requires α₂/H → 0.")

# Also check I(5)/H → 5/2?
r5 = [I_vals[5]/H for H, a1, I_vals in data if H > 0]
r6 = [I_vals[6]/H for H, a1, I_vals in data if H > 0]
r7 = [I_vals[7]/H for H, a1, I_vals in data if H > 0]
r8 = [I_vals[8]/H for H, a1, I_vals in data if H > 0]

print(f"\n  Ratio I(x)/H means at n=5:")
print(f"    I(3)/H = {np.mean(ratios):.4f}  (→ 3/2 = 1.500)")
print(f"    I(5)/H = {np.mean(r5):.4f}  (→ 5/2 = 2.500)")
print(f"    I(6)/H = {np.mean(r6):.4f}  (→ 6/2 = 3.000)")
print(f"    I(7)/H = {np.mean(r7):.4f}  (→ 7/2 = 3.500)")
print(f"    I(8)/H = {np.mean(r8):.4f}  (→ 8/2 = 4.000)")
print(f"\n  PATTERN: I(x)/H → x/2 as α₂/H → 0 (sparse regime)")
print(f"  This is because I(x) ≈ 1 + x·α₁ ≈ x·(H-1)/2 + 1 ≈ xH/2")

# ====================================================================
# PART 13: THE 5-6 BRIDGE — WHY (5,6) IS THE NATURAL COEFFICIENT PAIR
# ====================================================================
print("\n" + "=" * 70)
print("PART 13: WHY (5,6) = (2+3, 2·3) IS THE NATURAL PAIR")
print("=" * 70)

print("""
  The characteristic polynomial of the (2,3)-universe is:
    z² - (2+3)z + 2·3 = z² - 5z + 6

  By Vieta's formulas:
    sum of roots = 5 = 2+3
    product of roots = 6 = 2·3

  The coefficients (5,6) encode the universe completely:
    Given 5 and 6, you can recover 2 and 3 via the quadratic formula.
    z = (5 ± √(25-24))/2 = (5 ± 1)/2 ∈ {2, 3}

  The discriminant Δ = 5² - 4·6 = 25 - 24 = 1.
  A PERFECT SQUARE! This is why 2 and 3 are both integers.

  Δ = 1 is the smallest possible positive discriminant.
  So (2,3) is the UNIQUE pair of consecutive integers whose
  sum and product give discriminant 1.

  More: (n, n+1) has discriminant (2n+1)² - 4n(n+1) = 4n²+4n+1-4n²-4n = 1.
  So ALL consecutive integer pairs have discriminant 1!
  But (2,3) is the FIRST non-trivial pair (after (0,1) and (1,2)).

  And (1,2) gives z²-3z+2, which has sum=3, product=2.
  So the (1,2)-universe has coefficients (3,2) — the KEYS REVERSED!
""")

# The chain of universes
print("  CHAIN OF CONSECUTIVE-INTEGER UNIVERSES:")
print(f"  {'roots':>10}  {'sum':>4}  {'prod':>4}  {'Δ':>3}  {'char poly':>20}  notes")
for a in range(0, 8):
    b = a + 1
    s, p = a + b, a * b
    disc = s*s - 4*p
    notes = []
    if (a, b) == (2, 3):
        notes.append("← THE KEYS")
    if (a, b) == (1, 2):
        notes.append("← keys reversed")
    if (a, b) == (0, 1):
        notes.append("← trivial")
    if s in [7, 11, 15]:
        notes.append(f"sum={s} hierarchy")
    print(f"  ({a},{b}):  {s:>4}  {p:>4}  {disc:>3}  z²-{s}z+{p}  {', '.join(notes)}")

# ====================================================================
# PART 14: THE (5,6) → (7,8) SHIFT AS "ADDING REDEI"
# ====================================================================
print("\n" + "=" * 70)
print("PART 14: THE +2 SHIFT — ADDING THE FIRST KEY")
print("=" * 70)

print("""
  (5,6) + 2 → (7,8) [adding first key to both]

  What does "+2" mean in tournament theory?

  H = 1 + 2α₁ + 4α₂ + ...
  If we shift x → x+2:
    I(x+2) = 1 + (x+2)α₁ + (x+2)²α₂ + ...
           = 1 + xα₁ + 2α₁ + x²α₂ + 4xα₂ + 4α₂ + ...
           = I(x) + 2α₁ + 4xα₂ + 4α₂ + higher

  At n=5:
    I(x+2) = I(x) + 2α₁ + (4x+4)α₂

  In particular:
    I(7) = I(5) + 2α₁ + 24α₂
    I(8) = I(6) + 2α₁ + 28α₂

  The "+2 shift" adds 2α₁ plus correction terms involving α₂.
  When α₂ = 0: I(x+2) = I(x) + 2α₁ = I(x) + (H-1) exactly!

  So in the sparse regime (α₂=0):
    I(7) = I(5) + H - 1
    I(8) = I(6) + H - 1

  The (5,6)→(7,8) shift IS "adding one copy of (H-1)"!
""")

# Verify
print("  Verification at n=5:")
ok_7 = ok_8 = 0
total_sparse = 0
for H, a1, I_vals in data:
    a2 = (H - 1 - 2*a1) // 4
    if a2 == 0:
        total_sparse += 1
        if I_vals[7] == I_vals[5] + H - 1:
            ok_7 += 1
        if I_vals[8] == I_vals[6] + H - 1:
            ok_8 += 1

print(f"  Sparse (α₂=0): {total_sparse} tournaments")
print(f"  I(7) = I(5) + H-1: {ok_7}/{total_sparse} correct")
print(f"  I(8) = I(6) + H-1: {ok_8}/{total_sparse} correct")

# General case
print(f"\n  General case (with α₂):")
ok_7g = ok_8g = 0
for H, a1, I_vals in data:
    a2 = (H - 1 - 2*a1) // 4
    if I_vals[7] == I_vals[5] + 2*a1 + 24*a2:
        ok_7g += 1
    if I_vals[8] == I_vals[6] + 2*a1 + 28*a2:
        ok_8g += 1
print(f"  I(7) = I(5) + 2α₁ + 24α₂: {ok_7g}/{count} correct")
print(f"  I(8) = I(6) + 2α₁ + 28α₂: {ok_8g}/{count} correct")

# ====================================================================
# PART 15: SYNTHESIS — THE GRAND RECURRENCE WEB
# ====================================================================
print("\n" + "=" * 70)
print("PART 15: SYNTHESIS — THE GRAND RECURRENCE WEB")
print("=" * 70)

print("""
  THE CYCLOTOMIC-TOURNAMENT DICTIONARY:

  Level 0 (trivial):
    1 = Φ₁(2) — the empty tournament, I(x) = 1

  Level 1 (Rédei):
    3 = Φ₂(2) — parity (H odd), period-2 in OCF mod 3

  Level 2 (emergence):
    7 = Φ₃(2) — first barrier broken at n=7, period-3 cycling
    5 = Φ₄(2) — bridge number (F(5)=5), period-4 cycling

  Level 3 (Mersenne):
    31 = Φ₅(2) = 2⁵-1 — Mersenne prime, period-5 cycling
    11 = Φ₁₀(2) — L(5)=J(5), the Lucas-Jacobsthal bridge

  THE (5,6) ENCODING:
    z² - 5z + 6 = 0 encodes (2,3) as (sum, product).
    Discriminant = 1 (perfect square, consecutive integers).
    The recurrence x(n) = 5x(n-1) - 6x(n-2) generates ALL
    linear combinations A·2ⁿ + B·3ⁿ.

  THE (5,6) → (7,8) BRIDGE:
    Adding the first key (2) shifts (sum,product) to (7,8).
    In tournament space: I(7) = I(5) + 2α₁ + 24α₂.
    In sparse regime: I(x+2) = I(x) + H - 1.
    7 = Φ₃(2) is cyclotomic (irreducible, prime).
    8 = 2³ is exponential (composite, pure power).

  THE CONVERGENCE:
    k-nacci → 2 at rate (1/2)^k
    weighted k-nacci → 3 at rate (2/3)^k
    I(x)/H → x/2 in sparse regime
    I(3)/H → 3/2 = ratio of keys

  THE FIVE PRIMES WHERE 2 IS PRIMITIVE ROOT:
    3 (Φ₂), 5 (Φ₄), 11 (Φ₁₀), 13 (Φ₁₂), 19 (Φ₁₈)
    Period of 2 mod p = deg(Φ_d) = φ(d) = p-1.
    These are the primes where H mod p is MAXIMALLY informative.
""")

print("\nDone.")
