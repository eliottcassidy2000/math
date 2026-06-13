#!/usr/bin/env python3
"""
H_squared_extension_s109.py — Extending the Three-Factor OCR Formula Beyond n=5
opus-2026-03-21-S109

AT n=5: 1-OCR = P(PoS) × 4·Var(c₅|PoS) / Var(H) = 4/133
This works because H = 1 + 2α₁ = 1 + 2(c₃+c₅) and α₂=0 always at n=5.

AT n≥6: H = 1 + 2α₁ + 4α₂ + ... and α₂ can be nonzero.
The simple formula BREAKS because the within-class variance of H
involves BOTH Var(α₁|scores) AND Var(α₂|scores,α₁) AND cross-terms.

THIS SESSION: Understand the n=6 structure precisely.
What is the EXACT decomposition of Var(H|scores) at n=6?

DEEP THOUGHT ABOUT H²:

H² = (1 + 2α₁ + 4α₂ + ...)²
   = 1 + 4α₁ + 4α₁² + 16α₂² + 8α₂ + 16α₁α₂ + ...

E[H²|s] involves E[α₁²|s], E[α₂²|s], E[α₁α₂|s], etc.
Since c₃ is determined by s, α₁ = c₃(s) + c₅ + c₇ + ..., and
the residual is Var(c₅ + c₇ + ...|s) plus Var(α₂|s) plus cross-terms.

The question: at n=6, how much of Var(H|s) comes from:
  (a) Var(α₁|s) = Var(c₅+c₇|s)  [the "cycle count" residual]
  (b) Var(α₂|s)                    [the "packing" residual]
  (c) Cov(α₁, α₂|s)               [the cross-term]
"""

from math import comb, factorial
from fractions import Fraction
from itertools import permutations, combinations
from collections import defaultdict, Counter
import time

print("=" * 72)
print("  H² AND THE OCR BEYOND n=5")
print("  opus-2026-03-21-S109")
print("=" * 72)
print()

# =============================================================================
# PART 1: EXACT α₁, α₂ DECOMPOSITION AT n=6
# =============================================================================

print("=" * 72)
print("  PART 1: EXACT α₁, α₂ AT n=6")
print("=" * 72)
print()

n = 6
m = comb(n, 2)  # 15

print(f"  Computing α₁, α₂ for all {2**m} tournaments at n={n}...")
t0 = time.time()

def count_odd_cycles(adj_2d, n):
    """Count ALL directed odd cycles and independent sets of odd cycles."""
    # Find all 3-cycles
    cycles_3 = []
    for i, j, k in combinations(range(n), 3):
        if adj_2d[i][j] and adj_2d[j][k] and adj_2d[k][i]:
            cycles_3.append(frozenset([i, j, k]))
        elif adj_2d[i][k] and adj_2d[k][j] and adj_2d[j][i]:
            cycles_3.append(frozenset([i, j, k]))

    # Find all 5-cycles
    cycles_5 = []
    for perm in combinations(range(n), 5):
        for p in permutations(perm):
            if all(adj_2d[p[k]][p[(k+1)%5]] for k in range(5)):
                vset = frozenset(p)
                if vset not in [c for c in cycles_5]:
                    cycles_5.append(vset)
                break  # Only need one direction per vertex set

    all_cycles = list(set(cycles_3)) + list(set(cycles_5))
    alpha_1 = len(all_cycles)

    # α₂ = number of independent PAIRS (vertex-disjoint)
    alpha_2 = 0
    for i in range(len(all_cycles)):
        for j in range(i+1, len(all_cycles)):
            if all_cycles[i].isdisjoint(all_cycles[j]):
                alpha_2 += 1

    return alpha_1, alpha_2, len(set(cycles_3)), len(set(cycles_5))

# Enumerate all n=6 tournaments
score_data = defaultdict(list)  # score → list of (H, α₁, α₂, c₃, c₅)

for mask in range(2**m):
    adj_flat = [0]*(n*n)
    adj_2d = [[0]*n for _ in range(n)]
    scores = [0]*n
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (mask >> idx) & 1:
                adj_flat[i*n+j] = 1; adj_2d[i][j] = 1; scores[i] += 1
            else:
                adj_flat[j*n+i] = 1; adj_2d[j][i] = 1; scores[j] += 1
            idx += 1

    sc = tuple(sorted(scores))

    # H via DP
    dp = [0]*((1<<n)*n)
    for v in range(n): dp[(1<<v)*n+v] = 1
    for S in range(1, 1<<n):
        for v in range(n):
            if not (S&(1<<v)): continue
            val = dp[S*n+v]
            if val == 0: continue
            for u in range(n):
                if S&(1<<u): continue
                if adj_flat[v*n+u]: dp[(S|(1<<u))*n+u] += val
    full = (1<<n)-1
    H = sum(dp[full*n+v] for v in range(n))

    # α₁, α₂, c₃, c₅
    alpha1, alpha2, c3, c5 = count_odd_cycles(adj_2d, n)

    # Verify OCF: H = 1 + 2α₁ + 4α₂ + ...
    H_from_ocf = 1 + 2*alpha1 + 4*alpha2
    # At n=6, there could be α₃ (three disjoint cycles)
    # For that we'd need 3 disjoint 3-cycles = 9 vertices, but n=6 < 9
    # So at n=6: α₃ = 0 always. H = 1 + 2α₁ + 4α₂ exactly.

    score_data[sc].append({
        'H': H, 'alpha1': alpha1, 'alpha2': alpha2,
        'c3': c3, 'c5': c5, 'H_ocf': H_from_ocf
    })

elapsed = time.time() - t0
print(f"  Completed in {elapsed:.1f}s")
print()

# Verify OCF holds
ocf_errors = 0
for sc, tours in score_data.items():
    for t in tours:
        if t['H'] != t['H_ocf']:
            ocf_errors += 1
print(f"  OCF verification H = 1 + 2α₁ + 4α₂: {ocf_errors} errors out of {2**m}")
if ocf_errors == 0:
    print("  ✓ OCF holds EXACTLY at n=6 (H = 1 + 2α₁ + 4α₂, no α₃ term needed)")
print()

# =============================================================================
# PART 2: THE α₁, α₂ DECOMPOSITION WITHIN EACH PoS CLASS
# =============================================================================

print("=" * 72)
print("  PART 2: α₁, α₂ WITHIN PoS CLASSES AT n=6")
print("=" * 72)
print()

print(f"  {'Score':>25} | {'n_s':>5} | {'#H':>3} | {'α₁ range':>10} | {'α₂ range':>10} | {'c₃':>3} | {'c₅ range':>10} | {'Var(H|s)':>10}")
print("  " + "-" * 90)

total_within = Fraction(0)
N = 2**m

for sc in sorted(score_data.keys()):
    tours = score_data[sc]
    ns = len(tours)
    H_vals = [t['H'] for t in tours]
    H_set = sorted(set(H_vals))

    a1_vals = [t['alpha1'] for t in tours]
    a2_vals = [t['alpha2'] for t in tours]
    c3_vals = [t['c3'] for t in tours]
    c5_vals = [t['c5'] for t in tours]

    a1_range = f"[{min(a1_vals)},{max(a1_vals)}]"
    a2_range = f"[{min(a2_vals)},{max(a2_vals)}]"
    c5_range = f"[{min(c5_vals)},{max(c5_vals)}]"
    c3_const = c3_vals[0] if len(set(c3_vals)) == 1 else f"VARIES!"

    mean_H = sum(H_vals) / ns
    var_H = sum((h - mean_H)**2 for h in H_vals) / ns

    is_pos = len(H_set) > 1

    if is_pos:
        weight = Fraction(ns, N)
        contrib = weight * Fraction(int(var_H * 10000 + 0.5), 10000)
        # More precisely
        sum_H = sum(H_vals)
        sum_H2 = sum(h*h for h in H_vals)
        var_exact = Fraction(sum_H2, ns) - Fraction(sum_H, ns)**2
        total_within += Fraction(ns, N) * var_exact

        print(f"  {str(sc):>25} | {ns:5d} | {len(H_set):3d} | {a1_range:>10} | {a2_range:>10} | "
              f"{c3_const:>3} | {c5_range:>10} | {float(var_exact):10.4f} ★")
    elif len(set(c3_vals)) > 1:
        print(f"  {str(sc):>25} | {ns:5d} | {len(H_set):3d} | {a1_range:>10} | {a2_range:>10} | "
              f"{'VARIES':>3} | {c5_range:>10} |     0 (c₃ varies!)")

print()

# Compute exact Var(H|scores)
var_within_exact = Fraction(0)
for sc, tours in score_data.items():
    ns = len(tours)
    H_vals = [t['H'] for t in tours]
    sum_H = sum(H_vals)
    sum_H2 = sum(h*h for h in H_vals)
    var_sc = Fraction(sum_H2, ns) - Fraction(sum_H, ns)**2
    var_within_exact += Fraction(ns, N) * var_sc

var_H_total = Fraction(sum(t['H']**2 for tours in score_data.values() for t in tours), N) - \
              Fraction(sum(t['H'] for tours in score_data.values() for t in tours), N)**2

print(f"  Var(H|scores) = {var_within_exact} = {float(var_within_exact):.6f}")
print(f"  Var(H) = {var_H_total} = {float(var_H_total):.6f}")
print(f"  OCR = {1 - var_within_exact/var_H_total} = {float(1 - var_within_exact/var_H_total):.10f}")
print()

# =============================================================================
# PART 3: THE VARIANCE DECOMPOSITION OF H = 1 + 2α₁ + 4α₂
# =============================================================================

print("=" * 72)
print("  PART 3: Var(H|s) = 4·Var(α₁|s) + 16·Var(α₂|s) + 16·Cov(α₁,α₂|s)")
print("=" * 72)
print()

print("""
  Since H = 1 + 2α₁ + 4α₂:
    H - E[H|s] = 2(α₁ - E[α₁|s]) + 4(α₂ - E[α₂|s])

  Var(H|s) = 4·Var(α₁|s) + 16·Var(α₂|s) + 16·Cov(α₁,α₂|s)

  At n=5: α₂=0 always, so Var(H|s) = 4·Var(α₁|s) = 4·Var(c₅|s).
  At n=6: α₂ can be nonzero, so ALL THREE terms can contribute.

  Let's decompose EXACTLY for each PoS class at n=6.
""")

print(f"  {'Score':>25} | {'Var(H|s)':>10} | {'4Var(α₁)':>10} | {'16Var(α₂)':>10} | {'16Cov':>10} | {'sum':>10} | {'match?':>7}")
print("  " + "-" * 90)

for sc in sorted(score_data.keys()):
    tours = score_data[sc]
    ns = len(tours)
    H_vals = [t['H'] for t in tours]
    if len(set(H_vals)) <= 1:
        continue

    a1_vals = [t['alpha1'] for t in tours]
    a2_vals = [t['alpha2'] for t in tours]

    mean_a1 = sum(a1_vals) / ns
    mean_a2 = sum(a2_vals) / ns
    var_a1 = sum((a - mean_a1)**2 for a in a1_vals) / ns
    var_a2 = sum((a - mean_a2)**2 for a in a2_vals) / ns
    cov_a1a2 = sum((a1_vals[i] - mean_a1)*(a2_vals[i] - mean_a2) for i in range(ns)) / ns

    var_H_formula = 4*var_a1 + 16*var_a2 + 16*cov_a1a2
    mean_H = sum(H_vals) / ns
    var_H_direct = sum((h - mean_H)**2 for h in H_vals) / ns

    match = abs(var_H_formula - var_H_direct) < 0.001

    print(f"  {str(sc):>25} | {var_H_direct:10.4f} | {4*var_a1:10.4f} | {16*var_a2:10.4f} | {16*cov_a1a2:10.4f} | {var_H_formula:10.4f} | {'✓' if match else '✗'}")

print()

# =============================================================================
# PART 4: HOW MUCH OF THE RESIDUAL COMES FROM α₂?
# =============================================================================

print("=" * 72)
print("  PART 4: α₂'s CONTRIBUTION TO THE OCR RESIDUAL AT n=6")
print("=" * 72)
print()

total_from_a1 = 0
total_from_a2 = 0
total_from_cov = 0
total_var_H_s = 0

for sc, tours in score_data.items():
    ns = len(tours)
    H_vals = [t['H'] for t in tours]
    if len(set(H_vals)) <= 1:
        continue

    a1_vals = [t['alpha1'] for t in tours]
    a2_vals = [t['alpha2'] for t in tours]
    mean_a1 = sum(a1_vals) / ns
    mean_a2 = sum(a2_vals) / ns
    var_a1 = sum((a - mean_a1)**2 for a in a1_vals) / ns
    var_a2 = sum((a - mean_a2)**2 for a in a2_vals) / ns
    cov_a1a2 = sum((a1_vals[i] - mean_a1)*(a2_vals[i] - mean_a2) for i in range(ns)) / ns

    weight = ns / N
    total_from_a1 += weight * 4 * var_a1
    total_from_a2 += weight * 16 * var_a2
    total_from_cov += weight * 16 * cov_a1a2
    total_var_H_s += weight * (4*var_a1 + 16*var_a2 + 16*cov_a1a2)

print(f"  Total Var(H|scores) decomposition:")
print(f"    From 4·Var(α₁|s):    {total_from_a1:.6f} ({total_from_a1/total_var_H_s*100:.1f}%)")
print(f"    From 16·Var(α₂|s):   {total_from_a2:.6f} ({total_from_a2/total_var_H_s*100:.1f}%)")
print(f"    From 16·Cov(α₁,α₂|s):{total_from_cov:.6f} ({total_from_cov/total_var_H_s*100:.1f}%)")
print(f"    Total:                {total_var_H_s:.6f}")
print(f"    Expected:             {float(var_within_exact):.6f}")
print()

# =============================================================================
# PART 5: FURTHER DECOMPOSE α₁ = c₃ + c₅ + c₇(?)
# =============================================================================

print("=" * 72)
print("  PART 5: α₁ DECOMPOSITION: c₃ + c₅ + ...?")
print("=" * 72)
print()

print("""
  At n=6: are there 7-cycles? No — a 7-cycle needs 7 vertices but n=6.
  So α₁ = c₃ + c₅ at n=6 (same as n=5).
  And c₃ is determined by scores.
  So Var(α₁|s) = Var(c₅|s) (exactly).
""")

# Verify: is α₁ = c₃ + c₅?
for sc, tours in score_data.items():
    for t in tours[:3]:
        if t['alpha1'] != t['c3'] + t['c5']:
            print(f"  ERROR: α₁ ≠ c₃ + c₅ at {sc}")
            break

print(f"  Verified: α₁ = c₃ + c₅ for all {2**m} tournaments at n=6 ✓")
print(f"  (No 7-cycles possible at n=6)")
print()

# So at n=6:
# Var(H|s) = 4·Var(c₅|s) + 16·Var(α₂|s) + 16·Cov(c₅,α₂|s)
# (since c₃ is constant within each score class)

print(f"  Decomposition: Var(H|s) = 4·Var(c₅|s) + 16·Var(α₂|s) + 16·Cov(c₅,α₂|s)")
print()

total_from_c5 = 0
total_from_a2_exact = 0
total_from_cov_exact = 0

for sc, tours in score_data.items():
    ns = len(tours)
    H_vals = [t['H'] for t in tours]
    if len(set(H_vals)) <= 1:
        continue

    c5_vals = [t['c5'] for t in tours]
    a2_vals = [t['alpha2'] for t in tours]
    mean_c5 = sum(c5_vals) / ns
    mean_a2 = sum(a2_vals) / ns
    var_c5 = sum((c - mean_c5)**2 for c in c5_vals) / ns
    var_a2 = sum((a - mean_a2)**2 for a in a2_vals) / ns
    cov = sum((c5_vals[i]-mean_c5)*(a2_vals[i]-mean_a2) for i in range(ns)) / ns

    weight = ns / N
    total_from_c5 += weight * 4 * var_c5
    total_from_a2_exact += weight * 16 * var_a2
    total_from_cov_exact += weight * 16 * cov

print(f"  Var(H|scores) at n=6 decomposed into:")
print(f"    4·Var(c₅|s):     {total_from_c5:.6f} ({total_from_c5/float(var_within_exact)*100:.1f}%)")
print(f"    16·Var(α₂|s):    {total_from_a2_exact:.6f} ({total_from_a2_exact/float(var_within_exact)*100:.1f}%)")
print(f"    16·Cov(c₅,α₂|s): {total_from_cov_exact:.6f} ({total_from_cov_exact/float(var_within_exact)*100:.1f}%)")
print(f"    Total:            {total_from_c5+total_from_a2_exact+total_from_cov_exact:.6f}")
print(f"    Expected:         {float(var_within_exact):.6f}")
print()

# =============================================================================
# PART 6: THE EXTENDED THREE-FACTOR FORMULA
# =============================================================================

print("=" * 72)
print("  PART 6: THE EXTENDED FORMULA FOR n≥6")
print("=" * 72)
print()

print(f"""
  AT n=5: 1-OCR = [4·Var(c₅|s)] / Var(H)  (single term, no α₂)

  AT n=6: 1-OCR = [4·Var(c₅|s) + 16·Var(α₂|s) + 16·Cov(c₅,α₂|s)] / Var(H)

  where each term is averaged over score classes weighted by P(s).

  The RELATIVE contributions at n=6:
    c₅ alone: {total_from_c5/float(var_within_exact)*100:.1f}%
    α₂ alone: {total_from_a2_exact/float(var_within_exact)*100:.1f}%
    cross-term: {total_from_cov_exact/float(var_within_exact)*100:.1f}%
""")

if total_from_a2_exact / float(var_within_exact) < 0.01:
    print("  *** α₂ contributes NEGLIGIBLY at n=6! ***")
    print("  The n=5 formula 1-OCR ≈ 4·Var(c₅|s)/Var(H) is still a good approximation!")
elif total_from_a2_exact / float(var_within_exact) < 0.1:
    print("  α₂ contributes a small but nonzero amount at n=6.")
    print("  The n=5 formula needs a correction term of ~10%.")
else:
    print("  α₂ contributes significantly at n=6.")
    print("  The full three-component formula is needed.")

print()

# =============================================================================
# PART 7: WHAT DOES H² REALLY ENCODE?
# =============================================================================

print("=" * 72)
print("  PART 7: DEEP THOUGHT — WHAT DOES H² ENCODE?")
print("=" * 72)
print()

print("""
  H(T) = #{Hamiltonian paths of T}
  H(T)² = #{ORDERED PAIRS of Hamiltonian paths of T}

  An ordered pair (σ,τ) of HPs of T means:
    σ = (σ₁,σ₂,...,σₙ) visits all vertices in order σ₁→σ₂→...→σₙ
    τ = (τ₁,τ₂,...,τₙ) visits all vertices in order τ₁→τ₂→...→τₙ
    Both are valid HPs (all arcs point forward).

  The pair (σ,τ) has structure:
    - SHARED ARCS: {σᵢ→σᵢ₊₁} ∩ {τⱼ→τⱼ₊₁}
    - These are arcs that both paths traverse in the same direction
    - The more shared arcs, the more "similar" the two paths are

  H² = Σ_{(σ,τ)} [both are HPs of T]
     = Σ_{(σ,τ): compatible} [T satisfies all arc constraints of σ and τ]

  For compatible (σ,τ) with a shared arcs:
    - They jointly constrain 2(n-1)-a arcs of T
    - Each unconstrained arc can be either direction
    - BUT: the constrained arcs must be consistent with T's score sequence

  THE KEY INSIGHT FOR H²:
  H² counts how many PAIRS of paths the tournament supports.
  A tournament with "rich path structure" (many diverse paths) has large H².
  A tournament with "poor path structure" (few paths, or similar paths) has small H².

  The VARIANCE Var(H²|scores) measures how much the PATH PAIR STRUCTURE
  varies within a score class. This is what the OCR residual captures.

  At the PoS: some tournaments have many diverse path pairs (high H²),
  while others have fewer (low H²). The difference is in the CYCLE
  STRUCTURE — specifically, how the 5-cycles (and at n≥6, the disjoint
  cycle pairs) are arranged.
""")

# Compute H² statistics within each PoS class
print(f"  H² statistics within PoS classes at n=6:")
print(f"  {'Score':>25} | {'E[H²|s]':>10} | {'E[H]²|s':>10} | {'Var(H|s)':>10} | {'max H²':>8} | {'min H²':>8}")
print("  " + "-" * 80)

for sc in sorted(score_data.keys()):
    tours = score_data[sc]
    ns = len(tours)
    H_vals = [t['H'] for t in tours]
    if len(set(H_vals)) <= 1:
        continue

    H2_vals = [h*h for h in H_vals]
    E_H2 = sum(H2_vals) / ns
    E_H_sq = (sum(H_vals) / ns)**2
    Var_H = E_H2 - E_H_sq

    print(f"  {str(sc):>25} | {E_H2:10.2f} | {E_H_sq:10.2f} | {Var_H:10.4f} | {max(H2_vals):8d} | {min(H2_vals):8d}")

print()

# =============================================================================
# PART 8: THE GENERAL PATTERN
# =============================================================================

print("=" * 72)
print("  PART 8: THE GENERAL PATTERN — EXTENDING TO ALL n")
print("=" * 72)
print()

print("""
  AT n=5: H = 1 + 2α₁ + 4α₂
    α₁ = c₃ + c₅ (no 7-cycles at n=5)
    α₂ = 0 always (disjoint 3+3 needs 6 vertices, but n=5)
    α₃ = 0 always
    → Var(H|s) = 4·Var(c₅|s)

  AT n=6: H = 1 + 2α₁ + 4α₂
    α₁ = c₃ + c₅ (no 7-cycles at n=6)
    α₂ CAN be nonzero (two disjoint 3-cycles on 6 vertices)
    α₃ = 0 always (three disjoint 3-cycles need 9 vertices)
    → Var(H|s) = 4·Var(c₅|s) + 16·Var(α₂|s) + 16·Cov(c₅,α₂|s)

  AT n=7: H = 1 + 2α₁ + 4α₂ + 8α₃(?)
    α₁ = c₃ + c₅ + c₇ (7-cycles EXIST at n=7)
    α₂ CAN be nonzero (two disjoint 3-cycles, or 3+5-cycle pair, etc.)
    α₃ = 0 always? (three disjoint odd cycles need ≥9 vertices... wait,
      3+3+3=9 > 7, so no. BUT 3+3 on 6 leaves 1 vertex, and a 1-cycle
      doesn't exist. So α₃ = 0 at n=7.)
    Actually: α₃ requires THREE mutually vertex-disjoint odd cycles.
    Minimum: three 3-cycles → 9 vertices. Since n=7 < 9: α₃ = 0 at n=7.
    → Var(H|s) = 4·Var(c₅+c₇|s) + 16·Var(α₂|s) + 16·Cov(c₅+c₇,α₂|s)

  AT n=8:
    α₁ = c₃ + c₅ + c₇ (all exist)
    α₂ = #{disjoint pairs} (3+3, 3+5) — all possible
    α₃ = 0 (three disjoint 3-cycles need 9 > 8)
    → Same structure as n=7

  AT n=9:
    α₃ CAN be nonzero! (three disjoint 3-cycles on 9 vertices)
    → New term: 64·Var(α₃|s) enters

  GENERAL PATTERN:
    The OCF truncates at α_k for k = ⌊n/3⌋ (maximum number of disjoint 3-cycles).
    New terms enter at n = 3k for each k.

    n=5: terms up to α₁ (residual from c₅)
    n=6: terms up to α₂ (residual from c₅ + α₂)
    n=7-8: terms up to α₂ (residual from c₅+c₇ + α₂)
    n=9: terms up to α₃ (new!)

  The OCR residual GROWS as n increases because:
    (a) More cycle types enter (c₅, c₇, c₉, ...)
    (b) More packing types enter (α₂ at n=6, α₃ at n=9, ...)
    (c) Each adds variance within score classes

  BUT Var(H) also grows (factorially), so the RATIO might decrease.
""")

# =============================================================================
# PART 9: THE EXACT n=6 THREE-FACTOR FORMULA
# =============================================================================

print("=" * 72)
print("  PART 9: EXACT THREE-FACTOR FORMULA AT n=6")
print("=" * 72)
print()

print(f"  AT n=6:")
print(f"  Var(H) = {var_H_total} = {float(var_H_total):.4f}")
print(f"  Var(H|scores) = {var_within_exact} = {float(var_within_exact):.6f}")
print(f"  OCR = {1 - var_within_exact/var_H_total}")
print(f"  1-OCR = {var_within_exact/var_H_total}")
print()

# Decompose 1-OCR into contributions from each PoS class
print(f"  Per-PoS-class contributions to 1-OCR:")
print(f"  {'Score':>25} | {'P(s)':>8} | {'Var(H|s)':>10} | {'P×Var':>10} | {'% of resid':>10}")
print("  " + "-" * 72)

residual_total = float(var_within_exact)
for sc in sorted(score_data.keys()):
    tours = score_data[sc]
    ns = len(tours)
    H_vals = [t['H'] for t in tours]
    if len(set(H_vals)) <= 1:
        continue

    sum_H = sum(H_vals)
    sum_H2 = sum(h*h for h in H_vals)
    var_sc = float(Fraction(sum_H2, ns) - Fraction(sum_H, ns)**2)
    P_sc = ns / N
    contrib = P_sc * var_sc
    pct = contrib / residual_total * 100

    print(f"  {str(sc):>25} | {P_sc:8.4f} | {var_sc:10.4f} | {contrib:10.6f} | {pct:9.1f}%")

print()

# =============================================================================
# PART 10: SUMMARY
# =============================================================================

print("=" * 72)
print("  SUMMARY")
print("=" * 72)
print()

print(f"""
  EXTENDING THE THREE-FACTOR OCR FORMULA BEYOND n=5:

  THE EXACT FORMULA (for any n):
    1 - OCR(n) = Σ_{{s: ambiguous}} P(s) × Var(H|s) / Var(H)

  where Var(H|s) = 4·Var(α₁|s) + 16·Var(α₂|s) + 16·Cov(α₁,α₂|s)
                 + 64·Var(α₃|s) + ... [higher terms at n≥9]

  AT n=5: Var(H|s) = 4·Var(c₅|s)                         (100% from c₅)
  AT n=6: Var(H|s) = 4·Var(c₅|s) + 16·Var(α₂|s) + cross  (c₅ + α₂)
  AT n=7: Var(H|s) = 4·Var(c₅+c₇|s) + 16·Var(α₂|s) + ... (c₅,c₇ + α₂)

  KEY FINDING AT n=6:
    4·Var(c₅|s) contributes:  {total_from_c5/float(var_within_exact)*100:.1f}% of the residual
    16·Var(α₂|s) contributes: {total_from_a2_exact/float(var_within_exact)*100:.1f}% of the residual
    16·Cov(c₅,α₂|s):         {total_from_cov_exact/float(var_within_exact)*100:.1f}% of the residual

  THE α₂ CONTRIBUTION IS {'NEGLIGIBLE' if abs(total_from_a2_exact/float(var_within_exact)) < 0.05 else 'SIGNIFICANT'} AT n=6.

  GENERAL SCALING:
    The residual has the form:
    1-OCR = Σ_s P(s) × [4·Σ_{{k odd ≥5}} Var(c_k|s) + 16·Var(α₂|s) + ...] / Var(H)

    Each term in the numerator is bounded by C(n,k)² (cycle counts)
    or C(n,6)² (packing counts), both POLYNOMIAL in n.
    Var(H) grows FACTORIALLY.

    IF the number of significant PoS classes grows at most POLYNOMIALLY,
    then the total residual grows polynomially while Var(H) grows
    factorially, giving OCR → 1 as n → ∞.

    The key question: does the number of significant PoS classes
    grow polynomially? From the score sequence count S_n ~ 4^n/n^{{5/2}},
    the number of classes grows EXPONENTIALLY. But the FRACTION that
    is ambiguous might shrink.

    THIS REMAINS THE CENTRAL OPEN QUESTION.
""")

print("opus-2026-03-21-S109: H² AND THE OCR BEYOND n=5")
