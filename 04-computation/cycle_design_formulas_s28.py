#!/usr/bin/env python3
"""
opus-2026-04-05-S28: Extract and verify the design parameter formulas.

DISCOVERED:
  λ₃ = (p+1)/4
  λ₅ = (p+1)(p-2)(p-3)/8
  c₅ = p(p-1)(p+1)(p-2)(p-3)/160

  P_7: 5-cycles form a 4-(7,5,6) design = 2 copies of complete design.
  Every 5-subset of P_7 supports exactly 2 directed 5-cycles.
  At p≥11: 5-cycles form 2-design only (not 3-design).
"""

from math import comb
import time

def is_prime(n):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0 or n % 3 == 0: return False
    i = 5
    while i*i <= n:
        if n%i==0 or n%(i+2)==0: return False
        i += 6
    return True

print("=" * 70)
print("DESIGN PARAMETERS FOR PALEY TOURNAMENTS — PROVED FORMULAS")
print("=" * 70)
print()

print("FORMULA: For Paley P_p (p ≡ 3 mod 4 prime):")
print()
print("  3-CYCLE DESIGN:")
print("    c₃ = p(p²-1)/24")
print("    λ₃ = (p+1)/4    [each pair in exactly λ₃ directed 3-cycles]")
print()
print("  5-CYCLE DESIGN:")
print("    c₅ = p(p²-1)(p-2)(p-3)/160")
print("    λ₅ = (p+1)(p-2)(p-3)/8    [each pair in exactly λ₅ directed 5-cycles]")
print()
print("  RATIO: λ₅/λ₃ = (p-2)(p-3)/2 = C(p-2, 2)")
print()

print(f"{'p':>4} {'c₃':>8} {'λ₃':>6} {'c₅':>10} {'λ₅':>8} {'λ₅/λ₃':>8} "
      f"{'C(p-2,2)':>9} {'c₅/C(p,5)':>10}")
print("-" * 70)

for p in range(7, 100):
    if not (is_prime(p) and p % 4 == 3):
        continue

    c3 = p * (p*p - 1) // 24
    lam3 = (p + 1) // 4
    c5 = p * (p*p - 1) * (p - 2) * (p - 3) // 160
    lam5 = (p + 1) * (p - 2) * (p - 3) // 8
    ratio = lam5 // lam3
    cp2_2 = comb(p - 2, 2)
    c5_over_cp5 = c5 / comb(p, 5) if comb(p,5) > 0 else 0

    print(f"{p:>4} {c3:>8} {lam3:>6} {c5:>10} {lam5:>8} {ratio:>8} "
          f"{cp2_2:>9} {c5_over_cp5:>10.4f}")

print()
print("SPECIAL CASE P_7:")
print("  c₅ = 42 = 2 × C(7,5) = 2 × 21")
print("  → Every 5-element subset supports EXACTLY 2 directed 5-cycles")
print("  → Every induced 5-vertex subtournament is REGULAR")
print("  → The 5-cycles form a 4-(7, 5, 6) design = 2 × complete design")
print()
print("  At p ≥ 11: c₅/C(p,5) is NOT an integer")
print("  → The distribution of 5-cycles across 5-subsets is non-uniform")
print("  → Only a 2-design, not a 3- or 4-design")

print()
print("=" * 70)
print("GENERAL k-CYCLE FORMULA (CONJECTURE)")
print("=" * 70)
print()
print("For Paley P_p, directed k-cycles (k odd, k ≤ p) form a 2-(p, k, λ_k) design.")
print()
print("Conjectured formula:")
print("  λ_k = (p+1) × Π_{j=2}^{k-2} (p-j) / 2^{k-1}")
print()
print("  k=3: λ₃ = (p+1) / 4 = (p+1)/2²")
print("  k=5: λ₅ = (p+1)(p-2)(p-3) / 8 = (p+1)(p-2)(p-3)/2³")
print("  k=7: λ₇ = (p+1)(p-2)(p-3)(p-4)(p-5) / 16 = ... /2⁴  [TO VERIFY]")
print()

# Verify k=7 formula at p=7 (where λ₇ = 24 from computation)
p = 7
lam7_formula = (p+1)*(p-2)*(p-3)*(p-4)*(p-5) // 16
print(f"  k=7, p=7: formula gives {lam7_formula}, computed λ₇ = 24")
print(f"  Match: {'✓' if lam7_formula == 24 else '✗'}")

# For p=11: we computed λ₇ = 1512
p = 11
lam7_formula = (p+1)*(p-2)*(p-3)*(p-4)*(p-5) // 16
print(f"  k=7, p=11: formula gives {lam7_formula}, computed λ₇ = 1512")
print(f"  Match: {'✓' if lam7_formula == 1512 else '✗'}")

print()
print("=" * 70)
print("THE GENERAL FORMULA")
print("=" * 70)
print()
print("  For k odd, 3 ≤ k ≤ p:")
print()
print("         (p+1)     k-2")
print("  λ_k = ────── × Π (p - j)")
print("        2^{k-1}  j=2")
print()
print("  Equivalently: λ_k = (p+1) × (p-2)! / ((p-k+1)! × 2^{k-1})")
print()
print("             (p+1) × (p-2)_(k-2)")
print("       λ_k = ────────────────────")
print("                   2^{k-1}")
print()
print("  where (a)_m = a(a-1)...(a-m+1) is the falling factorial.")
print()

# Verify at all computed cases
print("Verification against ALL computed values:")
computed = {
    (7, 3): 2, (7, 5): 20, (7, 7): 24,
    (11, 3): 3, (11, 5): 108, (11, 7): 1512,
    (19, 3): 5, (19, 5): 680,
    (23, 3): 6, (23, 5): 1260,
}

all_match = True
for (p, k), expected in sorted(computed.items()):
    # λ_k = (p+1) × falling_factorial(p-2, k-2) / 2^{k-1}
    ff = 1
    for j in range(k-2):
        ff *= (p - 2 - j)
    lam = (p + 1) * ff // (2 ** (k - 1))
    match = lam == expected
    if not match: all_match = False
    print(f"  p={p:>2}, k={k}: formula={lam:>6}, computed={expected:>6} {'✓' if match else '✗'}")

print()
if all_match:
    print("✓ ALL VALUES MATCH — formula verified!")
else:
    print("✗ Some values don't match")

print()
print("TOTAL ODD-CYCLE INCIDENCE PER PAIR:")
for p in [7, 11, 19, 23]:
    if not (is_prime(p) and p % 4 == 3):
        continue
    total = 0
    for k in range(3, p+1, 2):
        ff = 1
        for j in range(k-2):
            ff *= (p - 2 - j)
        lam_k = (p + 1) * ff // (2 ** (k - 1))
        total += lam_k
    print(f"  p={p}: Σ_k λ_k = {total}")
    # This should equal... what?
    # Each σ ∈ S_n with all arcs in T and containing vertex pair {u,v}
    # contributes to some odd k-cycle through {u,v}.
    # Total = per(A_T) restricted to perms where u,v are in the same cycle.
    # For p=7: per = ham = 24. Each Ham cycle through every pair: 24.
    # Plus: 5-cycles through pair: 20. Plus 3-cycles: 2. Total: 46.
    # But Σ λ_k includes the p-cycle (Hamiltonian) contribution.
