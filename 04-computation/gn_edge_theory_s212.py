#!/usr/bin/env python3
"""
gn_edge_theory_s212.py — Exact formulas for T_n and relationships
opus-2026-03-22-S212

KEY DISCOVERY: T_n / V_n → m = C(n,2) as n → ∞.

This implies: asymptotically, each iso class contributes ~m transition orbits,
which means each class representative's m arc positions all produce DISTINCT
neighbor classes (no merging).

More precisely: T_n = V_n × m - correction_n.

Can we find the correction? And relate it to |E|?

THEOREM (exact): T_n = (1/n!) × Σ_{σ∈S_n} Fix(σ) × C(f(σ),2)
where Fix(σ) = # tournaments fixed by σ,
      f(σ) = # fixed points of σ,
and the Tb contribution is always 0 (even cycles kill Fix(σ)).
"""

import sys
from math import comb, factorial, gcd
from collections import Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  EXACT THEORY FOR T_n AND |E(G_n)|")
print("  opus-2026-03-22-S212")
print("=" * 80)

def partitions(n, max_part=None):
    if max_part is None:
        max_part = n
    if n == 0:
        yield []
        return
    for first in range(min(n, max_part), 0, -1):
        for rest in partitions(n - first, first):
            yield [first] + rest

def cycle_type_count(n, cycle_type):
    ct = Counter(cycle_type)
    result = factorial(n)
    for length, count in ct.items():
        result //= (length ** count) * factorial(count)
    return result

def fix_tournaments(cycle_type):
    for c in cycle_type:
        if c % 2 == 0:
            return 0
    k = len(cycle_type)
    exp = sum((c - 1) // 2 for c in cycle_type)
    for i in range(k):
        for j in range(i + 1, k):
            exp += gcd(cycle_type[i], cycle_type[j])
    return 2 ** exp

# ============================================================================
# EXACT COMPUTATION: T_n = V_n × m - D_n  where D_n is the "deficit"
# ============================================================================

print("\n  EXACT: T_n = V_n × m - D_n")
print("-" * 60)

for n in range(3, 14):
    m = comb(n, 2)
    nfact = factorial(n)

    V_sum = 0
    T_sum = 0

    for part in partitions(n):
        ct = cycle_type_count(n, part)
        ft = fix_tournaments(part)
        f = part.count(1)  # fixed points

        V_sum += ct * ft
        T_sum += ct * ft * comb(f, 2)

    V_n = V_sum // nfact
    T_n = T_sum // nfact
    D_n = V_n * m - T_n

    print(f"  n={n:2d}: V={V_n:>12d}, m={m:3d}, Vm={V_n*m:>14d}, "
          f"T={T_n:>14d}, D=Vm-T={D_n:>12d}")

# ============================================================================
# The deficit D_n
# ============================================================================

print("\n  DEFICIT D_n = V_n × m - T_n")
print("-" * 60)

D_vals = {}
V_vals = {}
T_vals = {}

for n in range(3, 14):
    m = comb(n, 2)
    nfact = factorial(n)
    V_sum = 0
    T_sum = 0
    for part in partitions(n):
        ct = cycle_type_count(n, part)
        ft = fix_tournaments(part)
        f = part.count(1)
        V_sum += ct * ft
        T_sum += ct * ft * comb(f, 2)
    V_n = V_sum // nfact
    T_n = T_sum // nfact
    D_n = V_n * m - T_n
    D_vals[n] = D_n
    V_vals[n] = V_n
    T_vals[n] = T_n

print(f"  D_n sequence: {[D_vals[n] for n in sorted(D_vals)]}")
print(f"  D_n / V_n: {[D_vals[n]/V_vals[n] for n in sorted(D_vals)]}")

# ============================================================================
# Can we express D_n in terms of V_n?
# ============================================================================

print("\n  D_n / V_n (should approach something?):")
for n in sorted(D_vals):
    ratio = D_vals[n] / V_vals[n]
    m = comb(n, 2)
    print(f"  n={n:2d}: D/V = {ratio:.6f}, m-D/V = {m - ratio:.6f}")

# ============================================================================
# EXACT FORMULA: D_n comes from non-identity permutations
# ============================================================================

print("\n  BREAKDOWN of D_n contribution by cycle type:")
print("-" * 60)

for n in range(3, 9):
    m = comb(n, 2)
    nfact = factorial(n)
    print(f"\n  n={n}:")

    for part in partitions(n):
        ct = cycle_type_count(n, part)
        ft = fix_tournaments(part)
        f = part.count(1)

        if ft == 0:
            continue

        # This cycle type's contribution to V: ct * ft / n!
        v_contrib = ct * ft / nfact
        # Contribution to T: ct * ft * C(f,2) / n!
        t_contrib = ct * ft * comb(f, 2) / nfact
        # Contribution to V*m: v_contrib * m
        vm_contrib = v_contrib * m
        # Contribution to deficit: vm_contrib - t_contrib
        d_contrib = vm_contrib - t_contrib
        # = v_contrib * (m - C(f,2))
        deficit_per = m - comb(f, 2)

        print(f"    {str(part):>20s}: ct={ct:>6d}, Fix={ft:>6d}, "
              f"f={f}, C(f,2)={comb(f,2):>3d}, "
              f"m-C(f,2)={deficit_per:>3d}, "
              f"v_contrib={v_contrib:.4f}, d_contrib={d_contrib:.4f}")

# ============================================================================
# KEY INSIGHT: For identity σ, f=n, C(f,2)=m, deficit=0.
# For non-identity σ with all-odd cycles, f < n, C(f,2) < m.
# The deficit comes entirely from non-identity automorphisms.
# ============================================================================

print("\n  KEY INSIGHT:")
print("  D_n = (1/n!) × Σ_{σ≠id, all-odd-cycles} ct(σ) × Fix(σ) × (m - C(f(σ),2))")
print("  The deficit measures how much non-trivial automorphisms reduce connectivity.")
print()

# ============================================================================
# EXACT FORMULA for D_n
# ============================================================================

print("  EXACT D_n formula:")
for n in range(3, 14):
    m = comb(n, 2)
    nfact = factorial(n)

    D_exact = 0
    for part in partitions(n):
        if part == [1] * n:
            continue  # skip identity
        ct = cycle_type_count(n, part)
        ft = fix_tournaments(part)
        if ft == 0:
            continue
        f = part.count(1)
        D_exact += ct * ft * (m - comb(f, 2))

    D_exact //= nfact

    print(f"  n={n:2d}: D_n = {D_exact:>12d}, check D_vals = {D_vals[n]:>12d}, "
          f"{'✓' if D_exact == D_vals[n] else '✗'}")

# ============================================================================
# RELATIONSHIP: T_n, D_n, V_n, and |E(G_n)|
# ============================================================================

E_known = {3: 1, 4: 5, 5: 30, 6: 290}
SL_w = {3: 12, 4: 144, 5: 1760, 6: 37920}

print("\n  RELATIONSHIP TABLE: T_n, D_n, V_n, |E|")
print("-" * 60)

print(f"  {'n':>3} {'V':>8} {'m':>4} {'T':>8} {'D':>8} {'|E|':>6} "
      f"{'T-2E':>8} {'D/(T-2E)':>10} {'SL_w/n!':>8}")

for n in sorted(E_known):
    V = V_vals[n]
    m = comb(n, 2)
    T = T_vals[n]
    D = D_vals[n]
    E = E_known[n]
    nfact = factorial(n)
    sl_orb = SL_w[n] / nfact

    t_minus_2e = T - 2*E
    ratio = D / t_minus_2e if t_minus_2e != 0 else float('inf')

    print(f"  {n:3d} {V:8d} {m:4d} {T:8d} {D:8d} {E:6d} "
          f"{t_minus_2e:8d} {ratio:10.4f} {sl_orb:8.2f}")

# T - 2E: this is the "excess orbits" = self-loop orbits + multi-weight edge orbits
# n=3: 4-2 = 2, n=4: 16-10 = 6, n=5: 88-60 = 28, n=6: 704-580 = 124

print(f"\n  T - 2|E| sequence: {[T_vals[n] - 2*E_known[n] for n in sorted(E_known)]}")
# 2, 6, 28, 124

# Is 2, 6, 28, 124 in OEIS? Or a known pattern?
# 2, 6, 28, 124...
# Diffs: 4, 22, 96
# Ratios: 3, 4.67, 4.43

# Actually, self-loop orbits:
# n=3: SL_w/n! = 2.0. T - 2E = 2. So self-loop orbits = 2 exactly (no multi-weight). ✓
# n=4: SL_w/n! = 6.0. T - 2E = 6. So self-loop orbits = 6 exactly. ✓
# n=5: SL_w/n! = 14.67. T - 2E = 28. So 28 - 14.67 = 13.33 multi-weight orbits?
# But orbits must be integers. The SL_w/n! = 14.67 means SL_orbits is NOT simply SL_w/n!.

# Actually, SL_w/n! = 1760/120 = 14.6667. This means SL_orbits ≠ 14.67.
# The issue is that some self-loop tournaments have |Aut| > 1, so their orbits
# have size < n!.

# Let me compute SL_orbits exactly from the data.
# At n=5: self-loop classes are 0,1,3,4,5,6,8.
# Class 0: sl_w=480, |Aut|=1. SL orbit size = n!/1 = 120. orbits = 480/120 = 4.
# Class 1: sl_w=40, |Aut|=3. But the orbit size for (T,e) is
#   n!/|Stab(T,e)| where Stab = {σ∈Aut(T): σ(e)=e}.
#   |Stab| divides |Aut|=3. If |Stab|=3, orbit = 120/3 = 40. orbits = 40/40 = 1.
#   If |Stab|=1, orbit = 120. orbits = 40/120 — not integer! So |Stab|=3.
# So class 1 has 1 self-loop orbit.
# Class 3: same as class 1: 1 orbit.
# Class 4: sl_w=360, |Aut|=1. orbits = 360/120 = 3.
# Class 5: sl_w=360, |Aut|=1. orbits = 3.
# Class 6: sl_w=240, |Aut|=1. orbits = 2.
# Class 8: sl_w=240, |Aut|=1. orbits = 2.

SL_orbits_exact = {
    3: 2,   # from class 0: 12/6 = 2
    4: 6,   # from class 0: 72/24=3, class 3: 72/24=3 → total 6
    5: 4+1+1+3+3+2+2,  # = 16
}
# n=5: SL_orbits = 16
# T - 2E = 88 - 60 = 28
# Multi-weight orbits = 28 - 16 = 12

# These 12 multi-weight orbits mean: 12 edges of G_5 have weight coming from
# 2 or more orbit types. Specifically: 12 extra orbits spread across 30 edges.
# Most edges have 2 orbits, some have 3 (contributing the extra 12 - some correction...)

# Actually: total cross-class orbits = T - SL_orbits = 88 - 16 = 72.
# And 72 orbits for 30 edges → 72/30 = 2.4 orbits per edge.
# If all edges had exactly 2 orbits: 60 orbits. But we have 72.
# Extra: 12 orbits → 12 edges have an extra orbit (3 instead of 2).
# Check: 18 edges × 2 orbits + 12 edges × 3 orbits = 36 + 36 = 72. ✓
# So 12 edges have 3 orbits and 18 edges have 2 orbits (at n=5).

print(f"\n  EXACT self-loop orbits:")
for n in sorted(SL_orbits_exact):
    T = T_vals[n]
    E = E_known[n]
    sl_orb = SL_orbits_exact[n]
    cross_orb = T - sl_orb
    ratio = cross_orb / E if E > 0 else 0
    extra = cross_orb - 2*E

    print(f"  n={n}: SL_orbits={sl_orb}, cross_orbits={cross_orb}, "
          f"|E|={E}, cross/E={ratio:.4f}, "
          f"extra_orbits={extra} (edges with multi-weight)")

# ============================================================================
# THEOREM: |E| = (T - SL_orbits - extra) / 2
# where extra = number of edges with >1 arc-position type per side
# ============================================================================

print(f"\n  EDGE FORMULA: |E| = (T - SL_orbits) / 2 - extra/2")
print(f"  where extra = edges with weight > n!")

# The formula is: 2|E| + extra = T - SL_orbits = cross_orbits
# |E| = (cross_orbits - extra) / 2

# The key unknown: what determines "extra"?

# extra/|E|: 0/1=0, 0/5=0, 12/30=0.4, ?/290
# At n=3,4: all edges have weight = n!, so extra = 0 → |E| = T/2 - SL/2
# At n=5: some edges have higher weight
# At n=6: need to compute

# ============================================================================
# The formula |E| = (T - SL_orbits) / 2 works EXACTLY at n=3,4
# ============================================================================

print(f"\n  VERIFICATION at n=3,4:")
for n in [3, 4]:
    T = T_vals[n]
    SL = SL_orbits_exact[n]
    E_pred = (T - SL) // 2
    E_actual = E_known[n]
    print(f"  n={n}: (T-SL)/2 = ({T}-{SL})/2 = {E_pred}, actual={E_actual} "
          f"{'✓' if E_pred == E_actual else '✗'}")

# n=3: (4-2)/2 = 1 ✓
# n=4: (16-6)/2 = 5 ✓
# n=5: (88-16)/2 = 36 ≠ 30 (because extra=12 → (88-16-12)/2=30 ✓)

print(f"\n  At n=5: (T-SL)/2 = (88-16)/2 = 36, but |E|=30. Deficit: 6 edges with 3 orbits.")
print(f"  (88-16-12)/2 = 30 ✓ (12 extra orbits from multi-weight edges)")

# ============================================================================
# SUMMARY OF WHAT WE KNOW
# ============================================================================

print(f"\n  SUMMARY:")
print(f"  T_n = (1/n!) Σ_σ Fix(σ) × C(f(σ), 2)  — computable for all n")
print(f"  V_n = (1/n!) Σ_σ Fix(σ)                 — computable for all n")
print(f"  D_n = V_n × m - T_n                     — the 'deficit'")
print(f"  SL_orbits: known at n≤5, need computation at n≥6")
print(f"  |E| = (T - SL_orbits) / 2 at n=3,4 (no multi-weight)")
print(f"  |E| = (T - SL_orbits - extra) / 2 at n≥5")
print(f"  The 'extra' term grows with n — this is the hard part.")

# Extended T_n values for reference:
print(f"\n  EXTENDED DATA:")
for n in range(3, 14):
    m = comb(n, 2)
    print(f"    n={n:2d}: V={V_vals[n]:>14d}, T={T_vals[n]:>16d}, "
          f"m={m:3d}, D={D_vals[n]:>14d}")

print(f"\n{'='*80}")
print(f"  DONE")
print(f"{'='*80}")
