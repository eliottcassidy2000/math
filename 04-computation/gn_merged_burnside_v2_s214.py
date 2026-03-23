#!/usr/bin/env python3
"""
gn_merged_burnside_v2_s214.py — CORRECTED Burnside for G_n/Z_2
opus-2026-03-22-S214

FIX: The anti-automorphism formula was wrong.
For Fix_anti(σ) = #{T : σ(T) = T^op}, the correct analysis is:

σ acts on undirected arcs {i,j} → {σi, σj}. The anti-aut condition
says: T(σi, σj) = 1 - T(i,j). So along each orbit of undirected arcs,
the value FLIPS each step.

For orbit length 1 (fixed arc {i,j}):
  - If σi=i, σj=j (both fixed): T(i,j) = 1-T(i,j) → IMPOSSIBLE
  - If σi=j, σj=i (swap): T(j,i) = 1-T(i,j), always true → 1 free choice

For orbit length L > 1:
  - Values alternate: v, 1-v, v, 1-v, ...
  - After L steps must return to v → need L EVEN
  - If L even: 1 free choice
  - If L odd: IMPOSSIBLE → Fix_anti = 0

So Fix_anti = 0 if:
  (a) σ has ≥ 2 fixed points (creates a both-fixed arc), OR
  (b) any undirected arc orbit has ODD length > 1

Otherwise: Fix_anti = 2^(# orbits with valid structure)
"""

import sys
from math import comb, factorial, gcd
from itertools import permutations
from collections import Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  CORRECTED BURNSIDE FOR G_n/Z_2")
print("  opus-2026-03-22-S214")
print("=" * 80)

# ============================================================================
# HELPERS
# ============================================================================

def gen_partitions(n, max_part=None):
    if max_part is None: max_part = n
    if n == 0: yield []; return
    for first in range(min(n, max_part), 0, -1):
        for rest in gen_partitions(n - first, first):
            yield [first] + rest

def cycle_type_count(n, ct):
    c = Counter(ct)
    r = factorial(n)
    for l, k in c.items(): r //= (l**k) * factorial(k)
    return r

def fix_tournaments(ct):
    for c in ct:
        if c % 2 == 0: return 0
    exp = sum((c-1)//2 for c in ct)
    for i in range(len(ct)):
        for j in range(i+1, len(ct)): exp += gcd(ct[i], ct[j])
    return 2**exp

def fix_arcs(ct):
    f = ct.count(1)
    t = ct.count(2)
    return comb(f, 2) + t

def fix_anti_tournaments_v2(ct, n):
    """CORRECTED: # tournaments T with σ(T) = T^op for σ with cycle type ct."""
    f = ct.count(1)

    # Rule (a): ≥ 2 fixed points → impossible (both-fixed arc exists)
    if f >= 2:
        return 0

    # Build representative permutation
    perm = [0] * n
    pos = 0
    for c in ct:
        for i in range(c):
            perm[pos + i] = pos + (i + 1) % c
        pos += c

    # Find orbits of σ on unordered arcs {i, j}
    visited = set()
    n_valid_orbits = 0

    for i in range(n):
        for j in range(i + 1, n):
            if (i, j) in visited:
                continue

            # Trace orbit
            orbit = []
            a, b = i, j
            while True:
                key = (min(a, b), max(a, b))
                if key in visited:
                    break
                visited.add(key)
                orbit.append(key)
                a, b = perm[a], perm[b]
                a, b = min(a, b), max(a, b)

            L = len(orbit)

            if L == 1:
                # Fixed arc: check if swap (valid) or both-fixed (invalid)
                a, b = orbit[0]
                if perm[a] == b and perm[b] == a:
                    # Swap → valid, 1 free choice
                    n_valid_orbits += 1
                elif perm[a] == a and perm[b] == b:
                    # Both fixed → impossible
                    return 0
                else:
                    # This shouldn't happen for a fixed arc
                    return 0
            elif L % 2 == 0:
                # Even orbit → valid, 1 free choice
                n_valid_orbits += 1
            else:
                # Odd orbit > 1 → impossible
                return 0

    return 2 ** n_valid_orbits

# ============================================================================
# VERIFICATION vs brute force
# ============================================================================

def fix_anti_brute(sigma, n):
    """Brute force: count T with sigma(T) = T^op."""
    m = comb(n, 2)
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    count = 0
    for bits in range(1 << m):
        adj = [[0] * n for _ in range(n)]
        for k, (i, j) in enumerate(pairs):
            if bits & (1 << k): adj[i][j] = 1
            else: adj[j][i] = 1
        ok = True
        for i in range(n):
            for j in range(n):
                if i == j: continue
                if adj[sigma[i]][sigma[j]] != 1 - adj[i][j]:
                    ok = False
                    break
            if not ok: break
        if ok: count += 1
    return count

print("\n  VERIFICATION: corrected formula vs brute force")
print("-" * 60)

all_match = True
for n in range(3, 7):
    print(f"\n  n={n}:")
    for ct in gen_partitions(n):
        ct_list = list(ct)
        perm = [0] * n
        pos = 0
        for c in ct_list:
            for i in range(c): perm[pos + i] = pos + (i + 1) % c
            pos += c

        formula_val = fix_anti_tournaments_v2(ct_list, n)

        if n <= 6:
            brute_val = fix_anti_brute(perm, n)
            match = "✓" if formula_val == brute_val else "✗"
            if formula_val != brute_val:
                all_match = False
            if formula_val > 0 or brute_val > 0:
                print(f"    {str(ct_list):>20s}: formula={formula_val:6d}, brute={brute_val:6d} {match}")
            else:
                print(f"    {str(ct_list):>20s}: 0 {match}")

print(f"\n  All match: {'YES' if all_match else 'NO — STILL BUGGY'}")

# ============================================================================
# MAIN COMPUTATION with corrected formula
# ============================================================================

print("\n" + "=" * 80)
print("  MAIN COMPUTATION: CORRECTED T_n, T_anti, T_merged")
print("=" * 80)

# Store results for analysis
results = {}

for n in range(3, 14):
    m = comb(n, 2)
    nfact = factorial(n)

    V_sum = 0; T_sum = 0; SC_sum = 0; T_anti_sum = 0

    for ct in gen_partitions(n):
        ct_list = list(ct)
        ctc = cycle_type_count(n, ct_list)
        ft = fix_tournaments(ct_list)
        fa = fix_arcs(ct_list)
        fta = fix_anti_tournaments_v2(ct_list, n)

        V_sum += ctc * ft
        T_sum += ctc * ft * fa
        SC_sum += ctc * fta
        T_anti_sum += ctc * fta * fa

    V_n = V_sum // nfact
    T_n = T_sum // nfact
    SC_n = SC_sum // nfact
    T_anti = T_anti_sum // nfact
    Vm = (V_n + SC_n) // 2
    Tm = (T_n + T_anti) // 2

    results[n] = {
        'V': V_n, 'T': T_n, 'SC': SC_n, 'T_anti': T_anti,
        'Vm': Vm, 'Tm': Tm, 'm': m
    }

    print(f"  n={n:2d}: V={V_n:>12d}, SC={SC_n:>8d}, V_merged={Vm:>10d}")
    print(f"         T={T_n:>12d}, T_anti={T_anti:>8d}, T_merged={Tm:>10d}")
    print(f"         m={m:3d}, T_m/V_m={Tm/Vm:.4f}")

# ============================================================================
# COMPARE WITH KNOWN DATA
# ============================================================================

print("\n" + "=" * 80)
print("  COMPARISON WITH KNOWN EDGE COUNTS")
print("=" * 80)

E_merged = {3: 1, 4: 3, 5: 21, 6: 143}
E_orig = {3: 1, 4: 5, 5: 30, 6: 290, 7: 4086}

print(f"\n  {'n':>3} {'V':>8} {'SC':>5} {'Vm':>8} {'Tm':>10} {'Em':>6} {'Tm/Em':>8} "
      f"{'Eo':>6} {'To/Eo':>8}")

for n in range(3, 8):
    r = results[n]
    em = E_merged.get(n, '?')
    eo = E_orig.get(n, '?')
    if isinstance(em, int):
        tm_em = f"{r['Tm']/em:.4f}"
        tm_2em = r['Tm'] - 2*em
    else:
        tm_em = '?'
        tm_2em = '?'
    if isinstance(eo, int):
        to_eo = f"{r['T']/eo:.4f}"
    else:
        to_eo = '?'
    print(f"  {n:3d} {r['V']:8d} {r['SC']:5d} {r['Vm']:8d} {r['Tm']:10d} "
          f"{str(em):>6s} {tm_em:>8s} {str(eo):>6s} {to_eo:>8s}")

# ============================================================================
# PREDICT MERGED EDGE COUNTS
# ============================================================================

print("\n" + "=" * 80)
print("  PREDICTIONS FOR E_merged AND E_original")
print("=" * 80)

# T_merged/E_merged ratios from known data
if all(n in E_merged for n in [3,4,5,6]):
    ratios = [results[n]['Tm'] / E_merged[n] for n in [3,4,5,6]]
    print(f"  T_m/E_m: {[f'{r:.4f}' for r in ratios]}")
    print(f"  T_m-2E_m: {[results[n]['Tm']-2*E_merged[n] for n in [3,4,5,6]]}")

# For n=7
r7 = results[7]
print(f"\n  n=7 predictions:")
print(f"    V_merged = {r7['Vm']}, T_merged = {r7['Tm']}, SC = {r7['SC']}")
for c in [2.0, 2.1, 2.2, 2.3, 2.5]:
    em_pred = r7['Tm'] / c
    print(f"    T_m/E_m={c:.1f}: E_merged ≈ {em_pred:.0f}")

# ============================================================================
# ALL KEY SEQUENCES
# ============================================================================

print("\n" + "=" * 80)
print("  ALL KEY SEQUENCES (n=3..13)")
print("=" * 80)

seq_names = ['V', 'SC', 'Vm', 'T', 'T_anti', 'Tm']
for name in seq_names:
    vals = [results[n][name] for n in range(3, 14)]
    print(f"  {name:>8s}: {vals}")

# Deficit in merged
print(f"  {'D_m':>8s}: {[results[n]['Vm']*comb(n,2) - results[n]['Tm'] for n in range(3, 14)]}")

# T_anti / SC
print(f"\n  T_anti/SC (should be like 'avg FA over anti-auts'):")
for n in range(3, 14):
    r = results[n]
    if r['SC'] > 0:
        print(f"    n={n}: T_anti/SC = {r['T_anti']/r['SC']:.4f}")

print(f"\n{'='*80}")
print(f"  DONE")
print(f"{'='*80}")
