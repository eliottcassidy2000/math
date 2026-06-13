#!/usr/bin/env python3
"""
gn_merged_burnside_v3_s214.py — CORRECT anti-aut formula via directed arc orbits
opus-2026-03-22-S214

The CORRECT analysis uses DIRECTED arc orbits:

1. σ acts on directed arcs (i,j) → (σi, σj).
2. Orbits come in conjugate pairs: orbit of (i,j) and orbit of (j,i).
3. Self-conjugate orbit: (j,i) appears IN the same orbit as (i,j), at position k.
4. Non-self-conjugate: (j,i) is in a separate orbit.

Rules:
- Self-conjugate orbit with k ODD: VALID (1 free choice)
- Self-conjugate orbit with k EVEN: INVALID → Fix_anti = 0
- Non-self-conjugate pair, even length: VALID (1 free choice per pair)
- Non-self-conjugate pair, odd length: INVALID → Fix_anti = 0
- Fixed directed arc (σi=i, σj=j): length 1, k=0 even → INVALID (or: f≥2 → 0)

Fix_anti = 2^(valid_count) if no invalid orbits, else 0.
"""

import sys
from math import comb, factorial, gcd
from collections import Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  CORRECT ANTI-AUTOMORPHISM FORMULA (DIRECTED ARC ORBITS)")
print("  opus-2026-03-22-S214")
print("=" * 80)

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

def fix_anti_correct(ct, n):
    """CORRECT Fix_anti via directed arc orbit analysis.
    FIX: non-self-conjugate pairs count as 1 free choice per PAIR, not per orbit."""
    f = ct.count(1)
    if f >= 2:
        return 0

    # Build permutation
    perm = [0] * n
    pos = 0
    for c in ct:
        for i in range(c):
            perm[pos + i] = pos + (i + 1) % c
        pos += c

    # Find orbits of σ on DIRECTED arcs (i,j), i≠j
    visited = set()
    free_choices = 0

    for i in range(n):
        for j in range(n):
            if i == j: continue
            if (i, j) in visited: continue

            # Trace orbit of directed arc (i,j)
            orbit = []
            a, b = i, j
            while (a, b) not in visited:
                visited.add((a, b))
                orbit.append((a, b))
                a, b = perm[a], perm[b]

            L = len(orbit)
            orbit_set = set(orbit)

            # Check if self-conjugate: does (j,i) appear in this orbit?
            reverse = (j, i)
            if reverse in orbit_set:
                # Self-conjugate. Find position k of (j,i).
                k = orbit.index(reverse)
                if k % 2 == 1:
                    # k odd → VALID, 1 free choice
                    free_choices += 1
                else:
                    # k even → INVALID
                    return 0
            else:
                # Non-self-conjugate. Paired with orbit of (j,i).
                # Check if the CONJUGATE orbit has already been visited.
                # If (j,i) is already visited, we already counted this pair.
                if reverse in visited:
                    # Already counted when we processed the conjugate orbit
                    pass
                else:
                    # First time seeing this pair
                    if L % 2 == 0:
                        # Even length → VALID (1 choice for the pair)
                        free_choices += 1
                    else:
                        # Odd length → INVALID
                        return 0

    return 2 ** free_choices

# ============================================================================
# VERIFICATION
# ============================================================================

def fix_anti_brute(sigma, n):
    m = comb(n, 2)
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    count = 0
    for bits in range(1 << m):
        adj = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(pairs):
            if bits & (1 << k): adj[i][j] = 1
            else: adj[j][i] = 1
        ok = True
        for i in range(n):
            for j in range(n):
                if i == j: continue
                if adj[sigma[i]][sigma[j]] != 1 - adj[i][j]:
                    ok = False; break
            if not ok: break
        if ok: count += 1
    return count

print("\n  VERIFICATION (all cycle types, n=3..6)")
print("-" * 60)

all_ok = True
for n in range(3, 7):
    print(f"\n  n={n}:")
    for ct in gen_partitions(n):
        ct_list = list(ct)
        perm = [0] * n
        pos = 0
        for c in ct_list:
            for i in range(c): perm[pos + i] = pos + (i + 1) % c
            pos += c

        formula = fix_anti_correct(ct_list, n)
        brute = fix_anti_brute(perm, n) if n <= 6 else None

        if brute is not None:
            match = "✓" if formula == brute else "✗"
            if formula != brute: all_ok = False
            if formula > 0 or brute > 0:
                print(f"    {str(ct_list):>24s}: formula={formula:6d}, brute={brute:6d} {match}")
        else:
            if formula > 0:
                print(f"    {str(ct_list):>24s}: formula={formula:6d}")

print(f"\n  ALL MATCH: {'YES ✓' if all_ok else 'NO ✗'}")

# ============================================================================
# MAIN COMPUTATION
# ============================================================================

print("\n" + "=" * 80)
print("  MAIN COMPUTATION: CORRECT T_n, SC_n, T_anti, T_merged")
print("=" * 80)

E_merged_known = {3: 1, 4: 3, 5: 21, 6: 143}
E_orig_known = {3: 1, 4: 5, 5: 30, 6: 290, 7: 4086}

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
        fta = fix_anti_correct(ct_list, n)

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

    results[n] = {'V': V_n, 'T': T_n, 'SC': SC_n, 'T_anti': T_anti, 'Vm': Vm, 'Tm': Tm, 'm': m}

    em_str = f", E_m={E_merged_known[n]}" if n in E_merged_known else ""
    eo_str = f", E_o={E_orig_known[n]}" if n in E_orig_known else ""
    print(f"  n={n:2d}: V={V_n:>12d}, SC={SC_n:>8d}, Vm={Vm:>10d}, "
          f"T={T_n:>12d}, T_anti={T_anti:>8d}, Tm={Tm:>10d}{em_str}{eo_str}")

# ============================================================================
# ANALYSIS
# ============================================================================

print("\n" + "=" * 80)
print("  KEY RATIOS AND PREDICTIONS")
print("=" * 80)

print(f"\n  {'n':>3} {'Vm':>10} {'Tm':>12} {'Em_known':>8} {'Tm/Em':>8} "
      f"{'Tm-2Em':>8} {'Eo_known':>8} {'To/Eo':>8}")

for n in range(3, 8):
    r = results[n]
    em = E_merged_known.get(n)
    eo = E_orig_known.get(n)
    tm_em = f"{r['Tm']/em:.4f}" if em else "?"
    tm_2em = f"{r['Tm']-2*em}" if em else "?"
    to_eo = f"{r['T']/eo:.4f}" if eo else "?"
    print(f"  {n:3d} {r['Vm']:10d} {r['Tm']:12d} {str(em or '?'):>8s} {tm_em:>8s} "
          f"{str(tm_2em):>8s} {str(eo or '?'):>8s} {to_eo:>8s}")

# Predict E_merged(7)
print(f"\n  PREDICTIONS for n=7:")
r7 = results[7]
# From ratios: Tm/Em ≈ 2.0-2.6 (range from n=3..6)
for c in [2.0, 2.1, 2.2, 2.3, 2.5]:
    em_pred = r7['Tm'] / c
    print(f"    Tm/Em={c:.1f}: E_merged(7) ≈ {em_pred:.0f}")

# ============================================================================
# ALL SEQUENCES
# ============================================================================

print("\n" + "=" * 80)
print("  ALL NEW SEQUENCES")
print("=" * 80)

for name in ['V', 'SC', 'Vm', 'T', 'T_anti', 'Tm']:
    vals = [results[n][name] for n in range(3, 14)]
    print(f"  {name:>8s}: {vals}")

# SC sequence should be known (self-complementary tournament iso classes)
print(f"\n  SC_n verification:")
print(f"    Known SC classes: n=3:2, n=4:2, n=5:8, n=6:12")
print(f"    Our SC_n: {[results[n]['SC'] for n in [3,4,5,6]]}")

# Vm verification
print(f"\n  V_merged verification:")
print(f"    Known: n=3:2, n=4:3, n=5:10, n=6:34")
print(f"    Formula (V+SC)/2: {[(results[n]['V']+results[n]['SC'])//2 for n in [3,4,5,6]]}")

print(f"\n{'='*80}")
print(f"  DONE")
print(f"{'='*80}")
