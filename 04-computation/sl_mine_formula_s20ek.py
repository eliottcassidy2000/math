#!/usr/bin/env python3
"""
sl_mine_formula_s20ek.py — Efficient D(n) formula (Burnside defect-1 bound on SL_mine)
kind-pasteur-2026-03-23-S20ek

PROVED formula:
  D(n) = (1/n!) * sum over cycle types ct with EXACTLY ONE even cycle 2k of:
         count(ct) * k * 2^{a(ct)}

where:
  a(ct) = sum(c//2 for c in ct) + sum(gcd(ci,cj) for i<j)
  k = half the length of the unique even cycle

D(n) is the Burnside defect-1 count:
  D(n) = sum_C neutral_labeled(C)
       = sum_C (sum over self-loop arc orbits O of |O|)

SL_mine ≤ D(n) with equality iff all self-loop arc orbits have size 1.
Verified: D = SL_mine at n=3,4,5. D = SL_mine + 2 at n=6.

The correction D - SL_mine comes from iso classes with |Aut|>1 that have
self-loop orbits of size > 1. This correction is small at all n.

Also computes: T_n, E_formula = (T - D)/2 (what E would be without multi-edges),
and the actual multi-edge-corrected values.
"""

import sys
from math import factorial, gcd
from collections import Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  EFFICIENT D(n) / SL_mine FORMULA")
print("  kind-pasteur-2026-03-23-S20ek")
print("=" * 80)

def partitions(n, max_val=None):
    if max_val is None: max_val = n
    if n == 0: yield (); return
    for i in range(min(n, max_val), 0, -1):
        for rest in partitions(n - i, i):
            yield (i,) + rest

def count_perms(ct, n):
    mult = Counter(ct)
    d = 1
    for c, m in mult.items(): d *= (c**m) * factorial(m)
    return factorial(n) // d

def arc_orbits_count(ct):
    """Number of arc orbits for cycle type ct."""
    cycles = list(ct)
    total = sum(c // 2 for c in cycles)
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            total += gcd(cycles[i], cycles[j])
    return total

def fixed_arcs_count(ct):
    """Number of arcs fixed by sigma with cycle type ct.
    An arc {u,v} is fixed iff both u,v are fixed points (1-cycles)."""
    f = sum(1 for c in ct if c == 1)
    return f * (f - 1) // 2

# Compute T_n via standard Burnside (on (tournament, arc) pairs)
def compute_T(n):
    """T_n = sum_C #{arc orbits of C} via Burnside."""
    total = 0
    for ct in partitions(n):
        if any(c % 2 == 0 for c in ct): continue  # even cycles => Fix(T)=0
        nperms = count_perms(ct, n)
        a = arc_orbits_count(ct)
        fa = fixed_arcs_count(ct)
        total += nperms * (2 ** a) * fa
    return total // factorial(n)

# Compute D(n) via defect-1 Burnside
def compute_D(n):
    """D(n) = sum_C neutral_labeled(C) via Burnside defect formula."""
    total = 0
    for ct in partitions(n):
        # Count even cycles
        even_cycles = [c for c in ct if c % 2 == 0]
        if len(even_cycles) != 1: continue

        k = even_cycles[0] // 2  # half the even cycle length
        nperms = count_perms(ct, n)
        a = arc_orbits_count(ct)
        total += nperms * k * (2 ** a)
    return total // factorial(n)

# Compute twin_SL (only 2-cycles contribute, k=1)
def compute_twin_SL(n):
    """twin_SL = defect formula restricted to 2-cycles only."""
    total = 0
    for ct in partitions(n):
        even_cycles = [c for c in ct if c % 2 == 0]
        if len(even_cycles) != 1: continue
        if even_cycles[0] != 2: continue  # only 2-cycles
        nperms = count_perms(ct, n)
        a = arc_orbits_count(ct)
        total += nperms * (2 ** a)  # k=1
    return total // factorial(n)

# Compute V (number of iso classes) via standard Burnside
def compute_V(n):
    total = 0
    for ct in partitions(n):
        if any(c % 2 == 0 for c in ct): continue
        total += count_perms(ct, n) * (2 ** arc_orbits_count(ct))
    return total // factorial(n)

print()
print(f"{'n':>3} {'V':>8} {'T':>10} {'D':>8} {'twin_SL':>8} {'D-twin':>8} {'V-1':>8} {'D/V':>8}")
print("-" * 80)

# Known exact SL_mine values (from direct computation)
SL_exact = {3: 2, 4: 6, 5: 16, 6: 58}
E_known = {3:1, 4:5, 5:30, 6:290, 7:4086, 8:91161, 9:3380751}

results = {}

for n in range(3, 25):
    t0 = time.time()
    V = compute_V(n)
    T = compute_T(n)
    D = compute_D(n)
    twin = compute_twin_SL(n)
    elapsed = time.time() - t0

    results[n] = {'V': V, 'T': T, 'D': D, 'twin': twin}

    print(f"{n:3d} {V:8d} {T:10d} {D:8d} {twin:8d} {D-twin:8d} {V-1:8d} {D/V:8.3f}", end="")

    # Known exact values check
    if n in SL_exact:
        print(f"  SL_exact={SL_exact[n]}, gap={D - SL_exact[n]}", end="")
    if n in E_known:
        E_from_D = (T - D) // 2  # E if no multi-edges AND D=SL
        print(f"  E_known={E_known[n]}", end="")

    if elapsed > 0.1:
        print(f"  ({elapsed:.1f}s)", end="")
    print()

    if elapsed > 30:
        break

# Analysis
print()
print("=" * 60)
print("ANALYSIS")
print("=" * 60)

print("\nD(n) = Burnside defect-1 count = sum_C neutral_labeled(C)")
print("SL_mine ≤ D(n), with D = SL_mine when all SL orbits have size 1")
print("The correction D - SL_mine comes from high-|Aut| classes")
print()

# Twin vs D ratio
print("twin_SL contribution to D:")
for n in sorted(results.keys()):
    r = results[n]
    if r['D'] > 0:
        ratio = r['twin'] / r['D']
        print(f"  n={n:2d}: twin={r['twin']:8d}, D={r['D']:8d}, twin/D={ratio:.4f}, non-twin/D={1-ratio:.4f}")

print()
print("D/V ratio (self-loop density per class):")
for n in sorted(results.keys()):
    r = results[n]
    print(f"  n={n:2d}: D/V = {r['D']}/{r['V']} = {r['D']/r['V']:.4f}")

# Multi-edge surplus analysis
print()
print("Multi-edge surplus = non_SL - 2E = (T - D) - 2E:")
for n in sorted(results.keys()):
    r = results[n]
    if n in E_known:
        nonSL_lower = r['T'] - r['D']  # using D as upper bound on SL
        surplus_lower = nonSL_lower - 2*E_known[n]
        print(f"  n={n:2d}: T-D={nonSL_lower:8d}, 2E={2*E_known[n]:8d}, surplus≥{surplus_lower:6d}")

print()
print("DONE.")
print("=" * 80)
