#!/usr/bin/env python3
"""
sl_mine_summary_s20el.py -- Final summary of SL_mine analysis
kind-pasteur-2026-03-23-S20el

KEY RESULTS:

1. PROVED FORMULA for D(n) (defect-1 Burnside count):
   D(n) = (1/n!) * sum_{ct with exactly one even cycle 2k} count(ct) * k * 2^{a(ct)}

   where a(ct) = sum(c//2) + sum(gcd(ci,cj) for i<j)

2. D(n) >= SL_mine(n) with small correction from |Aut|>1 classes

3. T - 2E != SL_mine at n >= 5 (multi-edge surplus in metagraph)

4. D(n) is Burnside-computable to arbitrary n (O(partitions) time)

5. Twin_SL/D -> 1 as n -> infinity (2-cycles dominate)
"""

import sys
from math import factorial, gcd, comb
from collections import Counter

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  SL_mine SUMMARY AND CORRECTED METAGRAPH EQUATIONS")
print("  kind-pasteur-2026-03-23-S20el")
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
    cycles = list(ct)
    total = sum(c // 2 for c in cycles)
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            total += gcd(cycles[i], cycles[j])
    return total

def fixed_arcs(ct):
    f = sum(1 for c in ct if c == 1)
    return f * (f-1) // 2

def compute_all(n):
    # V_n: iso classes
    V_raw = 0
    for ct in partitions(n):
        if any(c % 2 == 0 for c in ct): continue
        V_raw += count_perms(ct, n) * (2 ** arc_orbits_count(ct))
    V = V_raw // factorial(n)

    # T_n: total arc orbit transitions
    T_raw = 0
    for ct in partitions(n):
        if any(c % 2 == 0 for c in ct): continue
        T_raw += count_perms(ct, n) * (2 ** arc_orbits_count(ct)) * fixed_arcs(ct)
    T = T_raw // factorial(n)

    # D_n: defect-1 count (upper bound on SL_mine)
    D_raw = 0
    twin_raw = 0
    for ct in partitions(n):
        even_cycles = [c for c in ct if c % 2 == 0]
        if len(even_cycles) != 1: continue
        k = even_cycles[0] // 2
        nperms = count_perms(ct, n)
        a = arc_orbits_count(ct)
        contrib = nperms * k * (2 ** a)
        D_raw += contrib
        if even_cycles[0] == 2:
            twin_raw += contrib
    D = D_raw // factorial(n)
    twin = twin_raw // factorial(n)

    return V, T, D, twin

# Known exact values
SL_exact = {3: 2, 4: 6, 5: 16, 6: 58}
E_known = {3:1, 4:5, 5:30, 6:290, 7:4086, 8:91161, 9:3380751}

print()
print("CORRECTED METAGRAPH EQUATIONS:")
print("  T = SL_mine + non_SL")
print("  non_SL >= 2*E (multi-edges possible)")
print("  D >= SL_mine (orbit size correction)")
print("  E = (T - SL_mine)/2 ONLY IF no multi-edges")
print()

print(f"{'n':>3} {'V':>8} {'T':>10} {'D':>8} {'SL_mine':>8} {'D-SL':>6} {'non_SL':>8} {'2E':>8} {'surplus':>8} {'D/V':>7} {'tw/D':>6}")
print("-" * 100)

for n in range(3, 16):
    V, T, D, twin = compute_all(n)

    sl = SL_exact.get(n)
    if sl is not None:
        nonsl = T - sl
        gap = D - sl
    else:
        sl = D  # upper bound
        nonsl = T - D  # lower bound on non_SL
        gap = None

    E = E_known.get(n)
    if E:
        two_E = 2*E
        surplus = nonsl - two_E if gap is not None else f">={nonsl - two_E}"
    else:
        two_E = None
        surplus = None

    tw_ratio = twin/D if D > 0 else 0

    sl_str = f"{sl:8d}" if n in SL_exact else f"<={D:7d}"
    gap_str = f"{gap:6d}" if gap is not None else "  ?"
    nonsl_str = f"{nonsl:8d}" if n in SL_exact else f">={nonsl:7d}"
    twoE_str = f"{two_E:8d}" if two_E else "       ?"
    surplus_str = f"{surplus:8d}" if isinstance(surplus, int) else f"  {surplus}" if surplus else "       ?"

    print(f"{n:3d} {V:8d} {T:10d} {D:8d} {sl_str} {gap_str} {nonsl_str} {twoE_str} {surplus_str} {D/V:7.3f} {tw_ratio:6.3f}")

print()
print("LEGEND:")
print("  V = iso classes, T = total arc orbit transitions")
print("  D = Burnside defect-1 count (Burnside-computable upper bound on SL_mine)")
print("  SL_mine = actual self-loop orbit count")
print("  D-SL = orbit size correction (from |Aut|>1 classes)")
print("  non_SL = non-self-loop transitions = T - SL_mine")
print("  2E = twice the number of simple edges")
print("  surplus = non_SL - 2E = multi-edge surplus")
print("  tw/D = fraction of D from 2-cycle types (twin contribution)")

print()
print("=" * 80)
print("KEY FORMULAE (all Burnside-computable):")
print("=" * 80)
print()
print("V(n) = (1/n!) sum_{ct all-odd} count(ct) * 2^{a(ct)}")
print("T(n) = (1/n!) sum_{ct all-odd} count(ct) * 2^{a(ct)} * C(f(ct),2)")
print("D(n) = (1/n!) sum_{ct with 1 even cycle 2k} count(ct) * k * 2^{a(ct)}")
print()
print("where a(ct) = sum(c//2) + sum(gcd(ci,cj)) and f(ct) = #{1-cycles}")
print()
print("RELATIONS:")
print("  SL_mine <= D (equality iff all SL arc orbits have size 1)")
print("  E <= (T - SL_mine)/2 (equality iff no multi-edges)")
print("  E = (T - D + correction)/2 + multi_edge_correction")
print()
print("ASYMPTOTICS:")
print("  D/V -> 0 (most classes have no self-loops)")
print("  twin/D -> 1 (2-cycle types dominate)")
print("  D ~ twin_SL ~ sum over k-cycles of k*2^{a(2,k^*)} terms")
print()
print("DONE.")
print("=" * 80)
