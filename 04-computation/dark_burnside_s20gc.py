#!/usr/bin/env python3
"""
dark_burnside_s20gc.py — Burnside formula for dark tournaments to large n
kind-pasteur-2026-03-25-S20gc

D(n) = (1/n!) * sum_{sigma with at least one even cycle} 2^{pair_orbits(sigma)}

This is the COMPLEMENT of the tournament Burnside:
  V(n) = (1/n!) * sum_{all-odd sigma} 2^{a(sigma)}    [tournaments]
  D(n) = (1/n!) * sum_{even-cycle sigma} 2^{a(sigma)}  [dark tournaments]
  G(n) = V(n) + D(n) = (1/n!) * sum_{all sigma} 2^{a(sigma)}  [all graphs]

So D(n) = G(n) - V(n) = A000088(n) - A000568(n).

But we can also compute D(n) DIRECTLY via Burnside over even-cycle types.
This gives D(n) at all n without needing A000088 or A000568 separately.
"""

import sys
from math import factorial, gcd
from collections import Counter

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  DARK TOURNAMENT BURNSIDE FORMULA")
print("  kind-pasteur-2026-03-25-S20gc")
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

def pair_orbits(ct):
    """Number of pair orbits (for undirected graphs) = arc orbits for tournaments."""
    cycles = list(ct)
    total = sum(c // 2 for c in cycles)
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            total += gcd(cycles[i], cycles[j])
    return total

# Compute V(n), D(n), G(n) for each n
print(f"\n  {'n':>3} {'V(tour)':>10} {'D(dark)':>10} {'G(all)':>10} {'D/V':>8} {'D/G':>8}")
print(f"  {'-'*52}")

for n in range(1, 21):
    V_raw = 0
    D_raw = 0
    G_raw = 0

    for ct in partitions(n):
        nperms = count_perms(ct, n)
        a = pair_orbits(ct)
        contrib = nperms * (2 ** a)
        G_raw += contrib

        has_even = any(c % 2 == 0 for c in ct)
        if has_even:
            D_raw += contrib
        else:
            V_raw += contrib

    V = V_raw // factorial(n)
    D = D_raw // factorial(n)
    G = G_raw // factorial(n)

    assert G == V + D, f"G != V + D at n={n}"

    ratio_DV = D / V if V > 0 else 0
    ratio_DG = D / G if G > 0 else 0

    print(f"  {n:3d} {V:10d} {D:10d} {G:10d} {ratio_DV:8.4f} {ratio_DG:8.4f}")

print(f"""
ANALYSIS:

D/V ratio (dark/light):
  Peaks at n=5 (1.83) then DECREASES.
  At n=8: D/V = 0.79 (dark < light!)
  At n=10: D/V << 1 (dark becomes rare)
  At large n: D/V -> 0 (almost all graphs are even)

D/G ratio (dark/total):
  Starts at 0.50 (n=2,3) then decreases.
  At large n: D/G -> 0.

The dark tournament world SHRINKS relative to light as n grows.
This is because at large n, almost all graphs have trivial
automorphism group (|Aut|=1), and for trivial Aut, the sign
representation is always trivial (even). So generic graphs are even.

The DARK SEQUENCE (A000088 - A000568):
  0, 1, 2, 7, 22, 100, 588, 5466, ...
This is Burnside-computable to any n via the even-cycle sum.
""")

print("DONE.")
print("=" * 80)
