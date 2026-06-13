#!/usr/bin/env python3
"""
anti_defect_burnside_s20fk.py — Burnside formula for defect = M-1 (anti-D)
kind-pasteur-2026-03-24-S20fk

D(n) counts (T, sigma) with defect = 1 (almost-automorphism).
anti_D(n) counts (T, sigma) with defect = M-1 (almost-anti-automorphism).

The defect spectrum of (T, sigma) depends on the cycle type of sigma
and the arc orbit structure. For defect = M-1: all but 1 arc orbit
must be FULLY INCONSISTENT (all arcs disagree), and exactly 1 orbit
must have exactly 1 agreement.

This is the COMPLEMENT of the defect-1 condition.

For all-odd-cycle sigma: defect is always 0 (Fix(sigma) > 0) or
comes from arc orbits. All orbits have even reversals => defect even.
Max defect = M (all arcs disagree) iff the tournament has an
anti-automorphism (sigma maps T to T^c).

For sigma with cycle type having EXACTLY ONE even cycle 2k:
The odd-reversal orbit of size k forces at least 1 arc to agree (min 1).
So max defect = M - 1. This is achievable iff the even-reversal orbits
are ALL fully inconsistent and the odd orbit has exactly 1 agreement.

This is EXACTLY the defect-1 condition on the COMPLEMENT orbits!

CONJECTURE: anti_D(n) = D(n). The defect-1 and defect-(M-1) counts
are the same by a complement duality.

Let's verify this.
"""

import sys
from math import factorial, gcd, comb
from collections import Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  ANTI-DEFECT BURNSIDE: defect = M-1")
print("  kind-pasteur-2026-03-24-S20fk")
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

# D(n): defect-1 Burnside (our proved formula)
def compute_D(n):
    total = 0
    for ct in partitions(n):
        even_cycles = [c for c in ct if c % 2 == 0]
        if len(even_cycles) != 1: continue
        k = even_cycles[0] // 2
        nperms = count_perms(ct, n)
        a = arc_orbits_count(ct)
        total += nperms * k * (2 ** a)
    return total // factorial(n)

# anti_D(n): defect = M-1 Burnside
# By complement duality: if we replace T with T^c in the defect formula,
# defect(T^c, sigma) = #{arcs where sigma(T^c) != T^c}
#                    = #{arcs where 1 - sigma(T) != 1 - T}
#                    = #{arcs where sigma(T) != T}
#                    = defect(T, sigma)
#
# Wait: that says defect is COMPLEMENT-INVARIANT! defect(T, sigma) = defect(T^c, sigma).
# So: #{(T, sigma) : defect = d} = #{(T^c, sigma) : defect = d} = #{(T, sigma) : defect = d}.
# This is trivially true (complement is a bijection on T).
#
# But anti_D counts (T, sigma) with defect = M-1.
# And defect(T, sigma) + defect(T^c, sigma) is NOT generally M.
# Actually: defect(T^c, sigma) = defect(T, sigma) (as shown above).
# So the defect is complement-invariant.
#
# What about: defect(T, sigma) vs defect(T, sigma composed with complement)?
# sigma_comp = complement ∘ sigma: this maps T to complement(sigma(T)).
# defect(T, sigma_comp) = #{arcs where sigma_comp(T) != T}
#                       = #{arcs where complement(sigma(T)) != T}
#                       = #{arcs where sigma(T) = T} (complement flips)
#                       = M - defect(T, sigma)!
#
# YES! defect(T, sigma ∘ comp) = M - defect(T, sigma).
#
# So: the map sigma -> sigma ∘ comp bijection between permutations that
# maps defect d to defect M-d.
#
# THEREFORE: #{(T, sigma) : defect = d} = #{(T, sigma ∘ comp) : defect = M-d}
#          = #{(T, sigma') : defect = M-d} where sigma' = sigma ∘ comp ranges over...
# BUT: "comp" here is the complement operation on tournaments, which is NOT a permutation!
# It's the map A[i][j] -> 1-A[i][j]. This doesn't correspond to a sigma in S_n.
#
# Let me reconsider. The "defect of (T, sigma)" involves sigma in S_n acting on T.
# The complement is NOT in S_n. So the bijection sigma -> sigma ∘ comp doesn't stay in S_n.
#
# Hmm. Let me just compute anti_D by brute force for small n and compare with D.

print("\n  COMPARING D(n) and anti_D(n) (brute force):")
print(f"  {'n':>3} {'D':>8} {'anti_D':>8} {'Equal?':>8} {'M':>4}")

for n in range(3, 8):
    M = comb(n, 2)
    D = compute_D(n)

    # Brute force anti_D at small n
    if n <= 6:
        from itertools import permutations as perms
        tiles = []
        for y in range(1, n-1):
            for x in range(n, y+1, -1):
                tiles.append((x, y))
        m = len(tiles)
        VERTS = list(range(n, 0, -1))
        all_perms = list(perms(range(n)))

        tv = [(VERTS.index(x), VERTS.index(y)) for x, y in tiles]

        def b2a(bits):
            A = [[0]*n for _ in range(n)]
            for k in range(n-1): A[k][k+1] = 1
            for i in range(m):
                xi, yi = tv[i]
                if bits[i] == 0: A[xi][yi] = 1
                else: A[yi][xi] = 1
            return A

        anti_D_raw = 0
        for mask in range(1 << m):
            bits = [(mask >> k) & 1 for k in range(m)]
            A = b2a(bits)
            for sigma in all_perms:
                si = [0]*n
                for i in range(n): si[sigma[i]] = i
                d = 0
                for i in range(n):
                    for j in range(i+1, n):
                        if A[si[i]][si[j]] != A[i][j]:
                            d += 1
                if d == M - 1:
                    anti_D_raw += 1

        anti_D = anti_D_raw // factorial(n)
        print(f"  {n:3d} {D:8d} {anti_D:8d} {'YES' if D == anti_D else 'NO':>8} {M:4d}")
    else:
        print(f"  {n:3d} {D:8d} {'?':>8} {'?':>8} {M:4d}")

print("""
If D = anti_D: the defect spectrum is SYMMETRIC around M/2.
This would mean: the "almost-automorphism" structure perfectly mirrors
the "almost-anti-automorphism" structure. The D(n) formula would also
give the AAA count.

If D != anti_D: the spectrum is asymmetric, and the AAA count has a
different Burnside formula.
""")

print("DONE.")
print("=" * 80)
