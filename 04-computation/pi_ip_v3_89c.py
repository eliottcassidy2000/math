#!/usr/bin/env python3
"""
pi_ip_v3_89c.py — IP structure via smart enumeration
opus-2026-03-14-S89c

G(P_7) has 80 vertices — can't brute force 2^80.
But we only need IP coefficients up to degree = max packing ≤ floor(7/3) = 2.

Actually, the maximum number of vertex-disjoint odd cycles in P_7
is at most floor(7/3) = 2 (using 3-cycles) or floor(7/5) = 1 (using 5-cycles).
With mixed sizes: 3 + 3 = 6 vertices (max 2 cycles) or 3 + 5 = 8 > 7 (impossible).
Or one 7-cycle (max 1 cycle). So max packing = 2.

Therefore IP has degree ≤ 2:
  IP(G, x) = 1 + c₁x + c₂x²

And c₁ = 80 (number of cycles), c₂ = number of disjoint pairs.
IP(G, 2) = 1 + 80×2 + c₂×4 = 161 + 4c₂ = 189
So c₂ = (189 - 161)/4 = 7.

Let's verify this and extend to P_11.
"""

from itertools import combinations

def paley_tournament(p):
    qr = set()
    for a in range(1, p):
        qr.add((a*a) % p)
    adj = {i: set() for i in range(p)}
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in qr:
                adj[i].add(j)
    return adj, qr

def find_directed_cycles(adj, n, length):
    """Find directed cycles of given length (canonical: start is min vertex)."""
    cycles = []
    for start in range(n):
        stack = [(start, [start], 1 << start)]
        while stack:
            v, path, mask = stack.pop()
            if len(path) == length:
                if start in adj[v]:
                    cycles.append(tuple(path))
                continue
            for u in adj[v]:
                if u < start:
                    continue
                if mask & (1 << u):
                    continue
                stack.append((u, path + [u], mask | (1 << u)))
    return cycles

print("=" * 70)
print("PART 1: P_7 — IP coefficients")
print("=" * 70)

p = 7
adj, qr = paley_tournament(p)

# Enumerate all odd cycles
all_cycles = []
for length in range(3, p+1, 2):
    cycles = find_directed_cycles(adj, p, length)
    print(f"  t_{length}(P_7) = {len(cycles)}")
    all_cycles.extend(cycles)

c1 = len(all_cycles)
print(f"\n  Total odd cycles (IP_1): c₁ = {c1}")

# Max packing: can fit at most 2 vertex-disjoint odd cycles
# (since smallest cycle is length 3, and 3+3 = 6 ≤ 7, but 3+3+3 = 9 > 7)
# Also 3+5 = 8 > 7 and 5+5 = 10 > 7
# So max packing uses two 3-cycles.

# Count disjoint pairs
c2 = 0
disjoint_pairs = []
for i in range(len(all_cycles)):
    si = set(all_cycles[i])
    for j in range(i+1, len(all_cycles)):
        if not si & set(all_cycles[j]):
            c2 += 1
            disjoint_pairs.append((i, j))

print(f"  Disjoint pairs (IP_2): c₂ = {c2}")

# Verify IP(G, 2) = H(P_7)
ip_at_2 = 1 + c1 * 2 + c2 * 4
print(f"\n  IP(G, 2) = 1 + {c1}×2 + {c2}×4 = {ip_at_2}")
print(f"  H(P_7) = 189")
print(f"  Match: {'✓' if ip_at_2 == 189 else '✗'}")

# Full IP polynomial
print(f"\n  IP(G(P_7), x) = 1 + {c1}x + {c2}x²")

# Evaluations
for x in [-1, 0, 1, 2, 3]:
    val = 1 + c1*x + c2*x*x
    print(f"  IP(G, {x:2d}) = {val}")

# What are the 7 disjoint pairs? They must be pairs of 3-cycles.
print(f"\n  The {c2} disjoint pairs:")
for idx, (i, j) in enumerate(disjoint_pairs):
    ci = all_cycles[i]
    cj = all_cycles[j]
    print(f"    {set(ci)} ∥ {set(cj)} (lengths {len(ci)},{len(cj)})")

# Orbits of disjoint pairs under Z/7Z
print(f"\n  These {c2} pairs form {c2}//7 = {c2//7} orbit(s) under Z/7Z")

print()
print("=" * 70)
print("PART 2: P_11 — IP coefficients")
print("=" * 70)

p = 11
adj, qr = paley_tournament(p)

all_cycles_11 = []
for length in range(3, p+1, 2):
    print(f"  Finding {length}-cycles...", end=" ", flush=True)
    cycles = find_directed_cycles(adj, p, length)
    print(f"t_{length}(P_11) = {len(cycles)}")
    all_cycles_11.extend(cycles)

c1_11 = len(all_cycles_11)
print(f"\n  Total odd cycles (IP_1): c₁ = {c1_11}")

# Max packing for P_11:
# 3+3+3 = 9 ≤ 11: up to 3 disjoint 3-cycles (but need all to exist)
# 3+3+5 = 11: one 3-cycle pair + one 5-cycle (exact cover!)
# 5+5 = 10 ≤ 11: two disjoint 5-cycles
# 3+3+3+... nope, 3+3+3+3 = 12 > 11
# 7+3 = 10 ≤ 11
# 9+3 = 12 > 11
# So max packing ≤ 3 (three 3-cycles using 9 vertices)
# Or max packing could be 3 (3+3+5 = 11, which uses ALL vertices)

# Count disjoint pairs (IP_2)
print(f"\n  Counting disjoint pairs...")
# With ~16000+ cycles, we need to be efficient
# Use vertex sets

cycle_vsets = [frozenset(c) for c in all_cycles_11]

# Group by length
by_length_11 = {}
for i, c in enumerate(all_cycles_11):
    l = len(c)
    if l not in by_length_11:
        by_length_11[l] = []
    by_length_11[l].append(i)

print(f"  Cycles by length: {[(l, len(v)) for l, v in sorted(by_length_11.items())]}")

# Count disjoint pairs by type
c2_11 = 0
pair_types = {}

for l1 in sorted(by_length_11.keys()):
    for l2 in sorted(by_length_11.keys()):
        if l2 < l1:
            continue
        count = 0
        indices1 = by_length_11[l1]
        indices2 = by_length_11[l2] if l2 != l1 else indices1

        for ii, i in enumerate(indices1):
            si = cycle_vsets[i]
            start_j = ii + 1 if l2 == l1 else 0
            for jj in range(start_j, len(indices2)):
                j = indices2[jj]
                if not si & cycle_vsets[j]:
                    count += 1

        if count > 0:
            pair_types[(l1, l2)] = count
            c2_11 += count
            print(f"    ({l1},{l2})-disjoint pairs: {count}")

print(f"\n  Total disjoint pairs (IP_2): c₂ = {c2_11}")

# Count disjoint triples (IP_3)
print(f"\n  Counting disjoint triples...")
# For triples, the most likely type is (3,3,3) or (3,3,5)
# since 3+3+3=9 ≤ 11 and 3+3+5=11

c3_11 = 0

# (3,3,3) triples: three pairwise disjoint 3-cycles
three_idx = by_length_11.get(3, [])
print(f"  Checking (3,3,3) triples among {len(three_idx)} 3-cycles...")
for ii in range(len(three_idx)):
    i = three_idx[ii]
    si = cycle_vsets[i]
    for jj in range(ii+1, len(three_idx)):
        j = three_idx[jj]
        sj = cycle_vsets[j]
        if si & sj:
            continue
        sij = si | sj
        for kk in range(jj+1, len(three_idx)):
            k = three_idx[kk]
            if not sij & cycle_vsets[k]:
                c3_11 += 1

print(f"    (3,3,3) triples: {c3_11}")

# (3,3,5) triples
five_idx = by_length_11.get(5, [])
c3_335 = 0
for ii in range(len(three_idx)):
    i = three_idx[ii]
    si = cycle_vsets[i]
    for jj in range(ii+1, len(three_idx)):
        j = three_idx[jj]
        sj = cycle_vsets[j]
        if si & sj:
            continue
        sij = si | sj  # 6 vertices used
        for kk in range(len(five_idx)):
            k = five_idx[kk]
            if not sij & cycle_vsets[k]:
                c3_335 += 1

print(f"    (3,3,5) triples: {c3_335}")

# (3,5,5): 3+5+5 = 13 > 11, impossible
# (5,5,5): 15 > 11, impossible
# (3,3,7): 3+3+7 = 13 > 11, impossible

c3_total = c3_11 + c3_335
print(f"    Total triples (IP_3): c₃ = {c3_total}")

# Any quadruples? 3+3+3+3 = 12 > 11, impossible.
# So IP has degree ≤ 3.

# IP(G(P_11), 2) = 1 + c1×2 + c2×4 + c3×8
ip_at_2_11 = 1 + c1_11 * 2 + c2_11 * 4 + c3_total * 8
print(f"\n  IP(G(P_11), 2) = 1 + {c1_11}×2 + {c2_11}×4 + {c3_total}×8")
print(f"  = {ip_at_2_11}")
print(f"  H(P_11) = 95095")
print(f"  Match: {'✓' if ip_at_2_11 == 95095 else '✗'}")

if ip_at_2_11 != 95095:
    print(f"  Difference: {95095 - ip_at_2_11}")
    print(f"  Missing: need {(95095 - ip_at_2_11) // 8} more triples or higher")

print()
print("=" * 70)
print("PART 3: IP polynomial structure")
print("=" * 70)

# For P_7:
print(f"\n  P_7: IP(G, x) = 1 + 80x + 7x²")
print(f"  Zeros: x = (-80 ± √(6400-28))/14 = (-80 ± √6372)/14")
import math
disc = 80**2 - 4*7
print(f"  Discriminant = {disc} = {math.isqrt(disc) if disc >= 0 else 'complex'}²?")
sq = math.isqrt(disc)
print(f"  √{disc} = {sq} (exact: {sq*sq == disc})")
if sq*sq == disc:
    print(f"  Zeros: x = ({-80+sq}/{14}) = {(-80+sq)/14:.6f} and x = ({-80-sq}/{14}) = {(-80-sq)/14:.6f}")
else:
    print(f"  Zeros: x ≈ {(-80 + math.sqrt(disc))/14:.6f} and {(-80 - math.sqrt(disc))/14:.6f}")

# For P_11:
print(f"\n  P_11: IP(G, x) = 1 + {c1_11}x + {c2_11}x² + {c3_total}x³")
# Zeros of cubic
import numpy as np
coeffs = [c3_total, c2_11, c1_11, 1]  # highest degree first
roots = np.roots(coeffs)
print(f"  Zeros:")
for r in roots:
    if abs(r.imag) < 1e-6:
        print(f"    x = {r.real:.6f}")
    else:
        print(f"    x = {r.real:.6f} ± {abs(r.imag):.6f}i")

print()
print("=" * 70)
print("PART 4: IP decomposition of H = 1 + 2c₁ + 4c₂ + 8c₃")
print("=" * 70)

# H(P_7) = 189 = 1 + 2(80) + 4(7) = 1 + 160 + 28 = 189
# H(P_11) = 95095 = 1 + 2(c1) + 4(c2) + 8(c3)

# The "level" contributions:
print(f"\n  P_7: H = 189")
print(f"    Level 0 (empty set):     1")
print(f"    Level 1 (single cycles): 2×80 = 160  ({160/189*100:.1f}%)")
print(f"    Level 2 (disjoint pairs): 4×7 = 28   ({28/189*100:.1f}%)")
print(f"    Total: 189 ✓")

print(f"\n  P_11: H = 95095")
l0 = 1
l1 = 2 * c1_11
l2 = 4 * c2_11
l3 = 8 * c3_total
print(f"    Level 0: {l0}")
print(f"    Level 1: 2×{c1_11} = {l1}  ({l1/95095*100:.1f}%)")
print(f"    Level 2: 4×{c2_11} = {l2}  ({l2/95095*100:.1f}%)")
print(f"    Level 3: 8×{c3_total} = {l3}  ({l3/95095*100:.1f}%)")
print(f"    Total: {l0+l1+l2+l3}")
print(f"    Matches H: {'✓' if l0+l1+l2+l3 == 95095 else '✗ (missing higher-order packings)'}")

print()
print("=" * 70)
print("PART 5: The 7 disjoint pairs in P_7 — orbits under Z/7Z")
print("=" * 70)

# The 7 disjoint pairs should form 1 orbit under Z/7Z (since 7 is prime).
# Let's check what types they are.
p = 7
print(f"\n  Disjoint pairs:")
for i, j in disjoint_pairs:
    ci = all_cycles[i]
    cj = all_cycles[j]
    print(f"    {set(ci)} ∥ {set(cj)} (lengths {len(ci)},{len(cj)}), vertex union = {set(ci)|set(cj)}")

# They must be pairs of 3-cycles (since 5+3=8 > 7 and 5+5=10 > 7)
# And they use 6 of 7 vertices, leaving one vertex out.
# Under Z/7Z shift, each pair maps to another pair.

print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)

print(f"""
  IP DECOMPOSITION OF H(P_p):

  P_7:  IP = 1 + 80x + 7x²
        H = 1 + 160 + 28 = 189
        Cycle contributions: 84.7% level 1, 14.8% level 2

  P_11: IP = 1 + {c1_11}x + {c2_11}x² + {c3_total}x³
        H = 1 + {l1} + {l2} + {l3} = {l0+l1+l2+l3}
        {'✓ matches 95095' if l0+l1+l2+l3 == 95095 else '✗ MISSING higher terms'}

  The IP degree = max packing number ≤ floor(p/3).

  As p grows, higher-order packing terms become more important.
  The convergence of H/E[H] may relate to how many packing
  levels contribute.
""")

print("Done!")
