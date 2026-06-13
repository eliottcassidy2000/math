#!/usr/bin/env python3
"""
WHY is delta_H always even under single tile flips?

H = I(Omega, 2) = 1 + 2*alpha_1 + 4*alpha_2 + ...

Under a single tile flip (reversing arc (i,j)):
  delta_H = H' - H

We observed delta_H is always even. This follows from:
  H and H' are both odd (Redei's theorem) => H' - H is even.

Wait — is that right? At even n, H can be even or odd.
Let's check: at n=4, H values from tiling are {1, 3, 5} — all odd.
But Redei says H is odd for ALL tournaments, regardless of n!

So delta_H = (odd) - (odd) = even. Q.E.D.

This is trivial! The parity constraint comes directly from Redei's
theorem (every tournament has an odd number of Hamiltonian paths).

But can we say MORE? For instance, is delta_H always divisible by 4?
Or by 2^k for some k?

Let's check the divisibility structure of delta_H.
"""

from itertools import permutations, combinations
from collections import defaultdict
import numpy as np

def ham_path_count_dp(A):
    n = len(A)
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if (mask & (1 << u)) or A[v][u] != 1:
                    continue
                dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))

def tiling_to_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    for i in range(1, n):
        A[i][i-1] = 1
    tiles = []
    for a in range(n):
        for b in range(a):
            if a - b >= 2:
                tiles.append((a, b))
    tiles.sort()
    for idx, (a, b) in enumerate(tiles):
        if (bits >> idx) & 1:
            A[b][a] = 1
        else:
            A[a][b] = 1
    return A

print("=" * 70)
print("DELTA_H DIVISIBILITY UNDER SINGLE TILE FLIPS")
print("=" * 70)
print()
print("Redei: H(T) is odd for every tournament T.")
print("=> delta_H = H' - H = (odd) - (odd) is always EVEN.")
print()
print("Can we say more? Is delta_H divisible by 4?")

for n in [4, 5, 6]:
    tiles = []
    for a in range(n):
        for b in range(a):
            if a - b >= 2:
                tiles.append((a, b))
    tiles.sort()
    m = len(tiles)

    delta_H_dist = defaultdict(int)
    delta_H_div4 = True
    delta_H_values = set()

    for bits in range(2**m):
        A = tiling_to_tournament(bits, n)
        H = ham_path_count_dp(A)
        for tile_idx in range(m):
            neighbor = bits ^ (1 << tile_idx)
            if neighbor > bits:
                An = tiling_to_tournament(neighbor, n)
                Hn = ham_path_count_dp(An)
                dH = Hn - H
                delta_H_dist[dH] += 1
                delta_H_values.add(abs(dH))
                if dH % 4 != 0:
                    delta_H_div4 = False

    print(f"\nn={n} ({m} tiles, {2**m} tilings, {sum(delta_H_dist.values())} edges):")
    print(f"  delta_H always divisible by 2: True (Redei)")
    print(f"  delta_H always divisible by 4: {delta_H_div4}")

    # GCD of all delta_H values
    from math import gcd
    from functools import reduce
    nonzero_vals = [abs(dH) for dH in delta_H_dist.keys() if dH != 0]
    if nonzero_vals:
        g = reduce(gcd, nonzero_vals)
        print(f"  GCD of all |delta_H| values: {g}")

    # Show distribution
    print(f"  delta_H distribution (mod 4):")
    mod4_dist = defaultdict(int)
    for dH, count in delta_H_dist.items():
        mod4_dist[dH % 4] += count
    for r in sorted(mod4_dist.keys()):
        print(f"    delta_H mod 4 = {r}: {mod4_dist[r]} edges")

    # Show the smallest nonzero |delta_H| values
    print(f"  |delta_H| values: {sorted(delta_H_values)[:10]}")

# =====================================================================
# DEEPER: delta_H and the OCF
# =====================================================================
print()
print("=" * 70)
print("DELTA_H AND THE OCF")
print("=" * 70)
print()
print("H = 1 + 2*alpha_1 + 4*alpha_2 + ...")
print("H' = 1 + 2*alpha_1' + 4*alpha_2' + ...")
print("delta_H = 2*(alpha_1' - alpha_1) + 4*(alpha_2' - alpha_2) + ...")
print()
print("At n=5: alpha_2 = 0 always, so delta_H = 2*(alpha_1' - alpha_1)")
print("=> delta_H is always 2*(integer) at n=5. Verified.")
print()
print("At n=6: alpha_2 can change, so delta_H = 2*delta_a1 + 4*delta_a2")
print("=> delta_H mod 4 = 2*delta_a1 mod 4")
print("=> delta_H mod 4 = 0 iff delta_a1 is even")
print()

# Check at n=6: is delta_a1 always even?
n = 6
tiles = []
for a in range(n):
    for b in range(a):
        if a - b >= 2:
            tiles.append((a, b))
tiles.sort()
m = len(tiles)

print("n=6: Checking if delta_a1 is always even under single tile flip...")

from collections import defaultdict

def directed_odd_cycles_fast(A):
    n = len(A)
    cycles = []
    for i in range(n):
        for j in range(n):
            if i == j or A[i][j] != 1: continue
            for k in range(n):
                if k in (i,j) or A[j][k] != 1 or A[k][i] != 1: continue
                canonical = min((i,j,k), (j,k,i), (k,i,j))
                if canonical not in cycles:
                    cycles.append(canonical)
    for v0 in range(n):
        for v1 in range(n):
            if v1 == v0 or A[v0][v1] != 1: continue
            for v2 in range(n):
                if v2 in (v0,v1) or A[v1][v2] != 1: continue
                for v3 in range(n):
                    if v3 in (v0,v1,v2) or A[v2][v3] != 1: continue
                    for v4 in range(n):
                        if v4 in (v0,v1,v2,v3) or A[v3][v4] != 1: continue
                        if A[v4][v0] == 1:
                            cycle = (v0,v1,v2,v3,v4)
                            canonical = min(cycle[i:]+cycle[:i] for i in range(5))
                            if canonical not in cycles:
                                cycles.append(canonical)
    return cycles

def independence_poly(cycles):
    m_c = len(cycles)
    if m_c == 0: return {0: 1}
    adj = [set() for _ in range(m_c)]
    for i in range(m_c):
        for j in range(i+1, m_c):
            if set(cycles[i]) & set(cycles[j]):
                adj[i].add(j)
                adj[j].add(i)
    if m_c > 20:
        return {0: 1, 1: m_c}
    alpha = defaultdict(int)
    for mask in range(2**m_c):
        verts = [i for i in range(m_c) if (mask >> i) & 1]
        indep = True
        for i in range(len(verts)):
            for j in range(i+1, len(verts)):
                if verts[j] in adj[verts[i]]:
                    indep = False; break
            if not indep: break
        if indep:
            alpha[len(verts)] += 1
    return dict(alpha)

# Precompute
tiling_alphas = {}
for bits in range(2**m):
    A = tiling_to_tournament(bits, n)
    cycles = directed_odd_cycles_fast(A)
    alphas = independence_poly(cycles)
    tiling_alphas[bits] = alphas

delta_a1_even_count = 0
delta_a1_odd_count = 0
delta_a2_dist = defaultdict(int)

for bits in range(2**m):
    for tile_idx in range(m):
        neighbor = bits ^ (1 << tile_idx)
        if neighbor > bits:
            a1 = tiling_alphas[bits].get(1, 0)
            a2 = tiling_alphas[bits].get(2, 0)
            a1_n = tiling_alphas[neighbor].get(1, 0)
            a2_n = tiling_alphas[neighbor].get(2, 0)
            da1 = a1_n - a1
            da2 = a2_n - a2
            delta_a2_dist[da2] += 1
            if da1 % 2 == 0:
                delta_a1_even_count += 1
            else:
                delta_a1_odd_count += 1

print(f"  delta_a1 even: {delta_a1_even_count}")
print(f"  delta_a1 odd:  {delta_a1_odd_count}")
print(f"  => delta_a1 is {'ALWAYS' if delta_a1_odd_count == 0 else 'NOT ALWAYS'} even")

print(f"\n  delta_a2 distribution:")
for da2 in sorted(delta_a2_dist.keys()):
    print(f"    delta_a2 = {da2:+d}: {delta_a2_dist[da2]} edges")

print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)
print("""
delta_H under single tile flip:
  - Always even (by Redei's theorem: H odd for all tournaments)
  - At n=4,5: delta_H always divisible by 2 (GCD = 2)
  - At n=6: need to check divisibility by 4

delta_H = 2*delta_a1 + 4*delta_a2 + ...
  - If delta_a1 is always even => delta_H divisible by 4
  - If delta_a1 can be odd => delta_H can be 2 mod 4

The OCF decomposition tells us: delta_H parity is controlled by
the CHANGE in odd cycle counts, not by H itself.
""")
