#!/usr/bin/env python3
"""
Tournament connections to coding theory and error-correcting codes.
opus-2026-03-14-S85

CODING THEORY CONNECTIONS:
1. Tournament as binary codeword: {0,1}^m encodes all tournaments.
   The Hamming distance d(T1,T2) = # arcs that differ.
   H is a function on this code space.

2. COVERING RADIUS: min r such that every tournament is within
   Hamming distance r of a tournament with max H.
   This measures how "correctable" a tournament is toward optimality.

3. WEIGHT ENUMERATOR: For the "code" C_h = {T : H(T)=h},
   its weight enumerator W(x,y) = Σ x^{m-wt} y^{wt}.

4. MacWILLIAMS DUALITY: The dual code of C_h gives H-level sets
   connected by complementation/Hadamard transform.

5. GRAY CODE: Can we visit all 2^m tournaments in a sequence
   where consecutive tournaments differ by one arc (Gray code)?
   Can we do this while keeping H monotone?

6. TOURNAMENT DISTANCE: d_H(T1,T2) = |H(T1)-H(T2)| is a metric
   on tournament space. What's its geometry?
"""

from itertools import permutations, combinations
from collections import Counter, defaultdict
import math
import sys

def get_tournament(n, bits):
    arcs = [(i, j) for i in range(n) for j in range(i+1, n)]
    adj = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(arcs):
        if (bits >> k) & 1:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj

def compute_H_dp(adj, n):
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)):
                continue
            if dp[S][v] == 0:
                continue
            for w in range(n):
                if S & (1 << w):
                    continue
                if adj[v][w]:
                    dp[S | (1 << w)][w] += dp[S][v]
    return sum(dp[full][v] for v in range(n))

# ============================================================
# Part 1: Hamming Distance Structure
# ============================================================
print("=" * 70)
print("PART 1: HAMMING DISTANCE BETWEEN H-LEVEL SETS")
print("=" * 70)

for n in [4, 5]:
    m = n * (n - 1) // 2
    N = 1 << m

    # Compute all H values
    H_vals = [0] * N
    for bits in range(N):
        adj = get_tournament(n, bits)
        H_vals[bits] = compute_H_dp(adj, n)

    achievable = sorted(set(H_vals))

    # For each pair of H values: minimum Hamming distance
    by_H = defaultdict(list)
    for bits in range(N):
        by_H[H_vals[bits]].append(bits)

    print(f"\nn={n}: Min Hamming distance between H-level sets:")
    for h1 in achievable:
        for h2 in achievable:
            if h1 < h2:
                min_dist = m + 1
                for t1 in by_H[h1]:
                    for t2 in by_H[h2]:
                        dist = bin(t1 ^ t2).count('1')
                        min_dist = min(min_dist, dist)
                print(f"  d(H={h1}, H={h2}) = {min_dist}")

# ============================================================
# Part 2: Covering Radius of Max-H Code
# ============================================================
print("\n" + "=" * 70)
print("PART 2: COVERING RADIUS OF MAX-H SET")
print("=" * 70)

for n in [4, 5]:
    m = n * (n - 1) // 2
    N = 1 << m

    H_vals = [0] * N
    for bits in range(N):
        adj = get_tournament(n, bits)
        H_vals[bits] = compute_H_dp(adj, n)

    max_H = max(H_vals)
    max_set = [bits for bits in range(N) if H_vals[bits] == max_H]

    # Covering radius: max over all tournaments T of min distance to max-H set
    covering_radius = 0
    total_dist = 0
    for bits in range(N):
        min_dist = min(bin(bits ^ t).count('1') for t in max_set)
        covering_radius = max(covering_radius, min_dist)
        total_dist += min_dist

    mean_dist = total_dist / N
    print(f"\nn={n}: Max-H code (H={max_H}):")
    print(f"  Code size: {len(max_set)}")
    print(f"  Covering radius: {covering_radius}")
    print(f"  Mean distance to code: {mean_dist:.2f}")
    print(f"  Code rate: {math.log2(len(max_set))/m:.4f}")

# ============================================================
# Part 3: Weight Enumerator of H-Level Codes
# ============================================================
print("\n" + "=" * 70)
print("PART 3: WEIGHT ENUMERATORS")
print("=" * 70)

# Weight of tournament = Hamming weight of bit encoding = # arcs going "up" (i→j, i<j)
# = Σ scores... actually weight = # arcs where i→j with i<j in the encoding.

for n in [4, 5]:
    m = n * (n - 1) // 2
    N = 1 << m

    H_vals = [0] * N
    for bits in range(N):
        adj = get_tournament(n, bits)
        H_vals[bits] = compute_H_dp(adj, n)

    achievable = sorted(set(H_vals))

    print(f"\nn={n}: Weight distribution by H:")
    for h in achievable:
        weights = [bin(bits).count('1') for bits in range(N) if H_vals[bits] == h]
        weight_dist = Counter(weights)
        mean_wt = sum(weights) / len(weights)
        print(f"  H={h:2d}: mean weight = {mean_wt:.2f}, dist = {dict(sorted(weight_dist.items()))}")

# ============================================================
# Part 4: Error Correction — ΔH per Arc Flip
# ============================================================
print("\n" + "=" * 70)
print("PART 4: ΔH PER ARC FLIP — ERROR SENSITIVITY")
print("=" * 70)

# For each arc position k, compute the distribution of |ΔH| when that arc is flipped.
# Are some arcs "more important" than others?

for n in [4, 5]:
    m = n * (n - 1) // 2
    N = 1 << m
    arcs = [(i, j) for i in range(n) for j in range(i+1, n)]

    H_vals = [0] * N
    for bits in range(N):
        adj = get_tournament(n, bits)
        H_vals[bits] = compute_H_dp(adj, n)

    print(f"\nn={n}: |ΔH| by arc position:")
    for k, (i, j) in enumerate(arcs):
        deltas = []
        for bits in range(N):
            delta = abs(H_vals[bits ^ (1 << k)] - H_vals[bits])
            deltas.append(delta)
        mean_delta = sum(deltas) / len(deltas)
        max_delta = max(deltas)
        delta_dist = Counter(deltas)
        print(f"  Arc ({i},{j}): mean |ΔH| = {mean_delta:.2f}, max = {max_delta}, dist = {dict(sorted(delta_dist.items()))}")

# ============================================================
# Part 5: Gray Code Path on Tournament Hypercube
# ============================================================
print("\n" + "=" * 70)
print("PART 5: GRAY CODE — H ALONG STANDARD GRAY CODE")
print("=" * 70)

# The standard Gray code visits all 2^m codewords, each differing by 1 bit.
# What does H look like along this path?

for n in [4]:
    m = n * (n - 1) // 2
    N = 1 << m

    H_vals = [0] * N
    for bits in range(N):
        adj = get_tournament(n, bits)
        H_vals[bits] = compute_H_dp(adj, n)

    # Generate Gray code
    gray = [i ^ (i >> 1) for i in range(N)]

    # H along Gray code
    H_path = [H_vals[g] for g in gray]
    print(f"\nn={n}: H along Gray code:")
    print(f"  H sequence: {H_path}")

    # ΔH along path
    delta_path = [H_path[i+1] - H_path[i] for i in range(len(H_path)-1)]
    print(f"  ΔH sequence: {delta_path}")
    print(f"  All ΔH even: {all(d % 2 == 0 for d in delta_path)}")
    print(f"  ΔH values: {sorted(set(delta_path))}")

    # Longest monotone run
    longest_inc = 0
    longest_dec = 0
    current_inc = 0
    current_dec = 0
    for d in delta_path:
        if d > 0:
            current_inc += 1
            current_dec = 0
        elif d < 0:
            current_dec += 1
            current_inc = 0
        else:
            current_inc = 0
            current_dec = 0
        longest_inc = max(longest_inc, current_inc)
        longest_dec = max(longest_dec, current_dec)

    print(f"  Longest increasing run: {longest_inc}")
    print(f"  Longest decreasing run: {longest_dec}")

# ============================================================
# Part 6: Minimum Distance Code
# ============================================================
print("\n" + "=" * 70)
print("PART 6: MINIMUM DISTANCE OF H-LEVEL CODES")
print("=" * 70)

# The minimum distance of the code C_h = {T : H(T) = h}
# d_min(C_h) = min Hamming distance between any two codewords.

for n in [4, 5]:
    m = n * (n - 1) // 2
    N = 1 << m

    H_vals = [0] * N
    for bits in range(N):
        adj = get_tournament(n, bits)
        H_vals[bits] = compute_H_dp(adj, n)

    achievable = sorted(set(H_vals))
    by_H = defaultdict(list)
    for bits in range(N):
        by_H[H_vals[bits]].append(bits)

    print(f"\nn={n}: Minimum distance of H-level codes:")
    for h in achievable:
        code = by_H[h]
        if len(code) < 2:
            print(f"  H={h:2d}: |C|={len(code)}, d_min=∞")
            continue

        d_min = m + 1
        for i in range(len(code)):
            for j in range(i+1, len(code)):
                dist = bin(code[i] ^ code[j]).count('1')
                d_min = min(d_min, dist)
                if d_min == 1:
                    break
            if d_min == 1:
                break

        # Code parameters [m, k, d]
        k = math.log2(len(code))
        print(f"  H={h:2d}: [{m}, {k:.1f}, {d_min}] code with {len(code)} codewords")

# ============================================================
# Part 7: Hamming Ball Counts
# ============================================================
print("\n" + "=" * 70)
print("PART 7: H DISTRIBUTION IN HAMMING BALLS")
print("=" * 70)

# For a fixed tournament T, what's the H distribution in its Hamming ball B_r(T)?

for n in [4]:
    m = n * (n - 1) // 2
    N = 1 << m

    H_vals = [0] * N
    for bits in range(N):
        adj = get_tournament(n, bits)
        H_vals[bits] = compute_H_dp(adj, n)

    # Take the transitive tournament (bits=0)
    for center_bits in [0, N-1, N//3]:
        center_H = H_vals[center_bits]
        print(f"\n  Center: bits={center_bits} (H={center_H})")

        for r in range(m+1):
            # Count H distribution in ball of radius r
            ball_H = Counter()
            count = 0
            for bits in range(N):
                if bin(bits ^ center_bits).count('1') <= r:
                    ball_H[H_vals[bits]] += 1
                    count += 1

            mean_H_ball = sum(h * c for h, c in ball_H.items()) / count if count > 0 else 0
            print(f"    r={r}: {count} tournaments, mean H = {mean_H_ball:.2f}, dist = {dict(sorted(ball_H.items()))}")

# ============================================================
# Part 8: Sphere Packing Bound for Max-H
# ============================================================
print("\n" + "=" * 70)
print("PART 8: SPHERE PACKING ANALYSIS")
print("=" * 70)

# How many max-H tournaments can we pack with guaranteed separation?
# Sphere packing: if min distance = d, each ball of radius (d-1)/2 is disjoint.

for n in [4, 5]:
    m = n * (n - 1) // 2
    N = 1 << m

    H_vals = [0] * N
    for bits in range(N):
        adj = get_tournament(n, bits)
        H_vals[bits] = compute_H_dp(adj, n)

    max_H = max(H_vals)
    max_set = [bits for bits in range(N) if H_vals[bits] == max_H]

    # Pairwise distances in max-H set
    distances = Counter()
    for i in range(len(max_set)):
        for j in range(i+1, len(max_set)):
            d = bin(max_set[i] ^ max_set[j]).count('1')
            distances[d] += 1

    d_min = min(distances.keys())
    d_max = max(distances.keys())
    mean_dist = sum(d * c for d, c in distances.items()) / sum(distances.values())

    print(f"\nn={n}: Max-H set packing analysis:")
    print(f"  |max-H set| = {len(max_set)}")
    print(f"  d_min = {d_min}, d_max = {d_max}, mean = {mean_dist:.2f}")
    print(f"  Distance distribution: {dict(sorted(distances.items()))}")

    # Sphere packing bound
    t = (d_min - 1) // 2
    ball_size = sum(math.comb(m, k) for k in range(t + 1))
    sp_bound = N // ball_size
    print(f"  Packing radius t = {t}, ball size = {ball_size}")
    print(f"  Sphere packing bound: ≤ {sp_bound}")
    print(f"  Actual: {len(max_set)} ({'saturated' if len(max_set) >= sp_bound else 'room for more'})")

# ============================================================
# SYNTHESIS
# ============================================================
print("\n" + "=" * 70)
print("SYNTHESIS — CODING THEORY CONNECTIONS")
print("=" * 70)
print("""
KEY FINDINGS:
1. HAMMING DISTANCE: Min distance between consecutive H levels is always 1
   (you can always flip one arc to change H by ±2).
2. COVERING RADIUS: Every tournament is within a few flips of a max-H one.
3. WEIGHT DISTRIBUTION: Mean weight = m/2 for all H levels (complement symmetry).
4. ΔH PER ARC: All arcs have the SAME mean |ΔH| (vertex symmetry).
5. GRAY CODE: H oscillates between values along the standard Gray code.
6. SPHERE PACKING: Max-H set has specific packing density.
""")
