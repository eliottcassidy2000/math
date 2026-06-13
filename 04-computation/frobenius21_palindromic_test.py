#!/usr/bin/env python3
"""
Test palindromic N property for the NON-self-converse F_21 tournament.

If palindromic N FAILS for this VT tournament, then:
- M =/= (H/n)*I for this tournament
- THM-052 genuinely requires self-converse (or something beyond VT)

If palindromic N HOLDS despite not being self-converse, then:
- There's a different mechanism at work (not reflection-reversal)
- THM-052 might still hold for all VT at odd n

We compute N(0, b, j) = #{Ham paths with 0 at position j and b at position j+1,
                         OR 0 at position j+1 and b at position j}
via DP (bitmask DP, 2^21 * 21 states).

kind-pasteur-2026-03-06-S25e
"""

import time
import sys

POW2 = [1, 2, 4]

def mul(a, b):
    return ((a[0] + POW2[a[1]] * b[0]) % 7, (a[1] + b[1]) % 3)

def inv(a):
    j_inv = (-a[1]) % 3
    i_inv = (-POW2[j_inv] * a[0]) % 7
    return (i_inv, j_inv)

elements = [(i, j) for j in range(3) for i in range(7)]
idx = {e: k for k, e in enumerate(elements)}
n = 21

nonid = [e for e in elements if e != (0, 0)]
pairs = []
seen = set()
for g in nonid:
    if g not in seen:
        pairs.append((g, inv(g)))
        seen.add(g)
        seen.add(inv(g))

# Non-normal connection set (the one confirmed NOT self-converse)
S_nonnormal = set()
for k, (g, gi) in enumerate(pairs):
    if k < 5:
        S_nonnormal.add(g)
    else:
        S_nonnormal.add(gi)

# Build adjacency as compact array
adj = [[False]*n for _ in range(n)]
for a in elements:
    for b in elements:
        if a == b:
            continue
        diff = mul(inv(a), b)
        if diff in S_nonnormal:
            adj[idx[a]][idx[b]] = True

# Convert to bit adjacency for fast lookup
adj_bits = [0] * n
for i in range(n):
    for j in range(n):
        if adj[i][j]:
            adj_bits[i] |= (1 << j)

print(f"n={n}, tournament on F_21 (Frobenius)")
print(f"Connection set (non-normal): {sorted(S_nonnormal)}")

# Step 1: Compute total H(T) via DP
# dp[mask][last] = number of paths using vertices in mask, ending at last
print("\nStep 1: Computing H(T)...")
t0 = time.time()

# Use arrays for speed
full_mask = (1 << n) - 1
# dp[mask][last] — too much memory for 2^21 * 21 = 44M entries as dict
# Use arrays instead. 2^21 = 2M, * 21 = 44M ints. At 8 bytes each = 352 MB.
# That's a lot. Let me use a more memory-efficient approach.

# Actually, let's use the DP but store only what we need.
# For palindromic N test, we need to count paths where specific vertices
# appear at specific positions.

# For N(a, b, j): count paths P where {P[j], P[j+1]} = {a, b}
# (both orderings counted).

# More efficient: compute f(d, j) where d = idx of some non-identity element,
# using the circulant symmetry... wait, this is NOT circulant! F_21 is non-abelian.

# For VT, we have: N(g*a, g*b, j) = N(a, b, j) for all g in G.
# So N(a, b, j) depends on the "relative position" of a and b.
# For Cayley, this means N(a, b, j) = N(e, a^{-1}b, j) = f(a^{-1}b, j).

# But for non-abelian groups, we also need: N(a, b, j) = N(e, a^{-1}b, j).
# Let's verify this and just compute f(d, j) for a few values of d.

# The DP approach: we need dp[mask][last] for ALL masks and last vertices.
# Total memory: 2^21 * 21 entries. At 4 bytes per int = 176 MB.
# Let's try with numpy arrays for memory efficiency.

try:
    import numpy as np
    USE_NUMPY = True
except ImportError:
    USE_NUMPY = False

if USE_NUMPY:
    # Use numpy for memory efficiency
    # dp[last] is an array of size 2^21, indexed by mask
    # We process by popcount (number of bits set in mask)

    # Actually, we need dp[mask][last] but mask has 2^21 entries.
    # Let's use a dict of dicts: dp[(mask, last)] = count
    # With popcount-based iteration, we only need to store the current and next popcounts.
    pass

# For n=21, the full DP with 2^21 * 21 entries is ~44M entries.
# Using a flat array of int32 would be 176 MB. Let's try.

# Actually, for just counting H and N(0, b, j), we can be smarter:
# Use the forward DP and track the position of vertex 0 and the connection
# to b in each path.

# SIMPLIFIED APPROACH: Just compute H(T) first, then N(0, b, j) for a few b.

# For H(T): use the standard bitmask DP.
# dp[mask][v] = number of Ham paths through exactly the vertices in mask, ending at v.

print("Allocating DP array...")
sys.stdout.flush()

# Use dictionary-based DP (memory efficient, slower)
dp = {}
for v in range(n):
    dp[(1 << v, v)] = 1

print("Running DP...")
sys.stdout.flush()

for pop in range(2, n + 1):
    count = 0
    t_pop = time.time()
    for mask in range(1, 1 << n):
        if bin(mask).count('1') != pop:
            continue
        for last in range(n):
            if not ((mask >> last) & 1):
                continue
            prev_mask = mask ^ (1 << last)
            total = 0
            # Try all predecessors of last
            remaining = prev_mask
            while remaining:
                prev = (remaining & -remaining).bit_length() - 1
                if adj[prev][last] and (prev_mask, prev) in dp:
                    total += dp[(prev_mask, prev)]
                remaining &= remaining - 1
            if total > 0:
                dp[(mask, last)] = total
                count += 1

    elapsed = time.time() - t_pop
    # Only print every few popcounts to avoid spam
    if pop <= 5 or pop >= n - 2 or pop % 3 == 0:
        mem = len(dp) * 32  # rough estimate: 32 bytes per dict entry
        print(f"  pop={pop:2d}: {count:8d} nonzero entries, {elapsed:.1f}s, ~{mem // (1024*1024)}MB")
    sys.stdout.flush()

# Compute H
H = sum(dp.get((full_mask, v), 0) for v in range(n))
print(f"\nH(T) = {H}")
print(f"H/n = {H/n}" if H % n == 0 else f"H/n = {H}/{n} (NOT divisible!)")
print(f"Total DP time: {time.time() - t0:.1f}s")

# Step 2: Compute N(0, b, j) for a few values of b
# N(a, b, j) counts paths P where {P[j], P[j+1]} = {a, b}
# = C(a, b, j) + C(b, a, j) where C(a,b,j) = #{P: P[j]=a, P[j+1]=b, adj[a][b]=1}

# We need to extract from dp: for a path ending at vertex v through mask,
# what is the vertex at each position? We can't do this from the forward DP alone.

# Instead: use TWO DPs:
# prefix_dp[mask][last] = # paths through mask ending at last
# suffix_dp[mask][first] = # paths through mask starting at first

# N(a, b, j) = sum over (prefix ending at a of length j+1) * (suffix starting at b of length n-j-1)
# where the prefix uses a subset of vertices and suffix uses the complement.
# And adj[a][b] = 1 (or adj[b][a] = 1 for the other direction).

# We already have prefix_dp = dp. Now compute suffix_dp.
print("\nComputing suffix DP...")
t_suf = time.time()
sys.stdout.flush()

# suffix_dp[mask][first] = # paths through exactly the vertices in mask, starting at first
sdp = {}
for v in range(n):
    sdp[(1 << v, v)] = 1

for pop in range(2, n + 1):
    t_pop = time.time()
    for mask in range(1, 1 << n):
        if bin(mask).count('1') != pop:
            continue
        for first in range(n):
            if not ((mask >> first) & 1):
                continue
            rest = mask ^ (1 << first)
            total = 0
            remaining = rest
            while remaining:
                nxt = (remaining & -remaining).bit_length() - 1
                if adj[first][nxt] and (rest, nxt) in sdp:
                    total += sdp[(rest, nxt)]
                remaining &= remaining - 1
            if total > 0:
                sdp[(mask, first)] = total

    if pop <= 5 or pop >= n - 2 or pop % 3 == 0:
        elapsed = time.time() - t_pop
        print(f"  pop={pop:2d}: {elapsed:.1f}s")
    sys.stdout.flush()

print(f"Suffix DP time: {time.time() - t_suf:.1f}s")

# Now compute N(0, b, j) for b = 1 (index 1)
a = 0  # vertex (0,0), index 0
b_vertex = elements[1]  # vertex (1,0), index 1
b = 1

print(f"\nComputing N({elements[a]}, {b_vertex}, j)...")

N = [0] * (n - 1)  # N[j] for j = 0, ..., n-2

for j in range(n - 1):
    prefix_len = j + 1  # positions 0..j, so j+1 vertices
    suffix_len = n - j - 1  # positions j+1..n-1, so n-j-1 vertices

    # Case 1: P[j] = a, P[j+1] = b (need adj[a][b] = 1)
    if adj[a][b]:
        for prefix_mask in range(1, 1 << n):
            if bin(prefix_mask).count('1') != prefix_len:
                continue
            if not ((prefix_mask >> a) & 1):
                continue
            if (prefix_mask >> b) & 1:
                continue  # b should not be in prefix
            pcnt = dp.get((prefix_mask, a), 0)
            if pcnt == 0:
                continue
            suffix_mask = full_mask ^ prefix_mask
            scnt = sdp.get((suffix_mask, b), 0)
            N[j] += pcnt * scnt

    # Case 2: P[j] = b, P[j+1] = a (need adj[b][a] = 1)
    if adj[b][a]:
        for prefix_mask in range(1, 1 << n):
            if bin(prefix_mask).count('1') != prefix_len:
                continue
            if not ((prefix_mask >> b) & 1):
                continue
            if (prefix_mask >> a) & 1:
                continue  # a should not be in prefix
            pcnt = dp.get((prefix_mask, b), 0)
            if pcnt == 0:
                continue
            suffix_mask = full_mask ^ prefix_mask
            scnt = sdp.get((suffix_mask, a), 0)
            N[j] += pcnt * scnt

    if j <= 3 or j >= n - 4 or j % 4 == 0:
        print(f"  j={j:2d}: N = {N[j]}")
    sys.stdout.flush()

print(f"\nFull N sequence: {N}")

# Check palindromic
is_pal = all(N[j] == N[n-2-j] for j in range(n-1))
print(f"Palindromic? {is_pal}")

if not is_pal:
    print("N is NOT palindromic!")
    for j in range(n-1):
        if N[j] != N[n-2-j]:
            print(f"  N[{j}] = {N[j]} != N[{n-2-j}] = {N[n-2-j]}")
    alt_sum = sum((-1)**j * N[j] for j in range(n-1))
    print(f"Alternating sum = {alt_sum}")
    if alt_sum != 0:
        print("M[0,1] != 0 => M is NOT scalar!")
    else:
        print("Alternating sum is 0 despite non-palindromic (fortuitous cancellation?)")
else:
    alt_sum = sum((-1)**j * N[j] for j in range(n-1))
    print(f"Alternating sum = {alt_sum} (should be 0 for odd n)")

print(f"\nTotal time: {time.time() - t0:.1f}s")
print("DONE")
