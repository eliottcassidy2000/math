#!/usr/bin/env python3
"""
Compute N(0, b, j) for the non-SC F_21 tournament to test palindromic property.

From frobenius21_ham_count.py:
  H(T_nonnormal) = 123,522,430,238,361
  H / 21 = 5,882,020,487,541

We compute N(0, b, j) for a specific non-edge pair (0, b) using:
  N(a,b,j) = C(a,b,j) + C(b,a,j)
where C(a,b,j) = #{Ham paths with a at position j and b at position j+1, with edge a->b}

This requires two DPs:
  prefix_dp[mask][last] = # paths through mask ending at last
  suffix_dp[mask][first] = # paths through mask starting at first

Then: C(a,b,j) = sum over prefix_masks of size j+1 containing a but not b,
                  prefix_dp[mask][a] * suffix_dp[complement][b]
       (only if adj[a][b] = 1)

At n=21, each DP has ~7M entries (peak at pop~12). Total time ~12min per DP.

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
idx_map = {e: k for k, e in enumerate(elements)}
n = 21
full_mask = (1 << n) - 1

nonid = [e for e in elements if e != (0, 0)]
pairs = []
seen = set()
for g in nonid:
    if g not in seen:
        pairs.append((g, inv(g)))
        seen.add(g)
        seen.add(inv(g))

# Non-normal connection set
S_nonnormal = set()
for k, (g, gi) in enumerate(pairs):
    if k < 5:
        S_nonnormal.add(g)
    else:
        S_nonnormal.add(gi)

adj = [[False]*n for _ in range(n)]
for a in elements:
    for b in elements:
        if a == b:
            continue
        diff = mul(inv(a), b)
        if diff in S_nonnormal:
            adj[idx_map[a]][idx_map[b]] = True

print(f"n={n}, F_21 tournament (non-normal S)")
print(f"H = 123522430238361, H/21 = 5882020487541")
sys.stdout.flush()

# Step 1: Forward DP (prefix paths)
print("\nStep 1: Forward DP (prefix paths)...")
t0 = time.time()

# Store DP by popcount layers, keeping ALL layers for later use
prefix = [{} for _ in range(n + 1)]
for v in range(n):
    prefix[1][(1 << v, v)] = 1

for pop in range(2, n + 1):
    tp = time.time()
    for (mask, last), cnt in prefix[pop - 1].items():
        for nxt in range(n):
            if (mask >> nxt) & 1:
                continue
            if not adj[last][nxt]:
                continue
            new_mask = mask | (1 << nxt)
            key = (new_mask, nxt)
            if key in prefix[pop]:
                prefix[pop][key] += cnt
            else:
                prefix[pop][key] = cnt
    if pop <= 5 or pop >= n - 3 or pop % 4 == 0:
        print(f"  pop={pop:2d}: {len(prefix[pop]):8d} entries, {time.time()-t0:.1f}s")
    sys.stdout.flush()

# Verify H
H = sum(prefix[n].get((full_mask, v), 0) for v in range(n))
print(f"  H = {H}")
assert H == 123522430238361, f"H mismatch: {H}"

# Step 2: Backward DP (suffix paths)
print("\nStep 2: Backward DP (suffix paths)...")
t1 = time.time()

suffix = [{} for _ in range(n + 1)]
for v in range(n):
    suffix[1][(1 << v, v)] = 1

for pop in range(2, n + 1):
    for (mask, first), cnt in suffix[pop - 1].items():
        for prev in range(n):
            if (mask >> prev) & 1:
                continue
            if not adj[prev][first]:
                continue
            new_mask = mask | (1 << prev)
            key = (new_mask, prev)
            if key in suffix[pop]:
                suffix[pop][key] += cnt
            else:
                suffix[pop][key] = cnt
    if pop <= 5 or pop >= n - 3 or pop % 4 == 0:
        print(f"  pop={pop:2d}: {len(suffix[pop]):8d} entries, {time.time()-t1:.1f}s")
    sys.stdout.flush()

# Verify H from suffix
H2 = sum(suffix[n].get((full_mask, v), 0) for v in range(n))
print(f"  H (from suffix) = {H2}")
assert H2 == H, f"H mismatch in suffix: {H2}"

# Step 3: Compute N(0, b, j)
# Choose b = 1 (vertex (1,0))
a = 0  # vertex (0,0)
b = 1  # vertex (1,0)
print(f"\nStep 3: Computing N(v{a}, v{b}, j) = N({elements[a]}, {elements[b]}, j)")
print(f"  adj[{a}][{b}] = {adj[a][b]}, adj[{b}][{a}] = {adj[b][a]}")

N = [0] * (n - 1)

for j in range(n - 1):
    prefix_len = j + 1
    suffix_len = n - j - 1

    # Case 1: P[j] = a, P[j+1] = b (need adj[a][b])
    if adj[a][b]:
        for (mask, last), cnt in prefix[prefix_len].items():
            if last != a:
                continue
            if (mask >> b) & 1:
                continue
            suffix_mask = full_mask ^ mask
            scnt = suffix[suffix_len].get((suffix_mask, b), 0)
            if scnt > 0:
                N[j] += cnt * scnt

    # Case 2: P[j] = b, P[j+1] = a (need adj[b][a])
    if adj[b][a]:
        for (mask, last), cnt in prefix[prefix_len].items():
            if last != b:
                continue
            if (mask >> a) & 1:
                continue
            suffix_mask = full_mask ^ mask
            scnt = suffix[suffix_len].get((suffix_mask, a), 0)
            if scnt > 0:
                N[j] += cnt * scnt

    if j <= 4 or j >= n - 5 or j % 3 == 0:
        print(f"  j={j:2d}: N = {N[j]}")
    sys.stdout.flush()

print(f"\nFull N({elements[a]}, {elements[b]}, j) sequence:")
for j in range(n - 1):
    pal_j = n - 2 - j
    match = "=" if N[j] == N[pal_j] else "!="
    print(f"  j={j:2d}: N={N[j]:>18d}   (j'={pal_j:2d}: N={N[pal_j]:>18d}) {match}")

is_pal = all(N[j] == N[n-2-j] for j in range(n-1))
print(f"\nPalindromic? {is_pal}")

alt_sum = sum((-1)**j * N[j] for j in range(n-1))
print(f"Alternating sum = {alt_sum}")

if is_pal:
    print("PALINDROMIC => M[a,b] = 0 (scalar M holds for this pair!)")
else:
    if alt_sum == 0:
        print("NOT palindromic but alternating sum = 0 (M[a,b] = 0 anyway)")
    else:
        print(f"NOT palindromic and M[a,b] = {alt_sum} != 0")
        print("M is NOT scalar for this VT tournament!")

print(f"\nTotal time: {time.time()-t0:.1f}s")
print("DONE")
