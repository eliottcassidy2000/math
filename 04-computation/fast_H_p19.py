#!/usr/bin/env python3
"""
Fast Hamiltonian path count for p=19 using optimized DP with numpy.
The key optimization: vectorize the inner loop over vertices.

opus-2026-03-12-S62d
"""
import numpy as np
import time
import sys

def make_tournament(p, S):
    A = np.zeros((p, p), dtype=np.int8)
    for i in range(p):
        for s in S:
            A[i][(i + s) % p] = 1
    return A

def count_H_fast(A):
    """Optimized DP for Hamiltonian paths using numpy arrays."""
    n = len(A)
    # dp[mask][v] = number of Hamiltonian paths ending at v using vertices in mask
    # Use dict of arrays for sparse representation
    full = (1 << n) - 1

    # Initialize: single vertex paths
    dp = {}
    for i in range(n):
        mask = 1 << i
        if mask not in dp:
            dp[mask] = np.zeros(n, dtype=np.int64)
        dp[mask][i] = 1

    # Process masks in order of popcount
    for popcount in range(1, n):
        new_dp = {}
        for mask, counts in dp.items():
            if bin(mask).count('1') != popcount:
                continue
            for u in range(n):
                if counts[u] == 0:
                    continue
                if not (mask & (1 << u)):
                    continue
                for v in range(n):
                    if mask & (1 << v):
                        continue
                    if A[u][v]:
                        new_mask = mask | (1 << v)
                        if new_mask not in new_dp:
                            new_dp[new_mask] = np.zeros(n, dtype=np.int64)
                        new_dp[new_mask][v] += counts[u]

        # Merge new_dp into dp
        for mask, counts in new_dp.items():
            if mask in dp:
                dp[mask] += counts
            else:
                dp[mask] = counts

    if full in dp:
        return int(dp[full].sum())
    return 0

def count_H_array(A):
    """Fully array-based DP. Uses 2^n * n array (needs ~10GB for n=19)."""
    n = len(A)
    if n > 20:
        print(f"  n={n} too large for array-based DP (needs {2**n * n * 8 / 1e9:.1f} GB)")
        return -1

    full = (1 << n) - 1
    # Process by popcount layer to save memory
    # Layer k: all masks with popcount k
    # Only need current and next layer

    prev = {}  # mask -> array of counts per ending vertex
    for i in range(n):
        prev[1 << i] = np.zeros(n, dtype=np.int64)
        prev[1 << i][i] = 1

    for layer in range(2, n + 1):
        curr = {}
        for mask, counts in prev.items():
            for u in range(n):
                if counts[u] == 0:
                    continue
                # Extend from u to each neighbor v not in mask
                for v in range(n):
                    if (mask & (1 << v)) == 0 and A[u][v]:
                        new_mask = mask | (1 << v)
                        if new_mask not in curr:
                            curr[new_mask] = np.zeros(n, dtype=np.int64)
                        curr[new_mask][v] += counts[u]

        # Keep full mask if we're at the last layer
        if layer == n:
            if full in curr:
                return int(curr[full].sum())
            else:
                return 0
        prev = curr

    return 0

def get_QR(p):
    return sorted(set(pow(a, 2, p) for a in range(1, p)) - {0})

p = 19
m = (p - 1) // 2
QR = get_QR(p)
S_int = list(range(1, m + 1))

print(f"p = {p}, m = {m}")
print(f"QR = {QR}")
print(f"Interval = {S_int}")
print(f"2^{p} = {2**p} states")
print()

# Try the layer-by-layer approach
for name, S in [("Interval", S_int), ("Paley", QR)]:
    A = make_tournament(p, S)
    print(f"Computing H({name})...", flush=True)
    t0 = time.time()
    H = count_H_array(A)
    t1 = time.time()
    print(f"  H({name}) = {H} (computed in {t1-t0:.1f}s)")
    print()

# Now do single-swap analysis from Interval
print("=" * 72)
print("SINGLE-SWAP ANALYSIS FROM INTERVAL AT p=19")
print("=" * 72)
print()

A_int = make_tournament(p, S_int)
H_int = count_H_array(A_int)
print(f"H(Interval) = {H_int}")

pairs = [(j, p-j) for j in range(1, m+1)]
swaps = []

for a in S_int:
    for b in range(1, p):
        if b not in S_int:
            S_new = sorted([x for x in S_int if x != a] + [b])
            if len(S_new) == m and len(set(S_new)) == m:
                # Check valid tournament set
                valid = all((j in S_new) != ((p-j) in S_new) for j in range(1, m+1))
                if valid:
                    A_new = make_tournament(p, S_new)
                    t0 = time.time()
                    H_new = count_H_array(A_new)
                    t1 = time.time()
                    delta = H_new - H_int
                    swaps.append((a, b, H_new, delta, delta/H_int*100))
                    print(f"  Swap {a}→{b}: H={H_new}, ΔH={delta:+d} ({delta/H_int*100:+.4f}%) [{t1-t0:.1f}s]", flush=True)

swaps.sort(key=lambda x: -x[3])
print()
print(f"Best swap: {swaps[0][0]}→{swaps[0][1]} with ΔH = {swaps[0][3]:+d}")
all_neg = all(s[3] <= 0 for s in swaps)
print(f"ALL swaps negative: {all_neg}")
print(f"→ Interval is {'LOCAL MAX' if all_neg else 'NOT local max'} at p=19")
print()

# Also check Paley swaps
print("SINGLE-SWAP FROM PALEY:")
A_pal = make_tournament(p, QR)
H_pal = count_H_array(A_pal)
print(f"H(Paley) = {H_pal}")

pal_swaps = []
for a in QR:
    for b in range(1, p):
        if b not in QR:
            S_new = sorted([x for x in QR if x != a] + [b])
            if len(S_new) == m and len(set(S_new)) == m:
                valid = all((j in S_new) != ((p-j) in S_new) for j in range(1, m+1))
                if valid:
                    A_new = make_tournament(p, S_new)
                    H_new = count_H_array(A_new)
                    delta = H_new - H_pal
                    pal_swaps.append((a, b, H_new, delta))
                    print(f"  Swap {a}→{b}: ΔH={delta:+d} ({delta/H_pal*100:+.4f}%)", flush=True)

pal_swaps.sort(key=lambda x: -x[3])
pal_local = all(s[3] <= 0 for s in pal_swaps)
print(f"Paley local max: {pal_local}")
if not pal_local:
    print(f"Best swap from Paley: {pal_swaps[0][0]}→{pal_swaps[0][1]}, ΔH={pal_swaps[0][3]:+d}")

print()
print(f"SUMMARY p=19:")
print(f"  H(Interval) = {H_int}")
print(f"  H(Paley) = {H_pal}")
print(f"  Gap: {H_int - H_pal:+d} ({(H_int-H_pal)/H_pal*100:+.3f}%)")
print(f"  Interval local max: {all_neg}")
print(f"  Paley local max: {pal_local}")
print("\nDONE.")
