#!/usr/bin/env python3
"""Fast H=21 search using optimized DP."""
import random, time, sys
from collections import Counter

def H(adj_bits, n):
    """Count Hamiltonian paths. adj_bits[v] = bitmask of v's out-neighbors."""
    full = (1 << n) - 1
    # dp[mask][v] = # paths through mask ending at v
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1

    for mask in range(1, 1 << n):
        for v in range(n):
            c = dp.get((mask, v), 0)
            if c == 0:
                continue
            out = adj_bits[v] & ~mask  # out-neighbors not yet visited
            while out:
                u = (out & -out).bit_length() - 1
                key = (mask | (1 << u), u)
                dp[key] = dp.get(key, 0) + c
                out &= out - 1

    return sum(dp.get((full, v), 0) for v in range(n))


def bits_to_adj(n, bits):
    """Convert bit encoding to adjacency bitmasks."""
    adj = [0] * n
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << k):
                adj[i] |= (1 << j)
            else:
                adj[j] |= (1 << i)
            k += 1
    return adj


def rand_adj(n):
    """Random tournament as adjacency bitmasks."""
    adj = [0] * n
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                adj[i] |= (1 << j)
            else:
                adj[j] |= (1 << i)
    return adj


# ============================================================
# Exhaustive n <= 6
# ============================================================
print("=" * 60)
print("Exhaustive search n=3..6")
print("=" * 60)

all_h_values = set()

for n in range(3, 7):
    m = n*(n-1)//2
    total = 1 << m
    hc = Counter()
    t0 = time.time()

    for bits in range(total):
        adj = bits_to_adj(n, bits)
        h = H(adj, n)
        hc[h] += 1

    elapsed = time.time() - t0
    all_h = sorted(hc.keys())
    all_h_values.update(hc.keys())
    missing = [h for h in range(1, 52, 2) if h not in hc and h <= max(hc)]
    print(f"n={n} ({total} tours, {elapsed:.1f}s): H values = {all_h}")
    print(f"  Missing odd <= {max(hc)}: {missing}")
    print(f"  Near 21: ", end="")
    for hv in [17,19,21,23,25]:
        if hv in hc: print(f"H={hv}:{hc[hv]} ", end="")
    print()
    sys.stdout.flush()

# ============================================================
# Exhaustive n=7 (2M tournaments -- may take ~10 min in Python)
# Let's try with progress reporting
# ============================================================
print(f"\n{'='*60}")
print("Exhaustive search n=7 (2,097,152 tournaments)")
print("=" * 60)

n = 7
m = 21
total = 1 << m
hc = Counter()
t0 = time.time()

for bits in range(total):
    adj = bits_to_adj(n, bits)
    h = H(adj, n)
    hc[h] += 1
    if (bits + 1) % 200000 == 0:
        elapsed = time.time() - t0
        rate = (bits+1)/elapsed
        eta = (total - bits - 1)/rate
        print(f"  {bits+1}/{total} ({elapsed:.0f}s, {rate:.0f}/s, ETA {eta:.0f}s)", flush=True)

elapsed = time.time() - t0
all_h = sorted(hc.keys())
all_h_values.update(hc.keys())
missing = [h for h in range(1, 52, 2) if h not in hc and h <= max(hc)]
print(f"n=7 ({total} tours, {elapsed:.1f}s): all odd H = {[h for h in all_h if h%2==1]}")
print(f"  Missing odd <= 51: {missing}")
print(f"  Near 21: ", end="")
for hv in [15,17,19,21,23,25,27]:
    if hv in hc: print(f"H={hv}:{hc[hv]} ", end="")
print(flush=True)

# ============================================================
# Random sampling n=8, 9, 10
# ============================================================
for n, samples in [(8, 200000), (9, 100000), (10, 20000)]:
    print(f"\n{'='*60}")
    print(f"Random search n={n}: {samples} samples")
    print("=" * 60)

    hc = Counter()
    t0 = time.time()
    found_21 = None

    for trial in range(samples):
        adj = rand_adj(n)
        h = H(adj, n)
        hc[h] += 1
        if h == 21 and found_21 is None:
            found_21 = adj[:]
            print(f"  *** H=21 FOUND at trial {trial}! adj_bits={adj} ***")
        if (trial+1) % 50000 == 0:
            elapsed = time.time() - t0
            print(f"  {trial+1}/{samples} ({elapsed:.0f}s)", flush=True)

    elapsed = time.time() - t0
    all_h_values.update(hc.keys())

    print(f"  Time: {elapsed:.1f}s")
    print(f"  Near 21: ", end="")
    for hv in [15,17,19,21,23,25,27]:
        if hv in hc: print(f"H={hv}:{hc[hv]} ", end="")
    print()

    missing = [h for h in range(1, 52, 2) if h not in hc and h <= 51]
    print(f"  Missing odd <= 51: {missing}")

    if found_21:
        print(f"  H=21 FOUND! adj_bits = {found_21}")
    else:
        print(f"  H=21 NOT found")
    sys.stdout.flush()

# ============================================================
# Summary
# ============================================================
print(f"\n{'='*60}")
print("FINAL SUMMARY")
print("=" * 60)
missing_global = sorted(h for h in range(1, 52, 2) if h not in all_h_values)
print(f"Missing odd H values <= 51 across all searches: {missing_global}")
if 21 in all_h_values:
    print("CONCLUSION: H=21 IS achievable!")
else:
    print("CONCLUSION: H=21 NOT found. Evidence strongly suggests it is impossible.")
