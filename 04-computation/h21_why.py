#!/usr/bin/env python3
"""
Analyze WHY H=21 is impossible.

H(T) = I(Omega(T), 2) = 1 + 2*i_1 + 4*i_2 + 8*i_3 + ...

For H=21: need 2*i_1 + 4*i_2 + 8*i_3 + ... = 20

Key approach: Enumerate all achievable (i_1, i_2, i_3) vectors at n=6 and n=7,
and show none gives H=21.
"""

import sys, time, random
from itertools import combinations, permutations
from collections import Counter, defaultdict


def H_adj(adj, n):
    """Count Hamiltonian paths from adjacency bitmasks."""
    full = (1 << n) - 1
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            c = dp.get((mask, v), 0)
            if c == 0: continue
            out = adj[v] & ~mask
            while out:
                u = (out & -out).bit_length() - 1
                key = (mask | (1 << u), u)
                dp[key] = dp.get(key, 0) + c
                out &= out - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def bits_to_adj(n, bits):
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


def bits_to_matrix(n, bits):
    A = [[0]*n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << k):
                A[i][j] = 1
            else:
                A[j][i] = 1
            k += 1
    return A


def find_odd_cycle_vsets(A, n):
    """Find all directed odd cycle vertex sets."""
    result = []
    for size in range(3, n+1, 2):
        for subset in combinations(range(n), size):
            s = size
            sub = list(subset)
            idx = {v: i for i, v in enumerate(sub)}
            sub_A = [[0]*s for _ in range(s)]
            for a in sub:
                for b in sub:
                    if a != b and A[a][b]:
                        sub_A[idx[a]][idx[b]] = 1
            dp = [[0]*s for _ in range(1 << s)]
            dp[1][0] = 1
            for mask in range(1, 1 << s):
                if not (mask & 1): continue
                for v in range(s):
                    c = dp[mask][v]
                    if not (mask & (1 << v)) or c == 0: continue
                    for u in range(s):
                        if mask & (1 << u): continue
                        if sub_A[v][u]:
                            dp[mask | (1 << u)][u] += c
            full = (1 << s) - 1
            if any(dp[full][v] > 0 and sub_A[v][0] for v in range(1, s)):
                result.append(frozenset(subset))
    return result


def ocf_vector(A, n):
    """Return (i_0, i_1, i_2, ...) and I(Omega,2)."""
    cycles = find_odd_cycle_vsets(A, n)
    m = len(cycles)
    if m == 0:
        return (1,), 1

    adj_bits = [0] * m
    for a in range(m):
        for b in range(a+1, m):
            if cycles[a] & cycles[b]:
                adj_bits[a] |= 1 << b
                adj_bits[b] |= 1 << a

    size_counts = Counter()
    for mask in range(1 << m):
        ok = True
        seen = 0
        temp = mask
        while temp:
            v = (temp & -temp).bit_length() - 1
            if adj_bits[v] & seen:
                ok = False
                break
            seen |= 1 << v
            temp &= temp - 1
        if ok:
            size_counts[bin(mask).count('1')] += 1

    max_k = max(size_counts.keys()) if size_counts else 0
    ivec = tuple(size_counts.get(k, 0) for k in range(max_k + 1))
    I2 = sum(ivec[k] * (2**k) for k in range(len(ivec)))
    return ivec, I2


# ============================================================
# Exhaustive n=6: catalog ALL achieved (i_1, i_2, ...) vectors
# ============================================================
print("=" * 60)
print("Exhaustive OCF vectors at n=6")
print("=" * 60)

n = 6
m_bits = n*(n-1)//2
total = 1 << m_bits

# Group by H value
h_to_ivecs = defaultdict(set)  # H -> set of i-vectors
ivec_to_h = {}

t0 = time.time()
for bits in range(total):
    A = bits_to_matrix(n, bits)
    ivec, I2 = ocf_vector(A, n)
    h_to_ivecs[I2].add(ivec)
    ivec_to_h[ivec] = I2

elapsed = time.time() - t0
print(f"Time: {elapsed:.1f}s")

# All achieved H
all_h = sorted(h_to_ivecs.keys())
print(f"All H values at n=6: {all_h}")

# Show i-vectors for H near 21
for hv in range(15, 28, 2):
    if hv in h_to_ivecs:
        print(f"\nH={hv}: i-vectors = {sorted(h_to_ivecs[hv])}")
    else:
        print(f"\nH={hv}: NOT ACHIEVABLE")

# What i-vectors COULD give H=21?
print(f"\n{'='*60}")
print("Theoretical i-vectors for H=21:")
print("=" * 60)

for i1 in range(20):
    r1 = 20 - 2*i1
    if r1 < 0: break
    for i2 in range(r1//4 + 1):
        r2 = r1 - 4*i2
        if r2 < 0: break
        for i3 in range(r2//8 + 1):
            r3 = r2 - 8*i3
            if r3 == 0:
                ivec = (1, i1, i2, i3) if i3 > 0 else ((1, i1, i2) if i2 > 0 else (1, i1))
                achievable = ivec in ivec_to_h
                print(f"  (i_0=1, i_1={i1}, i_2={i2}, i_3={i3}): {'ACHIEVABLE' if achievable else 'NOT achievable'} at n=6")

# ============================================================
# Look at all i-vectors achieved at n=6
# ============================================================
print(f"\n{'='*60}")
print("All i-vectors achieved at n=6 (sorted by H):")
print("=" * 60)
for hv in sorted(h_to_ivecs.keys()):
    if hv <= 50:
        for ivec in sorted(h_to_ivecs[hv]):
            print(f"  H={hv}: i-vec={ivec}")

# ============================================================
# Key question: what i_1 values are achievable at n=6?
# ============================================================
print(f"\n{'='*60}")
print("Achievable i_1 (=alpha_1 = # directed odd cycle vertex sets) at n=6:")
print("=" * 60)
i1_vals = set()
for ivec in ivec_to_h:
    if len(ivec) > 1:
        i1_vals.add(ivec[1])
    else:
        i1_vals.add(0)
print(f"i_1 values: {sorted(i1_vals)}")

# ============================================================
# H mod 4 analysis
# ============================================================
print(f"\n{'='*60}")
print("H mod 4 analysis at n=6:")
print("=" * 60)
for hv in sorted(h_to_ivecs.keys()):
    if hv <= 50:
        print(f"  H={hv}: H mod 4 = {hv % 4}")

# Count by H mod 4
mod4_counts = Counter()
for bits in range(total):
    adj = bits_to_adj(n, bits)
    h = H_adj(adj, n)
    mod4_counts[h % 4] += 1
print(f"\nH mod 4 distribution: {dict(sorted(mod4_counts.items()))}")

sys.stdout.flush()

# ============================================================
# Try n=7 with sampling to find i-vectors near H=21
# ============================================================
print(f"\n{'='*60}")
print("n=7: Sample i-vectors near H=21 (10K samples)")
print("=" * 60)

n = 7
near_21 = defaultdict(list)
t0 = time.time()

for trial in range(10000):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    adj = [0]*n
    for i in range(n):
        for j in range(n):
            if A[i][j]:
                adj[i] |= (1 << j)
    h = H_adj(adj, n)
    if 15 <= h <= 27 and len(near_21[h]) < 5:
        near_21[h].append(A)

elapsed = time.time() - t0
print(f"Time: {elapsed:.1f}s")

for hv in sorted(near_21.keys()):
    print(f"\nH={hv} ({len(near_21[hv])} examples):")
    for A in near_21[hv][:2]:
        ivec, I2 = ocf_vector(A, n)
        print(f"  i-vec={ivec}, verify I(Omega,2)={I2}")

if 21 not in near_21:
    print(f"\nH=21: NOT FOUND at n=7")

# ============================================================
# Binary structure analysis
# ============================================================
print(f"\n{'='*60}")
print("Binary analysis of H values")
print("=" * 60)

print("H=7  = 0b111 = 1+2+4 = i_0 + 2*i_1 + 4*i_2, need i_1=1, i_2=1 (3 cycles, 1 pair disjoint)")
print("H=21 = 0b10101 = 1+4+16 = i_0 + 4*i_2 + 16*i_4, need i_1=0, i_2=1, i_4=1")
print("   OR H=21 = 1+20, and 20 = 2*10 => i_1=10")
print("   OR 20 = 2*8 + 4*1 => i_1=8, i_2=1")
print("   etc.")
print()
print("The question: is there a parity/divisibility constraint blocking 20?")
print()
print("Note H mod 4:")
print("  H = 1 + 2*i_1 + 4*i_2 + 8*i_3 + ...")
print("  H mod 4 = (1 + 2*i_1) mod 4")
print("  If i_1 is even: H mod 4 = 1")
print("  If i_1 is odd:  H mod 4 = 3")
print()
print("H=21: 21 mod 4 = 1, so i_1 must be even.")
print("H=7:   7 mod 4 = 3, so i_1 must be odd.")
print("H=19: 19 mod 4 = 3, so i_1 must be odd.")
print("H=23: 23 mod 4 = 3, so i_1 must be odd.")
print()

# Check: at n=6, what are the H mod 4 = 1 values?
h_mod4_1 = sorted(h for h in h_to_ivecs if h % 4 == 1)
print(f"H values with H mod 4 = 1 at n=6: {h_mod4_1}")
h_mod4_1_le50 = [h for h in h_mod4_1 if h <= 50]
print(f"  Up to 50: {h_mod4_1_le50}")
missing_mod4_1 = [h for h in range(1, 50, 4) if h not in set(h_mod4_1)]
print(f"  Missing (H=1 mod 4, <= 50): {missing_mod4_1}")

sys.stdout.flush()
