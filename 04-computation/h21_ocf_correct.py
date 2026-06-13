#!/usr/bin/env python3
"""
Correct OCF analysis for H=21.
Uses DIRECTED CYCLES as nodes in Omega (not vertex sets).

H(T) = I(Omega(T), 2) = 1 + 2*i_1 + 4*i_2 + ...
where i_k = # independent sets of size k in Omega.

For H=21: need 2*i_1 + 4*i_2 + 8*i_3 + ... = 20.
"""
import sys, time, random
from itertools import combinations, permutations
from collections import Counter, defaultdict

def H_matrix(A, n):
    full = (1 << n) - 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            c = dp[mask][v]
            if not (mask & (1 << v)) or c == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += c
    return sum(dp[full])


def find_directed_odd_cycles(A, n):
    """Find all directed odd cycles (canonical: min vertex first)."""
    cycles = []
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            first = verts[0]
            for perm in permutations(verts[1:]):
                path = (first,) + perm
                valid = True
                for i in range(length):
                    if not A[path[i]][path[(i+1) % length]]:
                        valid = False
                        break
                if valid:
                    cycles.append(path)
    return cycles


def ocf_ivec(A, n):
    """Compute (i_0, i_1, i_2, ...) using directed cycles."""
    cycles = find_directed_odd_cycles(A, n)
    m = len(cycles)
    if m == 0:
        return (1,), 1

    vsets = [frozenset(c) for c in cycles]
    adj_bits = [0] * m
    for a in range(m):
        for b in range(a+1, m):
            if vsets[a] & vsets[b]:
                adj_bits[a] |= 1 << b
                adj_bits[b] |= 1 << a

    if m > 22:
        # Too large for enumeration, use size-bounded approach
        return None, None

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
# Exhaustive n=5: all i-vectors
# ============================================================
print("=" * 60)
print("Exhaustive OCF i-vectors at n=5 (1024 tournaments)")
print("=" * 60)

n = 5
m_bits = 10
h_to_ivecs = defaultdict(set)
t0 = time.time()

for bits in range(1 << m_bits):
    A = [[0]*n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << k): A[i][j] = 1
            else: A[j][i] = 1
            k += 1
    ivec, I2 = ocf_ivec(A, n)
    h_to_ivecs[I2].add(ivec)

for hv in sorted(h_to_ivecs.keys()):
    print(f"  H={hv}: i-vecs = {sorted(h_to_ivecs[hv])}")

# ============================================================
# Exhaustive n=6: i-vectors (this may be slow due to cycle enumeration)
# ============================================================
print(f"\n{'='*60}")
print("Exhaustive OCF i-vectors at n=6 (32768 tournaments)")
print("=" * 60)

n = 6
m_bits = 15
h_to_ivecs = defaultdict(set)
all_ivecs = set()
t0 = time.time()

for bits in range(1 << m_bits):
    A = [[0]*n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << k): A[i][j] = 1
            else: A[j][i] = 1
            k += 1
    ivec, I2 = ocf_ivec(A, n)
    if ivec is not None:
        h_to_ivecs[I2].add(ivec)
        all_ivecs.add(ivec)
    if (bits + 1) % 5000 == 0:
        elapsed = time.time() - t0
        print(f"  {bits+1}/32768 ({elapsed:.0f}s)", flush=True)

elapsed = time.time() - t0
print(f"Total time: {elapsed:.0f}s")

print(f"\nAll H values: {sorted(h_to_ivecs.keys())}")
for hv in sorted(h_to_ivecs.keys()):
    if hv <= 50:
        print(f"  H={hv}: i-vecs = {sorted(h_to_ivecs[hv])}")

# What i_1 values (number of directed cycles) are achieved?
i1_vals = sorted(set(ivec[1] if len(ivec) > 1 else 0 for ivec in all_ivecs))
print(f"\nAchieved i_1 (# directed odd cycles) at n=6: {i1_vals}")

# Missing H values
all_h = sorted(h_to_ivecs.keys())
missing = [h for h in range(1, max(all_h)+1, 2) if h not in h_to_ivecs]
print(f"Missing odd H values: {missing}")

sys.stdout.flush()
