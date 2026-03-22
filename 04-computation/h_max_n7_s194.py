#!/usr/bin/env python3
"""
h_max_n7_s194.py — Find exact H_max at n=7
opus-2026-03-22-S194

Exhaustive search of all 2^21 = 2,097,152 tournaments at n=7
to find the true H_max. Also verify the OEIS sequence for H_max.
"""

from math import comb
import time

print("=" * 72)
print("  H_max AT n=7: EXHAUSTIVE SEARCH")
print("  opus-2026-03-22-S194")
print("=" * 72)
print()

def tournament_adj(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): adj[i][j] = 1
            else: adj[j][i] = 1
            idx += 1
    return adj

def H_dp(adj, n):
    dp = [0]*((1<<n)*n)
    for v in range(n): dp[(1<<v)*n+v] = 1
    for S in range(1, 1<<n):
        for v in range(n):
            if not (S&(1<<v)): continue
            val = dp[S*n+v]
            if val == 0: continue
            for u in range(n):
                if S&(1<<u): continue
                if adj[v][u]: dp[(S|(1<<u))*n+u] += val
    full = (1<<n)-1
    return sum(dp[full*n+v] for v in range(n))

def score_seq(adj, n):
    return tuple(sorted(sum(adj[i][j] for j in range(n) if j != i) for i in range(n)))

n = 7
m = comb(n, 2)  # 21
total = 1 << m  # 2,097,152

print(f"  n={n}, m={m}, total tournaments = {total:,}")
print(f"  Searching for H_max...")
print()

t0 = time.time()
best_h = 0
best_bits = []
count_checked = 0
milestone = total // 20

for bits in range(total):
    adj = tournament_adj(n, bits)
    h = H_dp(adj, n)

    if h > best_h:
        best_h = h
        best_bits = [bits]
    elif h == best_h:
        best_bits.append(bits)

    count_checked += 1
    if count_checked % milestone == 0:
        elapsed = time.time() - t0
        pct = 100 * count_checked / total
        rate = count_checked / elapsed
        remaining = (total - count_checked) / rate
        print(f"  {pct:.0f}% done ({count_checked:,}/{total:,}), "
              f"current best H = {best_h}, "
              f"elapsed = {elapsed:.0f}s, ETA = {remaining:.0f}s")

elapsed = time.time() - t0
print()
print(f"  DONE in {elapsed:.1f}s ({total/elapsed:.0f} tournaments/sec)")
print()
print(f"  H_max(7) = {best_h}")
print(f"  Achieved by {len(best_bits)} labeled tournaments")
print()

# Show some achievers
for bits in best_bits[:10]:
    adj = tournament_adj(n, bits)
    ss = score_seq(adj, n)
    print(f"    bits={bits}, score = {ss}")

# Count by score sequence
from collections import Counter
score_counter = Counter()
for bits in best_bits:
    adj = tournament_adj(n, bits)
    ss = score_seq(adj, n)
    score_counter[ss] += 1

print()
print(f"  Score sequences achieving H_max(7) = {best_h}:")
for ss, count in score_counter.most_common():
    print(f"    {ss}: {count} tournaments")

print()

# H_max sequence update
print(f"  H_max sequence (n=3..7): 3, 5, 15, 45, {best_h}")

# Check OEIS A000438 (which is the number of H-paths in complete graph Kn, different thing)
# Actually what OEIS sequence is this?
print(f"  Looking for OEIS match: 3, 5, 15, 45, {best_h}")
print()

# Is C_7 achieving H_max?
adj_c7 = [[0]*7 for _ in range(7)]
for i in range(7):
    for j in range(7):
        if i == j: continue
        if 0 < (j - i) % 7 <= 3:
            adj_c7[i][j] = 1
h_c7 = H_dp(adj_c7, 7)
print(f"  H(C_7) = {h_c7}")
if h_c7 == best_h:
    print(f"  C_7 achieves H_max(7) ✓")
else:
    print(f"  C_7 does NOT achieve H_max(7)! H(C_7) = {h_c7} < {best_h}")

print()
print("opus-2026-03-22-S194")
