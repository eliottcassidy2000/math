#!/usr/bin/env python3
"""
selfwt_H_s20cn.py -- kind-pasteur-2026-03-22-S20cn

Compute self_wt_7 using H (Hamiltonian path count) as discriminant.

For n<=7: (score, c3, H) almost certainly determines the iso class.
H(T) = H(T') iff T ~ T' for "most" pairs (with very rare exceptions).

We only compute H for arc flips that pass the (score) check.
H computation for n=7 is O(7 * 2^7) = 896 ops. Feasible for the subset.

Author: kind-pasteur-2026-03-22-S20cn
"""
import sys
from math import comb
import time
sys.stdout.reconfigure(line_buffering=True)

n = 7
m = comb(n, 2)  # 21
total = 1 << m

PAIRS = [(i,j) for i in range(n) for j in range(i+1,n)]
pair_idx = {}
for k, (p,q) in enumerate(PAIRS):
    pair_idx[(p,q)] = k

def compute_scores(bits):
    s = [0] * n
    for k, (i,j) in enumerate(PAIRS):
        if bits & (1 << k): s[i] += 1
        else: s[j] += 1
    return s

def H_count(bits):
    """Count Hamiltonian paths via bitmask DP."""
    dp = [0] * ((1 << n) * n)
    for v in range(n):
        dp[(1 << v) * n + v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)): continue
            val = dp[S * n + v]
            if val == 0: continue
            for u in range(n):
                if S & (1 << u): continue
                # Check if v->u in tournament
                if v < u:
                    has_edge = bits & (1 << pair_idx[(v,u)])
                else:
                    has_edge = not (bits & (1 << pair_idx[(u,v)]))
                if has_edge:
                    dp[(S | (1 << u)) * n + u] += val
    full = (1 << n) - 1
    return sum(dp[full * n + v] for v in range(n))

print("=" * 70)
print("  self_wt via H discriminant at n=7")
print("=" * 70)

t0 = time.time()
total_self_H = 0  # self-flips by (score, H) match
total_twins_check = 0  # just for verification

for bits in range(total):
    scores = compute_scores(bits)
    sorted_score = tuple(sorted(scores))

    # Compute H only once per tournament
    H_val = None

    for k in range(m):
        i, j = PAIRS[k]

        # Quick score check first
        if bits & (1 << k):
            ns_i, ns_j = scores[i] - 1, scores[j] + 1
        else:
            ns_i, ns_j = scores[i] + 1, scores[j] - 1

        new_s = list(scores)
        new_s[i] = ns_i
        new_s[j] = ns_j
        if tuple(sorted(new_s)) != sorted_score:
            continue

        # Score matches! Compute H for both tournaments
        if H_val is None:
            H_val = H_count(bits)

        flipped = bits ^ (1 << k)
        H_flip = H_count(flipped)

        if H_val == H_flip:
            total_self_H += 1

    if (bits + 1) % 200000 == 0:
        elapsed = time.time() - t0
        rate = (bits + 1) / elapsed
        eta = (total - bits - 1) / rate
        print(f"    {bits+1}/{total} ({rate:.0f}/s, ETA {eta:.0f}s) self_H={total_self_H}")

elapsed = time.time() - t0
print(f"\n  Done in {elapsed:.1f}s")

print(f"\n  RESULTS:")
print(f"  Total (score, H)-matched self-flips: {total_self_H}")
print(f"  This is an UPPER BOUND on self_wt (H not always unique)")
print(f"  Predicted self_wt: 1527456")
print(f"  Predicted NTE(7): 151200 = self_wt - 1376256")
print(f"  NTE from (score,H): {total_self_H - 1376256}")
