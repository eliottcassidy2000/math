#!/usr/bin/env python3
"""
selfwt_sample_n7_s20cn.py -- kind-pasteur-2026-03-22-S20cn

Verify self_wt_7 by sampling random tournaments.

Predicted: self_wt_7 = 1527456
         = avg self-flips per tournament * 2^21
         = 0.7285 self-flips per tournament (on average)

We sample 50000 random tournaments and count self-flip arcs.

A "self-flip" at arc e means T xor e is isomorphic to T.
We use (score, c3, deletion_H) fingerprint for fast iso check.

Author: kind-pasteur-2026-03-22-S20cn
"""
import sys
import random
from math import comb
from collections import Counter
sys.stdout.reconfigure(line_buffering=True)

n = 7
m = comb(n, 2)  # 21
PAIRS = [(i,j) for i in range(n) for j in range(i+1,n)]

def bits_to_adj(bits):
    adj = [[0]*n for _ in range(n)]
    for k, (i,j) in enumerate(PAIRS):
        if bits & (1 << k): adj[i][j] = 1
        else: adj[j][i] = 1
    return adj

def H_dp(adj, nn):
    dp = [0] * ((1 << nn) * nn)
    for v in range(nn):
        dp[(1 << v) * nn + v] = 1
    for S in range(1, 1 << nn):
        for v in range(nn):
            if not (S & (1 << v)): continue
            val = dp[S * nn + v]
            if val == 0: continue
            for u in range(nn):
                if S & (1 << u): continue
                if adj[v][u]:
                    dp[(S | (1 << u)) * nn + u] += val
    full = (1 << nn) - 1
    return sum(dp[full * nn + v] for v in range(nn))

def fingerprint(adj):
    """(score_seq, c3, deletion_H_fingerprint)"""
    scores = tuple(sorted(sum(adj[i][j] for j in range(n)) for i in range(n)))
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]: c3 += 1
                if adj[i][k] and adj[k][j] and adj[j][i]: c3 += 1
    del_h = []
    for v in range(n):
        verts = [i for i in range(n) if i != v]
        nn = len(verts)
        sub = [[adj[verts[i]][verts[j]] for j in range(nn)] for i in range(nn)]
        del_h.append(H_dp(sub, nn))
    return (scores, c3, tuple(sorted(del_h)))

print("=" * 70)
print("  SELF-FLIP SAMPLING AT n=7")
print("=" * 70)

N_SAMPLES = 50000
total_self_flips = 0
twin_flips = 0
nontwin_flips = 0
flip_histogram = Counter()

random.seed(42)

for trial in range(N_SAMPLES):
    bits = random.randint(0, (1 << m) - 1)
    adj = bits_to_adj(bits)
    fp_T = fingerprint(adj)

    self_count = 0
    twin_count = 0
    for k in range(m):
        flipped = bits ^ (1 << k)
        adj_flip = bits_to_adj(flipped)
        fp_flip = fingerprint(adj_flip)
        if fp_T == fp_flip:
            self_count += 1
            # Check if twin pair
            i, j = PAIRS[k]
            is_twin = all(adj[a][i] == adj[a][j] for a in range(n) if a != i and a != j)
            if is_twin:
                twin_count += 1

    total_self_flips += self_count
    twin_flips += twin_count
    nontwin_flips += (self_count - twin_count)
    flip_histogram[self_count] += 1

    if (trial + 1) % 10000 == 0:
        avg = total_self_flips / (trial + 1)
        print(f"  {trial+1}/{N_SAMPLES}: avg self-flips = {avg:.4f}, twin = {twin_flips/(trial+1):.4f}, non-twin = {nontwin_flips/(trial+1):.4f}")

avg_self_flips = total_self_flips / N_SAMPLES
avg_twin = twin_flips / N_SAMPLES
avg_nontwin = nontwin_flips / N_SAMPLES

print(f"\n  RESULTS ({N_SAMPLES} samples):")
print(f"  avg self-flips per tournament: {avg_self_flips:.4f}")
print(f"  avg twin-flips: {avg_twin:.4f}")
print(f"  avg non-twin flips: {avg_nontwin:.4f}")

# Extrapolate
est_self_wt = avg_self_flips * (1 << m)
est_twins = avg_twin * (1 << m)
est_nte = avg_nontwin * (1 << m)

print(f"\n  Estimated self_wt_7: {est_self_wt:.0f}")
print(f"  Predicted self_wt_7: 1527456")
print(f"  Ratio: {est_self_wt / 1527456:.4f}")

print(f"\n  Estimated twins: {est_twins:.0f}")
print(f"  Predicted twins: 1376256")
print(f"  Ratio: {est_twins / 1376256:.4f}")

print(f"\n  Estimated NTE: {est_nte:.0f}")
print(f"  Predicted NTE(7): 151200")
print(f"  Ratio: {est_nte / 151200:.4f}")

# Self-flip distribution
print(f"\n  Self-flip distribution:")
for k in sorted(flip_histogram.keys()):
    pct = 100 * flip_histogram[k] / N_SAMPLES
    print(f"    {k} self-flips: {flip_histogram[k]} ({pct:.1f}%)")
