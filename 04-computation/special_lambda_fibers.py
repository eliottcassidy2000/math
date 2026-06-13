#!/usr/bin/env python3
"""
special_lambda_fibers.py — opus-2026-03-13-S71c

The 31 ambiguous lambda fibers at n=7 cluster into just two H-value types:
{109,111} with (c3,c5)=(11,23) and {141,143} with (c3,c5)=(12,30).

QUESTIONS:
1. What do these special lambda graphs look like?
2. Are they related by any symmetry?
3. How many distinct labeled lambda patterns produce each type?
4. What is the score sequence of tournaments in these fibers?
5. Is there a simple characterization of "c7-ambiguous" lambda graphs?

Also: these are exactly the Vitali atom fibers — where the (1,1,2,2)
reversal changes c7 by ±1 while preserving lambda.
"""

import sys, time
import numpy as np
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                pm = mask ^ (1 << v)
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def lambda_graph(A, n):
    L = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            for w in range(n):
                if w == u or w == v: continue
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    L[u][v] += 1; L[v][u] += 1
    return L

def lambda_key(L, n):
    return tuple(L[i][j] for i in range(n) for j in range(i+1, n))

n = 7
tb = n*(n-1)//2
np.random.seed(42)

# Find examples from each of the two types
type1_examples = []  # H ∈ {109, 111}
type2_examples = []  # H ∈ {141, 143}
seen_lk = defaultdict(list)

print("Searching for examples of the two special fiber types...")
t0 = time.time()

for trial in range(200000):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    lk = lambda_key(L, n)

    seen_lk[lk].append(bits)

    if len(seen_lk[lk]) >= 2:
        H_vals = set()
        for b in seen_lk[lk]:
            A2 = bits_to_adj(b, n)
            H_vals.add(count_ham_paths(A2, n))

        if len(H_vals) > 1:
            H_sorted = sorted(H_vals)
            if H_sorted == [109, 111]:
                if len(type1_examples) < 5 and lk not in [e[0] for e in type1_examples]:
                    type1_examples.append((lk, list(seen_lk[lk])))
            elif H_sorted == [141, 143]:
                if len(type2_examples) < 5 and lk not in [e[0] for e in type2_examples]:
                    type2_examples.append((lk, list(seen_lk[lk])))

    if trial % 50000 == 0:
        dt = time.time() - t0
        print(f"  trial {trial}: {dt:.1f}s, type1={len(type1_examples)}, type2={len(type2_examples)}")

    if len(type1_examples) >= 3 and len(type2_examples) >= 3:
        break

dt = time.time() - t0
print(f"\nFound: type1={len(type1_examples)}, type2={len(type2_examples)}, {dt:.1f}s")

# Analyze each type
for label, examples in [("Type 1: H∈{109,111}, (c3,c5)=(11,23)", type1_examples),
                         ("Type 2: H∈{141,143}, (c3,c5)=(12,30)", type2_examples)]:
    print(f"\n{'='*70}")
    print(label)
    print(f"{'='*70}")

    if not examples:
        print("  No examples found")
        continue

    for idx, (lk, bits_list) in enumerate(examples[:2]):
        # Analyze the lambda graph
        A = bits_to_adj(bits_list[0], n)
        L = lambda_graph(A, n)
        scores = [int(sum(A[i])) for i in range(n)]

        print(f"\n  Example {idx+1}:")
        print(f"  Scores: {sorted(scores)}")
        print(f"  Lambda matrix:")
        for i in range(n):
            row = [L[i][j] for j in range(n)]
            print(f"    {row}")

        # Lambda degree sequence (sum of lambda values per vertex)
        lam_deg = [sum(L[i][j] for j in range(n) if j != i) for i in range(n)]
        print(f"  Lambda degrees: {sorted(lam_deg)}")

        # Lambda value histogram
        lam_vals = []
        for i in range(n):
            for j in range(i+1, n):
                lam_vals.append(L[i][j])
        lam_hist = defaultdict(int)
        for v in lam_vals:
            lam_hist[v] += 1
        print(f"  Lambda value histogram: {dict(sorted(lam_hist.items()))}")

        # Find tournaments with different H values in this fiber
        H_to_bits = defaultdict(list)
        for b in bits_list:
            A2 = bits_to_adj(b, n)
            H2 = count_ham_paths(A2, n)
            H_to_bits[H2].append(b)

        for H_val, bl in sorted(H_to_bits.items()):
            A2 = bits_to_adj(bl[0], n)
            c3 = int(np.trace(A2 @ A2 @ A2)) // 3
            c5 = int(np.trace(np.linalg.matrix_power(A2, 5))) // 5
            scores2 = sorted([int(sum(A2[i])) for i in range(n)])
            print(f"    H={H_val}: c3={c3}, c5={c5}, scores={scores2}")

# =====================================================================
print(f"\n{'='*70}")
print("SCORE SEQUENCES IN AMBIGUOUS FIBERS")
print("=" * 70)

# Collect score sequences for all ambiguous fibers
np.random.seed(42)
lam_to_data = defaultdict(list)
for trial in range(100000):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    lk = lambda_key(L, n)
    H = count_ham_paths(A, n)
    scores = tuple(sorted([int(sum(A[i])) for i in range(n)]))
    c3 = int(np.trace(A @ A @ A)) // 3
    lam_to_data[lk].append((H, scores, c3))

ambig_scores = defaultdict(int)
ambig_c3 = defaultdict(int)
total_ambig = 0
for lk, entries in lam_to_data.items():
    H_vals = set(e[0] for e in entries)
    if len(H_vals) > 1:
        total_ambig += 1
        for _, scores, c3 in entries:
            ambig_scores[scores] += 1
            ambig_c3[c3] += 1

print(f"Ambiguous lambda fibers found: {total_ambig}")
print(f"Score sequences in ambiguous fibers: {dict(sorted(ambig_scores.items()))}")
print(f"c3 values in ambiguous fibers: {dict(sorted(ambig_c3.items()))}")

# =====================================================================
print(f"\n{'='*70}")
print("ARE BOTH H VALUES ALWAYS ACHIEVED? (SATURATION)")
print("=" * 70)

# For each ambiguous fiber, what fraction of tournaments have each H value?
h_fraction = defaultdict(lambda: defaultdict(int))
for lk, entries in lam_to_data.items():
    H_vals = set(e[0] for e in entries)
    if len(H_vals) > 1:
        for H, _, _ in entries:
            h_fraction[sorted(H_vals)[0]][H] += 1

print("\nH-value distribution within ambiguous fibers:")
for H_low, dist in sorted(h_fraction.items()):
    total = sum(dist.values())
    print(f"  H∈{{{H_low},{H_low+2}}}: ", end="")
    for H_val, count in sorted(dist.items()):
        print(f"H={H_val}: {count}/{total} ({100*count/total:.1f}%), ", end="")
    print()

# =====================================================================
print(f"\n{'='*70}")
print("CONNECTION TO VITALI ATOM (1,1,2,2) REVERSAL")
print("=" * 70)

# The Vitali atom is the unique lambda-preserving reversal that can change c7.
# It reverses the edge between vertices with equal lambda-degree that differ
# in sigma.

print("\nThe 31 ambiguous fibers all correspond to Vitali atom pairs.")
print("THM-182: The atom is the UNIQUE lambda-preserving reversal that changes c7.")
print("Since c7 changes by ±1, H changes by ±2.")
print()

# Check: are the H=109 and H=111 tournaments always related by a single arc reversal?
if type1_examples:
    lk, bits_list = type1_examples[0]
    H_to_bits = defaultdict(list)
    for b in bits_list:
        A2 = bits_to_adj(b, n)
        H2 = count_ham_paths(A2, n)
        H_to_bits[H2].append(b)

    if 109 in H_to_bits and 111 in H_to_bits:
        A_low = bits_to_adj(H_to_bits[109][0], n)
        A_high = bits_to_adj(H_to_bits[111][0], n)

        # Count differing edges
        diff_count = 0
        diffs = []
        for i in range(n):
            for j in range(n):
                if i != j and A_low[i][j] != A_high[i][j]:
                    diff_count += 1
                    diffs.append((i, j))

        print(f"  Type 1 example: H=109 vs H=111")
        print(f"  Number of differing entries: {diff_count}")
        print(f"  Differing edges: {diffs[:10]}")

        # A single arc reversal changes 2 entries: A[i][j] and A[j][i]
        print(f"  This is {diff_count//2} arc reversal(s)")

print("\nDone.")
