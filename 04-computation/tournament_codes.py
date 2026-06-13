"""
tournament_codes.py -- kind-pasteur-2026-03-14-S74
Tournaments as error-correcting codes.

IDEA: Each tournament T on n vertices is a binary string of length m = C(n,2).
The set of tournaments with H(T) = k forms a "code" C_k of length m.
What are the error-correcting properties of these codes?

HAMMING DISTANCE: d(T, T') = #{arcs that differ} = Hamming distance in {0,1}^m.

QUESTIONS:
1. What is the minimum distance of C_k (the H-constant code)?
2. For the H-maximizers (C_max), is the minimum distance large?
3. Do the H-fiber codes have algebraic structure?

ALSO: The ARC-FLIP GRAPH is the Hamming graph on {0,1}^m restricted by H-value.
The "decoding" problem: given a tournament, find the closest H-maximizer.

AND: Connection to the Paley conference matrix — Paley tournaments have
connections to Hadamard codes.
"""

import numpy as np
from itertools import combinations
from collections import Counter, defaultdict
import sys

sys.stdout.reconfigure(encoding='utf-8')

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

def compute_H_dp(A, n):
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

def hamming_dist(a, b):
    return bin(a ^ b).count('1')

def main():
    print("=" * 70)
    print("TOURNAMENTS AS ERROR-CORRECTING CODES")
    print("kind-pasteur-2026-03-14-S74")
    print("=" * 70)

    for n in [4, 5, 6]:
        m = n * (n - 1) // 2
        N = 2**m
        print(f"\n{'='*70}")
        print(f"n = {n}, m = {m} (code length), N = {N} tournaments")
        print(f"{'='*70}")

        # Compute all H values
        H_map = {}
        count = 0
        for bits in range(N):
            count += 1
            if n >= 6 and count > 20000:
                break
            A = bits_to_adj(bits, n)
            H_map[bits] = compute_H_dp(A, n)

        # Group by H
        H_fibers = defaultdict(list)
        for bits, H in H_map.items():
            H_fibers[H].append(bits)

        # For each H-fiber, compute code parameters
        print(f"\n  H-FIBER CODE PARAMETERS [n, k, d]:")
        print(f"  {'H':>5} {'|C|':>8} {'min_d':>7} {'max_d':>7} {'avg_d':>8} {'rate':>8}")

        for H in sorted(H_fibers.keys()):
            codewords = H_fibers[H]
            size = len(codewords)

            if size <= 1:
                min_d = m
                max_d = 0
                avg_d = 0
            elif size <= 500:
                # Compute all pairwise distances
                dists = []
                for i in range(len(codewords)):
                    for j in range(i+1, len(codewords)):
                        dists.append(hamming_dist(codewords[i], codewords[j]))
                min_d = min(dists)
                max_d = max(dists)
                avg_d = np.mean(dists)
            else:
                # Sample pairwise distances
                import random
                random.seed(42)
                dists = []
                for _ in range(5000):
                    i, j = random.sample(range(len(codewords)), 2)
                    dists.append(hamming_dist(codewords[i], codewords[j]))
                min_d = min(dists)
                max_d = max(dists)
                avg_d = np.mean(dists)

            rate = np.log2(size) / m if size > 1 else 0
            print(f"  {H:5d} {size:8d} {min_d:7d} {max_d:7d} {avg_d:8.2f} {rate:8.4f}")

        # SPECIAL: H-maximizer code
        max_H = max(H_fibers.keys())
        max_cw = H_fibers[max_H]
        print(f"\n  H-MAXIMIZER CODE (H={max_H}):")
        print(f"    Codewords: {len(max_cw)}")
        if len(max_cw) <= 500:
            dists = [hamming_dist(max_cw[i], max_cw[j])
                     for i in range(len(max_cw)) for j in range(i+1, len(max_cw))]
            if dists:
                print(f"    Min distance: {min(dists)}")
                print(f"    Distance distribution: {dict(Counter(dists))}")

        # SPECIAL: Transitive tournament code (H=1)
        trans_cw = H_fibers.get(1, [])
        print(f"\n  TRANSITIVE CODE (H=1):")
        print(f"    Codewords: {len(trans_cw)}")
        if len(trans_cw) <= 500:
            dists = [hamming_dist(trans_cw[i], trans_cw[j])
                     for i in range(len(trans_cw)) for j in range(i+1, len(trans_cw))]
            if dists:
                dist_counter = Counter(dists)
                print(f"    Min distance: {min(dists)}")
                print(f"    Max distance: {max(dists)}")
                print(f"    Distance distribution: {dict(sorted(dist_counter.items()))}")

        # COVERING RADIUS: max distance from any tournament to nearest maximizer
        if len(max_cw) > 0 and n <= 5:
            max_cover = 0
            for bits in range(N):
                min_dist_to_max = min(hamming_dist(bits, cw) for cw in max_cw)
                max_cover = max(max_cover, min_dist_to_max)
            print(f"\n  COVERING RADIUS (max dist to nearest H-maximizer): {max_cover}")
            print(f"    This means: from ANY tournament, at most {max_cover} arc flips to reach H={max_H}")

        # WEIGHT DISTRIBUTION of maximizer code
        if max_cw and n <= 5:
            # Weight of each codeword = Hamming weight (number of 1-bits)
            weights = [bin(cw).count('1') for cw in max_cw]
            print(f"\n  WEIGHT DISTRIBUTION of H={max_H} code:")
            for w, cnt in sorted(Counter(weights).items()):
                print(f"    weight {w}: {cnt} codewords")

    print(f"\n{'='*70}")
    print("DONE")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
