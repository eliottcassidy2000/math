"""
e4_n6_verify.py -- kind-pasteur-2026-03-14-S105i
VERIFY level-4 Fourier coefficients at n=6 DIRECTLY.

At n=5: 60 nonzero, all magnitude 1/8, all cover 5 vertices.
At n=6: E_4 = 45/4. Need to check magnitudes and count.
"""

import sys, math
import numpy as np
from itertools import combinations

sys.stdout.reconfigure(encoding='utf-8')

def count_ham_paths(adj, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if (mask, v) not in dp:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def main():
    n = 6
    m = n*(n-1)//2  # 15
    N = 1 << m  # 32768

    print(f"COMPUTING LEVEL-4 FOURIER COEFFICIENTS AT n={n}")
    print(f"m = {m}, N = {N}")

    arcs = [(i,j) for i in range(n) for j in range(i+1,n)]

    # Compute H for all tournaments
    h_vals = np.zeros(N)
    for bits in range(N):
        adj = [[0]*n for _ in range(n)]
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if bits & (1 << idx):
                    adj[i][j] = 1
                else:
                    adj[j][i] = 1
                idx += 1
        h_vals[bits] = count_ham_paths(adj, n)

    mu = np.mean(h_vals)
    print(f"Mean(H) = {mu}")

    # Compute level-4 coefficients
    print(f"\nComputing level-4 Fourier coefficients (C({m},4) = {math.comb(m,4)} subsets)...")

    nonzero_coeffs = []
    e4_total = 0
    count = 0

    for s_idx in combinations(range(m), 4):
        # Compute chi_S for each tournament
        chi_vals = np.ones(N)
        for e in s_idx:
            signs = np.array([(1 if (bits & (1 << e)) else -1) for bits in range(N)])
            chi_vals *= signs
        coeff = np.dot(h_vals, chi_vals) / N
        e4_total += coeff**2
        if abs(coeff) > 1e-10:
            arc_set = [arcs[e] for e in s_idx]
            vertices = set()
            for arc in arc_set:
                vertices.update(arc)
            nonzero_coeffs.append((s_idx, coeff, len(vertices)))
        count += 1
        if count % 500 == 0:
            print(f"  {count}/{math.comb(m,4)} done, {len(nonzero_coeffs)} nonzero so far")

    print(f"\nTotal level-4 subsets: {count}")
    print(f"Nonzero: {len(nonzero_coeffs)}")
    print(f"E_4 = {e4_total}")
    print(f"E_4/E_0 = {e4_total/mu**2}")
    print(f"Expected E_4 = {45/4}")

    # Analyze nonzero coefficients
    from collections import Counter
    mags = Counter()
    vertex_coverage = Counter()
    for (s_idx, coeff, nv) in nonzero_coeffs:
        mags[round(abs(coeff), 8)] += 1
        vertex_coverage[nv] += 1

    print(f"\nUnique magnitudes: {sorted(mags.items())}")
    print(f"Vertex coverage: {sorted(vertex_coverage.items())}")

    # Print some examples of each type
    for mag in sorted(mags.keys()):
        examples = [(s, c, nv) for (s, c, nv) in nonzero_coeffs if abs(abs(c) - mag) < 1e-6]
        print(f"\n  |H_hat| = {mag}: {len(examples)} coefficients")
        for (s_idx, coeff, nv) in examples[:5]:
            arc_set = [arcs[e] for e in s_idx]
            degs = [0]*n
            for arc in arc_set:
                degs[arc[0]] += 1
                degs[arc[1]] += 1
            print(f"    {arc_set}, vertices={nv}, deg={sorted(degs)}, sign={'+' if coeff > 0 else '-'}")

    print(f"\n{'='*70}")
    print("DONE")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
