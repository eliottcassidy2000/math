"""
e4_formula_search.py -- kind-pasteur-2026-03-14-S105g
Search for a closed-form formula for E_4/E_0.

Known:
  n=5: E_4/E_0 = 1/60
  n=6: E_4/E_0 = 1/45

Key observations:
  60 = 5*4*3 = P(5,3) (falling factorial)
  45 = 9*5 = 3^2 * 5

  60 = C(5,2)*C(4,2) = 10*6? No, 10*6=60. YES!
  45 = C(6,2)*C(5,2)/(something)? C(6,2)=15, 15*3=45.

Let me think about what LEVEL-4 Fourier means.

Level-4 coefficients H_hat(S) where |S|=4 (four arcs).
From the exact formula at level 2: |H_hat({e1,e2})| = (n-2)!/2^(n-2)
when arcs e1,e2 share a vertex.

What is |H_hat({e1,e2,e3,e4})| for a 4-arc subset?
This depends on the intersection pattern of the 4 arcs.

The level-4 energy E_4 = sum_{|S|=4} H_hat(S)^2.

Possible 4-arc patterns (by vertex coverage):
  - All 4 arcs through a single vertex: "star pattern"
  - Two pairs sharing different vertices
  - Chain pattern
  - etc.

Let me compute E_4 directly from the data.
"""

import sys, math
import numpy as np
from fractions import Fraction
from collections import Counter
from itertools import combinations

sys.stdout.reconfigure(encoding='utf-8')

def count_ham_paths(adj, n):
    """Count Hamiltonian paths using Held-Karp DP."""
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
    n = 5
    m = n*(n-1)//2
    N = 1 << m

    print(f"COMPUTING LEVEL-4 FOURIER COEFFICIENTS AT n={n}")
    print(f"m = {m}, N = {N}")

    # Enumerate arcs
    arcs = []
    for i in range(n):
        for j in range(i+1, n):
            arcs.append((i, j))
    print(f"Arcs: {arcs}")

    # Compute H for all tournaments (in {-1,+1} basis)
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

    # Compute ALL Fourier coefficients
    # H_hat(S) = (1/N) * sum_T H(T) * chi_S(T)
    # where chi_S(T) = prod_{e in S} y_e(T), y_e = 2*x_e - 1

    # Convert to {-1,+1} basis
    # For each tournament T (bits), y_e = +1 if bit e is set, -1 if not
    # Actually: we need to be consistent. Let's say:
    # bit e = 1 means arc e oriented "forward" (i->j), y_e = +1
    # bit e = 0 means arc e oriented "backward" (j->i), y_e = -1

    # Compute Fourier coefficients for ALL subsets
    level_energy = {}  # level -> total energy

    # Level 0
    h_hat_0 = mu
    level_energy[0] = h_hat_0**2
    print(f"\nLevel 0: H_hat(empty) = {h_hat_0}, energy = {h_hat_0**2}")

    # Level 2
    e2_total = 0
    level2_coeffs = []
    for s_idx in combinations(range(m), 2):
        # chi_S(T) = y_{s1} * y_{s2}
        coeff = 0
        for bits in range(N):
            chi = 1
            for e in s_idx:
                if bits & (1 << e):
                    chi *= 1
                else:
                    chi *= -1
            coeff += h_vals[bits] * chi
        coeff /= N
        e2_total += coeff**2
        if abs(coeff) > 1e-10:
            level2_coeffs.append((s_idx, coeff))

    level_energy[2] = e2_total
    print(f"\nLevel 2: {len(level2_coeffs)} nonzero coefficients")
    print(f"  Energy = {e2_total}")
    print(f"  E_2/E_0 = {e2_total/h_hat_0**2}")

    # Print a few level-2 coefficients
    for (s_idx, coeff) in level2_coeffs[:10]:
        arc_pair = tuple(arcs[e] for e in s_idx)
        # Check if they share a vertex
        v1 = set(arc_pair[0])
        v2 = set(arc_pair[1])
        shared = v1 & v2
        print(f"    S={s_idx} = {arc_pair}, shared={shared}, H_hat = {coeff:.6f}")

    # Level 4
    e4_total = 0
    level4_coeffs = []
    for s_idx in combinations(range(m), 4):
        chi_vals = np.ones(N)
        for e in s_idx:
            signs = np.array([(1 if (bits & (1 << e)) else -1) for bits in range(N)])
            chi_vals *= signs
        coeff = np.dot(h_vals, chi_vals) / N
        e4_total += coeff**2
        if abs(coeff) > 1e-10:
            level4_coeffs.append((s_idx, coeff))

    level_energy[4] = e4_total
    print(f"\nLevel 4: {len(level4_coeffs)} nonzero coefficients")
    print(f"  Energy = {e4_total}")
    print(f"  E_4/E_0 = {e4_total/h_hat_0**2}")

    # Print ALL level-4 coefficients
    print(f"\n  ALL level-4 nonzero coefficients:")
    unique_mags = set()
    for (s_idx, coeff) in level4_coeffs:
        arc_set = tuple(arcs[e] for e in s_idx)
        # Find vertices covered
        vertices = set()
        for arc in arc_set:
            vertices.update(arc)
        print(f"    S={s_idx} arcs={arc_set} vertices={vertices} |V|={len(vertices)} H_hat={coeff:.6f}")
        unique_mags.add(round(abs(coeff), 8))

    print(f"\n  Unique magnitudes: {sorted(unique_mags)}")
    print(f"  Number of nonzero level-4 coefficients: {len(level4_coeffs)}")
    print(f"  E_4 = {e4_total}")
    print(f"  E_4/E_0 = {Fraction(int(round(e4_total * 900)), 900) if e4_total > 0 else 0}")

    # Expected: E_4 = 15/16
    print(f"\n  Expected E_4 = 15/16 = {15/16}")
    print(f"  Computed E_4 = {e4_total}")
    print(f"  Match: {abs(e4_total - 15/16) < 1e-10}")

    if len(level4_coeffs) > 0:
        mag = abs(level4_coeffs[0][1])
        print(f"\n  Each level-4 coefficient magnitude: {mag}")
        print(f"  Number of nonzero: {len(level4_coeffs)}")
        print(f"  Check: {len(level4_coeffs)} * {mag}^2 = {len(level4_coeffs) * mag**2}")
        print(f"  = {e4_total}")

    # Check what structure the nonzero 4-arc subsets have
    print(f"\n  STRUCTURE ANALYSIS of nonzero level-4 subsets:")
    vertex_counts = Counter()
    for (s_idx, coeff) in level4_coeffs:
        arc_set = [arcs[e] for e in s_idx]
        vertices = set()
        for arc in arc_set:
            vertices.update(arc)
        vertex_counts[len(vertices)] += 1

    for v_count, freq in sorted(vertex_counts.items()):
        print(f"    {v_count} vertices covered: {freq} subsets")

    # Level 6 (should not exist at n=5 since deg=4)
    print(f"\n  Checking level 6 (should be zero at n=5 by Degree Drop):")
    e6_sample = 0
    count = 0
    for s_idx in combinations(range(m), 6):
        chi_vals = np.ones(N)
        for e in s_idx:
            signs = np.array([(1 if (bits & (1 << e)) else -1) for bits in range(N)])
            chi_vals *= signs
        coeff = np.dot(h_vals, chi_vals) / N
        e6_sample += coeff**2
        count += 1
        if count >= 100:  # Just check first 100
            break
    print(f"    First {count} level-6 coefficients: total energy = {e6_sample:.10f}")
    print(f"    (Should be 0 by Degree Drop)")

    # Summary
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    total = sum(level_energy.values())
    for level, energy in sorted(level_energy.items()):
        pct = 100*energy/total
        print(f"  Level {level}: E_{level} = {energy:.6f} ({pct:.2f}%)")
    print(f"  Total: {total:.6f} (= Mean(H^2))")

if __name__ == '__main__':
    main()
