#!/usr/bin/env python3
"""
H_fourier_formula.py -- kind-pasteur-2026-03-13-S61

THEOREM CANDIDATE (THM-163): H as a Walsh-Hadamard expansion

Discovered structure:
  H(sigma) = E[H] + c_2 * sum_{edges sharing vertex} sign(e1,e2) * sigma_{e1} * sigma_{e2}
           + c_4 * (degree-4 terms) + ...

where:
  E[H] = n!/2^{n-1}
  c_2 = (n-2)!/2^{n-2}
  sign(e1,e2) depends on edge-vertex incidence

Degree-2 coefficient formula:
  For edges (a,b) and (c,d) sharing exactly one vertex v:
  hat(H)({e1,e2}) = (-1)^{pos_1 + pos_2 + 1} * (n-2)!/2^{n-2}
  where pos_i = 0 if v is the smaller vertex of edge i, 1 if larger.

For edges sharing NO vertex: hat(H)({e1,e2}) = 0.

This means the degree-2 part of H is:
  H_2(sigma) = (n-2)!/2^{n-2} * sum_{v in [n]} sum_{a<v<b} sigma_{av}*sigma_{vb}
             - (n-2)!/2^{n-2} * sum_{v in [n]} [sum_{a<b, both adj to v, v smaller in both} sigma_{va}*sigma_{vb}
                                                 + same for v larger in both]

This simplifies to a GRAPH LAPLACIAN structure!

Author: kind-pasteur-2026-03-13-S61
"""

import math
from itertools import combinations
from collections import defaultdict


def binary_to_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << pos):
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A


def tournament_to_sigma(A, n):
    vec = []
    for i in range(n):
        for j in range(i+1, n):
            vec.append(1 if A[i][j] else -1)
    return tuple(vec)


def count_ham_paths(A, n):
    if n <= 1:
        return 1
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + dp[key]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


print("=" * 70)
print("THM-163 CANDIDATE: H WALSH-HADAMARD COEFFICIENT FORMULA")
print("=" * 70)

# Verify degree-2 coefficient formula: |hat(H)| = (n-2)!/2^{n-2}
print("\n--- Degree-2 coefficient magnitude formula ---")
print(f"{'n':>3} {'(n-2)!/2^(n-2)':>20} {'observed |c_2|':>20} {'match':>8}")

for n in range(3, 8):
    m = n * (n - 1) // 2
    total = 1 << m

    predicted = math.factorial(n - 2) / (2 ** (n - 2))

    if n <= 6:
        # Compute exact Fourier coefficient for one edge pair
        H_vals = {}
        for bits in range(total):
            A = binary_to_tournament(bits, n)
            sigma = tournament_to_sigma(A, n)
            H = count_ham_paths(A, n)
            H_vals[sigma] = H

        # Pick edges (0,1) and (0,2) — share vertex 0
        edge_pairs = []
        for i in range(n):
            for j in range(i+1, n):
                edge_pairs.append((i, j))

        # Find indices of (0,1) and (0,2)
        idx_01 = edge_pairs.index((0, 1))
        idx_02 = edge_pairs.index((0, 2))

        # Compute hat(H)({01, 02})
        coeff = 0
        for sigma, H in H_vals.items():
            coeff += H * sigma[idx_01] * sigma[idx_02]
        coeff /= total

        observed = abs(coeff)
        match = abs(observed - predicted) < 1e-10
        print(f"{n:>3} {predicted:>20.6f} {observed:>20.6f} {str(match):>8}")

        # Also verify sign rule
        # (0,1) and (0,2): shared vertex 0, position 0 in both => sign = (-1)^{0+0+1} = -1
        expected_sign = -1
        actual_sign = 1 if coeff > 0 else -1
        sign_match = expected_sign == actual_sign
        print(f"    Sign: predicted {expected_sign}, actual {actual_sign}, match={sign_match}")
    else:
        print(f"{n:>3} {predicted:>20.6f} {'(too large)':>20}")


# ========================================================================
# VERIFY SIGN RULE FOR ALL EDGE PAIRS AT n=5
# ========================================================================
print(f"\n{'='*70}")
print("SIGN RULE VERIFICATION AT n=5")
print("=" * 70)

n = 5
m = n * (n - 1) // 2
total = 1 << m
predicted_mag = math.factorial(n - 2) / (2 ** (n - 2))

H_vals = {}
for bits in range(total):
    A = binary_to_tournament(bits, n)
    sigma = tournament_to_sigma(A, n)
    H = count_ham_paths(A, n)
    H_vals[sigma] = H

edge_pairs = []
for i in range(n):
    for j in range(i+1, n):
        edge_pairs.append((i, j))

print(f"\nEdge pair      | Shared v | pos1 | pos2 | predicted sign | actual coeff | match")
print("-" * 85)

all_match = True
for i in range(m):
    for j in range(i+1, m):
        e1 = edge_pairs[i]
        e2 = edge_pairs[j]
        shared = set(e1) & set(e2)

        if not shared:
            # No shared vertex => coefficient should be 0
            coeff = 0
            for sigma, H in H_vals.items():
                coeff += H * sigma[i] * sigma[j]
            coeff /= total
            if abs(coeff) > 1e-10:
                print(f"({e1[0]},{e1[1]})-({e2[0]},{e2[1]}) | none     | -    | -    | 0              | {coeff:.6f} | FAIL!")
                all_match = False
            continue

        v = shared.pop()
        pos1 = 0 if v == e1[0] else 1
        pos2 = 0 if v == e2[0] else 1
        predicted_sign = (-1) ** (pos1 + pos2 + 1)

        coeff = 0
        for sigma, H in H_vals.items():
            coeff += H * sigma[i] * sigma[j]
        coeff /= total

        actual_sign = 1 if coeff > 0 else -1
        match = (actual_sign == predicted_sign and abs(abs(coeff) - predicted_mag) < 1e-10)
        if not match:
            all_match = False

        print(f"({e1[0]},{e1[1]})-({e2[0]},{e2[1]}) | v={v}      | {pos1}    | {pos2}    | {'+' if predicted_sign > 0 else '-'}              | {coeff:>12.6f} | {'OK' if match else 'FAIL!'}")

print(f"\nAll degree-2 coefficients match formula: {all_match}")


# ========================================================================
# ALTERNATIVE FORMULA: H_2 AS SCORE REGULARITY MEASURE
# ========================================================================
print(f"\n{'='*70}")
print("H_2 AS SCORE REGULARITY MEASURE")
print("=" * 70)

print("""
The degree-2 part of H can be rewritten as:

  H_2(sigma) = c * sum_{v} [sum_{a<v} sigma_{av}] * [sum_{b>v} sigma_{vb}]
             - c * sum_{v} [C(out_v, 2) + C(in_v, 2) - ...terms...]

Actually, let's express it in terms of SCORES.
  sigma_{ij} = +1 if i->j, -1 if j->i (for i<j)
  out-degree of vertex v: d_v^+ = sum_{w != v} A[v][w]
  in-degree: d_v^- = n-1 - d_v^+
  Score: s_v = d_v^+

  sum_{w: w>v} sigma_{vw} = d_v^+(among w>v) - d_v^-(among w>v)
  sum_{w: w<v} sigma_{wv} = d_w^+(w->v among w<v) - d_w^-(v->w among w<v)
                           = (in from below) - (out to below)

Let me just compute H - E[H] and see if it correlates with score variance.
""")

n = 5
m = n * (n - 1) // 2
total = 1 << m
E_H = math.factorial(n) / (2 ** (n - 1))

score_var_vs_H = []
for bits in range(total):
    A = binary_to_tournament(bits, n)
    H = count_ham_paths(A, n)
    scores = [sum(A[v]) for v in range(n)]
    mean_s = sum(scores) / n
    var_s = sum((s - mean_s)**2 for s in scores) / n
    score_var_vs_H.append((var_s, H))

# Group by score variance
by_var = defaultdict(list)
for var_s, H in score_var_vs_H:
    by_var[round(var_s, 4)].append(H)

print(f"n={n}: E[H] = {E_H}")
print(f"\nScore variance | avg H   | min H | max H | count")
print("-" * 55)
for var_s in sorted(by_var.keys()):
    H_list = by_var[var_s]
    print(f"{var_s:>14.4f} | {sum(H_list)/len(H_list):>7.2f} | {min(H_list):>5d} | {max(H_list):>5d} | {len(H_list):>5d}")

# Compute correlation
import statistics
vars_list = [v for v, _ in score_var_vs_H]
H_list_all = [h for _, h in score_var_vs_H]
if len(set(vars_list)) > 1:
    # Pearson correlation
    mean_v = sum(vars_list) / len(vars_list)
    mean_h = sum(H_list_all) / len(H_list_all)
    cov = sum((v - mean_v) * (h - mean_h) for v, h in zip(vars_list, H_list_all)) / len(vars_list)
    std_v = (sum((v - mean_v)**2 for v in vars_list) / len(vars_list)) ** 0.5
    std_h = (sum((h - mean_h)**2 for h in H_list_all) / len(H_list_all)) ** 0.5
    corr = cov / (std_v * std_h) if std_v > 0 and std_h > 0 else 0
    print(f"\ncorr(score_variance, H) = {corr:.6f}")


# ========================================================================
# DEGREE-2 FORMULA IN TERMS OF VERTEX BALANCE
# ========================================================================
print(f"\n{'='*70}")
print("DEGREE-2 H AS VERTEX BALANCE MEASURE")
print("=" * 70)

# Define balance B_v = (out_v - in_v) = 2*s_v - (n-1)
# Then sigma_{vw} depends on whether v beats w
# H_2 involves products sigma_{va} * sigma_{vb} for edges sharing v

# The degree-2 contribution from vertex v:
# sum_{a<v, b>v} sigma_{av}*sigma_{vb} * sign = ?
# This equals a weighted cross product of "above" and "below" performance

print(f"\nDegree-2 part of H at n={n}:")
print(f"  c_2 = (n-2)!/2^(n-2) = {predicted_mag}")

# For each tournament, compute H_2 from formula and verify
for bits in [0, 1, 76, 341, 682, 1023]:
    if bits >= total:
        continue
    A = binary_to_tournament(bits, n)
    sigma = tournament_to_sigma(A, n)
    H = count_ham_paths(A, n)
    scores = [sum(A[v]) for v in range(n)]

    # Compute H_2 from Fourier formula
    H_2_fourier = 0
    for idx_i in range(m):
        for idx_j in range(idx_i + 1, m):
            e1 = edge_pairs[idx_i]
            e2 = edge_pairs[idx_j]
            shared = set(e1) & set(e2)
            if not shared:
                continue
            v = shared.pop()
            pos1 = 0 if v == e1[0] else 1
            pos2 = 0 if v == e2[0] else 1
            sign = (-1) ** (pos1 + pos2 + 1)
            H_2_fourier += predicted_mag * sign * sigma[idx_i] * sigma[idx_j]

    # Compute H_0 + H_2 and compare to H
    H_approx = E_H + H_2_fourier
    residual = H - H_approx

    print(f"\n  bits={bits}: H={H}, scores={scores}")
    print(f"    H_0 = {E_H:.1f}")
    print(f"    H_2 = {H_2_fourier:.4f}")
    print(f"    H_0 + H_2 = {H_approx:.4f}")
    print(f"    Residual (degree >= 4) = {residual:.4f}")
    print(f"    H_2 / (H - E_H) = {H_2_fourier / (H - E_H):.4f}" if H != E_H else "")


# ========================================================================
# THE DEEP FORMULA: H_2 = (n-2)!/2^{n-2} * sum_v Z_v
# ========================================================================
print(f"\n{'='*70}")
print("H_2 IN TERMS OF VERTEX CROSS-PRODUCTS Z_v")
print("=" * 70)

print("""
Define for each vertex v:
  Z_v = sum_{a<v<b} sigma_{av} * sigma_{vb}
      - sum_{a<b, both<v} sigma_{av} * sigma_{bv}
      - sum_{a<b, both>v} sigma_{va} * sigma_{vb}

This measures how much v "channels" Hamiltonian paths through it.
  Z_v > 0: v is a "bridge" (different above/below orientation)
  Z_v < 0: v is a "bottleneck" (consistent orientation)

Then: H_2 = (n-2)!/2^{n-2} * sum_v Z_v
""")

n = 5
m = n * (n - 1) // 2
total = 1 << m
c2 = math.factorial(n - 2) / (2 ** (n - 2))

# Verify the Z_v formula
for bits in [0, 76, 341, 1023]:
    if bits >= total:
        continue
    A = binary_to_tournament(bits, n)
    sigma = tournament_to_sigma(A, n)
    H = count_ham_paths(A, n)

    # Compute Z_v for each vertex
    # Need to convert sigma to a dict for easy access
    sigma_dict = {}
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            sigma_dict[(i, j)] = sigma[pos]
            sigma_dict[(j, i)] = -sigma[pos]
            pos += 1

    Z = [0] * n
    for v in range(n):
        # Cross terms: a < v < b
        for a in range(v):
            for b in range(v + 1, n):
                # sigma_{av} * sigma_{vb}
                # sigma_{av}: since a < v, this is sigma_dict[(a,v)]
                # sigma_{vb}: since v < b, this is sigma_dict[(v,b)]
                Z[v] += sigma_dict[(a, v)] * sigma_dict[(v, b)]

        # Same-side-below: a < b < v
        for a in range(v):
            for b in range(a + 1, v):
                # Both edges have v as the larger vertex
                # sigma_{av} * sigma_{bv}
                Z[v] -= sigma_dict[(a, v)] * sigma_dict[(b, v)]

        # Same-side-above: v < a < b
        for a in range(v + 1, n):
            for b in range(a + 1, n):
                # Both edges have v as the smaller vertex
                # sigma_{va} * sigma_{vb}
                Z[v] -= sigma_dict[(v, a)] * sigma_dict[(v, b)]

    total_Z = sum(Z)
    H_2_pred = c2 * total_Z
    H_2_actual = H - math.factorial(n) / (2 ** (n - 1))

    # Also need to subtract higher-degree terms for exact comparison
    # But if the degree-4 terms are small, H_2 should be close

    print(f"\n  bits={bits}: H={H}")
    print(f"    Z = {Z}")
    print(f"    sum(Z) = {total_Z}")
    print(f"    c2 * sum(Z) = {H_2_pred:.4f}")
    print(f"    H - E[H] = {H_2_actual:.4f}")
    print(f"    Match (deg-2 only): {'YES' if abs(H_2_pred - H_2_actual) < 0.01 else 'NO (higher terms)'}")

    # Verify: Z_v relates to score?
    scores = [sum(A[v]) for v in range(n)]
    print(f"    Scores: {scores}")


# ========================================================================
# BEAUTIFUL FORMULA: Z_v = (s_v - (n-1)/2) * (... something?)
# ========================================================================
print(f"\n{'='*70}")
print("Z_v vs SCORE ANALYSIS")
print("=" * 70)

n = 5
m = n * (n - 1) // 2
total = 1 << m
c2 = math.factorial(n - 2) / (2 ** (n - 2))

# Collect Z_v and score data across all tournaments
z_vs_score = defaultdict(list)

for bits in range(total):
    A = binary_to_tournament(bits, n)
    sigma = tournament_to_sigma(A, n)

    sigma_dict = {}
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            sigma_dict[(i, j)] = sigma[pos]
            sigma_dict[(j, i)] = -sigma[pos]
            pos += 1

    scores = [sum(A[v]) for v in range(n)]

    for v in range(n):
        Z = 0
        for a in range(v):
            for b in range(v + 1, n):
                Z += sigma_dict[(a, v)] * sigma_dict[(v, b)]
        for a in range(v):
            for b in range(a + 1, v):
                Z -= sigma_dict[(a, v)] * sigma_dict[(b, v)]
        for a in range(v + 1, n):
            for b in range(a + 1, n):
                Z -= sigma_dict[(v, a)] * sigma_dict[(v, b)]

        z_vs_score[scores[v]].append(Z)

print(f"  Score | avg Z_v  | min Z_v | max Z_v | count")
print(f"  ------+----------+---------+---------+------")
for s in sorted(z_vs_score):
    zs = z_vs_score[s]
    print(f"  {s:>5d} | {sum(zs)/len(zs):>8.4f} | {min(zs):>7d} | {max(zs):>7d} | {len(zs):>5d}")


# ========================================================================
# ENERGY FRACTION BY DEGREE
# ========================================================================
print(f"\n{'='*70}")
print("FOURIER ENERGY FRACTIONS BY n")
print("=" * 70)

for n in range(3, 7):
    m = n * (n - 1) // 2
    total = 1 << m

    H_vals = {}
    for bits in range(total):
        A = binary_to_tournament(bits, n)
        sigma = tournament_to_sigma(A, n)
        H = count_ham_paths(A, n)
        H_vals[sigma] = H

    # Compute energy by degree
    energy_by_deg = defaultdict(float)
    for k in range(0, m + 1, 2):  # only even degrees
        for S in combinations(range(m), k):
            coeff = 0
            for sigma, H in H_vals.items():
                chi = 1
                for idx in S:
                    chi *= sigma[idx]
                coeff += H * chi
            coeff /= total
            energy_by_deg[k] += coeff ** 2

    total_energy = sum(energy_by_deg.values())

    print(f"\n  n={n} (m={m}):")
    print(f"    E[H] = {math.factorial(n)/2**(n-1):.2f}")
    for deg in sorted(energy_by_deg):
        e = energy_by_deg[deg]
        pct = 100 * e / total_energy if total_energy > 0 else 0
        print(f"    Degree {deg}: energy = {e:.4f} ({pct:.2f}%)")

    # Verify: does the ratio degree-0/total stay at 75%?
    if 0 in energy_by_deg:
        ratio_0 = energy_by_deg[0] / total_energy
        print(f"    Degree-0 fraction = {ratio_0:.6f}")
        # E[H]^2 / E[H^2] = ?
        EH2 = sum(H**2 for H in H_vals.values()) / total
        EH = sum(H_vals.values()) / total
        print(f"    E[H]^2 / E[H^2] = {EH**2 / EH2:.6f}")

print("\n\nDONE.")
