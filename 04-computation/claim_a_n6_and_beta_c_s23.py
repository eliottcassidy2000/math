#!/usr/bin/env python3
"""
claim_a_n6_and_beta_c_s23.py — opus-2026-04-04-S23

1. Claim A decomposition at n=5→6 (32768 extensions)
2. β_c at n=7 (sampled)
3. The a_inter mechanism: μ correlation with H_sub
"""

import numpy as np
from itertools import combinations, permutations
from math import comb, factorial
from collections import defaultdict
import random
import time
import sys

random.seed(42)

def compute_H(A, n):
    if n <= 1: return 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if not (mask & (1 << u)) and A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[(1 << n) - 1])

def subtournament(A, n, vertices):
    verts = sorted(vertices)
    k = len(verts)
    B = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            B[i][j] = A[verts[i]][verts[j]]
    return B, k

def all_tournaments(n):
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i, j) in enumerate(edges):
            if (mask >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A

def extend(A_sub, n_sub, sigma):
    n = n_sub + 1
    A = [[0]*n for _ in range(n)]
    for i in range(n_sub):
        for j in range(n_sub):
            A[i][j] = A_sub[i][j]
    for i in range(n_sub):
        if sigma[i]:
            A[n_sub][i] = 1
        else:
            A[i][n_sub] = 1
    return A

def delta_c3(A_sub, n_sub, sigma):
    return sum(A_sub[i][j] for i in range(n_sub) if sigma[i] == 1
                           for j in range(n_sub) if sigma[j] == 0)


def claim_a_decomp_n6():
    """Decompose ΔH at n=5→6 by cycle length."""
    n_target = 6
    n_sub = 5
    new_v = 5

    print(f"\n{'='*70}")
    print(f"  CLAIM A DECOMPOSITION: n = {n_sub} → {n_target}")
    print(f"{'='*70}")

    total_data = []
    count = 0

    for A_sub in all_tournaments(n_sub):
        H_sub = compute_H(A_sub, n_sub)

        for sigma_int in range(1 << n_sub):
            sigma = tuple((sigma_int >> i) & 1 for i in range(n_sub))
            A_n = extend(A_sub, n_sub, sigma)
            H_n = compute_H(A_n, n_target)
            dH = H_n - H_sub
            dc3 = delta_c3(A_sub, n_sub, sigma)

            # Find 3-cycles through new_v and compute μ
            c3_contribution = 0
            for i in range(n_sub):
                for j in range(n_sub):
                    if i == j: continue
                    # Check cycle new_v → i → j → new_v
                    if A_n[new_v][i] and A_n[i][j] and A_n[j][new_v]:
                        # This is a directed 3-cycle (new_v, i, j)
                        # μ = H(T[V \ {new_v, i, j}])
                        complement = [u for u in range(n_target) if u not in (new_v, i, j)]
                        B, m = subtournament(A_n, n_target, complement)
                        mu = compute_H(B, m)
                        c3_contribution += mu
            # Each UNDIRECTED 3-cycle counted ONCE (i→j direction is the valid one)
            # Actually we're iterating over ordered (i,j) pairs so each directed cycle
            # is counted once (new_v→i→j→new_v). Each undirected triple gives exactly
            # one directed 3-cycle. So c3_contribution = Σ_C μ(C) for 3-cycles.

            # Find 5-cycles through new_v and compute μ
            c5_contribution = 0
            other_verts = list(range(n_sub))
            for combo in combinations(other_verts, 4):
                # 5-cycle on {new_v} ∪ combo
                for perm in permutations(combo):
                    cycle = [new_v] + list(perm)
                    if all(A_n[cycle[i]][cycle[(i+1) % 5]] for i in range(5)):
                        # Found a 5-cycle. V(C) = {new_v} ∪ combo, complement has 1 vertex.
                        complement = [u for u in range(n_target) if u not in set(cycle)]
                        if len(complement) == 0:
                            mu = 1
                        else:
                            B, m = subtournament(A_n, n_target, complement)
                            mu = compute_H(B, m)
                        c5_contribution += mu
                        break  # One directed 5-cycle per vertex set

            claim_a_check = 2 * (c3_contribution + c5_contribution)

            total_data.append({
                'dH': dH, 'dc3': dc3, 'H_sub': H_sub,
                'c3_contrib': 2 * c3_contribution,
                'c5_contrib': 2 * c5_contribution,
                'claim_a': claim_a_check,
            })

        count += 1
        if count % 100 == 0:
            print(f"  ... {count}/1024 parents", file=sys.stderr)

    N = len(total_data)

    # Verify Claim A
    mismatches = sum(1 for d in total_data if d['claim_a'] != d['dH'])
    print(f"\n  Claim A mismatches: {mismatches}/{N}")

    all_dc3 = np.array([d['dc3'] for d in total_data], dtype=float)
    all_dH = np.array([d['dH'] for d in total_data], dtype=float)
    C3 = np.array([d['c3_contrib'] for d in total_data], dtype=float)
    C5 = np.array([d['c5_contrib'] for d in total_data], dtype=float)

    # Regression of each contribution vs Δc₃
    print(f"\n  Contribution by cycle length:")
    print(f"  {'k':>3}  {'mean':>12}  {'fraction':>12}  {'coeff vs Δc₃':>15}")

    total_coeff = 0
    for k, Ck in [(3, C3), (5, C5)]:
        mean_Ck = Ck.mean()
        frac = mean_Ck / all_dH.mean() if all_dH.mean() > 0 else 0
        X = np.column_stack([all_dc3, np.ones(N)])
        c, _, _, _ = np.linalg.lstsq(X, Ck, rcond=None)
        total_coeff += c[0]
        print(f"  {k:3d}  {mean_Ck:12.4f}  {frac:12.6f}  {c[0]:15.6f}")

    print(f"  {'SUM':>3}  {all_dH.mean():12.4f}  {'1.000000':>12}  {total_coeff:15.6f}")

    predicted = factorial(n_target - 2) / 2**(n_target - 4)
    print(f"\n  Total coefficient: {total_coeff:.6f}")
    print(f"  Predicted a({n_target}) = {predicted:.6f}")
    print(f"  Ratio: {total_coeff/predicted:.6f}")

    # Now the KEY question: what is the average μ for 3-cycles?
    # μ = H(subtournament on 3 verts) for 3-cycles at n=6
    mu3_vals = []
    for d in total_data:
        if d['dc3'] > 0:
            avg_mu = d['c3_contrib'] / (2 * d['dc3'])  # c3_contrib = 2*Σμ, so Σμ/dc3
            mu3_vals.append(avg_mu)
    if mu3_vals:
        print(f"\n  Average μ for 3-cycles (= avg H(T_3 subtournament)):")
        print(f"    Mean: {np.mean(mu3_vals):.6f}")
        print(f"    This should be related to mean_H(3) = 7.5/4 = 1.875... no")
        print(f"    Mean H of random T_3: (1*6+3*2)/8 = 12/8 = 1.5")

    # Decompose μ further: how does μ depend on H_sub?
    mu3_by_Hsub = defaultdict(list)
    for d in total_data:
        if d['dc3'] > 0:
            mu3_by_Hsub[d['H_sub']].append(d['c3_contrib'] / (2 * d['dc3']))

    print(f"\n  Average μ(3-cycle) by parent H:")
    for H_sub in sorted(mu3_by_Hsub.keys()):
        vals = mu3_by_Hsub[H_sub]
        print(f"    H_sub={H_sub:3d}: mean μ = {np.mean(vals):.6f} (n={len(vals)})")

    return total_data


def beta_c_n7():
    """Compute phase transition at n=7 via sampling."""
    n = 7
    n_samples = 5000
    print(f"\n{'='*70}")
    print(f"  PHASE TRANSITION at n = {n} ({n_samples} samples)")
    print(f"{'='*70}")

    H_list = []
    for _ in range(n_samples):
        A = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
        H_list.append(compute_H(A, n))

    H_arr = np.array(H_list, dtype=float)
    print(f"  H range: [{H_arr.min():.0f}, {H_arr.max():.0f}], mean={H_arr.mean():.2f}")

    beta_range = np.linspace(0, 0.5, 501)
    max_C = 0
    beta_max = 0

    for beta in beta_range:
        log_w = beta * H_arr
        log_w -= log_w.max()
        w = np.exp(log_w)
        w /= w.sum()
        mean_H = np.sum(w * H_arr)
        var_H = np.sum(w * (H_arr - mean_H)**2)
        C = beta**2 * var_H
        if C > max_C:
            max_C = C
            beta_max = beta

    print(f"  Specific heat peak: C_max = {max_C:.4f} at β_c = {beta_max:.4f}")
    return beta_max


def a_inter_mechanism():
    """Investigate WHY a_inter ≈ 0.27.

    In Claim A: ΔH = 2Σ μ(C). For 3-cycles, μ = H(complement subtournament).
    The interaction term a_inter · Δc₃ · H_sub says μ ∝ H_sub.

    Let's test: is μ(3-cycle C) ≈ α + β·H_sub for some α, β?
    """
    n_sub = 5
    n_target = 6
    new_v = 5
    print(f"\n{'='*70}")
    print(f"  a_inter MECHANISM: μ vs H_sub at n = {n_sub} → {n_target}")
    print(f"{'='*70}")

    data = []
    count = 0
    for A_sub in all_tournaments(n_sub):
        H_sub = compute_H(A_sub, n_sub)

        for sigma_int in range(1 << n_sub):
            sigma = tuple((sigma_int >> i) & 1 for i in range(n_sub))
            A_n = extend(A_sub, n_sub, sigma)

            # For each 3-cycle through new_v, compute μ
            for i in range(n_sub):
                for j in range(n_sub):
                    if i == j: continue
                    if A_n[new_v][i] and A_n[i][j] and A_n[j][new_v]:
                        complement = [u for u in range(n_target) if u not in (new_v, i, j)]
                        B, m = subtournament(A_n, n_target, complement)
                        mu = compute_H(B, m)
                        data.append({
                            'mu': mu, 'H_sub': H_sub,
                            'complement_verts': tuple(complement),
                        })
        count += 1
        if count % 200 == 0:
            print(f"  ... {count}/1024 parents", file=sys.stderr)

    N = len(data)
    mu_arr = np.array([d['mu'] for d in data], dtype=float)
    H_sub_arr = np.array([d['H_sub'] for d in data], dtype=float)

    print(f"\n  {N} (3-cycle, parent) pairs")
    print(f"  μ range: [{mu_arr.min():.0f}, {mu_arr.max():.0f}], mean={mu_arr.mean():.4f}")
    print(f"  corr(μ, H_sub) = {np.corrcoef(mu_arr, H_sub_arr)[0,1]:.6f}")

    # Regression: μ ≈ a + b·H_sub
    X = np.column_stack([H_sub_arr, np.ones(N)])
    c, _, _, _ = np.linalg.lstsq(X, mu_arr, rcond=None)
    R2 = 1 - np.sum((mu_arr - X@c)**2) / np.sum((mu_arr - mu_arr.mean())**2)
    print(f"  μ ≈ {c[0]:.6f}·H_sub + {c[1]:.6f}, R²={R2:.6f}")
    print(f"  → Each unit of H_sub increases μ by {c[0]:.6f}")
    print(f"  → This maps to a_inter because:")
    print(f"     ΔH from 3-cycles = 2·Δc₃·mean_μ ≈ 2·Δc₃·({c[0]:.4f}·H_sub + {c[1]:.4f})")
    print(f"     = {2*c[0]:.4f}·Δc₃·H_sub + {2*c[1]:.4f}·Δc₃")
    print(f"  → Predicted a_inter from 3-cycles alone: {2*c[0]:.6f}")
    print(f"  → Measured a_inter at n=5→6: ~0.278")

    # By H_sub level
    by_H = defaultdict(list)
    for d in data:
        by_H[d['H_sub']].append(d['mu'])

    print(f"\n  Mean μ by H_sub:")
    for H in sorted(by_H.keys()):
        vals = by_H[H]
        print(f"    H_sub={H:3d}: mean μ = {np.mean(vals):.4f} ± {np.std(vals):.4f} (n={len(vals)})")

    return data


if __name__ == "__main__":
    t0 = time.time()

    # 1. Claim A decomposition at n=6
    claim_a_decomp_n6()
    t1 = time.time()
    print(f"\n  [Claim A n=6: {t1-t0:.1f}s]")

    # 2. a_inter mechanism
    a_inter_mechanism()
    t2 = time.time()
    print(f"\n  [a_inter: {t2-t1:.1f}s]")

    # 3. β_c at n=7
    beta_c_7 = beta_c_n7()
    t3 = time.time()
    print(f"\n  [β_c n=7: {t3-t2:.1f}s]")

    # Summary
    print(f"\n{'='*70}")
    print(f"  PHASE TRANSITION SCALING:")
    print(f"  n=5: β_c ≈ 0.700")
    print(f"  n=6: β_c ≈ 0.310")
    print(f"  n=7: β_c ≈ {beta_c_7:.4f}")
    print(f"  Pattern: β_c decreasing. Product β_c·max_H: {0.7*15:.1f}, {0.31*45:.1f}, {beta_c_7*189:.1f}")
    print(f"{'='*70}")
    print(f"\n  TOTAL: {t3-t0:.1f}s")
