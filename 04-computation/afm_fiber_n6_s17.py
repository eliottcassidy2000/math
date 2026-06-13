#!/usr/bin/env python3
"""
afm_fiber_n6_s17.py — opus-2026-04-04-S17

Extend fiber bundle AFM analysis to n=5→6.
Focus: H propagation formula and per-class magnon spectrum at n=6.
"""
import numpy as np
from itertools import combinations
from math import comb
from collections import defaultdict, Counter
import time
import sys
import random

random.seed(42)

def compute_H(A, n):
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

def scores(A, n):
    return tuple(sorted(sum(A[i][j] for j in range(n) if j != i) for i in range(n)))

def count_3cycles(A, n):
    c = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]:
                    c += 1
                elif A[i][k] and A[k][j] and A[j][i]:
                    c += 1
    return c

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

def extend_tournament(A_sub, n_sub, sigma):
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
    W = [i for i in range(n_sub) if sigma[i] == 1]
    L = [i for i in range(n_sub) if sigma[i] == 0]
    arcs = 0
    for i in W:
        for j in L:
            if A_sub[i][j]:
                arcs += 1
    return arcs


def h_propagation_n6():
    """H propagation from n=5 to n=6, sampled."""
    n_sub = 5
    n_target = 6
    print(f"\n{'='*70}")
    print(f"  H PROPAGATION: n=5 → n=6 (exhaustive over parents, full extensions)")
    print(f"{'='*70}")

    all_data = []
    count = 0
    for A_sub in all_tournaments(n_sub):
        H_sub = compute_H(A_sub, n_sub)
        c3_sub = count_3cycles(A_sub, n_sub)

        for sigma_int in range(1 << n_sub):
            sigma = tuple((sigma_int >> i) & 1 for i in range(n_sub))
            score_n = sum(sigma)
            A_n = extend_tournament(A_sub, n_sub, sigma)
            H_n = compute_H(A_n, n_target)
            dc3 = delta_c3(A_sub, n_sub, sigma)
            dH = H_n - H_sub

            all_data.append({
                'H_sub': H_sub, 'H_n': H_n, 'dH': dH,
                'c3_sub': c3_sub, 'dc3': dc3,
                'score_n': score_n,
            })

        count += 1
        if count % 200 == 0:
            print(f"  ... {count}/1024 parents done", file=sys.stderr)

    N = len(all_data)
    all_dH = np.array([d['dH'] for d in all_data])
    all_dc3 = np.array([d['dc3'] for d in all_data])
    all_score_n = np.array([d['score_n'] for d in all_data])
    all_H_sub = np.array([d['H_sub'] for d in all_data])
    all_c3_sub = np.array([d['c3_sub'] for d in all_data])

    print(f"\n  Total extensions: {N}")

    # ΔH by Δc₃
    dH_by_dc3 = defaultdict(list)
    for d in all_data:
        dH_by_dc3[d['dc3']].append(d['dH'])

    print(f"\n  ΔH by Δc₃:")
    print(f"  {'Δc₃':>5}  {'mean ΔH':>10}  {'std':>10}  {'min':>6}  {'max':>6}  {'count':>7}")
    for dc3 in sorted(dH_by_dc3.keys()):
        vals = dH_by_dc3[dc3]
        print(f"  {dc3:5d}  {np.mean(vals):10.4f}  {np.std(vals):10.4f}  "
              f"{min(vals):6d}  {max(vals):6d}  {len(vals):7d}")

    # Correlation
    corr = np.corrcoef(all_dc3, all_dH)[0,1]
    print(f"\n  Correlation(Δc₃, ΔH) = {corr:.6f}")

    # Full linear model
    X = np.column_stack([all_dc3, all_score_n, all_H_sub, all_c3_sub, np.ones(N)])
    coeffs, _, _, _ = np.linalg.lstsq(X, all_dH, rcond=None)
    pred = X @ coeffs
    R2 = 1 - np.sum((all_dH - pred)**2) / np.sum((all_dH - all_dH.mean())**2)
    print(f"\n  Full model: ΔH ≈ {coeffs[0]:.4f}*Δc₃ + {coeffs[1]:.4f}*score(n) "
          f"+ {coeffs[2]:.4f}*H_sub + {coeffs[3]:.4f}*c₃_sub + {coeffs[4]:.4f}")
    print(f"  R² = {R2:.6f}")

    # Simple: ΔH ~ a*Δc₃ + b
    X2 = np.column_stack([all_dc3, np.ones(N)])
    coeffs2, _, _, _ = np.linalg.lstsq(X2, all_dH, rcond=None)
    pred2 = X2 @ coeffs2
    R2_simple = 1 - np.sum((all_dH - pred2)**2) / np.sum((all_dH - all_dH.mean())**2)
    print(f"  Simple: ΔH ≈ {coeffs2[0]:.4f}*Δc₃ + {coeffs2[1]:.4f}")
    print(f"  R²(simple) = {R2_simple:.6f}")

    # ΔH ~ a*Δc₃ + b*H_sub + c
    X3 = np.column_stack([all_dc3, all_H_sub, np.ones(N)])
    coeffs3, _, _, _ = np.linalg.lstsq(X3, all_dH, rcond=None)
    pred3 = X3 @ coeffs3
    R2_3 = 1 - np.sum((all_dH - pred3)**2) / np.sum((all_dH - all_dH.mean())**2)
    print(f"  ΔH ≈ {coeffs3[0]:.4f}*Δc₃ + {coeffs3[1]:.4f}*H_sub + {coeffs3[2]:.4f}")
    print(f"  R²(Δc₃ + H_sub) = {R2_3:.6f}")

    # Test: ΔH ~ a*Δc₃*H_sub + b*Δc₃ + c (multiplicative interaction)
    X4 = np.column_stack([all_dc3 * all_H_sub, all_dc3, all_H_sub, np.ones(N)])
    coeffs4, _, _, _ = np.linalg.lstsq(X4, all_dH, rcond=None)
    pred4 = X4 @ coeffs4
    R2_4 = 1 - np.sum((all_dH - pred4)**2) / np.sum((all_dH - all_dH.mean())**2)
    print(f"  ΔH ≈ {coeffs4[0]:.6f}*Δc₃*H_sub + {coeffs4[1]:.4f}*Δc₃ + "
          f"{coeffs4[2]:.4f}*H_sub + {coeffs4[3]:.4f}")
    print(f"  R²(interaction) = {R2_4:.6f}")

    # Fiber connection per parent H
    dH_by_Hsub = defaultdict(list)
    dc3_by_Hsub = defaultdict(list)
    for d in all_data:
        dH_by_Hsub[d['H_sub']].append(d['dH'])
        dc3_by_Hsub[d['H_sub']].append(d['dc3'])

    print(f"\n  FIBER CONNECTION by parent H:")
    print(f"  {'H_sub':>6}  {'mean ΔH':>10}  {'mean Δc₃':>10}  {'ΔH/Δc₃':>8}  "
          f"{'corr':>8}  {'min ΔH':>7}  {'max ΔH':>7}")
    for H_sub in sorted(dH_by_Hsub.keys()):
        dH_vals = np.array(dH_by_Hsub[H_sub])
        dc3_vals = np.array(dc3_by_Hsub[H_sub])
        mean_dH = dH_vals.mean()
        mean_dc3 = dc3_vals.mean()
        ratio = mean_dH / mean_dc3 if mean_dc3 > 0 else float('nan')
        c = np.corrcoef(dc3_vals, dH_vals)[0,1] if np.std(dc3_vals) > 0 else 0
        print(f"  {H_sub:6d}  {mean_dH:10.4f}  {mean_dc3:10.4f}  {ratio:8.4f}  "
              f"{c:8.4f}  {dH_vals.min():7d}  {dH_vals.max():7d}")

    return all_data


def magnon_decomp_n6():
    """Inner vs boundary magnon decomposition at n=6 (sampled)."""
    n = 6
    new_v = n - 1  # vertex 5 = the "extension vertex"
    print(f"\n{'='*70}")
    print(f"  INNER vs BOUNDARY MAGNON DECOMPOSITION — n = {n} (sampled)")
    print(f"{'='*70}")

    inner_dH = []
    boundary_dH = []
    boundary_by_v = defaultdict(list)

    sample_size = 5000
    for _ in range(sample_size):
        A = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
        H0 = compute_H(A, n)

        for i in range(n):
            for j in range(i+1, n):
                A2 = [row[:] for row in A]
                A2[i][j] = 1 - A2[i][j]
                A2[j][i] = 1 - A2[j][i]
                H1 = compute_H(A2, n)
                dH = abs(H1 - H0)

                if i == new_v or j == new_v:
                    boundary_dH.append(dH)
                    other = j if i == new_v else i
                    boundary_by_v[other].append(dH)
                else:
                    inner_dH.append(dH)

    inner_arr = np.array(inner_dH)
    boundary_arr = np.array(boundary_dH)

    print(f"\n  Samples: {sample_size}")
    print(f"\n  INNER magnons (flip within T_5):")
    print(f"    Count: {len(inner_arr)}")
    print(f"    Mean |ΔH|: {inner_arr.mean():.4f}")
    print(f"    ΔH=0 rate: {100*np.sum(inner_arr==0)/len(inner_arr):.2f}%")

    print(f"\n  BOUNDARY magnons (flip arc to/from vertex 5):")
    print(f"    Count: {len(boundary_arr)}")
    print(f"    Mean |ΔH|: {boundary_arr.mean():.4f}")
    print(f"    ΔH=0 rate: {100*np.sum(boundary_arr==0)/len(boundary_arr):.2f}%")

    print(f"\n  RATIO: boundary/inner = {boundary_arr.mean()/inner_arr.mean():.4f}")

    # Boundary by vertex
    print(f"\n  Boundary by vertex:")
    for v in sorted(boundary_by_v.keys()):
        vals = np.array(boundary_by_v[v])
        print(f"    v={v}: mean|ΔH|={vals.mean():.4f}")

    return inner_dH, boundary_dH


if __name__ == "__main__":
    t0 = time.time()

    # H propagation n=5→6 (exhaustive: 1024 parents × 32 extensions = 32768)
    h_propagation_n6()

    t1 = time.time()
    print(f"\n  [H propagation n=5→6: {t1-t0:.1f}s]")

    # Magnon decomposition n=6 (sampled)
    magnon_decomp_n6()

    t2 = time.time()
    print(f"\n  [Magnon decomp n=6: {t2-t1:.1f}s]")
    print(f"\n  TOTAL: {t2-t0:.1f}s")
