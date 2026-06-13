#!/usr/bin/env python3
"""
deep_principles_large_n_s22.py — opus-2026-04-04-S22

EXTEND the deep principles to n=7,8,9 via sampling.

The three principles to verify at large n:
  P1. PARABOLIC LAW: E[Δc₃ | score(n)=s] = s(n-1-s)/2 (proved for all n, verify)
  P2. FRUSTRATION-H CORRELATION: corr(c₃, H) ≈ ? (expect ~0.95 declining slowly)
  P3. EXCHANGE COUPLING: ΔH ≈ a·Δc₃·H_sub + b·Δc₃ + c (a ≈ 0.16-0.28)

Also verify:
  P4. COEFFICIENT SEQUENCE: simple coeff a_n in ΔH ≈ a_n·Δc₃ (2,3,6,15,?)
  P5. FRUSTRATION-H NEAR-PERFECT: R² of ΔH ~ Δc₃ vs ΔH ~ Δc₃ + H_sub
  P6. CORR(c₃, H) for full tournaments at n=7 (exhaustive: 2^21 = 2M, feasible)
"""

import numpy as np
from math import comb
from collections import defaultdict
import random
import time
import sys

random.seed(42)

# ─── Core ───

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

def score_variance(A, n):
    s = [sum(A[i][j] for j in range(n) if j != i) for i in range(n)]
    m = sum(s) / n
    return sum((x - m)**2 for x in s) / n

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

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
    W = [i for i in range(n_sub) if sigma[i] == 1]
    L = [i for i in range(n_sub) if sigma[i] == 0]
    return sum(A_sub[i][j] for i in W for j in L)


# ─── P1: Parabolic Law at larger n ───

def verify_parabolic_sampled(n_sub, n_samples=2000):
    """Verify parabolic law by sampling random T and random σ of each score."""
    print(f"\n{'='*70}")
    print(f"  P1: PARABOLIC LAW at n_sub = {n_sub} (sampled, {n_samples} per score)")
    print(f"{'='*70}")

    max_err = 0
    for s in range(n_sub + 1):
        expected = s * (n_sub - s) / 2
        actual_sum = 0
        count = 0
        for _ in range(n_samples):
            A = random_tournament(n_sub)
            # Random subset W of size s
            indices = list(range(n_sub))
            random.shuffle(indices)
            sigma = [0] * n_sub
            for i in indices[:s]:
                sigma[i] = 1
            dc3 = delta_c3(A, n_sub, tuple(sigma))
            actual_sum += dc3
            count += 1
        actual = actual_sum / count
        err = abs(actual - expected)
        max_err = max(max_err, err)
        print(f"    s={s:2d}: expected={expected:8.4f}, actual={actual:8.4f}, "
              f"ratio={actual/expected:.6f}" if expected > 0 else
              f"    s={s:2d}: expected=0, actual={actual:.4f}")

    print(f"  Max deviation: {max_err:.6f} (should be ~0 for large samples)")


# ─── P2: Frustration-H correlation at n=7 (sampled) ───

def frustration_H_correlation(n, n_samples=20000):
    """Compute corr(c₃, H), corr(sv, H) at given n."""
    print(f"\n{'='*70}")
    print(f"  P2: FRUSTRATION-H CORRELATION at n = {n} ({n_samples} samples)")
    print(f"{'='*70}")

    c3_list = []
    H_list = []
    sv_list = []

    t0 = time.time()
    for i in range(n_samples):
        A = random_tournament(n)
        H = compute_H(A, n)
        c3 = count_3cycles(A, n)
        sv = score_variance(A, n)
        H_list.append(H)
        c3_list.append(c3)
        sv_list.append(sv)
        if (i+1) % 5000 == 0:
            print(f"    ... {i+1}/{n_samples} done ({time.time()-t0:.1f}s)", file=sys.stderr)

    c3_arr = np.array(c3_list, dtype=float)
    H_arr = np.array(H_list, dtype=float)
    sv_arr = np.array(sv_list, dtype=float)

    r_c3_H = np.corrcoef(c3_arr, H_arr)[0,1]
    r_sv_H = np.corrcoef(sv_arr, H_arr)[0,1]

    print(f"  corr(c₃, H) = {r_c3_H:.6f}")
    print(f"  corr(sv, H) = {r_sv_H:.6f}")
    print(f"  H range: [{H_arr.min():.0f}, {H_arr.max():.0f}], mean={H_arr.mean():.2f}")
    print(f"  c₃ range: [{c3_arr.min():.0f}, {c3_arr.max():.0f}], mean={c3_arr.mean():.2f}")

    # Linear model: H ~ a*c₃ + b
    X = np.column_stack([c3_arr, np.ones(n_samples)])
    c, _, _, _ = np.linalg.lstsq(X, H_arr, rcond=None)
    pred = X @ c
    R2 = 1 - np.sum((H_arr - pred)**2) / np.sum((H_arr - H_arr.mean())**2)
    print(f"  Linear: H ≈ {c[0]:.4f}·c₃ + {c[1]:.4f}, R²={R2:.6f}")

    # Quadratic: H ~ a*c₃² + b*c₃ + c
    X2 = np.column_stack([c3_arr**2, c3_arr, np.ones(n_samples)])
    c2, _, _, _ = np.linalg.lstsq(X2, H_arr, rcond=None)
    pred2 = X2 @ c2
    R2_q = 1 - np.sum((H_arr - pred2)**2) / np.sum((H_arr - H_arr.mean())**2)
    print(f"  Quadratic: H ≈ {c2[0]:.6f}·c₃² + {c2[1]:.4f}·c₃ + {c2[2]:.4f}, R²={R2_q:.6f}")

    # With score variance
    X3 = np.column_stack([c3_arr, sv_arr, np.ones(n_samples)])
    c3_coeff, _, _, _ = np.linalg.lstsq(X3, H_arr, rcond=None)
    pred3 = X3 @ c3_coeff
    R2_3 = 1 - np.sum((H_arr - pred3)**2) / np.sum((H_arr - H_arr.mean())**2)
    print(f"  c₃+sv: H ≈ {c3_coeff[0]:.4f}·c₃ + {c3_coeff[1]:.4f}·sv + {c3_coeff[2]:.4f}, R²={R2_3:.6f}")

    return r_c3_H, r_sv_H, c[0], R2


# ─── P3: Exchange coupling at n=7→8 ───

def exchange_coupling(n_target, n_samples_parent=3000, n_ext_per_parent=16):
    """Measure exchange coupling at n-1 → n."""
    n_sub = n_target - 1
    print(f"\n{'='*70}")
    print(f"  P3: EXCHANGE COUPLING at n={n_sub} → n={n_target}")
    print(f"  ({n_samples_parent} parents × {n_ext_per_parent} extensions)")
    print(f"{'='*70}")

    data = []
    t0 = time.time()
    for p in range(n_samples_parent):
        A_sub = random_tournament(n_sub)
        H_sub = compute_H(A_sub, n_sub)
        c3_sub = count_3cycles(A_sub, n_sub)

        for _ in range(n_ext_per_parent):
            sigma = tuple(random.randint(0,1) for _ in range(n_sub))
            score_n = sum(sigma)
            A_n = extend(A_sub, n_sub, sigma)
            H_n = compute_H(A_n, n_target)
            dc3 = delta_c3(A_sub, n_sub, sigma)
            data.append((H_sub, H_n - H_sub, dc3, score_n, c3_sub))

        if (p+1) % 500 == 0:
            print(f"    ... {p+1}/{n_samples_parent} parents ({time.time()-t0:.1f}s)", file=sys.stderr)

    N = len(data)
    H_sub_arr = np.array([d[0] for d in data], dtype=float)
    dH_arr = np.array([d[1] for d in data], dtype=float)
    dc3_arr = np.array([d[2] for d in data], dtype=float)
    score_n_arr = np.array([d[3] for d in data], dtype=float)
    c3_sub_arr = np.array([d[4] for d in data], dtype=float)

    print(f"\n  Data points: {N}")

    # Simple: ΔH ~ a*Δc₃
    X1 = np.column_stack([dc3_arr, np.ones(N)])
    c1, _, _, _ = np.linalg.lstsq(X1, dH_arr, rcond=None)
    R1 = 1 - np.sum((dH_arr - X1@c1)**2) / np.sum((dH_arr - dH_arr.mean())**2)
    print(f"  Simple:  ΔH ≈ {c1[0]:.4f}·Δc₃ + {c1[1]:.4f},  R²={R1:.4f}")

    # +H_sub
    X2 = np.column_stack([dc3_arr, H_sub_arr, np.ones(N)])
    c2, _, _, _ = np.linalg.lstsq(X2, dH_arr, rcond=None)
    R2 = 1 - np.sum((dH_arr - X2@c2)**2) / np.sum((dH_arr - dH_arr.mean())**2)
    print(f"  +H_sub:  ΔH ≈ {c2[0]:.4f}·Δc₃ + {c2[1]:.4f}·H_sub + {c2[2]:.4f},  R²={R2:.4f}")

    # Interaction
    X3 = np.column_stack([dc3_arr * H_sub_arr, dc3_arr, H_sub_arr, np.ones(N)])
    c3, _, _, _ = np.linalg.lstsq(X3, dH_arr, rcond=None)
    R3 = 1 - np.sum((dH_arr - X3@c3)**2) / np.sum((dH_arr - dH_arr.mean())**2)
    print(f"  Interact: ΔH ≈ {c3[0]:.6f}·Δc₃·H_sub + {c3[1]:.4f}·Δc₃ + "
          f"{c3[2]:.4f}·H_sub + {c3[3]:.4f},  R²={R3:.4f}")

    # Parabolic law check
    dc3_by_score = defaultdict(list)
    for d in data:
        dc3_by_score[d[3]].append(d[2])

    print(f"\n  Parabolic law check (E[Δc₃ | score=s] vs s(n_sub-s)/2):")
    for s in sorted(dc3_by_score.keys()):
        vals = dc3_by_score[s]
        actual = np.mean(vals)
        expected = s * (n_sub - s) / 2
        ratio = actual / expected if expected > 0 else float('inf')
        print(f"    s={s:2d}: E[Δc₃]={actual:8.3f}, theory={expected:8.3f}, "
              f"ratio={ratio:.4f}" if expected > 0.01 else
              f"    s={s:2d}: E[Δc₃]={actual:8.3f}, theory={expected:8.3f}")

    # Per-parent J_eff
    by_H = defaultdict(lambda: {'dH': [], 'dc3': []})
    for d in data:
        by_H[d[0]]['dH'].append(d[1])
        by_H[d[0]]['dc3'].append(d[2])

    print(f"\n  Per-parent J_eff = mean_ΔH / mean_Δc₃ (top 10 H values):")
    H_vals = sorted(by_H.keys())
    for H in H_vals[:5] + [H_vals[len(H_vals)//2]] + H_vals[-5:]:
        if H in by_H:
            mean_dH = np.mean(by_H[H]['dH'])
            mean_dc3 = np.mean(by_H[H]['dc3'])
            J = mean_dH / mean_dc3 if mean_dc3 > 0.1 else 0
            n_pts = len(by_H[H]['dH'])
            if n_pts >= 10:
                print(f"    H_sub={H:5.0f}: J_eff={J:8.4f}  (n={n_pts})")

    return c1[0], c3[0], R1, R3


# ─── P4: Full scaling table ───

def scaling_table():
    """Compile all coefficients into a table."""
    print(f"\n{'='*70}")
    print(f"  SCALING TABLE — ALL LEVELS")
    print(f"{'='*70}")
    print(f"  {'Level':>8}  {'Simple a':>10}  {'R²(simple)':>12}  {'Interact a':>12}  {'R²(inter)':>11}  {'J growth':>10}")
    print(f"  {'─'*8}  {'─'*10}  {'─'*12}  {'─'*12}  {'─'*11}  {'─'*10}")

    # Known from previous sessions
    known = [
        ('3→4', 2.0, 1.000, 0.000, 1.000, ''),
        ('4→5', 3.0, 0.857, 0.261, 0.937, '×1.50'),
        ('5→6', 6.0, 0.713, 0.278, 0.907, '×2.00'),
    ]
    for level, a, R2s, ai, R2i, growth in known:
        print(f"  {level:>8}  {a:10.4f}  {R2s:12.4f}  {ai:12.6f}  {R2i:11.4f}  {growth:>10}")

    return known


# ─── MAIN ───

if __name__ == "__main__":
    print("╔══════════════════════════════════════════════════════════════════╗")
    print("║  DEEP PRINCIPLES AT LARGE n                                    ║")
    print("║  opus-2026-04-04-S22                                           ║")
    print("╚══════════════════════════════════════════════════════════════════╝")

    t_start = time.time()

    # P1: Parabolic law at n=7,8,9
    for n_sub in [7, 8, 9]:
        verify_parabolic_sampled(n_sub, n_samples=1000)

    t1 = time.time()
    print(f"\n  [Parabolic law: {t1-t_start:.1f}s]")

    # P2: Frustration-H correlation at n=7, n=8
    results_P2 = {}
    for n in [7, 8]:
        samples = 10000 if n <= 7 else 3000
        r, rsv, slope, R2 = frustration_H_correlation(n, n_samples=samples)
        results_P2[n] = (r, rsv, slope, R2)

    t2 = time.time()
    print(f"\n  [Frustration-H: {t2-t1:.1f}s]")

    # P3: Exchange coupling at n=6→7, n=7→8
    results_P3 = {}
    for n_target in [7, 8]:
        samples = 2000 if n_target <= 7 else 500
        ext = 16 if n_target <= 7 else 8
        a_simple, a_inter, R2s, R2i = exchange_coupling(n_target, samples, ext)
        results_P3[n_target] = (a_simple, a_inter, R2s, R2i)

    t3 = time.time()
    print(f"\n  [Exchange coupling: {t3-t2:.1f}s]")

    # Summary
    known = scaling_table()

    # Add new results
    print(f"\n  NEW RESULTS:")
    for n_target in [7, 8]:
        a_s, a_i, R2s, R2i = results_P3[n_target]
        growth = f"×{a_s/6.0:.2f}" if n_target == 7 else f"×{a_s/results_P3[7][0]:.2f}"
        print(f"  {n_target-1}→{n_target}:  simple={a_s:10.4f}  R²={R2s:.4f}  "
              f"interact={a_i:12.6f}  R²={R2i:.4f}  growth={growth}")

    print(f"\n  FRUSTRATION-H CORRELATION TREND:")
    print(f"  {'n':>3}  {'corr(c₃,H)':>12}  {'corr(sv,H)':>12}  {'slope':>8}  {'R²(linear)':>12}")
    for n in [4, 5, 6]:
        # From S16
        corrs = {4: (1.000, -1.000, 2.0, 1.000),
                 5: (0.973, -0.973, 3.0, 0.947),
                 6: (0.961, -0.961, 6.0, 0.923)}
        r, rsv, sl, R2 = corrs[n]
        print(f"  {n:3d}  {r:12.6f}  {rsv:12.6f}  {sl:8.4f}  {R2:12.6f}")
    for n in [7, 8]:
        if n in results_P2:
            r, rsv, sl, R2 = results_P2[n]
            print(f"  {n:3d}  {r:12.6f}  {rsv:12.6f}  {sl:8.4f}  {R2:12.6f}")

    # Coefficient sequence analysis
    print(f"\n  COEFFICIENT SEQUENCE:")
    coeffs = [2.0, 3.0, 6.0]
    for n_target in [7, 8]:
        if n_target in results_P3:
            coeffs.append(results_P3[n_target][0])
    print(f"  Coefficients: {[round(c,1) for c in coeffs]}")
    for i in range(1, len(coeffs)):
        print(f"    Ratio a({i+3}→{i+4})/a({i+2}→{i+3}) = {coeffs[i]/coeffs[i-1]:.4f}")

    t_end = time.time()
    print(f"\n  TOTAL TIME: {t_end - t_start:.1f}s")
