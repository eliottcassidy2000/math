#!/usr/bin/env python3
"""
prove_coefficient_s23.py — opus-2026-04-04-S23

PROVE a(n) = (n-2)!/2^{n-4} via Claim A cycle-length decomposition.

CLAIM A: ΔH = H(T_n) - H(T_{n-1}) = 2 Σ_{C∋n, C odd} μ(C)
where μ(C) = I(Ω(T_{n-1})|_{avoid C\{n}}, 2) = H(T[V \ V(C)])

So: ΔH = 2 Σ_{k=3,5,7,...} Σ_{C∋n, |C|=k} H(T[V(T_n) \ V(C)])

For a k-cycle C through n: V(C) has k vertices, so V(T_n)\V(C) has n-k vertices.
  T[V(T_n)\V(C)] is a subtournament on n-k vertices.
  μ(C) = H(T[n-k vertices]).

DECOMPOSITION:
  ΔH = 2 Σ_{k=3} #{k-cycles thru n} · (average μ over those cycles)
     + 2 Σ_{k=5} #{k-cycles thru n} · (average μ)
     + ...

KEY OBSERVATIONS:
  - k=3: μ = H(T[n-3 verts]). These are subtournaments of T_{n-1} on n-3 vertices.
  - k=5: μ = H(T[n-5 verts]).
  - k=n (if n odd): μ = H(T[0 verts]) = 1 (empty tournament has 1 trivial path).
  - k=n-1 (if n-1 odd): μ = H(T[1 vert]) = 1.

So the cycle-length-k contribution to ΔH is:
  C_k = 2 · (#{k-cycles through n}) · E[H(T[n-k vertices]) | those cycles]

The TOTAL:
  ΔH = Σ_{k=3,5,...} C_k

We want to show: E[ΔH] / E[Δc₃] = (n-2)!/2^{n-4}
where the expectations are over uniform random (T_{n-1}, σ) extensions.

PLAN:
  1. Compute C_k exhaustively at n=5 (from n=4) and n=6 (from n=5)
  2. Find the contribution of each cycle length to the total coefficient
  3. Look for the recursive/factorial pattern
"""

import numpy as np
from itertools import combinations, permutations
from math import comb, factorial
from collections import defaultdict
import time
import sys

# ─── Core ───

def compute_H(A, n):
    if n == 0:
        return 1  # Convention: empty tournament has 1 path
    if n == 1:
        return 1
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
    """Extract subtournament on given vertices (0-indexed)."""
    verts = sorted(vertices)
    k = len(verts)
    B = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            B[i][j] = A[verts[i]][verts[j]]
    return B, k

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

def find_directed_odd_cycles_through_v(A, n, v):
    """Find all directed odd cycles through vertex v in tournament A.
    Returns list of (vertex_set, cycle_length)."""
    cycles = []
    # For each odd length k from 3 to n
    for k in range(3, n+1, 2):
        # Choose k-1 other vertices
        other_verts = [u for u in range(n) if u != v]
        for combo in combinations(other_verts, k-1):
            verts = [v] + list(combo)
            # Check all orderings starting from v
            for perm in permutations(combo):
                cycle = [v] + list(perm)
                if all(A[cycle[i]][cycle[(i+1) % k]] for i in range(k)):
                    cycles.append((frozenset(cycle), k))
                    break  # Only need one valid direction per vertex set
            # Actually, need to count distinct directed cycles
            # A directed cycle on vertex set S through v: fix v as starting point,
            # check all (k-1)! orderings of the others
            # But frozenset deduplicates by vertex set — need to be more careful

    # Better approach: enumerate by vertex set, count directed cycles
    result = []
    for k in range(3, n+1, 2):
        other_verts = [u for u in range(n) if u != v]
        for combo in combinations(other_verts, k-1):
            verts_set = frozenset([v] + list(combo))
            # Count directed k-cycles on this vertex set containing v
            # Fix v as position 0, try all permutations of the rest
            for perm in permutations(combo):
                cycle = [v] + list(perm)
                if all(A[cycle[i]][cycle[(i+1) % k]] for i in range(k)):
                    result.append((verts_set, k, tuple(cycle)))
                    # Each directed cycle is counted once (since we fix v as start)
    return result


def claim_a_decomposition(n_target):
    """Decompose ΔH = Σ_k C_k where C_k is the k-cycle contribution."""
    n_sub = n_target - 1
    new_v = n_target - 1  # 0-indexed

    print(f"\n{'='*70}")
    print(f"  CLAIM A DECOMPOSITION: n = {n_sub} → {n_target}")
    print(f"{'='*70}")

    # For each (T_{n-1}, σ), compute ΔH decomposed by cycle length
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

            # Find all directed odd cycles through new_v in A_n
            cycles = find_directed_odd_cycles_through_v(A_n, n_target, new_v)

            # For each cycle, compute μ(C) = H(T[V \ V(C)])
            contribution_by_k = defaultdict(int)
            for verts_set, k, cycle_tuple in cycles:
                complement_verts = [u for u in range(n_target) if u not in verts_set]
                if len(complement_verts) == 0:
                    mu = 1
                else:
                    B, m = subtournament(A_n, n_target, complement_verts)
                    mu = compute_H(B, m)
                contribution_by_k[k] += mu  # Will multiply by 2 below

            # Verify Claim A
            claim_a_sum = 2 * sum(contribution_by_k.values())
            if claim_a_sum != dH:
                # Could be cycle counting issue — let me debug
                pass

            total_data.append({
                'dH': dH, 'dc3': dc3, 'H_sub': H_sub,
                'by_k': dict(contribution_by_k),
                'claim_a_sum': claim_a_sum,
            })

        count += 1
        if count % 50 == 0:
            print(f"  ... {count} parents done", file=sys.stderr)

    N = len(total_data)

    # Check Claim A verification
    mismatches = sum(1 for d in total_data if d['claim_a_sum'] != d['dH'])
    print(f"\n  Claim A mismatches: {mismatches}/{N}")
    if mismatches > 0:
        for d in total_data[:5]:
            if d['claim_a_sum'] != d['dH']:
                print(f"    dH={d['dH']}, claim_a={d['claim_a_sum']}, by_k={d['by_k']}")

    # Decompose the coefficient by cycle length
    all_dc3 = np.array([d['dc3'] for d in total_data], dtype=float)
    all_dH = np.array([d['dH'] for d in total_data], dtype=float)

    # For each cycle length k, compute the regression coefficient
    max_k = n_target if n_target % 2 == 1 else n_target - 1
    print(f"\n  Contribution by cycle length k:")
    print(f"  {'k':>3}  {'mean 2·Σμ':>12}  {'as fraction of ΔH':>20}  {'coeff vs Δc₃':>15}")

    total_coeff = 0
    for k in range(3, max_k + 1, 2):
        # Extract k-cycle contribution: 2 * sum of μ values for k-cycles
        Ck = np.array([2 * d['by_k'].get(k, 0) for d in total_data], dtype=float)
        mean_Ck = Ck.mean()
        mean_dH = all_dH.mean()
        frac = mean_Ck / mean_dH if mean_dH != 0 else 0

        # Regression Ck vs Δc₃
        if np.std(all_dc3) > 0 and np.std(Ck) > 0:
            X = np.column_stack([all_dc3, np.ones(N)])
            c, _, _, _ = np.linalg.lstsq(X, Ck, rcond=None)
            coeff = c[0]
        else:
            coeff = 0

        total_coeff += coeff
        print(f"  {k:3d}  {mean_Ck:12.6f}  {frac:20.6f}  {coeff:15.6f}")

    print(f"  {'SUM':>3}  {all_dH.mean():12.6f}  {'1.000000':>20}  {total_coeff:15.6f}")

    # The predicted coefficient
    predicted = factorial(n_target - 2) / (2 ** (n_target - 4))
    print(f"\n  Total coefficient from decomposition: {total_coeff:.6f}")
    print(f"  Predicted (n-2)!/2^(n-4): {predicted:.6f}")
    print(f"  Ratio: {total_coeff/predicted:.6f}")

    # Key per-k analysis: coefficient of k-cycle contribution vs Δc₃
    print(f"\n  RECURSIVE STRUCTURE:")
    print(f"  k=3 contribution coefficient = coeff of 3-cycle μ vs Δc₃")
    print(f"  This should relate to H at (n-3) vertices times the 3-cycle count")

    return total_data


def variance_decomposition(n):
    """Decompose Var(H) into contributions from c₃, c₅, α₂, etc.
    to understand why R² ≈ 0.914."""
    print(f"\n{'='*70}")
    print(f"  VARIANCE DECOMPOSITION of H at n = {n}")
    print(f"{'='*70}")

    data = []
    for A in all_tournaments(n):
        H = compute_H(A, n)
        c3 = count_3cycles(A, n)

        # Count 5-cycles
        c5 = 0
        if n >= 5:
            for combo in combinations(range(n), 5):
                verts = list(combo)
                for perm in permutations(verts):
                    if perm[0] == min(perm):  # canonical start
                        if all(A[perm[i]][perm[(i+1) % 5]] for i in range(5)):
                            c5 += 1

        # Count independent pairs (α₂) — pairs of vertex-disjoint 3-cycles
        # This is expensive but doable at small n
        three_cycles = []
        for i in range(n):
            for j in range(i+1, n):
                for k in range(j+1, n):
                    if A[i][j] and A[j][k] and A[k][i]:
                        three_cycles.append(frozenset([i,j,k]))
                    elif A[i][k] and A[k][j] and A[j][i]:
                        three_cycles.append(frozenset([i,j,k]))

        alpha2 = 0
        for a in range(len(three_cycles)):
            for b in range(a+1, len(three_cycles)):
                if not (three_cycles[a] & three_cycles[b]):  # disjoint
                    alpha2 += 1

        data.append({'H': H, 'c3': c3, 'c5': c5, 'alpha2': alpha2})

    N = len(data)
    H_arr = np.array([d['H'] for d in data], dtype=float)
    c3_arr = np.array([d['c3'] for d in data], dtype=float)
    c5_arr = np.array([d['c5'] for d in data], dtype=float)
    a2_arr = np.array([d['alpha2'] for d in data], dtype=float)

    print(f"  {N} tournaments")

    # Correlations
    print(f"\n  Correlations:")
    for name, arr in [('c₃', c3_arr), ('c₅', c5_arr), ('α₂', a2_arr)]:
        if np.std(arr) > 0:
            r = np.corrcoef(arr, H_arr)[0,1]
            print(f"    corr({name}, H) = {r:.6f}")

    # Cross-correlations
    if np.std(c5_arr) > 0:
        print(f"    corr(c₃, c₅) = {np.corrcoef(c3_arr, c5_arr)[0,1]:.6f}")
    if np.std(a2_arr) > 0:
        print(f"    corr(c₃, α₂) = {np.corrcoef(c3_arr, a2_arr)[0,1]:.6f}")

    # R² from various models
    # c₃ only
    X1 = np.column_stack([c3_arr, np.ones(N)])
    c1, _, _, _ = np.linalg.lstsq(X1, H_arr, rcond=None)
    R2_1 = 1 - np.sum((H_arr - X1@c1)**2) / np.sum((H_arr - H_arr.mean())**2)
    print(f"\n  R²(c₃ only) = {R2_1:.6f}")

    # c₃ + c₅
    if np.std(c5_arr) > 0:
        X2 = np.column_stack([c3_arr, c5_arr, np.ones(N)])
        c2, _, _, _ = np.linalg.lstsq(X2, H_arr, rcond=None)
        R2_2 = 1 - np.sum((H_arr - X2@c2)**2) / np.sum((H_arr - H_arr.mean())**2)
        print(f"  R²(c₃ + c₅) = {R2_2:.6f} (ΔR² = {R2_2 - R2_1:.6f})")
        print(f"    H ≈ {c2[0]:.4f}·c₃ + {c2[1]:.4f}·c₅ + {c2[2]:.4f}")

    # c₃ + α₂
    if np.std(a2_arr) > 0:
        X3 = np.column_stack([c3_arr, a2_arr, np.ones(N)])
        c3_coeff, _, _, _ = np.linalg.lstsq(X3, H_arr, rcond=None)
        R2_3 = 1 - np.sum((H_arr - X3@c3_coeff)**2) / np.sum((H_arr - H_arr.mean())**2)
        print(f"  R²(c₃ + α₂) = {R2_3:.6f} (ΔR² = {R2_3 - R2_1:.6f})")
        print(f"    H ≈ {c3_coeff[0]:.4f}·c₃ + {c3_coeff[1]:.4f}·α₂ + {c3_coeff[2]:.4f}")

    # Full: c₃ + c₅ + α₂
    if np.std(c5_arr) > 0 and np.std(a2_arr) > 0:
        X4 = np.column_stack([c3_arr, c5_arr, a2_arr, np.ones(N)])
        c4, _, _, _ = np.linalg.lstsq(X4, H_arr, rcond=None)
        R2_4 = 1 - np.sum((H_arr - X4@c4)**2) / np.sum((H_arr - H_arr.mean())**2)
        print(f"  R²(c₃ + c₅ + α₂) = {R2_4:.6f} (ΔR² = {R2_4 - R2_1:.6f})")

    # OCF exact: H = 1 + 2(c₃+c₅+c₇+...) + 4α₂ + 8α₃ + ...
    # So H = 1 + 2·α₁ + 4·α₂ + ...
    alpha1_arr = c3_arr + c5_arr  # Only 3- and 5-cycles (ignoring 7-cycles at this n)
    X5 = np.column_stack([alpha1_arr, a2_arr, np.ones(N)])
    if np.std(alpha1_arr) > 0 and np.std(a2_arr) > 0:
        c5_coeff, _, _, _ = np.linalg.lstsq(X5, H_arr, rcond=None)
        R2_5 = 1 - np.sum((H_arr - X5@c5_coeff)**2) / np.sum((H_arr - H_arr.mean())**2)
        print(f"  R²(α₁ + α₂) = {R2_5:.6f}")
        print(f"    H ≈ {c5_coeff[0]:.4f}·α₁ + {c5_coeff[1]:.4f}·α₂ + {c5_coeff[2]:.4f}")
        print(f"    (OCF predicts: H = 2·α₁ + 4·α₂ + 1)")

    # What fraction of α₁ is c₃?
    if np.std(alpha1_arr) > 0:
        frac = np.mean(c3_arr) / np.mean(alpha1_arr)
        corr_c3_a1 = np.corrcoef(c3_arr, alpha1_arr)[0,1]
        print(f"\n  c₃ fraction of α₁: {frac:.4f}")
        print(f"  corr(c₃, α₁) = {corr_c3_a1:.6f}")
        print(f"  → c₃ is {frac*100:.1f}% of total odd cycles, correlated at {corr_c3_a1:.3f}")

    return data


def phase_transition_n6():
    """Compute β_c at n=6 (exhaustive)."""
    n = 6
    print(f"\n{'='*70}")
    print(f"  PHASE TRANSITION at n = {n} (exhaustive)")
    print(f"{'='*70}")

    H_list = []
    for A in all_tournaments(n):
        H_list.append(compute_H(A, n))

    H_arr = np.array(H_list, dtype=float)
    N = len(H_arr)
    print(f"  {N} tournaments, H in [{H_arr.min():.0f}, {H_arr.max():.0f}]")

    beta_range = np.linspace(0, 2, 201)
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


# ─── MAIN ───

if __name__ == "__main__":
    print("╔══════════════════════════════════════════════════════════════════╗")
    print("║  PROVING THE OPEN QUESTIONS                                    ║")
    print("║  opus-2026-04-04-S23                                           ║")
    print("╚══════════════════════════════════════════════════════════════════╝")

    t0 = time.time()

    # Task 1: Claim A decomposition at n=4→5
    claim_a_decomposition(5)

    t1 = time.time()
    print(f"\n  [Claim A decomp n=5: {t1-t0:.1f}s]")

    # Task 2: Variance decomposition at n=5 and n=6
    variance_decomposition(5)
    t2 = time.time()
    print(f"\n  [Variance decomp n=5: {t2-t1:.1f}s]")

    variance_decomposition(6)
    t3 = time.time()
    print(f"\n  [Variance decomp n=6: {t3-t2:.1f}s]")

    # Task 3: Phase transition at n=6
    beta_c_6 = phase_transition_n6()
    t4 = time.time()
    print(f"\n  [Phase transition n=6: {t4-t3:.1f}s]")

    # Summary
    print(f"\n{'='*70}")
    print(f"  PHASE TRANSITION SCALING:")
    print(f"  β_c(5) ≈ 0.7 (from S16)")
    print(f"  β_c(6) ≈ {beta_c_6:.4f}")
    print(f"{'='*70}")

    print(f"\n  TOTAL: {t4-t0:.1f}s")
