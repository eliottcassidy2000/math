#!/usr/bin/env python3
"""
afm_fiber_bundle_s17.py — opus-2026-04-04-S17

THE FIBER BUNDLE ANTIFERROMAGNET

Integration of the AFM framework (S16) with the fiber bundle recursion (S12).

KEY IDEA: T_n = (T_{n-1}, σ) where σ ∈ {0,1}^{n-1} is the extension signature.
Every AFM quantity decomposes:
  - c₃(T_n) = c₃(T_{n-1}) + Δc₃(T_{n-1}, σ)
  - H(T_n) = H(T_{n-1}) + ΔH(T_{n-1}, σ)   [Claim A: ΔH = 2Σ μ(C)]
  - score_var(T_n) = f(scores(T_{n-1}), σ)

The extension σ is the "fiber coordinate" — it determines how the new vertex
couples to the existing spin system. In AFM terms:
  - σ = the coupling of a new spin to the existing frustrated lattice
  - Δc₃ = the frustration the new spin generates with existing pairs
  - ΔH = the partition function contribution from new cycles through vertex n

MAGNON DECOMPOSITION:
  - Inner magnon: flip arc (i,j) with i,j < n (stays in same fiber)
  - Boundary magnon: flip arc (n,i) for some i (changes σ, crosses fibers)

This decomposition should reveal the magnon anisotropy that was invisible
in the S_n-averaged ensemble.
"""

import numpy as np
from itertools import combinations, permutations
from math import comb, factorial
from collections import defaultdict, Counter
import time
import sys

# ─────────────────────────────────────────────
#  Core utilities
# ─────────────────────────────────────────────

def compute_H(A, n):
    """H(T) via Held-Karp DP."""
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

def scores_raw(A, n):
    """Raw (unsorted) out-degrees."""
    return [sum(A[i][j] for j in range(n) if j != i) for i in range(n)]

def scores(A, n):
    return tuple(sorted(scores_raw(A, n)))

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
    s = scores_raw(A, n)
    mean = sum(s) / n
    return sum((x - mean)**2 for x in s) / n

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

# ─────────────────────────────────────────────
#  Fiber bundle structure
# ─────────────────────────────────────────────

def extend_tournament(A_sub, n_sub, sigma):
    """Extend T_{n-1} by adding vertex n-1 (0-indexed as n_sub) with signature sigma.
    sigma[i] = 1 means new vertex beats vertex i."""
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
    """Number of new 3-cycles created by adding vertex n with signature sigma.

    A new 3-cycle involves vertex n and two existing vertices i,j.
    The triple {n, i, j} forms a 3-cycle iff the arcs form a directed cycle:
      n→i→j→n or n→j→i→n

    For n→i: sigma[i]=1. For i→n: sigma[i]=0.
    For i→j: A_sub[i][j]=1. For j→i: A_sub[j][i]=1.

    The triple (n,i,j) is a 3-cycle iff:
      (sigma[i]=1 and A_sub[i][j]=1 and sigma[j]=0) — cycle n→i→j→n
      OR
      (sigma[j]=1 and A_sub[j][i]=1 and sigma[i]=0) — cycle n→j→i→n
    """
    dc3 = 0
    for i in range(n_sub):
        for j in range(i+1, n_sub):
            # Check n→i→j→n: sigma[i]=1, A[i][j]=1, sigma[j]=0
            if sigma[i] == 1 and A_sub[i][j] == 1 and sigma[j] == 0:
                dc3 += 1
            # Check n→j→i→n: sigma[j]=1, A[j][i]=1, sigma[i]=0
            elif sigma[j] == 1 and A_sub[j][i] == 1 and sigma[i] == 0:
                dc3 += 1
    return dc3

def delta_c3_formula(A_sub, n_sub, sigma):
    """Express Δc₃ in terms of score structure.

    Vertices split into W = {i : sigma[i]=1} (beaten by n) and
    L = {i : sigma[i]=0} (beat n). |W| = sum(sigma) = score(n).

    A new 3-cycle n→i→j→n requires i∈W, j∈L, i→j.
    So Δc₃ = #{(i,j) : i∈W, j∈L, i→j} = number of arcs from W to L.

    In the tournament T_{n-1}: total arcs from W to L = Σ_{i∈W} (# vertices in L that i beats).
    For each i∈W: # beaten vertices in L = out-degree of i restricted to L.
    """
    W = [i for i in range(n_sub) if sigma[i] == 1]
    L = [i for i in range(n_sub) if sigma[i] == 0]
    L_set = set(L)

    arcs_W_to_L = 0
    for i in W:
        for j in L:
            if A_sub[i][j]:
                arcs_W_to_L += 1

    return arcs_W_to_L


# ─────────────────────────────────────────────
#  1. FRUSTRATION PROPAGATION
# ─────────────────────────────────────────────

def frustration_propagation(n_target):
    """Verify: c₃(T_n) = c₃(T_{n-1}) + Δc₃(T_{n-1}, σ) for all tournaments."""
    n_sub = n_target - 1
    print(f"\n{'='*70}")
    print(f"  FRUSTRATION PROPAGATION: n={n_sub} → n={n_target}")
    print(f"{'='*70}")

    errors = 0
    total = 0

    # Statistics on Δc₃
    dc3_by_score_n = defaultdict(list)
    dc3_vs_parent_c3 = []

    for A_sub in all_tournaments(n_sub):
        c3_sub = count_3cycles(A_sub, n_sub)
        s_sub = scores_raw(A_sub, n_sub)

        for sigma_int in range(1 << n_sub):
            sigma = tuple((sigma_int >> i) & 1 for i in range(n_sub))
            score_n = sum(sigma)

            A_n = extend_tournament(A_sub, n_sub, sigma)
            c3_n = count_3cycles(A_n, n_target)

            dc3 = delta_c3(A_sub, n_sub, sigma)
            dc3_formula = delta_c3_formula(A_sub, n_sub, sigma)

            # Verify
            if c3_n != c3_sub + dc3:
                errors += 1
                if errors <= 3:
                    print(f"  ERROR: c3_sub={c3_sub}, dc3={dc3}, c3_n={c3_n}")

            if dc3 != dc3_formula:
                print(f"  FORMULA MISMATCH: dc3={dc3}, formula={dc3_formula}")

            dc3_by_score_n[score_n].append(dc3)
            dc3_vs_parent_c3.append((c3_sub, dc3, score_n))
            total += 1

    print(f"\n  Total extensions checked: {total}")
    print(f"  Errors in c₃ propagation: {errors}")
    nn = n_target
    print(f"  c3(T_{nn}) = c3(T_{nn-1}) + dc3 {'VERIFIED' if errors == 0 else 'FAILED'}")

    # Δc₃ statistics by score of new vertex
    print(f"\n  Δc₃ by score of new vertex (score(n) = Σσ):")
    print(f"  {'score(n)':>8}  {'mean Δc₃':>10}  {'std Δc₃':>10}  {'min':>5}  {'max':>5}  {'count':>7}")
    for s in sorted(dc3_by_score_n.keys()):
        vals = dc3_by_score_n[s]
        print(f"  {s:8d}  {np.mean(vals):10.4f}  {np.std(vals):10.4f}  "
              f"{min(vals):5d}  {max(vals):5d}  {len(vals):7d}")

    # Key formula test: E[Δc₃ | score(n)=s] = s(n-1-s)/?
    # For random tournament: arcs from W to L ~ Binomial(s*(n-1-s), 1/2)
    # So E[Δc₃] = s*(n-1-s)/2
    print(f"\n  Test: E[Δc₃ | score(n)=s] vs s*(n-1-s)/2 (random baseline):")
    for s in sorted(dc3_by_score_n.keys()):
        vals = dc3_by_score_n[s]
        expected = s * (n_sub - s) / 2
        actual = np.mean(vals)
        print(f"    score(n)={s}: actual={actual:.4f}, random={expected:.4f}, "
              f"ratio={actual/expected:.4f}" if expected > 0 else
              f"    score(n)={s}: actual={actual:.4f}, random=0")

    return dc3_by_score_n


# ─────────────────────────────────────────────
#  2. H PROPAGATION THROUGH FIBERS
# ─────────────────────────────────────────────

def h_propagation(n_target):
    """Compute ΔH = H(T_n) - H(T_{n-1}) for each extension.
    By Claim A: ΔH = 2Σ_{C∋n} μ(C).
    This is the 'extension energy' — the AFM energy cost of adding vertex n."""
    n_sub = n_target - 1
    print(f"\n{'='*70}")
    print(f"  H PROPAGATION (CLAIM A): n={n_sub} → n={n_target}")
    print(f"{'='*70}")

    dH_by_score_n = defaultdict(list)
    dH_by_dc3 = defaultdict(list)
    dH_by_parent_H = defaultdict(list)
    all_data = []

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

            dH_by_score_n[score_n].append(dH)
            dH_by_dc3[dc3].append(dH)
            dH_by_parent_H[H_sub].append(dH)

            all_data.append({
                'H_sub': H_sub, 'H_n': H_n, 'dH': dH,
                'c3_sub': c3_sub, 'dc3': dc3,
                'score_n': score_n, 'sigma': sigma
            })

    print(f"\n  ΔH by score of new vertex:")
    print(f"  {'score(n)':>8}  {'mean ΔH':>10}  {'std ΔH':>10}  {'min':>6}  {'max':>6}")
    for s in sorted(dH_by_score_n.keys()):
        vals = dH_by_score_n[s]
        print(f"  {s:8d}  {np.mean(vals):10.4f}  {np.std(vals):10.4f}  "
              f"{min(vals):6d}  {max(vals):6d}")

    print(f"\n  ΔH by Δc₃ (new frustration):")
    print(f"  {'Δc₃':>5}  {'mean ΔH':>10}  {'std ΔH':>10}  {'min':>6}  {'max':>6}  {'count':>7}")
    for dc3 in sorted(dH_by_dc3.keys()):
        vals = dH_by_dc3[dc3]
        print(f"  {dc3:5d}  {np.mean(vals):10.4f}  {np.std(vals):10.4f}  "
              f"{min(vals):6d}  {max(vals):6d}  {len(vals):7d}")

    # KEY TEST: is ΔH linearly related to Δc₃?
    all_dH = np.array([d['dH'] for d in all_data])
    all_dc3 = np.array([d['dc3'] for d in all_data])
    all_score_n = np.array([d['score_n'] for d in all_data])
    all_H_sub = np.array([d['H_sub'] for d in all_data])
    all_c3_sub = np.array([d['c3_sub'] for d in all_data])

    corr = np.corrcoef(all_dc3, all_dH)[0,1]
    print(f"\n  Correlation(Δc₃, ΔH) = {corr:.6f}")

    # Linear model: ΔH ~ a*Δc₃ + b*score_n + c*H_sub + d
    X = np.column_stack([all_dc3, all_score_n, all_H_sub, np.ones(len(all_data))])
    coeffs, _, _, _ = np.linalg.lstsq(X, all_dH, rcond=None)
    pred = X @ coeffs
    R2 = 1 - np.sum((all_dH - pred)**2) / np.sum((all_dH - all_dH.mean())**2)
    print(f"\n  Linear model: ΔH ≈ {coeffs[0]:.4f}*Δc₃ + {coeffs[1]:.4f}*score(n) "
          f"+ {coeffs[2]:.4f}*H_sub + {coeffs[3]:.4f}")
    print(f"  R² = {R2:.6f}")

    # Simpler model: ΔH ~ a*Δc₃ + b
    X2 = np.column_stack([all_dc3, np.ones(len(all_data))])
    coeffs2, _, _, _ = np.linalg.lstsq(X2, all_dH, rcond=None)
    pred2 = X2 @ coeffs2
    R2_simple = 1 - np.sum((all_dH - pred2)**2) / np.sum((all_dH - all_dH.mean())**2)
    print(f"  Simple model: ΔH ≈ {coeffs2[0]:.4f}*Δc₃ + {coeffs2[1]:.4f}")
    print(f"  R²(simple) = {R2_simple:.6f}")

    # ΔH by parent H
    print(f"\n  ΔH by parent H(T_{n_sub}):")
    print(f"  {'H_sub':>6}  {'mean ΔH':>10}  {'std ΔH':>10}")
    for H_sub in sorted(dH_by_parent_H.keys()):
        vals = dH_by_parent_H[H_sub]
        print(f"  {H_sub:6d}  {np.mean(vals):10.4f}  {np.std(vals):10.4f}")

    # The FIBER CONNECTION: for a fixed parent T_{n-1}, how does ΔH depend on σ?
    print(f"\n  FIBER CONNECTION: ΔH variation within single parent fiber")
    # Pick one parent at each H level
    parent_examples = {}
    for A_sub in all_tournaments(n_sub):
        H_sub = compute_H(A_sub, n_sub)
        if H_sub not in parent_examples:
            parent_examples[H_sub] = A_sub
        if len(parent_examples) >= 5:
            break

    for H_sub, A_sub in sorted(parent_examples.items()):
        dH_list = []
        dc3_list = []
        for sigma_int in range(1 << n_sub):
            sigma = tuple((sigma_int >> i) & 1 for i in range(n_sub))
            A_n = extend_tournament(A_sub, n_sub, sigma)
            H_n = compute_H(A_n, n_target)
            dc3 = delta_c3(A_sub, n_sub, sigma)
            dH_list.append(H_n - H_sub)
            dc3_list.append(dc3)
        dH_arr = np.array(dH_list)
        dc3_arr = np.array(dc3_list)
        corr = np.corrcoef(dc3_arr, dH_arr)[0,1] if np.std(dc3_arr) > 0 else float('nan')
        print(f"    Parent H={H_sub}: ΔH range [{dH_arr.min()}, {dH_arr.max()}], "
              f"mean={dH_arr.mean():.2f}, corr(Δc₃,ΔH)={corr:.4f}")

    return all_data


# ─────────────────────────────────────────────
#  3. INNER vs BOUNDARY MAGNON DECOMPOSITION
# ─────────────────────────────────────────────

def magnon_decomposition(n_target):
    """Decompose arc flips into inner magnons and boundary magnons.

    Inner magnon: flip arc (i,j) with i,j < n (both in T_{n-1}).
      - Changes T_{n-1} but preserves σ.
      - The ΔH comes from changes in the inherited structure.

    Boundary magnon: flip arc (n,i) for some i.
      - Preserves T_{n-1} but flips σ[i].
      - The ΔH comes from changing how vertex n couples to the lattice.
    """
    n_sub = n_target - 1
    print(f"\n{'='*70}")
    print(f"  INNER vs BOUNDARY MAGNON DECOMPOSITION — n = {n_target}")
    print(f"{'='*70}")

    inner_dH = []  # |ΔH| for inner magnons
    boundary_dH = []  # |ΔH| for boundary magnons
    inner_dH_signed = []
    boundary_dH_signed = []

    # Also track by position
    boundary_dH_by_vertex = defaultdict(list)  # vertex i → list of |ΔH|

    count = 0
    for A in all_tournaments(n_target):
        H0 = compute_H(A, n_target)

        # Inner magnons: flip (i,j) for i,j < n_sub (but we need i < n_target, j < n_target, both != n_sub)
        # Actually: vertex n is vertex n_target-1 (0-indexed).
        new_v = n_target - 1

        for i in range(n_target):
            for j in range(i+1, n_target):
                # Flip arc between i and j
                A2 = [row[:] for row in A]
                A2[i][j] = 1 - A2[i][j]
                A2[j][i] = 1 - A2[j][i]
                H1 = compute_H(A2, n_target)
                dH = H1 - H0

                if i == new_v or j == new_v:
                    # Boundary magnon
                    boundary_dH.append(abs(dH))
                    boundary_dH_signed.append(dH)
                    other_v = j if i == new_v else i
                    boundary_dH_by_vertex[other_v].append(abs(dH))
                else:
                    # Inner magnon
                    inner_dH.append(abs(dH))
                    inner_dH_signed.append(dH)

        count += 1
        if count % 2000 == 0:
            print(f"  ... {count} tournaments processed", file=sys.stderr)

    inner_arr = np.array(inner_dH)
    boundary_arr = np.array(boundary_dH)

    print(f"\n  INNER magnons (flip within T_{{n-1}}):")
    print(f"    Count: {len(inner_arr)}")
    print(f"    Mean |ΔH|: {inner_arr.mean():.4f}")
    print(f"    Std |ΔH|: {inner_arr.std():.4f}")
    print(f"    ΔH=0 rate: {100*np.sum(inner_arr==0)/len(inner_arr):.2f}%")

    print(f"\n  BOUNDARY magnons (flip arc to/from vertex n):")
    print(f"    Count: {len(boundary_arr)}")
    print(f"    Mean |ΔH|: {boundary_arr.mean():.4f}")
    print(f"    Std |ΔH|: {boundary_arr.std():.4f}")
    print(f"    ΔH=0 rate: {100*np.sum(boundary_arr==0)/len(boundary_arr):.2f}%")

    print(f"\n  RATIO: boundary/inner mean|ΔH| = {boundary_arr.mean()/inner_arr.mean():.4f}")

    # Boundary magnon by vertex
    print(f"\n  Boundary magnon spectrum by vertex:")
    for v in sorted(boundary_dH_by_vertex.keys()):
        vals = np.array(boundary_dH_by_vertex[v])
        print(f"    vertex {v}: mean|ΔH|={vals.mean():.4f}, ΔH=0 rate={100*np.sum(vals==0)/len(vals):.2f}%")

    # ΔH distributions
    inner_counter = Counter(inner_dH_signed)
    boundary_counter = Counter(boundary_dH_signed)
    print(f"\n  Inner ΔH distribution (top 10):")
    for dH, cnt in sorted(inner_counter.items(), key=lambda x: -x[1])[:10]:
        print(f"    ΔH={dH:+4d}: {cnt} ({100*cnt/len(inner_dH_signed):.2f}%)")
    print(f"\n  Boundary ΔH distribution (top 10):")
    for dH, cnt in sorted(boundary_counter.items(), key=lambda x: -x[1])[:10]:
        print(f"    ΔH={dH:+4d}: {cnt} ({100*cnt/len(boundary_dH_signed):.2f}%)")

    return inner_dH, boundary_dH


# ─────────────────────────────────────────────
#  4. FIBER-CONDITIONED MAGNON SPECTRUM (per iso class)
# ─────────────────────────────────────────────

def per_class_magnon(n_target):
    """Compute magnon spectrum conditioned on isomorphism class.
    This should show the ANISOTROPY hidden by the S_n average."""
    n_sub = n_target - 1
    print(f"\n{'='*70}")
    print(f"  PER-CLASS MAGNON SPECTRUM — n = {n_target}")
    print(f"{'='*70}")

    # Group tournaments by canonical form
    canon_map = {}  # tuple(upper triangle) → canonical key
    class_data = defaultdict(lambda: {'tournaments': [], 'H': None})

    for A in all_tournaments(n_target):
        key = tuple(A[i][j] for i in range(n_target) for j in range(i+1, n_target))
        # Simple canonical: sorted scores + H as proxy
        H = compute_H(A, n_target)
        c3 = count_3cycles(A, n_target)
        sc = scores(A, n_target)
        class_key = (sc, H, c3)  # Coarse grouping (not exact iso, but good proxy)
        class_data[class_key]['tournaments'].append(A)
        class_data[class_key]['H'] = H

    print(f"  {len(class_data)} coarse classes found")

    # For each class, compute magnon spectrum (inner vs boundary)
    print(f"\n  {'H':>4}  {'c₃':>4}  {'#tours':>6}  {'inner_mean':>11}  {'bound_mean':>11}  "
          f"{'ratio':>7}  {'inner_0%':>8}  {'bound_0%':>8}")

    class_results = []
    new_v = n_target - 1

    for class_key in sorted(class_data.keys(), key=lambda x: x[1]):
        tours = class_data[class_key]['tournaments']
        sc, H, c3 = class_key
        if len(tours) > 50:
            tours = tours[:50]  # Sample for speed

        inner_vals = []
        boundary_vals = []

        for A in tours:
            for i in range(n_target):
                for j in range(i+1, n_target):
                    A2 = [row[:] for row in A]
                    A2[i][j] = 1 - A2[i][j]
                    A2[j][i] = 1 - A2[j][i]
                    H1 = compute_H(A2, n_target)
                    dH = abs(H1 - H)

                    if i == new_v or j == new_v:
                        boundary_vals.append(dH)
                    else:
                        inner_vals.append(dH)

        if inner_vals and boundary_vals:
            inner_mean = np.mean(inner_vals)
            bound_mean = np.mean(boundary_vals)
            inner_zero = 100 * sum(1 for v in inner_vals if v == 0) / len(inner_vals)
            bound_zero = 100 * sum(1 for v in boundary_vals if v == 0) / len(boundary_vals)
            ratio = bound_mean / inner_mean if inner_mean > 0 else float('inf')

            class_results.append({
                'H': H, 'c3': c3, 'scores': sc,
                'inner_mean': inner_mean, 'bound_mean': bound_mean,
                'ratio': ratio, 'inner_zero': inner_zero, 'bound_zero': bound_zero
            })

            if len(class_results) <= 30:  # Print first 30
                print(f"  {H:4d}  {c3:4d}  {len(tours):6d}  {inner_mean:11.4f}  "
                      f"{bound_mean:11.4f}  {ratio:7.3f}  {inner_zero:7.2f}%  {bound_zero:7.2f}%")

    # Summary statistics
    if class_results:
        ratios = [r['ratio'] for r in class_results if r['ratio'] < 100]
        print(f"\n  Summary: boundary/inner ratio across classes:")
        print(f"    Mean: {np.mean(ratios):.4f}")
        print(f"    Std: {np.std(ratios):.4f}")
        print(f"    Range: [{min(ratios):.4f}, {max(ratios):.4f}]")
        print(f"    INTERPRETATION: ratio > 1 → boundary magnons cost MORE energy")
        print(f"                   ratio < 1 → inner magnons cost MORE energy")

    return class_results


# ─────────────────────────────────────────────
#  5. FIBER BUNDLE PARTITION FUNCTION
# ─────────────────────────────────────────────

def fiber_partition_function(n_target):
    """The partition function Z_n(β) = Σ exp(β·H(T_n)) decomposes over fibers:

    Z_n(β) = Σ_{T_{n-1}} Σ_σ exp(β·H(T_{n-1}, σ))
           = Σ_{T_{n-1}} exp(β·H(T_{n-1})) · Σ_σ exp(β·ΔH(T_{n-1}, σ))
           = Σ_{T_{n-1}} exp(β·H(T_{n-1})) · z(T_{n-1}, β)

    where z(T_{n-1}, β) = Σ_σ exp(β·ΔH(T_{n-1}, σ)) is the FIBER partition function.
    This measures how much the extension possibilities amplify the parent's weight.
    """
    n_sub = n_target - 1
    print(f"\n{'='*70}")
    print(f"  FIBER PARTITION FUNCTION — n={n_sub} → n={n_target}")
    print(f"{'='*70}")

    # For each parent tournament, compute the fiber PF
    fiber_pf_data = []

    for A_sub in all_tournaments(n_sub):
        H_sub = compute_H(A_sub, n_sub)
        c3_sub = count_3cycles(A_sub, n_sub)

        dH_list = []
        for sigma_int in range(1 << n_sub):
            sigma = tuple((sigma_int >> i) & 1 for i in range(n_sub))
            A_n = extend_tournament(A_sub, n_sub, sigma)
            H_n = compute_H(A_n, n_target)
            dH_list.append(H_n - H_sub)

        dH_arr = np.array(dH_list)

        fiber_pf_data.append({
            'H_sub': H_sub, 'c3_sub': c3_sub,
            'mean_dH': dH_arr.mean(),
            'max_dH': dH_arr.max(),
            'min_dH': dH_arr.min(),
            'std_dH': dH_arr.std(),
            'dH_values': dH_arr,
        })

    # Compute fiber PF at various β
    print(f"\n  Fiber partition function z(T, β) = Σ_σ exp(β·ΔH):")
    print(f"  {'H_sub':>6}  {'c₃':>4}  {'mean ΔH':>9}  {'z(β=0)':>8}  {'z(β=0.5)':>10}  "
          f"{'z(β=1)':>10}  {'z(β=2)':>10}")

    for d in sorted(fiber_pf_data, key=lambda x: x['H_sub']):
        dH = d['dH_values']
        z0 = len(dH)  # = 2^{n-1}
        z05 = np.sum(np.exp(0.5 * dH))
        z1 = np.sum(np.exp(1.0 * dH))
        z2 = np.sum(np.exp(2.0 * dH))
        d['z'] = {0: z0, 0.5: z05, 1.0: z1, 2.0: z2}

        if d['H_sub'] in [1, 3, 5] or d == fiber_pf_data[0] or d == fiber_pf_data[-1]:
            print(f"  {d['H_sub']:6d}  {d['c3_sub']:4d}  {d['mean_dH']:9.4f}  "
                  f"{z0:8.0f}  {z05:10.2f}  {z1:10.2f}  {z2:10.2f}")

    # The key question: does the fiber PF correlate with parent H?
    H_sub_arr = np.array([d['H_sub'] for d in fiber_pf_data])
    c3_sub_arr = np.array([d['c3_sub'] for d in fiber_pf_data])
    mean_dH_arr = np.array([d['mean_dH'] for d in fiber_pf_data])

    print(f"\n  Correlations:")
    print(f"    corr(H_sub, mean_ΔH)  = {np.corrcoef(H_sub_arr, mean_dH_arr)[0,1]:.6f}")
    print(f"    corr(c₃_sub, mean_ΔH) = {np.corrcoef(c3_sub_arr, mean_dH_arr)[0,1]:.6f}")

    # Compute β_c through recursion
    for beta in [0.5, 1.0, 2.0]:
        z_vals = np.array([d['z'][beta] for d in fiber_pf_data])
        mean_z = np.mean(z_vals)
        std_z = np.std(z_vals)
        # The "fiber amplification" at β
        print(f"\n    β={beta}: mean z = {mean_z:.2f}, std z = {std_z:.2f}, "
              f"CV = {std_z/mean_z:.4f}")
        # More frustrated parents amplify MORE at positive β
        corr_H_z = np.corrcoef(H_sub_arr, z_vals)[0,1]
        print(f"      corr(H_sub, z) = {corr_H_z:.6f}")
        print(f"      → {'Frustrated parents amplify MORE' if corr_H_z > 0 else 'Ordered parents amplify MORE'}")

    return fiber_pf_data


# ─────────────────────────────────────────────
#  MAIN
# ─────────────────────────────────────────────

if __name__ == "__main__":
    print("╔══════════════════════════════════════════════════════════════════╗")
    print("║  THE FIBER BUNDLE ANTIFERROMAGNET                              ║")
    print("║  opus-2026-04-04-S17                                           ║")
    print("╚══════════════════════════════════════════════════════════════════╝")

    t0 = time.time()

    # 1. Frustration propagation
    frustration_propagation(5)  # n=4→5
    frustration_propagation(6)  # n=5→6

    t1 = time.time()
    print(f"\n  [Frustration propagation: {t1-t0:.1f}s]")

    # 2. H propagation
    h_propagation(5)  # n=4→5

    t2 = time.time()
    print(f"\n  [H propagation: {t2-t1:.1f}s]")

    # 3. Magnon decomposition
    magnon_decomposition(5)  # n=5

    t3 = time.time()
    print(f"\n  [Magnon decomposition: {t3-t2:.1f}s]")

    # 4. Per-class magnon spectrum
    per_class_magnon(5)

    t4 = time.time()
    print(f"\n  [Per-class magnon: {t4-t3:.1f}s]")

    # 5. Fiber partition function
    fiber_partition_function(5)

    t5 = time.time()
    print(f"\n  [Fiber PF: {t5-t4:.1f}s]")

    print(f"\n\n{'='*70}")
    print(f"  TOTAL: {t5-t0:.1f}s")
    print(f"{'='*70}")
