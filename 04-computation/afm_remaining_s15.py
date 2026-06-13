#!/usr/bin/env python3
"""
afm_remaining_s15.py — opus-2026-04-04-S15

Run the remaining parts that didn't finish in the main script.
"""

import numpy as np
from itertools import combinations, permutations
from math import comb
from collections import defaultdict, Counter
import time

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
        yield mask, A

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

def score_variance(A, n):
    s = [sum(A[i]) for i in range(n)]
    mean = sum(s) / n
    return sum((x - mean)**2 for x in s) / n


def find_odd_cycles(A, n):
    """Find all directed odd cycles in tournament A (as frozensets of vertices)."""
    cycles = set()
    for length in range(3, n+1, 2):
        for combo in combinations(range(n), length):
            verts = list(combo)
            for perm in permutations(verts):
                if perm[0] == min(perm):
                    if all(A[perm[i]][perm[(i+1) % length]] for i in range(length)):
                        cycles.add(frozenset(perm))
    return list(cycles)


def count_independent_sets(cycles, n_cycles):
    """Count independent sets by size. Cycles share vertex = adjacent in conflict graph."""
    adj = defaultdict(set)
    for i in range(n_cycles):
        for j in range(i+1, n_cycles):
            if cycles[i] & cycles[j]:
                adj[i].add(j)
                adj[j].add(i)

    alphas = defaultdict(int)
    for mask in range(1 << n_cycles):
        subset = [i for i in range(n_cycles) if mask & (1 << i)]
        is_indep = True
        for i in range(len(subset)):
            for j in range(i+1, len(subset)):
                if subset[j] in adj.get(subset[i], set()):
                    is_indep = False
                    break
            if not is_indep:
                break
        if is_indep:
            alphas[len(subset)] += 1
    return dict(alphas)


# ─── 1. AFM GROUND STATE for n=6 ───

def ground_state_n6():
    n = 6
    print(f"\n{'='*80}")
    print(f"  AFM GROUND STATE — n = {n}")
    print(f"{'='*80}")

    max_c3 = 0
    max_H = 0

    data_by_H = defaultdict(list)
    data_by_c3 = defaultdict(list)

    for mask, A in all_tournaments(n):
        H = compute_H(A, n)
        c3 = count_3cycles(A, n)
        sc = scores(A, n)
        sv = score_variance(A, n)
        max_c3 = max(max_c3, c3)
        max_H = max(max_H, H)
        data_by_H[H].append({'c3': c3, 'sv': sv, 'scores': sc})
        data_by_c3[c3].append({'H': H, 'sv': sv, 'scores': sc})

    print(f"  Max c₃ = {max_c3}, Max H = {max_H}")

    H_max_data = data_by_H[max_H]
    c3_vals = [d['c3'] for d in H_max_data]
    print(f"\n  H-MAXIMIZERS (H={max_H}): {len(H_max_data)} tournaments")
    print(f"    c₃ values: {sorted(set(c3_vals))}")
    print(f"    Scores: {sorted(set(d['scores'] for d in H_max_data))}")

    c3_max_data = data_by_c3[max_c3]
    H_vals = [d['H'] for d in c3_max_data]
    print(f"\n  FRUSTRATION-MAX (c₃={max_c3}): {len(c3_max_data)} tournaments")
    print(f"    H values: {sorted(set(H_vals))}")
    print(f"    Scores: {sorted(set(d['scores'] for d in c3_max_data))}")

    # Does max c₃ always give max H?
    if max(H_vals) == max_H:
        print(f"\n  ✓ Frustration-maximizer ACHIEVES max H = {max_H}")
    else:
        print(f"\n  ✗ Frustration-maximizer only reaches H = {max(H_vals)} < {max_H}")


# ─── 2. THERMODYNAMIC FUNCTIONS at n=5 ───

def thermodynamics_n5():
    n = 5
    print(f"\n{'='*80}")
    print(f"  THERMODYNAMIC FUNCTIONS — n = {n}")
    print(f"{'='*80}")

    H_list, c3_list, sv_list = [], [], []
    for mask, A in all_tournaments(n):
        H_list.append(compute_H(A, n))
        c3_list.append(count_3cycles(A, n))
        sv_list.append(score_variance(A, n))

    H_arr = np.array(H_list, dtype=np.float64)
    c3_arr = np.array(c3_list, dtype=np.float64)
    sv_arr = np.array(sv_list, dtype=np.float64)
    N = len(H_arr)

    print(f"  {N} tournaments, H in [{H_arr.min():.0f}, {H_arr.max():.0f}]")

    beta_range = np.linspace(-3, 5, 81)

    max_C = 0
    beta_max_C = 0

    print(f"\n  {'β':>8}  {'<H>':>10}  {'Var(H)':>10}  {'C(β)':>10}  "
          f"{'<c₃>':>8}  {'<sv>':>8}")

    results = []
    for beta in beta_range:
        log_w = beta * H_arr
        log_w -= np.max(log_w)
        w = np.exp(log_w)
        Z = w.sum()
        w /= Z

        mean_H = np.sum(w * H_arr)
        var_H = np.sum(w * (H_arr - mean_H)**2)
        C = beta**2 * var_H
        mean_c3 = np.sum(w * c3_arr)
        mean_sv = np.sum(w * sv_arr)

        results.append((beta, mean_H, var_H, C, mean_c3, mean_sv))

        if C > max_C:
            max_C = C
            beta_max_C = beta

    # Print selected values
    for beta, mean_H, var_H, C, mean_c3, mean_sv in results:
        if abs(beta) < 0.05 or abs(beta - 1) < 0.1 or abs(beta - 2) < 0.1 or \
           abs(beta - beta_max_C) < 0.1 or abs(beta - beta_range[0]) < 0.05 or \
           abs(beta - beta_range[-1]) < 0.05 or abs(beta - 0.5) < 0.05:
            print(f"  {beta:8.3f}  {mean_H:10.4f}  {var_H:10.4f}  {C:10.4f}  "
                  f"{mean_c3:8.3f}  {mean_sv:8.4f}")

    print(f"\n  SPECIFIC HEAT PEAK: C_max = {max_C:.4f} at β = {beta_max_C:.3f}")
    print(f"  This is the 'PHASE TRANSITION' of the tournament ensemble")

    # Entropy: S(β) = log Z + β <H>
    print(f"\n  Phase diagram summary:")
    for beta in [-2, -1, 0, 0.5, 1, 2, 5]:
        log_w = beta * H_arr
        log_w -= np.max(log_w)
        w = np.exp(log_w)
        Z = w.sum()
        w /= Z
        mean_H = np.sum(w * H_arr)
        mean_c3 = np.sum(w * c3_arr)
        mean_sv = np.sum(w * sv_arr)
        S = np.log(Z) + beta * mean_H  # normalized entropy
        print(f"    β={beta:+5.1f}: <H>={mean_H:7.2f}, <c₃>={mean_c3:5.2f}, "
              f"<sv>={mean_sv:5.3f}, S={S:.3f}")
    print(f"    β → -∞: FERROMAGNETIC (transitive, H=1, c₃=0)")
    print(f"    β → +∞: ANTIFERROMAGNETIC (regular, max H, max c₃)")


# ─── 3. FORBIDDEN VALUE DEEP ANALYSIS at n=5 ───

def forbidden_deep_n5():
    n = 5
    print(f"\n{'='*80}")
    print(f"  FORBIDDEN H=7: OCF DECOMPOSITION at n = {n}")
    print(f"{'='*80}")

    ocf_data = []
    count = 0
    for mask, A in all_tournaments(n):
        H = compute_H(A, n)
        c3 = count_3cycles(A, n)

        cycles = find_odd_cycles(A, n)
        alpha_1 = len(cycles)
        alphas = count_independent_sets(cycles, len(cycles))

        H_ocf = sum(2**k * alphas.get(k, 0) for k in range(len(cycles)+1))
        if H_ocf != H:
            print(f"  OCF MISMATCH: H={H}, OCF={H_ocf}")

        ocf_data.append({
            'H': H, 'c3': c3, 'alpha_1': alpha_1,
            'alpha_2': alphas.get(2, 0),
            'alpha_3': alphas.get(3, 0),
            'alphas': alphas, 'scores': scores(A, n),
            'n_cycles_by_len': {}
        })

        # Count cycles by length
        for c in cycles:
            l = len(c)
            ocf_data[-1]['n_cycles_by_len'][l] = ocf_data[-1]['n_cycles_by_len'].get(l, 0) + 1

        count += 1

    by_H = defaultdict(list)
    for d in ocf_data:
        by_H[d['H']].append(d)

    print(f"\n  OCF decomposition for each H value:")
    print(f"  {'H':>4}  {'#':>5}  {'α₁':>4}  {'α₂':>4}  {'α₃':>4}  {'3-cyc':>5}  {'5-cyc':>5}  {'scores':>20}")
    for H in sorted(by_H.keys()):
        items = by_H[H]
        a1s = sorted(set(d['alpha_1'] for d in items))
        a2s = sorted(set(d['alpha_2'] for d in items))
        a3s = sorted(set(d['alpha_3'] for d in items))
        c3s = sorted(set(d['n_cycles_by_len'].get(3, 0) for d in items))
        c5s = sorted(set(d['n_cycles_by_len'].get(5, 0) for d in items))
        scs = sorted(set(d['scores'] for d in items))
        print(f"  {H:4d}  {len(items):5d}  {a1s}  {a2s}  {a3s}  {c3s}  {c5s}  {scs}")

    # The H=7 impossibility
    print(f"\n  ═══ WHY H=7 IS FORBIDDEN ═══")
    print(f"  H = 1 + 2α₁ + 4α₂ + 8α₃ + ...")
    print(f"  H = 7 requires α₁ + 2α₂ + 4α₃ + ... = 3")
    print(f"  Possible decompositions: (α₁,α₂,α₃) = (3,0,0) or (1,1,0)")
    print()

    # Check (3,0,0): exactly 3 odd cycles, none vertex-disjoint
    a1_3 = [d for d in ocf_data if d['alpha_1'] == 3]
    print(f"  Tournaments with α₁=3: {len(a1_3)}")
    if a1_3:
        a2_vals = sorted(set(d['alpha_2'] for d in a1_3))
        print(f"    Their α₂ values: {a2_vals}")
        for a2 in a2_vals:
            subset = [d for d in a1_3 if d['alpha_2'] == a2]
            Hs = sorted(set(d['H'] for d in subset))
            print(f"    α₁=3, α₂={a2}: H = {Hs} ({len(subset)} tournaments)")
            # Cycle decomposition
            cyc_lens = sorted(set(str(d['n_cycles_by_len']) for d in subset))
            print(f"      Cycle breakdown: {cyc_lens}")

    # Check (1,1,0): 1 cycle total??? No, α₁=1, α₂=1 means 1 odd cycle...
    # Wait: α₁ = # independent sets of size 1 = # odd cycles
    # α₂ = # independent sets of size 2 = # pairs of vertex-disjoint odd cycles
    # H = 1 + 2·1 + 4·1 = 7 needs (α₁, α₂) = (1, 1)
    # But if α₁ = 1 and α₂ = 1, that means 1 cycle and 1 disjoint pair...
    # No! α₂ = 1 means there's 1 pair of vertex-disjoint cycles. But with only 1 cycle, there can't be a pair.
    # So actually α₂ > 0 requires α₁ ≥ 2.
    # For (α₁, α₂) = (1, 1): impossible (need 2 cycles for a pair).
    # The only decomposition: (3, 0, 0).
    print(f"\n  The ONLY valid decomposition for H=7: α₁=3, α₂=0, α₃=0")
    print(f"  Meaning: exactly 3 odd cycles, ALL pairwise sharing vertices")

    a1_3_a2_0 = [d for d in ocf_data if d['alpha_1'] == 3 and d['alpha_2'] == 0]
    if a1_3_a2_0:
        print(f"  Found {len(a1_3_a2_0)} tournaments with this decomposition")
        print(f"  Their H values: {sorted(set(d['H'] for d in a1_3_a2_0))}")
        print(f"  BUT their H should be 1+6=7... checking OCF...")
        for d in a1_3_a2_0[:5]:
            print(f"    alphas={d['alphas']}, H={d['H']}")
    else:
        print(f"  ✓ NO tournament has (α₁=3, α₂=0)!")
        print(f"  → 3 odd cycles ALWAYS have at least one disjoint pair at n=5")
        print(f"  → H = 1+2·3+4·α₂ ≥ 1+6+4 = 11")
        print(f"  → H=7, H=9 with α₁=3 are impossible (9 needs (3,0) too)")

    # What about α₁=3 from 3-cycles + 5-cycles?
    print(f"\n  Can α₁=3 be from mixed 3+5 cycles?")
    for d in a1_3[:5] if a1_3 else []:
        print(f"    3-cycles: {d['n_cycles_by_len'].get(3, 0)}, "
              f"5-cycles: {d['n_cycles_by_len'].get(5, 0)}, "
              f"scores: {d['scores']}, α₂={d['alpha_2']}, H={d['H']}")


# ─── 4. THE KEY THEOREM: frustration forces disjoint pairs ───

def frustration_disjointness_theorem():
    """Test the core AFM theorem:
    'In a tournament on n ≥ 5 vertices with ≥ 3 odd cycles,
    at least two must be vertex-disjoint (α₂ ≥ 1).'

    If true, this gives H ≥ 1 + 2·3 + 4·1 = 11 whenever α₁ ≥ 3,
    which PROVES H=7 is impossible for n ≥ 5 and explains it for n ≤ 4 trivially."""

    for n in [4, 5, 6]:
        print(f"\n{'='*80}")
        print(f"  FRUSTRATION FORCES DISJOINTNESS — n = {n}")
        print(f"{'='*80}")

        violations = 0
        total = 0
        data = []

        for mask, A in all_tournaments(n):
            cycles = find_odd_cycles(A, n)
            alpha_1 = len(cycles)
            if alpha_1 >= 3:
                alphas = count_independent_sets(cycles, len(cycles))
                alpha_2 = alphas.get(2, 0)
                total += 1
                if alpha_2 == 0:
                    violations += 1
                    data.append({'alpha_1': alpha_1, 'scores': scores(A, n)})
            if n == 6 and total % 2000 == 0 and total > 0:
                print(f"    ... checked {total} with α₁≥3", end='\r')

        print(f"  Tournaments with α₁ ≥ 3: {total}")
        print(f"  Violations (α₁≥3 but α₂=0): {violations}")
        if violations == 0:
            print(f"  ✓ THEOREM HOLDS at n={n}: α₁ ≥ 3 → α₂ ≥ 1")
            print(f"    → H ≥ 1 + 2·3 + 4·1 = 11 whenever α₁ ≥ 3")
            print(f"    → H=7 = 1+2·3+0 is impossible (the 4α₂ term is always ≥ 4)")
        else:
            print(f"  ✗ VIOLATIONS found:")
            for d in data[:5]:
                print(f"    α₁={d['alpha_1']}, scores={d['scores']}")


if __name__ == "__main__":
    t0 = time.time()

    # 1. Ground state at n=6
    ground_state_n6()

    # 2. Thermodynamics at n=5
    thermodynamics_n5()

    # 3. Forbidden value deep analysis
    forbidden_deep_n5()

    # 4. The key theorem
    frustration_disjointness_theorem()

    t1 = time.time()
    print(f"\n\nTotal time: {t1-t0:.1f}s")
