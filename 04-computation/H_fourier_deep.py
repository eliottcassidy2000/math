#!/usr/bin/env python3
"""
H_fourier_deep.py -- kind-pasteur-2026-03-13-S61

Deep analysis of the H Walsh-Hadamard structure:
1. Z_v as function of score at all n
2. Degree-4 coefficient structure
3. Maximum Fourier degree
4. Connection to overlap weight via the Fourier decomposition
5. The "Vitali dimension" -- how the Fourier spectrum encodes the
   equivalence class structure

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


def compute_Zv(A, n):
    """Compute Z_v for each vertex using the Walsh-Hadamard degree-2 formula."""
    sigma_dict = {}
    for i in range(n):
        for j in range(i+1, n):
            sigma_dict[(i, j)] = 1 if A[i][j] else -1
            sigma_dict[(j, i)] = -(1 if A[i][j] else -1)

    Z = [0] * n
    for v in range(n):
        # Cross terms: a < v < b
        for a in range(v):
            for b in range(v + 1, n):
                Z[v] += sigma_dict[(a, v)] * sigma_dict[(v, b)]
        # Same-side-below: a < b < v (both edges have v as larger vertex)
        for a in range(v):
            for b in range(a + 1, v):
                Z[v] -= sigma_dict[(a, v)] * sigma_dict[(b, v)]
        # Same-side-above: v < a < b (both edges have v as smaller vertex)
        for a in range(v + 1, n):
            for b in range(a + 1, n):
                Z[v] -= sigma_dict[(v, a)] * sigma_dict[(v, b)]
    return Z


# ========================================================================
# ANALYSIS 1: Z_v AS FUNCTION OF SCORE
# ========================================================================
print("=" * 70)
print("ANALYSIS 1: Z_v AS DETERMINISTIC FUNCTION OF SCORE")
print("=" * 70)

for n in range(3, 8):
    m = n * (n - 1) // 2
    total = 1 << m
    half = (n - 1) / 2

    if total > 2**20:
        print(f"\nn={n}: (too large, sampling)")
        import random
        random.seed(42)
        sample_size = min(100000, total)
        z_vs_score = defaultdict(set)
        for _ in range(sample_size):
            bits = random.randint(0, total - 1)
            A = binary_to_tournament(bits, n)
            Z = compute_Zv(A, n)
            scores = [sum(A[v]) for v in range(n)]
            for v in range(n):
                z_vs_score[scores[v]].add(Z[v])
    else:
        z_vs_score = defaultdict(set)
        for bits in range(total):
            A = binary_to_tournament(bits, n)
            Z = compute_Zv(A, n)
            scores = [sum(A[v]) for v in range(n)]
            for v in range(n):
                z_vs_score[scores[v]].add(Z[v])

    print(f"\nn={n}: (n-1)/2 = {half}")
    print(f"  Score | Z_v values | -2(s-m)^2+2 | match?")
    print(f"  ------+------------+-------------+-------")

    all_deterministic = True
    z_formula = {}
    for s in sorted(z_vs_score.keys()):
        z_vals = z_vs_score[s]
        predicted = -2 * (s - half) ** 2 + 2
        is_det = len(z_vals) == 1
        if not is_det:
            all_deterministic = False
        z_formula[s] = z_vals

        z_str = str(sorted(z_vals)) if len(z_vals) <= 5 else f"{len(z_vals)} values"
        match = (len(z_vals) == 1 and abs(list(z_vals)[0] - predicted) < 0.01)
        print(f"  {s:>5d} | {z_str:>10s} | {predicted:>11.2f} | {'YES' if match else 'NO'}")

    print(f"  Z_v deterministic by score? {all_deterministic}")

    # Try to find the actual Z(s) formula
    if all_deterministic:
        z_actual = {s: list(v)[0] for s, v in z_formula.items()}
        scores_list = sorted(z_actual.keys())
        z_list = [z_actual[s] for s in scores_list]

        # Fit polynomial: Z(s) = a*s^2 + b*s + c
        # Using 3 points
        if len(scores_list) >= 3:
            s0, s1, s2 = scores_list[0], scores_list[1], scores_list[2]
            z0, z1, z2 = z_actual[s0], z_actual[s1], z_actual[s2]

            # Vandermonde system
            # z0 = a*s0^2 + b*s0 + c
            # z1 = a*s1^2 + b*s1 + c
            # z2 = a*s2^2 + b*s2 + c

            det = (s0**2 * (s1 - s2) - s1**2 * (s0 - s2) + s2**2 * (s0 - s1))
            if abs(det) > 1e-10:
                a = (z0 * (s1 - s2) - z1 * (s0 - s2) + z2 * (s0 - s1)) / det
                b = (s0**2 * (z1 - z2) - s1**2 * (z0 - z2) + s2**2 * (z0 - z1)) / det
                c = (s0**2 * (s1*z2 - s2*z1) - s1**2 * (s0*z2 - s2*z0) + s2**2 * (s0*z1 - s1*z0)) / det

                # Verify all points
                all_match = True
                for s in scores_list:
                    pred = a * s**2 + b * s + c
                    if abs(pred - z_actual[s]) > 0.01:
                        all_match = False
                        break

                if all_match:
                    # Express in centered form
                    vertex_s = -b / (2 * a) if abs(a) > 1e-10 else 0
                    vertex_z = c - b**2 / (4 * a) if abs(a) > 1e-10 else c
                    print(f"  Z(s) = {a:.4f} * s^2 + {b:.4f} * s + {c:.4f}")
                    print(f"       = {a:.4f} * (s - {vertex_s:.4f})^2 + {vertex_z:.4f}")
                    print(f"  Vertex at s = {vertex_s:.4f} = (n-1)/2 = {half}? "
                          f"{'YES' if abs(vertex_s - half) < 0.01 else 'NO'}")
                else:
                    print(f"  Z(s) is NOT quadratic in s!")
                    # Try cubic or higher
                    print(f"  Z values: {z_list}")


# ========================================================================
# ANALYSIS 2: DEGREE-4 STRUCTURE
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 2: DEGREE-4 FOURIER STRUCTURE")
print("=" * 70)

for n in [5, 6]:
    m = n * (n - 1) // 2
    total = 1 << m

    edge_pairs = []
    for i in range(n):
        for j in range(i+1, n):
            edge_pairs.append((i, j))

    H_vals = {}
    for bits in range(total):
        A = binary_to_tournament(bits, n)
        sigma = tournament_to_sigma(A, n)
        H = count_ham_paths(A, n)
        H_vals[sigma] = H

    print(f"\nn={n} (m={m}):")

    # Compute degree-4 coefficients
    deg4_coeffs = {}
    for S in combinations(range(m), 4):
        coeff = 0
        for sigma, H in H_vals.items():
            chi = sigma[S[0]] * sigma[S[1]] * sigma[S[2]] * sigma[S[3]]
            coeff += H * chi
        coeff /= total
        if abs(coeff) > 1e-10:
            deg4_coeffs[S] = coeff

    print(f"  Number of non-zero degree-4 coefficients: {len(deg4_coeffs)}")

    if deg4_coeffs:
        # Check: are all magnitudes the same?
        mags = set(round(abs(c), 8) for c in deg4_coeffs.values())
        print(f"  Distinct magnitudes: {sorted(mags)}")

        # Analyze the vertex-overlap structure of the 4-edge subsets
        # A 4-edge subset uses 4 edges from C(n,2) edges
        # The "vertex support" is the union of all endpoints
        # Possible support sizes: 4, 5, 6, 7, 8 (at n >= 8)

        support_dist = defaultdict(list)
        for S, c in deg4_coeffs.items():
            edges = [edge_pairs[idx] for idx in S]
            verts = set()
            for e in edges:
                verts.update(e)
            support_dist[len(verts)].append((S, c))

        print(f"\n  Support size distribution:")
        for sz in sorted(support_dist):
            items = support_dist[sz]
            mags = set(round(abs(c), 8) for _, c in items)
            print(f"    |support|={sz}: {len(items)} terms, magnitudes = {sorted(mags)}")

            # Show a few examples
            if len(items) <= 3:
                for S, c in items:
                    edges = [edge_pairs[idx] for idx in S]
                    print(f"      {edges}: coeff = {c:.6f}")
            else:
                for S, c in items[:2]:
                    edges = [edge_pairs[idx] for idx in S]
                    print(f"      {edges}: coeff = {c:.6f}")
                print(f"      ... and {len(items)-2} more")

        # CRITICAL CHECK: do degree-4 terms correspond to 5-vertex subgraphs?
        # A 4-cycle in K_n uses 4 edges and 4 vertices
        # A path-of-length-3 uses 4 edges and 5 vertices
        # etc.

        # Classify by graph structure
        print(f"\n  Graph structure of 4-edge subsets with non-zero coefficient:")
        graph_types = defaultdict(list)
        for S, c in deg4_coeffs.items():
            edges = [edge_pairs[idx] for idx in S]
            verts = set()
            for e in edges:
                verts.update(e)

            # Compute degree sequence within the 4-edge subgraph
            deg_seq = defaultdict(int)
            for e in edges:
                deg_seq[e[0]] += 1
                deg_seq[e[1]] += 1

            # Sort degree sequence for classification
            deg_type = tuple(sorted(deg_seq.values(), reverse=True))
            graph_types[deg_type].append((S, c, len(verts)))

        for dtype, items in sorted(graph_types.items()):
            support = items[0][2]
            print(f"    Degree type {dtype} (support={support}): "
                  f"{len(items)} terms, |coeff| = {abs(items[0][1]):.6f}")


# ========================================================================
# ANALYSIS 3: MAXIMUM FOURIER DEGREE
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 3: MAXIMUM NON-ZERO FOURIER DEGREE")
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

    max_nonzero_deg = 0
    for k in range(0, m + 1, 2):
        has_nonzero = False
        for S in combinations(range(m), k):
            coeff = 0
            for sigma, H in H_vals.items():
                chi = 1
                for idx in S:
                    chi *= sigma[idx]
                coeff += H * chi
            coeff /= total
            if abs(coeff) > 1e-10:
                has_nonzero = True
                break
        if has_nonzero:
            max_nonzero_deg = k

    predicted_max = 2 * ((n - 1) // 2)
    print(f"  n={n}: max non-zero even degree = {max_nonzero_deg}, "
          f"predicted 2*floor((n-1)/2) = {predicted_max}, "
          f"match = {max_nonzero_deg == predicted_max}")


# ========================================================================
# ANALYSIS 4: H_2 EXACT FORMULA IN TERMS OF SCORES
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 4: H_2 AS EXACT FUNCTION OF SCORE SEQUENCE")
print("=" * 70)

for n in range(3, 8):
    m = n * (n - 1) // 2
    total = 1 << m
    c2 = math.factorial(n - 2) / (2 ** (n - 2))
    EH = math.factorial(n) / (2 ** (n - 1))

    if total > 2**20:
        print(f"\nn={n}: (too large)")
        continue

    # For each tournament, compute H_2 and relate to score sequence
    score_seq_to_H2 = defaultdict(set)
    score_seq_to_sumZ = defaultdict(set)

    for bits in range(total):
        A = binary_to_tournament(bits, n)
        Z = compute_Zv(A, n)
        sumZ = sum(Z)
        H2 = c2 * sumZ

        scores = tuple(sorted([sum(A[v]) for v in range(n)]))
        score_seq_to_H2[scores].add(round(H2, 6))
        score_seq_to_sumZ[scores].add(sumZ)

    print(f"\nn={n}: c_2 = {c2:.6f}, E[H] = {EH:.2f}")
    print(f"  Score sequence -> H_2 deterministic? ", end="")
    deterministic = all(len(v) == 1 for v in score_seq_to_H2.values())
    print(f"{'YES' if deterministic else 'NO'}")

    if deterministic:
        print(f"  Score seq            | sum(Z) | H_2       | H_0+H_2")
        print(f"  ---------------------+--------+-----------+--------")
        for ss in sorted(score_seq_to_sumZ.keys()):
            sz = list(score_seq_to_sumZ[ss])[0]
            h2 = list(score_seq_to_H2[ss])[0]
            print(f"  {str(ss):>21s} | {sz:>6d} | {h2:>9.4f} | {EH + h2:>7.2f}")

        # Check: is sum(Z) = sum_v Z(s_v) where Z(s) = -2(s-m)^2 + 2?
        print(f"\n  Verify: sum(Z) = sum_v [-2(s_v - (n-1)/2)^2 + 2]?")
        half = (n - 1) / 2
        for ss in sorted(score_seq_to_sumZ.keys()):
            actual_sumZ = list(score_seq_to_sumZ[ss])[0]
            predicted_sumZ = sum(-2 * (s - half)**2 + 2 for s in ss)
            match = abs(actual_sumZ - predicted_sumZ) < 0.01
            print(f"    {ss}: actual={actual_sumZ}, predicted={predicted_sumZ:.0f}, match={match}")
    else:
        # H_2 is NOT determined by score sequence alone!
        print(f"  IMPORTANT: H_2 depends on MORE than just score sequence!")
        for ss in sorted(score_seq_to_H2.keys()):
            vals = score_seq_to_H2[ss]
            if len(vals) > 1:
                print(f"    {ss}: H_2 values = {sorted(vals)}")


# ========================================================================
# ANALYSIS 5: OVERLAP WEIGHT IN FOURIER DOMAIN
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 5: OVERLAP WEIGHT IN FOURIER DOMAIN")
print("=" * 70)

n = 5
m = n * (n - 1) // 2
total = 1 << m
c2 = math.factorial(n - 2) / (2 ** (n - 2))
EH = math.factorial(n) / (2 ** (n - 1))

print(f"\nn={n}: Understanding how cycle overlap connects to Fourier structure")

# For each tournament, compute:
# 1. H and its decomposition H = H_0 + H_2 + H_4
# 2. alpha decomposition H = 1 + 2*N + 4*alpha_2 + ...
# 3. Check: does H_2 correspond to the 2*N term?

for bits in [0, 76, 120, 341, 682, 1023]:
    if bits >= total:
        continue
    A = binary_to_tournament(bits, n)
    H = count_ham_paths(A, n)
    scores = sorted([sum(A[v]) for v in range(n)])

    Z = compute_Zv(A, n)
    sumZ = sum(Z)
    H_2 = c2 * sumZ
    H_4 = H - EH - H_2

    # Count odd cycles
    c3_sets = []
    c5_directed = 0
    for a, b, c in combinations(range(n), 3):
        if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
            c3_sets.append(frozenset([a, b, c]))

    c3_sets = list(set(c3_sets))
    c3 = len(c3_sets)

    # 5-cycles
    c5 = 0
    for subset in combinations(range(n), 5):
        verts = list(subset)
        # Count directed Ham cycles on this subset
        sub_A = [[A[i][j] for j in verts] for i in verts]
        # Simple DP
        k = 5
        dp = {}
        dp[(1, 0)] = 1
        for mask_val in range(1, 1 << k):
            for v in range(k):
                if not (mask_val & (1 << v)):
                    continue
                key = (mask_val, v)
                if key not in dp or dp[key] == 0:
                    continue
                for w in range(k):
                    if mask_val & (1 << w):
                        continue
                    if sub_A[v][w]:
                        nk = (mask_val | (1 << w), w)
                        dp[nk] = dp.get(nk, 0) + dp[key]
        full = (1 << k) - 1
        for v in range(1, k):
            if (full, v) in dp and sub_A[v][0]:
                c5 += dp[(full, v)]

    N = c3 + c5  # total directed odd cycles (c3 are undirected, each gives 2 directed)
    # Wait, c3 counts vertex sets, each vertex set has exactly 1 undirected cycle = 2 directed
    # For OCF: we count DIRECTED cycles. Each 3-vertex set with a cycle has 2 directed cycles.
    # Actually no: each 3-vertex set has exactly ONE directed 3-cycle (the other direction is the
    # SAME undirected cycle traversed backward). For tournament: each {a,b,c} with a cycle
    # has exactly 1 directed cycle (either a->b->c->a OR a->c->b->a, not both).
    # Actually for a tournament on 3 vertices: either it's transitive (0 cycles) or it has
    # exactly 2 directed 3-cycles (clockwise and counterclockwise). These are DIFFERENT directed cycles.
    # But for the independence polynomial / OCF, we count UNDIRECTED odd cycles (vertex sets).
    # alpha_1 = number of vertex sets that support a directed odd cycle.
    # So N = c3 (vertex sets) + c5 (vertex sets with 5-cycle)

    # Actually let me recount. At n=5:
    # c3 vertex sets: each contributes 1 to alpha_1
    # c5: at n=5, C(5,5)=1, and the full tournament either has Ham cycles or not
    # A 5-vertex set supports a 5-cycle iff the tournament restricted to those 5 vertices
    # has a Hamiltonian cycle. For our c5 count, c5 = number of directed Ham cycles / 2
    # (since each undirected cycle gives 2 directed ones).
    # Actually for 5-cycles: a 5-vertex tournament can have 0, 1, 2, or 3 directed 5-cycles.
    # Wait, it can have 0, 2, 4, 6, ... (even number, since each undirected cycle gives 2).
    # No wait, directed cycles on 5 vertices: each vertex ordering gives a cycle.
    # 5!/5 = 24 orderings up to rotation. Each can go CW or CCW = 12 undirected.
    # A tournament selects some of these. The count is even by complementation? Not necessarily.

    # Let me just count alpha_1 properly
    alpha_1_cycles = []  # vertex sets supporting at least one directed odd cycle
    for k in range(3, n + 1, 2):
        for subset in combinations(range(n), k):
            verts = list(subset)
            k_sub = len(verts)
            sub_A = [[0]*k_sub for _ in range(k_sub)]
            for ii in range(k_sub):
                for jj in range(k_sub):
                    sub_A[ii][jj] = A[verts[ii]][verts[jj]]

            # Check if there's at least one directed Ham cycle
            dp = {}
            dp[(1, 0)] = 1
            for mask_val in range(1, 1 << k_sub):
                for v in range(k_sub):
                    if not (mask_val & (1 << v)):
                        continue
                    key = (mask_val, v)
                    if key not in dp or dp[key] == 0:
                        continue
                    for w in range(k_sub):
                        if mask_val & (1 << w):
                            continue
                        if sub_A[v][w]:
                            nk = (mask_val | (1 << w), w)
                            dp[nk] = dp.get(nk, 0) + dp[key]

            full = (1 << k_sub) - 1
            has_cycle = False
            for v in range(1, k_sub):
                if (full, v) in dp and sub_A[v][0]:
                    has_cycle = True
                    break
            if has_cycle:
                alpha_1_cycles.append(frozenset(subset))

    N_alpha = len(alpha_1_cycles)

    # Disjoint pairs
    alpha_2 = 0
    for i in range(len(alpha_1_cycles)):
        for j in range(i+1, len(alpha_1_cycles)):
            if not (alpha_1_cycles[i] & alpha_1_cycles[j]):
                alpha_2 += 1

    H_ocf = 1 + 2 * N_alpha + 4 * alpha_2

    print(f"\n  bits={bits}: H={H}, scores={scores}")
    print(f"    OCF: alpha_1={N_alpha}, alpha_2={alpha_2}, H=1+{2*N_alpha}+{4*alpha_2}={H_ocf}")
    print(f"    Fourier: H_0={EH:.1f}, H_2={H_2:.2f}, H_4={H_4:.2f}")
    print(f"    Note: H_0+H_2 = {EH+H_2:.2f}, residual = {H_4:.2f}")

    # Compare: does H_2 relate to N linearly?
    print(f"    2*N = {2*N_alpha}, H_2 = {H_2:.2f}")
    print(f"    4*alpha_2 = {4*alpha_2}, H_4 = {H_4:.2f}")


print("\n" + "=" * 70)
print("DONE.")
