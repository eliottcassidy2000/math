#!/usr/bin/env python3
"""
vitali_210_bridge.py -- kind-pasteur-2026-03-13-S61

THE {2,1,0} BRIDGE: How overlap weights encode the Vitali structure

The key insight from THM-164: H(T) lives on a ~2-dimensional manifold
embedded in the high-dimensional tournament space {-1,+1}^m.

The {2,1,0} values arise in THREE distinct contexts:
  - OVERLAP WEIGHTS: |V(C_i) ∩ V(C_j)| in {0,1,2} for 3-cycles
  - FOURIER DEGREES: Only degrees 0, 2, 4 (mapped to indices 0, 1, 2) appear significantly
  - SCORE STRUCTURE: At n=3, scores are in {0,1,2}

THIS SCRIPT explores whether these three {2,1,0} structures are
THE SAME mathematical object seen from three different angles.

Specifically:
1. The overlap weight matrix W[i,j] for 3-cycles has eigenstructure
   related to the Fourier decomposition
2. The "average overlap weight" is connected to H_2 (score part)
3. The "overlap weight variance" is connected to H_4 (cycle part)
4. At n=3: everything collapses to a SINGLE number (overlap weight = 0 trivially)

Author: kind-pasteur-2026-03-13-S61
"""

import math
from itertools import combinations, permutations
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


def count_directed_ham_cycles_subset(A, subset):
    """Count directed Hamiltonian cycles on a vertex subset."""
    verts = list(subset)
    k = len(verts)
    if k < 3:
        return 0
    if k == 3:
        a, b, c = verts
        return (A[a][b] * A[b][c] * A[c][a]) + (A[a][c] * A[c][b] * A[b][a])

    dp = {}
    dp[(1, 0)] = 1
    for mask in range(1, 1 << k):
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nk = (mask | (1 << w), w)
                    dp[nk] = dp.get(nk, 0) + dp[key]
    full = (1 << k) - 1
    total = 0
    for v in range(1, k):
        if (full, v) in dp and A[verts[v]][verts[0]]:
            total += dp[(full, v)]
    return total


def compute_H2(A, n):
    c2 = math.factorial(n - 2) / (2 ** (n - 2))
    scores = [sum(A[v]) for v in range(n)]
    half = (n - 1) / 2
    Z = [-2 * (s - half)**2 + half for s in scores]
    return c2 * sum(Z)


def get_cycle_data(A, n, max_k=None):
    """Get all directed odd cycles with vertex sets and lengths."""
    if max_k is None:
        max_k = n
    # Collect: (vertex_set, length, count_of_directed_cycles)
    cycles = []
    for k in range(3, max_k + 1, 2):
        for subset in combinations(range(n), k):
            cnt = count_directed_ham_cycles_subset(A, subset)
            if cnt > 0:
                for _ in range(cnt):
                    cycles.append((frozenset(subset), k))
    return cycles


# ========================================================================
# ANALYSIS 1: THE {2,1,0} OVERLAP STRUCTURE AT n=3,4,5,6
# ========================================================================
print("=" * 70)
print("ANALYSIS 1: OVERLAP WEIGHT MATRIX SPECTRUM")
print("=" * 70)

for n in range(3, 7):
    m = n * (n - 1) // 2
    total = 1 << m
    EH = math.factorial(n) / (2 ** (n - 1))

    print(f"\n--- n={n} (m={m}, {total} tournaments) ---")

    # Collect statistics across all tournaments
    overlap_spectrum = defaultdict(lambda: defaultdict(int))  # H -> overlap_val -> count
    overlap_mean_by_H = defaultdict(list)
    overlap_var_by_H = defaultdict(list)
    alpha_by_H = defaultdict(list)  # H -> (alpha_1, alpha_2)

    for bits in range(total):
        A = binary_to_tournament(bits, n)
        H = count_ham_paths(A, n)

        # Get all odd-cycle vertex sets (not counting directed multiplicity)
        cycle_vsets = set()
        for k in range(3, n + 1, 2):
            for subset in combinations(range(n), k):
                cnt = count_directed_ham_cycles_subset(A, subset)
                if cnt > 0:
                    cycle_vsets.add(frozenset(subset))

        cycle_list = list(cycle_vsets)
        n_cyc = len(cycle_list)
        alpha_1 = n_cyc

        # Overlap weights (between distinct cycle vertex sets)
        ov_vals = []
        alpha_2 = 0
        for i in range(n_cyc):
            for j in range(i+1, n_cyc):
                ov = len(cycle_list[i] & cycle_list[j])
                ov_vals.append(ov)
                overlap_spectrum[H][ov] += 1
                if ov == 0:
                    alpha_2 += 1

        if ov_vals:
            mean_ov = sum(ov_vals) / len(ov_vals)
            var_ov = sum((x - mean_ov)**2 for x in ov_vals) / len(ov_vals)
        else:
            mean_ov = 0
            var_ov = 0

        overlap_mean_by_H[H].append(mean_ov)
        overlap_var_by_H[H].append(var_ov)
        alpha_by_H[H].append((alpha_1, alpha_2))

    # Print summary
    print(f"  {'H':>5s} | {'count':>6s} | {'alpha1':>7s} | {'alpha2':>7s} | {'avg_ov':>8s} | {'var_ov':>8s} | overlap dist")
    print(f"  {'':->5s}-+-{'':->6s}-+-{'':->7s}-+-{'':->7s}-+-{'':->8s}-+-{'':->8s}-+{'-'*20}")

    for H in sorted(overlap_spectrum.keys()):
        count = len(overlap_mean_by_H[H])
        a1s = set(a[0] for a in alpha_by_H[H])
        a2s = set(a[1] for a in alpha_by_H[H])
        avg_ov = sum(overlap_mean_by_H[H]) / count
        avg_var = sum(overlap_var_by_H[H]) / count

        ov_dist_str = " ".join(f"{k}:{v//count}" for k, v in sorted(overlap_spectrum[H].items()))

        a1_str = str(sorted(a1s)[0]) if len(a1s) == 1 else str(sorted(a1s))
        a2_str = str(sorted(a2s)[0]) if len(a2s) == 1 else str(sorted(a2s))

        print(f"  {H:>5d} | {count:>6d} | {a1_str:>7s} | {a2_str:>7s} | {avg_ov:>8.3f} | {avg_var:>8.3f} | {ov_dist_str}")


# ========================================================================
# ANALYSIS 2: CONNECTION TO FOURIER COMPONENTS
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 2: OVERLAP WEIGHT vs FOURIER COMPONENTS")
print("=" * 70)

for n in [5, 6]:
    m = n * (n - 1) // 2
    total = 1 << m
    EH = math.factorial(n) / (2 ** (n - 1))

    print(f"\n--- n={n} ---")

    # For each tournament, compute:
    # 1. H, H_2, H_4
    # 2. alpha_1, alpha_2 (from vertex-set overlap)
    # 3. Mean and variance of overlap weights

    data_rows = []
    for bits in range(total):
        A = binary_to_tournament(bits, n)
        H = count_ham_paths(A, n)
        H2 = compute_H2(A, n)
        H4 = H - EH - H2
        scores = tuple(sorted([sum(A[v]) for v in range(n)]))

        # Cycle vertex sets
        cycle_vsets = []
        for k in range(3, n + 1, 2):
            for subset in combinations(range(n), k):
                cnt = count_directed_ham_cycles_subset(A, subset)
                if cnt > 0:
                    cycle_vsets.append(frozenset(subset))

        alpha_1 = len(cycle_vsets)
        alpha_2 = 0
        total_ov = 0
        total_ov_sq = 0
        n_pairs = 0
        ov_dist = defaultdict(int)

        for i in range(len(cycle_vsets)):
            for j in range(i+1, len(cycle_vsets)):
                ov = len(cycle_vsets[i] & cycle_vsets[j])
                if ov == 0:
                    alpha_2 += 1
                total_ov += ov
                total_ov_sq += ov * ov
                n_pairs += 1
                ov_dist[ov] += 1

        mean_ov = total_ov / n_pairs if n_pairs > 0 else 0
        var_ov = (total_ov_sq / n_pairs - mean_ov**2) if n_pairs > 0 else 0

        data_rows.append({
            'bits': bits, 'H': H, 'H2': H2, 'H4': round(H4, 4),
            'scores': scores, 'alpha_1': alpha_1, 'alpha_2': alpha_2,
            'mean_ov': mean_ov, 'var_ov': var_ov, 'ov_dist': dict(ov_dist),
            'n_pairs': n_pairs
        })

    # Correlations
    # H4 vs mean_ov, H4 vs var_ov, H4 vs alpha_2, H4 vs (alpha_1 - alpha_2)
    def corr(xs, ys):
        n_d = len(xs)
        mx = sum(xs) / n_d
        my = sum(ys) / n_d
        cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / n_d
        vx = sum((x - mx)**2 for x in xs) / n_d
        vy = sum((y - my)**2 for y in ys) / n_d
        if vx > 0 and vy > 0:
            return cov / (vx**0.5 * vy**0.5)
        return 0

    H4s = [d['H4'] for d in data_rows]
    H2s = [d['H2'] for d in data_rows]
    a1s = [d['alpha_1'] for d in data_rows]
    a2s = [d['alpha_2'] for d in data_rows]
    mov = [d['mean_ov'] for d in data_rows]
    vov = [d['var_ov'] for d in data_rows]

    print(f"\n  Global correlations:")
    print(f"    corr(H_2, alpha_1) = {corr(H2s, a1s):.4f}")
    print(f"    corr(H_2, alpha_2) = {corr(H2s, a2s):.4f}")
    print(f"    corr(H_4, alpha_1) = {corr(H4s, a1s):.4f}")
    print(f"    corr(H_4, alpha_2) = {corr(H4s, a2s):.4f}")
    print(f"    corr(H_4, mean_overlap) = {corr(H4s, mov):.4f}")
    print(f"    corr(H_4, var_overlap) = {corr(H4s, vov):.4f}")

    # Within score classes
    print(f"\n  Within-score-class correlations (H_4 vs cycle metrics):")
    by_score = defaultdict(list)
    for d in data_rows:
        by_score[d['scores']].append(d)

    for sc in sorted(by_score.keys()):
        group = by_score[sc]
        if len(set(d['H4'] for d in group)) <= 1:
            continue  # No H4 variation

        h4 = [d['H4'] for d in group]
        a1 = [d['alpha_1'] for d in group]
        a2 = [d['alpha_2'] for d in group]
        mo = [d['mean_ov'] for d in group]
        vo = [d['var_ov'] for d in group]

        c_a1 = corr(h4, a1)
        c_a2 = corr(h4, a2)
        c_mo = corr(h4, mo)
        c_vo = corr(h4, vo)

        print(f"    Score {sc}: corr(H4,a1)={c_a1:.3f}, corr(H4,a2)={c_a2:.3f}, "
              f"corr(H4,mean_ov)={c_mo:.3f}, corr(H4,var_ov)={c_vo:.3f}")


# ========================================================================
# ANALYSIS 3: THE {2,1,0} ENCODING — IS OVERLAP WEIGHT A TERNARY CODE?
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 3: THE {2,1,0} ENCODING AS TERNARY CODE")
print("=" * 70)

n = 5
m = n * (n - 1) // 2
total = 1 << m
EH = math.factorial(n) / (2 ** (n - 1))

print(f"""
At n={n}: 3-cycle vertex sets have pairwise overlaps in {{0,1,2}}.
This naturally defines a TERNARY CODE on C(alpha_1, 2) symbols.

The code word for a tournament T is the vector of overlaps between
all pairs of 3-cycle vertex sets, taking values in {{0,1,2}}.

Question: Is this ternary code a complete invariant for H_4?
If so, the overlap weight distribution ALONE determines H.
""")

# For each tournament, compute the overlap "code word"
by_score_class = defaultdict(list)
for bits in range(total):
    A = binary_to_tournament(bits, n)
    scores = tuple(sorted([sum(A[v]) for v in range(n)]))

    # 3-cycle vertex sets
    c3_sets = []
    for a, b, c in combinations(range(n), 3):
        cnt = count_directed_ham_cycles_subset(A, [a, b, c])
        if cnt > 0:
            c3_sets.append(frozenset([a, b, c]))

    # Overlap code: sorted list of overlap values between all pairs
    # (order doesn't matter for now — just the multiset)
    ov_code = []
    for i in range(len(c3_sets)):
        for j in range(i+1, len(c3_sets)):
            ov_code.append(len(c3_sets[i] & c3_sets[j]))
    ov_code = tuple(sorted(ov_code))

    H = count_ham_paths(A, n)
    H2 = compute_H2(A, n)
    H4 = round(H - EH - H2, 4)

    by_score_class[scores].append({
        'bits': bits, 'H': H, 'H4': H4,
        'alpha_1': len(c3_sets), 'ov_code': ov_code,
        'c3_count': len(c3_sets)
    })

for sc in sorted(by_score_class.keys()):
    group = by_score_class[sc]
    if len(set(d['H4'] for d in group)) <= 1:
        continue

    print(f"\n  Score {sc}:")
    # Group by ov_code
    by_code = defaultdict(lambda: {'H4_vals': set(), 'count': 0})
    for d in group:
        by_code[d['ov_code']]['H4_vals'].add(d['H4'])
        by_code[d['ov_code']]['count'] += 1

    is_complete = all(len(v['H4_vals']) == 1 for v in by_code.values())
    print(f"    Overlap code determines H_4? {'YES' if is_complete else 'NO'}")
    print(f"    Number of distinct overlap codes: {len(by_code)}")

    for code in sorted(by_code.keys()):
        info = by_code[code]
        print(f"      code {code}: H4={sorted(info['H4_vals'])}, count={info['count']}")


# ========================================================================
# ANALYSIS 4: THE FUNDAMENTAL BRIDGE — OVERLAP EIGENVALUES
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 4: OVERLAP WEIGHT MATRIX EIGENSTRUCTURE")
print("=" * 70)

print("""
The overlap weight matrix W[i,j] = |V(C_i) cap V(C_j)| is a symmetric
matrix on the cycle space. Its eigenvalues encode the "spectral geometry"
of the cycle interaction network.

Conjecture: The eigenvalues of W relate to the Fourier degrees of H:
  - The largest eigenvalue relates to H_0 (global average)
  - The second eigenvalue relates to H_2 (score structure)
  - The third eigenvalue relates to H_4 (cycle structure)
""")

for n in [5, 6]:
    m = n * (n - 1) // 2
    total = 1 << m
    EH = math.factorial(n) / (2 ** (n - 1))

    print(f"\n--- n={n} ---")

    # Collect overlap matrices and their spectra for representative tournaments
    eig_data = defaultdict(list)

    for bits in range(total):
        A = binary_to_tournament(bits, n)
        H = count_ham_paths(A, n)
        scores = tuple(sorted([sum(A[v]) for v in range(n)]))

        # Cycle vertex sets (include all odd-length cycles)
        cycle_vsets = []
        for k in range(3, n + 1, 2):
            for subset in combinations(range(n), k):
                cnt = count_directed_ham_cycles_subset(A, subset)
                if cnt > 0:
                    cycle_vsets.append(frozenset(subset))

        if len(cycle_vsets) < 2:
            continue

        # Build overlap matrix
        nc = len(cycle_vsets)
        W = [[0]*nc for _ in range(nc)]
        for i in range(nc):
            for j in range(nc):
                if i != j:
                    W[i][j] = len(cycle_vsets[i] & cycle_vsets[j])

        # Compute eigenvalues (simple power iteration not needed — use characteristic poly for small matrices)
        # For small nc, compute trace and trace(W^2) to get spectral info
        tr_W = sum(W[i][i] for i in range(nc))  # Should be 0 (diagonal is 0)
        tr_W2 = sum(W[i][j] * W[j][i] for i in range(nc) for j in range(nc))
        # tr_W2 = sum of squared eigenvalues = 2 * sum of squared off-diagonal entries
        # (since W is symmetric)

        # Frobenius norm squared = sum of all W[i][j]^2
        frob_sq = sum(W[i][j]**2 for i in range(nc) for j in range(nc))

        # Sum of all entries
        entry_sum = sum(W[i][j] for i in range(nc) for j in range(nc))
        mean_entry = entry_sum / (nc * (nc - 1)) if nc > 1 else 0

        eig_data[H].append({
            'nc': nc, 'tr_W2': tr_W2, 'frob_sq': frob_sq,
            'entry_sum': entry_sum, 'mean_entry': round(mean_entry, 4),
            'scores': scores
        })

    print(f"  {'H':>5s} | {'nc':>4s} | {'sum_W':>8s} | {'mean_W':>8s} | {'Frob^2':>8s} | {'scores':>20s}")
    print(f"  {'':->5s}-+-{'':->4s}-+-{'':->8s}-+-{'':->8s}-+-{'':->8s}-+-{'':->20s}")

    for H in sorted(eig_data.keys()):
        group = eig_data[H]
        nc = set(d['nc'] for d in group)
        sums = set(d['entry_sum'] for d in group)
        means = set(d['mean_entry'] for d in group)
        frobs = set(d['frob_sq'] for d in group)
        scs = set(d['scores'] for d in group)

        nc_str = str(sorted(nc)) if len(nc) > 1 else str(sorted(nc)[0])
        sum_str = str(sorted(sums)) if len(sums) > 1 else str(sorted(sums)[0])
        mean_str = str(sorted(means)[:3]) if len(means) > 1 else str(sorted(means)[0])
        frob_str = str(sorted(frobs)[:3]) if len(frobs) > 1 else str(sorted(frobs)[0])

        print(f"  {H:>5d} | {nc_str:>4s} | {sum_str:>8s} | {mean_str:>8s} | {frob_str:>8s} | {str(sorted(scs)[:2])}")


# ========================================================================
# ANALYSIS 5: OVERLAP WEIGHT DETERMINES INDEPENDENCE POLYNOMIAL
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 5: FROM OVERLAP WEIGHT TO INDEPENDENCE POLYNOMIAL")
print("=" * 70)

print("""
The OCF gives H(T) = I(Omega(T), 2) = sum_{k>=0} alpha_k * 2^k.

The conflict graph Omega has adjacency defined by overlap > 0.
So the ternary structure {0,1,2} is COLLAPSED to binary {0, >0} = {independent, adjacent}.

Question: Does the FULL overlap weight (not just binary) carry MORE information?
Specifically, do two tournaments with the same Omega adjacency but different
overlap weight distributions have the same H?

Answer: YES, by definition — H depends only on Omega, not on the weights.
But the weight distribution may help PREDICT H or explain its structure.
""")

n = 5
m = n * (n - 1) // 2
total = 1 << m
EH = math.factorial(n) / (2 ** (n - 1))

# Check: do different overlap distributions give same H?
# Group by H and check if all have same Omega adjacency
omega_by_H = defaultdict(set)
weight_by_H = defaultdict(set)

for bits in range(total):
    A = binary_to_tournament(bits, n)
    H = count_ham_paths(A, n)

    cycle_vsets = []
    for k in range(3, n + 1, 2):
        for subset in combinations(range(n), k):
            cnt = count_directed_ham_cycles_subset(A, subset)
            if cnt > 0:
                cycle_vsets.append(frozenset(subset))

    nc = len(cycle_vsets)

    # Binary adjacency (Omega)
    omega_edges = set()
    for i in range(nc):
        for j in range(i+1, nc):
            if cycle_vsets[i] & cycle_vsets[j]:
                omega_edges.add((i, j))

    # Weight distribution (multiset of overlap values)
    weight_dist = []
    for i in range(nc):
        for j in range(i+1, nc):
            weight_dist.append(len(cycle_vsets[i] & cycle_vsets[j]))
    weight_dist = tuple(sorted(weight_dist))

    omega_by_H[H].add(frozenset(omega_edges))
    weight_by_H[H].add(weight_dist)

print(f"\nn={n}: Omega and weight diversity by H value:")
print(f"  {'H':>5s} | {'#Omega types':>12s} | {'#weight types':>14s}")
print(f"  {'':->5s}-+-{'':->12s}-+-{'':->14s}")
for H in sorted(omega_by_H.keys()):
    print(f"  {H:>5d} | {len(omega_by_H[H]):>12d} | {len(weight_by_H[H]):>14d}")


# ========================================================================
# ANALYSIS 6: THE DEEP {2,1,0} — TERNARY INFORMATION CONTENT
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 6: TERNARY INFORMATION CONTENT")
print("=" * 70)

print("""
The TERNARY information (overlap=0, 1, or 2) is RICHER than the binary
Omega graph (overlap=0 vs >0). How much richer?

For each tournament, define:
  - I_binary = # distinct Omega adjacency classes (up to relabeling)
  - I_ternary = # distinct overlap weight distributions
  - I_H = # distinct H values

The hierarchy: I_H <= I_ternary ?? I_binary

If I_ternary > I_binary, then the ternary structure carries information
BEYOND what Omega captures — but that information doesn't affect H!
""")

for n in [5, 6]:
    m = n * (n - 1) // 2
    total_t = 1 << m
    EH_loc = math.factorial(n) / (2 ** (n - 1))

    scores_to_info = defaultdict(lambda: {
        'omega_types': set(), 'weight_types': set(), 'H_vals': set()
    })

    for bits in range(total_t):
        A = binary_to_tournament(bits, n)
        H = count_ham_paths(A, n)
        scores = tuple(sorted([sum(A[v]) for v in range(n)]))

        # Cycle vertex sets
        cycle_vsets = []
        for k in range(3, n + 1, 2):
            for subset in combinations(range(n), k):
                cnt = count_directed_ham_cycles_subset(A, subset)
                if cnt > 0:
                    cycle_vsets.append(frozenset(subset))

        nc = len(cycle_vsets)

        # Binary adjacency signature
        adj_sig = []
        weight_sig = []
        for i in range(nc):
            for j in range(i+1, nc):
                ov = len(cycle_vsets[i] & cycle_vsets[j])
                adj_sig.append(1 if ov > 0 else 0)
                weight_sig.append(ov)

        scores_to_info[scores]['omega_types'].add(tuple(sorted(adj_sig)))
        scores_to_info[scores]['weight_types'].add(tuple(sorted(weight_sig)))
        scores_to_info[scores]['H_vals'].add(H)

    print(f"\nn={n}: Ternary vs binary information by score class:")
    print(f"  {'Score':>25s} | {'#H':>4s} | {'#binary':>8s} | {'#ternary':>9s} | {'ternary>binary?':>15s}")
    print(f"  {'':->25s}-+-{'':->4s}-+-{'':->8s}-+-{'':->9s}-+-{'':->15s}")

    any_richer = False
    for sc in sorted(scores_to_info.keys()):
        info = scores_to_info[sc]
        n_H = len(info['H_vals'])
        n_bin = len(info['omega_types'])
        n_ter = len(info['weight_types'])
        richer = n_ter > n_bin
        if richer:
            any_richer = True
        print(f"  {str(sc):>25s} | {n_H:>4d} | {n_bin:>8d} | {n_ter:>9d} | {'YES!' if richer else 'no'}")

    if any_richer:
        print(f"\n  RESULT: Ternary carries MORE information than binary at n={n}!")
        print(f"  This means the overlap WEIGHT (not just adjacency) distinguishes")
        print(f"  tournaments that have the same conflict graph Omega.")
    else:
        print(f"\n  RESULT: Ternary and binary carry SAME information at n={n}.")
        print(f"  At this n, knowing which cycles conflict is enough.")


print("\n" + "=" * 70)
print("DONE.")
