#!/usr/bin/env python3
"""
information_theory_deep.py — opus-2026-03-13-S67j

DEEP INFORMATION-THEORETIC ANALYSIS OF TOURNAMENT H(T)

8 connections explored:

1. RANKING ENTROPY: H(T)/n! as a probability, S(T) = log(H(T)) as entropy
2. CHANNEL CAPACITY: Tournament as noisy pairwise comparison channel
3. MINIMUM DESCRIPTION LENGTH: Bits needed to specify H(T)
4. CAUSAL ENTROPY: Entropy production in greedy ascent trajectories
5. DIRECTED INFORMATION: Tournament as a causal DAG, flow of information
6. SOURCE CODING: Tournament source at inverse temperature β
7. KOLMOGOROV COMPLEXITY vs H: Algorithmic simplicity vs ranking ambiguity
8. FISHER INFORMATION: Score sequence as sufficient statistic

Inspired by arXiv:2409.01006 (hypergraph rewriting & causality).
"""

import numpy as np
from itertools import permutations
import math
import random

def all_tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)
    ts = []
    for bits in range(2**m):
        A = np.zeros((n,n), dtype=int)
        for k, (i,j) in enumerate(edges):
            if (bits >> k) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        ts.append((bits, A))
    return ts, edges

def count_ham_paths(A):
    n = A.shape[0]
    count = 0
    for perm in permutations(range(n)):
        valid = all(A[perm[i]][perm[i+1]] == 1 for i in range(n-1))
        if valid:
            count += 1
    return count

def score_seq(A):
    return tuple(sorted(A.sum(axis=1).astype(int)))

# =====================================================================
# 1. RANKING ENTROPY
# =====================================================================
print("=" * 70)
print("1. RANKING ENTROPY: S(T) = log₂(H(T))")
print("=" * 70)
print("  H(T)/n! = fraction of total orderings consistent with T")
print("  S(T) = log₂(H(T)) = 'ranking entropy' in bits")
print("  Maximum: S_max = log₂(n!) for transitive (all orderings work) — WRONG")
print("  Actually: H(transitive) = 1, so S_min = 0")
print("  H(regular) = max → S_max = log₂(H_max)")

for n in [3, 4, 5]:
    ts, edges = all_tournaments(n)
    m = len(edges)

    H_vals = {bits: count_ham_paths(A) for bits, A in ts}
    max_H = max(H_vals.values())
    min_H = min(H_vals.values())

    print(f"\n  n={n}:")
    print(f"    S_min = log₂({min_H}) = {math.log2(min_H):.4f} bits")
    print(f"    S_max = log₂({max_H}) = {math.log2(max_H):.4f} bits")
    print(f"    S_mean = log₂(E[H]) = log₂({math.factorial(n)/2**(n-1):.2f}) = "
          f"{math.log2(math.factorial(n)/2**(n-1)):.4f} bits")
    print(f"    S_max/log₂(n!) = {math.log2(max_H)/math.log2(math.factorial(n)):.4f}")

    # Entropy of H-distribution (treating H values as symbols)
    h_counts = {}
    for h in H_vals.values():
        h_counts[h] = h_counts.get(h, 0) + 1
    p_vals = [c / len(H_vals) for c in h_counts.values()]
    H_entropy = -sum(p * math.log2(p) for p in p_vals if p > 0)
    print(f"    Shannon entropy of H-value distribution: {H_entropy:.4f} bits")
    print(f"    Max possible (log₂ of distinct H values): {math.log2(len(h_counts)):.4f}")

# =====================================================================
# 2. TOURNAMENT AS NOISY CHANNEL
# =====================================================================
print("\n" + "=" * 70)
print("2. NOISY CHANNEL: TOURNAMENT AS PAIRWISE COMPARISON CHANNEL")
print("=" * 70)
print("  Model: n items have a true ranking π ∈ S_n")
print("  Channel: each pair (i,j) compared independently")
print("    Correct comparison w.p. p > 1/2, error w.p. 1-p")
print("  Output: tournament T")
print("  H(T) = # rankings consistent with T")
print("  Channel capacity: how many bits of ranking info per comparison?")

n = 5
m = n*(n-1)//2
ts, edges = all_tournaments(n)
H_vals = {bits: count_ham_paths(A) for bits, A in ts}

# For each accuracy level p, compute expected H
for p_acc in [0.51, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0]:
    # Expected H when true ranking is observed through noisy channel
    # For a true ranking π, each edge is correct w.p. p_acc
    # The expected number of consistent rankings = ?
    # This is related to the Mallows model

    # For the identity ranking, P(T | π=id) = p_acc^(agreements) * (1-p_acc)^(disagreements)
    # where agreements = # edges consistent with id

    # Expected H = E_T[H(T)] under this noise model
    # = sum_T P(T|π) * H(T)

    # For identity ranking: P(T|π=id) depends on # of "forward" edges
    expected_H = 0
    for bits, A in ts:
        # Count forward edges (consistent with natural ordering 0<1<...<n-1)
        forward = sum(1 for k, (i,j) in enumerate(edges)
                     if ((bits >> k) & 1) == 1)  # i->j matches i<j
        p_T = p_acc**forward * (1-p_acc)**(m-forward)
        expected_H += p_T * H_vals[bits]

    # Mutual information: I(π; T) = H(π) - H(π|T)
    # H(π) = log₂(n!)
    # H(π|T) = E_T[log₂(H(T))]  (since each consistent ranking equally likely)
    H_pi = math.log2(math.factorial(n))

    expected_log_H = 0
    for bits, A in ts:
        forward = sum(1 for k, (i,j) in enumerate(edges)
                     if ((bits >> k) & 1) == 1)
        p_T = p_acc**forward * (1-p_acc)**(m-forward)
        if p_T > 0 and H_vals[bits] > 0:
            expected_log_H += p_T * math.log2(H_vals[bits])

    MI = H_pi - expected_log_H
    efficiency = MI / m  # bits per comparison

    print(f"  p={p_acc:.2f}: E[H]={expected_H:.2f}, "
          f"I(π;T)={MI:.4f} bits, efficiency={efficiency:.4f} bits/comparison")

print(f"\n  Maximum I(π;T) = log₂({n}!) = {math.log2(math.factorial(n)):.4f} bits")
print(f"  With {m} comparisons → max {math.log2(math.factorial(n))/m:.4f} bits/comparison")

# =====================================================================
# 3. MINIMUM DESCRIPTION LENGTH
# =====================================================================
print("\n" + "=" * 70)
print("3. MINIMUM DESCRIPTION LENGTH FOR H(T)")
print("=" * 70)

for n in [3, 4, 5]:
    ts, edges = all_tournaments(n)
    m = len(edges)

    H_vals = {bits: count_ham_paths(A) for bits, A in ts}

    # Full description: m bits
    full_bits = m

    # Score description: n integers, each in {0,...,n-1}
    # Actually: sorted score sequence is a partition of m into n parts
    score_seqs = set()
    for bits, A in ts:
        score_seqs.add(score_seq(A))
    score_bits = math.log2(len(score_seqs)) if score_seqs else 0

    # How much of H is determined by scores?
    score_to_H = {}
    for bits, A in ts:
        ss = score_seq(A)
        if ss not in score_to_H:
            score_to_H[ss] = set()
        score_to_H[ss].add(H_vals[bits])

    # Conditional entropy H(H_value | score)
    cond_ent = 0
    for ss, h_set in score_to_H.items():
        # Count tournaments with this score
        count_ss = sum(1 for b, A in ts if score_seq(A) == ss)
        # H-value distribution within this score class
        h_counts = {}
        for b, A in ts:
            if score_seq(A) == ss:
                h = H_vals[b]
                h_counts[h] = h_counts.get(h, 0) + 1
        for h, c in h_counts.items():
            p_h = c / count_ss
            cond_ent -= (count_ss / len(ts)) * p_h * math.log2(p_h) if p_h > 0 else 0

    # MDL: score_bits + cond_ent
    h_distinct = len(set(H_vals.values()))
    h_entropy = math.log2(h_distinct)

    print(f"\n  n={n}, m={m}:")
    print(f"    Full tournament: {full_bits} bits")
    print(f"    Score sequences: {len(score_seqs)} ({score_bits:.2f} bits)")
    print(f"    H distinct values: {h_distinct} ({h_entropy:.2f} bits)")
    print(f"    H(H_value | score): {cond_ent:.4f} bits")
    print(f"    MDL for H: score ({score_bits:.2f}) + residual ({cond_ent:.4f}) "
          f"= {score_bits + cond_ent:.4f} bits")
    print(f"    Compression ratio: {(score_bits + cond_ent)/full_bits:.4f}")
    print(f"    Score sufficiency: {1 - cond_ent/h_entropy:.4f}" if h_entropy > 0 else "")

# =====================================================================
# 4. CAUSAL ENTROPY OF GREEDY ASCENT
# =====================================================================
print("\n" + "=" * 70)
print("4. CAUSAL ENTROPY: UNCERTAINTY IN OPTIMIZATION TRAJECTORIES")
print("=" * 70)
print("  Inspired by Bajaj 2024 (arXiv:2409.01006): events have causal structure")
print("  In our setting: each arc reversal is an 'event'")
print("  Greedy ascent = sequence of events with causal dependencies")
print("  Causal entropy = log(# distinct greedy paths to optimum)")

for n in [4, 5]:
    ts, edges = all_tournaments(n)
    m = len(edges)

    H_vals = {bits: count_ham_paths(A) for bits, A in ts}
    max_H = max(H_vals.values())

    # For each tournament, count distinct steepest ascent paths
    # (there may be ties in gradient, creating branching)
    path_counts = {}
    path_lengths = {}

    for start_bits in range(2**m):
        # BFS/DFS to enumerate all greedy paths
        def count_paths(bits):
            h_current = H_vals[bits]
            # Find all maximally improving neighbors
            max_delta = 0
            for k in range(m):
                neighbor = bits ^ (1 << k)
                delta = H_vals[neighbor] - h_current
                if delta > max_delta:
                    max_delta = delta

            if max_delta <= 0:
                return 1, 0  # local max, 1 path of length 0

            # All neighbors achieving max_delta
            best_neighbors = []
            for k in range(m):
                neighbor = bits ^ (1 << k)
                if H_vals[neighbor] - h_current == max_delta:
                    best_neighbors.append(neighbor)

            total_paths = 0
            max_len = 0
            for nb in best_neighbors:
                p, l = count_paths(nb)
                total_paths += p
                max_len = max(max_len, l + 1)
            return total_paths, max_len

        n_paths, path_len = count_paths(start_bits)
        final_H = H_vals[start_bits]
        # Follow one path to find final H
        current = start_bits
        while True:
            h_current = H_vals[current]
            max_delta = 0
            best_k = -1
            for k in range(m):
                delta = H_vals[current ^ (1 << k)] - h_current
                if delta > max_delta:
                    max_delta = delta
                    best_k = k
            if best_k == -1:
                break
            current = current ^ (1 << best_k)
        final_H = H_vals[current]

        if final_H not in path_counts:
            path_counts[final_H] = []
        path_counts[final_H].append(n_paths)
        path_lengths[final_H] = path_lengths.get(final_H, [])
        path_lengths[final_H].append(path_len)

    print(f"\n  n={n}, m={m}:")
    for final_h in sorted(path_counts.keys()):
        counts = path_counts[final_h]
        lengths = path_lengths[final_h]
        avg_paths = np.mean(counts)
        max_paths = max(counts)
        avg_len = np.mean(lengths)
        causal_entropy = np.mean([math.log2(c) if c > 0 else 0 for c in counts])
        print(f"    Final H={final_h}: {len(counts)} starts, "
              f"avg_paths={avg_paths:.1f}, max_paths={max_paths}, "
              f"causal_entropy={causal_entropy:.3f} bits, "
              f"avg_steps={avg_len:.1f}")

# =====================================================================
# 5. DIRECTED INFORMATION FLOW
# =====================================================================
print("\n" + "=" * 70)
print("5. DIRECTED INFORMATION FLOW IN TOURNAMENTS")
print("=" * 70)
print("  Each arc i→j represents 'i beats j' = information flows from i to j")
print("  The REACHABILITY MATRIX R[i][j] = 1 iff there's a directed path from i to j")
print("  For tournaments: R[i][j] = 1 for all i≠j (Redei's theorem: HP exists)")
print("  So tournaments are 'maximally connected' information networks")
print("")
print("  More refined: TRANSFER ENTROPY T(i→j) = mutual info beyond j's past")
print("  In tournaments: vertex j's 'past' = {k : k→j}, 'future' = {k : j→k}")
print("  Transfer entropy captures indirect influence beyond pairwise comparison")

for n in [5]:
    ts, edges = all_tournaments(n)
    m = len(edges)

    # For each tournament, compute vertex "influence"
    # Influence of v = how much flipping v's outgoing arcs changes H
    print(f"\n  n={n}: vertex influence analysis")

    for bits, A in ts[:5]:  # just a few examples
        h0 = count_ham_paths(A)
        scores = A.sum(axis=1).astype(int)

        # Influence: flip ALL outgoing arcs of vertex v
        influences = []
        for v in range(n):
            A_flipped = A.copy()
            for j in range(n):
                if j != v:
                    A_flipped[v][j] = 1 - A[v][j]
                    A_flipped[j][v] = 1 - A[j][v]
            h_flipped = count_ham_paths(A_flipped)
            influences.append(abs(h_flipped - h0))

        # Correlation between score and influence
        corr = np.corrcoef(scores, influences)[0,1] if np.std(scores) > 0 else 0
        print(f"    T(bits={bits}): H={h0}, scores={tuple(scores)}, "
              f"influences={influences}, corr(s,inf)={corr:.3f}")

# =====================================================================
# 6. SOURCE CODING AT TEMPERATURE β
# =====================================================================
print("\n" + "=" * 70)
print("6. SOURCE CODING: ENTROPY vs TEMPERATURE")
print("=" * 70)

n = 5
ts, edges = all_tournaments(n)
m = len(edges)
H_vals = np.array([count_ham_paths(A) for _, A in ts])

print(f"  n={n}: Tournament source at temperature T = 1/β")
print(f"  P(T) ∝ exp(β·H(T))")
print(f"  Source entropy S(β) = -Σ P(T) log₂ P(T)")

for beta in [0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0]:
    log_weights = beta * H_vals
    log_weights -= np.max(log_weights)  # numerical stability
    weights = np.exp(log_weights)
    Z = np.sum(weights)
    probs = weights / Z

    entropy = -np.sum(probs * np.log2(probs + 1e-300))
    E_H = np.sum(probs * H_vals)
    max_entropy = math.log2(2**m)

    # Number of "effective" configurations
    n_eff = 2**entropy

    print(f"  β={beta:7.3f}: S={entropy:6.2f} bits, <H>={E_H:6.2f}, "
          f"n_eff={n_eff:8.1f}, S/m={entropy/m:.4f}")

# Phase transition: where does entropy drop fastest?
print(f"\n  Phase transition analysis:")
betas = np.linspace(0.01, 5.0, 200)
prev_S = None
max_dS = 0
critical_beta = 0
for beta in betas:
    log_w = beta * H_vals - np.max(beta * H_vals)
    w = np.exp(log_w)
    Z = np.sum(w)
    p = w / Z
    S = -np.sum(p * np.log2(p + 1e-300))
    if prev_S is not None:
        dS = abs(S - prev_S) / (betas[1] - betas[0])
        if dS > max_dS:
            max_dS = dS
            critical_beta = beta
    prev_S = S

print(f"    Steepest entropy drop at β_c ≈ {critical_beta:.3f} (T_c ≈ {1/critical_beta:.3f})")

# =====================================================================
# 7. KOLMOGOROV COMPLEXITY vs H
# =====================================================================
print("\n" + "=" * 70)
print("7. KOLMOGOROV COMPLEXITY vs RANKING AMBIGUITY")
print("=" * 70)
print("  Conjecture: K(T) and H(T) are ANTI-correlated")
print("  Simple tournaments (low K) tend to have extreme H")
print("  Complex tournaments (high K) tend to have middling H")

n = 5
ts, edges = all_tournaments(n)
m = len(edges)

# Approximate K(T) by compression ratio
# Compute compressibility: length of run-length encoding of bit string
for bits, A in ts[:10]:
    h = count_ham_paths(A)
    bit_str = format(bits, f'0{m}b')

    # Run-length encoding size
    runs = 1
    for i in range(1, len(bit_str)):
        if bit_str[i] != bit_str[i-1]:
            runs += 1
    rle_size = runs * (1 + math.ceil(math.log2(m + 1)))

    # Symmetry: is T isomorphic to T^op?
    # T^op has all arcs reversed
    bits_op = 0
    for k in range(m):
        if not ((bits >> k) & 1):
            bits_op |= (1 << k)

    # Is T = T^op (self-complementary)?
    is_sc = (bits == bits_op)

    ss = score_seq(A)
    score_var = np.var([s for s in ss])

    print(f"  bits={bits:4d}: H={h:3d}, runs={runs}, RLE={rle_size}, "
          f"SC={is_sc}, scores={ss}, Var={score_var:.2f}")

# =====================================================================
# 8. FISHER INFORMATION OF SCORE STATISTIC
# =====================================================================
print("\n" + "=" * 70)
print("8. FISHER INFORMATION: SCORE AS SUFFICIENT STATISTIC")
print("=" * 70)
print("  Score sequence s = (s_1,...,s_n) is a statistic of T")
print("  By THM-163: H_2 = c_2 · n · (m - 2·Var(s))")
print("  Score is 'almost sufficient' for H (85-100% info captured)")
print("")
print("  Fisher information I_F(θ) = Var(∂logP/∂θ) measures info content")
print("  For exponential family P(T) ∝ exp(β·H(T)):")
print("    I_F(β) = Var_β[H(T)] = <H²> - <H>²")

n = 5
ts, edges = all_tournaments(n)
m = len(edges)
H_vals = np.array([count_ham_paths(A) for _, A in ts])

print(f"\n  n={n}: Fisher information I_F(β) = Var_β[H]")
for beta in [0.01, 0.1, 0.5, 1.0, 2.0, 5.0]:
    log_w = beta * H_vals - np.max(beta * H_vals)
    w = np.exp(log_w)
    Z = np.sum(w)
    p = w / Z

    E_H = np.sum(p * H_vals)
    E_H2 = np.sum(p * H_vals**2)
    var_H = E_H2 - E_H**2
    fisher = beta**2 * var_H  # = specific heat C_v

    # Cramer-Rao bound: Var(θ_hat) ≥ 1/I_F
    # Minimum number of samples needed to estimate β to precision ε:
    # n_samples ≥ 1/(ε² · I_F)

    print(f"  β={beta:.2f}: <H>={E_H:.2f}, Var[H]={var_H:.4f}, "
          f"I_F={fisher:.4f}")

    # Score-based Fisher information (using only H_2 ≈ f(score_var))
    score_vars = np.array([np.var(A.sum(axis=1)) for _, A in ts])
    E_sv = np.sum(p * score_vars)
    E_sv2 = np.sum(p * score_vars**2)
    var_sv = E_sv2 - E_sv**2
    print(f"        Var[score_var]={var_sv:.4f}, "
          f"Fisher ratio: {beta**2*var_sv/(fisher+1e-10):.4f} (=1 if score sufficient)")

print("\n\nDONE — information_theory_deep.py complete")
