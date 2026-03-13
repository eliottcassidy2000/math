#!/usr/bin/env python3
"""
cross_field_bridges.py — opus-2026-03-13-S67j

CREATIVE CROSS-FIELD CONNECTIONS for tournament spectral theory.

Explores 6 bridges:
1. ISING MODEL: H as spin glass partition function on {-1,+1}^m
2. BOOLEAN COMPLEXITY: Fourier influence, noise sensitivity, junta proximity
3. ARC REVERSAL LANDSCAPE: DPO rewriting → gradient flow on tournament hypercube
4. CODING THEORY: Regular tournaments as codewords, minimum distance
5. QUANTUM ANALOGY: H as quasi-classical observable, entanglement entropy
6. MARKOV MIXING: Arc reversal chain, spectral gap, H as eigenfunction

Each connection is investigated computationally at small n.
"""

import numpy as np
from itertools import combinations, permutations
import math

def all_tournaments(n):
    """Generate all tournaments on n vertices as adjacency matrices."""
    edges = [(i,j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)
    tournaments = []
    for bits in range(2**m):
        A = np.zeros((n,n), dtype=int)
        for k, (i,j) in enumerate(edges):
            if (bits >> k) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        tournaments.append((bits, A))
    return tournaments, edges

def count_hamiltonian_paths(A):
    """Count Hamiltonian paths in tournament A."""
    n = A.shape[0]
    count = 0
    for perm in permutations(range(n)):
        valid = True
        for i in range(n-1):
            if A[perm[i]][perm[i+1]] != 1:
                valid = False
                break
        if valid:
            count += 1
    return count

def sigma_encoding(bits, m):
    """Convert bits to sigma in {-1,+1}^m."""
    return np.array([1 if (bits >> k) & 1 else -1 for k in range(m)])

def walsh_hadamard_transform(f_values, m):
    """Full Walsh-Hadamard transform of function on {-1,+1}^m."""
    N = 2**m
    coeffs = {}
    for S_mask in range(N):
        val = 0
        for x in range(N):
            sigma = sigma_encoding(x, m)
            # chi_S(sigma) = prod_{i in S} sigma_i
            chi = 1
            for k in range(m):
                if (S_mask >> k) & 1:
                    chi *= sigma[k]
            val += f_values[x] * chi
        coeffs[S_mask] = val / N
    return coeffs

def score_sequence(A):
    return tuple(sorted(A.sum(axis=1).astype(int)))

# =====================================================================
# BRIDGE 1: ISING MODEL INTERPRETATION
# =====================================================================
print("=" * 70)
print("BRIDGE 1: ISING MODEL — H AS SPIN GLASS ON BOOLEAN HYPERCUBE")
print("=" * 70)

for n in [3, 4, 5]:
    tournaments, edges = all_tournaments(n)
    m = len(edges)

    # Compute H for all tournaments
    H_values = {}
    for bits, A in tournaments:
        H_values[bits] = count_hamiltonian_paths(A)

    # Walsh-Hadamard decomposition
    f_array = np.array([H_values[bits] for bits in range(2**m)])
    coeffs = walsh_hadamard_transform(f_array, m)

    print(f"\n  n={n}, m={m}:")

    # Ising interpretation: H(σ) = Σ_S ĥ(S) χ_S(σ)
    # Degree 0 = external field (constant)
    # Degree 2 = pairwise interactions J_{ij}
    # Degree 4 = 4-body interactions

    # Extract coupling constants
    J_matrix = np.zeros((m, m))  # pairwise couplings
    for S_mask in range(2**m):
        S_bits = [k for k in range(m) if (S_mask >> k) & 1]
        if len(S_bits) == 2:
            i, j = S_bits
            J_matrix[i][j] = coeffs[S_mask]
            J_matrix[j][i] = coeffs[S_mask]

    # Statistics of Ising couplings
    J_upper = [J_matrix[i][j] for i in range(m) for j in range(i+1,m)]
    if J_upper:
        J_upper = np.array(J_upper)
        J_nonzero = J_upper[np.abs(J_upper) > 1e-12]
        print(f"    Pairwise couplings: {len(J_nonzero)} non-zero out of {len(J_upper)}")
        if len(J_nonzero) > 0:
            print(f"    |J| range: [{min(abs(J_nonzero)):.6f}, {max(abs(J_nonzero)):.6f}]")
            # Frustration: fraction of triangles with odd number of negative couplings
            # (In Ising physics, frustrated systems have complex ground states)

    # Ground state = tournament maximizing H
    max_H = max(f_array)
    ground_states = [bits for bits in range(2**m) if H_values[bits] == max_H]
    print(f"    Ground state H = {max_H}, degeneracy = {len(ground_states)}")

    # Partition function Z(β) = Σ exp(β·H)
    # At β=0: Z = 2^m (infinite temperature)
    # As β→∞: Z → degeneracy * exp(β·max_H)
    for beta in [0.1, 0.5, 1.0, 2.0]:
        Z = sum(math.exp(beta * H_values[bits]) for bits in range(2**m))
        E_H = sum(H_values[bits] * math.exp(beta * H_values[bits]) for bits in range(2**m)) / Z
        print(f"    β={beta:.1f}: Z={Z:.2f}, <H>={E_H:.4f}, <H>/max={E_H/max_H:.4f}")

# =====================================================================
# BRIDGE 2: BOOLEAN FUNCTION COMPLEXITY
# =====================================================================
print("\n" + "=" * 70)
print("BRIDGE 2: BOOLEAN COMPLEXITY — INFLUENCE AND NOISE SENSITIVITY")
print("=" * 70)

for n in [3, 4, 5]:
    tournaments, edges = all_tournaments(n)
    m = len(edges)

    f_array = np.array([count_hamiltonian_paths(A) for _, A in tournaments])
    coeffs = walsh_hadamard_transform(f_array, m)

    print(f"\n  n={n}, m={m}:")

    # Total influence I[f] = Σ_i Inf_i(f) where Inf_i(f) = Σ_{S∋i} ĥ(S)²
    influences = np.zeros(m)
    for S_mask in range(2**m):
        S_bits = [k for k in range(m) if (S_mask >> k) & 1]
        for k in S_bits:
            influences[k] += coeffs[S_mask]**2

    total_influence = sum(influences)
    total_variance = sum(coeffs[S]**2 for S in range(1, 2**m))  # exclude constant

    print(f"    Total influence: I[H] = {total_influence:.6f}")
    print(f"    Total variance: Var[H] = {total_variance:.6f}")
    print(f"    I[H] / Var[H] = {total_influence/total_variance:.6f}" if total_variance > 0 else "")
    print(f"    Mean influence per edge: {total_influence/m:.6f}")
    print(f"    Max influence: {max(influences):.6f}")
    print(f"    Min influence: {min(influences):.6f}")
    print(f"    Influence range ratio: {max(influences)/min(influences):.6f}" if min(influences) > 1e-12 else "")

    # Noise sensitivity: NS_ρ(f) = Σ_{|S|>0} (1-ρ^2)^|S| ρ^{2(m-|S|)} ĥ(S)²
    # Simplified: NS_ε(f) = Σ_{S≠∅} (1-(1-2ε)^{2|S|}) ĥ(S)²  / (4*Var)
    # For boolean: fraction that changes when each bit flipped with prob ε
    for eps in [0.01, 0.05, 0.1]:
        rho = 1 - 2*eps
        ns = sum((1 - rho**(2*bin(S).count('1'))) * coeffs[S]**2
                 for S in range(1, 2**m))
        print(f"    NS_{eps:.2f} = {ns:.6f} (noise sensitivity)")

    # Degree distribution of energy
    degree_energy = {}
    for S_mask in range(2**m):
        deg = bin(S_mask).count('1')
        degree_energy[deg] = degree_energy.get(deg, 0) + coeffs[S_mask]**2
    print(f"    Fourier energy by degree:")
    for deg in sorted(degree_energy.keys()):
        if degree_energy[deg] > 1e-15:
            pct = 100 * degree_energy[deg] / sum(degree_energy.values())
            print(f"      deg {deg}: {degree_energy[deg]:.6f} ({pct:.2f}%)")

    # Friedgut-Kalai-Naor proximity to junta
    # A function with total influence I is ε-close to a (I/ε)-junta
    # For H: low influence suggests it's close to depending on few variables
    print(f"    FKN bound: H is ε-close to {total_influence:.1f}/ε-junta")

# =====================================================================
# BRIDGE 3: ARC REVERSAL LANDSCAPE (DPO Rewriting)
# =====================================================================
print("\n" + "=" * 70)
print("BRIDGE 3: ARC REVERSAL LANDSCAPE — DPO REWRITING ON TOURNAMENTS")
print("=" * 70)
print("  (Inspired by Bajaj 2024, arXiv:2409.01006 — hypergraph rewriting)")
print("  Each arc reversal = a DPO rewrite rule on the tournament digraph")

for n in [3, 4, 5]:
    tournaments, edges = all_tournaments(n)
    m = len(edges)

    H_values = {bits: count_hamiltonian_paths(A) for bits, A in tournaments}

    print(f"\n  n={n}, m={m}:")

    # Build flip graph: tournaments connected by single arc reversal
    # In σ-space, this is the hypercube graph on {-1,+1}^m
    # Gradient of H at each tournament: ΔH for each possible arc flip

    gradients = {}
    for bits in range(2**m):
        grad = np.zeros(m)
        for k in range(m):
            # Flip bit k
            flipped = bits ^ (1 << k)
            grad[k] = H_values[flipped] - H_values[bits]
        gradients[bits] = grad

    # Local maxima: all gradients ≤ 0
    local_max = [bits for bits in range(2**m)
                 if all(gradients[bits][k] <= 0 for k in range(m))]
    local_min = [bits for bits in range(2**m)
                 if all(gradients[bits][k] >= 0 for k in range(m))]

    global_max = max(H_values.values())
    global_min = min(H_values.values())

    print(f"    Local maxima: {len(local_max)} (H values: {sorted(set(H_values[b] for b in local_max))})")
    print(f"    Local minima: {len(local_min)} (H values: {sorted(set(H_values[b] for b in local_min))})")
    print(f"    Global max: {global_max}, Global min: {global_min}")

    # Gradient magnitude distribution
    grad_mags = [np.linalg.norm(gradients[bits]) for bits in range(2**m)]
    print(f"    |∇H| range: [{min(grad_mags):.4f}, {max(grad_mags):.4f}]")
    print(f"    Mean |∇H|: {np.mean(grad_mags):.4f}")

    # Connection to score sequence: gradient at regular tournaments
    for bits, A in tournaments:
        scores = tuple(sorted(A.sum(axis=1).astype(int)))
        if n % 2 == 1 and scores == tuple([(n-1)//2]*n):
            # Regular tournament
            grad = gradients[bits]
            print(f"    Regular T (bits={bits}): H={H_values[bits]}, |∇H|={np.linalg.norm(grad):.4f}")
            # ΔH for each arc is related to the Z_v balance at endpoints
            break

    # Steepest ascent from transitive tournament
    # Follow greedy path to local maximum
    trans_bits = None
    for bits, A in tournaments:
        if score_sequence(A) == tuple(range(n)):
            trans_bits = bits
            break

    if trans_bits is not None:
        path = [trans_bits]
        current = trans_bits
        while True:
            grad = gradients[current]
            best_k = np.argmax(grad)
            if grad[best_k] <= 0:
                break
            current = current ^ (1 << best_k)
            path.append(current)

        H_path = [H_values[b] for b in path]
        print(f"    Steepest ascent from transitive: {len(path)} steps")
        print(f"    H along path: {H_path}")
        if path:
            final_scores = score_sequence(dict(tournaments)[path[-1]])
            print(f"    Arrives at score seq: {final_scores}, H={H_values[path[-1]]}")

    # Basin of attraction sizes (greedy ascent from each tournament)
    basins = {}
    for bits in range(2**m):
        current = bits
        while True:
            grad = gradients[current]
            best_k = np.argmax(grad)
            if grad[best_k] <= 0:
                break
            current = current ^ (1 << best_k)
        basins[current] = basins.get(current, 0) + 1

    print(f"    Basin structure ({len(basins)} attractors):")
    for attractor, size in sorted(basins.items(), key=lambda x: -x[1])[:5]:
        print(f"      H={H_values[attractor]}, basin={size}, "
              f"scores={score_sequence(dict(tournaments)[attractor])}")

# =====================================================================
# BRIDGE 4: CODING THEORY — REGULAR TOURNAMENTS AS CODEWORDS
# =====================================================================
print("\n" + "=" * 70)
print("BRIDGE 4: CODING THEORY — TOURNAMENT CODEWORDS")
print("=" * 70)

for n in [3, 5]:
    tournaments, edges = all_tournaments(n)
    m = len(edges)

    H_values = {bits: count_hamiltonian_paths(A) for bits, A in tournaments}
    max_H = max(H_values.values())

    # "Codewords" = tournaments achieving maximum H
    codewords = [bits for bits in range(2**m) if H_values[bits] == max_H]

    print(f"\n  n={n}, m={m}, max H = {max_H}:")
    print(f"    Codewords (H=max): {len(codewords)}")

    # Minimum Hamming distance between codewords
    if len(codewords) > 1:
        min_dist = m
        for i in range(len(codewords)):
            for j in range(i+1, len(codewords)):
                d = bin(codewords[i] ^ codewords[j]).count('1')
                min_dist = min(min_dist, d)
        print(f"    Minimum Hamming distance: {min_dist}")
        print(f"    Rate: R = log2({len(codewords)})/{m} = {math.log2(len(codewords))/m:.4f}")

        # Distance distribution
        distances = {}
        for i in range(len(codewords)):
            for j in range(i+1, len(codewords)):
                d = bin(codewords[i] ^ codewords[j]).count('1')
                distances[d] = distances.get(d, 0) + 1
        print(f"    Distance distribution: {dict(sorted(distances.items()))}")

    # Weight distribution (number of 1-bits in each codeword)
    weights = [bin(c).count('1') for c in codewords]
    print(f"    Weight distribution: {sorted(set(weights))}")

    # Dual distance (smallest non-zero Fourier coefficient of indicator)
    indicator = np.zeros(2**m)
    for c in codewords:
        indicator[c] = 1
    ind_coeffs = walsh_hadamard_transform(indicator, m)
    nonzero_degrees = set()
    for S_mask in range(1, 2**m):
        if abs(ind_coeffs[S_mask]) > 1e-12:
            nonzero_degrees.add(bin(S_mask).count('1'))
    if nonzero_degrees:
        print(f"    Dual distance (min nonzero Fourier deg): {min(nonzero_degrees)}")

# =====================================================================
# BRIDGE 5: MARKOV CHAIN MIXING
# =====================================================================
print("\n" + "=" * 70)
print("BRIDGE 5: MARKOV CHAIN — ARC REVERSAL MIXING TIME")
print("=" * 70)

for n in [3, 4, 5]:
    tournaments, edges = all_tournaments(n)
    m = len(edges)
    N = 2**m

    H_values = {bits: count_hamiltonian_paths(A) for bits, A in tournaments}

    print(f"\n  n={n}, m={m}, {N} tournaments:")

    # Uniform arc reversal chain: at each step, pick random edge, flip
    # This is a random walk on hypercube {0,1}^m
    # Transition matrix: P(x,y) = 1/m if Hamming distance 1, else 0
    # (lazy version: P(x,x) = 1/2, P(x,y) = 1/(2m) if dist 1)

    # Eigenvalues of hypercube random walk are known: 1-2k/m for k=0,...,m
    # with multiplicity C(m,k)

    # H-weighted chain (Metropolis): accept flip if H increases or with prob exp(-β·ΔH)
    # This has stationary distribution proportional to exp(β·H)

    # Instead of building the full transition matrix, use the known structure:
    # The hypercube adjacency eigenvalues are m-2k, k=0,...,m
    # So the lazy walk eigenvalues are 1-k/m

    print(f"    Hypercube walk spectral gap: 1/m = {1/m:.4f}")
    print(f"    Mixing time (relaxation): m·log(2^m) = {m * m * math.log(2):.1f}")

    # H as an eigenfunction: since H is a polynomial in σ of degree ≤ 2floor((n-1)/2),
    # it lives in the eigenspaces with eigenvalue 1-k/m for k ≤ degree

    # More interesting: Metropolis chain that samples tournaments proportional to exp(β·H)
    # Build it for n=3,4 where feasible
    if N <= 256:
        for beta in [0.5, 1.0, 2.0]:
            P = np.zeros((N, N))
            for x in range(N):
                for k in range(m):
                    y = x ^ (1 << k)
                    delta_H = H_values[y] - H_values[x]
                    accept = min(1, math.exp(beta * delta_H))
                    P[x][y] = accept / m
                P[x][x] = 1 - sum(P[x][y] for y in range(N) if y != x)

            # Find stationary distribution
            eigenvalues = np.linalg.eigvals(P.T)
            eigenvalues_sorted = sorted(abs(eigenvalues), reverse=True)
            spectral_gap = 1 - eigenvalues_sorted[1]

            # Stationary distribution should be proportional to exp(β·H)
            Z_beta = sum(math.exp(beta * H_values[b]) for b in range(N))
            pi_expected = np.array([math.exp(beta * H_values[b]) / Z_beta for b in range(N)])

            # Check mixing by power iteration
            pi_uniform = np.ones(N) / N
            pi_t = pi_uniform.copy()
            for _ in range(100):
                pi_t = pi_t @ P

            tv_distance = 0.5 * np.sum(np.abs(pi_t - pi_expected))

            print(f"    Metropolis β={beta:.1f}: gap={spectral_gap:.4f}, "
                  f"τ_mix≈{1/spectral_gap:.1f}, TV(100 steps)={tv_distance:.6f}")

# =====================================================================
# BRIDGE 6: ENTROPY AND INFORMATION
# =====================================================================
print("\n" + "=" * 70)
print("BRIDGE 6: INFORMATION THEORY — TOURNAMENT ENTROPY")
print("=" * 70)

for n in [3, 4, 5]:
    tournaments, edges = all_tournaments(n)
    m = len(edges)

    H_values = {bits: count_hamiltonian_paths(A) for bits, A in tournaments}
    total_H = sum(H_values.values())

    print(f"\n  n={n}, m={m}:")

    # Treat H/total as probability distribution on tournaments
    probs = np.array([H_values[b] / total_H for b in range(2**m)])
    entropy = -sum(p * math.log2(p) for p in probs if p > 0)
    max_entropy = math.log2(2**m)

    print(f"    Shannon entropy of H-distribution: {entropy:.4f} bits")
    print(f"    Maximum entropy: {max_entropy:.4f} bits")
    print(f"    Efficiency: {entropy/max_entropy:.4f}")

    # KL divergence from uniform
    kl = sum(p * math.log2(p * 2**m) for p in probs if p > 0)
    print(f"    KL(H-dist || uniform) = {kl:.4f} bits")

    # Mutual information between H and score sequence
    # I(H; score_seq) = H(H_values) - H(H_values | score_seq)
    score_to_H = {}
    for bits, A in tournaments:
        ss = score_sequence(A)
        if ss not in score_to_H:
            score_to_H[ss] = []
        score_to_H[ss].append(H_values[bits])

    # H(H_values | score_seq)
    cond_entropy = 0
    for ss, h_list in score_to_H.items():
        p_ss = len(h_list) / (2**m)
        h_counts = {}
        for h in h_list:
            h_counts[h] = h_counts.get(h, 0) + 1
        for count in h_counts.values():
            p_h_given_ss = count / len(h_list)
            cond_entropy -= p_ss * p_h_given_ss * math.log2(p_h_given_ss) if p_h_given_ss > 0 else 0

    # H(H_values)
    h_counts_total = {}
    for h in H_values.values():
        h_counts_total[h] = h_counts_total.get(h, 0) + 1
    h_entropy = -sum((c / 2**m) * math.log2(c / 2**m) for c in h_counts_total.values())

    mi = h_entropy - cond_entropy
    print(f"    H(H-values) = {h_entropy:.4f} bits")
    print(f"    H(H | score_seq) = {cond_entropy:.4f} bits")
    print(f"    I(H; score_seq) = {mi:.4f} bits")
    print(f"    H explained by scores: {mi/h_entropy*100:.1f}%" if h_entropy > 0 else "")

print("\n\nDONE — cross_field_bridges.py complete")
