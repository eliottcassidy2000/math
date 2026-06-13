#!/usr/bin/env python3
"""
tournament_probe_design.py -- kind-pasteur-2026-03-21-S12

Design document / proof of concept for TournamentProbe:
A parameter-free LLM analyzer using tournament theory.

PRODUCT CONCEPT:
Instead of training behavioral probes (Napolitano's approach),
extract tournament structure from attention patterns and compute
invariants that are MATHEMATICALLY GROUNDED in our proved theorems.

KEY ADVANTAGE: No training required. All invariants are functorial
(respect structure-preserving maps), so they generalize across
models, tasks, and contexts without retraining.

INVARIANT SUITE:
1. H(T_attn) — Hamiltonian path count (OCF gives this = I(Omega,2))
2. Score sequence — degree distribution of attention dominance
3. alpha_k — independence number profile (from OCF decomposition)
4. beta_k — path homology Betti numbers
5. Cartan energy — anti/sym/scalar fractions
6. Transitivity index — fraction of transitive triples
7. Regularity — how close to doubly-regular (Paley-like)

IMPLEMENTATION SKETCH (for HuggingFace models):
  model = AutoModel.from_pretrained("...")
  outputs = model(**inputs, output_attentions=True)
  for layer, attn in enumerate(outputs.attentions):
      for head in range(attn.shape[1]):
          A = attn[0, head].numpy()  # (seq_len, seq_len)
          T = attention_to_tournament(A)
          report = tournament_analysis(T)

Author: kind-pasteur-2026-03-21-S12
"""

import numpy as np
from collections import defaultdict
from itertools import combinations, permutations


# ============================================================
# CORE: Soft tournament (differentiable for training)
# ============================================================

def soft_tournament(A, temperature=0.1):
    """
    Differentiable approximation to tournament extraction.

    Hard version: T[i,j] = 1 if A[i,j] > A[j,i]
    Soft version: T[i,j] = sigmoid((A[i,j] - A[j,i]) / tau)

    As tau -> 0, approaches hard tournament.
    For tau > 0, differentiable w.r.t. A.

    Property: T[i,j] + T[j,i] = 1 for all i,j (soft tournament axiom).
    """
    n = A.shape[0]
    diff = A - A.T  # antisymmetric
    T_soft = 1 / (1 + np.exp(-diff / temperature))
    # Set diagonal to 0
    np.fill_diagonal(T_soft, 0)
    return T_soft


def soft_hamiltonian_path_count(T_soft):
    """
    Approximate H(T) for soft tournament using matrix permanent relaxation.

    For hard {0,1} tournaments: H(T) = sum over permutations sigma of
    product T[sigma(1)][sigma(2)] * T[sigma(2)][sigma(3)] * ...

    This is related to but NOT equal to the permanent of T.
    It's the "path permanent" — sum over directed Hamiltonian path labelings.

    For soft [0,1] tournaments: the same formula gives a continuous relaxation.
    """
    n = T_soft.shape[0]
    if n > 12:
        print(f"  Warning: exact soft H for n={n} requires O(2^n * n^2) time")
        return None

    # DP: dp[mask][v] = sum of products over all paths ending at v
    #     visiting exactly the vertices in mask
    dp = np.zeros((1 << n, n))
    for v in range(n):
        dp[1 << v][v] = 1.0

    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] < 1e-15:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                dp[mask | (1 << w)][w] += dp[mask][v] * T_soft[v][w]

    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))


# ============================================================
# TOURNAMENT INVARIANT PROFILES
# ============================================================

def tournament_profile(T, name=""):
    """
    Compute a comprehensive profile of tournament T.
    Returns dict of invariants suitable for comparison.
    """
    n = T.shape[0]
    profile = {'name': name, 'n': n}

    # Score sequence
    scores = sorted(T.sum(axis=1))
    profile['scores'] = list(scores)
    profile['score_variance'] = np.var(scores)
    profile['is_regular'] = all(s == (n-1)/2 for s in scores)

    # Transitivity
    trans = 0
    total = 0
    for triple in combinations(range(n), 3):
        total += 1
        a, b, c = triple
        s = [T[a][b]+T[a][c], T[b][a]+T[b][c], T[c][a]+T[c][b]]
        if 0 in s:
            trans += 1
    profile['transitivity'] = trans / total if total > 0 else 1.0

    # H(T) (exact, if small enough)
    if n <= 14:
        dp = defaultdict(int)
        for v in range(n):
            dp[(1 << v, v)] = 1
        for mask in range(1, 1 << n):
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                if dp[(mask, v)] == 0:
                    continue
                for w in range(n):
                    if mask & (1 << w):
                        continue
                    if T[v][w]:
                        dp[(mask | (1 << w), w)] += dp[(mask, v)]
        full = (1 << n) - 1
        H = sum(dp[(full, v)] for v in range(n))
        profile['H'] = H

        # H relative to maximum (Paley)
        # Known maxima: n=3:3, n=5:15, n=7:189
        known_max = {3: 3, 5: 15, 7: 189}
        if n in known_max:
            profile['H_fraction'] = H / known_max[n]

    # Cycle counts (if small enough)
    if n <= 9:
        c3 = 0
        for triple in combinations(range(n), 3):
            a, b, c = triple
            # Check for directed 3-cycle
            if T[a][b] and T[b][c] and T[c][a]:
                c3 += 1
            elif T[a][c] and T[c][b] and T[b][a]:
                c3 += 1
        profile['c3'] = c3

    # Cartan energy
    A_float = T.astype(float)
    A_anti = (A_float - A_float.T) / 2
    A_sym = (A_float + A_float.T) / 2
    norm2 = np.sum(A_float**2)
    if norm2 > 0:
        profile['anti_energy'] = np.sum(A_anti**2) / norm2
        profile['sym_energy'] = np.sum(A_sym**2) / norm2
    else:
        profile['anti_energy'] = 0
        profile['sym_energy'] = 0

    return profile


def demonstrate_soft_tournament_bridge():
    """
    Show the bridge between hard tournaments (our theory)
    and soft attention (transformer reality).

    Key finding: as temperature -> 0, soft H -> hard H.
    At intermediate temperatures, soft H interpolates smoothly.
    This makes tournament invariants differentiable probes.
    """
    print("SOFT TOURNAMENT BRIDGE")
    print("=" * 60)

    # Create a Paley T_7 attention matrix
    p = 7
    qr = set()
    for k in range(1, p):
        qr.add((k*k) % p)

    A = np.zeros((p, p))
    for i in range(p):
        for j in range(p):
            if i == j:
                continue
            if (j - i) % p in qr:
                A[i][j] = 1.0
    row_sums = A.sum(axis=1, keepdims=True)
    A = A / np.maximum(row_sums, 1e-10)

    print("\nPaley T_7 attention matrix (row-normalized):")

    # Hard tournament
    T_hard = np.zeros((p, p), dtype=int)
    for i in range(p):
        for j in range(i+1, p):
            if A[i, j] > A[j, i]:
                T_hard[i][j] = 1
            elif A[j, i] > A[i, j]:
                T_hard[j][i] = 1
            else:
                T_hard[i][j] = 1

    dp = defaultdict(int)
    for v in range(p):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << p):
        for v in range(p):
            if not (mask & (1 << v)):
                continue
            if dp[(mask, v)] == 0:
                continue
            for w in range(p):
                if mask & (1 << w):
                    continue
                if T_hard[v][w]:
                    dp[(mask | (1 << w), w)] += dp[(mask, v)]
    full = (1 << p) - 1
    H_hard = sum(dp[(full, v)] for v in range(p))
    print(f"  Hard H(T_7) = {H_hard}")

    # Soft tournament at various temperatures
    print("\nSoft H at varying temperature:")
    for tau in [10.0, 1.0, 0.5, 0.1, 0.05, 0.01, 0.001]:
        T_soft = soft_tournament(A, temperature=tau)
        H_soft = soft_hamiltonian_path_count(T_soft)
        if H_soft is not None:
            print(f"  tau={tau:6.3f}: soft H = {H_soft:.2f} "
                  f"(ratio to hard: {H_soft/H_hard:.4f})")

    print("\nAs tau -> 0, soft H converges to hard H = 189")
    print("This makes tournament invariants DIFFERENTIABLE probes!")


def profile_comparison():
    """
    Compare tournament profiles across different attention pattern types.
    This is what TournamentProbe would display.
    """
    print("\n\nTOURNAMENT PROFILE COMPARISON")
    print("=" * 60)

    np.random.seed(42)
    profiles = []

    # 1. Paley T_7
    p = 7
    qr = set()
    for k in range(1, p):
        qr.add((k*k) % p)
    T = np.zeros((p, p), dtype=int)
    for i in range(p):
        for j in range(p):
            if i != j and (j-i) % p in qr:
                T[i][j] = 1
    profiles.append(tournament_profile(T, "Paley T_7"))

    # 2. Transitive T_7
    T = np.zeros((7, 7), dtype=int)
    for i in range(7):
        for j in range(i+1, 7):
            T[i][j] = 1
    profiles.append(tournament_profile(T, "Transitive T_7"))

    # 3. Random tournaments
    for trial in range(5):
        T = np.zeros((7, 7), dtype=int)
        for i in range(7):
            for j in range(i+1, 7):
                if np.random.random() < 0.5:
                    T[i][j] = 1
                else:
                    T[j][i] = 1
        profiles.append(tournament_profile(T, f"Random T_7 #{trial+1}"))

    # 4. Random attention -> tournament
    for trial in range(3):
        logits = np.random.randn(7, 7)
        exp_l = np.exp(logits - logits.max(axis=1, keepdims=True))
        A = exp_l / exp_l.sum(axis=1, keepdims=True)
        T = np.zeros((7, 7), dtype=int)
        for i in range(7):
            for j in range(i+1, 7):
                if A[i, j] > A[j, i]:
                    T[i][j] = 1
                else:
                    T[j][i] = 1
        profiles.append(tournament_profile(T, f"Attn-Tournament #{trial+1}"))

    # Display
    print(f"\n{'Name':25s} {'H':>6s} {'c3':>4s} {'trans':>6s} "
          f"{'score_var':>9s} {'reg':>5s} {'anti_E':>7s}")
    print("-" * 75)
    for p in profiles:
        H = p.get('H', '?')
        c3 = p.get('c3', '?')
        trans = f"{p['transitivity']:.3f}"
        sv = f"{p['score_variance']:.3f}"
        reg = "Y" if p['is_regular'] else "N"
        ae = f"{p['anti_energy']:.3f}"
        print(f"{p['name']:25s} {str(H):>6s} {str(c3):>4s} {trans:>6s} "
              f"{sv:>9s} {reg:>5s} {ae:>7s}")

    print("\nKey: H=Hamiltonian paths, c3=directed 3-cycles, "
          "trans=transitivity, anti_E=Cartan antisymmetric energy")
    print("Paley T_7 has MAXIMUM H=189 among ALL 7-vertex tournaments")


if __name__ == "__main__":
    demonstrate_soft_tournament_bridge()
    profile_comparison()
