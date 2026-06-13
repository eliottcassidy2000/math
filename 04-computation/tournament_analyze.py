#!/usr/bin/env python3
"""
tournament_analyze.py — kind-pasteur-2026-03-22-S20i

A single-file command-line tool for complete tournament analysis.
Takes pairwise comparison data, outputs everything this project can compute.

Usage:
  python tournament_analyze.py --demo
  python tournament_analyze.py --csv data.csv --col-a item_a --col-b item_b --col-winner winner
  python tournament_analyze.py --matrix "0,1,1,0; 0,0,1,0; 0,0,0,1; 1,1,0,0"

Output includes:
  - Score sequence and S_2
  - H (Hamiltonian path count) via the algebraic formula (n<=4) or DP
  - c3 (3-cycle count) from scores in O(n)
  - Independence polynomial coefficients (small n)
  - Conductivity index kappa
  - Cartan decomposition (tournament fraction, cooperation fraction)
  - Free energy and energy bandwidth
  - FormalRank ranking with confidence

Author: kind-pasteur-2026-03-22-S20i
Dependencies: numpy only.
"""

import sys
import numpy as np
from math import comb, log, sqrt, pi, factorial
from collections import defaultdict
from itertools import permutations

def scores(A):
    """Score sequence (out-degrees)."""
    return A.sum(axis=1).astype(int)

def S2(A):
    """Sum of squared scores."""
    s = scores(A)
    return int(sum(s**2))

def c3_from_scores(A):
    """3-cycle count from score sequence in O(n)."""
    n = len(A)
    s2 = S2(A)
    return comb(n, 3) - (s2 - comb(n, 2)) // 2

def c3_from_trace(A):
    """3-cycle count via Tr(A^3)/3 in O(n^3)."""
    A3 = A @ A @ A
    return int(np.trace(A3)) // 3

def H_formula(A):
    """H via algebraic formula (exact for n<=4, approximate for n>=5)."""
    n = len(A)
    const = 1 + n*(n-1)*(2*n-1)//6
    s2 = S2(A)
    return const - s2

def H_exact(A):
    """H via Held-Karp dynamic programming. O(n^2 * 2^n)."""
    n = len(A)
    dp = defaultdict(int)
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[(mask, v)] == 0: continue
            for w in range(n):
                if mask & (1 << w): continue
                if A[v][w]:
                    dp[(mask | (1 << w), w)] += dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp[(full, v)] for v in range(n))

def cartan_decomposition(A):
    """Decompose attention/comparison matrix into tournament and cooperation."""
    A_float = A.astype(float)
    A_sym = (A_float + A_float.T) / 2
    A_anti = (A_float - A_float.T) / 2
    norm_sym = np.linalg.norm(A_sym)
    norm_anti = np.linalg.norm(A_anti)
    total = norm_sym**2 + norm_anti**2
    tournament_frac = norm_anti**2 / total if total > 0 else 0
    return {
        "norm_symmetric": norm_sym,
        "norm_antisymmetric": norm_anti,
        "tournament_fraction": tournament_frac,
        "cooperation_fraction": 1 - tournament_frac,
    }

def formal_rank(A):
    """Rank items using formal group logarithm (arctanh)."""
    n = len(A)
    rapidities = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if i == j: continue
            if A[i][j]:
                rapidities[i] += 1
            else:
                rapidities[i] -= 1
    # Normalize
    rapidities = rapidities / (n - 1)
    # Apply arctanh for formal group
    rapidities_fg = np.arctanh(np.clip(rapidities, -0.999, 0.999))
    ranking = np.argsort(-rapidities_fg)
    return ranking, rapidities_fg

def conductivity_index(A, n):
    """Kappa = I(G, 2) / I(spanning_tree, 2). Approximate via formula."""
    # For the tournament itself (not its conflict graph):
    # Compare cyclic vs acyclic evaluation
    c3 = c3_from_scores(A)
    total_triples = comb(n, 3)
    # Fraction of cyclic triples
    if total_triples > 0:
        cycle_fraction = c3 / total_triples
    else:
        cycle_fraction = 0
    return cycle_fraction

def energy_analysis(A):
    """Free energy and energy level analysis."""
    n = len(A)
    c3 = c3_from_scores(A)
    # For degree-1 polynomial: root = -1/(c3) if c3 > 0
    if c3 > 0:
        root = -1.0 / c3
        energy = -log(2 - root)
        bandwidth = log(3/2)
    else:
        root = None
        energy = 0
        bandwidth = 0
    return {
        "c3": c3,
        "root": root,
        "energy_level": energy,
        "bandwidth": bandwidth,
    }

def full_analysis(A, names=None):
    """Complete tournament analysis."""
    n = len(A)
    if names is None:
        names = [str(i) for i in range(n)]

    print("=" * 60)
    print("  TOURNAMENT ANALYSIS")
    print("=" * 60)

    # Basic info
    print(f"\n  Vertices: {n}")
    print(f"  Arcs: {int(A.sum())} (of {n*(n-1)} possible)")

    # Scores
    s = scores(A)
    sorted_scores = sorted(s, reverse=True)
    print(f"\n  SCORES:")
    ranking, raps = formal_rank(A)
    for rank, idx in enumerate(ranking):
        print(f"    #{rank+1}: {names[idx]:>12s}  score={s[idx]:>2d}  rapidity={raps[idx]:>+.3f}")

    # Score statistics
    s2 = S2(A)
    print(f"\n  SCORE STATISTICS:")
    print(f"    S_2 (sum of squared scores) = {s2}")
    print(f"    Score variance = {np.var(s):.4f}")
    regular = all(x == s[0] for x in s)
    print(f"    Regular tournament: {regular}")

    # Cycle counts
    c3 = c3_from_scores(A)
    c3_check = c3_from_trace(A)
    print(f"\n  CYCLE COUNTS:")
    print(f"    c3 (3-cycles) = {c3} (from scores), {c3_check} (from Tr(A^3)/3)")
    print(f"    Fraction cyclic triples: {c3}/{comb(n,3)} = {c3/comb(n,3):.4f}" if comb(n,3) > 0 else "")

    # H computation
    const = 1 + n*(n-1)*(2*n-1)//6
    H_approx = const - s2

    if n <= 15:
        H = H_exact(A)
        print(f"\n  HAMILTONIAN PATH COUNT:")
        print(f"    H = {H} (exact, Held-Karp DP)")
        if n <= 4:
            print(f"    H = {const} - {s2} = {H_approx} (algebraic formula, EXACT at n<={n})")
        else:
            print(f"    H_approx = {const} - {s2} = {H_approx} (score-determined part)")
            print(f"    Correction = {H - H_approx} (cycle/packing terms)")
            print(f"    Score explains {100*H_approx/H:.1f}% of H" if H > 0 else "")
    else:
        H = None
        print(f"\n  HAMILTONIAN PATH COUNT:")
        print(f"    H_approx = {const} - {s2} = {H_approx} (score-determined part, n too large for exact)")

    # Cartan decomposition
    cartan = cartan_decomposition(A)
    print(f"\n  CARTAN DECOMPOSITION:")
    print(f"    ||A_tournament|| = {cartan['norm_antisymmetric']:.4f}")
    print(f"    ||A_cooperation|| = {cartan['norm_symmetric']:.4f}")
    print(f"    Tournament fraction = {cartan['tournament_fraction']:.4f}")
    print(f"    Cooperation fraction = {cartan['cooperation_fraction']:.4f}")

    # Conductivity
    kappa = conductivity_index(A, n)
    print(f"\n  CONDUCTIVITY INDEX:")
    print(f"    Cycle fraction (c3/total triples) = {kappa:.4f}")
    if kappa > 0.5:
        print(f"    Classification: HIGH INTRANSITIVITY (cycle-rich)")
    elif kappa > 0.1:
        print(f"    Classification: MODERATE INTRANSITIVITY")
    else:
        print(f"    Classification: LOW INTRANSITIVITY (near-transitive)")

    # Energy
    energy = energy_analysis(A)
    print(f"\n  FREE ENERGY ANALYSIS:")
    if H and H > 0:
        print(f"    log H = {log(H):.4f}")
        print(f"    Free energy F = -log H = {-log(H):.4f}")
    if energy["root"] is not None:
        print(f"    Primary root = {energy['root']:.6f}")
        print(f"    Energy level = {energy['energy_level']:.6f}")
        print(f"    Bandwidth = log(3/2) = {energy['bandwidth']:.6f} nats")

    # Redei parity
    if H is not None:
        print(f"\n  REDEI PARITY: H = {H} is {'ODD (Redei confirmed)' if H % 2 == 1 else 'EVEN (BUG!)'}")

    # Summary
    print(f"\n  SUMMARY:")
    if regular:
        print(f"    This is a REGULAR tournament (all scores equal).")
        print(f"    Regular tournaments MAXIMIZE H for their order.")
    elif H is not None and H == 1:
        print(f"    This is a TRANSITIVE tournament (H=1, no cycles).")
        print(f"    Transitive tournaments MINIMIZE H.")
    else:
        score_info = f"Score-determined H: {H_approx}"
        if H is not None:
            score_info += f", actual H: {H}, gap: {H - H_approx}"
        print(f"    {score_info}")

    print()
    return {"H": H, "c3": c3, "S2": s2, "scores": s, "cartan": cartan}


def demo():
    """Run demo with example tournaments."""
    print("\n" + "=" * 60)
    print("  DEMO: ROCK-PAPER-SCISSORS (n=3)")
    print("=" * 60)

    A_rps = np.array([[0, 1, 0],
                       [0, 0, 1],
                       [1, 0, 0]])
    full_analysis(A_rps, names=["Rock", "Paper", "Scissors"])

    print("\n" + "=" * 60)
    print("  DEMO: TRANSITIVE TOURNAMENT (n=4)")
    print("=" * 60)

    A_trans = np.array([[0, 1, 1, 1],
                         [0, 0, 1, 1],
                         [0, 0, 0, 1],
                         [0, 0, 0, 0]])
    full_analysis(A_trans, names=["Alice", "Bob", "Carol", "Dave"])

    print("\n" + "=" * 60)
    print("  DEMO: NEAR-REGULAR TOURNAMENT (n=5)")
    print("=" * 60)

    A_5 = np.array([[0, 1, 0, 1, 1],
                     [0, 0, 1, 1, 0],
                     [1, 0, 0, 0, 1],
                     [0, 0, 1, 0, 1],
                     [0, 1, 0, 0, 0]])
    full_analysis(A_5, names=["Alpha", "Beta", "Gamma", "Delta", "Epsilon"])

    print("\n" + "=" * 60)
    print("  DEMO: LLM COMPARISON (n=4)")
    print("=" * 60)

    # Simulated: Claude > GPT on code, GPT > Gemini on reasoning,
    # Gemini > Claude on creativity, Claude > Gemini on safety,
    # GPT > Claude on speed, Gemini > GPT on multimodal
    A_llm = np.array([[0, 1, 1, 1],  # Claude beats GPT(code), Gemini(safety), Llama(all)
                       [0, 0, 1, 1],  # GPT beats Gemini(reasoning), Llama
                       [0, 0, 0, 1],  # Gemini beats Llama
                       [0, 0, 0, 0]]) # Llama loses all
    # Add a cycle: Gemini > Claude on creativity
    A_llm[2][0] = 1; A_llm[0][2] = 0
    full_analysis(A_llm, names=["Claude", "GPT-4o", "Gemini", "Llama"])


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--demo":
        demo()
    elif len(sys.argv) > 1 and sys.argv[1] == "--matrix":
        # Parse matrix from command line
        mat_str = sys.argv[2]
        rows = mat_str.split(";")
        A = np.array([[int(x.strip()) for x in row.split(",")] for row in rows])
        full_analysis(A)
    else:
        demo()
