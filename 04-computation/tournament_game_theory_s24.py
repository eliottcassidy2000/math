#!/usr/bin/env python3
"""
tournament_game_theory_s24.py — opus-2026-04-05-S24

Connect tournaments to n-person game theory (von Neumann, Morgenstern, Nash).

A tournament on n vertices IS:
- A 2-player zero-sum game (adjacency = payoff matrix)
- An n-person cooperative game with v(S) = H(T|_S)
- A social choice profile (Condorcet structure)

We compute:
1. Mixed Nash equilibrium of the zero-sum game
2. Shapley value of the "Hamiltonian game" v(S) = H(T|_S)
3. Core of the H-game (when non-empty)
4. Tournament solution concepts: top cycle, uncovered set, Copeland score
5. Markov stationary distribution (random out-neighbor walk)

Then relate all of these to H(T), c₃, and the metagraph structure.

KEY INSIGHT: Claim A says H(T) - H(T-v) = 2*Σ_{C∋v} μ(C).
This IS the marginal contribution of v to the grand coalition H(T).
The Shapley value weights this across ALL coalitions, giving a
"fair" allocation of Hamiltonian path count to vertices.
"""

import numpy as np
from scipy.optimize import linprog
from itertools import combinations, permutations
from collections import defaultdict, Counter
from math import factorial
import time
import sys

sys.stdout.reconfigure(encoding='utf-8')

# ================================================================
# Tournament primitives
# ================================================================

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_ham_paths(A, n):
    full = (1 << n) - 1
    dp = np.zeros((1 << n, n), dtype=int)
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if dp[mask][v] == 0 or not (mask & (1 << v)):
                continue
            for u in range(n):
                if not (mask & (1 << u)) and A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return int(sum(dp[full]))

def score_c3(A, n):
    total = n * (n-1) * (n-2) // 6
    for v in range(n):
        s = int(A[v].sum())
        total -= s * (s-1) // 2
    return total

def canonical_form(A, n):
    best = None
    for perm in permutations(range(n)):
        bits = 0
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if A[perm[i]][perm[j]] == 1:
                    bits |= (1 << idx)
                idx += 1
        if best is None or bits < best:
            best = bits
    return best

# ================================================================
# 1. Nash equilibrium of the zero-sum tournament game
# ================================================================

def nash_equilibrium(A, n):
    """
    The tournament zero-sum game: Player 1 picks row i, Player 2 picks col j.
    Payoff to Player 1: M[i][j] = 1 if A[i][j]=1 (i beats j), -1 if A[j][i]=1, 0 if i=j.

    For a symmetric zero-sum game, the Nash equilibrium is a SINGLE mixed strategy p
    solving: min_p max_j (Mp)_j, or equivalently the LP:
        max v s.t. Mp >= v*1, p >= 0, sum(p) = 1

    Returns (p, value) where p is the equilibrium strategy and value is the game value.
    """
    M = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            if i == j:
                M[i][j] = 0
            elif A[i][j] == 1:
                M[i][j] = 1
            else:
                M[i][j] = -1

    # LP: maximize v
    # Variables: [p_0, p_1, ..., p_{n-1}, v]
    # Constraints:
    #   M^T p >= v*1  =>  -M^T p + v*1 <= 0
    #   sum(p) = 1
    #   p >= 0
    c = np.zeros(n + 1)
    c[n] = -1  # maximize v = minimize -v

    # Inequality: -M^T @ p + v <= 0 for each column j
    A_ub = np.zeros((n, n + 1))
    for j in range(n):
        for i in range(n):
            A_ub[j][i] = -M[i][j]  # -M^T p
        A_ub[j][n] = 1  # + v
    b_ub = np.zeros(n)

    # Equality: sum(p) = 1
    A_eq = np.zeros((1, n + 1))
    A_eq[0, :n] = 1
    b_eq = np.array([1.0])

    bounds = [(0, None)] * n + [(None, None)]  # p >= 0, v unrestricted

    result = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method='highs')

    if result.success:
        p = result.x[:n]
        v = result.x[n]
        # Clean up near-zero entries
        p = np.maximum(p, 0)
        p /= p.sum()
        return p, v
    else:
        return None, None

def nash_support_size(p, tol=1e-6):
    return int(np.sum(p > tol))

def nash_entropy(p, tol=1e-10):
    h = 0
    for pi in p:
        if pi > tol:
            h -= pi * np.log2(pi)
    return h

# ================================================================
# 2. Shapley value of the Hamiltonian game
# ================================================================

def shapley_value_H_game(A, n):
    """
    Cooperative game where v(S) = H(T|_S) for each coalition S.
    The Shapley value of player i is:
    φ_i = Σ_{S ⊆ N\{i}} [|S|!(n-|S|-1)!/n!] * [v(S∪{i}) - v(S)]

    This measures the "fair" contribution of vertex i to Hamiltonian path count.
    """
    # Precompute H for all subsets of vertices
    H_cache = {}
    for mask in range(1, 1 << n):
        verts = [v for v in range(n) if mask & (1 << v)]
        k = len(verts)
        if k == 1:
            H_cache[mask] = 1
        else:
            # Extract subtournament
            sub = np.zeros((k, k), dtype=int)
            for a in range(k):
                for b in range(k):
                    sub[a][b] = A[verts[a]][verts[b]]
            H_cache[mask] = count_ham_paths(sub, k)

    # Compute Shapley value for each vertex
    phi = np.zeros(n)
    fact_n = factorial(n)
    for i in range(n):
        for mask in range(0, 1 << n):
            if mask & (1 << i):
                continue  # i must not be in S
            verts_in_S = [v for v in range(n) if mask & (1 << v)]
            s = len(verts_in_S)
            weight = factorial(s) * factorial(n - s - 1) / fact_n

            mask_with_i = mask | (1 << i)
            v_S = H_cache.get(mask, 0) if mask > 0 else 0
            v_Si = H_cache[mask_with_i]

            phi[i] += weight * (v_Si - v_S)

    return phi

# ================================================================
# 3. Core of the H-game
# ================================================================

def h_game_core_check(A, n, H_cache=None):
    """
    The core of the cooperative game v(S) = H(T|_S) is:
    {x : sum(x) = v(N) = H(T), sum(x_i for i in S) >= v(S) for all S}

    Check if core is non-empty via LP.
    Returns (is_nonempty, x_core) where x_core is a core allocation (if exists).
    """
    if H_cache is None:
        H_cache = {}
        for mask in range(1, 1 << n):
            verts = [v for v in range(n) if mask & (1 << v)]
            k = len(verts)
            if k == 1:
                H_cache[mask] = 1
            else:
                sub = np.zeros((k, k), dtype=int)
                for a in range(k):
                    for b in range(k):
                        sub[a][b] = A[verts[a]][verts[b]]
                H_cache[mask] = count_ham_paths(sub, k)

    full_mask = (1 << n) - 1
    H_total = H_cache[full_mask]

    # LP: minimize 0 (feasibility check)
    # Variables: x_0, ..., x_{n-1}
    # Constraints:
    #   sum(x) = H_total
    #   sum(x_i for i in S) >= v(S) for all S
    #   (we can also add x_i >= 0 for individual rationality)

    num_constraints = 0
    A_ub_list = []
    b_ub_list = []

    for mask in range(1, full_mask):  # all proper non-empty subsets
        verts_in_S = [v for v in range(n) if mask & (1 << v)]
        v_S = H_cache[mask]
        # Constraint: -sum(x_i for i in S) <= -v(S)
        row = np.zeros(n)
        for v in verts_in_S:
            row[v] = -1
        A_ub_list.append(row)
        b_ub_list.append(-v_S)

    A_ub = np.array(A_ub_list) if A_ub_list else None
    b_ub = np.array(b_ub_list) if b_ub_list else None

    A_eq = np.ones((1, n))
    b_eq = np.array([float(H_total)])

    c = np.zeros(n)  # feasibility

    result = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=[(None, None)]*n, method='highs')

    if result.success:
        return True, result.x
    else:
        return False, None

# ================================================================
# 4. Tournament solution concepts
# ================================================================

def copeland_scores(A, n):
    """Copeland score = out-degree. Higher = more dominant."""
    return [int(A[i].sum()) for i in range(n)]

def top_cycle(A, n):
    """Top cycle (Smith set): smallest set that beats all outsiders.
    Equivalently: strongly connected components, take the top one."""
    # Find SCCs via DFS
    visited = [False] * n
    finish = []

    def dfs1(v):
        visited[v] = True
        for u in range(n):
            if A[v][u] and not visited[u]:
                dfs1(u)
        finish.append(v)

    for v in range(n):
        if not visited[v]:
            dfs1(v)

    # Reverse graph
    AT = A.T

    visited2 = [False] * n
    sccs = []

    def dfs2(v, scc):
        visited2[v] = True
        scc.append(v)
        for u in range(n):
            if AT[v][u] and not visited2[u]:
                dfs2(u, scc)

    for v in reversed(finish):
        if not visited2[v]:
            scc = []
            dfs2(v, scc)
            sccs.append(scc)

    # Top SCC: the one with no incoming edges from other SCCs
    # In a tournament, the top cycle is the first SCC in topological order
    # which is the first SCC found in the second pass
    return sorted(sccs[0])

def uncovered_set(A, n):
    """
    Vertex i covers vertex j if: i beats j AND for every k that j beats, i also beats k.
    The uncovered set = {i : no j covers i}.
    """
    uncovered = []
    for i in range(n):
        is_covered = False
        for j in range(n):
            if j == i:
                continue
            if A[j][i] == 0:
                continue  # j doesn't beat i
            # Check if j covers i: j beats i, and for all k, if i beats k then j beats k
            covers = True
            for k in range(n):
                if k == i or k == j:
                    continue
                if A[i][k] and not A[j][k]:
                    covers = False
                    break
            if covers:
                is_covered = True
                break
        if not is_covered:
            uncovered.append(i)
    return uncovered

def markov_stationary(A, n):
    """
    Random walk: from vertex i, move to a random out-neighbor.
    Transition matrix P[i][j] = A[i][j] / out_degree(i).
    Stationary distribution: pi P = pi.
    """
    scores = [int(A[i].sum()) for i in range(n)]
    P = np.zeros((n, n))
    for i in range(n):
        if scores[i] > 0:
            for j in range(n):
                P[i][j] = A[i][j] / scores[i]

    # Find stationary: left eigenvector of P with eigenvalue 1
    evals, evecs = np.linalg.eig(P.T)
    idx = np.argmin(np.abs(evals - 1.0))
    pi = np.real(evecs[:, idx])
    pi = np.abs(pi)
    pi /= pi.sum()
    return pi

# ================================================================
# Main analysis
# ================================================================

def analyze_tournament(A, n, label=""):
    """Full game-theoretic analysis of one tournament."""
    H = count_ham_paths(A, n)
    c3 = score_c3(A, n)
    scores = copeland_scores(A, n)

    result = {'H': H, 'c3': c3, 'scores': tuple(sorted(scores))}

    # Nash equilibrium
    p_nash, v_nash = nash_equilibrium(A, n)
    if p_nash is not None:
        result['nash_p'] = p_nash
        result['nash_value'] = v_nash
        result['nash_support'] = nash_support_size(p_nash)
        result['nash_entropy'] = nash_entropy(p_nash)
    else:
        result['nash_support'] = -1

    # Shapley value
    phi = shapley_value_H_game(A, n)
    result['shapley'] = phi
    result['shapley_max'] = float(np.max(phi))
    result['shapley_min'] = float(np.min(phi))
    result['shapley_spread'] = float(np.max(phi) - np.min(phi))
    result['shapley_sorted'] = tuple(sorted(phi))

    # Core
    core_exists, x_core = h_game_core_check(A, n)
    result['core_exists'] = core_exists
    if core_exists:
        result['core_point'] = x_core

    # Solution concepts
    result['top_cycle_size'] = len(top_cycle(A, n))
    result['uncovered_size'] = len(uncovered_set(A, n))

    # Markov
    pi_markov = markov_stationary(A, n)
    result['markov_entropy'] = nash_entropy(pi_markov)
    result['markov_max'] = float(np.max(pi_markov))

    return result


print("=" * 70)
print("TOURNAMENT GAME THEORY ANALYSIS")
print("von Neumann, Morgenstern, Nash meets Redei")
print("=" * 70)

for n in [4, 5]:
    m = n * (n - 1) // 2
    total = 1 << m
    print(f"\n{'='*60}")
    print(f"n = {n}")
    print(f"{'='*60}")

    # Group by iso class
    class_results = {}
    canon_map = {}

    t0 = time.time()
    for bits in range(total):
        A = bits_to_adj(bits, n)
        cf = canonical_form(A, n)
        if cf not in canon_map:
            canon_map[cf] = bits
            result = analyze_tournament(A, n)
            class_results[cf] = result

    elapsed = time.time() - t0
    print(f"\n  {len(class_results)} iso classes analyzed in {elapsed:.1f}s\n")

    # Sort by H
    sorted_classes = sorted(class_results.items(), key=lambda x: (x[1]['H'], x[1]['c3']))

    print(f"  {'H':>3} {'c3':>3} {'score':>15} | {'Nash':>4} {'val':>5} {'ent':>5} | "
          f"{'Shap_spread':>11} {'Shap_sorted':>25} | {'Core':>4} {'TC':>2} {'UC':>2} {'M_ent':>5}")
    print(f"  {'-'*3} {'-'*3} {'-'*15}-|-{'-'*4}-{'-'*5}-{'-'*5}-|-"
          f"{'-'*11}-{'-'*25}-|-{'-'*4}-{'-'*2}-{'-'*2}-{'-'*5}")

    for cf, r in sorted_classes:
        score_str = str(r['scores'])
        nash_sup = r['nash_support']
        nash_val = r.get('nash_value', 0)
        nash_ent = r.get('nash_entropy', 0)
        shap_spread = r['shapley_spread']
        shap_sorted = tuple(round(x, 2) for x in r['shapley_sorted'])
        core = "YES" if r['core_exists'] else "no"
        tc = r['top_cycle_size']
        uc = r['uncovered_size']
        m_ent = r['markov_entropy']

        print(f"  {r['H']:>3} {r['c3']:>3} {score_str:>15} | {nash_sup:>4} {nash_val:>5.2f} {nash_ent:>5.2f} | "
              f"{shap_spread:>11.3f} {str(shap_sorted):>25} | {core:>4} {tc:>2} {uc:>2} {m_ent:>5.2f}")

    # Key correlations
    print(f"\n  KEY RELATIONSHIPS:")
    Hs = [r['H'] for _, r in sorted_classes]
    c3s = [r['c3'] for _, r in sorted_classes]
    nash_sups = [r['nash_support'] for _, r in sorted_classes]
    nash_ents = [r.get('nash_entropy', 0) for _, r in sorted_classes]
    shap_spreads = [r['shapley_spread'] for _, r in sorted_classes]
    tc_sizes = [r['top_cycle_size'] for _, r in sorted_classes]
    uc_sizes = [r['uncovered_size'] for _, r in sorted_classes]
    m_ents = [r['markov_entropy'] for _, r in sorted_classes]

    if len(Hs) > 2:
        H_arr = np.array(Hs, dtype=float)
        for name, vals in [
            ("Nash support", nash_sups),
            ("Nash entropy", nash_ents),
            ("Shapley spread", shap_spreads),
            ("Top cycle size", tc_sizes),
            ("Uncovered set size", uc_sizes),
            ("Markov entropy", m_ents),
        ]:
            v = np.array(vals, dtype=float)
            if np.std(v) > 0 and np.std(H_arr) > 0:
                corr = np.corrcoef(H_arr, v)[0, 1]
                print(f"    corr(H, {name:20s}) = {corr:+.4f}")

    # Does every vertex's Shapley value sum to H?
    print(f"\n  SHAPLEY VALUE PROPERTIES:")
    for cf, r in sorted_classes[:3]:
        phi = r['shapley']
        print(f"    H={r['H']}: sum(phi)={sum(phi):.2f} (should = H), phi={[round(x,3) for x in phi]}")

    # Core analysis
    core_yes = sum(1 for _, r in sorted_classes if r['core_exists'])
    core_no = sum(1 for _, r in sorted_classes if not r['core_exists'])
    print(f"\n  CORE: {core_yes} classes with non-empty core, {core_no} with empty core")
    for cf, r in sorted_classes:
        if r['core_exists']:
            print(f"    H={r['H']}, score={r['scores']}: core point = {[round(x,2) for x in r['core_point']]}")

    # What makes the Shapley spread large?
    print(f"\n  SHAPLEY SPREAD (inequality among vertices):")
    for cf, r in sorted(class_results.items(), key=lambda x: -x[1]['shapley_spread']):
        print(f"    H={r['H']}, c3={r['c3']}, score={r['scores']}: spread={r['shapley_spread']:.3f}, "
              f"sorted_phi={tuple(round(x,2) for x in r['shapley_sorted'])}")
