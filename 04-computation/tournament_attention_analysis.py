#!/usr/bin/env python3
"""
tournament_attention_analysis.py -- kind-pasteur-2026-03-21-S12

Apply tournament theory to transformer attention patterns.

Key ideas:
1. Threshold attention matrices to get tournaments on token sequences
2. Compute tournament invariants: H(T), I(Omega(T), x), Betti numbers
3. Cartan decomposition: A_sym = (A+A^T)/2, A_anti = (A-A^T)/2
4. Measure which component carries more task-relevant information

This script works WITHOUT a GPU — it uses pre-computed attention patterns
or small models. For actual LLM analysis, use a GPU machine.

Connections to Napolitano (2026) "Mathematics Is All You Need":
- Paper claims gl(4,R) fiber bundle in hidden states
- Our approach: extract TOURNAMENT structure from attention, use OCF
- Testing: does A_sym (Napolitano's "dark") or A_anti (our tournaments)
  carry more information?

Author: kind-pasteur-2026-03-21-S12
"""

import numpy as np
from itertools import combinations
from collections import defaultdict


# ============================================================
# CORE TOURNAMENT FUNCTIONS (from our existing codebase)
# ============================================================

def attention_to_tournament(A, threshold='median'):
    """
    Convert attention matrix A (n×n, non-negative) to tournament.

    Method: For each pair (i,j), i->j if A[i,j] > A[j,i].
    Ties broken by index (i->j if i < j). This gives a tournament.

    The 'threshold' parameter is unused for pairwise comparison method.
    We use pairwise comparison because it's the natural tournament construction:
    i "dominates" j in the attention pattern iff i attends to j more than j attends to i.

    Returns: T[i][j] = 1 if i->j, 0 otherwise
    """
    n = A.shape[0]
    T = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            if A[i, j] > A[j, i]:
                T[i][j] = 1
            elif A[j, i] > A[i, j]:
                T[j][i] = 1
            else:
                # tie-break by index
                T[i][j] = 1
    return T


def count_ham_paths(T):
    """Count Hamiltonian paths in tournament T using DP. O(2^n * n^2)."""
    n = T.shape[0]
    if n > 20:
        raise ValueError(f"n={n} too large for exact HP count")
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
    return sum(dp[(full, v)] for v in range(n))


def find_odd_cycles(T, max_length=None):
    """
    Find all directed odd cycles in tournament T.
    Returns list of vertex tuples, each a directed cycle.
    """
    n = T.shape[0]
    if max_length is None:
        max_length = n
    cycles = []
    for length in range(3, max_length + 1, 2):  # odd lengths only
        for verts in combinations(range(n), length):
            # Check all directed Hamiltonian cycles on this vertex set
            from itertools import permutations as perms
            vlist = list(verts)
            # Fix first vertex to avoid counting rotations
            first = vlist[0]
            rest = vlist[1:]
            for perm in perms(rest):
                path = [first] + list(perm)
                # Check if this forms a directed cycle
                is_cycle = True
                for k in range(length):
                    if not T[path[k]][path[(k+1) % length]]:
                        is_cycle = False
                        break
                if is_cycle:
                    cycles.append(tuple(path))
    return cycles


def build_conflict_graph(cycles):
    """
    Build conflict graph Omega: vertices = cycles, edges = shared vertices.
    Returns adjacency matrix of Omega.
    """
    k = len(cycles)
    adj = np.zeros((k, k), dtype=int)
    for i in range(k):
        set_i = set(cycles[i])
        for j in range(i+1, k):
            set_j = set(cycles[j])
            if set_i & set_j:  # shared vertices
                adj[i][j] = adj[j][i] = 1
    return adj


def independence_polynomial(adj, x_val=2):
    """
    Compute I(G, x) at x = x_val by enumerating independent sets.
    Works for small graphs only (|V| < 25).
    """
    k = adj.shape[0]
    if k > 24:
        raise ValueError(f"Graph too large: {k} vertices")

    total = 0
    for mask in range(1 << k):
        # Check if mask is an independent set
        is_indep = True
        verts = [v for v in range(k) if mask & (1 << v)]
        for a in range(len(verts)):
            for b in range(a+1, len(verts)):
                if adj[verts[a]][verts[b]]:
                    is_indep = False
                    break
            if not is_indep:
                break
        if is_indep:
            total += x_val ** len(verts)
    return total


def compute_scores(T):
    """Score sequence of tournament (out-degrees)."""
    return sorted(T.sum(axis=1))


# ============================================================
# CARTAN DECOMPOSITION OF ATTENTION
# ============================================================

def cartan_decompose(A):
    """
    Decompose matrix A into:
    - A_anti = (A - A^T)/2  (antisymmetric = "active" = tournament structure)
    - A_sym = (A + A^T)/2   (symmetric = "dark" = similarity structure)
    - A_scalar = tr(A)/n * I (scalar = center of gl(n))
    - A_traceless_sym = A_sym - A_scalar (traceless symmetric = non-compact)

    This is the Cartan decomposition gl(n,R) = so(n) + p + R*I
    where:
    - so(n) has dim n(n-1)/2 ["active" in Napolitano's terminology]
    - p has dim n(n+1)/2 - 1 ["dark" in Napolitano's terminology]
    - R*I has dim 1 ["null" in Killing form]

    Returns dict with components and their Frobenius norms.
    """
    n = A.shape[0]
    A_anti = (A - A.T) / 2
    A_sym = (A + A.T) / 2
    scalar = np.trace(A) / n
    A_scalar = scalar * np.eye(n)
    A_traceless_sym = A_sym - A_scalar

    norm_anti = np.linalg.norm(A_anti, 'fro')
    norm_sym = np.linalg.norm(A_traceless_sym, 'fro')
    norm_scalar = abs(scalar) * np.sqrt(n)
    norm_total = np.linalg.norm(A, 'fro')

    return {
        'antisymmetric': A_anti,       # so(n) part = tournament
        'symmetric_traceless': A_traceless_sym,  # p part = "dark"
        'scalar': A_scalar,            # center
        'norm_anti': norm_anti,
        'norm_sym': norm_sym,
        'norm_scalar': norm_scalar,
        'norm_total': norm_total,
        'anti_fraction': norm_anti**2 / norm_total**2 if norm_total > 0 else 0,
        'dark_fraction': norm_sym**2 / norm_total**2 if norm_total > 0 else 0,
        'dimensions': {
            'active': n*(n-1)//2,
            'dark': n*(n+1)//2 - 1,
            'null': 1,
            'total': n*n
        }
    }


# ============================================================
# SYNTHETIC ATTENTION PATTERNS FOR TESTING
# ============================================================

def random_attention_matrix(n, temperature=1.0, seed=None):
    """
    Generate a random attention matrix (each row sums to 1).
    Uses softmax of random logits.
    """
    if seed is not None:
        np.random.seed(seed)
    logits = np.random.randn(n, n) / temperature
    # Apply softmax row-wise
    exp_logits = np.exp(logits - logits.max(axis=1, keepdims=True))
    A = exp_logits / exp_logits.sum(axis=1, keepdims=True)
    return A


def paley_attention(p):
    """
    Create attention matrix corresponding to Paley tournament T_p.
    A[i,j] = 1/p if i->j in Paley, 0 otherwise, then row-normalized.

    This tests whether Paley structure (our H-maximizer) gives
    distinctive Cartan decomposition properties.
    """
    if p < 3 or not all(p % i != 0 for i in range(2, int(p**0.5)+1)):
        raise ValueError(f"{p} is not prime")

    # Quadratic residues mod p
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

    # Row-normalize
    row_sums = A.sum(axis=1, keepdims=True)
    A = A / np.maximum(row_sums, 1e-10)
    return A


def transitive_attention(n):
    """
    Attention matrix for transitive tournament: i->j iff i < j.
    All weight on dominated tokens.
    """
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1.0
    row_sums = A.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1  # vertex 0 has no outgoing
    A = A / row_sums
    return A


# ============================================================
# ANALYSIS FUNCTIONS
# ============================================================

def analyze_attention_tournament(A, name="", max_n_for_cycles=9):
    """
    Full tournament analysis of an attention matrix.

    Computes:
    1. Tournament T from pairwise comparison
    2. H(T) = number of Hamiltonian paths
    3. Score sequence
    4. Odd cycles and conflict graph (if n small enough)
    5. I(Omega(T), 2) = should equal H(T) by OCF
    6. Cartan decomposition
    """
    n = A.shape[0]
    print(f"\n{'='*60}")
    print(f"ATTENTION TOURNAMENT ANALYSIS: {name}")
    print(f"{'='*60}")
    print(f"Token count (tournament size): n = {n}")

    # 1. Convert to tournament
    T = attention_to_tournament(A)
    scores = compute_scores(T)
    print(f"Score sequence: {list(scores)}")
    is_regular = all(s == (n-1)/2 for s in scores)
    print(f"Regular: {is_regular}")

    # 2. Count Hamiltonian paths
    if n <= 16:
        H = count_ham_paths(T)
        print(f"H(T) = {H}")
    else:
        print(f"H(T): skipped (n={n} too large for exact computation)")
        H = None

    # 3. Cartan decomposition
    cartan = cartan_decompose(A)
    print(f"\nCartan decomposition of attention matrix:")
    print(f"  ||A_anti|| / ||A|| = {cartan['anti_fraction']:.4f}  "
          f"(tournament/'active' fraction)")
    print(f"  ||A_sym||  / ||A|| = {cartan['dark_fraction']:.4f}  "
          f"(similarity/'dark' fraction)")
    print(f"  Dimensions: active={cartan['dimensions']['active']}, "
          f"dark={cartan['dimensions']['dark']}, null=1")

    # 4. Odd cycles and OCF verification
    if n <= max_n_for_cycles:
        cycles = find_odd_cycles(T)
        print(f"\nOdd cycle analysis:")
        print(f"  Total directed odd cycles: {len(cycles)}")

        # Count by length
        by_length = defaultdict(int)
        for c in cycles:
            by_length[len(c)] += 1
        for length in sorted(by_length):
            print(f"    {length}-cycles: {by_length[length]}")

        if len(cycles) <= 24:
            omega = build_conflict_graph(cycles)
            I_val = independence_polynomial(omega, 2)
            print(f"  I(Omega(T), 2) = {I_val}")
            if H is not None:
                if I_val == H:
                    print(f"  OCF VERIFIED: H(T) = I(Omega(T), 2) = {H} [OK]")
                else:
                    print(f"  OCF MISMATCH: H={H} != I={I_val} [FAIL]")
        else:
            print(f"  Conflict graph too large ({len(cycles)} vertices) "
                  f"for independence polynomial")

    # 5. Additional tournament properties
    # Transitivity index = fraction of transitive triples
    trans_count = 0
    total_triples = 0
    for triple in combinations(range(n), 3):
        total_triples += 1
        a, b, c = triple
        scores_triple = [T[a][b] + T[a][c], T[b][a] + T[b][c], T[c][a] + T[c][b]]
        if 0 in scores_triple:  # one vertex beaten by both others
            trans_count += 1
    trans_index = trans_count / total_triples if total_triples > 0 else 0
    print(f"\nTransitivity index: {trans_index:.4f} "
          f"({trans_count}/{total_triples} transitive triples)")

    return {
        'T': T, 'H': H, 'scores': scores, 'cartan': cartan,
        'is_regular': is_regular, 'transitivity': trans_index
    }


def comparative_analysis(n=7, num_random=50, seed=42):
    """
    Compare Cartan decomposition across different attention pattern types.
    Tests whether tournament (antisymmetric) or similarity (symmetric)
    component dominates.
    """
    print("\n" + "="*60)
    print("COMPARATIVE CARTAN DECOMPOSITION ANALYSIS")
    print("="*60)
    print(f"n = {n}, {num_random} random attention matrices")

    results = {
        'random': [],
        'paley': None,
        'transitive': None
    }

    # Random attention matrices
    np.random.seed(seed)
    for i in range(num_random):
        A = random_attention_matrix(n, temperature=1.0, seed=seed+i)
        cartan = cartan_decompose(A)
        results['random'].append(cartan['anti_fraction'])

    # Paley attention (only if n is prime >= 3)
    if n >= 3 and all(n % i != 0 for i in range(2, int(n**0.5)+1)):
        A_paley = paley_attention(n)
        results['paley'] = cartan_decompose(A_paley)

    # Transitive attention
    A_trans = transitive_attention(n)
    results['transitive'] = cartan_decompose(A_trans)

    # Report
    random_anti = results['random']
    print(f"\nRandom attention anti-fraction:")
    print(f"  Mean:   {np.mean(random_anti):.4f}")
    print(f"  Std:    {np.std(random_anti):.4f}")
    print(f"  Min:    {np.min(random_anti):.4f}")
    print(f"  Max:    {np.max(random_anti):.4f}")

    if results['paley'] is not None:
        print(f"\nPaley T_{n} attention anti-fraction: "
              f"{results['paley']['anti_fraction']:.4f}")
        print(f"  (dark fraction: {results['paley']['dark_fraction']:.4f})")

    print(f"\nTransitive attention anti-fraction: "
          f"{results['transitive']['anti_fraction']:.4f}")
    print(f"  (dark fraction: {results['transitive']['dark_fraction']:.4f})")

    # Key comparison
    print(f"\n--- KEY FINDING ---")
    if results['paley'] is not None:
        paley_anti = results['paley']['anti_fraction']
        if paley_anti > np.mean(random_anti) + 2*np.std(random_anti):
            print(f"Paley has SIGNIFICANTLY MORE tournament structure "
                  f"than random")
        elif paley_anti < np.mean(random_anti) - 2*np.std(random_anti):
            print(f"Paley has SIGNIFICANTLY LESS tournament structure "
                  f"than random")
        else:
            print(f"Paley tournament structure is within normal range")

    return results


def killing_form_analysis(n_values=[3, 5, 7, 11]):
    """
    Verify the Killing form computation for our tournament Lie algebra
    and compare with gl(n) Killing form structure.

    Our result: K = -2(n-2)*I for so(n)
    Napolitano: gl(4,R) has Killing form with signature (9,6,1)

    General gl(n,R): Killing form K(X,Y) = 2n*tr(XY) - 2*tr(X)*tr(Y)
    Restricted to sl(n,R): signature ((n^2+n-2)/2, (n^2-n)/2)
    """
    print("\n" + "="*60)
    print("KILLING FORM COMPARISON: so(n) vs gl(n,R)")
    print("="*60)

    for n in n_values:
        dim_so = n*(n-1)//2
        dim_gl = n*n
        dim_sl = n*n - 1

        # so(n) Killing form: K = -2(n-2)*I
        # Eigenvalues: all equal to -2(n-2), so signature = (0, dim_so, 0)
        # (all negative!)
        killing_so = -2*(n-2)

        # gl(n,R) Killing form on sl(n,R):
        # signature = ((n^2+n-2)/2, (n^2-n)/2)
        pos_gl = (n*n + n - 2)//2
        neg_gl = (n*n - n)//2
        null_gl = 1  # scalar part

        # The "active" = neg_gl = n(n-1)/2 = dim(so(n))!
        # The "dark" = pos_gl = (n^2+n-2)/2 = dim(symmetric traceless) + check

        print(f"\nn = {n}:")
        print(f"  so({n}): dim = {dim_so}, "
              f"Killing eigenvalue = {killing_so}, signature = (0, {dim_so}, 0)")
        print(f"  gl({n},R): dim = {dim_gl}, "
              f"Killing signature on sl = ({pos_gl}, {neg_gl}, {null_gl})")
        print(f"  'Active' (negative) = {neg_gl} = dim(so({n})) = {dim_so}: "
              f"{'MATCH' if neg_gl == dim_so else 'MISMATCH'}")
        print(f"  'Dark' (positive) = {pos_gl}, null = {null_gl}")
        print(f"  Napolitano ratio: dark/active = {pos_gl}/{neg_gl} = "
              f"{pos_gl/neg_gl:.3f}")

        # For n=4 specifically (Napolitano's case)
        if n == 4:
            print(f"\n  *** n=4 is Napolitano's case: ***")
            print(f"  gl(4,R): 16 dims = 6 active + 9 dark + 1 null = "
                  f"{neg_gl} + {pos_gl} + {null_gl}")
            print(f"  Active = so(4) ≅ su(2) ⊕ su(2)")
            print(f"  Dark = symmetric traceless 4×4 matrices")


# ============================================================
# MAIN EXECUTION
# ============================================================

if __name__ == "__main__":
    print("Tournament Attention Analysis")
    print("Connecting tournament theory to transformer attention patterns")
    print("Based on analysis of Napolitano (2026) paper")
    print()

    # Part 1: Killing form comparison
    killing_form_analysis()

    # Part 2: Analyze specific attention patterns
    print("\n\n" + "#"*60)
    print("# TOURNAMENT ANALYSIS OF ATTENTION PATTERNS")
    print("#"*60)

    # Paley T_7 attention
    A_paley = paley_attention(7)
    r_paley = analyze_attention_tournament(A_paley, "Paley T_7 attention")

    # Transitive attention
    A_trans = transitive_attention(7)
    r_trans = analyze_attention_tournament(A_trans, "Transitive attention (n=7)")

    # Random attention
    A_rand = random_attention_matrix(7, seed=42)
    r_rand = analyze_attention_tournament(A_rand, "Random attention (n=7)")

    # Part 3: Comparative analysis
    comparative_analysis(n=7, num_random=100)

    # Part 4: Scale test
    print("\n\n" + "#"*60)
    print("# SCALING ANALYSIS: Cartan decomposition vs n")
    print("#"*60)

    for n in [3, 5, 7, 11, 13, 17]:
        np.random.seed(42)
        fracs = []
        for _ in range(50):
            A = random_attention_matrix(n)
            c = cartan_decompose(A)
            fracs.append(c['anti_fraction'])
        print(f"n={n:2d}: anti_fraction = {np.mean(fracs):.4f} ± {np.std(fracs):.4f}"
              f"  (active_dims={n*(n-1)//2}, dark_dims={(n*n+n-2)//2},"
              f" ratio={((n*n+n-2)//2)/(n*(n-1)//2):.2f})")

    print("\n\nDone. Results should be saved to "
          "05-knowledge/results/tournament_attention_analysis.out")
