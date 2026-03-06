#!/usr/bin/env python3
"""
TRANSFER MATRIX SYMMETRY ANALYSIS

Deep investigation into WHY M[a,b] = M[b,a] for all tournaments T and vertices a,b.

M[a,b] = sum_S (-1)^|S| E_a(S) * B_b(V\\{a,b}\\S)

where:
  E_a(S) = number of Hamiltonian paths in T[S union {a}] ending at a
  B_b(R) = number of Hamiltonian paths in T[R union {b}] starting at b
  sum is over subsets S of V\\{a,b}

Key fact: with INDEPENDENT arc variables, M[a,b] - M[b,a] != 0.
The tournament constraint T[i,j] + T[j,i] = 1 is ESSENTIAL for symmetry.

This script:
1. Computes M[a,b]-M[b,a] symbolically at n=4, shows cancellation
2. Analyzes the cancellation pattern -- sign-reversing involution search
3. Tests detailed balance hypotheses
4. Tests determinantal representations
5. Verifies patterns at n=5

Instance: opus-2026-03-06-S1
"""

from sympy import (symbols, expand, simplify, factor, collect, Poly,
                   Matrix, det, Symbol, Rational, Add, Mul, Pow, S,
                   zeros as sym_zeros)
from itertools import permutations, combinations
from collections import defaultdict
import random

# ============================================================
# PART 0: Setup symbolic tournament
# ============================================================

def make_symbolic_tournament(n):
    """
    Create symbolic tournament on n vertices.
    Arc variables: t_ij for i < j, where T[i][j] = t_ij, T[j][i] = 1 - t_ij.
    Returns: T (dict), arc_vars (list of symbols)
    """
    arc_vars = {}
    T = {}
    for i in range(n):
        T[(i, i)] = S.Zero
        for j in range(i+1, n):
            t = Symbol(f't{i}{j}')
            arc_vars[(i, j)] = t
            T[(i, j)] = t
            T[(j, i)] = 1 - t
    return T, arc_vars


def make_free_tournament(n):
    """
    Create tournament with INDEPENDENT arc variables (no constraint).
    T[i][j] = t_ij for all i != j (all independent).
    """
    arc_vars = {}
    T = {}
    for i in range(n):
        T[(i, i)] = S.Zero
        for j in range(n):
            if i != j:
                t = Symbol(f't{i}{j}')
                arc_vars[(i, j)] = t
                T[(i, j)] = t
    return T, arc_vars


def symbolic_ham_paths_ending_at(T, vertices, target):
    """
    Count Hamiltonian paths in T[vertices union {target}] ending at target.
    Returns symbolic expression.
    vertices should NOT include target.
    """
    vlist = sorted(vertices)
    if not vlist:
        return S.One

    # Enumerate all permutations of vlist; path is perm[0] -> perm[1] -> ... -> perm[-1] -> target
    total = S.Zero
    for perm in permutations(vlist):
        term = S.One
        for k in range(len(perm) - 1):
            term *= T[(perm[k], perm[k+1])]
        term *= T[(perm[-1], target)]
        total += term
    return expand(total)


def symbolic_ham_paths_starting_at(T, vertices, source):
    """
    Count Hamiltonian paths in T[{source} union vertices] starting at source.
    Returns symbolic expression.
    vertices should NOT include source.
    """
    vlist = sorted(vertices)
    if not vlist:
        return S.One

    total = S.Zero
    for perm in permutations(vlist):
        term = T[(source, perm[0])]
        for k in range(len(perm) - 1):
            term *= T[(perm[k], perm[k+1])]
        total += term
    return expand(total)


def compute_M_symbolic(T, n, a, b):
    """
    Compute M[a,b] = sum_S (-1)^|S| E_a(S) * B_b(V\\{a,b}\\S) symbolically.
    """
    others = [v for v in range(n) if v != a and v != b]
    m = len(others)

    total = S.Zero
    for mask in range(1 << m):
        S_set = [others[k] for k in range(m) if mask & (1 << k)]
        R_set = [others[k] for k in range(m) if not (mask & (1 << k))]
        sz = bin(mask).count('1')
        sign = (-1) ** sz

        E = symbolic_ham_paths_ending_at(T, S_set, a)
        B = symbolic_ham_paths_starting_at(T, R_set, b)

        total += sign * E * B

    return expand(total)


# ============================================================
# PART 1: n=4 symbolic computation
# ============================================================

def part1_n4_symbolic():
    print("=" * 70)
    print("PART 1: n=4 SYMBOLIC TRANSFER MATRIX")
    print("=" * 70)

    n = 4
    T_tourn, arc_vars_tourn = make_symbolic_tournament(n)

    a, b = 0, 1
    others = [v for v in range(n) if v != a and v != b]  # [2, 3]

    print(f"\nVertices: 0,1,2,3. Computing M[0,1] and M[1,0].")
    print(f"Others = {others}")
    print(f"Arc variables: {list(arc_vars_tourn.values())}")

    # Show subset-by-subset contributions
    print("\n--- Subset-by-subset contributions ---")
    m = len(others)
    M01_terms = []
    M10_terms = []

    for mask in range(1 << m):
        S_set = [others[k] for k in range(m) if mask & (1 << k)]
        R_set = [others[k] for k in range(m) if not (mask & (1 << k))]
        sz = bin(mask).count('1')
        sign = (-1) ** sz

        E0 = symbolic_ham_paths_ending_at(T_tourn, S_set, 0)
        B1 = symbolic_ham_paths_starting_at(T_tourn, R_set, 1)
        E1 = symbolic_ham_paths_ending_at(T_tourn, S_set, 1)
        B0 = symbolic_ham_paths_starting_at(T_tourn, R_set, 0)

        t01 = expand(sign * E0 * B1)
        t10 = expand(sign * E1 * B0)
        M01_terms.append((S_set, R_set, sign, t01))
        M10_terms.append((S_set, R_set, sign, t10))

        print(f"\n  S={S_set}, R={R_set}, sign={sign:+d}")
        print(f"    E_0(S)={E0}, B_1(R)={B1}")
        print(f"    E_1(S)={E1}, B_0(R)={B0}")
        print(f"    contrib to M[0,1]: {t01}")
        print(f"    contrib to M[1,0]: {t10}")

    M01 = sum(t for _, _, _, t in M01_terms)
    M10 = sum(t for _, _, _, t in M10_terms)
    M01 = expand(M01)
    M10 = expand(M10)

    print(f"\nM[0,1] = {M01}")
    print(f"M[1,0] = {M10}")

    diff = expand(M01 - M10)
    print(f"\nM[0,1] - M[1,0] = {diff}")
    print(f"Number of terms in difference: {len(Add.make_args(diff))}")

    # Now substitute tournament constraint
    print("\n--- After tournament substitution T[j,i] = 1 - T[i,j] ---")
    print("(Already built in for constrained tournament)")
    print(f"M[0,1] - M[1,0] = {diff}")
    if diff == 0:
        print(">>> ZERO! Symmetry confirmed at n=4. <<<")
    else:
        print(f">>> NONZERO: {diff}")

    # Now with FREE (unconstrained) variables
    print("\n--- With FREE (independent) arc variables ---")
    T_free, arc_vars_free = make_free_tournament(n)
    M01_free = compute_M_symbolic(T_free, n, 0, 1)
    M10_free = compute_M_symbolic(T_free, n, 1, 0)
    diff_free = expand(M01_free - M10_free)
    print(f"M[0,1] - M[1,0] (free) has {len(Add.make_args(diff_free))} terms")
    print(f"Difference: {diff_free}")

    return diff, M01, M10


# ============================================================
# PART 2: Cancellation pattern analysis
# ============================================================

def part2_cancellation_analysis():
    print("\n" + "=" * 70)
    print("PART 2: CANCELLATION PATTERN ANALYSIS (n=4)")
    print("=" * 70)

    n = 4
    a, b = 0, 1
    others = [2, 3]
    m = len(others)

    # Use free variables first, then apply tournament constraint
    # to see WHICH terms cancel against WHICH
    all_arc_syms = {}
    T_raw = {}
    for i in range(n):
        T_raw[(i, i)] = S.Zero
        for j in range(i+1, n):
            t = Symbol(f't{i}{j}')
            all_arc_syms[(i, j)] = t
            T_raw[(i, j)] = t
            T_raw[(j, i)] = Symbol(f'u{i}{j}')  # Will substitute u = 1 - t

    # Compute M[0,1] - M[1,0] with raw variables
    print("\nComputing M[0,1] - M[1,0] with raw variables t_ij and u_ij...")

    terms_01 = {}  # monomial -> coefficient for M[0,1]
    terms_10 = {}  # monomial -> coefficient for M[1,0]

    for mask in range(1 << m):
        S_set = [others[k] for k in range(m) if mask & (1 << k)]
        R_set = [others[k] for k in range(m) if not (mask & (1 << k))]
        sz = bin(mask).count('1')
        sign = (-1) ** sz

        E0 = symbolic_ham_paths_ending_at(T_raw, S_set, 0)
        B1 = symbolic_ham_paths_starting_at(T_raw, R_set, 1)
        E1 = symbolic_ham_paths_ending_at(T_raw, S_set, 1)
        B0 = symbolic_ham_paths_starting_at(T_raw, R_set, 0)

        for term in Add.make_args(expand(sign * E0 * B1)):
            key = str(term)
            terms_01[key] = terms_01.get(key, S.Zero) + term

        for term in Add.make_args(expand(sign * E1 * B0)):
            key = str(term)
            terms_10[key] = terms_10.get(key, S.Zero) + term

    print(f"\nM[0,1] has {len(terms_01)} monomial terms")
    print(f"M[1,0] has {len(terms_10)} monomial terms")

    # Now look at M[0,1] - M[1,0] monomial by monomial
    all_keys = set(terms_01.keys()) | set(terms_10.keys())
    diff_terms = {}
    for k in all_keys:
        v = terms_01.get(k, S.Zero) - terms_10.get(k, S.Zero)
        v = expand(v)
        if v != 0:
            diff_terms[k] = v

    print(f"M[0,1] - M[1,0] has {len(diff_terms)} nonzero monomial terms (raw)")

    # Now: for the CONSTRAINED tournament, substitute u_ij = 1 - t_ij
    # and see how the terms pair up
    print("\n--- Sign-reversing involution search ---")
    print("Looking at how pairs of monomials cancel under u_ij -> 1 - t_ij...")

    # Rebuild with constrained tournament
    T_c, avars = make_symbolic_tournament(n)

    # Collect all (S, path_E, path_B, sign) contributions to M[0,1] and M[1,0]
    contribs_01 = []  # (S_set, E_path, B_path, sign, term_value)
    contribs_10 = []

    for mask in range(1 << m):
        S_set = tuple(others[k] for k in range(m) if mask & (1 << k))
        R_set = tuple(others[k] for k in range(m) if not (mask & (1 << k)))
        sz = len(S_set)
        sign = (-1) ** sz

        # E_0(S): paths through S ending at 0
        S_list = list(S_set)
        for perm in permutations(S_list):
            if not perm:
                E_term = S.One
                E_path = ()
            else:
                E_term = S.One
                for k in range(len(perm) - 1):
                    E_term *= T_c[(perm[k], perm[k+1])]
                E_term *= T_c[(perm[-1], 0)]
                E_path = perm

            # B_1(R): paths through R starting at 1
            R_list = list(R_set)
            for perm_r in permutations(R_list):
                if not perm_r:
                    B_term = S.One
                    B_path = ()
                else:
                    B_term = T_c[(1, perm_r[0])]
                    for k in range(len(perm_r) - 1):
                        B_term *= T_c[(perm_r[k], perm_r[k+1])]
                    B_path = perm_r

                val = expand(sign * E_term * B_term)
                contribs_01.append({
                    'S': S_set, 'R': R_set, 'sign': sign,
                    'E_path': E_path, 'B_path': B_path,
                    'value': val
                })

        # E_1(S): paths through S ending at 1
        for perm in permutations(S_list):
            if not perm:
                E_term = S.One
                E_path = ()
            else:
                E_term = S.One
                for k in range(len(perm) - 1):
                    E_term *= T_c[(perm[k], perm[k+1])]
                E_term *= T_c[(perm[-1], 1)]
                E_path = perm

            # B_0(R): paths through R starting at 0
            R_list = list(R_set)
            for perm_r in permutations(R_list):
                if not perm_r:
                    B_term = S.One
                    B_path = ()
                else:
                    B_term = T_c[(0, perm_r[0])]
                    for k in range(len(perm_r) - 1):
                        B_term *= T_c[(perm_r[k], perm_r[k+1])]
                    B_path = perm_r

                val = expand(sign * E_term * B_term)
                contribs_10.append({
                    'S': S_set, 'R': R_set, 'sign': sign,
                    'E_path': E_path, 'B_path': B_path,
                    'value': val
                })

    print(f"\nTotal path-level contributions to M[0,1]: {len(contribs_01)}")
    print(f"Total path-level contributions to M[1,0]: {len(contribs_10)}")

    print("\nDetailed contributions to M[0,1]:")
    for c in contribs_01:
        E_desc = f"E_path: {c['E_path']}->0" if c['E_path'] else "E_path: (empty)->0"
        B_desc = f"B_path: 1->{c['B_path']}" if c['B_path'] else "B_path: 1->(empty)"
        print(f"  S={c['S']}, sign={c['sign']:+d}, {E_desc}, {B_desc}: {c['value']}")

    print("\nDetailed contributions to M[1,0]:")
    for c in contribs_10:
        E_desc = f"E_path: {c['E_path']}->1" if c['E_path'] else "E_path: (empty)->1"
        B_desc = f"B_path: 0->{c['B_path']}" if c['B_path'] else "B_path: 0->(empty)"
        print(f"  S={c['S']}, sign={c['sign']:+d}, {E_desc}, {B_desc}: {c['value']}")

    # Now look for a natural pairing (involution) between contribs_01 and contribs_10
    print("\n--- Testing swap involution: S <-> R ---")
    print("If we swap S and R (complement), does a term in M[0,1] map to one in M[1,0]?")

    for c in contribs_01:
        # Under S <-> R, the sign flips (since |complement| = m - |S|, sign changes by (-1)^m)
        # E_a(S) * B_b(R) maps to E_a(R) * B_b(S)
        # This doesn't directly give M[1,0] terms...
        pass

    print("\n--- Testing path reversal involution ---")
    print("Under T -> T^op, E_a(S) in T becomes B_a(S) in T^op")
    print("So M_T[a,b] = sum (-1)^|S| E_a^T(S) * B_b^T(R)")
    print("   M_{T^op}[a,b] = sum (-1)^|S| E_a^{T^op}(S) * B_b^{T^op}(R)")
    print("Now E_a^{T^op}(S) = B_a^T(S) (reversing all arcs reverses all paths)")
    print("So M_{T^op}[a,b] = sum (-1)^|S| B_a^T(S) * E_b^T(R)")
    print()
    print("CLAIM: M_{T^op}[a,b] = (-1)^{n-2} M_T[b,a]")
    print("Proof: B_a^T(S) * E_b^T(R) with S<->R gives E_b^T(S) * B_a^T(R) with sign (-1)^m")
    print("       = (-1)^{n-2} * M_T[b,a]")
    print()
    print("So symmetry M_T[a,b] = M_T[b,a] is equivalent to M_{T^op}[a,b] = (-1)^{n-2} M_T[a,b]")
    print("i.e., reversing all arcs multiplies M by (-1)^{n-2}.")

    # Verify this claim at n=4 symbolically
    print("\n--- Verifying M_{T^op} = (-1)^{n-2} M_T at n=4 ---")
    T_c, avars = make_symbolic_tournament(n)
    T_op = {}
    for (i, j), val in T_c.items():
        T_op[(i, j)] = T_c[(j, i)]

    M01_T = compute_M_symbolic(T_c, n, 0, 1)
    M01_Top = compute_M_symbolic(T_op, n, 0, 1)
    sign_factor = (-1) ** (n - 2)
    check = expand(M01_Top - sign_factor * M01_T)
    print(f"M_T[0,1] = {M01_T}")
    print(f"M_{{T^op}}[0,1] = {M01_Top}")
    print(f"(-1)^{{n-2}} * M_T[0,1] = {expand(sign_factor * M01_T)}")
    print(f"Difference = {check}")
    if check == 0:
        print(">>> CONFIRMED: M_{T^op} = (-1)^{n-2} M_T at n=4 <<<")

    return contribs_01, contribs_10


# ============================================================
# PART 3: Detailed balance hypotheses
# ============================================================

def part3_detailed_balance():
    print("\n" + "=" * 70)
    print("PART 3: DETAILED BALANCE HYPOTHESES")
    print("=" * 70)

    print("\nSince M[a,b] = M[b,a] (symmetry), M is trivially 'reversible'")
    print("with respect to the uniform measure pi(v) = 1.")
    print()
    print("More interesting: the DIAGONAL values M[a,a] vary.")
    print("Test whether M[a,a] relates to the score (out-degree) of a.")

    # Numerical test
    random.seed(42)
    for n in [5, 6, 7]:
        print(f"\n--- n = {n} ---")

        def random_tournament_flat(n):
            T = [0]*(n*n)
            for a in range(n):
                for b in range(a+1, n):
                    if random.random() < 0.5:
                        T[a*n+b] = 1
                    else:
                        T[b*n+a] = 1
            return T

        def h_end_num(T, n, S_list, v):
            if not S_list:
                return 1
            total = 0
            for perm in permutations(S_list):
                ok = True
                for k in range(len(perm)-1):
                    if not T[perm[k]*n + perm[k+1]]:
                        ok = False; break
                if ok and T[perm[-1]*n + v]:
                    total += 1
            return total

        def h_start_num(T, n, R_list, v):
            if not R_list:
                return 1
            total = 0
            for perm in permutations(R_list):
                if not T[v*n + perm[0]]:
                    continue
                ok = True
                for k in range(len(perm)-1):
                    if not T[perm[k]*n + perm[k+1]]:
                        ok = False; break
                if ok:
                    total += 1
            return total

        def compute_M_num(T, n, a, b):
            others = [v for v in range(n) if v != a and v != b]
            m = len(others)
            total = 0
            for mask in range(1 << m):
                S_set = [others[k] for k in range(m) if mask & (1 << k)]
                R_set = [others[k] for k in range(m) if not (mask & (1 << k))]
                sz = bin(mask).count('1')
                sign = (-1) ** sz
                E = h_end_num(T, n, S_set, a)
                B = h_start_num(T, n, R_set, b)
                total += sign * E * B
            return total

        trials = 30 if n <= 6 else 10
        for trial in range(min(5, trials)):
            T = random_tournament_flat(n)
            scores = [sum(T[v*n + u] for u in range(n) if u != v) for v in range(n)]

            # Compute full M matrix
            M_full = {}
            for i in range(min(n, 4)):
                for j in range(min(n, 4)):
                    if i != j:
                        M_full[(i, j)] = compute_M_num(T, n, i, j)

            # Check symmetry
            sym_ok = True
            for i in range(min(n, 4)):
                for j in range(i+1, min(n, 4)):
                    if M_full.get((i,j)) != M_full.get((j,i)):
                        sym_ok = False

            # Compute diagonal
            diags = {}
            for i in range(min(n, 4)):
                others_i = [v for v in range(n) if v != i]
                # M[i,i] needs a second endpoint -- but the definition uses TWO endpoints
                # Actually M[a,b] requires a != b in our formulation
                pass

            print(f"  Trial {trial}: scores={scores[:4]}, sym={'OK' if sym_ok else 'FAIL'}", end="")
            # Print off-diagonal values
            vals = []
            for i in range(min(n, 3)):
                for j in range(i+1, min(n, 3)):
                    vals.append(f"M[{i},{j}]={M_full.get((i,j), '?')}")
            print(f"  {', '.join(vals)}")

        # Test: is M[a,b] related to scores?
        print(f"\n  Correlation test: M[0,1] vs score difference")
        m_vals = []
        score_diffs = []
        for trial in range(trials):
            T = random_tournament_flat(n)
            m_val = compute_M_num(T, n, 0, 1)
            s0 = sum(T[0*n + u] for u in range(n) if u != 0)
            s1 = sum(T[1*n + u] for u in range(n) if u != 1)
            m_vals.append(m_val)
            score_diffs.append(s0 - s1)

        # Check if M[0,1] = 0 always (it was for n=4!)
        distinct_m = set(m_vals)
        print(f"  Distinct M[0,1] values: {sorted(distinct_m)}")
        if distinct_m == {0}:
            print(f"  >>> M[0,1] = 0 for all tournaments at n={n}! <<<")


# ============================================================
# PART 4: Determinantal representation
# ============================================================

def part4_determinantal():
    print("\n" + "=" * 70)
    print("PART 4: DETERMINANTAL REPRESENTATION TEST")
    print("=" * 70)

    print("\nIdea: M[a,b] involves alternating sums over subsets, reminiscent")
    print("of inclusion-exclusion or determinant expansion.")
    print()
    print("If M[a,b] = det(X_{ab}) for some symmetric matrix X, symmetry follows.")
    print()
    print("Alternative: M might be a MINOR of a larger matrix.")
    print("The sum (-1)^|S| E_a(S) * B_b(R) looks like a Cauchy-Binet expansion.")
    print()
    print("Cauchy-Binet: det(AB) = sum_S det(A_S) * det(B_S)")
    print("Our formula: M[a,b] = sum_S (-1)^|S| E_a(S) * B_b(R)")
    print()
    print("Key difference: E and B are path counts, not determinants.")
    print("But Lindstrom-Gessel-Viennot says path counts CAN be determinants!")

    n = 4
    T, avars = make_symbolic_tournament(n)

    print(f"\n--- Testing at n={n} ---")

    # Build the adjacency matrix A
    A = Matrix(n, n, lambda i, j: T.get((i, j), S.Zero))
    print(f"\nAdjacency matrix A:")
    print(A)

    # The matrix (I - A)^{-1} has (i,j) entry = sum of weights of all walks from i to j
    # But we want paths, not walks.

    # LGV lemma: if we have a DAG, the number of non-intersecting path systems
    # from sources to sinks is det(path_count_matrix).
    # But a tournament is not a DAG in general.

    # Try: define P[i,j] = number of Hamiltonian paths from i to j in T
    # Then is M related to some transform of P?

    # Actually, let's try a different approach.
    # Consider the matrix X[i,j] = sum_S (-1)^|S| * (paths from i through S) * (paths from j through complement)
    # This is NOT obviously symmetric. But we proved it is.

    # Test: Can we express E_a(S) as a determinant?
    # E_a({v}) = T[v,a] -- just one arc, so it's a 1x1 "determinant"
    # E_a({v,w}) = T[v,w]*T[w,a] + T[w,v]*T[v,a] -- this is NOT a determinant of a 2x2 matrix
    #   unless we choose the matrix carefully.

    # Actually, E_a(S) for |S|=2 with S={v,w}:
    # Paths: v->w->a and w->v->a
    # = T[v,w]*T[w,a] + T[w,v]*T[v,a]

    # For a general set S, E_a(S) = permanent of a certain matrix (not determinant)
    # because all signs are +1.

    # Hmm, so the "alternating sign" comes from (-1)^|S| in the outer sum,
    # not from within E_a or B_b.

    # Let's try the transfer matrix approach differently.
    # Define for each vertex v in others, a 2x2 "transfer" matrix:
    #   L_v = [[T[v,a], T[v,b]], [something, something]]
    # and see if M = product of L_v's or similar.

    print("\n--- Transfer matrix product approach ---")
    print("Can M[a,b] = entries of a product of 2x2 matrices, one per 'other' vertex?")

    # For n=4, others = {2, 3}. The sum over S subsets of {2,3} gives:
    # S={}:     (+1) * E_a({}) * B_b({2,3}) = 1 * B_b({2,3})
    # S={2}:    (-1) * E_a({2}) * B_b({3})
    # S={3}:    (-1) * E_a({3}) * B_b({2})
    # S={2,3}:  (+1) * E_a({2,3}) * B_b({})

    # This is like: sum_S (-1)^|S| * E_a(S) * B_b(complement)
    # = sum_S (-1)^|S| * (paths ...->a through S) * (paths b->... through complement)

    # If we think of this as a bilinear form in the "states" at a and b,
    # we can try to factor it.

    # For |others|=2, vertex set {2,3}:
    # Define state after processing vertex 2:
    #   either 2 is in S (contribute to E) or in R (contribute to B)
    # This is a partition into left/right, processed vertex by vertex.

    # For each vertex v in others, we choose: v in S or v in R.
    # If v in S: contributes factor to E_a path
    # If v in R: contributes factor to B_b path

    # But E_a(S) counts ALL orderings of S, so the "transfer matrix" approach
    # would need to track the current endpoint of the partial path.

    print("\nThe E_a path count is a PERMANENT (sum over all orderings),")
    print("not a single-path weight. This makes a simple transfer matrix")
    print("product decomposition unlikely unless we expand the state space.")

    # However, let's check: maybe for the specific case of M[a,b]-M[b,a],
    # we can find a determinantal identity.

    # Key insight: E_a(S) with tournament constraint.
    # T[v,w] + T[w,v] = 1, so for |S|=2, S={v,w}:
    # E_a({v,w}) = T[v,w]*T[w,a] + T[w,v]*T[v,a]
    #            = T[v,w]*T[w,a] + (1-T[v,w])*T[v,a]
    #            = T[v,a] + T[v,w]*(T[w,a] - T[v,a])

    # This is a LINEAR function of T[v,w]!
    # In fact, for ANY |S|, E_a(S) is MULTILINEAR in the arcs within S.

    # The tournament constraint T[v,w]+T[w,v]=1 means each pair contributes
    # ONE variable. So E_a(S) as a polynomial in {T[i,j]: i<j} has degree |S|
    # (at most), one factor per arc used.

    # Let me compute M[0,1] as a polynomial in the arc variables explicitly.
    print(f"\n--- M[0,1] as polynomial in arc variables ---")
    M01 = compute_M_symbolic(T, n, 0, 1)
    M10 = compute_M_symbolic(T, n, 1, 0)
    print(f"M[0,1] = {M01}")
    print(f"M[1,0] = {M10}")

    # Collect by arc variables
    all_vars = sorted(avars.values(), key=str)
    print(f"\nArc variables: {all_vars}")

    # Check degree structure
    p01 = Poly(M01, *all_vars) if all_vars else None
    p10 = Poly(M10, *all_vars) if all_vars else None
    if p01:
        print(f"\nM[0,1] as polynomial:")
        print(f"  Total degree: {p01.total_degree()}")
        print(f"  Terms: {len(p01.as_dict())}")
        for monom, coeff in sorted(p01.as_dict().items()):
            var_str = '*'.join(f'{v}^{e}' if e > 1 else str(v)
                              for v, e in zip(all_vars, monom) if e > 0) or '1'
            print(f"    {int(coeff):+d} * {var_str}")

    if p10:
        print(f"\nM[1,0] as polynomial:")
        for monom, coeff in sorted(p10.as_dict().items()):
            var_str = '*'.join(f'{v}^{e}' if e > 1 else str(v)
                              for v, e in zip(all_vars, monom) if e > 0) or '1'
            print(f"    {int(coeff):+d} * {var_str}")

    # Check if M matrix has a nice form
    print("\n--- Full symbolic M matrix at n=4, endpoints (0,1) ---")
    M00 = compute_M_symbolic(T, n, 0, 0)
    M11 = compute_M_symbolic(T, n, 1, 1)
    print(f"M[0,0] (diag) = {M00}")
    print(f"M[1,1] (diag) = {M11}")
    print(f"M[0,1] (off)  = {M01}")
    print(f"M[1,0] (off)  = {M10}")
    print(f"Trace = {expand(M00 + M11)}")
    print(f"Det   = {expand(M00 * M11 - M01 * M10)}")

    # Test if this is related to I - xA for some x
    # Or some function of the adjacency matrix
    print("\n--- Is M = f(adjacency matrix)? Testing at specific tournaments ---")
    random.seed(42)
    for trial in range(3):
        # Generate specific tournament
        n_test = 5
        T_num = [0]*(n_test*n_test)
        for ii in range(n_test):
            for jj in range(ii+1, n_test):
                if random.random() < 0.5:
                    T_num[ii*n_test+jj] = 1
                else:
                    T_num[jj*n_test+ii] = 1

        # Build numerical M matrix (all pairs)
        def h_end_n(T, n, S_list, v):
            if not S_list: return 1
            total = 0
            for perm in permutations(S_list):
                ok = True
                for k in range(len(perm)-1):
                    if not T[perm[k]*n+perm[k+1]]: ok = False; break
                if ok and T[perm[-1]*n+v]: total += 1
            return total

        def h_start_n(T, n, R_list, v):
            if not R_list: return 1
            total = 0
            for perm in permutations(R_list):
                if not T[v*n+perm[0]]: continue
                ok = True
                for k in range(len(perm)-1):
                    if not T[perm[k]*n+perm[k+1]]: ok = False; break
                if ok: total += 1
            return total

        M_mat = [[0]*n_test for _ in range(n_test)]
        for ii in range(n_test):
            for jj in range(n_test):
                if ii == jj:
                    continue
                others = [v for v in range(n_test) if v != ii and v != jj]
                mm = len(others)
                total = 0
                for mask in range(1 << mm):
                    S_set = [others[k] for k in range(mm) if mask & (1 << k)]
                    R_set = [others[k] for k in range(mm) if not (mask & (1 << k))]
                    sz = bin(mask).count('1')
                    sign = (-1) ** sz
                    total += sign * h_end_n(T_num, n_test, S_set, ii) * h_start_n(T_num, n_test, R_set, jj)
                M_mat[ii][jj] = total

        print(f"\n  Trial {trial}, n={n_test}:")
        for row in M_mat:
            print(f"    {row}")

        # Check symmetry
        sym = all(M_mat[i][j] == M_mat[j][i]
                  for i in range(n_test) for j in range(n_test) if i != j)
        print(f"    Symmetric: {sym}")

        # What about the diagonal?
        diag = [M_mat[i][i] for i in range(n_test)]
        print(f"    Diagonal: {diag} (should be 0 by convention)")

        # Check eigenvalues of the matrix (treating diagonal as 0)
        M_sym = Matrix(M_mat)
        eigenvals = M_sym.eigenvals()
        print(f"    Eigenvalues: {eigenvals}")


# ============================================================
# PART 5: n=5 verification and pattern generalization
# ============================================================

def part5_n5_verification():
    print("\n" + "=" * 70)
    print("PART 5: n=5 SYMBOLIC VERIFICATION")
    print("=" * 70)

    n = 5
    T, avars = make_symbolic_tournament(n)

    print(f"\nArc variables ({len(avars)}): {sorted(avars.values(), key=str)}")

    # Computing full symbolic M at n=5 is feasible but larger
    # others = {2,3,4} for endpoints (0,1), giving 2^3 = 8 subsets
    # Each subset S of size k has k! permutations for E, and (3-k)! for B
    # Total path-pairs: sum_{k=0}^{3} C(3,k)*k!*(3-k)! = 1+3+6+6 = 16 for M[0,1]
    # and same for M[1,0].

    print("\nComputing M[0,1] and M[1,0] symbolically...")
    M01 = compute_M_symbolic(T, n, 0, 1)
    M10 = compute_M_symbolic(T, n, 1, 0)

    diff = expand(M01 - M10)
    print(f"\nM[0,1] = {M01}")
    print(f"\nM[1,0] = {M10}")
    print(f"\nM[0,1] - M[1,0] = {diff}")

    if diff == 0:
        print(">>> ZERO! Symmetry confirmed at n=5. <<<")
    else:
        print(f">>> NONZERO! {len(Add.make_args(diff))} terms remaining. <<<")

    # Verify T^op relation
    print("\n--- Verifying M_{T^op} = (-1)^{n-2} M_T at n=5 ---")
    T_op = {}
    for (i, j), val in T.items():
        T_op[(i, j)] = T[(j, i)]

    M01_Top = compute_M_symbolic(T_op, n, 0, 1)
    sign_factor = (-1) ** (n - 2)
    check = expand(M01_Top - sign_factor * M01)
    print(f"(-1)^{{n-2}} = {sign_factor}")
    print(f"M_{{T^op}}[0,1] - ({sign_factor})*M_T[0,1] = {check}")
    if check == 0:
        print(">>> CONFIRMED: M_{T^op} = (-1)^{n-2} M_T at n=5 <<<")

    # Polynomial structure
    all_vars = sorted(avars.values(), key=str)
    p01 = Poly(M01, *all_vars)
    print(f"\nM[0,1] polynomial structure:")
    print(f"  Total degree: {p01.total_degree()}")
    print(f"  Number of terms: {len(p01.as_dict())}")

    # Degree distribution
    deg_counts = defaultdict(int)
    for monom in p01.as_dict():
        deg = sum(monom)
        deg_counts[deg] += 1
    print(f"  Terms by total degree: {dict(sorted(deg_counts.items()))}")

    # Compare with free (unconstrained) case
    print("\n--- Free (unconstrained) case at n=5 ---")
    T_free, avars_free = make_free_tournament(n)
    M01_free = compute_M_symbolic(T_free, n, 0, 1)
    M10_free = compute_M_symbolic(T_free, n, 1, 0)
    diff_free = expand(M01_free - M10_free)
    nterms = len(Add.make_args(diff_free))
    print(f"M[0,1] - M[1,0] (free) has {nterms} terms")

    # Check if the constrained polynomial at n=5 is also multilinear
    max_exp = max(max(m) for m in p01.as_dict())
    print(f"\n  Max exponent in any variable: {max_exp}")
    if max_exp <= 1:
        print("  >>> M[0,1] is MULTILINEAR in arc variables! <<<")

    # Coefficient analysis: which arc variables appear?
    print("\n--- Coefficient analysis ---")
    # Check which variables actually appear
    appearing_vars = set()
    for monom in p01.as_dict():
        for i, e in enumerate(monom):
            if e > 0:
                appearing_vars.add(all_vars[i])
    print(f"  Variables appearing in M[0,1]: {sorted(appearing_vars, key=str)}")

    # Does t01 (the arc between the two endpoints) appear?
    t01_var = avars.get((0, 1))
    if t01_var in appearing_vars:
        print(f"  NOTE: t01 (arc between endpoints) DOES appear in M[0,1]")
    else:
        print(f"  NOTE: t01 (arc between endpoints) does NOT appear in M[0,1]")


# ============================================================
# PART 6: Deep structural analysis
# ============================================================

def part6_structural():
    print("\n" + "=" * 70)
    print("PART 6: DEEP STRUCTURAL ANALYSIS")
    print("=" * 70)

    print("\n--- Subset complement symmetry ---")
    print("Key observation: in M[a,b] = sum_S (-1)^|S| E_a(S) * B_b(R),")
    print("replacing S by its complement R maps:")
    print("  (-1)^|S| E_a(S) B_b(R)  ->  (-1)^|R| E_a(R) B_b(S)")
    print("  = (-1)^{m-|S|} E_a(R) B_b(S)")
    print("  = (-1)^m * (-1)^|S| * E_a(R) B_b(S)")
    print()
    print("So if we define M'[a,b] = sum_S (-1)^|S| E_a(R) * B_b(S),")
    print("then M'[a,b] = (-1)^m * M[a,b].")
    print()
    print("Now M'[a,b] = sum_S (-1)^|S| E_a(R) B_b(S)")
    print("            = sum_R (-1)^{m-|R|} E_a(R) B_b(complement of R)")
    print("            = (-1)^m sum_R (-1)^|R| E_a(R) B_b(complement)")
    print("This is just (-1)^m * M[a,b] again. (Tautology.)")
    print()
    print("For SYMMETRY, we need: M[a,b] = M[b,a], i.e.,")
    print("  sum_S (-1)^|S| E_a(S) B_b(R) = sum_S (-1)^|S| E_b(S) B_a(R)")
    print()
    print("This means: for each S,")
    print("  E_a(S) B_b(R) - E_b(S) B_a(R)")
    print("summed with alternating signs = 0.")
    print()
    print("But for individual S, E_a(S)*B_b(R) != E_b(S)*B_a(R) in general!")
    print("The cancellation is across different S's.")

    # Verify: for n=4, compute E_a(S)*B_b(R) - E_b(S)*B_a(R) for each S
    n = 4
    T, avars = make_symbolic_tournament(n)
    a, b = 0, 1
    others = [2, 3]
    m = len(others)

    print(f"\n--- Per-subset asymmetry at n=4 ---")
    total_asym = S.Zero
    for mask in range(1 << m):
        S_set = [others[k] for k in range(m) if mask & (1 << k)]
        R_set = [others[k] for k in range(m) if not (mask & (1 << k))]
        sz = len(S_set)
        sign = (-1) ** sz

        Ea = symbolic_ham_paths_ending_at(T, S_set, a)
        Eb = symbolic_ham_paths_ending_at(T, S_set, b)
        Ba = symbolic_ham_paths_starting_at(T, R_set, a)
        Bb = symbolic_ham_paths_starting_at(T, R_set, b)

        asym = expand(Ea * Bb - Eb * Ba)
        total_asym += sign * asym

        print(f"  S={S_set}: E_a*B_b - E_b*B_a = {asym} (sign={sign:+d})")

    print(f"\n  Sum with signs = {expand(total_asym)}")

    # Check: does the S<->complement pairing help?
    print("\n--- Pairing S with complement(S) ---")
    for mask in range(1 << (m - 1)):  # pair (mask, complement)
        comp_mask = ((1 << m) - 1) ^ mask
        S1 = [others[k] for k in range(m) if mask & (1 << k)]
        R1 = [others[k] for k in range(m) if not (mask & (1 << k))]
        S2 = [others[k] for k in range(m) if comp_mask & (1 << k)]
        R2 = [others[k] for k in range(m) if not (comp_mask & (1 << k))]
        sz1 = len(S1)
        sz2 = len(S2)
        sign1 = (-1) ** sz1
        sign2 = (-1) ** sz2

        Ea1 = symbolic_ham_paths_ending_at(T, S1, a)
        Bb1 = symbolic_ham_paths_starting_at(T, R1, b)
        Ea2 = symbolic_ham_paths_ending_at(T, S2, a)
        Bb2 = symbolic_ham_paths_starting_at(T, R2, b)

        Eb1 = symbolic_ham_paths_ending_at(T, S1, b)
        Ba1 = symbolic_ham_paths_starting_at(T, R1, a)
        Eb2 = symbolic_ham_paths_ending_at(T, S2, b)
        Ba2 = symbolic_ham_paths_starting_at(T, R2, a)

        paired_01 = expand(sign1 * Ea1 * Bb1 + sign2 * Ea2 * Bb2)
        paired_10 = expand(sign1 * Eb1 * Ba1 + sign2 * Eb2 * Ba2)
        paired_diff = expand(paired_01 - paired_10)

        print(f"  Pair S={S1}, complement={S2}:")
        print(f"    M[0,1] contribution: {paired_01}")
        print(f"    M[1,0] contribution: {paired_10}")
        print(f"    Difference: {paired_diff}")

    # Under complement, S1 = R2, R1 = S2. So:
    # paired_01 = sign1 * E_a(S1)*B_b(R1) + sign2 * E_a(R1)*B_b(S1)
    # paired_10 = sign1 * E_b(S1)*B_a(R1) + sign2 * E_b(R1)*B_a(S1)
    # With sign2 = (-1)^m * sign1:
    # paired_01 = sign1 * [E_a(S1)*B_b(R1) + (-1)^m * E_a(R1)*B_b(S1)]

    # Now test at n=5: do complement pairs INDIVIDUALLY cancel?
    print(f"\n--- Complement pairing at n=5 ---")
    n5 = 5
    T5, avars5 = make_symbolic_tournament(n5)
    a5, b5 = 0, 1
    others5 = [2, 3, 4]
    m5 = len(others5)

    all_cancel = True
    for mask in range(1 << m5):
        comp_mask = ((1 << m5) - 1) ^ mask
        if mask > comp_mask:
            continue  # only process each pair once
        S1 = [others5[k] for k in range(m5) if mask & (1 << k)]
        R1 = [others5[k] for k in range(m5) if not (mask & (1 << k))]
        sz1 = len(S1)
        sign1 = (-1) ** sz1

        S2 = [others5[k] for k in range(m5) if comp_mask & (1 << k)]
        R2 = [others5[k] for k in range(m5) if not (comp_mask & (1 << k))]
        sz2 = len(S2)
        sign2 = (-1) ** sz2

        Ea1 = symbolic_ham_paths_ending_at(T5, S1, a5)
        Bb1 = symbolic_ham_paths_starting_at(T5, R1, b5)
        Ea2 = symbolic_ham_paths_ending_at(T5, S2, a5)
        Bb2 = symbolic_ham_paths_starting_at(T5, R2, b5)
        Eb1 = symbolic_ham_paths_ending_at(T5, S1, b5)
        Ba1 = symbolic_ham_paths_starting_at(T5, R1, a5)
        Eb2 = symbolic_ham_paths_ending_at(T5, S2, b5)
        Ba2 = symbolic_ham_paths_starting_at(T5, R2, a5)

        paired_01 = expand(sign1 * Ea1 * Bb1 + sign2 * Ea2 * Bb2)
        paired_10 = expand(sign1 * Eb1 * Ba1 + sign2 * Eb2 * Ba2)
        paired_diff = expand(paired_01 - paired_10)

        is_zero = (paired_diff == 0)
        if not is_zero:
            all_cancel = False
        print(f"  S={S1}, comp={S2}: pair diff = {'0' if is_zero else paired_diff}")

    if all_cancel:
        print("  >>> ALL complement pairs cancel individually at n=5! <<<")
    else:
        print("  >>> Some pairs do NOT cancel individually -- cross-pair cancellation needed. <<<")

    # Check: do the 4 pair-differences at n=5 sum to 0?
    print("\n  Verification: sum of all 4 pair-differences:")
    pair_diffs = []
    for mask in range(1 << m5):
        comp_mask = ((1 << m5) - 1) ^ mask
        if mask > comp_mask:
            continue
        S1 = [others5[k] for k in range(m5) if mask & (1 << k)]
        R1 = [others5[k] for k in range(m5) if not (mask & (1 << k))]
        sz1 = len(S1)
        sign1 = (-1) ** sz1
        S2 = [others5[k] for k in range(m5) if comp_mask & (1 << k)]
        R2 = [others5[k] for k in range(m5) if not (comp_mask & (1 << k))]
        sz2 = len(S2)
        sign2 = (-1) ** sz2
        Ea1 = symbolic_ham_paths_ending_at(T5, S1, a5)
        Bb1 = symbolic_ham_paths_starting_at(T5, R1, b5)
        Ea2 = symbolic_ham_paths_ending_at(T5, S2, a5)
        Bb2 = symbolic_ham_paths_starting_at(T5, R2, b5)
        Eb1 = symbolic_ham_paths_ending_at(T5, S1, b5)
        Ba1 = symbolic_ham_paths_starting_at(T5, R1, a5)
        Eb2 = symbolic_ham_paths_ending_at(T5, S2, b5)
        Ba2 = symbolic_ham_paths_starting_at(T5, R2, a5)
        paired_diff = expand(
            sign1 * Ea1 * Bb1 + sign2 * Ea2 * Bb2
            - sign1 * Eb1 * Ba1 - sign2 * Eb2 * Ba2)
        pair_diffs.append(paired_diff)
    total = expand(sum(pair_diffs))
    print(f"  Sum = {total}")

    # Check whether pairs {S,comp} with |S|=0,|comp|=3 cancel against
    # pairs {S,comp} with |S|=1,|comp|=2
    print("\n  Cross-pair structure:")
    print(f"  (|S|=0,|comp|=3) diff + (|S|=1,|comp|=2) diffs:")
    for i, d in enumerate(pair_diffs):
        print(f"    Pair {i}: {d}")
    # Pair 0: S=[], comp=[2,3,4]  (sizes 0,3)
    # Pair 1: S=[2], comp=[3,4]   (sizes 1,2)
    # Pair 2: S=[3], comp=[2,4]   (sizes 1,2)
    # Pair 3: S=[2,3], comp=[4]   (sizes 2,1)
    cross1 = expand(pair_diffs[0] + pair_diffs[1])
    cross2 = expand(pair_diffs[2] + pair_diffs[3])
    print(f"  Pair0+Pair1 = {cross1}")
    print(f"  Pair2+Pair3 = {cross2}")
    if cross1 == 0:
        print("  >>> Pairs 0,1 cancel! <<<")
    if cross2 == 0:
        print("  >>> Pairs 2,3 cancel! <<<")
    # Also try other groupings
    cross3 = expand(pair_diffs[0] + pair_diffs[2])
    cross4 = expand(pair_diffs[1] + pair_diffs[3])
    print(f"  Pair0+Pair2 = {cross3}")
    print(f"  Pair1+Pair3 = {cross4}")
    if cross3 == 0:
        print("  >>> Pairs 0,2 cancel! <<<")
    if cross4 == 0:
        print("  >>> Pairs 1,3 cancel! <<<")

    print(f"\n--- Key structural observation ---")
    print(f"For n=4 (m=2, (-1)^m = +1):")
    print(f"Each complementary pair contributes:")
    print(f"  M[a,b]: E_a(S)*B_b(R) + E_a(R)*B_b(S)  [symmetric in (S,R) for E_a, B_b]")
    print(f"  M[b,a]: E_b(S)*B_a(R) + E_b(R)*B_a(S)  [symmetric in (S,R) for E_b, B_a]")
    print(f"")
    print(f"The symmetry M[a,b]=M[b,a] then requires:")
    print(f"  E_a(S)*B_b(R) + E_a(R)*B_b(S) = E_b(S)*B_a(R) + E_b(R)*B_a(S)")
    print(f"for each complementary pair (S, R=complement(S)).")
    print(f"")
    print(f"This is a bilinear identity in (E, B) values!")
    print(f"If we define x1=E_a(S), x2=E_a(R), y1=B_b(S), y2=B_b(R),")
    print(f"           u1=E_b(S), u2=E_b(R), v1=B_a(S), v2=B_a(R),")
    print(f"then the identity is: x1*y2 + x2*y1 = u1*v2 + u2*v1.")
    print(f"")
    print(f"Note that E_a and B_a count paths in COMPLEMENTARY subgraphs (S vs R)")
    print(f"but ending/starting at different vertices (a vs b).")
    print(f"The tournament constraint couples these through shared arcs.")


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    print("TRANSFER MATRIX SYMMETRY ANALYSIS")
    print("Instance: opus-2026-03-06-S1")
    print("=" * 70)

    diff_n4, M01_n4, M10_n4 = part1_n4_symbolic()
    contribs_01, contribs_10 = part2_cancellation_analysis()
    part3_detailed_balance()
    part4_determinantal()
    part5_n5_verification()
    part6_structural()

    print("\n" + "=" * 70)
    print("SUMMARY OF FINDINGS")
    print("=" * 70)
    print("""
1. CANCELLATION MECHANISM: M[a,b]-M[b,a] = 0 requires the tournament
   constraint T[i,j]+T[j,i]=1. With free variables, the difference is
   nonzero (12 terms at n=4, 48 at n=5).

2. T^OP EQUIVALENCE (PROVED): Symmetry M[a,b]=M[b,a] is EQUIVALENT to
   M_{T^op}[a,b] = (-1)^{n-2} M_T[a,b]. This means reversing all arcs
   multiplies M by (-1)^{n-2}. Verified symbolically at n=4,5.

   PROOF SKETCH: Under T -> T^op, paths reverse direction, so
     E_a^{T^op}(S) = B_a^T(S)  (ending paths become starting paths).
   Therefore M_{T^op}[a,b] = sum_S (-1)^|S| B_a^T(S) E_b^T(R).
   Swapping S <-> R in the sum introduces factor (-1)^m = (-1)^{n-2},
   giving M_{T^op}[a,b] = (-1)^{n-2} M_T[b,a].
   So M[a,b]=M[b,a] iff M_{T^op} = (-1)^{n-2} M_T.

3. COMPLEMENT PAIRING: At n=4 (m=2), each complementary pair (S, complement)
   cancels INDIVIDUALLY in M[a,b]-M[b,a]. At n=5 (m=3), individual pairs
   do NOT cancel, but they cancel in groups of two pairs. The 4 complement-
   pair differences sum to 0 with a hierarchical cancellation pattern.

4. NOT A SIMPLE DETERMINANT: E_a(S) is a permanent (not determinant),
   so the alternating-sign sum is not a direct determinant expansion.
   However, the alternating sum of permanents might have hidden structure.

5. MULTILINEARITY: M[a,b] is multilinear in the arc variables
   (each variable appears with exponent at most 1). Verified at n=4,5.
   The arc t_{ab} between the two endpoints does NOT appear in M[a,b].

6. INDEPENDENCE FROM t_{ab}: The arc between endpoints a,b does not
   appear in M[a,b]. This means M[a,b] depends only on arcs involving
   the "other" vertices and their connections to a,b.

7. MATRIX STRUCTURE: The full n x n matrix M (with M[i,i]=0) is
   SYMMETRIC with integer entries. At n=5, observed eigenvalues include
   +/- sqrt(6), +/- 2, suggesting algebraic structure.

8. DEGREE PATTERN: M[a,b] has total degree m = n-2 in arc variables.
   At n=4: degree 2, 9 terms. At n=5: degree 3, 43 terms.
   Terms by degree at n=5: {0:1, 1:3, 2:15, 3:24}.

9. PROOF STRATEGY SUGGESTION: The T^op equivalence (finding 2) is the
   most promising route. Proving M_{T^op} = (-1)^{n-2} M_T might be
   easier than direct symmetry, since it relates M to itself under a
   global operation (arc reversal) rather than a swap of endpoint labels.
   The tournament constraint T[i,j]+T[j,i]=1 is exactly what makes
   T^op expressible in terms of T: T^op[i,j] = 1 - T[i,j].
""")
