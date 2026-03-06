#!/usr/bin/env python3
"""
Even r-powers via involution — attempt to prove M[a,b] is even in r.

Instance: opus-2026-03-06-S22

STRATEGY: The even cycle vanishing theorem (T148) uses a cycle-reversal
involution to show p_mu = 0 for mu with even parts. Can we adapt this
to show odd r-powers vanish in M[a,b]?

KEY INSIGHT from S10: Both properties stem from the T <-> T^op involution.
- Even cycle vanishing: reversing an even k-cycle changes sign by (-1)^{k-1} = -1
- Even r-powers: under r -> -r, path weights pick up (-1)^length

M[a,b] = sum_{S subset U} (-1)^|S| E_a(S+a) B_b(R+b)

where each E_a, B_b is a sum over Hamiltonian paths with weights
product of t_ij = r + s_ij.

The r^k coefficient of M[a,b] selects terms where exactly k edges
contribute their r part and the rest contribute their s part.

ATTEMPT 1: Direct involution on the (S, path_a, path_b) triples.
ATTEMPT 2: Generating function approach — M as trace/permanent.
ATTEMPT 3: Induction on n using vertex deletion.
"""

from itertools import permutations, combinations
from sympy import symbols, expand, Poly, collect, Rational, sympify, Symbol
from sympy import zeros as sympy_zeros
from collections import defaultdict
import random

def make_symbols(n):
    """Create symbolic edge weights for n vertices."""
    r = Symbol('r')
    s = {}
    for i in range(n):
        for j in range(n):
            if i < j:
                s[(i,j)] = Symbol(f's{i}{j}')
                s[(j,i)] = -s[(i,j)]
    return r, s

def edge_weight(i, j, r, s):
    """t(i,j) = r + s_ij"""
    return r + s[(i,j)]

def ham_paths_ending_at(vertex_set, target, r, s):
    """Sum of weights of Hamiltonian paths in T[vertex_set] ending at target."""
    vs = list(vertex_set)
    if len(vs) == 1:
        return 1  # trivial path
    total = 0
    others = [v for v in vs if v != target]
    for perm in permutations(others):
        path = list(perm) + [target]
        w = 1
        for k in range(len(path)-1):
            w *= edge_weight(path[k], path[k+1], r, s)
        total += w
    return expand(total)

def ham_paths_beginning_at(vertex_set, source, r, s):
    """Sum of weights of Hamiltonian paths in T[vertex_set] beginning at source."""
    vs = list(vertex_set)
    if len(vs) == 1:
        return 1
    total = 0
    others = [v for v in vs if v != source]
    for perm in permutations(others):
        path = [source] + list(perm)
        w = 1
        for k in range(len(path)-1):
            w *= edge_weight(path[k], path[k+1], r, s)
        total += w
    return expand(total)

def compute_M(a, b, n, r, s):
    """Compute M[a,b] = sum_S (-1)^|S| E_a(S+a) B_b(R+b)."""
    U = [v for v in range(n) if v != a and v != b]
    total = 0
    for mask in range(2**len(U)):
        S = [U[i] for i in range(len(U)) if (mask >> i) & 1]
        R = [u for u in U if u not in S]
        sign = (-1)**len(S)
        Sa = set(S + [a])
        Rb = set(R + [b])
        Ea = ham_paths_ending_at(Sa, a, r, s)
        Bb = ham_paths_beginning_at(Rb, b, r, s)
        total += sign * Ea * Bb
    return expand(total)


def analyze_r_coefficients(M_expr, r, n):
    """Extract coefficients of each r-power."""
    p = Poly(M_expr, r)
    coeffs = {}
    for monom, coeff in p.as_dict().items():
        deg = monom[0]
        coeffs[deg] = coeff
    return coeffs


# =============================================================
# ATTEMPT 1: Enumerate (S, path_a, path_b) triples and find involution
# =============================================================

def enumerate_triples_n4():
    """
    At n=4 with a=0, b=1, U={2,3}:
    List all (S, path_a, path_b) triples with their r-contributions.
    """
    print("=== ATTEMPT 1: Enumerate triples at n=4 ===")
    print("a=0, b=1, U={2,3}\n")

    n = 4; a = 0; b = 1
    r_sym, s_sym = make_symbols(n)

    U = [2, 3]

    triples = []
    for mask in range(4):
        S = [U[i] for i in range(2) if (mask >> i) & 1]
        R = [u for u in U if u not in S]
        sign = (-1)**len(S)

        Sa = sorted(S + [a])
        Rb = sorted(R + [b])

        # Enumerate paths ending at a in Sa
        others_a = [v for v in Sa if v != a]
        paths_a = list(permutations(others_a))

        # Enumerate paths beginning at b in Rb
        others_b = [v for v in Rb if v != b]
        paths_b = list(permutations(others_b))

        for pa in paths_a:
            path_a = list(pa) + [a]
            for pb in paths_b:
                path_b = [b] + list(pb)

                # Compute weight
                w = sign
                edges_a = [(path_a[k], path_a[k+1]) for k in range(len(path_a)-1)]
                edges_b = [(path_b[k], path_b[k+1]) for k in range(len(path_b)-1)]

                all_edges = edges_a + edges_b
                # Each edge contributes (r + s_ij)
                # r^1 coefficient: choose 1 edge for r, rest for s

                r1_coeff = 0
                for idx in range(len(all_edges)):
                    prod = sign
                    for j, (u,v) in enumerate(all_edges):
                        if j == idx:
                            prod *= 1  # r contribution
                        else:
                            prod *= s_sym[(u,v)]
                    r1_coeff += prod

                total_weight = sign
                for u,v in all_edges:
                    total_weight *= edge_weight(u, v, r_sym, s_sym)
                total_weight = expand(total_weight)

                triples.append({
                    'S': S, 'R': R, 'sign': sign,
                    'path_a': path_a, 'path_b': path_b,
                    'edges_a': edges_a, 'edges_b': edges_b,
                    'weight': total_weight,
                })

    print(f"Total triples: {len(triples)}")

    # Group by r-degree contribution
    for t in triples:
        p = Poly(t['weight'], r_sym)
        t['r_coeffs'] = {m[0]: c for m, c in p.as_dict().items()}

    # Sum all r^1 contributions
    r1_total = sum(t['r_coeffs'].get(1, 0) for t in triples)
    print(f"Total r^1 coefficient: {expand(r1_total)}")

    # Print each triple's r^1 contribution
    print("\nTriple-by-triple r^1 contributions:")
    for i, t in enumerate(triples):
        r1 = t['r_coeffs'].get(1, 0)
        if r1 != 0:
            print(f"  [{i}] S={t['S']}, path_a={t['path_a']}, path_b={t['path_b']}, "
                  f"sign={t['sign']}, r^1={r1}")

    # Now try to find a PAIRING involution
    print("\n--- Looking for involution on triples ---")
    # For each triple, try swapping an edge's direction (corresponding to
    # moving a vertex between S and R)

    # Key idea: for vertex u in U, consider triples where u ∈ S vs u ∈ R.
    # The "toggle u" operation moves u between S and R, changing sign.
    # But it also changes which paths are available.

    for u in U:
        print(f"\nToggle vertex {u}:")
        with_u = [t for t in triples if u in t['S']]
        without_u = [t for t in triples if u not in t['S']]

        r1_with = sum(t['r_coeffs'].get(1, 0) for t in with_u)
        r1_without = sum(t['r_coeffs'].get(1, 0) for t in without_u)

        print(f"  r^1 from S∋{u}: {expand(r1_with)}")
        print(f"  r^1 from S∌{u}: {expand(r1_without)}")
        print(f"  Sum: {expand(r1_with + r1_without)}")
        print(f"  Negation? {expand(r1_with + r1_without) == 0}")


# =============================================================
# ATTEMPT 2: Connection to permanent / determinant
# =============================================================

def attempt_permanent_connection():
    """
    The inclusion-exclusion formula for permanents:
    perm(A) = sum_{S⊆[n]} (-1)^|S| prod_i (sum_{j not in S} A[i,j])

    But M[a,b] has a different structure. Let's see if M can be written
    as a permanent or determinant of some matrix.
    """
    print("\n\n=== ATTEMPT 2: Permanent / determinant connection ===")

    n = 4; a = 0; b = 1
    r, s = make_symbols(n)

    M = compute_M(a, b, n, r, s)
    print(f"M[0,1] at n=4 = {M}")

    # Try: M[a,b] = permanent of some 2x2 matrix?
    # At n=4, the "path" from a to b through U has length 2.
    # There are 2 intermediate vertices (2 and 3).
    # Could M[a,b] = sum over orderings of U: product of appropriate weights?

    # Actually, at n=4:
    # M[0,1] = t02*t21 * t03*t31 ... no, it's more complex

    # Let me just check: is M the permanent of the 2x2 matrix
    # [[t(2,0)*t(0,3), t(2,1)*t(1,3)], [...]] ?

    # Actually, a cleaner formulation:
    # At n=4, M[0,1] should be a polynomial of degree 2 in r.
    # Let me factor it.

    coeffs = analyze_r_coefficients(M, r, n)
    print(f"\nr-coefficients:")
    for deg in sorted(coeffs.keys()):
        print(f"  r^{deg}: {coeffs[deg]}")

    # The r^0 term is the s-only part
    # The r^2 term should be a constant (number of path-pair covers)

    # Now try the "weighted Laplacian" approach
    # Kirchhoff's theorem: number of spanning trees = any cofactor of Laplacian
    # For directed graphs: arborescences = det of reduced Laplacian

    # Define the skew-adjacency matrix S_mat
    from sympy import Matrix
    S_mat = sympy_zeros(n, n)
    for i in range(n):
        for j in range(n):
            if i < j:
                S_mat[i,j] = s[(i,j)]
                S_mat[j,i] = -s[(i,j)]

    print(f"\nSkew matrix S = {S_mat}")

    # A = rJ' + S where J' = all-ones minus identity
    A = sympy_zeros(n, n)
    for i in range(n):
        for j in range(n):
            if i != j:
                A[i,j] = r + s[(i,j)]

    # Try: cofactor of (I - A)?  No, A has huge entries...
    # Try: adj(rI - S)?
    rI_S = r * Matrix.eye(n) - S_mat
    adj_rIS = rI_S.adjugate()

    print(f"\nadj(rI - S)[0,1] = {expand(adj_rIS[0,1])}")
    print(f"M[0,1]           = {M}")
    print(f"Match: {expand(adj_rIS[0,1] - M) == 0}")

    # What about cofactor of rI-S with sign adjustment?
    # cofactor(i,j) = (-1)^{i+j} det(minor)
    minor = rI_S.minor_submatrix(0, 1)
    cof01 = (-1)**(0+1) * minor.det()
    print(f"\ncofactor(rI-S, 0,1) = {expand(cof01)}")
    print(f"Match: {expand(cof01 - M) == 0}")

    # Let me try (rI + S)
    rI_pS = r * Matrix.eye(n) + S_mat
    adj_rIpS = rI_pS.adjugate()
    print(f"\nadj(rI + S)[0,1] = {expand(adj_rIpS[0,1])}")
    print(f"Match: {expand(adj_rIpS[0,1] - M) == 0}")


# =============================================================
# ATTEMPT 3: Induction via vertex deletion
# =============================================================

def attempt_induction():
    """
    Try: M_n[a,b] = sum_v (something involving M_{n-1}).

    If M[a,b] is even in r at n-1, can we show it's even at n?

    Approach: fix an intermediate vertex u ∈ U. Split the sum
    based on whether u is in S or R.
    """
    print("\n\n=== ATTEMPT 3: Induction via vertex deletion ===")

    # At n=5, a=0, b=1, U={2,3,4}
    # Pick u=2.
    # If u ∈ S: M contribution involves E_a(S+a) where 2∈S, B_b(R+b) where 2∉R
    # If u ∈ R: M contribution involves E_a(S+a) where 2∉S, B_b(R+b) where 2∈R

    n = 5; a = 0; b = 1; u = 2
    r_sym, s_sym = make_symbols(n)

    U = [v for v in range(n) if v != a and v != b]
    U_minus_u = [v for v in U if v != u]

    # When u ∈ S: S = S' ∪ {u} where S' ⊆ U\{u}
    # E_a(S'+{u}+a) = paths ending at a through S'∪{u}∪{a}
    # B_b(R'+b) where R' = (U\{u})\S' = U_minus_u \ S'
    # sign = (-1)^{|S'|+1}

    M_u_in_S = 0
    M_u_in_R = 0

    for mask in range(2**len(U_minus_u)):
        Sp = [U_minus_u[i] for i in range(len(U_minus_u)) if (mask >> i) & 1]
        Rp = [v for v in U_minus_u if v not in Sp]

        # u in S case
        S_full = Sp + [u]
        R_full = Rp
        sign_S = (-1)**len(S_full)
        Sa = set(S_full + [a])
        Rb = set(R_full + [b])
        Ea_S = ham_paths_ending_at(Sa, a, r_sym, s_sym)
        Bb_S = ham_paths_beginning_at(Rb, b, r_sym, s_sym)
        M_u_in_S += sign_S * Ea_S * Bb_S

        # u in R case
        S_full2 = Sp
        R_full2 = Rp + [u]
        sign_R = (-1)**len(S_full2)
        Sa2 = set(S_full2 + [a])
        Rb2 = set(R_full2 + [b])
        Ea_R = ham_paths_ending_at(Sa2, a, r_sym, s_sym)
        Bb_R = ham_paths_beginning_at(Rb2, b, r_sym, s_sym)
        M_u_in_R += sign_R * Ea_R * Bb_R

    M_u_in_S = expand(M_u_in_S)
    M_u_in_R = expand(M_u_in_R)
    M_total = expand(M_u_in_S + M_u_in_R)

    print(f"n=5, a=0, b=1, splitting on u={u}")

    # Check r-parity of each part
    p_S = Poly(M_u_in_S, r_sym)
    p_R = Poly(M_u_in_R, r_sym)

    print(f"\nM (u∈S part):")
    for deg in sorted(set(m[0] for m in p_S.as_dict())):
        coeff = sum(c for m, c in p_S.as_dict().items() if m[0] == deg)
        print(f"  r^{deg}: {coeff}")

    print(f"\nM (u∈R part):")
    for deg in sorted(set(m[0] for m in p_R.as_dict())):
        coeff = sum(c for m, c in p_R.as_dict().items() if m[0] == deg)
        print(f"  r^{deg}: {coeff}")

    # Do the parts individually have even r-powers?
    odd_S = any(m[0] % 2 == 1 for m in p_S.as_dict())
    odd_R = any(m[0] % 2 == 1 for m in p_R.as_dict())
    print(f"\nu∈S part has odd r-powers: {odd_S}")
    print(f"u∈R part has odd r-powers: {odd_R}")

    if odd_S and odd_R:
        print("=> Cancellation between u∈S and u∈R needed (not individual property)")

    # Now check: does M(u∈S) + M(u∈R) have only even r-powers?
    p_total = Poly(M_total, r_sym)
    odd_total = any(m[0] % 2 == 1 for m, c in p_total.as_dict().items() if c != 0)
    print(f"\nTotal has odd r-powers: {odd_total}")

    # Check the toggle pairing structure
    print(f"\n--- Toggle u={u} pairing ---")
    print(f"r^1 in u∈S: {sum(c for m,c in p_S.as_dict().items() if m[0]==1)}")
    print(f"r^1 in u∈R: {sum(c for m,c in p_R.as_dict().items() if m[0]==1)}")
    print(f"Sum: {expand(sum(c for m,c in p_S.as_dict().items() if m[0]==1) + sum(c for m,c in p_R.as_dict().items() if m[0]==1))}")


# =============================================================
# ATTEMPT 4: The "Euler characteristic" interpretation
# =============================================================

def attempt_euler_char():
    """
    M[a,b] = sum_S (-1)^|S| f(S) g(U\S)

    This is a Euler characteristic computation on the Boolean lattice.
    The alternating sum can be interpreted as an inclusion-exclusion
    or as a simplicial homology computation.

    If f and g have specific symmetry under complement (S <-> U\S),
    this constrains which "Fourier modes" survive.

    On the Boolean lattice 2^U, functions have a Fourier expansion:
    f(S) = sum_{T⊆U} f_hat(T) * chi_T(S) where chi_T(S) = (-1)^|S∩T|

    The (-1)^|S| factor is chi_U(S). So:
    M[a,b] = sum_S chi_U(S) f(S) g(U\S)

    This is a convolution in the Fourier basis. Let's compute what
    constraints "even in r" imposes on the Fourier coefficients.
    """
    print("\n\n=== ATTEMPT 4: Boolean Fourier analysis ===")

    n = 4; a = 0; b = 1
    r_sym, s_sym = make_symbols(n)
    U = [2, 3]

    # Compute f(S) = E_a(S+a) and g(S) = B_b(S+b) for all S ⊆ U
    f_vals = {}
    g_vals = {}

    for mask in range(2**len(U)):
        S = frozenset(U[i] for i in range(len(U)) if (mask >> i) & 1)
        f_vals[S] = ham_paths_ending_at(S | {a}, a, r_sym, s_sym)
        g_vals[S] = ham_paths_beginning_at(S | {b}, b, r_sym, s_sym)

    print("f(S) = E_a(S+a):")
    for S in sorted(f_vals.keys(), key=len):
        print(f"  f({set(S)}) = {f_vals[S]}")

    print("\ng(S) = B_b(S+b):")
    for S in sorted(g_vals.keys(), key=len):
        print(f"  g({set(S)}) = {g_vals[S]}")

    # M[a,b] = sum_S (-1)^|S| f(S) g(U\S)
    # where U\S = complement of S in U
    U_set = frozenset(U)
    M_check = 0
    for S in f_vals:
        R = U_set - S
        M_check += (-1)**len(S) * f_vals[S] * g_vals[R]
    M_check = expand(M_check)

    M_direct = compute_M(a, b, n, r_sym, s_sym)
    print(f"\nM[0,1] (convolution) = {M_check}")
    print(f"M[0,1] (direct)      = {M_direct}")
    print(f"Match: {expand(M_check - M_direct) == 0}")

    # Now compute Fourier transforms
    # f_hat(T) = 2^{-|U|} * sum_S (-1)^{|S∩T|} f(S)
    print("\n--- Fourier transforms ---")

    def fourier_transform(vals, U_list):
        U_set = frozenset(U_list)
        n_U = len(U_list)
        hat = {}
        for mask_T in range(2**n_U):
            T = frozenset(U_list[i] for i in range(n_U) if (mask_T >> i) & 1)
            total = 0
            for mask_S in range(2**n_U):
                S = frozenset(U_list[i] for i in range(n_U) if (mask_S >> i) & 1)
                total += (-1)**len(S & T) * vals[S]
            hat[T] = expand(total * Rational(1, 2**n_U))
        return hat

    f_hat = fourier_transform(f_vals, U)
    g_hat = fourier_transform(g_vals, U)

    print("f_hat(T):")
    for T in sorted(f_hat.keys(), key=len):
        print(f"  f_hat({set(T)}) = {f_hat[T]}")

    print("\ng_hat(T):")
    for T in sorted(g_hat.keys(), key=len):
        print(f"  g_hat({set(T)}) = {g_hat[T]}")

    # In the Fourier basis, the convolution with (-1)^|S| becomes:
    # M[a,b] = 2^|U| * sum_T f_hat(T) * g_hat(U\T) * (-1)^{|U\T|}
    # Wait, need to be more careful.

    # M = sum_S (-1)^|S| f(S) g(U\S) = sum_S chi_U(S) f(S) g(U\S)
    # Using Fourier: f(S) = sum_T f_hat(T) chi_T(S)
    #                g(S) = sum_V g_hat(V) chi_V(S)
    # g(U\S) = sum_V g_hat(V) chi_V(U\S) = sum_V g_hat(V) (-1)^{|V\S|}
    #        = sum_V g_hat(V) (-1)^{|V|} (-1)^{|V∩S|} ... hmm

    # Actually chi_V(U\S) = (-1)^{|(U\S)∩V|} = (-1)^{|V| - |S∩V|} = (-1)^|V| * (-1)^{|S∩V|}
    # So chi_V(U\S) = (-1)^|V| chi_V(S)

    # Therefore:
    # M = sum_S chi_U(S) * [sum_T f_hat(T) chi_T(S)] * [sum_V g_hat(V) (-1)^|V| chi_V(S)]
    # = sum_{T,V} f_hat(T) g_hat(V) (-1)^|V| * sum_S chi_U(S) chi_T(S) chi_V(S)
    # = sum_{T,V} f_hat(T) g_hat(V) (-1)^|V| * sum_S (-1)^{|S∩U| + |S∩T| + |S∩V|}
    # = sum_{T,V} f_hat(T) g_hat(V) (-1)^|V| * sum_S (-1)^{|S∩(U Δ T Δ V)|}

    # The inner sum = 2^|U| if U Δ T Δ V = ∅, else 0.
    # U Δ T Δ V = ∅ iff T Δ V = U (since all are subsets of U)
    # iff V = U Δ T = U\T ∪ T\U... since T,V ⊆ U, V = U\T.

    # So M = 2^|U| sum_T f_hat(T) g_hat(U\T) (-1)^|U\T|

    M_fourier = 0
    for T in f_hat:
        UT = U_set - T
        M_fourier += f_hat[T] * g_hat[UT] * (-1)**len(UT)
    M_fourier = expand(4 * M_fourier)  # 2^|U| = 4

    print(f"\nM from Fourier: {M_fourier}")
    print(f"Match: {expand(M_fourier - M_direct) == 0}")

    # So M[a,b] = 2^|U| * sum_T f_hat(T) g_hat(U\T) (-1)^{|U\T|}

    # For M to be even in r, we need each term in this Fourier expression
    # to be even in r. Let's check:
    print("\n--- Fourier term r-parity ---")
    for T in sorted(f_hat.keys(), key=len):
        UT = U_set - T
        term = f_hat[T] * g_hat[UT] * (-1)**len(UT)
        term = expand(term)
        if term == 0:
            print(f"  T={set(T)}: 0")
            continue
        p = Poly(term, r_sym)
        degs = sorted(set(m[0] for m in p.as_dict()))
        has_odd = any(d % 2 == 1 for d in degs)
        print(f"  T={set(T)}: degrees={degs}, has_odd={has_odd}, value={term}")


# =============================================================
# ATTEMPT 5: Connection to the ADJUGATE of (I - tA) where A is adjacency
# =============================================================

def attempt_transfer_matrix_formulation():
    """
    In many graph theory contexts, the transfer matrix arises as
    adj(I - zA) or (I - zA)^{-1} where A is the adjacency matrix.

    For a tournament, A[i,j] = t(i,j) = r + s_ij.

    (I - zA)^{-1} = sum_k z^k A^k counts weighted walks.

    The (a,b) entry of A^{n-1} counts weighted walks of length n-1 from a to b.

    The Hamiltonian path count is related to this by inclusion-exclusion:
    M[a,b] = sum_S (-1)^|S| [A_{T[V\S]}^{|V\S|-1}]_{a,b}

    Wait, that's not quite right either. Let me think...

    Actually, the MFMC (matrix-forest theorem) gives:
    adj(I - A)[i,j] = sum over spanning forests with root structure...

    For directed graphs, (I-A)^{-1} = sum A^k, and
    adj(I-A) = det(I-A) * (I-A)^{-1}

    The key: det(I-A) for our weighted tournament is a polynomial in r.
    If it's even in r, and adj(I-A)[a,b] is the product of det and (I-A)^{-1}[a,b],
    then we'd need (I-A)^{-1}[a,b] to also be even...

    This probably doesn't lead anywhere directly. Let me just check.
    """
    print("\n\n=== ATTEMPT 5: Transfer matrix via (I - zA) ===")

    n = 4; a = 0; b = 1
    r_sym, s_sym = make_symbols(n)

    from sympy import Matrix, eye
    z = Symbol('z')

    A = sympy_zeros(n, n)
    for i in range(n):
        for j in range(n):
            if i != j:
                A[i,j] = r_sym + s_sym[(i,j)]

    # A^{n-1} at n=4 is A^3
    A3 = A**3
    walk_ab = expand(A3[a, b])
    print(f"A^3[0,1] (walks of length 3 from 0 to 1) = {walk_ab}")

    M_direct = compute_M(a, b, n, r_sym, s_sym)
    print(f"M[0,1] (Hamiltonian paths) = {M_direct}")
    print(f"Match: {expand(walk_ab - M_direct) == 0}")

    # Of course they don't match — A^k counts ALL walks, not just Hamiltonian
    # The difference is walks that revisit vertices.
    diff = expand(walk_ab - M_direct)
    print(f"Difference (non-Hamiltonian walks): {diff}")

    # Is the difference also even in r?
    p_diff = Poly(diff, r_sym)
    degs_diff = sorted(set(m[0] for m in p_diff.as_dict()))
    print(f"Difference r-degrees: {degs_diff}")

    # A^k[a,b] is a product of k edge weights, each = r + s_ij.
    # Total degree in r is k. For k=n-1, max r degree is n-1.

    # A is a matrix of degree 1 in r. A(-r,s) = -A^T(r,s) because:
    # A[i,j](-r) = -r + s_ij = -(r - s_ij) = -(r + s_ji) = -A[j,i](r)
    # So A(-r) = -A^T(r).

    # Therefore A^k(-r) = (-1)^k (A^T)^k(r) = (-1)^k (A^k)^T(r)
    # So A^k(-r)[a,b] = (-1)^k A^k(r)[b,a]

    # For k=n-1 (odd for even n, even for odd n):
    # A^{n-1}(-r)[a,b] = (-1)^{n-1} A^{n-1}(r)[b,a]

    # If A^{n-1}[a,b] is even in r, then A^{n-1}(-r)[a,b] = A^{n-1}(r)[a,b]
    # So A^{n-1}(r)[a,b] = (-1)^{n-1} A^{n-1}(r)[b,a]
    # For n even: A^{n-1}[a,b] = -A^{n-1}[b,a] (antisymmetric!)
    # For n odd: A^{n-1}[a,b] = A^{n-1}[b,a] (symmetric!)

    # But A^{n-1}[a,b] counts ALL walks, including non-Hamiltonian ones.
    # Is A^{n-1}[a,b] even in r?
    p_walk = Poly(walk_ab, r_sym)
    degs_walk = sorted(set(m[0] for m in p_walk.as_dict()))
    print(f"\nA^3[0,1] r-degrees: {degs_walk}")

    # Check A^3[0,1](-r) vs A^3[1,0](r)
    A3_ba = expand(A3[b, a])
    walk_neg = expand(walk_ab.subs(r_sym, -r_sym))
    print(f"A^3[0,1](-r) = {walk_neg}")
    print(f"(-1)^3 * A^3[1,0](r) = {expand(-A3_ba)}")
    print(f"Match: {expand(walk_neg - (-A3_ba)) == 0}")

    # Now the KEY question: A^k[a,b] counts walks with repeats.
    # M[a,b] = sum_S (-1)^|S| E_a B_b is the Hamiltonian part.
    # The non-Hamiltonian part is the difference.
    # If BOTH the total walks AND the non-Hamiltonian walks are even in r,
    # then M is even in r.

    # But A^{n-1}[a,b] is NOT even in r (it has all degrees 0..n-1).
    # So evenness of M comes from specific cancellation in inclusion-exclusion.

    print(f"\nA^3 is NOT even in r. Evenness of M is a property of the I-E formula.")


def main():
    print("=" * 70)
    print("EVEN R-POWERS VIA INVOLUTION — PROOF ATTEMPTS")
    print("=" * 70)

    enumerate_triples_n4()
    attempt_permanent_connection()
    attempt_induction()
    attempt_euler_char()
    attempt_transfer_matrix_formulation()

    print("\n\n" + "=" * 70)
    print("SUMMARY OF FINDINGS")
    print("=" * 70)
    print("""
The key structural facts:
1. A(-r,s) = -A^T(r,s) where A is the weighted adjacency matrix.
2. This gives A^k(-r)[a,b] = (-1)^k A^k[b,a](r) for walk counts.
3. For Hamiltonian paths specifically: M[a,b](-r) = M[b,a](r).
4. Symmetry M[a,b] = M[b,a] <=> M even in r <=> only even r-powers.

The involution approach from even cycle vanishing doesn't directly
adapt because:
- ECF works on permutations (cycle structure)
- Transfer matrix works on path-pair covers (subset structure)
- The natural involution (path reversal) gives M(-r)=M^T, not M=M^T

REMAINING PROMISING DIRECTIONS:
- Boolean Fourier analysis: M expressed as Fourier convolution on 2^U
- The Fourier terms individually may NOT be even, but their sum is
- Need a "parity selection" argument for why odd-frequency Fourier
  modes cancel in the convolution
""")


if __name__ == '__main__':
    main()
