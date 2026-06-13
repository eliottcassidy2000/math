#!/usr/bin/env python3
"""
Deep analysis of transfer matrix symmetry: F(a,b) - F(b,a) = 0.

Goal: Understand the cancellation mechanism well enough to prove it for all n.

Approach:
1. For each subset S, compute D(S) = E_a(S)*B_b(R) - E_b(S)*B_a(R)
2. Track individual monomial contributions to see HOW they cancel
3. Look for an involution on (S, path-pair) terms
4. Test if the cancellation has a "local" structure (pairs of subsets that cancel)

Instance: opus-2026-03-06-S6
"""
from sympy import Symbol, expand, Poly, degree
from itertools import permutations, combinations
from collections import defaultdict

def make_arc_var(i, j, arc_vars):
    if i < j:
        return arc_vars[(i, j)]
    elif i > j:
        return 1 - arc_vars[(j, i)]
    else:
        return 0

def sym_h_end(arc_vars, S_list, v):
    """Symbolic count of Ham paths through S ending at v."""
    if not S_list:
        return 1
    S = list(S_list)
    m = len(S)
    total = 0
    for perm in permutations(S):
        prod = 1
        for k in range(m - 1):
            prod = prod * make_arc_var(perm[k], perm[k + 1], arc_vars)
        prod = prod * make_arc_var(perm[m - 1], v, arc_vars)
        total = total + prod
    return expand(total)

def sym_h_start(arc_vars, R_list, v):
    """Symbolic count of Ham paths through R starting from v."""
    if not R_list:
        return 1
    R = list(R_list)
    m = len(R)
    total = 0
    for perm in permutations(R):
        prod = make_arc_var(v, perm[0], arc_vars)
        for k in range(m - 1):
            prod = prod * make_arc_var(perm[k], perm[k + 1], arc_vars)
        total = total + prod
    return expand(total)

def analyze_per_subset(n):
    """Analyze D(S) = E_a(S)*B_b(R) - E_b(S)*B_a(R) for each subset S."""
    arc_vars = {}
    for i in range(n):
        for j in range(i + 1, n):
            arc_vars[(i, j)] = Symbol(f't{i}{j}')

    a, b = 0, 1
    U = [v for v in range(n) if v != a and v != b]
    m = len(U)

    total_diff = 0
    subset_diffs = {}

    for smask in range(1 << m):
        S = [U[k] for k in range(m) if smask & (1 << k)]
        R = [U[k] for k in range(m) if not (smask & (1 << k))]
        sign = (-1) ** len(S)

        Ea = sym_h_end(arc_vars, S, a)
        Eb = sym_h_end(arc_vars, S, b)
        Ba = sym_h_start(arc_vars, R, a)
        Bb = sym_h_start(arc_vars, R, b)

        D = expand(Ea * Bb - Eb * Ba)
        signed_D = expand(sign * D)

        subset_diffs[tuple(S)] = {
            'S': S, 'R': R, '|S|': len(S), 'sign': sign,
            'D(S)': D, 'signed_D': signed_D,
            'Ea': Ea, 'Eb': Eb, 'Ba': Ba, 'Bb': Bb
        }
        total_diff += signed_D

    total_diff = expand(total_diff)
    return subset_diffs, total_diff, arc_vars

def analyze_complement_pairing(n):
    """Check if D(S) and D(U\S) have a simple relationship."""
    diffs, total, arc_vars = analyze_per_subset(n)

    print(f"\n{'='*60}")
    print(f"COMPLEMENT PAIRING ANALYSIS at n={n}")
    print(f"{'='*60}")

    a, b = 0, 1
    U = [v for v in range(n) if v != a and v != b]
    m = len(U)

    for smask in range(1 << (m - 1)):  # Only half, pair with complement
        S = tuple(U[k] for k in range(m) if smask & (1 << k))
        R = tuple(U[k] for k in range(m) if not (smask & (1 << k)))

        D_S = diffs[S]['D(S)']
        D_R = diffs[R]['D(S)']
        sign_S = diffs[S]['sign']
        sign_R = diffs[R]['sign']

        print(f"\nS={S}, R={R}")
        print(f"  sign(S)={sign_S:+d}, sign(R)={sign_R:+d}")
        print(f"  D(S) = {D_S}")
        print(f"  D(R) = {D_R}")

        # Check relationship
        summ = expand(D_S + D_R)
        diff = expand(D_S - D_R)
        print(f"  D(S)+D(R) = {summ}")
        print(f"  D(S)-D(R) = {diff}")

        # Check signed sum
        signed_sum = expand(sign_S * D_S + sign_R * D_R)
        print(f"  sign(S)*D(S) + sign(R)*D(R) = {signed_sum}")

def analyze_swapping_structure(n):
    """Analyze what happens when we swap vertices between S and R."""
    diffs, total, arc_vars = analyze_per_subset(n)

    print(f"\n{'='*60}")
    print(f"SWAPPING STRUCTURE at n={n}")
    print(f"{'='*60}")

    a, b = 0, 1
    U = [v for v in range(n) if v != a and v != b]
    m = len(U)

    # For each pair of subsets differing by one element
    for smask in range(1 << m):
        S = tuple(U[k] for k in range(m) if smask & (1 << k))
        for k in range(m):
            if smask & (1 << k):  # vertex U[k] is in S
                # Move U[k] from S to R
                smask2 = smask ^ (1 << k)
                S2 = tuple(U[j] for j in range(m) if smask2 & (1 << j))
                if S < S2:  # avoid double counting
                    v_moved = U[k]
                    D_S = diffs[S]['signed_D']
                    D_S2 = diffs[S2]['signed_D']
                    pair_sum = expand(D_S + D_S2)
                    if pair_sum != 0:
                        print(f"  S={S} <-> S2={S2} (moving {v_moved}): sum = {pair_sum}")

def analyze_determinantal_structure(n):
    """
    Key idea: D(S) = det([[E_a(S), B_a(R)], [E_b(S), B_b(R)]]).
    Can we express the sum as a single determinant?
    """
    diffs, total, arc_vars = analyze_per_subset(n)

    print(f"\n{'='*60}")
    print(f"DETERMINANTAL STRUCTURE at n={n}")
    print(f"{'='*60}")

    a, b = 0, 1
    U = [v for v in range(n) if v != a and v != b]
    m = len(U)

    # Define the "inclusion-exclusion path counts"
    # F_v = sum_S (-1)^|S| E_v(S) x^S  (generating function over subsets)
    # G_v = sum_S (-1)^|S| B_v(R) x^R

    # The off-diagonal of M is sum_S (-1)^|S| E_a(S)*B_b(R)
    # = sum_S (-1)^|S| E_a(S)*B_b(U\S)

    # What if we write E_a(S) as a "signed permanent" or something?

    # Let's look at whether each D(S) factors
    for S_key in sorted(diffs.keys(), key=len):
        info = diffs[S_key]
        D = info['D(S)']
        if D != 0:
            print(f"  S={S_key}: D(S) = {D}")

    print(f"\n  Total: {total}")

def test_row_linearity(n):
    """
    Test if E_a(S) and B_a(R) satisfy any linear relationship
    that could explain the determinantal cancellation.
    """
    arc_vars = {}
    for i in range(n):
        for j in range(i + 1, n):
            arc_vars[(i, j)] = Symbol(f't{i}{j}')

    a, b = 0, 1
    U = [v for v in range(n) if v != a and v != b]
    m = len(U)

    print(f"\n{'='*60}")
    print(f"ROW LINEARITY TEST at n={n}")
    print(f"{'='*60}")

    # For each vertex v in U, check if E_a({v}) + E_b({v}) has a pattern
    for v in U:
        Ea_v = sym_h_end(arc_vars, [v], a)
        Eb_v = sym_h_end(arc_vars, [v], b)
        Ba_v = sym_h_start(arc_vars, [v], a)
        Bb_v = sym_h_start(arc_vars, [v], b)

        print(f"\n  v={v}:")
        print(f"    E_a({{v}}) = {Ea_v} = T[{v},{a}]")
        print(f"    E_b({{v}}) = {Eb_v} = T[{v},{b}]")
        print(f"    B_a({{v}}) = {Ba_v} = T[{a},{v}]")
        print(f"    B_b({{v}}) = {Bb_v} = T[{b},{v}]")
        print(f"    E_a + E_b = {expand(Ea_v + Eb_v)}")
        print(f"    B_a + B_b = {expand(Ba_v + Bb_v)}")

        # Key: T[v,a] + T[v,b] and T[a,v] + T[b,v]
        # Using tournament: T[v,a] = 1-T[a,v], T[v,b] = 1-T[b,v]
        # So E_a({v}) + E_b({v}) = (1-T[a,v]) + (1-T[b,v]) = 2 - T[a,v] - T[b,v]
        # And B_a({v}) + B_b({v}) = T[a,v] + T[b,v]
        # Therefore E_sum + B_sum = 2!

def test_cauchy_binet(n):
    """
    Test if the transfer matrix can be written as C^T * D for some matrices C, D,
    which would make symmetry follow from C = D (or similar).

    M[a,b] = sum_S (-1)^|S| E_a(S) * B_b(R)

    This looks like a matrix product if we index by subsets S:
    M = E^T * Diag((-1)^|S|) * B

    where E[S,v] = E_v(S) and B[S,v] = B_v(U\S).
    """
    arc_vars = {}
    for i in range(n):
        for j in range(i + 1, n):
            arc_vars[(i, j)] = Symbol(f't{i}{j}')

    a, b = 0, 1
    U = [v for v in range(n) if v != a and v != b]
    m = len(U)

    print(f"\n{'='*60}")
    print(f"CAUCHY-BINET DECOMPOSITION at n={n}")
    print(f"{'='*60}")

    # Build the E and B matrices (indexed by subset masks)
    E_mat = {}  # E_mat[smask, v] = E_v(S(smask))
    B_mat = {}  # B_mat[smask, v] = B_v(R(smask))

    for smask in range(1 << m):
        S = [U[k] for k in range(m) if smask & (1 << k)]
        R = [U[k] for k in range(m) if not (smask & (1 << k))]

        E_mat[smask, 0] = sym_h_end(arc_vars, S, a)
        E_mat[smask, 1] = sym_h_end(arc_vars, S, b)
        B_mat[smask, 0] = sym_h_start(arc_vars, R, a)
        B_mat[smask, 1] = sym_h_start(arc_vars, R, b)

    # M = sum_S (-1)^|S| * [E_a(S), E_b(S)]^T * [B_a(R), B_b(R)]
    # This is a rank-1 update per subset S.
    # Equivalently: M = E^T * Lambda * B where Lambda[S,S'] = delta_{S,S'} * (-1)^|S|

    # For symmetry M = M^T, need E^T * Lambda * B = B^T * Lambda * E
    # i.e., Lambda * B * E^{-T} is symmetric, or E * Lambda * B = (E * Lambda * B)^T

    # Alternative: define C[S] = (-1)^|S| * E[S], then M = C^T * B.
    # Symmetry: C^T * B = B^T * C, i.e., sum_S C[S,a]*B[S,b] = sum_S B[S,a]*C[S,b]

    # Let's check if E[S,a]*B[S,b] = E[U\S,b]*B[U\S,a] (with appropriate signs)
    # This would give a PAIRING between S and U\S that proves symmetry

    print("Testing E(S,a)*B(S,b) vs E(R,b)*B(R,a) (complement pairing):")
    for smask in range(1 << m):
        S = tuple(U[k] for k in range(m) if smask & (1 << k))
        R = tuple(U[k] for k in range(m) if not (smask & (1 << k)))
        rmask = ((1 << m) - 1) ^ smask

        prod_ab = expand(E_mat[smask, 0] * B_mat[smask, 1])
        prod_ba_comp = expand(E_mat[rmask, 1] * B_mat[rmask, 0])

        diff = expand(prod_ab - prod_ba_comp)
        sign_s = (-1)**len(S)
        sign_r = (-1)**len(R)

        if diff != 0 or smask <= rmask:
            print(f"  S={S}: E_a(S)*B_b(R)={prod_ab}")
            print(f"  R={R}: E_b(R)*B_a(S)={prod_ba_comp}")
            print(f"  signs: (-1)^|S|={sign_s}, (-1)^|R|={sign_r}")
            print(f"  diff={diff}")
            print()

if __name__ == '__main__':
    n = 4
    print(f"{'#'*60}")
    print(f"# TRANSFER MATRIX SYMMETRY ANALYSIS at n={n}")
    print(f"{'#'*60}")

    # 1. Per-subset analysis
    diffs, total, arc_vars = analyze_per_subset(n)
    print(f"\nTotal difference F(a,b) - F(b,a) = {total}")

    # 2. Complement pairing
    analyze_complement_pairing(n)

    # 3. Row linearity
    test_row_linearity(n)

    # 4. Cauchy-Binet
    test_cauchy_binet(n)

    # Do the same at n=5
    if True:
        n = 5
        print(f"\n\n{'#'*60}")
        print(f"# n={n}")
        print(f"{'#'*60}")
        diffs5, total5, _ = analyze_per_subset(n)
        print(f"\nTotal diff = {total5}")
        analyze_complement_pairing(n)
        test_cauchy_binet(n)
