#!/usr/bin/env python3
"""
Compare M[a,b] with cofactors explicitly to find the relationship.

Instance: opus-2026-03-06-S22
"""

from itertools import permutations
from sympy import symbols, expand, Symbol, Matrix, det, zeros, ones, Poly, factor

def make_symbols(n):
    r = Symbol('r')
    s = {}
    for i in range(n):
        for j in range(n):
            if i < j:
                s[(i,j)] = Symbol(f's{i}{j}')
                s[(j,i)] = -s[(i,j)]
    return r, s

def edge_weight(i, j, r, s):
    return r + s[(i,j)]

def ham_paths_ending_at(vertex_set, target, r, s):
    vs = list(vertex_set)
    if len(vs) == 1:
        return 1
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

def compute_M(a, b, vertices, r, s):
    V = list(vertices)
    U = [v for v in V if v != a and v != b]
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


def explicit_comparison(n):
    print(f"\n{'='*60}")
    print(f"EXPLICIT COMPARISON: n={n}")
    print(f"{'='*60}")

    r, s = make_symbols(n)
    V = list(range(n))

    A = zeros(n)
    for i in range(n):
        for j in range(n):
            if i != j:
                A[i,j] = r + s[(i,j)]

    a, b = 0, 1
    M_ab = compute_M(a, b, V, r, s)

    print(f"\n  M[{a},{b}] = {M_ab}")

    # Print some cofactors
    for i in range(n):
        for j in range(n):
            cof = expand((-1)**(i+j) * det(A.minor_submatrix(i, j)))
            if n <= 4:
                print(f"  cofactor({i},{j}) = {cof}")

    # Also try: permanent-like objects, or traces of adj(A)
    # Or: what if we use a DIFFERENT matrix?
    # The "Hamiltonian path matrix" might involve (I - A)^{-1} type expressions.

    # Try: the matrix B where B[i,j] = permanent of A with row i, col j deleted
    # (This is the "permanent minor" -- too expensive for large n but try n=3)

    if n == 3:
        # Compute all permanents of 2x2 submatrices
        from sympy import Rational
        print(f"\n  --- Permanent minors ---")
        for i in range(n):
            for j in range(n):
                sub = A.minor_submatrix(i, j)
                perm_ij = expand(sub[0,0]*sub[1,1] + sub[0,1]*sub[1,0])
                print(f"  perm_minor({i},{j}) = {perm_ij}")
                if expand(M_ab - perm_ij) == 0:
                    print(f"  *** MATCH: M[{a},{b}] = perm_minor({i},{j})")
                if expand(M_ab + perm_ij) == 0:
                    print(f"  *** MATCH: M[{a},{b}] = -perm_minor({i},{j})")

    # Key idea: M[a,b] involves HAMILTONIAN paths, which are permanents, not determinants.
    # The inclusion-exclusion in M converts a permanent into... what exactly?
    #
    # For a 2-vertex partition {S+a, R+b}:
    #   E_a(S+a) = permanent contribution from paths ending at a through S+a
    #   B_b(R+b) = permanent contribution from paths starting at b through R+b
    # And we sum over all partitions with alternating signs.
    #
    # This is the inclusion-exclusion formula for the permanent of a certain matrix!
    # Specifically, M[a,b] should relate to the permanent of a matrix derived from A.

    # Actually, let me think about what the "transfer matrix" really is.
    # In the Grinberg-Stanley / Redei-Berge framework:
    # H(T) = number of Hamiltonian paths = permanent-like sum
    # But with the inclusion-exclusion, M is more like the "immanant" with the sign character.

    # Let me try: does M[a,b] = sum over permutations sigma of V\{a,b}
    # of some signed product?

    # Actually, M[a,b] = sum_S (-1)^|S| E_a(S+a) B_b(R+b)
    # This is a signed sum over 2-path-covers of V.
    # A 2-path-cover is: a path p1 ending at a, a path p2 starting at b,
    # with V(p1) ∪ V(p2) = V, V(p1) ∩ V(p2) = ∅.

    # The sign is (-1)^{|V(p1)|-1} = (-1)^{|S|}.

    # This is reminiscent of the formula for det in terms of cycle covers,
    # but for paths instead of cycles.

    # NEW IDEA: Think of M[a,b] as a HAFNIAN or PFAFFIAN of a modified matrix.
    # Or: relate it to the permanent of the "reduced" matrix.

    print(f"\n  --- Difference M[{a},{b}] - cofactor(1,0) ---")
    cof_10 = expand((-1)**(1+0) * det(A.minor_submatrix(1, 0)))
    diff = expand(M_ab - cof_10)
    print(f"  M[{a},{b}] - cof(1,0) = {diff}")
    print(f"  M[{a},{b}] = {M_ab}")
    print(f"  cof(1,0) = {cof_10}")

    # What about the PERMANENT of A minor?
    # perm(A with row b, col a deleted)?
    if n <= 4:
        sub = A.minor_submatrix(b, a)  # remove row b=1, col a=0
        # compute permanent of sub
        m = sub.shape[0]
        from itertools import permutations as perms
        perm_val = 0
        for p in perms(range(m)):
            prod = 1
            for i in range(m):
                prod *= sub[i, p[i]]
            perm_val += prod
        perm_val = expand(perm_val)
        print(f"\n  perm(A_{{ba}}) = {perm_val}")
        diff_perm = expand(M_ab - perm_val)
        print(f"  M[{a},{b}] - perm(A_{{ba}}) = {diff_perm}")

        # det(A_ba)
        det_ba = expand(det(sub))
        print(f"  det(A_{{ba}}) = {det_ba}")
        print(f"  M[{a},{b}] - det(A_{{ba}}) = {expand(M_ab - det_ba)}")

        # (perm + det) / 2?
        half_sum = expand((perm_val + det_ba) / 2)
        print(f"  (perm + det)/2 = {half_sum}")
        print(f"  M[{a},{b}] = (perm+det)/2? {expand(M_ab - half_sum) == 0}")

        # (perm - det) / 2?
        half_diff = expand((perm_val - det_ba) / 2)
        print(f"  (perm - det)/2 = {half_diff}")
        print(f"  M[{a},{b}] = (perm-det)/2? {expand(M_ab - half_diff) == 0}")


def explore_permanent_structure(n):
    """Check if M is a permanent of some modified matrix."""
    print(f"\n{'='*60}")
    print(f"PERMANENT STRUCTURE: n={n}")
    print(f"{'='*60}")

    r, s = make_symbols(n)
    V = list(range(n))

    A = zeros(n)
    for i in range(n):
        for j in range(n):
            if i != j:
                A[i,j] = r + s[(i,j)]

    # Build the full M matrix
    M = zeros(n)
    for a in range(n):
        for b in range(n):
            if a != b:
                M[a,b] = compute_M(a, b, V, r, s)

    print(f"\n  Transfer matrix M:")
    for a in range(n):
        for b in range(n):
            if a != b:
                print(f"  M[{a},{b}] = {M[a,b]}")

    # Check: is M symmetric?
    print(f"\n  Symmetry check:")
    for a in range(n):
        for b in range(a+1, n):
            d = expand(M[a,b] - M[b,a])
            print(f"  M[{a},{b}] - M[{b},{a}] = {d}")

    # Check: is M even in r?
    print(f"\n  Even-r check:")
    for a in range(n):
        for b in range(n):
            if a != b:
                M_neg = expand(M[a,b].subs(r, -r))
                even = expand(M_neg - M[a,b]) == 0
                print(f"  M[{a},{b}](-r) = M[{a},{b}]? {even}")

    # Check: row sums and column sums
    print(f"\n  Row/column sums:")
    for a in range(n):
        row_sum = sum(M[a,b] for b in range(n) if b != a)
        col_sum = sum(M[b,a] for b in range(n) if b != a)
        print(f"  Row {a}: {expand(row_sum)}")
        print(f"  Col {a}: {expand(col_sum)}")


def main():
    for n in [3, 4]:
        explicit_comparison(n)

    for n in [3, 4]:
        explore_permanent_structure(n)


if __name__ == '__main__':
    main()
