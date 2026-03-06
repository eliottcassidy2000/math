#!/usr/bin/env python3
"""
Test whether M[a,b] equals a cofactor or minor of the adjacency matrix A = rJ' + S.

KEY INSIGHT: A(-r) = -A^T (since J' is symmetric, S is skew-symmetric).
If M[a,b] = cofactor_{ba}(A) or similar, then:
  M[a,b](-r) = cofactor_{ba}(-A^T) = (-1)^{n-3} cofactor_{ab}(A^T)
             = (-1)^{n-3} cofactor_{ab}(A)  [since det is transpose-invariant]
Combined with M[b,a] = cofactor_{ab}(A), this would prove both symmetry and even r-powers.

Instance: opus-2026-03-06-S22
"""

from itertools import permutations
from sympy import symbols, expand, Symbol, Matrix, det, zeros, ones

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
    """Compute M[a,b] on given vertex set."""
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

def build_adjacency(n, r, s):
    """Build A = rJ' + S where J' = J - I, S skew-symmetric."""
    A = zeros(n)
    for i in range(n):
        for j in range(n):
            if i != j:
                A[i,j] = r + s[(i,j)]
    return A

def test_cofactor_identity(n):
    """Test if M[a,b] equals some cofactor/minor of A."""
    print(f"\n{'='*60}")
    print(f"COFACTOR IDENTITY TEST: n={n}")
    print(f"{'='*60}")

    r, s = make_symbols(n)
    V = list(range(n))
    A = build_adjacency(n, r, s)

    print(f"\n  A = rJ' + S, size {n}x{n}")

    # Test all (a,b) pairs for small n
    pairs = [(0,1), (1,0)]
    if n <= 4:
        pairs = [(i,j) for i in range(n) for j in range(n) if i != j]

    for a, b in pairs:
        M_ab = compute_M(a, b, V, r, s)

        # Test 1: M[a,b] = (-1)^{a+b} det(A with row b, col a deleted)?
        # This is the standard cofactor C_{ba}
        A_del = A.minor_submatrix(b, a)
        cof_ba = expand((-1)**(a+b) * det(A_del))
        diff1 = expand(M_ab - cof_ba)

        # Test 2: M[a,b] = (-1)^{a+b} det(A with row a, col b deleted)?
        A_del2 = A.minor_submatrix(a, b)
        cof_ab = expand((-1)**(a+b) * det(A_del2))
        diff2 = expand(M_ab - cof_ab)

        # Test 3: Just the minor (no sign)
        minor_ba = expand(det(A.minor_submatrix(b, a)))
        diff3 = expand(M_ab - minor_ba)

        minor_ab = expand(det(A.minor_submatrix(a, b)))
        diff4 = expand(M_ab - minor_ab)

        # Test 5: M[a,b] = det(A with row a, col a deleted)? (diagonal minor)
        if a == b:
            continue
        minor_aa = expand(det(A.minor_submatrix(a, a)))
        diff5 = expand(M_ab - minor_aa)

        print(f"\n  ({a},{b}):")
        print(f"    M[{a},{b}] = cofactor_ba? {diff1 == 0}")
        print(f"    M[{a},{b}] = cofactor_ab? {diff2 == 0}")
        print(f"    M[{a},{b}] = minor_ba?    {diff3 == 0}")
        print(f"    M[{a},{b}] = minor_ab?    {diff4 == 0}")
        print(f"    M[{a},{b}] = minor_aa?    {diff5 == 0}")

        if diff1 != 0 and diff2 != 0 and diff3 != 0 and diff4 != 0:
            # Try scalar multiples
            if M_ab != 0:
                for c in [1, -1, 2, -2]:
                    if expand(c * M_ab - minor_ba) == 0:
                        print(f"    minor_ba = {c} * M[{a},{b}]!")
                    if expand(c * M_ab - minor_ab) == 0:
                        print(f"    minor_ab = {c} * M[{a},{b}]!")
                    if expand(c * M_ab - cof_ba) == 0:
                        print(f"    cof_ba = {c} * M[{a},{b}]!")
                    if expand(c * M_ab - cof_ab) == 0:
                        print(f"    cof_ab = {c} * M[{a},{b}]!")


def test_modified_matrices(n):
    """Test cofactors of modified matrices: I-A, I+A, lambda*I - A, etc."""
    print(f"\n{'='*60}")
    print(f"MODIFIED MATRIX TEST: n={n}")
    print(f"{'='*60}")

    r, s = make_symbols(n)
    V = list(range(n))
    A = build_adjacency(n, r, s)
    I_n = Matrix.eye(n)
    J = ones(n, n)

    # Various candidate matrices
    candidates = {
        "I - A": I_n - A,
        "I + A": I_n + A,
        "A - I": A - I_n,
        "J - A": J - A,
        "A": A,
    }

    a, b = 0, 1
    M_ab = compute_M(a, b, V, r, s)
    M_ba = compute_M(b, a, V, r, s)

    print(f"\n  Reference: M[{a},{b}] computed directly")
    print(f"  M[{a},{b}] = M[{b},{a}]? {expand(M_ab - M_ba) == 0}")

    for name, mat in candidates.items():
        minor_ba = expand(det(mat.minor_submatrix(b, a)))
        minor_ab = expand(det(mat.minor_submatrix(a, b)))
        cof_ba = expand((-1)**(a+b) * minor_ba)
        cof_ab = expand((-1)**(a+b) * minor_ab)

        match_minor_ba = expand(M_ab - minor_ba) == 0
        match_minor_ab = expand(M_ab - minor_ab) == 0
        match_cof_ba = expand(M_ab - cof_ba) == 0
        match_cof_ab = expand(M_ab - cof_ab) == 0

        if any([match_minor_ba, match_minor_ab, match_cof_ba, match_cof_ab]):
            print(f"\n  *** MATCH with {name}:")
            if match_minor_ba: print(f"      M[{a},{b}] = minor({b},{a})")
            if match_minor_ab: print(f"      M[{a},{b}] = minor({a},{b})")
            if match_cof_ba:   print(f"      M[{a},{b}] = cofactor({b},{a})")
            if match_cof_ab:   print(f"      M[{a},{b}] = cofactor({a},{b})")

        # Also try ratios
        if minor_ba != 0 and M_ab != 0:
            ratio = expand(minor_ba / M_ab)
            # Check if ratio is a simple constant
            try:
                from sympy import Rational
                if ratio.is_number:
                    print(f"  {name}: minor_ba / M = {ratio}")
            except:
                pass


def test_det_full(n):
    """Test if det(A) relates to something useful."""
    print(f"\n{'='*60}")
    print(f"FULL DETERMINANT: n={n}")
    print(f"{'='*60}")

    r, s = make_symbols(n)
    A = build_adjacency(n, r, s)

    d = expand(det(A))
    print(f"  det(A) = {d}")

    # Check: det(A) should satisfy det(A(-r)) = (-1)^n det(A^T) = (-1)^n det(A)
    d_neg = expand(d.subs(r, -r))
    print(f"  det(A(-r)) = {d_neg}")
    print(f"  (-1)^n det(A) = {expand((-1)**n * d)}")
    print(f"  det(A(-r)) = (-1)^n det(A)? {expand(d_neg - (-1)**n * d) == 0}")


def test_adjugate_identity(n):
    """Test if M[a,b] relates to entries of adj(A) or (A^{-1} * det(A))."""
    print(f"\n{'='*60}")
    print(f"ADJUGATE IDENTITY TEST: n={n}")
    print(f"{'='*60}")

    r, s = make_symbols(n)
    V = list(range(n))
    A = build_adjacency(n, r, s)

    # adj(A)[i,j] = (-1)^{i+j} det(A with row j, col i deleted) = cofactor(j,i)
    # So adj(A)[i,j] = cofactor_{ji}(A)

    a, b = 0, 1
    M_ab = compute_M(a, b, V, r, s)

    # adj(A)[a,b] = cofactor_{ba}(A)
    adj_ab = expand((-1)**(a+b) * det(A.minor_submatrix(b, a)))
    # adj(A)[b,a] = cofactor_{ab}(A)
    adj_ba = expand((-1)**(a+b) * det(A.minor_submatrix(a, b)))

    print(f"  M[{a},{b}] = adj(A)[{a},{b}]? {expand(M_ab - adj_ab) == 0}")
    print(f"  M[{a},{b}] = adj(A)[{b},{a}]? {expand(M_ab - adj_ba) == 0}")

    # If M[a,b] = adj(A)[b,a], then M[a,b] = cofactor_{ab}(A)
    # And symmetry of M follows from: A(-r) = -A^T implies
    # cofactor_{ab}(A(-r)) = cofactor_{ab}(-A^T) = (-1)^{n-3} cofactor_{ba}(A)
    # So M[a,b](-r) = (-1)^{n-3} M[b,a]
    # If n is odd: M[a,b](-r) = M[b,a], combined with M[b,a](-r) = M[a,b]
    #   gives M[a,b](-r,-r) = M[a,b], i.e. M[a,b](r) = M[a,b](r) trivially
    # Need more...


def main():
    for n in [3, 4]:
        test_cofactor_identity(n)

    for n in [3, 4]:
        test_modified_matrices(n)

    for n in [3, 4]:
        test_det_full(n)

    test_adjugate_identity(3)
    test_adjugate_identity(4)

    print(f"\n\n{'='*60}")
    print("ANALYSIS")
    print(f"{'='*60}")
    print("""
If M[a,b] = cofactor_{ba}(A) for A = rJ' + S:
  - A(-r) = -A^T  (since J' symmetric, S skew-symmetric)
  - cofactor_{ba}(-A^T) = (-1)^{n-3} cofactor_{ab}(A^T) = (-1)^{n-3} cofactor_{ab}(A)
  - So M[a,b](-r) = (-1)^{n-3} M[b,a]
  - This means M is even in r IFF M[a,b] = (-1)^{n-3} M[b,a]
  - For n odd: M[a,b] = M[b,a] (symmetry) AND M even in r
  - For n even: M[a,b] = -M[b,a] (antisymmetry)... but we KNOW M[a,b] = M[b,a]!

  So simple cofactors of A won't work for even n.

  Need a different matrix or a different relationship.
""")


if __name__ == '__main__':
    main()
