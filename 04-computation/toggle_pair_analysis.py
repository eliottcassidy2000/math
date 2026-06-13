#!/usr/bin/env python3
"""
Toggle pair analysis — does the toggle on vertex u cancel odd r-powers
at the INDIVIDUAL pair level (S, S Δ {u}), or only in aggregate?

If individual pairs cancel, that gives a direct involution proof.
If only aggregate cancels, need a different structure.

Also: homogeneity analysis. M[a,b] is homogeneous of degree n-2 in
(r, s_ij) jointly. What constraints does this impose?

Instance: opus-2026-03-06-S22
"""

from itertools import permutations, combinations
from sympy import symbols, expand, Poly, Rational, Symbol, collect
from collections import defaultdict

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


def test_individual_toggle_pairs(n, a, b):
    """Check if individual (S, S Δ {u}) pairs cancel for odd r-powers."""
    print(f"\n=== Individual toggle pair analysis: n={n}, a={a}, b={b} ===")

    r, s = make_symbols(n)
    U = [v for v in range(n) if v != a and v != b]

    for u in U:
        print(f"\nToggle vertex u={u}:")
        pair_cancels = True

        # Enumerate pairs (S, S Δ {u}) where u ∈ S
        U_minus_u = [v for v in U if v != u]

        for mask in range(2**len(U_minus_u)):
            S_base = [U_minus_u[i] for i in range(len(U_minus_u)) if (mask >> i) & 1]

            # S contains u
            S1 = set(S_base + [u])
            R1 = set(v for v in U if v not in S1)
            sign1 = (-1)**len(S1)

            # S does not contain u
            S2 = set(S_base)
            R2 = set(v for v in U if v not in S2)
            sign2 = (-1)**len(S2)

            Ea1 = ham_paths_ending_at(S1 | {a}, a, r, s)
            Bb1 = ham_paths_beginning_at(R1 | {b}, b, r, s)
            term1 = expand(sign1 * Ea1 * Bb1)

            Ea2 = ham_paths_ending_at(S2 | {a}, a, r, s)
            Bb2 = ham_paths_beginning_at(R2 | {b}, b, r, s)
            term2 = expand(sign2 * Ea2 * Bb2)

            pair_sum = expand(term1 + term2)

            # Check r-parity
            if pair_sum == 0:
                parity_ok = True
                odd_coeffs = {}
            else:
                p = Poly(pair_sum, r)
                odd_coeffs = {m[0]: c for m, c in p.as_dict().items() if m[0] % 2 == 1 and c != 0}
                parity_ok = len(odd_coeffs) == 0

            if not parity_ok:
                pair_cancels = False
                print(f"  FAIL: S1={sorted(S1)}, S2={sorted(S2)}")
                print(f"    pair_sum r-degrees: {sorted(set(m[0] for m in Poly(pair_sum, r).as_dict()))}")
                for deg, c in sorted(odd_coeffs.items()):
                    print(f"    r^{deg}: {c}")
            else:
                print(f"  OK: S1={sorted(S1)}, S2={sorted(S2)} — pair has even r-powers only")

        if pair_cancels:
            print(f"  ALL pairs cancel for u={u}!")
        else:
            print(f"  Some pairs do NOT individually cancel for u={u}")


def test_homogeneity(n, a, b):
    """Verify M[a,b] is homogeneous of degree n-2."""
    print(f"\n=== Homogeneity check: n={n}, a={a}, b={b} ===")

    r, s = make_symbols(n)
    U = [v for v in range(n) if v != a and v != b]

    M = 0
    for mask in range(2**len(U)):
        S = [U[i] for i in range(len(U)) if (mask >> i) & 1]
        R = [u for u in U if u not in S]
        sign = (-1)**len(S)
        Sa = set(S + [a])
        Rb = set(R + [b])
        Ea = ham_paths_ending_at(Sa, a, r, s)
        Bb = ham_paths_beginning_at(Rb, b, r, s)
        M += sign * Ea * Bb
    M = expand(M)

    # Check: M should be homogeneous of degree n-2 in (r, s_ij)
    # where deg(r) = 1 and deg(s_ij) = 1
    p = Poly(M, r)
    homogeneous = True
    for monom, coeff in p.as_dict().items():
        r_deg = monom[0]
        # coeff is a polynomial in s_ij's
        # Each s_ij has degree 1, so coeff should have s-degree n-2-r_deg
        # Check by substituting s_ij -> t*s_ij and checking coefficient of t^(n-2-r_deg)
        t = Symbol('t')
        coeff_scaled = coeff
        for key in s:
            if key[0] < key[1]:
                coeff_scaled = coeff_scaled.subs(s[key], t * s[key])
        coeff_poly = Poly(expand(coeff_scaled), t)
        degrees = [m[0] for m, c in coeff_poly.as_dict().items() if c != 0]
        if len(set(degrees)) > 1:
            print(f"  NOT homogeneous! r^{r_deg} has s-degrees {sorted(set(degrees))}")
            homogeneous = False
        elif degrees and degrees[0] != n - 2 - r_deg:
            print(f"  WRONG degree! r^{r_deg} has s-degree {degrees[0]}, expected {n-2-r_deg}")
            homogeneous = False
        else:
            print(f"  r^{r_deg}: s-degree = {degrees[0] if degrees else 0}, expected {n-2-r_deg} ✓")

    if homogeneous:
        print(f"  M is homogeneous of degree {n-2} ✓")

    return M


def test_symmetry_direct(n, a, b):
    """Directly compute and compare M[a,b] and M[b,a]."""
    print(f"\n=== Direct symmetry check: n={n}, a={a}, b={b} ===")

    r, s = make_symbols(n)

    def compute_M_entry(aa, bb):
        U = [v for v in range(n) if v != aa and v != bb]
        M = 0
        for mask in range(2**len(U)):
            S = [U[i] for i in range(len(U)) if (mask >> i) & 1]
            R = [u for u in U if u not in S]
            sign = (-1)**len(S)
            Sa = set(S + [aa])
            Rb = set(R + [bb])
            Ea = ham_paths_ending_at(Sa, aa, r, s)
            Bb = ham_paths_beginning_at(Rb, bb, r, s)
            M += sign * Ea * Bb
        return expand(M)

    Mab = compute_M_entry(a, b)
    Mba = compute_M_entry(b, a)

    D = expand(Mab - Mba)
    print(f"  M[{a},{b}] = {Mab}")
    print(f"  M[{b},{a}] = {Mba}")
    print(f"  D = M[{a},{b}] - M[{b},{a}] = {D}")

    if D == 0:
        print(f"  SYMMETRIC ✓")
    else:
        print(f"  NOT symmetric!")
        # Analyze D
        p = Poly(D, r)
        for monom, coeff in sorted(p.as_dict().items()):
            print(f"    r^{monom[0]}: {coeff}")


def test_derivative_at_zero():
    """
    Compute dM[a,b]/dr at r=0 symbolically.
    If M is even in r, this should be 0.

    dM/dr = sum_S (-1)^|S| d/dr[E_a(S+a) B_b(R+b)]
    = sum_S (-1)^|S| [E'_a B_b + E_a B'_b]

    where E'_a = dE_a/dr, etc.

    At r=0: E_a(S+a; 0, s) = sum over Ham paths in T[S+a] of product of s_ij.
    E'_a(S+a; 0, s) = sum over paths of (sum over edges e of product_{e'≠e} s_{e'}).
    """
    print("\n\n=== Derivative dM/dr at r=0 ===")

    n = 5; a = 0; b = 1
    r, s = make_symbols(n)
    U = [v for v in range(n) if v != a and v != b]

    # Compute M at symbolic r, then differentiate
    M = 0
    for mask in range(2**len(U)):
        S = [U[i] for i in range(len(U)) if (mask >> i) & 1]
        R = [u for u in U if u not in S]
        sign = (-1)**len(S)
        Sa = set(S + [a])
        Rb = set(R + [b])
        Ea = ham_paths_ending_at(Sa, a, r, s)
        Bb = ham_paths_beginning_at(Rb, b, r, s)
        M += sign * Ea * Bb
    M = expand(M)

    from sympy import diff
    dMdr = expand(diff(M, r))
    dMdr_at_0 = expand(dMdr.subs(r, 0))

    print(f"  n={n}, dM[{a},{b}]/dr|_{{r=0}} = {dMdr_at_0}")

    # Also check second derivative (should be nonzero if r^2 is present)
    d2Mdr2 = expand(diff(M, r, 2))
    d2Mdr2_at_0 = expand(d2Mdr2.subs(r, 0))
    print(f"  d²M/dr²|_{{r=0}} = {d2Mdr2_at_0}")

    # And third
    d3Mdr3 = expand(diff(M, r, 3))
    d3Mdr3_at_0 = expand(d3Mdr3.subs(r, 0))
    print(f"  d³M/dr³|_{{r=0}} = {d3Mdr3_at_0}")

    return M


def main():
    print("=" * 70)
    print("TOGGLE PAIR ANALYSIS — INDIVIDUAL vs AGGREGATE CANCELLATION")
    print("=" * 70)

    # Test at n=4
    test_individual_toggle_pairs(4, 0, 1)

    # Test at n=5
    test_individual_toggle_pairs(5, 0, 1)

    # Homogeneity
    test_homogeneity(4, 0, 1)
    test_homogeneity(5, 0, 1)

    # Direct symmetry
    test_symmetry_direct(4, 0, 1)
    test_symmetry_direct(5, 0, 1)

    # Derivative
    test_derivative_at_zero()

    print("\n" + "=" * 70)
    print("CONCLUSIONS")
    print("=" * 70)

if __name__ == '__main__':
    main()
