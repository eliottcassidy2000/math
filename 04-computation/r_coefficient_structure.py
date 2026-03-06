#!/usr/bin/env python3
"""
Systematic analysis of the r^k coefficients of M[a,b].

Key observations so far:
- [r^{n-2}] M = (n-2)! if n even, 0 if n odd (from counting argument)
- [r^{n-3}] M = 0 (first odd power vanishes)
- [r^{n-4}] = structured function of s_{uv} with u∈U, v∈{a,b}
- All odd r-powers vanish

This script computes ALL coefficients systematically to find the pattern.

Instance: opus-2026-03-06-S22
"""

from itertools import permutations
from sympy import symbols, expand, Symbol, Poly, diff, factorial, binomial
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

def compute_M(a, b, n, r, s):
    V = list(range(n))
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


def analyze_coefficients(n, a=0, b=1):
    """Full coefficient analysis."""
    print(f"\n{'='*60}")
    print(f"COEFFICIENT STRUCTURE: n={n}, a={a}, b={b}")
    print(f"{'='*60}")

    r, s = make_symbols(n)
    M_ab = compute_M(a, b, n, r, s)

    p = Poly(M_ab, r)
    coeffs = {}
    for monom, coeff in p.as_dict().items():
        coeffs[monom[0]] = coeff

    for k in range(n-1):
        c = coeffs.get(k, 0)
        parity = "ODD" if k % 2 == 1 else "even"
        print(f"\n  [r^{k}] M = {expand(c)}  [{parity}]")

        if c != 0 and k % 2 == 0:
            # Analyze the s-degree structure
            if hasattr(c, 'as_ordered_terms'):
                terms = c.as_ordered_terms()
                s_degrees = set()
                for term in terms:
                    # Count number of s variables
                    s_vars = [v for v in term.free_symbols if str(v).startswith('s')]
                    # Degree in s
                    deg = sum(term.as_powers_dict().get(v, 0) for v in s_vars)
                    s_degrees.add(deg)
                if len(s_degrees) == 1:
                    print(f"    s-degree: {s_degrees.pop()} (homogeneous)")
                else:
                    print(f"    s-degrees: {sorted(s_degrees)}")
            else:
                print(f"    (constant)")

    # Verify the top-degree formula
    m = n - 2
    if m % 2 == 0:
        expected_top = factorial(m)
    else:
        expected_top = 0
    actual_top = coeffs.get(m, 0)
    print(f"\n  [r^{m}] = {actual_top}, expected {expected_top}: {'✓' if expand(actual_top - expected_top) == 0 else '✗'}")


def check_r2_formula(n, a=0, b=1):
    """Check if [r^2] M = 2 * sum_{u in U} (s_{au} + s_{bu}) for general n."""
    print(f"\n{'='*60}")
    print(f"r^2 FORMULA CHECK: n={n}")
    print(f"{'='*60}")

    r, s = make_symbols(n)
    M_ab = compute_M(a, b, n, r, s)

    # Extract r^2 coefficient
    d2M = diff(M_ab, r, 2)
    r2_coeff = expand(d2M.subs(r, 0) / 2)

    U = [v for v in range(n) if v != a and v != b]

    # Candidate: 2 * sum_{u in U} (s_{au} + s_{bu})
    candidate = 2 * sum(s[(a, u)] + s[(b, u)] for u in U)

    # For n=4: U={2,3}, candidate = 2*(s02+s03+s12+s13)
    # But [r^2] = 2 (constant), so this can't be right for n=4.

    # Wait: for n=4, degree n-2=2, so r^2 has s-degree 0.
    # The formula should depend on n.

    # Let me just compare
    diff_val = expand(r2_coeff - candidate)
    print(f"  [r^2] M = {r2_coeff}")
    print(f"  2*sum(s_au+s_bu) = {expand(candidate)}")
    print(f"  difference = {diff_val}")


def count_top_coefficient():
    """Verify: [r^{n-2}] M = (n-2)! * chi(n even) by counting."""
    print(f"\n{'='*60}")
    print("TOP COEFFICIENT VERIFICATION")
    print(f"{'='*60}")

    for n in range(3, 7):
        r, s = make_symbols(n)
        m = n - 2

        # [r^m] of E_a(S+a) = |S|! (number of Ham paths ending at a through S+a)
        # [r^m] of M = sum_k C(m,k) (-1)^k k! (m-k)!
        # = m! * sum_k (-1)^k = m! * (1+(-1)^m)/2 if m >= 0

        formula = factorial(m) * (1 + (-1)**m) // 2

        # Actually compute
        M_ab = compute_M(0, 1, n, r, s)
        p = Poly(M_ab, r)
        actual = p.as_dict().get((m,), 0)

        print(f"  n={n}: [r^{m}] = {actual}, formula = {formula}, match = {expand(actual - formula) == 0}")

        # Also the next coefficient [r^{m-1}]:
        if m >= 1:
            r_m1 = p.as_dict().get((m-1,), 0)
            print(f"         [r^{m-1}] = {r_m1}")


def per_subset_structure(n, a=0, b=1):
    """Show the per-subset [r^k] contributions for all k."""
    print(f"\n{'='*60}")
    print(f"PER-SUBSET COEFFICIENTS: n={n}")
    print(f"{'='*60}")

    r, s = make_symbols(n)
    U = [v for v in range(n) if v != a and v != b]

    for mask in range(2**len(U)):
        S = [U[i] for i in range(len(U)) if (mask >> i) & 1]
        R = [u for u in U if u not in S]
        sign = (-1)**len(S)
        Sa = set(S + [a])
        Rb = set(R + [b])
        Ea = ham_paths_ending_at(Sa, a, r, s)
        Bb = ham_paths_beginning_at(Rb, b, r, s)
        term = expand(sign * Ea * Bb)

        p = Poly(term, r) if term != 0 else None
        print(f"\n  S={sorted(S)} (sign={sign:+d}):")
        if p:
            for monom, coeff in sorted(p.as_dict().items()):
                print(f"    [r^{monom[0]}]: {coeff}")
        else:
            print(f"    0")


def main():
    for n in [3, 4, 5]:
        analyze_coefficients(n)

    count_top_coefficient()

    for n in [4, 5]:
        check_r2_formula(n)

    for n in [4, 5]:
        per_subset_structure(n)


if __name__ == '__main__':
    main()
