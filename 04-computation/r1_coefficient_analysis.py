#!/usr/bin/env python3
"""
Analyze the coefficient of r^1 in M[a,b] — why does it vanish?

For M[a,b] = sum_S (-1)^|S| E_a(S+a) B_b(R+b), the coefficient of r^1
comes from terms where exactly one arc contributes r and the rest contribute s.

The coefficient of r^1 in M[a,b] is:
[r^1] M = sum_S (-1)^|S| sum_{edge e} [product of s on other edges]

where the sum is over all 2-path-covers and all arcs e in the cover.

If we can show this vanishes, that's the base case (first odd power = 0).
Then r^3 = 0 would follow from a similar but higher-order argument.

Instance: opus-2026-03-06-S22
"""

from itertools import permutations
from sympy import symbols, expand, Symbol, Poly, diff

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


def analyze_r1_coefficient(n, a, b):
    """Extract and analyze the coefficient of r^1 in M[a,b]."""
    print(f"\n{'='*60}")
    print(f"r^1 COEFFICIENT ANALYSIS: n={n}, a={a}, b={b}")
    print(f"{'='*60}")

    r, s = make_symbols(n)
    M_ab = compute_M(a, b, n, r, s)

    # Extract r^1 coefficient
    dM = diff(M_ab, r)
    r1_coeff = expand(dM.subs(r, 0))
    print(f"\n  [r^1] M[{a},{b}] = {r1_coeff}")

    if r1_coeff == 0:
        print(f"  ZERO — r^1 vanishes ✓")
    else:
        print(f"  NONZERO! r^1 does NOT vanish")

    # Now decompose by subset S to see where the cancellation happens
    U = [v for v in range(n) if v != a and v != b]

    print(f"\n  Per-subset contributions to [r^1]:")
    for mask in range(2**len(U)):
        S = [U[i] for i in range(len(U)) if (mask >> i) & 1]
        R = [u for u in U if u not in S]
        sign = (-1)**len(S)
        Sa = set(S + [a])
        Rb = set(R + [b])
        Ea = ham_paths_ending_at(Sa, a, r, s)
        Bb = ham_paths_beginning_at(Rb, b, r, s)
        term = expand(sign * Ea * Bb)

        # r^1 coefficient of this term
        d_term = diff(term, r)
        r1_term = expand(d_term.subs(r, 0))

        print(f"  S={sorted(S)}: [r^1] = {r1_term}")

    # Now: what IS the r^1 coefficient structurally?
    # For each 2-path-cover (P_a on S+a, P_b on R+b), the weight is
    # product_{e in P_a ∪ P_b} (r + s_e).
    # The r^1 coefficient is: sum_e product_{e' != e} s_{e'} summed over covers.
    # = (d/dr product)|_{r=0}
    # = sum_e prod_{e'!=e} s_{e'}
    #
    # This is the "sum of partial products" — for each arc in the cover,
    # replace that arc's s_e by 1 (since it contributes r, not s).

    # At r=0: E_a(S+a; 0, s) = sum of products of s along ham paths ending at a
    # d/dr E_a(S+a; r, s)|_{r=0} = sum over paths of (sum over edge positions k of
    #   product_{j != k} s_{e_j})

    # This is hard to simplify directly. Let me instead look at the STRUCTURE
    # of the r^1 coefficient for small n.

    if n <= 5:
        # Substitute specific s values to check if r1 is "accidentally" zero
        # or structurally zero
        test_vals = {}
        for i in range(n):
            for j in range(i+1, n):
                test_vals[s[(i,j)]] = i + j  # arbitrary distinct values

        r1_num = r1_coeff.subs(test_vals)
        print(f"\n  Numerical test: [r^1] with s_ij = i+j: {r1_num}")


def analyze_all_odd_coefficients(n, a, b):
    """Extract all odd r-power coefficients."""
    print(f"\n{'='*60}")
    print(f"ALL ODD COEFFICIENTS: n={n}, a={a}, b={b}")
    print(f"{'='*60}")

    r, s = make_symbols(n)
    M_ab = compute_M(a, b, n, r, s)

    p = Poly(M_ab, r)
    print(f"\n  M[{a},{b}] as polynomial in r:")
    for monom, coeff in sorted(p.as_dict().items()):
        parity = "ODD" if monom[0] % 2 == 1 else "even"
        print(f"    r^{monom[0]}: {coeff}  [{parity}]")


def trace_r1_cancellation_n4():
    """For n=4, trace exactly how r^1 cancels between subsets."""
    print(f"\n{'='*60}")
    print(f"DETAILED r^1 CANCELLATION: n=4, a=0, b=1")
    print(f"{'='*60}")

    n, a, b = 4, 0, 1
    r, s = make_symbols(n)
    U = [2, 3]

    # For each subset S, compute [r^1] of the contribution
    # and break it into path-level contributions
    for mask in range(4):
        S = [U[i] for i in range(2) if (mask >> i) & 1]
        R = [u for u in U if u not in S]
        sign = (-1)**len(S)
        Sa = list(set(S + [a]))
        Rb = list(set(R + [b]))

        print(f"\n  S={S}, R={R}, sign={sign}")
        print(f"    Paths ending at {a} through {sorted(Sa)}:")

        # Enumerate Hamiltonian paths ending at a in Sa
        others_a = [v for v in Sa if v != a]
        if len(others_a) == 0:
            paths_a = [([a], 1)]
        else:
            paths_a = []
            for perm in permutations(others_a):
                path = list(perm) + [a]
                w_at_0 = 1
                for k in range(len(path)-1):
                    w_at_0 *= s[(path[k], path[k+1])]
                dw = 0
                for pos in range(len(path)-1):
                    partial = 1
                    for k in range(len(path)-1):
                        if k == pos:
                            partial *= 1  # the r contribution
                        else:
                            partial *= s[(path[k], path[k+1])]
                    dw += partial
                paths_a.append((path, expand(w_at_0), expand(dw)))
                print(f"      {' -> '.join(map(str, path))}: w(0)={w_at_0}, dw/dr(0)={dw}")

        print(f"    Paths starting at {b} through {sorted(Rb)}:")
        others_b = [v for v in Rb if v != b]
        if len(others_b) == 0:
            paths_b = [([b], 1)]
        else:
            paths_b = []
            for perm in permutations(others_b):
                path = [b] + list(perm)
                w_at_0 = 1
                for k in range(len(path)-1):
                    w_at_0 *= s[(path[k], path[k+1])]
                dw = 0
                for pos in range(len(path)-1):
                    partial = 1
                    for k in range(len(path)-1):
                        if k == pos:
                            partial *= 1
                        else:
                            partial *= s[(path[k], path[k+1])]
                    dw += partial
                paths_b.append((path, expand(w_at_0), expand(dw)))
                print(f"      {' -> '.join(map(str, path))}: w(0)={w_at_0}, dw/dr(0)={dw}")

    # The r^1 coefficient of sign * E_a * B_b is:
    # sign * [E'_a(0) * B_b(0) + E_a(0) * B'_b(0)]
    # where primes are d/dr.

    print(f"\n  [r^1] contributions by Leibniz rule:")
    total_r1 = 0
    for mask in range(4):
        S = [U[i] for i in range(2) if (mask >> i) & 1]
        R = [u for u in U if u not in S]
        sign = (-1)**len(S)
        Sa = set(S + [a])
        Rb = set(R + [b])

        Ea = ham_paths_ending_at(Sa, a, r, s)
        Bb = ham_paths_beginning_at(Rb, b, r, s)

        Ea_0 = expand(Ea.subs(r, 0))
        Bb_0 = expand(Bb.subs(r, 0))
        dEa = expand(diff(Ea, r).subs(r, 0))
        dBb = expand(diff(Bb, r).subs(r, 0))

        r1_contrib = expand(sign * (dEa * Bb_0 + Ea_0 * dBb))
        total_r1 = expand(total_r1 + r1_contrib)

        print(f"  S={S}: sign*[E'B + EB'] = {sign}*[{dEa}*{Bb_0} + {Ea_0}*{dBb}] = {r1_contrib}")

    print(f"\n  Total [r^1] = {total_r1}")


def main():
    for n in [3, 4, 5]:
        analyze_r1_coefficient(n, 0, 1)

    for n in [3, 4, 5]:
        analyze_all_odd_coefficients(n, 0, 1)

    trace_r1_cancellation_n4()


if __name__ == '__main__':
    main()
