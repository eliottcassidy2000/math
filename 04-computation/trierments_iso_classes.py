#!/usr/bin/env python3
"""
Count unlabeled isomorphism classes of t(r)ienerments, refined by number of
dual-headed (↔) edges.

A t(r)ienerment on n vertices assigns each pair {i,j} one of three states:
  - i→j  (directed)
  - i←j  (directed other way)
  - i↔j  (dual-headed, a "tie")

This generalizes tournaments (which only allow → and ←).

The generating polynomial f(n,x) has:
  [x^k] f(n,x)  =  # unlabeled iso classes with exactly k ↔ edges

Key identities:
  f(n, 0)            = A000568(n)        # pure tournaments
  sum(f(n, x) at x=1) = A007107(n)      # all unlabeled oriented graphs
  [x^C(n,2)] f(n,x) = 1                 # all-↔ graph
  [x^{C(n,2)-1}] f(n,x) = 1            # one directed edge, rest ↔

Burnside framework:
  For σ ∈ S_n, the unordered pairs {i,j} form orbits under σ of two types:
  - Type-A (size d): orbit closes with σ^d(i)=i, σ^d(j)=j.
                     Any of →,←,↔ can be assigned. Factor: (2 + x^d).
  - Type-B (size d): orbit closes with σ^d(i)=j, σ^d(j)=i (swap!).
                     Only ↔ is consistent. Factor: x^d.
                     Arises ONLY from antipodal pairs in even-length cycles.

The "negative t(r)ienerment" space (ties = absent edge rather than ↔) is
isomorphic: the bijection ↔ <--> absent commutes with S_n, so iso-class
counts match. Hence total iso classes = A007107 (unlabeled oriented graphs).

Instance: claude/trierments-dual-arrows-6XsqY, 2026-05-01
"""

from itertools import permutations
from math import factorial


def get_pair_orbits(perm, n):
    """
    Decompose all C(n,2) pairs into type-A and type-B orbits under perm.

    Returns (list of type-A orbit sizes, list of type-B orbit sizes).
    """
    visited = set()
    type_a = []
    type_b = []

    for si in range(n):
        for sj in range(si + 1, n):
            if (si, sj) in visited:
                continue

            # Trace orbit of unordered pair {si, sj}
            orbit_size = 0
            ci, cj = si, sj
            while True:
                key = (min(ci, cj), max(ci, cj))
                if key in visited:
                    break
                visited.add(key)
                orbit_size += 1
                ci, cj = perm[ci], perm[cj]

            # After orbit_size steps from (si, sj), where do we land?
            ci, cj = si, sj
            for _ in range(orbit_size):
                ci, cj = perm[ci], perm[cj]

            if (ci, cj) == (si, sj):
                type_a.append(orbit_size)
            else:  # must be (sj, si)
                type_b.append(orbit_size)

    return type_a, type_b


def poly_multiply_factor(poly, factor_coeffs, max_deg):
    """Multiply polynomial (list of coeffs) by factor_coeffs, cap at max_deg."""
    result = [0] * (max_deg + 1)
    for i, pi in enumerate(poly):
        if pi == 0:
            continue
        for j, fj in enumerate(factor_coeffs):
            if i + j <= max_deg:
                result[i + j] += pi * fj
    return result


def burnside_poly(n):
    """
    Compute f(n, x) via Burnside's lemma over all n! permutations.

    Returns list of integer coefficients [c_0, c_1, ..., c_{C(n,2)}].
    """
    m = n * (n - 1) // 2  # C(n,2) = total edges
    total = [0] * (m + 1)

    for perm in permutations(range(n)):
        type_a, type_b = get_pair_orbits(list(perm), n)

        # Build contribution polynomial: Π_A(2+x^d) × Π_B(x^e)
        contrib = [0] * (m + 1)
        contrib[0] = 1

        for d in type_a:
            factor = [0] * min(d + 1, m + 1)
            factor[0] = 2
            if d <= m:
                factor[d] = 1
            contrib = poly_multiply_factor(contrib, factor, m)

        for e in type_b:
            # multiply by x^e = shift right by e
            new = [0] * (m + 1)
            for k in range(m + 1 - e):
                new[k + e] += contrib[k]
            contrib = new

        for k in range(m + 1):
            total[k] += contrib[k]

    nfact = factorial(n)
    return [c // nfact for c in total]


def format_poly(coeffs):
    terms = []
    for k, c in enumerate(coeffs):
        if c == 0:
            continue
        if k == 0:
            terms.append(str(c))
        elif k == 1:
            terms.append(f"{c}x" if c != 1 else "x")
        else:
            terms.append(f"{c}x^{k}" if c != 1 else f"x^{k}")
    return " + ".join(terms) if terms else "0"


if __name__ == "__main__":
    # A000568: unlabeled tournaments
    A000568 = [1, 1, 1, 2, 4, 12, 56, 456, 6880, 191536]
    # A007107: unlabeled oriented graphs
    A007107 = [1, 1, 2, 7, 42, 582, 21480, 2142288]

    print("T(r)ienerment iso classes by number of ↔ (dual-headed) edges")
    print("=" * 65)
    print()

    for n in range(1, 8):
        poly = burnside_poly(n)
        total = sum(poly)
        m = n * (n - 1) // 2

        print(f"n = {n}  (C({n},2) = {m} edges)")
        print(f"  f(n,x) = {format_poly(poly)}")
        print(f"  Coefficients [x^0 .. x^{m}]: {poly}")
        print(f"  Total iso classes: {total}  (A007107 check: {A007107[n] if n < len(A007107) else '?'})")
        t0 = poly[0]
        expected_t0 = A000568[n] if n < len(A000568) else "?"
        print(f"  f(n,0) = {t0}  (A000568 check: {expected_t0})")
        print(f"  Last two coeffs [x^{m-1}, x^{m}]: {poly[m-1] if m >= 1 else 'N/A'}, {poly[m]}")
        print()
