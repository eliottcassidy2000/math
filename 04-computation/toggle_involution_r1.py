#!/usr/bin/env python3
"""
TOGGLE INVOLUTION FOR r^1 VANISHING

At n=4, the r^1 coefficient of M[a,b] vanishes because each s_ij term
appears exactly twice with opposite signs. The two appearances differ
by toggling a single element of S.

This script checks whether this pattern holds for larger n:
- For each s-monomial at r^1, find the toggle partner
- Check that the toggle is a perfect matching (involution)

kind-pasteur-2026-03-06-S23
"""
from itertools import permutations
from sympy import symbols, expand, Poly
from collections import defaultdict

def setup(n):
    r = symbols('r')
    sv = {}
    for i in range(n):
        for j in range(i+1, n):
            sv[(i,j)] = symbols(f's{i}{j}')
    def s(i, j):
        if i == j: return 0
        if i < j: return sv[(i,j)]
        return -sv[(j,i)]
    def t(i, j):
        if i == j: return 0
        return r + s(i, j)
    return r, sv, s, t

def transfer_M(t_fn, n, a, b):
    U = [v for v in range(n) if v != a and v != b]
    result = 0
    for mask in range(1 << len(U)):
        S = [U[i] for i in range(len(U)) if mask & (1 << i)]
        R = [U[i] for i in range(len(U)) if not (mask & (1 << i))]
        sign = (-1)**len(S)

        S_set = set(S) | {a}
        R_set = set(R) | {b}

        ea = 0
        for p in permutations(sorted(S_set)):
            if p[-1] != a: continue
            prod = 1
            for i in range(len(p)-1):
                prod *= t_fn(p[i], p[i+1])
            ea += prod

        if len(S_set) == 1:
            ea = 1

        bb = 0
        for p in permutations(sorted(R_set)):
            if p[0] != b: continue
            prod = 1
            for i in range(len(p)-1):
                prod *= t_fn(p[i], p[i+1])
            bb += prod

        if len(R_set) == 1:
            bb = 1

        result += sign * ea * bb
    return expand(result)

print("=" * 70)
print("TOGGLE INVOLUTION ANALYSIS FOR r^1 VANISHING")
print("=" * 70)

for n in [4, 5, 6]:
    r, sv, s, t = setup(n)
    a, b = 0, 1

    M_ab = transfer_M(t, n, a, b)
    p_r = Poly(M_ab, r)

    print(f"\nn={n}:")

    for k in range(p_r.degree() + 1):
        coeff = expand(p_r.nth(k))
        parity = "EVEN" if k % 2 == 0 else "ODD"
        if coeff == 0:
            print(f"  r^{k} [{parity}]: 0")
        else:
            s_vars = list(sv.values())
            if s_vars:
                try:
                    ps = Poly(coeff, *s_vars)
                    nterms = len(ps.as_dict())
                    print(f"  r^{k} [{parity}]: {nterms} s-monomials")
                except:
                    print(f"  r^{k} [{parity}]: {coeff}")
            else:
                print(f"  r^{k} [{parity}]: {coeff}")

    # Detailed analysis of r^1 contributions grouped by S
    print(f"\n  --- Detailed r^1 analysis ---")
    U = [v for v in range(n) if v != a and v != b]

    r1_by_S = {}
    for mask in range(1 << len(U)):
        S = tuple(sorted([U[i] for i in range(len(U)) if mask & (1 << i)]))
        R = tuple(sorted([U[i] for i in range(len(U)) if not (mask & (1 << i))]))
        sign = (-1)**len(S)

        S_set = set(S) | {a}
        R_set = set(R) | {b}

        ea = 0
        for p in permutations(sorted(S_set)):
            if p[-1] != a: continue
            prod = 1
            for i in range(len(p)-1):
                prod *= t(p[i], p[i+1])
            ea += prod
        if len(S_set) == 1: ea = 1

        bb = 0
        for p in permutations(sorted(R_set)):
            if p[0] != b: continue
            prod = 1
            for i in range(len(p)-1):
                prod *= t(p[i], p[i+1])
            bb += prod
        if len(R_set) == 1: bb = 1

        product = expand(sign * ea * bb)
        pr = Poly(product, r) if product != 0 else None
        r1 = expand(pr.nth(1)) if pr and pr.degree() >= 1 else 0

        if r1 != 0:
            r1_by_S[S] = r1

    # Print contributions
    for S in sorted(r1_by_S.keys(), key=lambda x: (len(x), x)):
        val = r1_by_S[S]
        print(f"    S={set(S)}: {val}")

    total = expand(sum(r1_by_S.values()))
    print(f"    TOTAL r^1 = {total}")

    # Look for toggle pairs: S and S delta {v} with canceling contributions
    print(f"\n  --- Toggle pair search ---")
    S_list = list(r1_by_S.keys())
    used = set()
    pairs_found = 0

    for S in S_list:
        if S in used: continue
        S_set = set(S)
        found_toggle = False
        for v in U:
            S2_set = S_set.symmetric_difference({v})
            S2 = tuple(sorted(S2_set))
            if S2 in r1_by_S and S2 not in used and S2 != S:
                # Check if they cancel
                if expand(r1_by_S[S] + r1_by_S[S2]) == 0:
                    print(f"    PAIR: S={S_set} <-> S={S2_set} (toggle vertex {v})")
                    used.add(S)
                    used.add(S2)
                    pairs_found += 1
                    found_toggle = True
                    break
        if not found_toggle and S not in used:
            # Check if contribution is itself 0
            if r1_by_S[S] == 0:
                used.add(S)
            else:
                print(f"    UNPAIRED: S={S_set}, r^1 = {r1_by_S[S]}")

    unpaired = [S for S in S_list if S not in used]
    if unpaired:
        print(f"    {len(unpaired)} unpaired sets remaining")
    else:
        print(f"    ALL PAIRED via single-vertex toggles! ({pairs_found} pairs)")

    # Also check r^3 for n >= 6
    if n >= 6:
        r3 = expand(p_r.nth(3)) if p_r.degree() >= 3 else 0
        print(f"\n  r^3 = {r3}")

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print("""
If ALL r^1 contributions pair via single-vertex toggles, this gives
a SIGN-REVERSING INVOLUTION proof that the r^1 coefficient vanishes.

The toggle involution: for a term from subset S contributing to r^1,
there exists a unique vertex v such that the term from S delta {v}
has the opposite contribution, and these two cancel.

OPEN: Does this extend to r^3 and higher odd powers?
""")
