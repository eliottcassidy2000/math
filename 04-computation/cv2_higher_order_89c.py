#!/usr/bin/env python3
"""
cv2_higher_order_89c.py — Exact |S|=2k contributions to E[∏(1+Z_j)]
opus-2026-03-14-S89c

Found: |S|=2 contrib = 2(n-2)/(n)_2
       |S|=4 contrib = 2(n-4)²/(n)_4

Goal: find general formula for |S|=2k contribution.
"""

from fractions import Fraction
from itertools import permutations, combinations
from math import factorial
from functools import reduce

def falling_factorial(n, k):
    """n(n-1)(n-2)...(n-k+1)"""
    return reduce(lambda a,b: a*b, range(n, n-k, -1), 1)

def compute_contributions_by_size(n):
    """Compute |S|=k contribution for each k."""
    N = factorial(n)
    m = n - 1  # positions 0..m-1

    # Precompute Z values for all permutations
    all_Z = []
    for perm in permutations(range(n)):
        Z = []
        for j in range(m):
            xj = 1 if perm[j+1] == perm[j] + 1 else 0
            yj = 1 if perm[j+1] == perm[j] - 1 else 0
            Z.append(xj - yj)
        all_Z.append(Z)

    contributions = {}
    for size in range(0, m+1, 2):  # Only even sizes (odd vanish)
        total = Fraction(0)
        count = 0
        for S in combinations(range(m), size):
            moment_sum = 0
            for Z in all_Z:
                prod = 1
                for j in S:
                    prod *= Z[j]
                moment_sum += prod
            total += Fraction(moment_sum, N)
            count += 1
        contributions[size] = total

    return contributions

print("="*70)
print("EXACT |S|=2k CONTRIBUTIONS AND CLOSED FORMS")
print("="*70)

all_contribs = {}
for n in range(3, 11):
    print(f"\nn={n}:", flush=True)
    contribs = compute_contributions_by_size(n)
    all_contribs[n] = contribs

    for k in sorted(contribs):
        if contribs[k] != 0:
            ff = falling_factorial(n, k)
            scaled = contribs[k] * ff / 2
            print(f"  |S|={k}: {contribs[k]} = {float(contribs[k]):.12f}")
            print(f"         × (n)_{k}/2 = {scaled}")

# Now extract the "g" values: |S|=2k = 2·g_k(n) / (n)_{2k}
print("\n" + "="*70)
print("SCALED CONTRIBUTIONS: g_k(n) = |S|=2k × (n)_{2k} / 2")
print("="*70)

for k in range(1, 5):
    print(f"\nk={k} (|S|={2*k}):")
    vals = []
    for n in range(2*k+1, 11):
        if 2*k in all_contribs.get(n, {}):
            c = all_contribs[n][2*k]
            ff = falling_factorial(n, 2*k)
            g = c * ff / 2
            m = n - 2*k  # offset
            vals.append((n, m, g))
            print(f"  n={n}, m=n-{2*k}={m}: g = {g}")

    # Try to find pattern in g values
    if len(vals) >= 3:
        print(f"\n  Pattern analysis for g_{k}(m):")
        gs = [(m, g) for n, m, g in vals]

        # Check if it's a polynomial in m
        for deg in range(1, len(gs)):
            # Fit polynomial of degree deg through first deg+1 points
            # Using Lagrange interpolation with exact fractions
            from sympy import symbols, interpolate, simplify, factor, Rational

            x = symbols('x')
            points = [(Rational(m), Rational(g.numerator, g.denominator)) for m, g in gs[:deg+1]]
            try:
                poly = interpolate(points, x)
                poly_simplified = simplify(poly)
                poly_factored = factor(poly)

                # Check remaining points
                all_match = True
                for m, g in gs[deg+1:]:
                    pred = poly.subs(x, m)
                    if pred != Rational(g.numerator, g.denominator):
                        all_match = False
                        break

                if all_match and len(gs) > deg+1:
                    print(f"    Degree {deg} polynomial EXACT: g_{k}(m) = {poly_factored}")
                    break
                elif len(gs) == deg+1:
                    print(f"    Degree {deg} through all points: g_{k}(m) = {poly_factored}")
            except:
                pass

print("\n" + "="*70)
print("SUMMARY TABLE")
print("="*70)

print(f"\n{'n':>3}", end="")
for k in range(1, 5):
    print(f"  {'|S|='+str(2*k):>15}", end="")
print(f"  {'Sum':>15}")

for n in range(3, 11):
    print(f"{n:>3}", end="")
    total = Fraction(0)
    for k in range(1, 5):
        c = all_contribs.get(n, {}).get(2*k, Fraction(0))
        total += c
        print(f"  {float(c):>15.10f}", end="")
    total += 1  # |S|=0
    print(f"  {float(total):>15.10f}")

# Check the closed form: |S|=2 = 2(n-2)/(n)_2, |S|=4 = 2(n-4)^2/(n)_4
print("\nVerification of closed forms:")
for n in range(3, 11):
    # |S|=2
    pred2 = Fraction(2*(n-2), falling_factorial(n, 2))
    act2 = all_contribs[n].get(2, Fraction(0))
    print(f"  n={n}: |S|=2: pred={pred2}, actual={act2}, {'✓' if pred2==act2 else '✗'}")

    # |S|=4
    if n >= 5:
        pred4 = Fraction(2*(n-4)**2, falling_factorial(n, 4))
        act4 = all_contribs[n].get(4, Fraction(0))
        print(f"         |S|=4: pred={pred4}, actual={act4}, {'✓' if pred4==act4 else '✗'}")

print("\nDone!")
