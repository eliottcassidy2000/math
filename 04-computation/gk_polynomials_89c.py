#!/usr/bin/env python3
"""
gk_polynomials_89c.py — Study the g_k polynomials in CV² expansion
opus-2026-03-14-S89c

|S|=2k contribution = 2·g_k(n-2k) / (n)_{2k}

Known: g_1(m) = m, g_2(m) = m², g_3(m) = m(2m²+1)/3

Goal: find g_4, g_5, and identify the general pattern.
Need n=11 or 12 for g_4 data. n=11 has 11! = 40M perms — feasible but slow.
"""

from fractions import Fraction
from itertools import permutations, combinations
from math import factorial
from functools import reduce
from sympy import symbols, interpolate, Rational, factor, simplify, expand

def falling_factorial(n, k):
    return reduce(lambda a,b: a*b, range(n, n-k, -1), 1)

def compute_gk(n_max):
    """Compute g_k values for all accessible n."""
    all_gk = {}

    for n in range(3, n_max + 1):
        print(f"Computing n={n} ({factorial(n)} perms)...", flush=True)
        N = factorial(n)
        m = n - 1

        # Precompute Z for all permutations
        all_Z = []
        for perm in permutations(range(n)):
            Z = []
            for j in range(m):
                xj = 1 if perm[j+1] == perm[j] + 1 else 0
                yj = 1 if perm[j+1] == perm[j] - 1 else 0
                Z.append(xj - yj)
            all_Z.append(Z)

        for size in range(0, m+1, 2):
            if size == 0:
                continue
            total = Fraction(0)
            for S in combinations(range(m), size):
                moment_sum = 0
                for Z in all_Z:
                    prod = 1
                    for j in S:
                        prod *= Z[j]
                    moment_sum += prod
                total += Fraction(moment_sum, N)

            if total != 0:
                k = size // 2
                ff = falling_factorial(n, size)
                g = total * ff / 2
                if k not in all_gk:
                    all_gk[k] = []
                all_gk[k].append((n, n - size, g))  # (n, m=n-2k, g_k(m))

        print(f"  Done n={n}", flush=True)

    return all_gk

# Compute up to n=10 (fast)
all_gk = compute_gk(10)

print("\n" + "="*70)
print("g_k POLYNOMIAL ANALYSIS")
print("="*70)

x = symbols('x')

for k in sorted(all_gk):
    data = all_gk[k]
    print(f"\nk={k}: g_{k}(m)")
    for n, m, g in data:
        print(f"  m={m}: g = {g}")

    ms = [(Rational(m), Rational(g.numerator, g.denominator)) for n, m, g in data]

    # Find minimal degree polynomial fitting all points
    for deg in range(1, len(ms)):
        poly = interpolate(ms[:deg+1], x)
        all_match = True
        for m_val, g_val in ms[deg+1:]:
            if poly.subs(x, m_val) != g_val:
                all_match = False
                break
        if all_match:
            poly_f = factor(poly)
            poly_e = expand(poly)
            print(f"  → g_{k}(m) = {poly_f} = {poly_e}")
            break
    else:
        # Just show the interpolating polynomial through all points
        poly = interpolate(ms, x)
        print(f"  → g_{k}(m) = {factor(poly)} (through all {len(ms)} points)")

# Now try to extend to n=11 for one more g_4 data point
print("\n" + "="*70)
print("COMPUTING n=11 (this may take a few minutes)")
print("="*70)

# n=11 has 39916800 perms — too slow to enumerate
# Instead, use sampling or exploit the domino structure

# Key insight: only domino subsets contribute.
# For |S|=8 at n=11: domino tilings of positions [0..9] with 4 dominoes.
# Number of such tilings: C(10-4, 4) = C(6,4) = 15.
# For each tiling, compute E[Z_{d1}Z_{d1+1}...Z_{d4}Z_{d4+1}].

# But computing each expectation requires iterating over all 11! perms...
# unless we can find a pattern.

# Instead, let's compute g_4 by computing only the |S|=8 contribution at n=11.
# But we need ALL 11! permutations for the exact value.

# Alternative: compute W(n=11) from C code output and use the known
# |S|=2,4,6 contributions to extract |S|=8.

print("Using W(n) values to extract g_4:")
# From nud_weight.c: W(11) = 46963358
# W/n! = 46963358/39916800
W11 = 46963358
nf11 = factorial(11)
cv2_11 = Fraction(W11, nf11) - 1

# |S|=2 at n=11 = 2(n-2)/(n(n-1)) = 2·9/(11·10) = 18/110 = 9/55
s2_11 = Fraction(2*(11-2), 11*10)

# |S|=4 at n=11 = 2(n-4)²/(n)_4 = 2·49/(11·10·9·8) = 98/7920 = 49/3960
s4_11 = Fraction(2*(11-4)**2, falling_factorial(11, 4))

# |S|=6 at n=11: g_3(m) = m(2m²+1)/3 with m = 11-6 = 5
# g_3(5) = 5(2·25+1)/3 = 5·51/3 = 85
g3_5 = 5 * (2*25 + 1) // 3
s6_11 = Fraction(2*g3_5, falling_factorial(11, 6))

# |S|=8 = cv2_11 - s2_11 - s4_11 - s6_11 - |S|=10
# |S|=10 should be tiny (only 1 tiling of 10 positions with 5 dominoes: {0,1,2,...,9})
# Actually: positions 0..9, 5 dominoes covering all 10 positions.
# C(10-5, 5) = C(5,5) = 1. So exactly 1 tiling.
# E[Z_0Z_1Z_2Z_3Z_4Z_5Z_6Z_7Z_8Z_9] = ?
# This is E[∏_{j=0}^{9} Z_j]. Since |S|=10 is even, it can be nonzero.
# But all 10 positions must have Z_j ≠ 0, meaning every consecutive pair
# in the permutation must be a unit ascent or descent.
# That means σ(1) = σ(0) ± 1, σ(2) = σ(1) ± 1, etc.
# The permutation must be a path on the integer graph!
# There are exactly 2 such paths: 0,1,2,...,10 and 10,9,...,0
# Each has Z_j = +1 or -1 respectively, product = 1.
# So E[∏Z_j] = 2/11! (both paths give product 1) = 2/39916800.

s10_11 = Fraction(2, nf11) * falling_factorial(11, 10) / 2
# Wait: |S|=10 = 2·g_5(1)/(11)_10. And g_5(1) should be 1 (based on pattern).
# (11)_10 = 11!/1! = 39916800
# |S|=10 = 2·1/39916800 = 1/19958400

# Actually let me just compute: the only |S|=10 subset is {0,1,...,9}.
# E[Z_0·Z_1·...·Z_9] = (# perms with all Z nonzero and product positive) - (product negative) / N
# All ascending: 0,1,...,10 → product = 1. All descending: 10,9,...,0 → Z_j all -1, product = (-1)^10 = 1.
# Any other configuration with all Z nonzero? It would need to be a path on integers
# that doesn't repeat values. E.g., 0,1,2,1,... repeats. So only the two monotone paths.
# So E = (1+1)/11! = 2/39916800.
# And |S|=10 = 2/39916800. Multiply by (11)_10/2 = 39916800/2 = 19958400.
# g_5(1) = (2/39916800) × 39916800/2 = 1. ✓

s10_11 = Fraction(2, nf11)

s8_11 = cv2_11 - s2_11 - s4_11 - s6_11 - s10_11
print(f"\nn=11:")
print(f"  CV² = {cv2_11} = {float(cv2_11):.12f}")
print(f"  |S|=2 = {s2_11} = {float(s2_11):.12f}")
print(f"  |S|=4 = {s4_11} = {float(s4_11):.12f}")
print(f"  |S|=6 = {s6_11} = {float(s6_11):.12f}")
print(f"  |S|=10 = {s10_11} = {float(s10_11):.15f}")
print(f"  |S|=8 = {s8_11} = {float(s8_11):.12f}")

# Extract g_4 at m = n-8 = 3
g4_val = s8_11 * falling_factorial(11, 8) / 2
print(f"  g_4(3) = {g4_val}")

# We know g_4(1) = 1, g_4(2) = 8
# Now g_4(3) = ?
print(f"\ng_4 values: g_4(1)=1, g_4(2)=8, g_4(3)={g4_val}")

# Fit degree 2 polynomial through (1,1), (2,8), (3,g4_val)
pts = [(Rational(1), Rational(1)), (Rational(2), Rational(8)),
       (Rational(3), Rational(g4_val.numerator, g4_val.denominator))]
poly4 = interpolate(pts, x)
print(f"  g_4(m) = {factor(poly4)} = {expand(poly4)}")

# Verify pattern: g_1 = m (deg 1), g_2 = m² (deg 2),
# g_3 = m(2m²+1)/3 (deg 3), g_4 should be deg 4?
# With only 3 points we can fit deg 2, not deg 4.
# Need more data.

# Try to get g_4(4) from n=12 using W(12)
# W(12) is not available from C code. Let me see if we can compute it.
# Actually, let me try computing W(12) directly in Python.
# n=12 has 479M perms — too slow for brute force.

# Instead, use the bitmask DP approach in Python for n=11 or 12
print("\n" + "="*70)
print("COMPUTING W(n) via bitmask DP for n=11,12")
print("="*70)

def compute_W_dp(n):
    """Bitmask DP for W(n) = Σ_{σ∈NUD(n)} 2^{adj1(σ)}."""
    full = (1 << n) - 1
    # dp[mask][v] = weighted count
    dp = {}

    # Init: single vertices
    for v in range(n):
        dp[((1 << v), v)] = 1

    for mask in range(1, full + 1):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            cnt = dp[key]

            for u in range(n):
                if mask & (1 << u):
                    continue
                if u == v - 1:  # unit descent, forbidden
                    continue
                weight = 2 * cnt if u == v + 1 else cnt
                new_key = (mask | (1 << u), u)
                dp[new_key] = dp.get(new_key, 0) + weight

    # Sum over full mask
    W = 0
    for v in range(n):
        key = (full, v)
        W += dp.get(key, 0)
    return W

for n_test in [11, 12]:
    print(f"\nComputing W({n_test}) via DP...", flush=True)
    import time
    t0 = time.time()
    W_val = compute_W_dp(n_test)
    t1 = time.time()
    nf = factorial(n_test)
    cv2 = Fraction(W_val, nf) - 1
    print(f"  W({n_test}) = {W_val} (took {t1-t0:.1f}s)")
    print(f"  CV² = {float(cv2):.12f}")
    print(f"  n×CV² = {float(n_test*cv2):.12f}")

    # Extract g_4
    s2 = Fraction(2*(n_test-2), n_test*(n_test-1))
    s4 = Fraction(2*(n_test-4)**2, falling_factorial(n_test, 4))
    m6 = n_test - 6
    g3_m = m6 * (2*m6**2 + 1) // 3 if m6 > 0 else 0
    s6 = Fraction(2*g3_m, falling_factorial(n_test, 6)) if m6 > 0 else 0

    remaining = cv2 - s2 - s4 - s6
    # This includes |S|=8 + |S|=10 + ...

    m8 = n_test - 8
    if m8 > 0:
        # |S|=10 is tiny, compute separately
        m10 = n_test - 10
        if m10 > 0:
            # |S|=10: the only full-tiling is 0..n-2 with (n-1)/2 dominoes if n-1 even
            # For n=12: positions 0..10, 11 positions, can't tile fully with dominoes (odd!)
            # So |S|=10 at n=12 needs 5 dominoes covering 10 of 11 positions. Not a full tiling.
            # Actually |S|=10 ≠ full tiling. |S|=10 means we pick 10 positions from n-1 = 11.
            # The 10 positions must be decomposable as 5 dominos.
            # This requires the 10 positions to be exactly 5 non-overlapping adjacent pairs.
            # At n=12 (m=11, positions 0..10): C(11-5,5) = C(6,5) = 6 tilings.
            # Each E[Z∏] involves 11 values...  complex to compute exactly.
            pass
        print(f"  |S|≥8 for n={n_test}: {float(remaining - s6 + s6):.15f}")
        print(f"    |S|=2={float(s2):.12f}, |S|=4={float(s4):.12f}, |S|=6={float(s6):.12f}")
        print(f"    Remaining (|S|≥8) = {float(remaining):.15f}")

        if m8 >= 1:
            # Approximate: assume |S|≥10 is negligible
            g4_approx = remaining * falling_factorial(n_test, 8) / 2
            print(f"    g_4({m8}) ≈ {float(g4_approx):.6f} (exact: {g4_approx})")

print("\nDone!")
