#!/usr/bin/env python3
"""
gk_degree_analysis_s112b.py — Follow-up: investigate the small j>=4 correction
coefficients in the shifted binomial basis C(m-1, j).

Key finding from s112: In the C(m-1, j) basis, the corrections E_k(m) have
large j=2,3 coefficients (related to the true g_k structure) but TINY j>=4
coefficients. Let's see if these tiny coefficients have structure.

kind-pasteur-2026-03-15-S112
"""

from fractions import Fraction
from math import comb, factorial
from functools import reduce
from sympy import symbols, Rational, Poly, expand, factor, cancel

x = symbols('x')

# True g_k polynomials (verified, from THM-216)
gk_coeffs_3 = {
    3: (2, 0, 1, 0),
    4: (10, -33, 50, -24),
    5: (388, -2040, 3431, -1776),
    6: (69660, -380445, 653748, -342960),
    7: (19826270, -109486152, 189674605, -100014720),
    8: (7309726742, -40641958545, 70757788486, -37425556680),
    9: (3262687720240, -18232387983408, 31858349908595, -16888649645424),
}

def g_true_poly(k):
    if k == 1: return x
    if k == 2: return x**2
    if k in gk_coeffs_3:
        a, b, c, d = gk_coeffs_3[k]
        return Rational(a, 3)*x**3 + Rational(b, 3)*x**2 + Rational(c, 3)*x + Rational(d, 3)
    return None

def g_naive_poly(k):
    result = Rational(0)
    for r in range(1, k + 1):
        coeff = comb(k - 1, r - 1) * 2**(r - 1)
        binom_poly = Rational(1)
        for j in range(r):
            binom_poly *= (x - j)
        binom_poly /= factorial(r)
        result += coeff * binom_poly
    return expand(result)

# =================================================================
# PART A: j>=4 coefficients in shifted basis for corrections
# =================================================================
print("=" * 78)
print("j>=4 CORRECTION COEFFICIENTS IN SHIFTED BASIS C(m-1, j)")
print("=" * 78)
print()
print("From Part 14 of s112, the corrections E_k in C(m-1,j) basis are:")
print("k=4: [-8]*C(m-1,4)")
print("k=5: [160]*C(m-1,2) + [720]*C(m-1,3) + [-48]*C(m-1,4) + [-16]*C(m-1,5)")
print("k=6: large*C(m-1,2) + large*C(m-1,3) + [-160]*C(m-1,4) + [-112]*C(m-1,5) + [-32]*C(m-1,6)")
print("etc.")
print()
print("The j>=4 coefficients are SMALL. Let's tabulate them:")
print()

all_shifted_high = {}  # k -> {j: coeff}

for k in range(4, 10):
    true_p = g_true_poly(k)
    naive_p = g_naive_poly(k)
    if true_p is None:
        break
    E = expand(true_p - naive_p)
    deg_E = Poly(E, x).degree()

    # Shifted binomial basis C(m-1, j): evaluate E at m=1+t for t=0,1,...
    shifted_vals = [E.subs(x, 1 + t) for t in range(deg_E + 2)]
    curr = list(shifted_vals)
    shifted_coeffs = []
    for j in range(len(curr)):
        shifted_coeffs.append(int(curr[0]))
        curr = [curr[i+1] - curr[i] for i in range(len(curr) - 1)]
        if not curr:
            break

    high_coeffs = {j: c for j, c in enumerate(shifted_coeffs) if j >= 4 and c != 0}
    all_shifted_high[k] = high_coeffs

    print(f"k={k}: j>=4 coeffs: {high_coeffs}")

# Compare with naive j>=4 coefficients (which are C(k-1,j-1)*2^{j-1})
print()
print("Naive j>=4 coefficients for comparison:")
for k in range(4, 10):
    naive_high = {j: comb(k-1, j-1) * 2**(j-1) for j in range(4, k+1)}
    print(f"k={k}: {naive_high}")

# =================================================================
# PART B: The j>=4 shifted correction coefficients are just
# -C(k-1,j-1)*2^{j-1} PLUS the contribution from lower terms?
# No -- they are different. Let me check against the standard
# binomial correction coefficients (which ARE -C(k-1,j-1)*2^{j-1}).
# =================================================================
print()
print("=" * 78)
print("COMPARISON: shifted vs standard binomial correction at j>=4")
print("=" * 78)

for k in range(4, 10):
    true_p = g_true_poly(k)
    naive_p = g_naive_poly(k)
    if true_p is None:
        break
    E = expand(true_p - naive_p)
    deg_E = Poly(E, x).degree()

    # Standard binomial basis: evaluate E at m=0,1,...
    std_vals = [E.subs(x, m) for m in range(deg_E + 2)]
    curr = list(std_vals)
    std_coeffs = []
    for j in range(len(curr)):
        std_coeffs.append(int(curr[0]))
        curr = [curr[i+1] - curr[i] for i in range(len(curr) - 1)]
        if not curr:
            break

    # Shifted: evaluate at m=1,2,...
    shifted_vals = [E.subs(x, 1 + t) for t in range(deg_E + 2)]
    curr = list(shifted_vals)
    shifted_coeffs = []
    for j in range(len(curr)):
        shifted_coeffs.append(int(curr[0]))
        curr = [curr[i+1] - curr[i] for i in range(len(curr) - 1)]
        if not curr:
            break

    print(f"\nk={k}:")
    print(f"  Standard C(m,j) [j>=4]: {[std_coeffs[j] if j < len(std_coeffs) else 0 for j in range(4, k+1)]}")
    print(f"  Shifted  C(m-1,j) [j>=4]: {[shifted_coeffs[j] if j < len(shifted_coeffs) else 0 for j in range(4, k+1)]}")
    print(f"  -Naive C(k-1,j-1)*2^(j-1) [j>=4]: {[-comb(k-1,j-1)*2**(j-1) for j in range(4, k+1)]}")

# =================================================================
# PART C: Factor the shifted j>=4 coefficients
# =================================================================
print()
print("=" * 78)
print("FACTORING THE SHIFTED j>=4 CORRECTION COEFFICIENTS")
print("=" * 78)
print()

# The shifted coefficients at j>=4 are small. Let me see if they
# relate to something like C(2k-2, 2j-2) or similar.

print(f"{'k':>3} {'j':>3} {'shifted coeff':>16} {'=-C(k-1,j-1)*2^(j-1)?':>25} {'factor':>20}")
print("-" * 70)

for k in range(4, 10):
    high = all_shifted_high.get(k, {})
    for j in sorted(high):
        c = high[j]
        naive_neg = -comb(k-1, j-1) * 2**(j-1)
        match = "YES" if c == naive_neg else f"NO (naive_neg={naive_neg})"

        # Try to factor as something nice
        # Check c / 2^(j-1)
        ratio = Fraction(c, 2**(j-1))
        # Check c / C(k-1, j-1)
        ratio2 = Fraction(c, comb(k-1, j-1)) if comb(k-1, j-1) != 0 else "N/A"

        print(f"{k:>3} {j:>3} {c:>16} {match:>25} c/2^(j-1)={ratio}, c/C={ratio2}")

# =================================================================
# PART D: Check if the shifted j>=4 coefficients follow a SIMPLE pattern
# =================================================================
print()
print("=" * 78)
print("SHIFTED j>=4 COEFFICIENTS: looking for 2^j * polynomial(k,j) pattern")
print("=" * 78)

# For each j=4,5,..., collect the coefficient across k values
for j in range(4, 10):
    print(f"\nj={j}:")
    coeffs_at_j = []
    for k in range(j, 10):
        high = all_shifted_high.get(k, {})
        c = high.get(j, 0)
        c_over_2j = Fraction(c, 2**j) if c != 0 else 0
        coeffs_at_j.append((k, c, c_over_2j))
        print(f"  k={k}: coeff = {c}, coeff/2^{j} = {c_over_2j}")

# =================================================================
# PART E: Look at the shifted coefficients divided by -2^{j-1}
# These might be falling factorials in k
# =================================================================
print()
print("=" * 78)
print("SHIFTED COEFFICIENTS / (-2^{j-1}): pattern in k?")
print("=" * 78)

for j in range(4, 10):
    print(f"\nj={j}: c_{j}^(k) / (-2^{j-1}) =")
    vals = []
    for k in range(j, 10):
        high = all_shifted_high.get(k, {})
        c = high.get(j, 0)
        ratio = Fraction(c, -(2**(j-1))) if c != 0 else 0
        vals.append((k, ratio))
        print(f"  k={k}: {ratio}")

    # Try to fit as polynomial in k
    if len(vals) >= 3:
        # Forward differences
        int_vals = [v[1] for v in vals]
        diffs = [int_vals]
        while len(diffs[-1]) > 1:
            new = [diffs[-1][i+1] - diffs[-1][i] for i in range(len(diffs[-1]) - 1)]
            diffs.append(new)

        print(f"  Forward differences: {[[float(x) if isinstance(x, Fraction) else x for x in d] for d in diffs]}")

# =================================================================
# PART F: The KEY structural observation: rewrite E_k in terms of
# the naive high-order piece and the "redistribution"
# =================================================================
print()
print("=" * 78)
print("PART F: REDISTRIBUTION POLYNOMIALS")
print("=" * 78)
print()
print("E_k(m) = [true g_k] - [naive g_k]")
print("       = [degree-3 polynomial] - [degree-k polynomial]")
print()
print("The high-order (j>=4) part of naive = sum C(k-1,j-1)*2^{j-1}*C(m,j)")
print("always factors as: x*(x-1)*(x-2)*(x-3) * Q_{k-4}(x) / denominator")
print("where Q is a polynomial of degree k-4.")
print()
print("This means: E_k has (x-1)(x-2) as a factor (confirmed: E_k(1)=E_k(2)=0).")
print("For k=4: E_4 = -(x-1)(x-2)(x-3)(x-4)/3")
print("  has EXTRA roots at x=3 and x=4!")
print()

for k in range(4, 10):
    true_p = g_true_poly(k)
    naive_p = g_naive_poly(k)
    if true_p is None:
        break
    E = expand(true_p - naive_p)

    # Find the factorization
    f = factor(E)
    print(f"k={k}: E_k factored = {f}")

    # Check E_k(3), E_k(4):
    e3 = E.subs(x, 3)
    e4 = E.subs(x, 4)
    print(f"       E_k(3) = {e3}")
    print(f"       E_k(4) = {e4}")
    print()

# =================================================================
# PART G: The pairwise rho is NOT constant! It depends on gap.
# But from the data: for gap>=1, rho stabilizes.
# This suggests the correlations involve only "nearest-neighbor" corrections.
# =================================================================
print()
print("=" * 78)
print("PART G: WHY alpha_{k,4} = 0 — the 4-body weight identity")
print("=" * 78)
print()
print("The effective weight for r clusters (all singletons) is:")
print("  w_r^eff(k) = alpha_{k,r} / C(k-1,r-1)")
print()
print("For r=1,2,3 the weights are nonzero and k-dependent.")
print("For r=4, the weight is ZERO for ALL k>=4.")
print()
print("This is equivalent to saying: the coefficient of C(m,4) in g_k(m)")
print("is zero for all k>=4.")
print()
print("To see WHY, note that the 4-cluster contribution to g_k involves:")
print("  Sum over (s1,...,s4) with s_i>=1, sum=k of:")
print("    C(m,4) * weight(s1,...,s4; n)")
print()
print("For the alpha_{k,4}=0 identity, we need:")
print("  Sum over compositions of k into 4 parts: weight(s1,...,s4; n) = 0")
print()
print("For k=4: the only composition is (1,1,1,1).")
print("  weight(1,1,1,1; n) = 0")
print("  This is the 4-body EXACT CANCELLATION.")
print()
print("For k=5: compositions are (2,1,1,1) (4 permutations) + (1,1,1,1) impossible.")
print("  Actually, compositions of 5 into 4 parts:")
print("  (2,1,1,1) and permutations: C(3,0)=1 way to pick which part is 2? No,")
print("  C(5-1,4-1) = C(4,3) = 4 total ordered compositions.")

# Let me enumerate compositions
from itertools import combinations_with_replacement

def compositions(n, k):
    """Generate all compositions of n into k positive parts."""
    if k == 1:
        yield (n,)
        return
    for first in range(1, n - k + 2):
        for rest in compositions(n - first, k - 1):
            yield (first,) + rest

for k in range(4, 8):
    print(f"\nk={k}: compositions of k into 4 parts:")
    comps = list(compositions(k, 4))
    print(f"  {len(comps)} compositions")
    for c in comps[:10]:
        print(f"    {c}")
    if len(comps) > 10:
        print(f"    ... ({len(comps) - 10} more)")

# =================================================================
# PART H: The r-cluster weight involves cluster sizes and CROSS-correlations
# =================================================================
print()
print("=" * 78)
print("PART H: WHAT THE ZERO EFFECTIVE WEIGHT MEANS PHYSICALLY")
print("=" * 78)
print("""
The effective r-cluster weight alpha_{k,r} / C(k-1,r-1) represents the
AVERAGE weight of an r-cluster configuration, summed over all compositions
of k into r parts, when placed with C(m,r) spacings.

For r >= 4 clusters, this average is ZERO. This means the weights
of different compositions (with their correlations) CANCEL each other.

For r=4, k=4: the only composition is (1,1,1,1).
  The weight is h_1^4 * [correlation correction] = 0.
  This means the 4-body correlation for 4 separated single-dominos
  is EXACTLY -h_1^4 (complete cancellation).

For r=4, k=5: compositions are (2,1,1,1) in 4 orderings.
  The total weight sum includes:
  - 4 terms with one cluster of size 2 and three singletons
  Their weights must sum to zero.

For r=4, k=6: compositions (3,1,1,1), (2,2,1,1), etc.
  Again, the weighted sum over all compositions must vanish.

The PATTERN: for r >= 4, the r-cluster weight polynomial
  W_r(k) = sum over compositions of k into r parts: weight(s1,...,sr; n=m+2k)
vanishes as a polynomial in m at degree r.

This is a NON-TRIVIAL identity about the multi-body correlations of the
domino process {Z_j Z_{j+1}} in random permutations.

PROOF APPROACH: This could potentially be proved by showing that the
moment-generating function of the domino process has a special form
(e.g., determinantal or Pfaffian structure) that forces the cumulants
beyond order 3 to contribute only at bounded degree.
""")

# =================================================================
# PART I: Check the pair-separated weight EXACTLY
# =================================================================
print("=" * 78)
print("PART I: VERIFYING 4-DOMINO WEIGHT IS ZERO")
print("=" * 78)
print()
print("For 4 separated single-dominos at positions p_1, p_1+1, p_2, p_2+1,")
print("p_3, p_3+1, p_4, p_4+1 (all gaps >= 1) in a permutation of n elements:")
print()
print("E[Z_{p1}Z_{p1+1} * Z_{p2}Z_{p2+1} * Z_{p3}Z_{p3+1} * Z_{p4}Z_{p4+1}]")
print()
print("The naive (product) approximation gives h_1^4 = 16/(n)_2^4.")
print()
print("If the pairwise correlations are the ONLY correlations (no 3-body+),")
print("then the weight would be h_1^4 * rho^6 = 16/(n)_2^4 * [C(n,2)/C(n-2,2)]^6.")
print()
print("The effective weight in g_k is:")
print("  alpha_{k,4}/C(k-1,3) * C(m,4) / [(n)_{2k}/2 * h_1^4 * C(m,4)]")
print("  = 0.")
print()
print("This means EITHER:")
print("(a) The pairwise factorization is wrong (there are genuine 3-body and 4-body cumulants), OR")
print("(b) The pairwise rho^6 correction already gives zero.")
print()
print("Let's check (b): does rho^6 * h_1^4 * (n)_8/2 * C(m,4) give the right g_4 coefficient?")

# For k=4, r=4: contribution = C(m,4) * (n)_8/2 * h_1^4 * rho^6
# where rho = n(n-1)/((n-2)(n-3))
# h_1 = 2/(n(n-1))
# h_1^4 = 16/(n(n-1))^4
# rho^6 = [n(n-1)]^6 / [(n-2)(n-3)]^6
# h_1^4 * rho^6 = 16 * [n(n-1)]^2 / [(n-2)(n-3)]^6
# (n)_8 = n(n-1)(n-2)(n-3)(n-4)(n-5)(n-6)(n-7)
# Contribution = C(m,4) * (n)_8 / 2 * 16 * [n(n-1)]^2 / [(n-2)(n-3)]^6
#              = C(m,4) * 8 * (n)_8 * [n(n-1)]^2 / [(n-2)(n-3)]^6

# For n = m + 8 (since k=4):
# Let's compute this symbolically

n_sym = symbols('n')
m_sym = symbols('m')

# (n)_8 = n!/(n-8)!
ff_8 = 1
for i in range(8):
    ff_8 *= (n_sym - i)
ff_8 = expand(ff_8)

# h_1^4 * rho^6 * (n)_8 / 2
# = 8 * ff_8 * [n(n-1)]^2 / [(n-2)(n-3)]^6

numerator = 8 * ff_8 * (n_sym * (n_sym - 1))**2
denominator = ((n_sym - 2) * (n_sym - 3))**6

ratio = cancel(numerator / denominator)
print(f"\nPairwise rho^6 contribution (before C(m,4)):")
print(f"  = {ratio}")
print(f"  = {factor(ratio)}")

# Substitute n = m + 8:
ratio_m = ratio.subs(n_sym, m_sym + 8)
print(f"\nSubstituting n = m + 8:")
print(f"  = {factor(ratio_m)}")
print(f"  This is a rational function of m, NOT a polynomial.")
print(f"  And NOT zero. So pairwise rho alone does NOT explain alpha_4 = 0.")

# The naive contribution (without rho) is C(m,4) * 8 = 8 * C(m,4)
# (from the naive formula: C(3,3)*2^3 = 8)
# So the pairwise-rho-corrected contribution is C(m,4) * ratio(n=m+8)
# And alpha_{4,4} = ratio(n=m+8) SHOULD be... a constant?
# No, ratio is a rational function of n, not a polynomial.

# Actually the issue is that rho^6 as a correction is for the ALL-SEPARATED
# case, but the true factorization might not be pairwise.
print()
print("The pairwise rho model gives a non-constant rational function of n,")
print("not zero. This means the pairwise factorization is INSUFFICIENT.")
print("There MUST be genuine higher-body cumulants (3-body, 4-body, etc.)")
print("that contribute to making alpha_{k,4} = 0.")
print()
print("CONCLUSION: The degree-3 property requires a FULL cumulant analysis,")
print("not just pairwise correlations. The vanishing of alpha_{k,r} for r>=4")
print("is a deep identity involving all orders of cumulants simultaneously.")

# =================================================================
# PART J: The correction coefficients at j=0,1,2 have a striking pattern
# =================================================================
print()
print("=" * 78)
print("PART J: PATTERN IN LOW-ORDER CORRECTION COEFFICIENTS")
print("=" * 78)
print()
print("From Part 6 of s112, the low-order corrections are:")
print()

for k in range(4, 10):
    true_p = g_true_poly(k)
    naive_p = g_naive_poly(k)
    if true_p is None:
        break
    E = expand(true_p - naive_p)
    deg_E = Poly(E, x).degree()

    # Standard binomial basis
    std_vals = [E.subs(x, m) for m in range(deg_E + 2)]
    curr = list(std_vals)
    std_coeffs = []
    for j in range(len(curr)):
        std_coeffs.append(int(curr[0]))
        curr = [curr[i+1] - curr[i] for i in range(len(curr) - 1)]
        if not curr:
            break

    c0, c1, c2, c3 = std_coeffs[0], std_coeffs[1], std_coeffs[2], std_coeffs[3]

    print(f"k={k}: c0={c0}, c1={c1}, c2={c2}, c3={c3}")
    print(f"       c0 = -c1 = c2? {c0 == -c1 == c2}")
    print(f"       |c0| = {abs(c0)}")

    # Factor |c0|
    from sympy import factorint
    if abs(c0) > 1:
        print(f"       |c0| = {dict(factorint(abs(c0)))}")

    # Ratio c3 / c0
    if c0 != 0:
        ratio_30 = Fraction(c3, c0)
        print(f"       c3/c0 = {ratio_30} = {float(ratio_30):.6f}")
    print()

print()
print("STRIKING PATTERN: c0 = -c1 = c2 for ALL k >= 4!")
print("And c3 differs, absorbing the remaining mass.")
print()
print("c0 values: ", end="")
c0_vals = []
for k in range(4, 10):
    true_p = g_true_poly(k)
    naive_p = g_naive_poly(k)
    if true_p is None:
        break
    E = expand(true_p - naive_p)
    c0 = int(E.subs(x, 0))
    c0_vals.append(c0)
    print(f"{c0}", end=", ")
print()

print()
print("These are: -g_k(0) for k>=4 (since naive(0) = 0).")
print("g_k(0) values:")
for k in range(3, 10):
    p = g_true_poly(k)
    if p is None:
        break
    v = p.subs(x, 0)
    print(f"  g_{k}(0) = {v}")

print()
print("g_k(0) = d_k/3 where d_k is from the coefficient table.")
print("g_k(0) is the 'boundary value' of the polynomial at m=0.")
print("It's NOT a meaningful value (m>=1 in practice), but it determines")
print("the constant term in the binomial expansion.")

# =================================================================
# PART K: Summary of key identities
# =================================================================
print()
print("=" * 78)
print("KEY IDENTITIES DISCOVERED")
print("=" * 78)
print("""
1. NAIVE FORMULA EXACTNESS: g_k^naive(m) = sum C(k-1,r-1)*C(m,r)*2^{r-1}
   is EXACT for k = 1, 2, 3. This is the cluster-independence formula.

2. CORRECTION STRUCTURE: For k >= 4, E_k(m) = g_k(m) - g_k^naive(m):
   (a) E_k always has (m-1)(m-2) as a factor (E_k(1) = E_k(2) = 0)
   (b) For k=4: E_4 = -(m-1)(m-2)(m-3)(m-4)/3 (complete falling factorial)
   (c) For k>=5: E_k = (m-1)(m-2) * Q_{k-2}(m) where Q has degree k-2

3. BINOMIAL COEFFICIENTS:
   (a) For j >= 4: correction coeff = -C(k-1,j-1)*2^{j-1} (trivially, since true=0)
   (b) For j = 0,1,2: correction coeff c_0 = -c_1 = c_2 = -g_k(0) for all k >= 4
   (c) For j = 3: correction coeff absorbs remaining mass

4. EFFECTIVE CLUSTER WEIGHTS:
   - For r=1: weight = alpha_{k,1} (grows explosively with k)
   - For r=2: weight = alpha_{k,2}/C(k-1,1) (grows, negative)
   - For r=3: weight = alpha_{k,3}/C(k-1,2) (grows)
   - For r>=4: weight = 0 for ALL k (the central mystery)

5. INTERPRETATION: The vanishing of r>=4 effective weights means that
   in the cumulant expansion of the r-domino joint moment, the sum of
   all partition contributions at "4+ clusters" level cancels exactly.
   This is a COMBINATORIAL IDENTITY about random permutation statistics.

6. OPEN: Prove alpha_{k,r} = 0 for r >= 4 analytically. This requires
   understanding the cumulant structure of the domino process
   {Z_j * Z_{j+1}} over random permutations.
""")

print("Done!")
