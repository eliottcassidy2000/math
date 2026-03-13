#!/usr/bin/env python3
"""
Deep investigation of H(P_p) formula.

DISCOVERED: H(P_3) = 3·1^1, H(P_7) = 7·3^3 = p·m^m (exact!)
But H(P_11) = 95095 ≠ 11·5^5.

H(P_11) = 5·7·11·13·19 = 11 · 5·7·13·19 = 11 · 8645
H(P_19) = 3²·5·7·11·19·23·774463

Note 5·7·13·19 = 8645. And 8645/5^5 = 2.7664.

Let's look at this differently:
  H(P_3) = 3 = (2)! · 3/2 = (p-1)! · p/(p-1)!... hmm
  H(P_7) = 189 = 6! · 189/720... nah

Let me look at the combinatorial structure more carefully.

Author: opus-2026-03-13-S67k
"""

from math import factorial, comb, prod
import numpy as np

def paley_adjacency(p):
    qr = set()
    for x in range(1, p):
        qr.add((x * x) % p)
    A = np.zeros((p, p), dtype=int)
    for i in range(p):
        for j in range(p):
            if i != j and ((j - i) % p) in qr:
                A[i][j] = 1
    return A

# Known values
known = {3: 3, 7: 189, 11: 95095, 19: 1172695746915}

print("=" * 70)
print("DEEP ANALYSIS OF H(P_p) FORMULA")
print("=" * 70)

# Check double factorial and related
print("\n--- Various formula candidates ---\n")
for p, H in known.items():
    m = (p - 1) // 2
    print(f"p={p}, m={m}, H={H}:")

    # (2m)!! = 2^m · m!
    double_fact = prod(range(2, 2*m+1, 2))
    print(f"  (2m)!! = {double_fact}")
    print(f"  H/(2m)!! = {H/double_fact:.6f}")

    # p · C(2m, m)
    candidate = p * comb(2*m, m)
    print(f"  p·C(2m,m) = {candidate}, ratio = {H/candidate:.6f}")

    # p · (2m)! / (m!)²
    candidate2 = p * factorial(2*m) // (factorial(m)**2)
    print(f"  p·(2m)!/(m!)² = {candidate2}, ratio = {H/candidate2:.6f}")

    # Check if H divides (p-1)!·p = p!
    print(f"  p!/H = {factorial(p)/H:.6f}")
    print(f"  (p-1)!/H = {factorial(p-1)/H:.6f}")
    print(f"  2^m·(p-1)!/(p·H) = {2**m * factorial(p-1)/(p*H):.6f}")

    # H·2^m / p!
    print(f"  H·2^m / p! = {H * 2**m / factorial(p):.6f}")

    # Does H = (product of first m odd primes)?
    print(f"  (2m+1)!/2^m/m! = {factorial(2*m+1)/(2**m * factorial(m)):.2f}, H/this = {H/(factorial(2*m+1)/(2**m * factorial(m))):.6f}")

    # Number of different regular tournaments
    # For P_p, the cycle structure is key
    print()

# The H·2^m/p! ratio is interesting:
# p=3: 1.0
# p=7: 0.3
# p=11: 0.076235
# p=19: ...

print("--- H · 2^{p-1} / p! = H / E[H] pattern ---\n")
for p, H in known.items():
    ratio = H * 2**(p-1) / factorial(p)
    log_ratio = np.log(ratio)
    print(f"p={p}: H/E[H] = {ratio:.6f}, ln = {log_ratio:.6f}, /ln(p) = {ratio/np.log(p):.6f}")

# Look for a PRODUCT formula
# H(P_3) = 3
# H(P_7) = 3^3 · 7 = 27 · 7
# H(P_11) = 5 · 7 · 11 · 13 · 19
# H(P_19) = 9 · 5 · 7 · 11 · 19 · 23 · 774463

# The factorization of H(P_11) = 5·7·11·13·19
# These are the primes 5, 7, 11, 13, 19.
# The primes between m+1=6 and 2m+1=11 are: 7, 11
# Actually let me look at this differently.
# 5·7·11·13·19 = ?
# 5·19 = 95, 7·13 = 91, 11: product = 95·91·11 = 95095. Checks out.

print("\n--- Checking if H relates to C(p-1, m) ---\n")
for p, H in known.items():
    m = (p - 1) // 2
    # C(p-1, m) = C(2m, m)
    central = comb(2*m, m)
    print(f"p={p}: C({2*m},{m}) = {central}, H/C = {H/central:.6f}")
    print(f"  H / (p · C(2m,m)) = {H/(p*central):.6f}")

# For P_7: H/C(6,3) = 189/20 = 9.45
# For P_11: H/C(10,5) = 95095/252 = 377.36

# Let me try: H(P_p) = product_{k=1}^{m} (something about k)
print("\n--- Checking product formulas ---\n")
for p, H in known.items():
    m = (p - 1) // 2
    # Product of (p - 2k) for k = 0..m-1? That gives p, p-2, p-4, ..., 1
    prod_odd = prod(range(p, 0, -2))
    print(f"p={p}: p!! = {prod_odd}, H/p!! = {H/prod_odd:.6f}")

    # Product of (2k+1) for k = 0..m
    prod_odd2 = prod(2*k+1 for k in range(m+1))
    print(f"  (2m+1)!! = {prod_odd2}, H/(2m+1)!! = {H/prod_odd2:.6f}")

# VERY INTERESTING:
# p=3: H/p!! = 3/3 = 1
# p=7: H/p!! = 189/105 = 1.8
# p=11: H/p!! = 95095/10395 = 9.148...
# p=19: H/p!! = ...

print("\n--- RATIO: H / prod(2k+1, k=0..m) = H/(2m+1)!! ---\n")
for p, H in known.items():
    m = (p - 1) // 2
    odd_prod = prod(2*k+1 for k in range(m+1))
    ratio = H / odd_prod
    print(f"p={p}, m={m}: (2m+1)!! = {odd_prod}, H/(2m+1)!! = {ratio:.6f}")

# Let me check if H(P_p) = permanent of some matrix related to Paley
print("\n--- Permanent analysis ---\n")
for p in [3, 7, 11]:
    A = paley_adjacency(p)
    # Compute permanent by inclusion-exclusion (Ryser's formula)
    # For small n this is feasible
    n = p
    perm = 0
    for S_mask in range(1 << n):
        # S = columns included
        col_sums = np.zeros(n, dtype=int)
        sign = (-1) ** (n - bin(S_mask).count('1'))
        for j in range(n):
            if S_mask & (1 << j):
                col_sums += A[:, j]
        perm += sign * np.prod(col_sums)

    H = known[p]
    print(f"P_{p}: perm(A) = {perm}, H = {H}")
    print(f"  H/perm = {H/perm:.6f}" if perm != 0 else "  perm = 0")

    # Also check perm of (I + A)
    IpA = (np.eye(n) + A).astype(int)
    perm_IpA = 0
    for S_mask in range(1 << n):
        col_sums = np.zeros(n, dtype=int)
        sign = (-1) ** (n - bin(S_mask).count('1'))
        for j in range(n):
            if S_mask & (1 << j):
                col_sums += IpA[:, j]
        perm_IpA += sign * np.prod(col_sums)

    print(f"  perm(I+A) = {perm_IpA}")
    print(f"  H/perm(I+A) = {H/perm_IpA:.6f}" if perm_IpA != 0 else "  perm(I+A) = 0")

    # Check det(I+A)
    det_IpA = int(round(np.linalg.det(IpA.astype(float))))
    print(f"  det(I+A) = {det_IpA}")
    print()

# KEY OBSERVATION: H(P_p) values
# p=3: 3
# p=7: 189 = 3³ · 7
# p=11: 95095 = 5 · 7 · 11 · 13 · 19
# p=19: 1172695746915 = 3² · 5 · 7 · 11 · 19 · 23 · 774463

# Let me check: is H(P_11) = 11!! · something?
# 11!! = 10395. 95095/10395 = 9.148...
# Is 95095/11 = 8645 = 5·7·13·19?
# 5·7 = 35, 13·19 = 247, 35·247 = 8645. Yes.
# These are (m, m+2, p+2, 2p-1) = (5, 7, 13, 19)?
# m=5, p=11. 5,7,13,19 — primes, but what's the pattern?

# For p=7: H/p = 189/7 = 27 = 3^3 = m^m
# For p=11: H/p = 95095/11 = 8645. Is this m^m = 5^5 = 3125? No.
#   8645 / 3125 = 2.7664. Hmm.

# What if we look at H as multinomial-like?
# H(P_7) = 7!/(2!·2!·2!) = 5040/8 = 630? No, H=189.
# 189 = C(7,3) · 27/5 ... messy.

# OCF formula: H = 1 + 2α₁ + 4α₂ + ...
# For Paley P_p: α₁ = c3 + c5 + c7 = all odd cycles
# α₂ = number of disjoint cycle pairs
# For P_7: α₁ = 14 + 42 + 24 = 80, α₂ = 7
# H = 1 + 2·80 + 4·7 = 1 + 160 + 28 = 189 ✓

print("--- OCF verification ---\n")
for p in [7, 11]:
    A = paley_adjacency(p)
    n = p
    # Count 3-cycles
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] and A[j][k] and A[k][i]) or \
                   (A[i][k] and A[k][j] and A[j][i]):
                    c3 += 1

    # Count 5-cycles (for P_7 only — too expensive for P_11 naively)
    if p <= 7:
        from itertools import permutations
        c5 = 0
        for perm in combinations(range(n), 5):
            for cycle_order in permutations(perm):
                is_cycle = True
                for idx in range(5):
                    if not A[cycle_order[idx]][cycle_order[(idx+1) % 5]]:
                        is_cycle = False
                        break
                if is_cycle:
                    c5 += 1
        c5 //= 5  # each cycle counted 5 times (cyclic symmetry)
        # But also counted in both directions... for directed odd cycles,
        # each cycle has exactly one direction
        # Actually for tournaments, an odd cycle on 5 vertices can only go one way
        # (determined by arc orientations). So c5 counted correctly.

        # Actually c5 should be counted more carefully.
        # A directed 5-cycle on vertices v_0,...,v_4 means
        # v_0→v_1→v_2→v_3→v_4→v_0. Count = (1/5) · Σ over 5-element
        # subsets of 5-permutations that form a cycle.
        # The above counts each directed cycle 5 times (rotations), so /5 is correct.

        # For P_7: expected c5 = 42 (from cross-field synthesis)
        print(f"P_{p}: c3 = {c3}, c5 = {c5}")

        # Alpha_1: total odd cycles
        # For P_7: also need c7
        if p == 7:
            c7 = 0
            # Only one 7-cycle vertex set: all vertices
            for perm in permutations(range(7)):
                is_cycle = True
                for idx in range(7):
                    if not A[perm[idx]][perm[(idx+1) % 7]]:
                        is_cycle = False
                        break
                if is_cycle:
                    c7 += 1
            c7 //= 7  # cyclic rotations
            print(f"  c7 = {c7}")
            alpha1 = c3 + c5 + c7
            print(f"  alpha1 = {alpha1}")

            # Alpha_2: disjoint cycle pairs
            # At n=7: only 3+3 disjoint pairs possible (need 6 vertices, leave 1 out)
            # For each vertex v not in a 3-cycle pair:
            #   count pairs of vertex-disjoint 3-cycles in T-v
            # This is disj_33 from the cycle cascade analysis
            # For Paley P_7: disj_33 = 7 (from cycle_cascade_n7.py output)
            alpha2 = 7  # from known result
            print(f"  alpha2 = {alpha2} (from known)")
            H_ocf = 1 + 2*alpha1 + 4*alpha2
            print(f"  H = 1 + 2·{alpha1} + 4·{alpha2} = {H_ocf}")
            print(f"  Known H = {known[p]}")
            print(f"  Match: {H_ocf == known[p]}")
    else:
        print(f"P_{p}: c3 = {c3} (c5, c7 too expensive to compute here)")

    print()

# The real insight: H(P_p)/E[H] ~ ln(p) · constant
# p=3: 2.0, ln(3) = 1.099
# p=7: 2.4, ln(7) = 1.946
# p=11: 2.44, ln(11) = 2.398
# p=19: 2.527, ln(19) = 2.944

# Ratio / ln(p): 1.82, 1.23, 1.02, 0.86
# This is DECREASING. So it's not ~ ln(p).

# Maybe ~ (ln p)^{1/2}?
# sqrt(ln 3) = 1.048, 2/1.048 = 1.91
# sqrt(ln 7) = 1.395, 2.4/1.395 = 1.72
# sqrt(ln 11) = 1.548, 2.44/1.548 = 1.58
# sqrt(ln 19) = 1.716, 2.527/1.716 = 1.47
# Still decreasing but more slowly.

# Try ratio / (ln p / ln ln p):
print("--- Asymptotic analysis ---\n")
for p, H in known.items():
    ratio = H * 2**(p-1) / factorial(p)
    lp = np.log(p)
    llp = np.log(np.log(p)) if lp > 1 else 1
    print(f"p={p}: ratio={ratio:.4f}, ln(p)={lp:.4f}, ratio/sqrt(ln p)={ratio/np.sqrt(lp):.4f}, ratio·ln(ln p)/ln p = {ratio*llp/lp:.4f}")

print("\nDone.")
