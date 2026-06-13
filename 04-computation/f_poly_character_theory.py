#!/usr/bin/env python3
"""
f_poly_character_theory.py — F(T,x) as a character polynomial.

KEY INSIGHT: F(T,x) = sum_{sigma} x^{des_T(sigma)} where des_T(sigma)
counts the "T-descents" of the permutation sigma.

For the transitive tournament, des_T = standard descents.
The Eulerian polynomial A_n(x) = sum x^{des(sigma)} is the character
of the "descent representation" of S_n.

For general tournaments, des_T(sigma) = # positions i where
A[sigma(i)][sigma(i+1)] = 1 (forward edge).

QUESTION: Does F(T,x) encode the character of some S_n representation?

SOLOMON'S DESCENT ALGEBRA: For any subset S of {1,...,n-1},
  D_S = sum_{des(sigma)=S} sigma
is an element of the group algebra. The D_S span Solomon's descent algebra.

For tournaments, we have a SINGLE STATISTIC (number of T-descents),
not a descent SET. So F(T,x) is a "type B" version — it only counts
the size of the descent set, not which positions are descents.

NOVEL: The Ehrhart-Macdonald reciprocity for F(T,x)?
F(T,x) = x^{n-1} F(T, 1/x) (palindrome). This is exactly Ehrhart
reciprocity for a reflexive polytope! So F(T,x)/H might be the
h*-polynomial of some polytope.

Author: opus-2026-03-07-S44
"""
from itertools import permutations, combinations
import math
import random
from functools import reduce
from math import gcd

def tournament_from_bits(bits, n):
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> pos) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A

def compute_F_dp(A, n):
    full = (1 << n) - 1
    dp = [[None]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = {0: 1}
    for mask in range(1, 1 << n):
        for v in range(n):
            if dp[mask][v] is None:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                new_mask = mask | (1 << u)
                fwd_edge = A[v][u]
                if dp[new_mask][u] is None:
                    dp[new_mask][u] = {}
                for fwd_count, cnt in dp[mask][v].items():
                    new_fwd = fwd_count + fwd_edge
                    dp[new_mask][u][new_fwd] = dp[new_mask][u].get(new_fwd, 0) + cnt
    F = [0] * n
    for v in range(n):
        if dp[full][v] is not None:
            for fwd_count, cnt in dp[full][v].items():
                F[fwd_count] += cnt
    return F

# ============================================================
# DESCENT SET REFINEMENT
# ============================================================
print("=" * 60)
print("DESCENT SET DISTRIBUTION (not just count)")
print("=" * 60)
# For each permutation, record WHICH positions are T-forward.
# The descent set is a subset of {0, 1, ..., n-2}.
# F_k counts perms with |descent_set| = k.
# But the FULL descent set gives more info.

for n in [5]:
    m = n*(n-1)//2
    print(f"\nn={n}:")

    seen_F = set()
    for bits in range(1 << m):
        A = tournament_from_bits(bits, n)
        F = compute_F_dp(A, n)
        key = tuple(F)
        if key in seen_F:
            continue
        seen_F.add(key)
        H = F[n-1]

        # Count by descent set
        desc_sets = {}
        for P in permutations(range(n)):
            desc = frozenset(i for i in range(n-1) if A[P[i]][P[i+1]])
            desc_sets[desc] = desc_sets.get(desc, 0) + 1

        # How many distinct descent set counts?
        distinct_sets = len(desc_sets)
        # For comparison: there are 2^{n-1} possible descent sets
        total_possible = 2**(n-1)

        # beta coefficients: beta_S = # perms with descent set exactly S
        # These are the "refined" coefficients.
        # Check: sum of beta_S for |S|=k should = F_k
        for k in range(n):
            computed = sum(cnt for S, cnt in desc_sets.items() if len(S) == k)
            assert computed == F[k], f"Mismatch at k={k}"

        # Check symmetry: does beta_{S} = beta_{complement(S)}?
        # S = {0,...,n-2} \ S' maps to complement via palindrome
        complement = lambda S: frozenset(range(n-1)) - S
        symm_ok = all(
            desc_sets.get(S, 0) == desc_sets.get(complement(S), 0)
            for S in desc_sets
        )

        # Which sets are achieved?
        empty_sets = total_possible - distinct_sets

        print(f"  H={H:3d}: {distinct_sets}/{total_possible} desc sets used, "
              f"complement symmetry: {symm_ok}")

        # Show the beta_S for small n
        if n == 5:
            for S in sorted(desc_sets.keys(), key=lambda s: (len(s), sorted(s))):
                comp = complement(S)
                print(f"    S={str(set(S)):15s} beta={desc_sets[S]:3d}, comp={str(set(comp)):15s} beta_comp={desc_sets.get(comp,0):3d}")

# ============================================================
# NOVEL: h*-POLYNOMIAL INTERPRETATION
# ============================================================
print("\n" + "=" * 60)
print("F(T,x)/H AS h*-POLYNOMIAL OF REFLEXIVE POLYTOPE?")
print("=" * 60)
# If F(T,x)/H is the h*-polynomial of a reflexive polytope P,
# then P has volume H/n! * something.
# The h*-poly satisfies: h*(x) = x^d h*(1/x) (palindromic)
# and h*(1) = d! * Vol(P) / some normalization.
# h*(1) = F(1)/H = n!/H.

for n in [5]:
    m = n*(n-1)//2
    print(f"\nn={n}:")

    seen = set()
    for bits in range(1 << m):
        A = tournament_from_bits(bits, n)
        F = compute_F_dp(A, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)
        H = F[n-1]

        h_star = [F[k]/H for k in range(n)]
        h_star_1 = sum(h_star)  # = n!/H

        # Check: all h*_k > 0? (necessary for Ehrhart theory)
        all_pos = all(h >= 0 for h in h_star)

        # Check: h*_k integer? (necessary for lattice polytope)
        all_int = all(abs(h - round(h)) < 1e-10 for h in h_star)

        print(f"  H={H:3d}: h*={[f'{h:.2f}' for h in h_star]}, "
              f"h*(1)={h_star_1:.2f}=n!/H, int={all_int}, pos={all_pos}")

# ============================================================
# NOVEL: F(T,x) AND THE EULERIAN POLYNOMIAL RATIO
# ============================================================
print("\n" + "=" * 60)
print("F(T,x) / Eul_n(x) RATIO AT SPECIFIC POINTS")
print("=" * 60)
# The Eulerian polynomial E_n(x) = sum A(n,k) x^k
# F(transitive, x) = E_n(x) (up to reversal)
# For other tournaments, how does F relate to E?

def eulerian_poly(n):
    """Compute Eulerian polynomial coefficients A(n,0), A(n,1), ..., A(n,n-1)."""
    A = [0] * n
    A[0] = 1
    for i in range(1, n):
        new_A = [0] * n
        for k in range(n):
            new_A[k] = (k+1) * A[k] + (i - k) * (A[k-1] if k > 0 else 0)
        A = new_A
    return A

for n in [5, 7]:
    eul = eulerian_poly(n)
    eul_rev = list(reversed(eul))
    print(f"\n  n={n}: Eulerian (reversed) = {eul_rev}")

    # For each x, compute ratio F(T,x)/E_n(x)
    m = n*(n-1)//2
    if n == 5:
        iterator = range(1 << m)
    else:
        random.seed(42)
        iterator = [random.getrandbits(m) for _ in range(100)]

    seen = set()
    for bits in iterator:
        A = tournament_from_bits(bits, n)
        F = compute_F_dp(A, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)

        H = F[n-1]
        # Ratio at x=2
        F_2 = sum(F[k] * 2**k for k in range(n))
        E_2 = sum(eul_rev[k] * 2**k for k in range(n))
        ratio_2 = F_2 / E_2

        # Ratio at x=1 (should be 1 since both = n!)
        # Ratio at x=-1
        F_neg1 = sum(F[k] * (-1)**k for k in range(n))
        E_neg1 = sum(eul_rev[k] * (-1)**k for k in range(n))

        if n <= 5:
            print(f"  H={H:3d}: F(2)/E(2)={ratio_2:.6f}, "
                  f"F(-1)={F_neg1:5d}, E(-1)={E_neg1:5d}")
