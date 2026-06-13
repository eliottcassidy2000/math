#!/usr/bin/env python3
"""
H_moment_formula.py -- kind-pasteur-2026-03-13-S60

At p=11, H is an EXACT linear function of (S4, S6, S8):
  H = a*S4 + b*S6 + c*S8 + d

The coefficient of S10 is essentially 0 (< 6 in the fit).
This means the S10 contribution cancels in the H formula.

But H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ... and
alpha_j depends on the full cycle interaction structure.

This script:
1. Determines the EXACT integer coefficients in H = f(S4,S6,S8) at p=11
2. Explains why S10 cancels via the OCF structure
3. Tests whether this extends to p=13 (where H computation is infeasible
   but we can check whether the moment structure predicts correctly)
4. Derives the theoretical formula from Re(z^k) expansions
"""

import cmath
import numpy as np
from itertools import combinations
from collections import defaultdict
from fractions import Fraction


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i+s)%p] = 1
    return A


def compute_H(A, p):
    n = p
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n + 1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(n):
                    if (prev_mask & (1 << u)) and A[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def compute_eigenvalue_moments(p, S, max_moment=12):
    m = (p - 1) // 2
    omega = cmath.exp(2j * cmath.pi / p)
    D_vals = []
    for t in range(1, m + 1):
        lam = sum(omega ** (s * t) for s in S)
        D_vals.append(lam.imag)
    moments = {}
    for k in range(2, max_moment + 1, 2):
        moments[k] = sum(d**k for d in D_vals)
    return moments


def exact_coefficients(p):
    """Find EXACT rational coefficients of H = a*S4 + b*S6 + c*S8 + d."""
    m = (p - 1) // 2
    N = 1 << m

    if p > 11:
        print(f"  p={p} too large for exact H computation")
        return

    print(f"\n{'='*70}")
    print(f"EXACT H FORMULA at p={p}")
    print(f"{'='*70}")

    # Collect data
    data = []
    for bits in range(N):
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))
        A = build_adj(p, S)
        H = compute_H(A, p)
        moments = compute_eigenvalue_moments(p, S)
        data.append((H, moments))

    # Extract unique (S4, S6, S8, ...) tuples
    seen = {}
    for H, moments in data:
        key = (round(moments[4], 4),
               round(moments.get(6, 0), 4),
               round(moments.get(8, 0), 4))
        if key not in seen:
            seen[key] = H
        else:
            assert seen[key] == H, f"Non-unique: {key} maps to {seen[key]} and {H}"

    print(f"  {len(seen)} unique (S4,S6,S8) tuples, all mapping to unique H")

    # Use linear algebra with Fraction for exact coefficients
    # H = a*S4 + b*S6 + c*S8 + d
    # Need at least 4 equations
    keys = sorted(seen.keys())
    n_eq = len(keys)

    if n_eq >= 4:
        # Build system
        A_mat = []
        b_vec = []
        for s4, s6, s8 in keys:
            A_mat.append([Fraction(s4).limit_denominator(100000),
                          Fraction(s6).limit_denominator(100000),
                          Fraction(s8).limit_denominator(100000),
                          Fraction(1)])
            b_vec.append(Fraction(seen[(s4, s6, s8)]))

        # Gaussian elimination on Fraction matrix
        n_var = 4
        aug = [A_mat[i] + [b_vec[i]] for i in range(n_eq)]

        # Forward elimination
        for col in range(n_var):
            # Find pivot
            pivot_row = None
            for row in range(col, n_eq):
                if aug[row][col] != 0:
                    pivot_row = row
                    break
            if pivot_row is None:
                print(f"  Singular system at column {col}")
                return
            aug[col], aug[pivot_row] = aug[pivot_row], aug[col]

            # Eliminate below
            for row in range(col + 1, n_eq):
                if aug[row][col] != 0:
                    factor = aug[row][col] / aug[col][col]
                    for j in range(n_var + 1):
                        aug[row][j] -= factor * aug[col][j]

        # Back substitution
        coeffs = [Fraction(0)] * n_var
        for col in range(n_var - 1, -1, -1):
            s = aug[col][n_var]
            for j in range(col + 1, n_var):
                s -= aug[col][j] * coeffs[j]
            coeffs[col] = s / aug[col][col]

        print(f"\n  EXACT coefficients:")
        print(f"    H = ({coeffs[0]})*S4 + ({coeffs[1]})*S6 + ({coeffs[2]})*S8 + ({coeffs[3]})")
        print(f"    = ({float(coeffs[0]):.6f})*S4 + ({float(coeffs[1]):.6f})*S6 + ({float(coeffs[2]):.6f})*S8 + ({float(coeffs[3]):.2f})")

        # Verify on all data
        max_err = 0
        for s4, s6, s8 in keys:
            s4f = Fraction(s4).limit_denominator(10000)
            s6f = Fraction(s6).limit_denominator(10000)
            s8f = Fraction(s8).limit_denominator(10000)
            pred = coeffs[0]*s4f + coeffs[1]*s6f + coeffs[2]*s8f + coeffs[3]
            err = abs(float(pred) - seen[(s4, s6, s8)])
            max_err = max(max_err, err)
        print(f"\n  Max error = {max_err:.6f}")

        # Try to simplify: common denominator
        denoms = [c.denominator for c in coeffs]
        from math import gcd
        from functools import reduce
        lcm = reduce(lambda a, b: a * b // gcd(a, b), denoms)
        print(f"\n  Common denominator: {lcm}")
        print(f"  {lcm}*H = {lcm*coeffs[0]}*S4 + {lcm*coeffs[1]}*S6 + {lcm*coeffs[2]}*S8 + {lcm*coeffs[3]}")

        return coeffs


def theoretical_derivation(p):
    """Derive the formula theoretically.

    H = 1 + 2*sum_k c_k + 4*alpha_2 + 8*alpha_3 + ...

    Since c_k = (1/k)[m^k + 2*sum Re(z^k)], we can express sum_k c_k
    in terms of moments S_{2j}.

    For regular tournaments on Z_p:
      sum_{k odd, 3<=k<=p} c_k = alpha_1

    And:
      alpha_2 = disj3 = (p/4)*E + C = (p/4)(S4/p + const) + C' = S4/4 + const

    For alpha_3: involves triple-disjointness of cycles, which depends on
    the full interaction structure of Omega(T).

    The key insight: at p=7, all cycles are 3-cycles or 7-cycles.
    3-cycles are determined by c3 (constant).
    7-cycles are determined by c7 (depends on S4 only since p=7, floor((7-3)/2)=2
    parameters but m=3 so S6 is the last independent moment... hmm).

    Actually at p=7: m=3, so D has only 3 values. S2, S4, S6 are functions
    of (D1^2, D2^2, D3^2) but S2 is constant. So effectively 2 free parameters.
    But c5 depends on S4 only, and c7 on (S4, S6). Since c7 = p - sum_{k<7} c_k
    (each permutation contributes to exactly one cycle length... wait, that's
    not true for directed cycles in tournaments).

    Actually, sum_k k*c_k != n! or anything simple. Let me just compute.
    """
    m = (p - 1) // 2
    print(f"\n{'='*70}")
    print(f"THEORETICAL STRUCTURE at p={p}")
    print(f"{'='*70}")

    # Re(z^k) = sum_{j even, 0<=j<=k} C(k,j) (-1/2)^{k-j} (-1)^{j/2} D^j
    # The coefficient of D^{2r} in Re(z^k) is:
    # C(k, 2r) * (-1/2)^{k-2r} * (-1)^r

    print(f"\n  Coefficients of D^{{2r}} in Re(z^k):")
    print(f"  {'k':>3} {'D^0':>12} {'D^2':>12} {'D^4':>12} {'D^6':>12} {'D^8':>12} {'D^10':>12}")

    max_k = min(p, 13)
    coeff_table = {}  # coeff_table[k][r] = coefficient of D^{2r} in Re(z^k)
    for k in range(3, max_k + 1, 2):
        coeff_table[k] = {}
        vals = []
        for r in range(k // 2 + 1):
            j = 2 * r
            from math import comb
            c = Fraction(comb(k, j)) * Fraction(-1, 2)**(k - j) * (-1)**r
            coeff_table[k][r] = c
            vals.append(f"{float(c):12.6f}")
        while len(vals) < 6:
            vals.append(f"{'---':>12}")
        print(f"  {k:3d} {'  '.join(vals)}")

    # c_k = (1/k) * [m^k + 2 * sum_{t=1}^m sum_r coeff[k][r] * D_t^{2r}]
    #      = (1/k) * [m^k + 2 * sum_r coeff[k][r] * S_{2r}]
    # where S_0 = m (number of half-frequencies)

    print(f"\n  c_k as function of moments:")
    for k in range(3, max_k + 1, 2):
        terms = []
        for r in sorted(coeff_table[k].keys()):
            c = coeff_table[k][r]
            if c != 0:
                if r == 0:
                    terms.append(f"{float(2*c/k):.6f}*m")
                else:
                    terms.append(f"{float(2*c/k):.6f}*S{2*r}")
        const = Fraction(1, k) * Fraction(m)**k
        print(f"    c_{k} = {float(const):.2f} + {' + '.join(terms)}")

    # Now sum: alpha_1 = sum c_k
    print(f"\n  alpha_1 = sum c_k:")
    total_const = Fraction(0)
    total_S = defaultdict(Fraction)
    for k in range(3, p + 1, 2):
        total_const += Fraction(1, k) * Fraction(m)**k
        for r in range(k // 2 + 1):
            j = 2 * r
            from math import comb
            c = Fraction(comb(k, j)) * Fraction(-1, 2)**(k - j) * (-1)**r
            if r == 0:
                total_const += Fraction(2, k) * c * m
            else:
                total_S[2*r] += Fraction(2, k) * c

    print(f"    alpha_1 = {float(total_const):.4f}", end='')
    for sr in sorted(total_S.keys()):
        print(f" + {float(total_S[sr]):.6f}*S{sr}", end='')
    print()

    # H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...
    # alpha_2 = disj3 = (p/4)*E + C = S4/4 + const (since E = S4/p + const')

    # Wait: E*p = m^4 + m/8 + mp/4 + 2*S4  => S4 = (E*p - const)/2
    # So alpha_2 = (p/4)*E + C = (p/4)*(S4/p + const) + C
    #            = S4/4 + const' + C
    # But actually S4 here is sum_{t=1}^m D_t^4, and E includes the S_hat(0) term.
    # Let me just use: alpha_2 depends on S4 only (from THM-156).

    print(f"\n  H structure:")
    print(f"    H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...")
    print(f"    alpha_1 depends on: S4, S6, ..., S{{p-1}} (up to (p-3)/2 moments)")
    print(f"    alpha_2 depends on: S4 only")
    print(f"    alpha_3 depends on: ??? (cycle interaction structure)")

    # The fact that H = f(S4, S6, S8) at p=11 with only 3 parameters
    # means alpha_3 (and higher) contribute at most S8-level moments.

    # At p=11: alpha_1 depends on S4, S6, S8, S10 (from c_5, c_7, c_9, c_11)
    # alpha_2 depends on S4 only
    # If H = f(S4, S6, S8) then alpha_3 must combine with alpha_1 to cancel S10.

    # Specifically: 2*alpha_1 + 8*alpha_3 = 2*(sum c_k) + 8*alpha_3
    # The S10 coefficient in 2*alpha_1 is 2*coeff_S10_in_alpha1
    # For this to cancel: 8*d(alpha_3)/d(S10) = -2*coeff_S10_in_alpha1

    if p == 11:
        s10_coeff_alpha1 = float(total_S.get(10, 0))
        print(f"\n  At p=11:")
        print(f"    S10 coefficient in alpha_1 = {s10_coeff_alpha1:.6f}")
        print(f"    2 * S10 coeff in alpha_1 = {2*s10_coeff_alpha1:.6f}")
        print(f"    This must be cancelled by 8*(d alpha_3/d S10)")
        print(f"    => d alpha_3/d S10 = {-2*s10_coeff_alpha1/8:.6f}")


# ================================================================
# MAIN
# ================================================================

print("=" * 70)
print("H MOMENT FORMULA: Exact coefficients and theoretical derivation")
print("=" * 70)

for p in [7, 11]:
    exact_coefficients(p)

for p in [7, 11, 13]:
    theoretical_derivation(p)

print("\nDONE.")
