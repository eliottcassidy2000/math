#!/usr/bin/env python3
"""
zero_sum_walsh_analysis.py -- Zero-sum structure in Walsh coefficients

OPUS FINDING (S67b): Degree-4 Walsh coefficients split by
  W = #{(eps) in {+1,-1}^4 : sum eps_i * pair_vals[i] = 0 mod p}
  - W > 0 (additive relation): large positive net contribution
  - W = 0 (generic): small negative net contribution

This script:
1. Classifies ALL Walsh coefficients by zero-sum structure
2. Connects zero-sum count W to representation profile r_S
3. Shows WHY zero-sum structure drives H-maximization
4. Tests a formula: h_hat[S] = C * (-1)^{...} * W(S)?
5. Investigates the degree-2 analogue

KEY INSIGHT: If pair_vals = {a_1, ..., a_d} with a_i in [1,m],
then W = #{eps : sum eps_i * a_i = 0 mod p}.
This counts the number of representations of 0 in the sumset
{+a_1, -a_1} + ... + {+a_d, -a_d}.
For d=2: W > 0 iff a_i + a_j = 0 or a_i - a_j = 0 mod p.
Since a_i in [1,m], a_i - a_j = 0 only if i=j.
And a_i + a_j = 0 mod p iff a_i + a_j = p, which happens iff a_i = p - a_j.
But a_i, a_j in [1, m], so a_i + a_j <= 2m = p-1 < p. Thus W = 0 for d=2
UNLESS a_i = a_j (impossible for a set).

Wait, that's not right. The sum includes +/- signs:
{+a_1, -a_1} + {+a_2, -a_2} = {a_1+a_2, a_1-a_2, -a_1+a_2, -a_1-a_2}
We need this to contain 0 mod p.
- a_1 + a_2 = 0 mod p: impossible (both positive, < p)
- a_1 - a_2 = 0: iff a_1 = a_2 (impossible for distinct)
- -(a_1 + a_2) = 0 mod p: iff a_1 + a_2 = p, iff a_1 = p - a_2
  But a_1 in [1,m], a_2 in [1,m], so a_1 + a_2 in [2, 2m] = [2, p-1] < p.
  So a_1 + a_2 = p is impossible.
- -a_1 + a_2 = 0: iff a_1 = a_2.

So for degree-2 with pair indices i,j (pair_vals a_i, a_j in [1,m]):
W = 0 always! (since a_i != a_j and a_i + a_j < p)

But opus found that the SIGN of degree-2 coefficients depends on
whether pair_sum = a_i + a_j falls in QR or NQR...

Actually, the zero-sum analysis uses FULL pairs {+a, -a} = {a, p-a}.
Let me reconsider. The pair (i, p-i) contributes the FULL element,
and the orientation sigma_i chooses which element to include.
The Walsh coefficient h_hat[{i,j}] involves flipping sigma_i and sigma_j.

For degree-4: four pair indices {i,j,k,l} with pair_vals {a_i, a_j, a_k, a_l}.
W counts #{eps : eps_1*a_i + eps_2*a_j + eps_3*a_k + eps_4*a_l = 0 mod p}
where each eps in {+1, -1}. Total 2^4 = 16 possibilities.

This is related to the number of ways to form 0 from {a_i, p-a_i} choices.

Author: kind-pasteur-2026-03-12-S60
"""

import cmath
import math
from itertools import combinations
from collections import defaultdict


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i + s) % p] = 1
    return A


def count_ham_paths(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    dp[mask | (1 << w)][w] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))


def walsh_decomposition(p, pairs):
    m = len(pairs)
    n_orient = 1 << m

    H_vals = {}
    for bits in range(n_orient):
        sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
        S = sorted(pairs[i][0] if sigma[i] == 1 else pairs[i][1] for i in range(m))
        A = build_adj(p, S)
        H = count_ham_paths(A, p)
        H_vals[bits] = H

    h_hat = {}
    for S_bits in range(n_orient):
        S_set = frozenset(i for i in range(m) if S_bits & (1 << i))
        val = 0
        for bits in range(n_orient):
            chi = 1
            for i in S_set:
                if not (bits & (1 << i)):
                    chi *= -1
            val += H_vals[bits] * chi
        h_hat[S_set] = val / n_orient

    return h_hat, H_vals


def zero_sum_count(vals, p):
    """Count #{eps in {+1,-1}^d : sum eps_i * vals[i] = 0 mod p}."""
    d = len(vals)
    W = 0
    for bits in range(1 << d):
        eps = [(1 if bits & (1 << j) else -1) for j in range(d)]
        s = sum(e * v for e, v in zip(eps, vals)) % p
        if s == 0:
            W += 1
    return W


def main():
    print("=" * 70)
    print("ZERO-SUM WALSH ANALYSIS")
    print("=" * 70)

    # PART 1: Full degree-2 and degree-4 classification
    for p in [7, 11, 13]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]
        QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)

        print(f"\n{'='*60}")
        print(f"p={p}, m={m}")
        print(f"{'='*60}")

        h_hat, H_vals = walsh_decomposition(p, pairs)

        # DEGREE 2 ANALYSIS
        print(f"\n  --- DEGREE 2 ---")
        print(f"  {'pair':>10} {'pair_vals':>12} {'sum':>5} {'diff':>5} {'W':>3} {'h_hat':>10} {'|h_hat|':>10}")

        d2_by_W = defaultdict(list)
        d2_by_sum = defaultdict(list)
        for S_set, val in sorted(h_hat.items(), key=lambda x: tuple(sorted(x[0]))):
            if len(S_set) != 2:
                continue
            indices = sorted(S_set)
            pair_vals = [pairs[i][0] for i in indices]
            pv_sum = sum(pair_vals) % p
            pv_diff = abs(pair_vals[0] - pair_vals[1]) % p

            # Zero-sum: #{eps : eps1*a + eps2*b = 0 mod p}
            W = zero_sum_count(pair_vals, p)

            # Also check: a+b mod p, a-b mod p, -a+b mod p, -a-b mod p
            sums_mod_p = [(pair_vals[0] + pair_vals[1]) % p,
                          (pair_vals[0] - pair_vals[1]) % p,
                          (-pair_vals[0] + pair_vals[1]) % p,
                          (-pair_vals[0] - pair_vals[1]) % p]

            d2_by_W[W].append(val)
            d2_by_sum[pv_sum].append(val)

            print(f"  ({indices[0]+1},{indices[1]+1}){'':<5} {str(pair_vals):>12} "
                  f"{pv_sum:>5} {pv_diff:>5} {W:>3} {val:>10.1f} {abs(val):>10.1f}")

        print(f"\n  Degree-2 by W:")
        for W in sorted(d2_by_W):
            vals = d2_by_W[W]
            print(f"    W={W}: count={len(vals)}, net={sum(vals):.1f}, "
                  f"magnitudes={sorted(set(round(abs(v), 1) for v in vals))}")

        # DEGREE 4 ANALYSIS
        if m >= 4:
            print(f"\n  --- DEGREE 4 ---")
            print(f"  {'set':>15} {'pair_vals':>16} {'sum':>5} {'W':>3} {'h_hat':>10}")

            d4_by_W = defaultdict(list)
            for S_set, val in sorted(h_hat.items(), key=lambda x: tuple(sorted(x[0]))):
                if len(S_set) != 4:
                    continue
                indices = sorted(S_set)
                pair_vals = [pairs[i][0] for i in indices]
                pv_sum = sum(pair_vals) % p
                W = zero_sum_count(pair_vals, p)
                d4_by_W[W].append(val)

                print(f"  {set(i+1 for i in indices)!s:>15} {str(pair_vals):>16} "
                      f"{pv_sum:>5} {W:>3} {val:>10.1f}")

            print(f"\n  Degree-4 by W:")
            for W in sorted(d4_by_W):
                vals = d4_by_W[W]
                print(f"    W={W}: count={len(vals)}, net={sum(vals):.1f}, "
                      f"magnitudes={sorted(set(round(abs(v), 1) for v in vals))}")

        # DEGREE 2: Test sign law
        print(f"\n  --- DEGREE-2 SIGN LAW ANALYSIS ---")
        print(f"  Testing: sign(h_hat[{{i,j}}]) depends on gap structure")

        # For degree 2, the pair_sum = a_i + a_j
        # The key question: what determines sign(h_hat)?
        # Hypothesis: sign related to whether pair_sum is "close to m" or "far from m"
        for S_set, val in sorted(h_hat.items(), key=lambda x: tuple(sorted(x[0]))):
            if len(S_set) != 2 or abs(val) < 0.01:
                continue
            indices = sorted(S_set)
            pair_vals = [pairs[i][0] for i in indices]
            pv_sum = sum(pair_vals)
            # Distance from m+1 (= (p+1)/2)
            dist_from_m1 = abs(pv_sum - (m + 1))
            # Is pair_sum in [2, m+1] or [m+2, 2m]?
            half = "lower" if pv_sum <= m + 1 else "upper"
            sign = "+" if val > 0 else "-"

            # Check: pair_sum mod p, and its Legendre symbol
            chi_sum = pow(pv_sum % p, (p-1)//2, p) if (pv_sum % p) != 0 else 0
            chi_str = "+1" if chi_sum == 1 else "-1" if chi_sum == p - 1 else "0"

        # NEW: For degree 2, check product of pair indices
        print(f"\n  --- DEGREE-2: PRODUCT LAW ---")
        print(f"  Testing: h_hat[{{i,j}}] ~ chi(i*j) * something?")

        for S_set, val in sorted(h_hat.items(), key=lambda x: tuple(sorted(x[0]))):
            if len(S_set) != 2 or abs(val) < 0.01:
                continue
            indices = sorted(S_set)
            a, b = pairs[indices[0]][0], pairs[indices[1]][0]
            prod_ab = (a * b) % p
            chi_prod = pow(prod_ab, (p-1)//2, p)
            chi_prod = 1 if chi_prod == 1 else -1
            sign = 1 if val > 0 else -1
            match = "MATCH" if sign == chi_prod else "FAIL"
            print(f"    ({indices[0]+1},{indices[1]+1}): a*b={a}*{b}={prod_ab} mod {p}, "
                  f"chi(a*b)={chi_prod:+d}, sign(h)={sign:+d}, {match}")

    # PART 2: The number of magnitudes and their ratios
    print(f"\n{'='*60}")
    print("PART 2: MAGNITUDE STRUCTURE")
    print("=" * 60)

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]

        h_hat, _ = walsh_decomposition(p, pairs)

        print(f"\n  p={p}:")
        for d in [2, 4]:
            mags = sorted(set(round(abs(v), 2) for S_set, v in h_hat.items()
                              if len(S_set) == d and abs(v) > 0.01))
            print(f"    Degree {d}: distinct |h_hat| = {mags}")
            if len(mags) >= 2:
                ratios = [round(mags[i]/mags[0], 4) for i in range(len(mags))]
                print(f"    Ratios to smallest: {ratios}")
                # Check if ratios are rational with small denominator
                for r in ratios[1:]:
                    for denom in range(1, 20):
                        if abs(r * denom - round(r * denom)) < 0.01:
                            print(f"      {r:.4f} = {round(r * denom)}/{denom}")
                            break

    # PART 3: The Gauss sum connection to degree-2 signs
    print(f"\n{'='*60}")
    print("PART 3: GAUSS SUM CONNECTION")
    print("=" * 60)
    print("For circulant tournaments, S_hat(k) = sum_{s in S} omega^{ks}")
    print("The imaginary part D(k;sigma) = sum_{i=1}^m sigma_i * sin(2*pi*k*i/p)")
    print("h_hat[{i,j}] involves the product D(k;i)*D(k;j) summed over k")

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]
        omega = cmath.exp(2j * cmath.pi / p)

        print(f"\n  p={p}:")

        # Compute D_k(i) = sin(2*pi*k*i/p) for k=1..m, i=1..m
        D = [[0.0] * m for _ in range(m)]  # D[k][i], k=0..m-1 -> frequency k+1
        for ki in range(m):
            k = ki + 1
            for ii in range(m):
                i_val = pairs[ii][0]  # = ii + 1
                D[ki][ii] = math.sin(2 * math.pi * k * i_val / p)

        # The "J matrix" in Walsh space should relate to D^T D
        # h_hat[{i,j}] ~ sum_k f(Q_k) * D[k][i] * D[k][j] for some function f
        # Let's compute sum_k D[k][i]*D[k][j] and compare with h_hat[{i,j}]

        DD = [[0.0] * m for _ in range(m)]
        for i in range(m):
            for j in range(m):
                DD[i][j] = sum(D[ki][i] * D[ki][j] for ki in range(m))

        print(f"    D^T D matrix (sum_k sin(2pi*k*a/p)*sin(2pi*k*b/p)):")
        for i in range(m):
            row = "      "
            for j in range(m):
                row += f"{DD[i][j]:>8.3f}"
            print(row)

        # Compare with h_hat degree-2
        print(f"\n    Comparison h_hat[{{i,j}}] vs D^T D [i,j]:")
        h_hat_computed, _ = walsh_decomposition(p, pairs)
        for S_set, val in sorted(h_hat_computed.items(), key=lambda x: tuple(sorted(x[0]))):
            if len(S_set) != 2:
                continue
            i, j = sorted(S_set)
            dd_val = DD[i][j]
            if abs(dd_val) > 0.001:
                ratio = val / dd_val
                print(f"    ({i+1},{j+1}): h_hat={val:>10.2f}, D^TD={dd_val:>8.4f}, ratio={ratio:>10.2f}")

        # Now try weighted: sum_k Q_k * D[k][i] * D[k][j]
        S_int = list(range(1, m + 1))
        Q_vals = []
        for k in range(1, m + 1):
            val = sum(omega ** (k * s) for s in S_int)
            Q_vals.append(abs(val)**2)

        QDD = [[0.0] * m for _ in range(m)]
        for i in range(m):
            for j in range(m):
                QDD[i][j] = sum(Q_vals[ki] * D[ki][i] * D[ki][j] for ki in range(m))

        print(f"\n    Comparison h_hat[{{i,j}}] vs sum_k Q_k*D_k(i)*D_k(j):")
        for S_set, val in sorted(h_hat_computed.items(), key=lambda x: tuple(sorted(x[0]))):
            if len(S_set) != 2:
                continue
            i, j = sorted(S_set)
            qdd_val = QDD[i][j]
            if abs(qdd_val) > 0.001:
                ratio = val / qdd_val
                print(f"    ({i+1},{j+1}): h_hat={val:>10.2f}, Q*D^TD={qdd_val:>8.4f}, ratio={ratio:>10.2f}")

    # PART 4: The interval evaluation formula
    print(f"\n{'='*60}")
    print("PART 4: INTERVAL EVALUATION OF WALSH SUMS")
    print("=" * 60)
    print("At Interval (sigma = all +1):")
    print("  chi_S(sigma_Int) = 1 for all S")
    print("  So f_d(Int) = sum_{|S|=d} h_hat[S]")
    print("  = sum of ALL degree-d Walsh coefficients")

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]

        h_hat, H_vals = walsh_decomposition(p, pairs)

        print(f"\n  p={p}:")
        for d in range(m + 1):
            f_d = sum(val for S_set, val in h_hat.items() if len(S_set) == d)
            if abs(f_d) > 0.01:
                n_pos = sum(1 for S_set, val in h_hat.items() if len(S_set) == d and val > 0.01)
                n_neg = sum(1 for S_set, val in h_hat.items() if len(S_set) == d and val < -0.01)
                n_zero = sum(1 for S_set, val in h_hat.items() if len(S_set) == d and abs(val) <= 0.01)
                print(f"    deg {d}: f_d(Int) = {f_d:>12.1f}  (+:{n_pos}, -:{n_neg}, 0:{n_zero})")

        # Check: H(Int) = sum f_d(Int)
        H_int = sum(val for S_set, val in h_hat.items())
        sigma_int_bits = (1 << m) - 1
        print(f"    sum f_d = {H_int:.1f} vs H(Int) = {H_vals[sigma_int_bits]}")

    # PART 5: Explicit formula attempt
    print(f"\n{'='*60}")
    print("PART 5: EXPLICIT FORMULA FOR |h_hat|")
    print("=" * 60)
    print("At p=11: all degree-2 |h_hat| take 2 values: ~8.25 and ~272.25")
    print("At p=13: all degree-2 |h_hat| take 2 values: ~778.38 and ~2632.50")
    print("At p=11: all degree-4 |h_hat| take 1 value: ~118.25")
    print("At p=13: all degree-4 |h_hat| take 3 values: ~219.38, ~3508.38, ~5209.75")

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]
        h_hat, _ = walsh_decomposition(p, pairs)

        print(f"\n  p={p}, m={m}:")
        for d in [2, 4]:
            vals = [abs(val) for S_set, val in h_hat.items() if len(S_set) == d and abs(val) > 0.01]
            if not vals:
                continue
            vals_unique = sorted(set(round(v, 2) for v in vals))
            print(f"    deg {d}: |h_hat| values = {vals_unique}")

            # Check: are these related to p, m, or Fibonacci?
            for v in vals_unique:
                print(f"      {v:.2f} / p = {v/p:.4f}")
                print(f"      {v:.2f} / m = {v/m:.4f}")
                print(f"      {v:.2f} * 2^m / p! = {v * (2**m) / math.factorial(p):.8f}")
                # Check if v = p^k * something simple
                for k in range(1, 5):
                    r = v / p**k
                    if abs(r - round(r)) < 0.01:
                        print(f"      {v:.2f} = {round(r)} * p^{k}")
                    r = v / (p**k / 2)
                    if abs(r - round(r)) < 0.01:
                        print(f"      {v:.2f} = {round(r)} * p^{k}/2")


if __name__ == '__main__':
    main()
