#!/usr/bin/env python3
"""
true_deg2_walsh.py -- CORRECT degree-2 Walsh of D^{2n} via full 2^m transform

BUG FOUND in d2n_sign_analysis.py: the extraction formula
  coeff_xy = (1/4)[f(+,+) - f(+,-) - f(-,+) + f(-,-)]
with other sigmas fixed at +1 conflates degree-2 with degree-4+ Walsh.

At p=7 (m=3): no degree-4 exists, so extraction was correct.
At p>=11 (m>=5): degree-4 contamination made ALL D^{2n} appear positive.

The TRUE degree-2 Walsh of D^4 is -2*T_pure (algebraically proved),
which is nonzero ONLY for 3-resonant pairs (3a=+/-b mod p).

This script uses full 2^m Walsh transform for exact results.

Author: kind-pasteur-2026-03-12-S60
"""

import math
from collections import defaultdict


def legendre(a, p):
    if a % p == 0:
        return 0
    return 1 if pow(a, (p-1)//2, p) == 1 else -1


def compute_all_D2n(p, max_n=5):
    """For each sigma in {+/-1}^m, compute sum_{k=1}^{p-1} D_k^{2n} for n=1..max_n.

    Returns dict: sigma_bits -> {n: value}
    """
    m = (p - 1) // 2
    results = {}

    for bits in range(1 << m):
        sigma = [(1 if bits & (1 << i) else -1) for i in range(m)]

        # Compute D_k for each k
        D_vals = {}
        for k in range(1, p):
            D_vals[k] = sum(
                sigma[l] * math.sin(2 * math.pi * k * (l+1) / p)
                for l in range(m)
            )

        # Sum D_k^{2n} for each n
        power_sums = {}
        for n in range(1, max_n + 1):
            power_sums[n] = sum(D_vals[k] ** (2*n) for k in range(1, p))

        results[bits] = power_sums

    return results


def walsh_deg2(values, m, i, j):
    """Compute degree-2 Walsh coefficient at {i,j} from function values.

    h_hat[{i,j}] = (1/2^m) sum_sigma f(sigma) * sigma_i * sigma_j
    """
    total = 0
    for bits in range(1 << m):
        si = 1 if bits & (1 << i) else -1
        sj = 1 if bits & (1 << j) else -1
        total += values[bits] * si * sj
    return total / (1 << m)


def main():
    print("=" * 70)
    print("TRUE DEGREE-2 WALSH OF D^{2n} (FULL 2^m TRANSFORM)")
    print("=" * 70)

    # PART 1: Correct degree-2 Walsh for each D^{2n}
    for p in [7, 11, 13]:
        m = (p - 1) // 2
        max_n = min(5, (p-1)//2 + 1)
        print(f"\n{'='*60}")
        print(f"p={p} ({p%4} mod 4), m={m}")
        print(f"{'='*60}")

        # Compute D^{2n} sums for all orientations
        all_D2n = compute_all_D2n(p, max_n)

        # For each n, extract the degree-2 Walsh
        for n in range(1, max_n + 1):
            print(f"\n  D^{2*n} degree-2 Walsh:")
            # Build function values for this n
            fn_vals = {bits: all_D2n[bits][n] for bits in range(1 << m)}

            match_chi = 0
            total_nonzero = 0
            total_pairs = 0

            for i in range(m):
                for j in range(i+1, m):
                    a, b = i+1, j+1  # chord values (1-indexed)
                    chi_ab = legendre(a * b, p)
                    w = walsh_deg2(fn_vals, m, i, j)
                    sign_w = 1 if w > 0.01 else (-1 if w < -0.01 else 0)

                    # 3-resonance check
                    res = []
                    if (3*a - b) % p == 0: res.append("3a=b")
                    if (3*a + b) % p == 0: res.append("3a=-b")
                    if (a - 3*b) % p == 0: res.append("a=3b")
                    if (a + 3*b) % p == 0: res.append("a=-3b")
                    res_str = ",".join(res) if res else "-"

                    total_pairs += 1
                    if sign_w != 0:
                        total_nonzero += 1
                        if sign_w == chi_ab:
                            match_chi += 1

                    if abs(w) > 0.01 or p <= 11:  # print all for small p
                        print(f"    ({a},{b}): W={w:>12.4f}, sign={sign_w:+d}, "
                              f"chi={chi_ab:+d}, match={sign_w==chi_ab if sign_w!=0 else 'N/A':>5}, "
                              f"res={res_str}")

            print(f"    SUMMARY: sign=chi in {match_chi}/{total_nonzero} nonzero "
                  f"(of {total_pairs} pairs)")

    # PART 2: Verify D^4 = -2*T_pure algebraically
    print(f"\n{'='*60}")
    print("PART 2: VERIFY D^4 TRUE DEG-2 = -2*T_PURE")
    print(f"{'='*60}")

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        print(f"\n  p={p}:")

        all_D2n = compute_all_D2n(p, 2)
        fn_D4 = {bits: all_D2n[bits][2] for bits in range(1 << m)}

        for i in range(m):
            for j in range(i+1, m):
                a, b = i+1, j+1
                w_true = walsh_deg2(fn_D4, m, i, j)

                # T_pure = sum_k 4*(Sa*Sb^3 + Sa^3*Sb)
                T_pure = 0
                for k in range(1, p):
                    Sa = math.sin(2*math.pi*k*a/p)
                    Sb = math.sin(2*math.pi*k*b/p)
                    T_pure += 4 * (Sa * Sb**3 + Sa**3 * Sb)

                predicted = -2 * T_pure
                match = abs(w_true - predicted) < 0.01

                if abs(w_true) > 0.01 or p <= 11:
                    print(f"    ({a},{b}): true_deg2={w_true:>10.4f}, "
                          f"-2*T_pure={predicted:>10.4f}, match={match}")

    # PART 3: D^6 true degree-2 via algebraic decomposition
    print(f"\n{'='*60}")
    print("PART 3: D^6 TRUE DEGREE-2 WALSH")
    print(f"{'='*60}")
    print("D^6 = (sum S_l sigma_l)^6. Degree-2 at {a,b}:")
    print("Patterns with n_a odd, n_b odd, sum=6:")
    print("  (5,1,0): 6*Sa^5*Sb")
    print("  (1,5,0): 6*Sa*Sb^5")
    print("  (3,3,0): 20*Sa^3*Sb^3")
    print("  (3,1,2c): 60*Sa^3*Sb*Sc^2")
    print("  (1,3,2c): 60*Sa*Sb^3*Sc^2")
    print("  (1,1,4c): 30*Sa*Sb*Sc^4")
    print("  (1,1,2c,2d): 180*Sa*Sb*Sc^2*Sd^2")
    print()
    print("But (1,1,2c,2d) with c!=d is degree-4 Walsh! (4 distinct sigma vars)")
    print("True degree-2 excludes (1,1,2c,2d) cross-terms.")

    for p in [7, 11]:
        m = (p - 1) // 2
        print(f"\n  p={p}:")

        all_D2n = compute_all_D2n(p, 3)
        fn_D6 = {bits: all_D2n[bits][3] for bits in range(1 << m)}

        for i in range(m):
            for j in range(i+1, m):
                a, b = i+1, j+1
                w_true = walsh_deg2(fn_D6, m, i, j)

                # Compute the "pure + single-c" contributions (true degree-2)
                T_d6_pure = 0
                T_d6_single_c = 0
                for k in range(1, p):
                    Sa = math.sin(2*math.pi*k*a/p)
                    Sb = math.sin(2*math.pi*k*b/p)

                    # (5,1) + (1,5) + (3,3)
                    T_d6_pure += (6*Sa**5*Sb + 6*Sa*Sb**5 + 20*Sa**3*Sb**3)

                    # (3,1,2c) + (1,3,2c) + (1,1,4c) — single extra var
                    for cl in range(1, m+1):
                        if cl == a or cl == b:
                            continue
                        Sc = math.sin(2*math.pi*k*cl/p)
                        T_d6_single_c += (60*Sa**3*Sb*Sc**2 + 60*Sa*Sb**3*Sc**2
                                          + 30*Sa*Sb*Sc**4)

                T_d6_total = T_d6_pure + T_d6_single_c
                chi_ab = legendre(a*b, p)

                print(f"    ({a},{b}): true_deg2={w_true:>12.4f}, "
                      f"pure+single={T_d6_total:>12.4f}, "
                      f"diff={w_true-T_d6_total:>10.4f}, chi={chi_ab:+d}")

    # PART 4: Gauss sum factorization with CORRECT degree-2
    print(f"\n{'='*60}")
    print("PART 4: GAUSS SUM FACTORIZATION (CORRECTED)")
    print(f"{'='*60}")
    print("g(a) = chi(a)*sqrt(p) at p=3 mod 4")
    print("g(a)*g(b) = chi(ab)*p")
    print()
    print("TRUE_deg2(D^{2n}) / (chi(ab)*p) should reveal the sign mechanism")

    for p in [7, 11, 19]:
        if p > 13:
            continue  # 2^9 = 512 is feasible
        m = (p - 1) // 2
        chi_neg1 = legendre(-1, p)
        max_n = min(5, m+1) if p <= 11 else 3

        print(f"\n  p={p} ({p%4} mod 4), chi(-1)={chi_neg1:+d}:")

        all_D2n = compute_all_D2n(p, max_n)

        for n in range(1, max_n + 1):
            fn = {bits: all_D2n[bits][n] for bits in range(1 << m)}
            print(f"\n    D^{2*n} true degree-2 / (chi(ab)*p):")

            ratios = []
            for i in range(m):
                for j in range(i+1, m):
                    a, b = i+1, j+1
                    chi_ab = legendre(a*b, p)
                    w = walsh_deg2(fn, m, i, j)

                    if abs(w) > 0.01 and chi_ab != 0:
                        ratio = w / (chi_ab * p)
                        ratios.append(ratio)
                        print(f"      ({a},{b}): W={w:>10.4f}, chi={chi_ab:+d}, "
                              f"ratio={ratio:>10.6f}")
                    elif abs(w) <= 0.01:
                        print(f"      ({a},{b}): W={w:>10.4f} (~ 0), chi={chi_ab:+d}")

            if ratios:
                r0 = ratios[0]
                is_const = all(abs(r - r0) < 0.001 for r in ratios)
                all_same_sign = all(r > 0 for r in ratios) or all(r < 0 for r in ratios)
                print(f"      Constant ratio? {is_const} (first={r0:.6f})")
                print(f"      All same sign? {all_same_sign}")

    # PART 5: Full h_hat decomposition by trace power
    print(f"\n{'='*60}")
    print("PART 5: FULL H_HAT DECOMPOSITION BY TRACE POWER")
    print(f"{'='*60}")
    print("h_hat[{a,b}] = sum_{k even} weight_k * (1/k) * deg2(tr(A^k))")
    print("where deg2(tr(A^k)) = sum_n coeff_{k,n} * deg2(sum_freq D^{2n})")
    print()
    print("For circulant tournament, this simplifies to understanding")
    print("deg2(sum_freq D^{2n}) for each n.")

    for p in [7, 11]:
        m = (p - 1) // 2
        print(f"\n  p={p}, m={m}:")

        # Compute H for all orientations (brute force)
        pairs = [(s, p - s) for s in range(1, m + 1)]

        def build_adj_from_bits(bits):
            S = sorted(pairs[i][0] if bits & (1 << i) else pairs[i][1] for i in range(m))
            A = [[0]*p for _ in range(p)]
            for v in range(p):
                for s in S:
                    A[v][(v + s) % p] = 1
            return A

        def count_hp(A, n):
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

        H_vals = {}
        for bits in range(1 << m):
            A = build_adj_from_bits(bits)
            H_vals[bits] = count_hp(A, p)

        # Degree-2 Walsh of H
        for i in range(m):
            for j in range(i+1, m):
                a, b = i+1, j+1
                h_hat = walsh_deg2(H_vals, m, i, j)
                chi_ab = legendre(a*b, p)
                sign_h = 1 if h_hat > 0 else (-1 if h_hat < 0 else 0)

                # Now decompose: which D^{2n} deg-2 contributes most?
                all_D2n = compute_all_D2n(p, min(5, m+1))
                contributions = {}
                for n in range(1, min(5, m+1) + 1):
                    fn = {bits: all_D2n[bits][n] for bits in range(1 << m)}
                    contributions[n] = walsh_deg2(fn, m, i, j)

                print(f"    ({a},{b}): h_hat={h_hat:>10.2f}, sign={sign_h:+d}, "
                      f"chi={chi_ab:+d}, match={sign_h==chi_ab}")
                for n in sorted(contributions):
                    if abs(contributions[n]) > 0.001:
                        print(f"        D^{2*n} true deg-2: {contributions[n]:>12.4f}")

    # PART 6: The big question at p=19
    print(f"\n{'='*60}")
    print("PART 6: p=19 (m=9, 2^9=512 orientations)")
    print(f"{'='*60}")

    p = 19
    m = (p - 1) // 2
    max_n = 4

    all_D2n = compute_all_D2n(p, max_n)

    print(f"\n  D^{{2n}} true degree-2 at p={p}:")
    for n in range(1, max_n + 1):
        fn = {bits: all_D2n[bits][n] for bits in range(1 << m)}
        print(f"\n  D^{2*n}:")

        match_chi = 0
        total_nz = 0
        total = 0
        for i in range(m):
            for j in range(i+1, m):
                a, b = i+1, j+1
                chi_ab = legendre(a*b, p)
                w = walsh_deg2(fn, m, i, j)
                sign_w = 1 if w > 0.01 else (-1 if w < -0.01 else 0)

                total += 1
                if sign_w != 0:
                    total_nz += 1
                    if sign_w == chi_ab:
                        match_chi += 1

                # Print for nonzero or first few
                if abs(w) > 0.01:
                    res = []
                    if (3*a - b) % p == 0: res.append("3a=b")
                    if (3*a + b) % p == 0: res.append("3a=-b")
                    if (a - 3*b) % p == 0: res.append("a=3b")
                    if (a + 3*b) % p == 0: res.append("a=-3b")
                    print(f"    ({a},{b}): W={w:>12.4f}, chi={chi_ab:+d}, "
                          f"sign={sign_w:+d}, res={','.join(res) if res else '-'}")

        print(f"    SUMMARY: sign=chi in {match_chi}/{total_nz} nonzero "
              f"({total} pairs, {total-total_nz} zero)")


if __name__ == '__main__':
    main()
