#!/usr/bin/env python3
"""
thm143_formula_verify.py -- Verify THM-143 closed-form formula against brute force

THM-143 (PROVED): For Interval tournament on Z_p with S={1,...,m}, m=(p-1)/2:

    co_occ_k(d) = C(2m-1, k-2) - C(m-1, k-2) - (m-d) * C(m-2, k-3)

    Slope:     b_k = C(m-2, k-3)
    Intercept: a_k = C(2m-1, k-2) - C(m-1, k-2) - m * C(m-2, k-3)

Proof uses: composition counting + Vandermonde identity + absorption identity.

This script verifies the formula against ALL brute-force data from previous runs.

Author: kind-pasteur-2026-03-12-S59c
"""

from math import comb


def co_occ_formula(p, k, d):
    """THM-143: exact co_occ_k(d) for Interval tournament."""
    m = (p - 1) // 2
    d_eff = min(d, p - d)  # symmetry
    return comb(2*m - 1, k - 2) - comb(m - 1, k - 2) - (m - d_eff) * comb(m - 2, k - 3)


def slope_formula(p, k):
    m = (p - 1) // 2
    return comb(m - 2, k - 3)


def intercept_formula(p, k):
    m = (p - 1) // 2
    return comb(2*m - 1, k - 2) - comb(m - 1, k - 2) - m * comb(m - 2, k - 3)


def main():
    print("=" * 70)
    print("THM-143 FORMULA VERIFICATION")
    print("=" * 70)

    # All brute-force data collected from previous runs
    # Format: (p, k, profile_values_for_d_1_to_m)
    data = [
        # p=7
        (7, 3, [1, 2, 3]),
        (7, 5, [10, 10, 10]),
        # p=11
        (11, 3, [1, 2, 3, 4, 5]),
        (11, 5, [68, 71, 74, 77, 80]),
        (11, 7, [126, 126, 126, 126, 126]),
        (11, 9, [36, 36, 36, 36, 36]),
        # p=13
        (13, 3, [1, 2, 3, 4, 5, 6]),
        (13, 5, [125, 131, 137, 143, 149, 155]),
        (13, 7, [456, 457, 458, 459, 460, 461]),
        (13, 9, [330, 330, 330, 330, 330, 330]),
        # p=17
        (17, 3, [1, 2, 3, 4, 5, 6, 7, 8]),
        (17, 5, [315, 330, 345, 360, 375, 390, 405, 420]),
        (17, 7, [2877, 2892, 2907, 2922, 2937, 2952, 2967, 2982]),
        (17, 9, [6427, 6428, 6429, 6430, 6431, 6432, 6433, 6434]),
        (17, 11, [5005, 5005, 5005, 5005, 5005, 5005, 5005, 5005]),
        (17, 13, [1365, 1365, 1365, 1365, 1365, 1365, 1365, 1365]),
        # p=19
        (19, 3, [1, 2, 3, 4, 5, 6, 7, 8, 9]),
        (19, 5, [456, 477, 498, 519, 540, 561, 582, 603, 624]),
        (19, 7, [5852, 5887, 5922, 5957, 5992, 6027, 6062, 6097, 6132]),
        (19, 9, [19384, 19391, 19398, 19405, 19412, 19419, 19426, 19433, 19440]),
        (19, 11, [24310, 24310, 24310, 24310, 24310, 24310, 24310, 24310, 24310]),
        (19, 13, [12376, 12376, 12376, 12376, 12376, 12376, 12376, 12376, 12376]),
    ]

    all_pass = True
    print(f"\n{'p':>3} {'k':>3} | {'slope':>8} {'interc':>8} | {'formula':>10} {'actual':>10} | match")
    print("-" * 70)

    for p, k, profile in data:
        m = (p - 1) // 2
        b = slope_formula(p, k)
        a = intercept_formula(p, k)

        match_all = True
        for i, d in enumerate(range(1, m + 1)):
            formula_val = co_occ_formula(p, k, d)
            actual_val = profile[i]
            if formula_val != actual_val:
                match_all = False
                all_pass = False
                print(f"  MISMATCH at p={p}, k={k}, d={d}: "
                      f"formula={formula_val}, actual={actual_val}")

        status = "PASS" if match_all else "FAIL"
        print(f"{p:>3} {k:>3} | {b:>8} {a:>8} | "
              f"{co_occ_formula(p, k, 1):>10} {profile[0]:>10} | {status}")

    print(f"\n{'='*70}")
    if all_pass:
        print("ALL 22 TEST CASES PASS -- THM-143 FULLY VERIFIED")
    else:
        print("SOME TESTS FAILED")

    # ====== PREDICTIONS FOR LARGER PRIMES ======
    print(f"\n{'='*70}")
    print("PREDICTIONS FOR LARGER PRIMES")
    print("=" * 70)

    for p in [23, 29, 31, 37, 41, 43, 47]:
        m = (p - 1) // 2
        print(f"\np={p}, m={m}:")
        for k in range(3, min(p, 20), 2):
            b = slope_formula(p, k)
            a = intercept_formula(p, k)
            if b > 0 or k <= 7:
                print(f"  k={k:>2}: slope = C({m-2},{k-3}) = {b:>10}, "
                      f"intercept = {a:>10}")
        # When does slope become 0?
        k_max = m + 1  # C(m-2, k-3) = 0 when k-3 > m-2, i.e., k > m+1
        # Actually C(m-2, k-3) = 0 when k-3 > m-2, i.e., k > m+1
        # But k must be odd, so largest odd k with nonzero slope is...
        k_crit = m + 1 if (m + 1) % 2 == 1 else m
        print(f"  Slope nonzero for k <= {k_crit} (p >= 2k-1)")

    # ====== DISJOINTNESS EXCESS FORMULA ======
    print(f"\n{'='*70}")
    print("DISJOINTNESS EXCESS FROM CO_OCC VARIANCE")
    print("=" * 70)

    print("\nFor k=3 (THM-142 check):")
    for p in [7, 11, 13, 17, 19, 23, 29]:
        m = (p - 1) // 2
        # #{ov >= 2}_Int = p * sum_{d=1}^{m} C(co_occ_3(d), 2)
        #                = p * sum_{d=1}^m C(d, 2)
        #                = p * m(m-1)(m+1)/6  [telescoping sum of C(d,2)]
        ov2_int = p * sum(comb(d, 2) for d in range(1, m + 1))
        # Paley: co_occ = (p-1)/4, so #{ov>=2} = p * m * C((p-1)/4, 2)
        if p % 4 == 3:
            paley_co = (p - 1) // 4
            ov2_pal = p * m * comb(paley_co, 2)
            excess = ov2_int - ov2_pal
            thm142 = p * (p-1) * (p+1) * (p-3) // 192
            print(f"  p={p}: excess = {excess}, THM-142 = {thm142}, match = {excess == thm142}")
        else:
            print(f"  p={p}: no Paley (p != 3 mod 4)")

    print("\nFor k=5 (NEW):")
    for p in [11, 13, 19, 23]:
        m = (p - 1) // 2
        b5 = slope_formula(p, 5)
        a5 = intercept_formula(p, 5)
        # #{ov >= 2} for Interval 5-cycles
        ov2_int = p * sum(comb(a5 + b5 * d, 2) for d in range(1, m + 1))
        # For Paley (if exists): co_occ is constant
        if p % 4 == 3:
            # Paley co_occ_5: constant = n_0 * (k-1) / (p-1) ... need to compute
            # Actually, for Paley: co_occ_5(d) = c for all d,
            # so sum of co_occ = c * (p-1) = n_0 * (k-1)
            # But we don't know n_0 for Paley without computing.
            # Let's just report Interval's #{ov>=2}
            print(f"  p={p}: Interval ov2 = {ov2_int}, slope = {b5}")
        else:
            print(f"  p={p}: Interval ov2 = {ov2_int}, slope = {b5}")

    # ====== TOTAL VERTEX SETS FROM FORMULA ======
    print(f"\n{'='*70}")
    print("TOTAL VERTEX SETS n_0 FROM FORMULA")
    print("=" * 70)
    print("n_0 = sum_{d=1}^{p-1} co_occ_k(d) / (k-1)")
    print("    = [2 * sum_{d=1}^m (a_k + b_k * d)] / (k-1)")
    print("    = 2 * [m * a_k + b_k * m*(m+1)/2] / (k-1)")
    print()

    for p in [7, 11, 13, 17, 19]:
        m = (p - 1) // 2
        print(f"p={p}, m={m}:")
        for k in range(3, min(p, 14), 2):
            a = intercept_formula(p, k)
            b = slope_formula(p, k)
            n0 = 2 * (m * a + b * m * (m + 1) // 2) // (k - 1)
            n_total = n0 * p // k
            print(f"  k={k}: a={a}, b={b}, n_0={n0}, n_total={n_total}")


if __name__ == '__main__':
    main()
