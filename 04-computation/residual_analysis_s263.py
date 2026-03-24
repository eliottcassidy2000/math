#!/usr/bin/env python3
"""
residual_analysis_s263.py — Understanding the residual through factorization
opus-2026-03-23-S263

The residual R(n) = T_n - twin_SL(n) - 2E(n) = complex_SL + MW.
Sequence: 0, 2, 16, 76, 444, 2926, 31442 at n=3..9.

The factorization tells us:
  V_n = (1/n!) Σ Fix_even × 2^{k-1}
  twin_SL = (1/n!) Σ 2C(f,2) × Fix_even
  T_n = (1/n!) Σ C(f,2) × Fix_even × 2^{k-1}

So: T_n - twin_SL = (1/n!) Σ [C(f,2) × 2^{k-1} - 2C(f,2)] × Fix_even
                   = (1/n!) Σ C(f,2) × [2^{k-1} - 2] × Fix_even

And 2E = T_n - twin_SL - R, so:
  2E = (1/n!) Σ C(f,2) × [2^{k-1} - 2] × Fix_even - R

If R = 0: E = (T-twin_SL)/2 = (1/(2n!)) Σ C(f,2) × [2^{k-1} - 2] × Fix_even

The residual measures the deviation from this.

Let's analyze:
1. The ratio R/T
2. The ratio R/V
3. Whether R has its own Burnside-like structure
4. The sequence R/2: 0, 1, 8, 38, 222, 1463, 15721 — look for patterns
"""

from math import comb, factorial, gcd
from collections import Counter

def partitions(n, mx=None):
    if mx is None: mx = n
    if n == 0: yield []; return
    for f in range(min(n, mx), 0, -1):
        for r in partitions(n - f, f): yield [f] + r

def ccs_val(n, ct):
    c = Counter(ct); r = factorial(n)
    for l, k in c.items(): r //= (l ** k) * factorial(k)
    return r

def arc_orbits(ct):
    total = 0
    for i in range(len(ct)):
        total += (ct[i] - 1) // 2
        for j in range(i + 1, len(ct)):
            total += gcd(ct[i], ct[j])
    return total

if __name__ == '__main__':
    print("=" * 70)
    print("  RESIDUAL ANALYSIS")
    print("  opus-2026-03-23-S263")
    print("=" * 70)

    known_E = {3:1, 4:5, 5:30, 6:290, 7:4086, 8:91161, 9:3380751}

    res_seq = []

    print(f"\n{'n':>3} {'T':>12} {'twin_SL':>10} {'T-twin':>12} {'2E':>12} "
          f"{'R':>10} {'R/2':>8} {'R/T':>8} {'R/V':>8}")

    for n in range(3, 13):
        nf = factorial(n)
        m = comb(n, 2)

        # Compute T and twin_SL via factored Burnside
        T_sum = 0
        twin_sum = 0
        V_sum = 0

        for ct in partitions(n):
            ct_list = list(ct)
            if any(c % 2 == 0 for c in ct_list): continue
            cc = ccs_val(n, ct_list)
            ao = arc_orbits(ct_list)
            k = len(ct_list)
            f = ct_list.count(1)

            fix_even = 2 ** (ao - k + 1)
            fix_tourn = fix_even * (2 ** (k - 1))

            V_sum += cc * fix_tourn
            T_sum += cc * fix_tourn * comb(f, 2)
            twin_sum += cc * 2 * comb(f, 2) * fix_even if f >= 2 else 0

        V = V_sum // nf
        T = T_sum // nf
        twin_SL = twin_sum // nf

        T_minus_twin = T - twin_SL
        E = known_E.get(n)

        if E is not None:
            R = T_minus_twin - 2 * E
            R_half = R // 2
            res_seq.append(R_half)
            r_over_t = R / T if T > 0 else 0
            r_over_v = R / V if V > 0 else 0

            print(f"{n:>3} {T:>12} {twin_SL:>10} {T_minus_twin:>12} {2*E:>12} "
                  f"{R:>10} {R_half:>8} {r_over_t:>8.5f} {r_over_v:>8.4f}")
        else:
            E_approx = T_minus_twin // 2
            print(f"{n:>3} {T:>12} {twin_SL:>10} {T_minus_twin:>12} {'~'+str(2*E_approx):>12} "
                  f"{'?':>10} {'?':>8} {'?':>8} {'?':>8}")

    print(f"\n  R/2 sequence: {res_seq}")

    # Pattern analysis
    print(f"\n{'='*70}")
    print("  PATTERN ANALYSIS OF R/2: 0, 1, 8, 38, 222, 1463, 15721")
    print(f"{'='*70}")

    half_res = [0, 1, 8, 38, 222, 1463, 15721]
    print(f"  Ratios: ", end="")
    for i in range(1, len(half_res)):
        if half_res[i-1] > 0:
            print(f"{half_res[i]/half_res[i-1]:.3f}, ", end="")
        else:
            print(f"inf, ", end="")
    print()

    # Differences
    print(f"  First differences: {[half_res[i]-half_res[i-1] for i in range(1, len(half_res))]}")
    print(f"  Ratios of diffs: ", end="")
    diffs = [half_res[i]-half_res[i-1] for i in range(1, len(half_res))]
    for i in range(1, len(diffs)):
        if diffs[i-1] > 0:
            print(f"{diffs[i]/diffs[i-1]:.3f}, ", end="")
    print()

    # Check: is R/2 = a(n) * (n-2)! for some a(n)?
    print(f"\n  R/2 divided by (n-2)!:")
    for i, n in enumerate(range(3, 10)):
        if half_res[i] > 0:
            ratio = half_res[i] / factorial(n-2)
            print(f"    n={n}: {half_res[i]} / {factorial(n-2)} = {ratio:.4f}")

    # Check: is R/2 divisible by derangement D(n-2)?
    def derangement(nn):
        if nn == 0: return 1
        if nn == 1: return 0
        d = [1, 0]
        for kk in range(2, nn+1): d.append((kk-1)*(d[-1]+d[-2]))
        return d[nn]

    print(f"\n  R/2 divided by D(n-2):")
    for i, n in enumerate(range(3, 10)):
        d = derangement(n-2)
        if d > 0 and half_res[i] > 0:
            ratio = half_res[i] / d
            print(f"    n={n}: {half_res[i]} / D({n-2})={d} = {ratio:.4f}")

    # Is R related to the edge-orbit formula gap?
    # gap = EO - E = (n-2)! + twin_SL/2 + residual/2
    print(f"\n  Gap = EO - E = (n-2)! + twin_SL/2 + R/2:")
    for i, n in enumerate(range(3, 10)):
        E = known_E[n]
        gap_vals = [
            factorial(n-2),   # the (n-2)! from edge_orbits
            half_res[i],      # R/2
        ]
        # twin_SL/2 is what we need to add to get gap
        # gap = EO - E = T/2 + (n-2)! - E
        #     = (T - 2E)/2 + (n-2)!
        #     = (twin_SL + R)/2 + (n-2)!
        # So gap = twin_SL/2 + R/2 + (n-2)!

    # Key insight: R(n) = T(n) - twin_SL(n) - 2E(n)
    # Through factorization:
    # T(n) = (1/n!) Σ C(f,2) × Fix_even × 2^{k-1}
    # twin_SL(n) = (1/n!) Σ 2C(f,2) × Fix_even
    # T - twin_SL = (1/n!) Σ C(f,2) × Fix_even × [2^{k-1} - 2]
    # 2E = T - twin_SL - R
    # So R = (1/n!) Σ C(f,2) × Fix_even × [2^{k-1} - 2] - 2E

    # This means R captures the part of the "transition-minus-twin" that
    # doesn't correspond to new edges. It's the SELF-LOOPS that aren't twins
    # plus the MULTI-WIRE excess.

    # The factored view of T - twin_SL:
    print(f"\n{'='*70}")
    print("  T - twin_SL DECOMPOSITION BY CYCLE TYPE")
    print(f"{'='*70}")

    for n in range(3, 9):
        nf = factorial(n)
        print(f"\n  n={n}:")
        total_diff = 0
        for ct in partitions(n):
            ct_list = list(ct)
            if any(c % 2 == 0 for c in ct_list): continue
            cc = ccs_val(n, ct_list)
            ao = arc_orbits(ct_list)
            k = len(ct_list)
            f = ct_list.count(1)
            fix_even = 2 ** (ao - k + 1)

            # Contribution to T: C(f,2) × Fix_even × 2^{k-1}
            t_contrib = cc * comb(f, 2) * fix_even * (2 ** (k-1))
            # Contribution to twin_SL: 2C(f,2) × Fix_even
            twin_contrib = cc * 2 * comb(f, 2) * fix_even if f >= 2 else 0
            # Difference
            diff = t_contrib - twin_contrib

            if diff != 0:
                factor = comb(f, 2) * fix_even * (2**(k-1) - 2)
                total_diff += cc * factor
                print(f"    {str(ct_list):>15}: cc={cc}, diff_per_sigma = "
                      f"C({f},2) × {fix_even} × ({2**(k-1)} - 2) = "
                      f"{comb(f,2)} × {fix_even} × {2**(k-1)-2} = {factor}")

        T_minus_twin = total_diff // nf
        E_known = known_E.get(n, '?')
        R_val = total_diff // nf - 2 * known_E.get(n, 0)
        print(f"    Total: T-twin = {T_minus_twin}, 2E = {2*known_E.get(n,0)}, R = {R_val}")

    print("\nDONE.")
