#!/usr/bin/env python3
"""
twin_sl_through_factorization_s263.py
opus-2026-03-23-S263

The factorization theorem says:
  Fix_tournament(σ) = Fix_even(σ) × 2^{k(σ)-1}

This separates tournament enumeration into cycle (even graph) and cut (score) parts.

KEY QUESTION: How does twin_SL decompose under this factorization?

Twin SL: twin_SL(n) = (1/n!) Σ_{σ odd} ccs(σ) × C(f(σ),2) × 2^{ao(σ) - extra}

The twin mechanism requires TWO fixed points that are "twins" (identical neighborhoods).
In the factorization:
- The CYCLE part: the two fixed points have the same cycle-space projection neighbors
- The CUT part: the two fixed points have the same score

So twin_SL decomposes into:
  twin_SL = (1/n!) Σ_{σ odd} ccs(σ) × [cycle_twin(σ) × cut_twin(σ)]

Can we compute cycle_twin and cut_twin separately?

Also explore: the edge formula E ≈ (T - twin_SL)/2 through the factorization.
If T = Σ Fix_even × 2^{k-1} and twin_SL involves C(f,2), what's the structure?
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

def edge_orbits(ct):
    total = 0
    for i in range(len(ct)):
        total += ct[i] // 2
        for j in range(i + 1, len(ct)):
            total += gcd(ct[i], ct[j])
    return total

if __name__ == '__main__':
    print("=" * 70)
    print("  TWIN SL THROUGH THE FACTORIZATION LENS")
    print("  opus-2026-03-23-S263")
    print("=" * 70)

    for n in range(3, 12):
        nf = factorial(n)

        # Standard quantities
        T_sum = 0      # total transitions
        twin_sum = 0   # twin SL

        # Factored quantities
        cycle_T = 0    # Σ ccs × Fix_even
        cut_T = 0      # weighted average of 2^{k-1}

        # Per cycle type analysis
        for ct in partitions(n):
            ct_list = list(ct)
            if any(c % 2 == 0 for c in ct_list): continue

            cc = ccs_val(n, ct_list)
            ao = arc_orbits(ct_list)
            k = len(ct_list)
            f = ct_list.count(1)

            fix_tourn = 2 ** ao
            fix_even = 2 ** (ao - k + 1)  # = fix_tourn / 2^{k-1}
            fix_cut = 2 ** (k - 1)

            # Twin SL contribution from this cycle type:
            # Twin requires 2 fixed points → C(f,2) twin pairs
            # Each twin pair: the arc between them is "free" in the orientation
            # For the rest: fix count = 2^{ao - 1} (one fewer free arc)
            # But this is the kind-pasteur formula's inner structure

            # Standard twin_SL:
            # twin_SL(σ) = C(f,2) × 2^{ao - something}
            # From S20dj: twin_SL = C(f,2) × 2^{arc_orbits - #other_cycles}
            # where #other_cycles is the number of non-fixed-point cycles
            nfc = k - f  # non-fixed-point cycles
            twin_per_sigma = comb(f, 2) * (2 ** (ao - nfc)) if f >= 2 else 0

            # But wait: let me reconsider. The formula from kind-pasteur S20dj:
            # twin_SL(n) = (1/n!) Σ_{odd σ} ccs × C(f,2) × 2^{arc_orbits - #other_cycles}
            # "other_cycles" = cycles of length > 1 = nfc

            T_sum += cc * fix_tourn
            twin_sum += cc * twin_per_sigma

        T = T_sum // nf
        twin_SL = twin_sum // nf

        # E approximation
        E_approx = (T - twin_SL) // 2

        # Also compute the factored version
        # T = (1/n!) Σ_{odd σ} ccs × Fix_even × 2^{k-1}
        # twin_SL = (1/n!) Σ_{odd σ, f≥2} ccs × C(f,2) × 2^{ao - nfc}

        # Analyze the twin fraction
        twin_frac = twin_SL / T if T > 0 else 0

        # The twin fraction should be ~ 2^{2-n} = 4/2^n
        predicted_frac = 4.0 / (2 ** n)

        print(f"\n  n={n}: T={T}, twin_SL={twin_SL}, E≈{E_approx}")
        print(f"    twin_frac = {twin_frac:.6f}, predicted = {predicted_frac:.6f}, "
              f"ratio = {twin_frac/predicted_frac:.4f}")

        # Decompose twin_SL contribution by cycle type
        if n <= 8:
            print(f"    {'ct':>15} {'ccs':>6} {'f':>3} {'k':>3} {'nfc':>4} "
                  f"{'C(f,2)':>6} {'2^(ao-nfc)':>12} {'twin_contrib':>14}")
            for ct in partitions(n):
                ct_list = list(ct)
                if any(c % 2 == 0 for c in ct_list): continue
                cc = ccs_val(n, ct_list)
                ao = arc_orbits(ct_list)
                k = len(ct_list)
                f = ct_list.count(1)
                nfc = k - f
                twin_c = comb(f,2) * (2 ** (ao - nfc)) if f >= 2 else 0
                if twin_c > 0:
                    print(f"    {str(ct_list):>15} {cc:>6} {f:>3} {k:>3} {nfc:>4} "
                          f"{comb(f,2):>6} {2**(ao-nfc):>12} {cc*twin_c:>14}")

    # Now: the connection between twin_SL and the factorization
    print(f"\n{'='*70}")
    print("  TWIN_SL IN FACTORED FORM")
    print(f"{'='*70}")

    print("""
  For odd-cycle σ with k cycles, f fixed points, nfc non-fixed cycles:
    k = f + nfc
    ao = Σ (c_i-1)/2 + Σ gcd(c_i,c_j)

  twin_SL contribution = C(f,2) × 2^{ao - nfc}

  Factored: 2^{ao - nfc} = 2^{ao - k + f} = 2^{ao - k + 1} × 2^{f - 1}
                          = Fix_even × 2^{f - 1}

  So: twin_SL(σ) = C(f,2) × Fix_even(σ) × 2^{f - 1}

  The twin mechanism decomposes as:
  - Fix_even(σ): the cycle-space (even graph) part
  - C(f,2) × 2^{f-1}: the cut-space (score) part with twin constraint

  Note: C(f,2) × 2^{f-1} = f(f-1)/2 × 2^{f-1} = f(f-1) × 2^{f-2}

  The "extra" 2^{k-1} in Fix_tourn comes from cut freedom.
  The "twin reduction" from 2^{k-1} to C(f,2) × 2^{f-1} is:
    2^{k-1} / (C(f,2) × 2^{f-1}) = 2^{nfc} / C(f,2)

  This is the RATIO of total cut-space orientations to twin-compatible ones.
""")

    # Verify the factored twin formula
    print("  VERIFICATION: twin_SL(σ) = C(f,2) × Fix_even(σ) × 2^{f-1}")
    for n in range(3, 10):
        nf = factorial(n)
        twin_direct = 0
        twin_factored = 0

        for ct in partitions(n):
            ct_list = list(ct)
            if any(c % 2 == 0 for c in ct_list): continue
            cc = ccs_val(n, ct_list)
            ao = arc_orbits(ct_list)
            k = len(ct_list)
            f = ct_list.count(1)
            nfc = k - f
            fix_even = 2 ** (ao - k + 1)

            # Direct formula
            td = comb(f,2) * (2 ** (ao - nfc)) if f >= 2 else 0
            # Factored formula
            tf = comb(f,2) * fix_even * (2 ** (f - 1)) if f >= 2 else 0

            twin_direct += cc * td
            twin_factored += cc * tf

        td_val = twin_direct // nf
        tf_val = twin_factored // nf
        print(f"  n={n}: direct={td_val}, factored={tf_val}, match={'✓' if td_val == tf_val else '✗'}")

    print("\nDONE.")
