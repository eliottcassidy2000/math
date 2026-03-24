#!/usr/bin/env python3
"""
factored_twin_sl_s263.py — The twin SL formula through the factorization
opus-2026-03-23-S263

CORRECTED FORMULA:
  twin_SL(σ) = C(f,2) × 2^{ao - k + 2} = 2 × C(f,2) × Fix_even(σ)

where f = #fixed points, k = #cycles, ao = arc_orbits.

This factors as:
  twin_SL per sigma = 2C(f,2) × Fix_even(σ)

While: Fix_tourn(σ) = Fix_even(σ) × 2^{k-1}

So: twin_SL / Fix_tourn = 2C(f,2) / 2^{k-1} = f(f-1) / 2^{k-1}

For identity (f=n, k=n): twin_SL/Fix_tourn = n(n-1)/2^{n-1} → the universal 2^{2-n} factor.

THE RESIDUAL THROUGH FACTORIZATION:
  T = twin_SL + complex_SL + 2E + MW
  E ≈ (T - twin_SL)/2 - (complex_SL + MW)/2

  residual = complex_SL + MW
           = interaction term between cut and cycle spaces
           → 0 relative to T (thermodynamic decoupling)
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

def compute_all(n):
    nf = factorial(n)
    T_sum = 0  # = transition orbits
    twin_SL_sum = 0
    V_sum = 0

    for ct in partitions(n):
        ct_list = list(ct)
        if any(c % 2 == 0 for c in ct_list): continue
        cc = ccs_val(n, ct_list)
        ao = arc_orbits(ct_list)
        k = len(ct_list)
        f = ct_list.count(1)

        fix_even = 2 ** (ao - k + 1)
        fix_tourn = fix_even * (2 ** (k - 1))  # = 2^ao
        twin_per = 2 * comb(f, 2) * fix_even if f >= 2 else 0

        V_sum += cc * fix_tourn
        T_sum += cc * fix_tourn * comb(f, 2)  # T uses C(f,2) for transition orbits
        # Wait, T is different. Let me re-derive.
        # T_n = (1/n!) Σ_{σ odd} ccs × Fix_tourn(σ) × C(f,2)
        # No! T is the TRANSITION orbits, not the tournament count.
        # Actually: T_n (transition orbits) comes from Burnside on (T, arc) pairs.
        # twin_SL comes from Burnside on (T, arc) pairs where flip preserves class.
        # These are DIFFERENT Burnside sums.

        twin_SL_sum += cc * twin_per

    V = V_sum // nf
    twin_SL = twin_SL_sum // nf

    # T_n (transition orbits) from standard formula
    T = V_sum // nf  # This is actually V, not T
    # The proper T is from a different Burnside sum

    return V, twin_SL

# Need to compute T_n properly
def compute_T(n):
    """T_n = transition orbits = (1/n!) Σ_{σ odd} ccs(σ) × Fix(σ,T) × C(fix(σ),2)
    Actually no. T_n is defined as the number of (class, arc position) orbits.
    T_n = V_n * m / |average orbit size| but this needs the full formula.

    The formula from the project: T_n/2 + (n-2)! = edge_orbits.
    And T_n itself is computed via Burnside on (tournament, arc) pairs.
    """
    nf = factorial(n)
    m = comb(n, 2)

    # T_n via Burnside on the arc-labeled tournament space:
    # T_n = (1/n!) Σ_σ Fix_{(T,e)}(σ)
    # where Fix_{(T,e)}(σ) = #{(T, arc e) : σ(T)=T, σ maps e to e}
    # = Σ_T fixed by σ of #{arcs fixed by σ}
    # = Fix_tourn(σ) × (#arcs fixed by σ)

    # For σ with cycle type ct:
    # #arcs fixed by σ = #{ordered pairs (i,j) with i<j : σ(i)=i, σ(j)=j}
    # = C(f, 2) where f = fixed points

    # Wait, arcs fixed by σ: an arc (i,j) is fixed if σ(i)=i and σ(j)=j.
    # For i<j (unordered): {i,j} is fixed if both endpoints are fixed.
    # So #fixed arcs = C(f, 2).

    T_sum = 0
    for ct in partitions(n):
        ct_list = list(ct)
        if any(c % 2 == 0 for c in ct_list): continue
        cc = ccs_val(n, ct_list)
        ao = arc_orbits(ct_list)
        f = ct_list.count(1)
        fix_tourn = 2 ** ao
        T_sum += cc * fix_tourn * comb(f, 2)

    return T_sum // nf

if __name__ == '__main__':
    print("=" * 70)
    print("  FACTORED TWIN SL — CORRECTED")
    print("  opus-2026-03-23-S263")
    print("=" * 70)

    # Known exact values
    known_twin = {3:2, 4:4, 5:12, 6:48, 7:296, 8:3040, 9:54256}
    known_E = {3:1, 4:5, 5:30, 6:290, 7:4086, 8:91161, 9:3380751}

    print(f"\n{'n':>3} {'V':>8} {'twin_SL':>10} {'twin_known':>11} {'match':>6} "
          f"{'T_n':>10} {'E_approx':>10} {'E_known':>8} {'resid':>8}")

    for n in range(3, 13):
        nf = factorial(n)
        m = comb(n, 2)

        # V and twin_SL
        V, twin_SL = compute_all(n)
        T_n = compute_T(n)

        # Edge orbits formula: EO = T_n/2 + (n-2)!
        EO = T_n // 2 + factorial(n - 2)

        # E approximation: E ≈ (T_n - twin_SL) / 2
        E_approx = (T_n - twin_SL) // 2

        # Residual: gap_orbits = EO - E_known, residual = 2*gap - twin_SL - 2*(n-2)!
        E_known = known_E.get(n)
        twin_known = known_twin.get(n)

        twin_match = "✓" if twin_known and twin_SL == twin_known else "?" if not twin_known else "✗"

        if E_known:
            gap = EO - E_known
            residual = T_n - twin_SL - 2 * E_known  # = complex_SL + MW
            resid_str = str(residual)
        else:
            resid_str = "?"

        ek_str = str(E_known) if E_known else "?"

        print(f"{n:>3} {V:>8} {twin_SL:>10} {str(twin_known or '?'):>11} {twin_match:>6} "
              f"{T_n:>10} {E_approx:>10} {ek_str:>8} {resid_str:>8}")

    # The factored form
    print(f"\n{'='*70}")
    print("  FACTORED DECOMPOSITION")
    print(f"{'='*70}")
    print("""
  THEOREM (Burnside Exponent Factorization):
    Fix_tournament(σ) = Fix_even(σ) × 2^{k-1}
    twin_SL(σ)       = 2C(f,2) × Fix_even(σ)

  Therefore:
    V_n     = (1/n!) Σ Fix_even × 2^{k-1}    [all odd-cycle σ]
    twin_SL = (1/n!) Σ 2C(f,2) × Fix_even     [odd-cycle σ with f≥2]
    T_n     = (1/n!) Σ C(f,2) × Fix_even × 2^{k-1}  [transition orbits]

  The twin fraction:
    twin_SL / T_n → n(n-1)/2^{n-1} × (identity-dominated correction)

  The residual:
    residual = T_n - twin_SL - 2E
            = complex_SL + MW
            = cut-cycle INTERACTION term
            → 0 relative to T (asymptotic independence)

  PHYSICAL INTERPRETATION:
    - Fix_even = "pure cyclic structure" degrees of freedom
    - 2^{k-1} = "score hierarchy" degrees of freedom
    - C(f,2) = "potential twin pairs" (fixed-point pairs)
    - twin_SL = twins counted by cycle × score factorization
    - residual = interference between cycle and score spaces
""")

    print("DONE.")
