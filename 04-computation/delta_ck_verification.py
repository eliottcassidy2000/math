#!/usr/bin/env python3
"""
delta_ck_verification.py -- kind-pasteur-2026-03-13-S60

Verify that the Delta(c_k) between the two E=118 groups at p=13 is
EXACTLY predicted by the eigenvalue expansion:

  c_k = (1/k)[m^k + 2*sum_{t=1}^m Re(z_t^k)]
  Re(z^k) = sum_{r=0}^{floor(k/2)} C(k,2r)*(-1/2)^{k-2r}*(-1)^r * D^{2r}

So Delta(c_k) = (2/k) * sum_r coeff_r * Delta(S_{2r})

where coeff_r = C(k,2r)*(-1/2)^{k-2r}*(-1)^r and Delta(S_{2r}) is the
moment difference between the two groups.

At p=13, E=118 split: Delta(S8)=156, Delta(S10)=2925, Delta(S12)=36884.25
"""

from fractions import Fraction
from math import comb


def eigenvalue_coeff(k, r):
    """Coefficient of D^{2r} in Re(z^k)."""
    return Fraction(comb(k, 2*r)) * Fraction(-1, 2)**(k - 2*r) * (-1)**r


def predict_delta_ck(k, deltas):
    """Predict Delta(c_k) from moment deltas.
    deltas = {2r: Delta(S_{2r})} for all relevant r.
    """
    total = Fraction(0)
    for r in range(1, k // 2 + 1):
        coeff = eigenvalue_coeff(k, r)
        sr = 2 * r
        if sr in deltas:
            total += coeff * Fraction(deltas[sr]).limit_denominator(1000000)
    return Fraction(2, k) * total


def main():
    print("=" * 70)
    print("DELTA(c_k) VERIFICATION via eigenvalue expansion")
    print("=" * 70)

    # Known deltas from p=13 E=118 split
    # From the output: S8 diff = 156.000, S10 diff = 2925.000, S12 diff = 36884.250
    deltas = {
        8: Fraction(156),
        10: Fraction(2925),
        12: Fraction(36884250, 1000)  # 36884.250
    }

    # Known actual deltas from computation
    actual = {
        3: 0,     # c3 same
        5: 0,     # c5 same
        7: 0,     # c7 same
        9: -156,
        11: 2340,
        13: -21762,
    }

    print(f"\n  Known moment deltas (H_hi - H_lo):")
    for sr, dv in sorted(deltas.items()):
        print(f"    Delta(S{sr}) = {float(dv):.4f}")

    print(f"\n  {'k':>3} {'Predicted':>14} {'Actual':>10} {'Match':>6}")
    print(f"  {'-'*36}")

    for k in range(3, 14, 2):
        pred = predict_delta_ck(k, deltas)
        act = actual.get(k, '?')
        match = abs(float(pred) - act) < 0.01 if isinstance(act, (int, float)) else '?'
        print(f"  {k:3d} {float(pred):14.4f} {act:>10} {str(match):>6}")

        # Show the coefficient breakdown
        terms = []
        for r in range(1, k // 2 + 1):
            sr = 2 * r
            if sr in deltas:
                coeff = eigenvalue_coeff(k, r)
                frac_coeff = Fraction(2, k) * coeff
                val = float(frac_coeff * deltas[sr])
                terms.append(f"({float(frac_coeff):.6f})*Delta(S{sr})={val:.2f}")
        if terms:
            print(f"      = {' + '.join(terms)}")

    # Why c3, c5, c7 have zero delta:
    print(f"\n  --- Why c3, c5, c7 have Delta = 0 ---")
    for k in [3, 5, 7]:
        max_D_power = 2 * (k // 2)  # maximum even power of D in Re(z^k)
        print(f"    c_{k}: Re(z^{k}) has D^{{2r}} for r=0,...,{k//2}, "
              f"max D power = {max_D_power}")
        if max_D_power < 8:
            print(f"      Max D power {max_D_power} < 8, so no S8 (or higher) dependence")
            print(f"      Since Delta(S4)=0, Delta(S6)=0: Delta(c_{k})=0")

    # S10 cancellation at p=11
    print(f"\n  --- S10 CANCELLATION at p=11 ---")
    print(f"  At p=11, m=5, H = f(S4,S6,S8) exactly (S10 not needed)")
    print(f"  This means: sum_{{k odd}} (2/k) * coeff_S10_in_Re(z^k) = 0")
    print(f"  when summed with appropriate alpha_j weights")

    # Compute the S10 coefficient in alpha_1 = sum c_k
    print(f"\n  S_{10} coefficient in alpha_1 = sum c_k for p=11:")
    total_s10_coeff = Fraction(0)
    for k in range(3, 12, 2):
        r = 5  # S10 = S_{2*5}
        if 2*r <= k:
            coeff = eigenvalue_coeff(k, r)
            total_s10_coeff += Fraction(2, k) * coeff
            print(f"    c_{k}: (2/{k})*C({k},{2*r})*(-1/2)^{k-2*r}*(-1)^{r} "
                  f"= {float(Fraction(2,k)*coeff):.8f}")
    print(f"    Total S10 coeff in alpha_1 = {float(total_s10_coeff):.8f} = {total_s10_coeff}")

    # At p=13, check if any cancellation occurs
    print(f"\n  S_{12} coefficient in alpha_1 = sum c_k for p=13:")
    total_s12_coeff = Fraction(0)
    for k in range(3, 14, 2):
        r = 6  # S12 = S_{2*6}
        if 2*r <= k:
            coeff = eigenvalue_coeff(k, r)
            total_s12_coeff += Fraction(2, k) * coeff
            print(f"    c_{k}: (2/{k})*C({k},{2*r})*(-1/2)^{k-2*r}*(-1)^{r} "
                  f"= {float(Fraction(2,k)*coeff):.8f}")
    print(f"    Total S12 coeff in alpha_1 = {float(total_s12_coeff):.8f} = {total_s12_coeff}")

    # Generalise: compute S_{p-1} coeff in alpha_1 for various primes
    print(f"\n  --- S_{{p-1}} coefficient in alpha_1 for various primes ---")
    for p in [7, 11, 13, 17, 19, 23]:
        m = (p - 1) // 2
        r = m  # S_{2m} = S_{p-1}
        total = Fraction(0)
        for k in range(3, p + 1, 2):
            if 2*r <= k:
                coeff = eigenvalue_coeff(k, r)
                total += Fraction(2, k) * coeff
        print(f"    p={p:>2} (mod 4 = {p%4}): S_{p-1} coeff in alpha_1 = "
              f"{float(total):.8f} = {total}")

    # Check: is this coefficient always 0 for p=3 mod 4?
    print(f"\n  PATTERN: Does S_{{p-1}} coeff vanish for p = 3 mod 4?")


if __name__ == '__main__':
    main()
    print("\nDONE.")
