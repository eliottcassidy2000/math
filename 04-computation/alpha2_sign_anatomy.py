#!/usr/bin/env python3
"""
alpha2_sign_anatomy.py -- Understand the alpha_2 sign anomaly at p=11

alpha_1 and alpha_3 follow the product law: sign(h_hat_alpha_j[{a,b}]) = chi(ab)
alpha_2 does NOT: it has two magnitude classes with mixed signs.

The key question: what determines the alpha_2 sign?

At p=7 (from alpha2_walsh_anatomy.py):
  h_hat_alpha_2[{0,1}] = -1.75, chi(1*2)=chi(2)=+1 -> anti-product-law
  h_hat_alpha_2[{0,2}] = -1.75, chi(1*3)=chi(3)=-1 -> product-law
  h_hat_alpha_2[{1,2}] = -1.75, chi(2*3)=chi(6)=-1 -> product-law

  So at p=7: alpha_2 has UNIFORM sign (-) and UNIFORM magnitude (1.75)
  Product law holds for 2/3 pairs, anti for 1/3.
  But actually ALL have sign -, and chi(ab) values are +,-,-, so:
  sign(-) = -chi(ab)? For pair (0,1): chi=+1, sign=-, so sign=-chi(ab)
  For pair (0,2): chi=-1, sign=-, so sign=chi(ab)
  For pair (1,2): chi=-1, sign=-, so sign=chi(ab)

  So at p=7 alpha_2 has constant NEGATIVE sign, which agrees with product-law
  for QR pairs and anti-product-law for NQR pairs.

At p=11: TWO magnitude classes with different chi-agreement patterns.

Hypothesis: the alpha_2 sign depends on the resonance level q of the pair.

Author: kind-pasteur-2026-03-12-S60
"""


def legendre(a, p):
    if a % p == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def resonance_level(a, b, p):
    """Minimum odd q such that qa = +/-b or a = +/-qb mod p."""
    for q in range(1, p, 2):
        if (q * a - b) % p == 0 or (q * a + b) % p == 0:
            return q
        if q > 1 and ((a - q * b) % p == 0 or (a + q * b) % p == 0):
            return q
    return p


def main():
    for p in [7, 11]:
        m = (p - 1) // 2
        print("=" * 60)
        print(f"ALPHA_2 SIGN ANATOMY at p={p}")
        print("=" * 60)

        # alpha_2 Walsh coefficients (from previous computations)
        if p == 7:
            # All have h_hat = -1.75
            alpha2_walsh = {
                (0, 1): -1.75,
                (0, 2): -1.75,
                (1, 2): -1.75,
            }
        elif p == 11:
            alpha2_walsh = {
                (0, 1): 47.4375,
                (0, 2): 29.5625,
                (0, 3): 29.5625,
                (0, 4): -47.4375,
                (1, 2): -29.5625,
                (1, 3): 47.4375,
                (1, 4): -29.5625,
                (2, 3): -47.4375,
                (2, 4): -47.4375,
                (3, 4): 29.5625,
            }

        print(f"\n  {'Pair':>8s} {'gaps':>8s} {'q':>3s} {'chi(ab)':>8s} "
              f"{'h_hat':>12s} {'sign':>5s} {'pattern':>12s}")

        for (a, b), h in sorted(alpha2_walsh.items()):
            gap_a, gap_b = a + 1, b + 1
            chi_ab = legendre(gap_a * gap_b, p)
            q = resonance_level(gap_a, gap_b, p)
            sign = '+' if h > 0 else '-'

            # Check various patterns
            product_law = (h > 0) == (chi_ab > 0)
            anti_product = (h > 0) != (chi_ab > 0)

            # What about chi(q) * chi(ab)?
            chi_q = legendre(q, p)

            pattern = "product" if product_law else "anti-product"

            print(f"  ({a},{b})  ({gap_a},{gap_b})  {q:>3d}  {chi_ab:>+8d}  "
                  f"{h:>12.4f}  {sign:>5s}  {pattern:>12s}  "
                  f"chi(q)={chi_q:+d}")

        # Test: is sign(h_hat_alpha2) = -chi(q) * chi(ab)?
        print(f"\n  Testing: sign(h_hat_alpha2) = -chi(q) * chi(ab)?")
        for (a, b), h in sorted(alpha2_walsh.items()):
            gap_a, gap_b = a + 1, b + 1
            chi_ab = legendre(gap_a * gap_b, p)
            q = resonance_level(gap_a, gap_b, p)
            chi_q = legendre(q, p)
            predicted = -chi_q * chi_ab
            actual = 1 if h > 0 else -1
            match = predicted == actual
            print(f"    ({a},{b}): predicted={predicted:+d}, actual={actual:+d}, "
                  f"{'OK' if match else 'FAIL'}")

        # Test: is sign(h_hat_alpha2) = chi(ab) * (-1)^{something}?
        # More generally, classify by q
        print(f"\n  Grouped by resonance level q:")
        from collections import defaultdict
        by_q = defaultdict(list)
        for (a, b), h in sorted(alpha2_walsh.items()):
            gap_a, gap_b = a + 1, b + 1
            q = resonance_level(gap_a, gap_b, p)
            by_q[q].append((a, b, h))

        for q in sorted(by_q.keys()):
            entries = by_q[q]
            magnitudes = [abs(h) for _, _, h in entries]
            signs_vs_chi = []
            for a, b, h in entries:
                gap_a, gap_b = a + 1, b + 1
                chi_ab = legendre(gap_a * gap_b, p)
                s = 1 if h > 0 else -1
                signs_vs_chi.append(s * chi_ab)  # +1 if same, -1 if opposite

            print(f"    q={q}: |h_hat|={set(magnitudes)}, "
                  f"sign*chi(ab)={set(signs_vs_chi)}")

        # H comparison: how does the cancellation work per q-class?
        print(f"\n  OCF cancellation per resonance class:")

        if p == 7:
            a1_walsh = {(0,1): 5.25, (0,2): 5.25, (1,2): 5.25}
        elif p == 11:
            a1_walsh = {
                (0,1): -35.75, (0,2): 13.75, (0,3): 13.75, (0,4): 35.75,
                (1,2): -13.75, (1,3): -35.75, (1,4): -13.75,
                (2,3): 35.75, (2,4): 35.75, (3,4): 13.75,
            }
            a3_walsh = {
                (0,1): -15.8125, (0,2): 15.8125, (0,3): 15.8125, (0,4): 15.8125,
                (1,2): -15.8125, (1,3): -15.8125, (1,4): -15.8125,
                (2,3): 15.8125, (2,4): 15.8125, (3,4): 15.8125,
            }

        for q in sorted(by_q.keys()):
            entries = by_q[q]
            a, b, h_a2 = entries[0]
            gap_a, gap_b = a + 1, b + 1
            chi_ab = legendre(gap_a * gap_b, p)

            h_a1 = a1_walsh.get((a, b), 0)
            if p == 11:
                h_a3 = a3_walsh.get((a, b), 0)
            else:
                h_a3 = 0

            h_H = 2*h_a1 + 4*h_a2 + (8*h_a3 if p >= 11 else 0)

            print(f"\n    q={q}, pair ({a},{b}), chi(ab)={chi_ab:+d}:")
            print(f"      2*h_a1  = {2*h_a1:>12.4f} (sign vs chi: {'same' if (h_a1>0)==(chi_ab>0) else 'anti'})")
            print(f"      4*h_a2  = {4*h_a2:>12.4f} (sign vs chi: {'same' if (h_a2>0)==(chi_ab>0) else 'anti'})")
            if p >= 11:
                print(f"      8*h_a3  = {8*h_a3:>12.4f} (sign vs chi: {'same' if (h_a3>0)==(chi_ab>0) else 'anti'})")
            print(f"      h_hat_H = {h_H:>12.4f} (sign vs chi: {'same' if (h_H>0)==(chi_ab>0) else 'anti'})")

    print("\nDONE.")


if __name__ == '__main__':
    main()
