#!/usr/bin/env python3
"""
copart_consecutive_qr.py -- Connect co-participation theorem to classical
number theory of consecutive quadratic residues.

The Chord Co-Participation Theorem (HYP-601) says:
  copart > 0 iff:
    chi(q)=+1 AND chi(q+1)=-1, OR
    chi(q)=-1 AND chi(q-1)=+1.

This is equivalent to: q is adjacent to a QR/QNR TRANSITION.
Specifically: (q, q+1) or (q-1, q) is a (QR, QNR) or (QNR, QR) pair.

CLASSICAL RESULT (Davenport 1931, Perron): The number of consecutive pairs
(n, n+1) with chi(n)=+1, chi(n+1)=-1 is asymptotically (p-1)/4.

QUESTION: What fraction of odd q values in [3, p-2] give copart > 0?
Is there a closed form? Does it relate to known number theory?

Author: kind-pasteur-2026-03-12-S60
"""

import math


def legendre(a, p):
    if a % p == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def copart_fraction(p):
    """Fraction of odd q in [3, p-2] that give copart > 0."""
    m = (p - 1) // 2
    odd_q = [q for q in range(3, p, 2)]
    copart_count = 0
    for q in odd_q:
        chi_q = legendre(q, p)
        if chi_q == 1:
            copart = (legendre(q + 1, p) == -1)
        else:
            copart = (legendre(q - 1, p) == 1)
        if copart:
            copart_count += 1
    return copart_count, len(odd_q), copart_count / len(odd_q) if odd_q else 0


def consecutive_qr_analysis(p):
    """Count all consecutive (QR, QNR) and (QNR, QR) transitions."""
    m = (p - 1) // 2

    # Type A: chi(n)=+1, chi(n+1)=-1 (QR -> QNR)
    # Type B: chi(n)=-1, chi(n+1)=+1 (QNR -> QR)
    # Type C: chi(n)=+1, chi(n+1)=+1 (QR -> QR)
    # Type D: chi(n)=-1, chi(n+1)=-1 (QNR -> QNR)

    counts = {'A': 0, 'B': 0, 'C': 0, 'D': 0}
    for n in range(1, p - 1):
        chi_n = legendre(n, p)
        chi_n1 = legendre(n + 1, p)
        if chi_n == 1 and chi_n1 == -1:
            counts['A'] += 1
        elif chi_n == -1 and chi_n1 == 1:
            counts['B'] += 1
        elif chi_n == 1 and chi_n1 == 1:
            counts['C'] += 1
        else:
            counts['D'] += 1
    return counts


def main():
    print("=" * 70)
    print("COPART vs CONSECUTIVE QR ANALYSIS")
    print("=" * 70)

    # PART 1: Copart fractions at various primes
    print("\n--- PART 1: COPART FRACTION FOR ODD q ---")
    print(f"{'p':>5} {'p mod 4':>7} {'copart_Y':>9} {'total_q':>8} {'fraction':>10} {'(p-3)/4':>8}")

    for p in [7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83]:
        if p % 4 != 3:
            continue
        cy, total, frac = copart_fraction(p)
        print(f"{p:>5} {p%4:>7} {cy:>9} {total:>8} {frac:>10.4f} {(p-3)/4:>8.1f}")

    # PART 2: Consecutive QR pair analysis
    print("\n\n--- PART 2: CONSECUTIVE QR PAIR TYPES ---")
    print(f"{'p':>5} {'A(QR->QNR)':>12} {'B(QNR->QR)':>12} {'C(QR->QR)':>12} {'D(QNR->QNR)':>12} {'A+B':>6} {'C+D':>6}")

    for p in [7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83]:
        if p % 4 != 3:
            continue
        counts = consecutive_qr_analysis(p)
        print(f"{p:>5} {counts['A']:>12} {counts['B']:>12} {counts['C']:>12} {counts['D']:>12} "
              f"{counts['A']+counts['B']:>6} {counts['C']+counts['D']:>6}")

    # PART 3: Connection to copart
    # Copart > 0 when:
    #   chi(q)=+1 AND chi(q+1)=-1: this is a TYPE A transition at n=q (odd)
    #   chi(q)=-1 AND chi(q-1)=+1: this is a TYPE A transition at n=q-1 (even)
    # So copart counts = (# Type A transitions at odd n) + (# Type A transitions at even n, shifted)
    # Actually: copart at odd q counts transitions NEAR q, not all transitions.

    print("\n\n--- PART 3: DETAILED TRANSITION ANALYSIS ---")
    print("Copart > 0 iff (chi(q)=+1, chi(q+1)=-1) or (chi(q)=-1, chi(q-1)=+1)")
    print("First condition = Type A transition at ODD q")
    print("Second condition = Type B transition at ODD q (since chi(q-1)=+1, chi(q)=-1)")
    print("NO! Second is: chi(q)=-1 and chi(q-1)=+1, i.e., Type A at EVEN n=q-1")
    print()

    for p in [7, 11, 19, 23, 31]:
        if p % 4 != 3:
            continue
        m = (p - 1) // 2
        print(f"  p={p}:")

        # Count by detailed mechanism
        type_a_odd = 0   # chi(q)=+1, chi(q+1)=-1, q odd
        type_rev_even = 0  # chi(q)=-1, chi(q-1)=+1, q odd => A at even n=q-1

        for q in range(3, p, 2):
            chi_q = legendre(q, p)
            if chi_q == 1:
                if legendre(q + 1, p) == -1:
                    type_a_odd += 1
                    print(f"    q={q:>2}: chi(q)=+1, chi(q+1)=-1 -> copart=Y [A-type at q]")
                else:
                    print(f"    q={q:>2}: chi(q)=+1, chi(q+1)=+1 -> copart=N [C-type at q]")
            else:
                if legendre(q - 1, p) == 1:
                    type_rev_even += 1
                    print(f"    q={q:>2}: chi(q)=-1, chi(q-1)=+1 -> copart=Y [A-type at q-1]")
                else:
                    print(f"    q={q:>2}: chi(q)=-1, chi(q-1)=-1 -> copart=N [D-type at q-1]")

        print(f"    Total: {type_a_odd} from A-at-odd + {type_rev_even} from A-at-even = {type_a_odd + type_rev_even}")
        print()

    # PART 4: The q^2-1 connection (from copart_resonance_proof.py)
    print("\n--- PART 4: chi(q^2-1) CONNECTION ---")
    print("From algebra: copart > 0 iff chi((a+b)(b-a)) = chi(b^2-a^2) = chi(q^2-1) = ...")
    print("For q-resonance (qa=b): b^2-a^2 = a^2(q^2-1)")
    print("chi(b^2-a^2) = chi(a^2)*chi(q^2-1) = chi(q^2-1)")
    print()
    print("So: copart > 0 iff chi(q^2-1) = +1")
    print("    q^2-1 = (q-1)(q+1)")
    print("    chi(q^2-1) = chi(q-1)*chi(q+1)")
    print()
    print("CLAIM: copart > 0 iff chi(q-1)*chi(q+1) = +1")
    print("Let's verify this is equivalent to HYP-601:")
    print()

    for p in [7, 11, 19, 23, 31, 43]:
        if p % 4 != 3:
            continue
        m = (p - 1) // 2
        match = 0
        total = 0
        for q in range(3, p, 2):
            chi_q = legendre(q, p)
            chi_q1 = legendre(q + 1, p)
            chi_qm1 = legendre(q - 1, p)

            # HYP-601 prediction
            if chi_q == 1:
                hyp601 = (chi_q1 == -1)
            else:
                hyp601 = (chi_qm1 == 1)

            # chi(q^2-1) prediction
            chi_q2m1 = legendre(q * q - 1, p)
            q2m1_pred = (chi_q2m1 == 1)

            # Alternative: chi(q-1)*chi(q+1)
            product_pred = (chi_qm1 * chi_q1 == 1)

            total += 1
            if hyp601 == q2m1_pred:
                match += 1

        print(f"  p={p}: HYP-601 vs chi(q^2-1): {match}/{total} match")

    # PART 5: The PRODUCT of chi(q-1)*chi(q+1)
    print("\n\n--- PART 5: VERIFY chi(q-1)*chi(q+1) = chi(q^2-1) ---")
    print("This is just multiplicativity of chi. Always true.")
    print()
    print("So the copart theorem simplifies to:")
    print("  copart(q) > 0 iff chi(q^2-1) = +1")
    print("  i.e., (q-1)(q+1) is a QR mod p")
    print()
    print("But wait -- is this EQUIVALENT to HYP-601?")
    print("HYP-601: chi(q)=+1 => copart iff chi(q+1)=-1")
    print("         chi(q)=-1 => copart iff chi(q-1)=+1")
    print()
    print("Case 1: chi(q)=+1")
    print("  chi(q-1)*chi(q+1)=+1 means chi(q-1)=chi(q+1)")
    print("  q is odd, q-1 is even. For copart: chi(q+1)=-1")
    print("  If chi(q+1)=-1, then chi(q-1)=-1, product=+1. OK.")
    print("  If chi(q+1)=+1, then chi(q-1)=+1, product=+1. PROBLEM!")
    print("  Wait, product=+1 in BOTH cases if chi(q-1)=chi(q+1)!")
    print()
    print("So chi(q^2-1)=+1 does NOT reduce to HYP-601 directly.")
    print("Let me recheck...")

    # Actually verify if they match
    print("\n\n--- PART 6: CAREFUL COMPARISON ---")
    for p in [7, 11, 19, 23, 31, 43]:
        if p % 4 != 3:
            continue
        mismatches = []
        for q in range(3, p, 2):
            chi_q = legendre(q, p)
            chi_q1 = legendre(q + 1, p)
            chi_qm1 = legendre(q - 1, p)

            if chi_q == 1:
                hyp601 = (chi_q1 == -1)
            else:
                hyp601 = (chi_qm1 == 1)

            chi_q2m1 = chi_qm1 * chi_q1  # = chi(q^2-1)
            q2m1_pred = (chi_q2m1 == 1)

            if hyp601 != q2m1_pred:
                mismatches.append((q, chi_q, chi_qm1, chi_q1, hyp601, q2m1_pred))

        if mismatches:
            print(f"  p={p}: {len(mismatches)} MISMATCHES:")
            for q, cq, cqm1, cq1, h, qp in mismatches:
                print(f"    q={q}: chi(q)={cq:+d}, chi(q-1)={cqm1:+d}, chi(q+1)={cq1:+d}, "
                      f"HYP-601={h}, chi(q^2-1)={qp}")
        else:
            print(f"  p={p}: PERFECT MATCH")

    # PART 7: Connection to sign law
    print("\n\n--- PART 7: SIGN LAW CONNECTION ---")
    print("The sign formula from resonance cascade:")
    print("  sign(D^{2n}_onset) * chi(ab) = chi(q)")
    print()
    print("For the product law sign(h_hat) = chi(ab):")
    print("  Need sign(h_hat) = chi(ab)")
    print("  If h_hat dominated by D^{2n}_onset:")
    print("    sign(h_hat) = sign(D^{2n}_onset) = chi(q)*chi(ab)")
    print("    Product law requires chi(q)=+1")
    print()
    print("At p=3 mod 4 with chi(q)=-1 (e.g., q=3 at p=19):")
    print("  Dominant D^4 has wrong sign")
    print("  ALSO: copart(q=3) depends on chi(8)=chi(2)^3")
    print()

    for p in [7, 11, 19, 23, 31, 43]:
        if p % 4 != 3:
            continue
        chi2 = legendre(2, p)
        chi3 = legendre(3, p)
        # q=3: copart iff chi(3)=-1 => chi(3-1)=chi(2)=+1
        #       or chi(3)=+1 => chi(3+1)=chi(4)=chi(2)^2=+1... need =-1, never true
        if chi3 == -1:
            copart_q3 = (chi2 == 1)
            mech = f"chi(3)=-1, chi(2)={chi2:+d}, copart={'Y' if copart_q3 else 'N'}"
        else:
            chi4 = legendre(4, p)
            copart_q3 = (chi4 == -1)
            mech = f"chi(3)=+1, chi(4)={chi4:+d}, copart={'Y' if copart_q3 else 'N'}"

        # q=5
        chi5 = legendre(5, p)
        if chi5 == 1:
            copart_q5 = (legendre(6, p) == -1)
            mech5 = f"chi(5)=+1, chi(6)={legendre(6,p):+d}, copart={'Y' if copart_q5 else 'N'}"
        else:
            copart_q5 = (legendre(4, p) == 1)
            mech5 = f"chi(5)=-1, chi(4)={legendre(4,p):+d}, copart={'Y' if copart_q5 else 'N'}"

        print(f"  p={p}: chi(2)={chi2:+d}, chi(3)={chi3:+d}")
        print(f"    q=3: {mech}")
        print(f"    q=5: {mech5}")
        print(f"    Sign law OK for q=3? chi(3)={'=+1 YES' if chi3==1 else '=-1 NO (D^4 wrong sign)'}")
        print()

    # PART 8: The fraction of copart=Y converges to what?
    print("\n--- PART 8: ASYMPTOTIC COPART FRACTION ---")
    print("Heuristically, chi(q) and chi(q+/-1) are 'independent' for large p.")
    print("Then P(copart=Y) = P(chi(q)=+1)*P(chi(q+1)=-1) + P(chi(q)=-1)*P(chi(q-1)=+1)")
    print("                  = (1/2)(1/2) + (1/2)(1/2) = 1/2")
    print()
    print(f"{'p':>5} {'fraction':>10} {'1/2':>6}")
    for p in [7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83, 103, 107, 127, 131, 151, 163, 167, 179, 191, 199]:
        if p % 4 != 3:
            continue
        cy, total, frac = copart_fraction(p)
        print(f"{p:>5} {frac:>10.4f} {'0.5':>6}")


if __name__ == '__main__':
    main()
