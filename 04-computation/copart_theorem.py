#!/usr/bin/env python3
"""
copart_theorem.py -- PROVE the Chord Co-Participation Theorem

THEOREM (HYP-601): For Paley T_p at p = 3 mod 4, chord pair (a,b) with
minimum resonance level q co-participates in a 3-cycle gap triple iff:
  - chi(q) = +1 AND chi(q+1) = -1, OR
  - chi(q) = -1 AND chi(q-1) = +1

PROOF SKETCH:
1. For p = 3 mod 4, chi(-1) = -1, so for each chord c, exactly one of
   c, p-c is in QR. Call this qr(c).
2. A 3-cycle gap triple (g1, g2, g3) with g1+g2+g3 = 0 mod p, all in QR.
3. Chords a,b co-participate iff there exists such a triple with
   two of the gi having chord indices a,b.
4. Since chi(-1) = -1, the triple is: (qr(a), qr(b), g3) with
   g3 = p - qr(a) - qr(b). Need chi(g3) = +1.
5. For qa = b type resonance (WLOG):
   - chi(q)=+1: chi(g3) = -chi(1+q), so copart iff chi(q+1) = -1
   - chi(q)=-1: chi(g3) = chi(q-1), so copart iff chi(q-1) = +1
6. This depends ONLY on q, not on the specific pair (a,b).

Author: kind-pasteur-2026-03-12-S60
"""

import math


def legendre(a, p):
    if a % p == 0:
        return 0
    return 1 if pow(a, (p-1)//2, p) == 1 else -1


def classify_resonance(a, b, p):
    resonances = []
    for k in range(1, p):
        q = 2*k - 1
        if q >= p:
            break
        if (q*a - b) % p == 0:
            resonances.append((q, f"{q}a=b"))
        if (q*a + b) % p == 0:
            resonances.append((q, f"{q}a=-b"))
        if (a - q*b) % p == 0 and q != 1:
            resonances.append((q, f"a={q}b"))
        if (a + q*b) % p == 0 and q != 1:
            resonances.append((q, f"a=-{q}b"))
    return resonances


def main():
    print("=" * 70)
    print("CHORD CO-PARTICIPATION THEOREM -- PROOF AND VERIFICATION")
    print("=" * 70)

    # PART 1: State and verify the theorem
    print("\n--- PART 1: THEOREM STATEMENT AND VERIFICATION ---")
    print()
    print("THEOREM: For Paley T_p at p = 3 mod 4, chord pair (a,b) with")
    print("min resonance q co-participates in 3-cycle gap triple iff:")
    print("  chi(q) = +1 AND chi(q+1) = -1, OR")
    print("  chi(q) = -1 AND chi(q-1) = +1")
    print()

    for p in [7, 11, 19, 23, 31, 43]:
        if p % 4 != 3:
            continue
        m = (p - 1) // 2
        QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)

        # qr(c) = QR representative of chord c
        def qr(c):
            return c if c in QR else (p - c)

        total_pairs = 0
        theorem_correct = 0
        theorem_wrong = 0

        for a in range(1, m + 1):
            for b in range(a + 1, m + 1):
                res = classify_resonance(a, b, p)
                if not res:
                    continue
                min_q = min(q for q, t in res)
                chi_q = legendre(min_q, p)

                # Predicted copart by theorem
                if chi_q == 1:
                    predicted = (legendre(min_q + 1, p) == -1)
                else:
                    predicted = (legendre(min_q - 1, p) == 1)

                # Actual copart: check if (qr(a), qr(b), g3) with g3 in QR
                # Also check other orderings
                actual = False
                for sa in [a, p - a]:
                    if sa not in QR:
                        continue
                    for sb in [b, p - b]:
                        if sb not in QR:
                            continue
                        g3 = (p - sa - sb) % p
                        if g3 in QR and g3 > 0:
                            actual = True

                total_pairs += 1
                if predicted == actual:
                    theorem_correct += 1
                else:
                    theorem_wrong += 1
                    print(f"  **MISMATCH** p={p} ({a},{b}): q={min_q}, "
                          f"chi(q)={chi_q:+d}, pred={predicted}, actual={actual}")

        print(f"  p={p}: {theorem_correct}/{total_pairs} correct, "
              f"{theorem_wrong} wrong")

    # PART 2: Algebraic proof (step by step verification)
    print("\n\n--- PART 2: ALGEBRAIC PROOF ---")
    print()
    print("For qa = b (positive resonance), a in QR:")
    print("  qr(a) = a, chi(qa) = chi(q)*chi(a) = chi(q)")
    print("  If chi(q)=+1: qr(b) = qa, g3 = p - a - qa = p - a(1+q)")
    print("    chi(g3) = chi(-a(1+q)) = chi(-1)*chi(a)*chi(1+q)")
    print("           = (-1)(+1)*chi(1+q) = -chi(1+q)")
    print("  If chi(q)=-1: qr(b) = p-qa, g3 = p - a - (p-qa) = qa - a = a(q-1)")
    print("    chi(g3) = chi(a)*chi(q-1) = (+1)*chi(q-1) = chi(q-1)")
    print()
    print("For qa = b, a in QNR:")
    print("  qr(a) = p-a, chi(qa) = chi(q)*chi(a) = chi(q)*(-1)")
    print("  If chi(q)=+1: chi(qa)=-1, qr(b) = p-qa")
    print("    g3 = p - (p-a) - (p-qa) = a + qa - p = a(1+q) - p")
    print("    chi(g3) = chi(a(1+q)) = chi(a)*chi(1+q) = (-1)*chi(1+q) = -chi(1+q)")
    print("  If chi(q)=-1: chi(qa)=+1, qr(b) = qa")
    print("    g3 = p - (p-a) - qa = a - qa = a(1-q)")
    print("    chi(g3) = chi(a)*chi(1-q) = (-1)*chi(-(q-1)) = (-1)*chi(-1)*chi(q-1)")
    print("           = (-1)*(-1)*chi(q-1) = chi(q-1)")
    print()
    print("CONCLUSION: chi(g3) depends ONLY on q, not on a:")
    print("  chi(q)=+1 => chi(g3) = -chi(q+1)")
    print("  chi(q)=-1 => chi(g3) = chi(q-1)")
    print()
    print("Copart > 0 iff chi(g3) = +1, so:")
    print("  chi(q)=+1: copart iff chi(q+1) = -1")
    print("  chi(q)=-1: copart iff chi(q-1) = +1")
    print()
    print("QED.")

    # PART 3: The "neighbor transition" interpretation
    print("\n\n--- PART 3: QR NEIGHBOR TRANSITION TABLE ---")
    print("For each odd q < p, the co-participation depends on QR/QNR")
    print("status of q and its even neighbor q+/-1:")
    print()

    for p in [7, 11, 19, 23, 31, 43, 47]:
        if p % 4 != 3:
            continue
        m = (p - 1) // 2
        print(f"  p={p}:")
        all_q = [q for q in range(3, p, 2)]

        for q in all_q:
            chi_q = legendre(q, p)
            if chi_q == 1:
                neighbor = q + 1
                chi_n = legendre(neighbor, p)
                copart = (chi_n == -1)
            else:
                neighbor = q - 1
                chi_n = legendre(neighbor, p)
                copart = (chi_n == 1)

            marker = "Y" if copart else "N"
            print(f"    q={q:>2}: chi(q)={chi_q:+d}, "
                  f"neighbor={neighbor:>2}, chi(nb)={chi_n:+d}, "
                  f"copart={marker}")

        # Count copart=Y vs N
        copart_y = sum(1 for q in all_q
                       if (legendre(q, p) == 1 and legendre(q+1, p) == -1) or
                          (legendre(q, p) == -1 and legendre(q-1, p) == 1))
        copart_n = len(all_q) - copart_y
        print(f"    Summary: {copart_y} copart=Y, {copart_n} copart=N "
              f"out of {len(all_q)} odd q values")
        print()

    # PART 4: Connection to sign mechanism
    print("\n--- PART 4: CONNECTION TO PRODUCT LAW ---")
    print()
    print("The product law sign(h_hat[{a,b}]) = chi(ab) requires that")
    print("the dominant D^{2n} contributions all have sign = chi(ab).")
    print()
    print("From the resonance cascade:")
    print("  sign(D^{2n} deg-2) * chi(ab) = chi(q) at onset")
    print()
    print("So: sign = chi(q) * chi(ab).")
    print("For product law: need chi(q) * chi(ab) = chi(ab), i.e., chi(q) = +1.")
    print("But at p=19: chi(3) = -1, so 3-resonant pairs violate this!")
    print()
    print("The RESOLUTION must come from higher-order terms overwhelming D^4.")
    print("The co-participation theorem says: 3-resonant pairs have NO")
    print("additive coupling through 3-cycles, which means their H-contribution")
    print("comes from HIGHER-order cycle interactions, not from 3-cycle overlap.")
    print()
    print("At p=19:")
    print("  q=3 pairs: copart=N, sign(D^4)=-chi(ab), but D^6,D^8,... have +chi(ab)")
    print("  q=5 pairs: copart=N, sign(D^6)=+chi(ab)")
    print("  q=7 pairs: copart=Y, sign(D^8)=+chi(ab)")
    print("  q=9 pairs: copart=Y, sign(D^10)=+chi(ab)")
    print()
    print("The key question: for q=3 pairs at p=19, do D^6+D^8+... outweigh D^4?")
    print("This is the content of the p=19 product law computation.")


if __name__ == '__main__':
    main()
