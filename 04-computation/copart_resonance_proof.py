#!/usr/bin/env python3
"""
copart_resonance_proof.py -- Prove the chord co-participation vs resonance theorem

DISCOVERY from overlap_resonance_bridge.py:
- At p=11 Paley: copart(a,b) = 0 for q=3 resonant pairs, = p for q=5 resonant pairs
- This is EXACT (not approximate)

CONJECTURE (HYP-600): For Paley T_p at p=3 mod 4:
  chord_copart(a,b) in 3-cycles = { p, if a and b appear in same zero-sum triple in QR
                                    { 0, if not

Since zero-sum triples in QR are determined by additive structure, and resonance
is multiplicative structure, this creates a BRIDGE:

  copart = p * 1_{additive relationship exists}
  resonance level q = multiplicative relationship

PROOF STRATEGY:
For Paley at p=3 mod 4, QR is "additively flat" (|QR_hat(t)|^2 = (p+1)/4 for all t!=0).
A zero-sum triple (a,b,c) in QR means a+b+c=0 with a,b,c in QR.
The chord indices are min(a,p-a), min(b,p-b), min(c,p-c).

For chords i,j: copart(i,j) > 0 iff there exists a zero-sum triple
{g1,g2,g3} in QR with {min(g1,p-g1), min(g2,p-g2)} = {i,j}.

This is equivalent to: there exists s in QR such that
{i, j, s} or {i, p-j, s} or {p-i, j, s} or {p-i, p-j, s} sums to 0 mod p.

Author: kind-pasteur-2026-03-12-S60
"""

import math
from collections import defaultdict
from itertools import combinations


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


def zero_sum_triples_in_QR(p):
    """Find all zero-sum triples (a,b,c) with a,b,c in QR, a+b+c=0 mod p.
    Return as sorted tuples to avoid duplicates."""
    QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
    triples = set()
    for a in QR:
        for b in QR:
            c = (p - a - b) % p
            if c in QR and c != 0:
                triples.add(tuple(sorted([a, b, c])))
    return sorted(triples), QR


def main():
    print("=" * 70)
    print("CHORD CO-PARTICIPATION vs RESONANCE: PROOF ATTEMPT")
    print("=" * 70)

    # PART 1: Verify the pattern at multiple primes
    print("\n--- PART 1: VERIFY PATTERN AT ALL p=3 mod 4 ---")

    for p in [7, 11, 19, 23]:
        if p % 4 != 3:
            continue

        m = (p - 1) // 2
        QR = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        QR_set = set(QR)

        print(f"\n  p={p}, m={m}, QR={QR}")

        # Find zero-sum triples in QR
        triples, _ = zero_sum_triples_in_QR(p)
        print(f"  Zero-sum triples in QR: {len(triples)}")

        # For each triple, extract chord pair
        chord_pairs_from_triples = defaultdict(int)
        for t in triples:
            # Each triple (a,b,c) has 3 "gap values"; chord index = min(g, p-g)
            chords = set()
            for g in t:
                chords.add(min(g, p - g))
            chords = sorted(chords)
            # Add all pairs from this chord set
            for i in range(len(chords)):
                for j in range(i + 1, len(chords)):
                    chord_pairs_from_triples[(chords[i], chords[j])] += 1

        # Compare with resonance classification
        print(f"\n  Chord pair analysis:")
        for a in range(1, m + 1):
            for b in range(a + 1, m + 1):
                res = classify_resonance(a, b, p)
                min_q = min(q for q, t in res) if res else 'inf'
                chi_ab = legendre(a * b, p)

                cp_from_zst = chord_pairs_from_triples.get((a, b), 0)

                # Expected copart: p * (number of gap triple types using this chord pair)
                # Actually: each zero-sum triple contributes p 3-cycles (by circulant symmetry)
                # And each 3-cycle is a DIRECTED cycle, so we count each direction separately

                print(f"    ({a},{b}): min_q={min_q:>3}, chi={chi_ab:+d}, "
                      f"zst_chord_pairs={cp_from_zst:>3}, "
                      f"expected_copart={cp_from_zst * p:>4}")

    # PART 2: Algebraic proof of the connection
    print("\n\n--- PART 2: ALGEBRAIC CONNECTION ---")
    print("For chords a,b (with 1 <= a < b <= m):")
    print("Zero-sum triple in QR using chords a,b exists iff")
    print("there is a QR element c such that one of:")
    print("  +a + b + c = 0,  +a + b - c = 0 (i.e., c = a+b or c = -(a+b))")
    print("  +a - b + c = 0,  +a - b - c = 0 (i.e., c = b-a or c = a-b)")
    print("  -a + b + c = 0,  -a + b - c = 0 (same as above by symmetry)")
    print("  -a - b + c = 0,  -a - b - c = 0 (i.e., c = a+b, same as first)")
    print()
    print("So the conditions reduce to:")
    print("  (a+b) mod p in QR, OR (b-a) mod p in QR, OR")
    print("  (p-a-b) mod p in QR, OR (a-b) mod p in QR")
    print()
    print("Since QR is closed under p-x at p=3 mod 4 (chi(-1)=-1 means x in QR iff p-x in QNR),")
    print("wait - chi(-1) = (-1)^{(p-1)/2} = -1 at p=3 mod 4. So -x is QNR when x is QR.")
    print("Therefore (p-a-b) in QR iff -(a+b) in QR iff (a+b) in QNR.")
    print("And (a-b) in QR iff (b-a) in QNR.")
    print()
    print("Conditions for zero-sum triple using chords a,b:")
    print("  (a+b) in QR  OR  (b-a) in QR  (the other two are complementary)")

    for p in [7, 11, 19, 23]:
        if p % 4 != 3:
            continue
        m = (p - 1) // 2
        QR_set = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)

        print(f"\n  p={p}:")
        for a in range(1, m + 1):
            for b in range(a + 1, m + 1):
                sum_ab = (a + b) % p
                diff_ab = (b - a) % p
                sum_in_QR = sum_ab in QR_set
                diff_in_QR = diff_ab in QR_set

                res = classify_resonance(a, b, p)
                min_q = min(q for q, t in res) if res else 'inf'

                has_zst = sum_in_QR or diff_in_QR

                # Also check for the "gap triple" version (need DIRECTED sums)
                # A 3-cycle gap triple (g1,g2,g3) with g1+g2+g3 = p (since gaps are positive)
                # Wait, gaps in a directed 3-cycle: g1+g2+g3 = 0 mod p, but gi > 0
                # So g1+g2+g3 = p. Each gi in S = QR.
                # Chord indices: min(gi, p-gi) for each gi
                # For chords a,b: need gi such that min(gi, p-gi) = a and another gj with min(gj, p-gj) = b
                # gi = a or gi = p-a; gj = b or gj = p-b
                # Third gap = p - gi - gj must also be in QR

                actual_zst = False
                for sa in [a, p - a]:
                    if sa not in QR_set:
                        continue
                    for sb in [b, p - b]:
                        if sb not in QR_set:
                            continue
                        sc = (p - sa - sb) % p
                        if sc in QR_set and sc > 0:
                            actual_zst = True

                print(f"    ({a},{b}): q={min_q:>3}, "
                      f"a+b={sum_ab:>2}({'QR' if sum_in_QR else 'QNR':>3}), "
                      f"b-a={diff_ab:>2}({'QR' if diff_in_QR else 'QNR':>3}), "
                      f"has_zst={'Y' if has_zst else 'N'}, "
                      f"actual_zst={'Y' if actual_zst else 'N'}, "
                      f"{'MATCH' if has_zst == actual_zst else '**MISMATCH**'}")

    # PART 3: What determines the q-level from the additive structure?
    print("\n\n--- PART 3: RESONANCE LEVEL FROM ADDITIVE STRUCTURE ---")
    print("If copart(a,b)=0, then (a+b) and (b-a) are both QNR (or one is 0).")
    print("What is the multiplicative meaning of this condition?")
    print()
    print("chi(a+b) = -1 and chi(b-a) = -1")
    print("chi(a+b)*chi(b-a) = 1")
    print("chi((a+b)(b-a)) = 1, i.e., chi(b^2 - a^2) = 1")
    print()
    print("For q-resonance: qa = +/-b, so b = qa or b = -qa = (p-q)a.")
    print("Then b^2 - a^2 = (q^2 - 1)*a^2, so chi(b^2-a^2) = chi(q^2-1).")
    print()
    print("For copart=0: chi(b^2-a^2) = 1, so chi(q^2-1) = 1.")
    print("For q=3: q^2-1 = 8. chi(8) = chi(2)^3.")

    for p in [7, 11, 19, 23]:
        if p % 4 != 3:
            continue
        m = (p - 1) // 2

        print(f"\n  p={p}:")
        for q in [3, 5, 7, 9, 11, 13]:
            if q >= p:
                break
            chi_q21 = legendre(q*q - 1, p)
            chi_q = legendre(q, p)
            print(f"    q={q}: chi(q^2-1)={chi_q21:+d}, chi(q)={chi_q:+d}, "
                  f"q^2-1={q*q-1}")

    # PART 4: The sign mechanism at each resonance level
    print("\n\n--- PART 4: COPART vs D^{2n} SIGN ---")
    print("Connection: copart=0 (no additive coupling) correlates with")
    print("the D^4 contribution being nonzero (3-resonance = q=3)")
    print("copart=p (additive coupling exists) correlates with")
    print("the D^4 contribution being zero (higher-q resonance)")
    print()
    print("This is because 3-resonance means a/b = 3 (multiplicatively),")
    print("which implies a+b and b-a are BOTH QNR (no zero-sum triple),")
    print("while 5-resonance or higher allows zero-sum triples to exist.")

    for p in [7, 11, 19]:
        m = (p - 1) // 2
        QR_set = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)

        print(f"\n  p={p}:")
        for a in range(1, m + 1):
            for b in range(a + 1, m + 1):
                res = classify_resonance(a, b, p)
                min_q = min(q for q, t in res) if res else 'inf'
                chi_ab = legendre(a * b, p)

                # Check if a+b or b-a in QR (simplified copart test)
                sum_qr = (a + b) % p in QR_set
                diff_qr = (b - a) % p in QR_set

                # For actual copart: need to check all sign combinations
                actual_copart = False
                for sa in [a, p - a]:
                    if sa not in QR_set:
                        continue
                    for sb in [b, p - b]:
                        if sb not in QR_set:
                            continue
                        sc = (p - sa - sb) % p
                        if sc in QR_set and sc > 0:
                            actual_copart = True

                print(f"    ({a},{b}): q={min_q:>3}, chi={chi_ab:+d}, "
                      f"copart={'Y' if actual_copart else 'N':>1}, "
                      f"chi(q^2-1)={legendre(min_q**2-1, p) if min_q != float('inf') else 'N/A':>3}")


if __name__ == '__main__':
    main()
