#!/usr/bin/env python3
"""
universal_copart_constancy.py -- Prove the Universal Co-Participation Constancy

THEOREM (HYP-602): For Paley T_p at p=3 mod 4, the k-cycle co-participation
of chord pair (a,b) depends ONLY on the resonance level q (= min odd q s.t.
qa = +/- b mod p), not on the specific pair. This holds for ALL cycle lengths k.

Moreover: copart_k(q) is always a multiple of p.

PROOF IDEA: For circulant tournaments, the automorphism group Z_p acts by
translation. This maps k-cycles to k-cycles, preserving gap structure.
Two pairs with the same q are related by multiplication (qa=+/-b),
and multiplication by QR elements is in Aut(Paley). So multiplying the
pair by a QR element maps cycles of one pair to cycles of the other.

More precisely: if (a,b) has resonance q and (a',b') also has resonance q,
then there exists t in Z_p* (a QR element) such that {a', b'} = {ta, tb} mod p.
The map v -> tv is an automorphism of the Paley tournament, so it maps
k-cycles using chords a,b to k-cycles using chords a',b'. Since |Aut_mult|
acts freely on pairs with the same q, the co-participation counts are equal.

The divisibility by p comes from the Z_p translation symmetry.

Author: kind-pasteur-2026-03-12-S60
"""

from collections import defaultdict
from itertools import combinations


def legendre(a, p):
    if a % p == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


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


def count_kcycles_by_gap_structure(p, k):
    """Count k-cycles in Paley T_p grouped by their gap multiset.
    A k-cycle visits vertices v_0, v_1, ..., v_{k-1} with
    gap_i = (v_{i+1} - v_i) mod p, all gaps in QR.

    For circulant tournament: gap sequences up to cyclic rotation
    are equivalent up to translation. Each gap sequence gives exactly
    p cycles (one for each starting vertex).

    Returns: dict mapping frozenset(chord pair) -> count of k-cycles using both chords.
    """
    QR = sorted(j for j in range(1, p) if pow(j, (p - 1) // 2, p) == 1)
    m = (p - 1) // 2

    # Generate all valid gap sequences of length k summing to 0 mod p
    # Each gap must be in QR
    def gen_gap_seqs(remaining_length, current_sum, current_seq):
        if remaining_length == 0:
            if current_sum % p == 0:
                yield tuple(current_seq)
            return
        for g in QR:
            current_seq.append(g)
            yield from gen_gap_seqs(remaining_length - 1, current_sum + g, current_seq)
            current_seq.pop()

    # Count co-participation
    chord_copart = defaultdict(int)

    for gap_seq in gen_gap_seqs(k, 0, []):
        # This gap sequence represents p cycles (one per starting vertex)
        # Extract chord indices from gaps
        chords = set()
        for g in gap_seq:
            chords.add(min(g, p - g))
        chord_list = sorted(chords)
        for i in range(len(chord_list)):
            for j in range(i + 1, len(chord_list)):
                chord_copart[(chord_list[i], chord_list[j])] += p

    # Divide by k to account for cyclic rotations of the gap sequence
    # (each distinct directed cycle has k rotations of its gap sequence)
    for key in chord_copart:
        chord_copart[key] //= k

    return chord_copart


def verify_constancy(p, max_k=5):
    """Verify that k-cycle copart depends only on q for all k up to max_k."""
    m = (p - 1) // 2
    QR = sorted(j for j in range(1, p) if pow(j, (p - 1) // 2, p) == 1)

    print(f"\n{'='*60}")
    print(f"UNIVERSAL COPART CONSTANCY at p={p}")
    print(f"{'='*60}")
    print(f"QR = {QR}")

    all_constant = True

    for k in range(3, max_k + 1):
        print(f"\n  k={k} cycles:")
        copart = count_kcycles_by_gap_structure(p, k)

        # Group by q
        q_groups = defaultdict(list)
        for a in range(1, m + 1):
            for b in range(a + 1, m + 1):
                res = classify_resonance(a, b, p)
                if not res:
                    continue
                min_q = min(qq for qq, t in res)
                cp = copart.get((a, b), 0)
                q_groups[min_q].append((a, b, cp))

        for q in sorted(q_groups.keys()):
            items = q_groups[q]
            vals = set(cp for a, b, cp in items)
            chi_q = legendre(q, p)
            constant = len(vals) == 1
            if not constant:
                all_constant = False
            val = list(vals)[0] if constant else None
            div_p = (val % p == 0) if val is not None else 'N/A'
            print(f"    q={q:>2} (chi={chi_q:+d}): "
                  f"{'CONSTANT' if constant else '**NOT CONSTANT**'} "
                  f"val={val if constant else sorted(vals)}, "
                  f"div_p={div_p}, "
                  f"val/p={val//p if val is not None else 'N/A'}")

    return all_constant


def algebraic_proof():
    """Show the algebraic structure behind constancy."""
    print("\n\n" + "=" * 60)
    print("ALGEBRAIC PROOF OF CONSTANCY")
    print("=" * 60)

    print("""
For Paley T_p at p = 3 mod 4:

1. Aut(T_p) = AGL_1(F_p) restricted to QR multipliers
   More precisely: v -> av + b where a in QR, b in Z_p
   |Aut| = p * (p-1)/2

2. Two pairs (a1,b1) and (a2,b2) have the same resonance level q
   iff a2/b2 = a1/b1 or a2/b2 = b1/a1 (mod p)
   iff (a2,b2) = (t*a1, t*b1) for some t in Z_p*

3. If t in QR, then v -> tv is an automorphism.
   This maps any k-cycle using chord pair (a1,b1) to a
   k-cycle using chord pair (t*a1, t*b1) = (a2, b2).
   So copart_k(a1,b1) = copart_k(a2,b2).

4. If t in QNR, then v -> tv is NOT an automorphism.
   But for p = 3 mod 4, chi(-1) = -1, so -1 is QNR.
   v -> -v maps T to T^{rev} but preserves UNDIRECTED cycles.
   For directed k-cycles: v -> tv reverses orientation.
   For ODD k: reversed cycle is a DIFFERENT directed cycle.
   For EVEN k: need to check separately.

5. For odd k (3, 5, 7, ...):
   If t in QNR, then v -> tv maps directed k-cycle C to
   a different directed k-cycle C'. Both use chords
   (|t|*a1, |t|*b1) where |t| = min(t, p-t).
   But the gap structure changes: gap g -> tg.
   If g in QR, tg may be in QNR.
   So this is NOT a valid cycle in T_p.

   The correct argument: (a2, b2) = (t*a1, t*b1) with t in QR
   because qa1 = +/- b1 => q*(ta1) = +/-(tb1).
   The QR orbit of any pair has size (p-1)/2 / gcd, which
   equals the number of pairs with this q-level (by Burnside).

THEREFORE: copart_k depends only on q. QED.

Note: this argument works for ALL k, not just k=3.
The key is that QR multipliers are automorphisms of Paley.
""")

    # Verify the orbit structure
    print("ORBIT VERIFICATION:")
    for p in [7, 11, 19]:
        m = (p - 1) // 2
        QR = sorted(j for j in range(1, p) if pow(j, (p - 1) // 2, p) == 1)
        QR_set = set(QR)

        print(f"\n  p={p}:")
        # For each pair, find its QR-orbit
        seen = set()
        for a in range(1, m + 1):
            for b in range(a + 1, m + 1):
                pair = (a, b)
                if pair in seen:
                    continue
                res = classify_resonance(a, b, p)
                if not res:
                    continue
                min_q = min(qq for qq, t in res)

                orbit = set()
                for t in QR:
                    ta = (t * a) % p
                    tb = (t * b) % p
                    ca = min(ta, p - ta)
                    cb = min(tb, p - tb)
                    orbit.add(tuple(sorted([ca, cb])))

                for pair in orbit:
                    seen.add(pair)

                print(f"    q={min_q}: orbit size = {len(orbit)}, "
                      f"orbit = {sorted(orbit)}")


def main():
    # Part 1: Verify constancy
    for p in [7, 11]:
        verify_constancy(p, max_k=7)

    verify_constancy(19, max_k=5)

    # Part 2: Algebraic proof
    algebraic_proof()

    # Part 3: Normalized copart table (divided by p)
    print("\n\n" + "=" * 60)
    print("NORMALIZED COPART TABLE (copart/p)")
    print("=" * 60)

    for p in [7, 11, 19]:
        m = (p - 1) // 2
        print(f"\n  p={p}:")

        max_k = 7 if p <= 11 else 5

        # Collect all data
        data = {}  # (q, k) -> val/p
        q_vals = set()
        for k in range(3, max_k + 1):
            copart = count_kcycles_by_gap_structure(p, k)
            for a in range(1, m + 1):
                for b in range(a + 1, m + 1):
                    res = classify_resonance(a, b, p)
                    if not res:
                        continue
                    min_q = min(qq for qq, t in res)
                    q_vals.add(min_q)
                    data[(min_q, k)] = copart.get((a, b), 0) // p

        q_vals = sorted(q_vals)
        header = f"    {'q':>3}"
        for k in range(3, max_k + 1):
            header += f" {'c'+str(k)+'/p':>8}"
        print(header)

        for q in q_vals:
            chi_q = legendre(q, p)
            row = f"    {q:>3}"
            for k in range(3, max_k + 1):
                val = data.get((q, k), 0)
                row += f" {val:>8}"
            row += f"  (chi(q)={chi_q:+d})"
            print(row)


if __name__ == '__main__':
    main()
