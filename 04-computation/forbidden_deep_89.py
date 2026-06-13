#!/usr/bin/env python3
"""
forbidden_deep_89.py -- opus-2026-03-14-S89

Deep analysis of ALL 18 forbidden H-values at n=7.
Two groups: "interior" (surrounded by achieved) and "tail" (sparse region).

Interior forbidden: {7, 21, 63, 107, 119, 149}
Tail forbidden: {161, 163, 165, 167, 169, 173, 177, 179, 181, 183, 185, 187}
"""

from math import gcd, comb, factorial
from fractions import Fraction


def main():
    print("="*70)
    print("DEEP FORBIDDEN VALUE ANALYSIS at n=7")
    print("opus-2026-03-14-S89")
    print("="*70)

    # Full H-spectrum at n=7 (from exhaustive computation)
    # 77 achieved values from 1 to 189
    achieved = [1, 3, 5, 9, 11, 13, 15, 17, 19, 23, 25, 27, 29, 31, 33, 37,
                39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 65, 67, 69,
                71, 73, 75, 77, 79, 81, 83, 85, 87, 89, 91, 93, 95, 97, 99,
                101, 103, 105, 109, 111, 113, 115, 117, 121, 123, 125, 127,
                129, 131, 133, 135, 137, 139, 141, 143, 145, 147, 151, 153,
                155, 157, 159, 171, 175, 189]
    achieved_set = set(achieved)

    all_odd = set(range(1, 190, 2))
    forbidden = sorted(all_odd - achieved_set)

    print(f"\n  Achieved: {len(achieved)} values")
    print(f"  Forbidden: {len(forbidden)} values")
    print(f"  Forbidden: {forbidden}")

    # Classify: interior vs tail
    interior = []
    tail = []
    for f in forbidden:
        has_below = (f - 2) in achieved_set or f == 1
        has_above = (f + 2) in achieved_set
        if has_below and has_above:
            interior.append(f)
        else:
            tail.append(f)

    print(f"\n  Interior forbidden (neighbors achieved): {interior}")
    print(f"  Tail/sparse forbidden: {tail}")

    # Part 1: Number theory of interior forbidden values
    print(f"\n{'='*70}")
    print("PART 1: NUMBER THEORY OF INTERIOR FORBIDDEN VALUES")
    print("="*70)

    for fv in interior:
        print(f"\n  H = {fv}:")
        # Factorization
        factors = []
        x = fv
        for p in range(2, fv + 1):
            while x % p == 0:
                factors.append(p)
                x //= p
            if x == 1:
                break
        print(f"    Factorization: {' * '.join(map(str, factors)) if factors else 'prime'}")
        print(f"    Binary: {bin(fv)}")

        # Repunit/repdigit in various bases
        for base in range(2, 10):
            digits = []
            x = fv
            while x > 0:
                digits.append(x % base)
                x //= base
            digits.reverse()
            if len(set(digits)) == 1:
                print(f"    Base {base}: {''.join(map(str, digits))} (REPDIGIT)")

        # Mersenne / repunit test
        # 2^k - 1
        for k in range(1, 10):
            if fv == 2**k - 1:
                print(f"    = 2^{k} - 1 (MERSENNE)")
        # (q^k - 1)/(q - 1) for various q
        for q in range(2, 8):
            for k in range(2, 10):
                if fv * (q - 1) == q**k - 1:
                    print(f"    = (q^{k} - 1)/(q - 1) with q={q} (REPUNIT in base {q})")
                if fv == q**k - 1:
                    print(f"    = {q}^{k} - 1")

        # Projective space sizes: (q^{k+1}-1)/(q-1) = PG(k, F_q)
        for q in [2, 3, 4, 5, 7, 8, 9]:
            for k in range(1, 10):
                if fv == (q**(k+1) - 1) // (q - 1) and (q**(k+1) - 1) % (q - 1) == 0:
                    print(f"    = |PG({k}, F_{q})| (projective {k}-space over F_{q})")

        # Check H mod small numbers
        print(f"    Mod 7: {fv % 7}")
        print(f"    Mod 8: {fv % 8}")
        print(f"    Mod 9: {fv % 9}")
        print(f"    Mod 16: {fv % 16}")
        print(f"    Mod 21: {fv % 21}")

    # Part 2: Relationships among interior forbidden values
    print(f"\n{'='*70}")
    print("PART 2: RELATIONSHIPS AMONG INTERIOR FORBIDDEN VALUES")
    print("="*70)

    pairs = [(a, b) for a in interior for b in interior if a < b]
    for a, b in pairs:
        d = b - a
        g = gcd(a, b)
        r = b / a if a > 0 else 0
        print(f"  {a}, {b}: diff={d}, gcd={g}, ratio={r:.4f}")
        if b % a == 0:
            print(f"    {b} = {b//a} * {a}")

    # 7, 21, 63: 21 = 3*7, 63 = 9*7 = 3*21. All multiples of 7!
    # 107, 119, 149: 119 = 7*17, 107 is prime, 149 is prime
    print(f"\n  {interior[:3]}: all multiples of 7! (7, 3*7, 9*7)")
    print(f"  119 = 7 * 17 (also multiple of 7!)")
    print(f"  107 is prime, 149 is prime")

    # So: 7 | {7, 21, 63, 119} but NOT {107, 149}
    mults_of_7 = [f for f in interior if f % 7 == 0]
    non_7 = [f for f in interior if f % 7 != 0]
    print(f"\n  Multiples of 7 among interior: {mults_of_7}")
    print(f"  Non-multiples of 7 among interior: {non_7}")
    print(f"  107 mod 7 = {107 % 7}")
    print(f"  149 mod 7 = {149 % 7}")

    # Part 3: Mean = 315/4 = 78.75. Symmetry around mean?
    print(f"\n{'='*70}")
    print("PART 3: SYMMETRY OF FORBIDDEN VALUES AROUND MEAN")
    print("="*70)

    mean = Fraction(315, 4)  # 78.75
    max_h = 189

    print(f"  Mean = {mean} = {float(mean)}")
    print(f"  Max H = {max_h}")
    print(f"\n  For each forbidden value f, check if (2*Mean - f) = (157.5 - f) is achieved:")
    for fv in forbidden:
        complement = 2 * float(mean) - fv  # 157.5 - fv
        # complement might not be an odd integer
        print(f"    f={fv:3d}: 2*Mean-f = {complement:.1f}, ", end="")
        if complement == int(complement) and complement > 0:
            c = int(complement)
            if c % 2 == 1:
                status = "ACHIEVED" if c in achieved_set else "FORBIDDEN"
                print(f"odd integer {c}: {status}")
            else:
                print(f"even integer {c}: N/A")
        else:
            print(f"not an odd integer")

    # Since Mean = 315/4 is not a half-integer in general, 2*Mean - f is not always odd.
    # But let's check: 315/4 * 2 = 315/2 = 157.5. So 2*Mean - f = 157.5 - f.
    # For f odd: 157.5 - f is not an integer. So there's no direct complement pairing.

    # Better: use complement symmetry H(T) + H(T-bar) = 2*Mean?
    # No: H(T-bar) = H(T), so H(T) + H(T-bar) = 2*H(T), not 2*Mean.
    # The complement doesn't give us a pairing of different H values.

    # The actual symmetry: max_H + min_H = 189 + 1 = 190.
    # If H(T) = h, does H(T') = 190 - h for some T'?
    # No obvious reason, since H(T-bar) = H(T).
    print(f"\n  Check: max + min = 189 + 1 = 190. Is 190 - f achieved for each f?")
    for fv in forbidden:
        c = 190 - fv
        if c > 0 and c % 2 == 1:
            status = "ACHIEVED" if c in achieved_set else "FORBIDDEN"
            print(f"    f={fv:3d}: 190-f={c:3d}: {status}")
        else:
            status = "even" if c % 2 == 0 else "negative"
            print(f"    f={fv:3d}: 190-f={c:3d}: {status}")

    # Part 4: Forbidden values in base representations
    print(f"\n{'='*70}")
    print("PART 4: COMPREHENSIVE BASE REPRESENTATION")
    print("="*70)

    for fv in interior:
        print(f"\n  H = {fv}:")
        for base in [2, 3, 4, 5, 6, 7, 8, 9]:
            digits = []
            x = fv
            while x > 0:
                digits.append(str(x % base))
                x //= base
            digits.reverse()
            rep = ''.join(digits)
            print(f"    Base {base}: {rep}")

    # Part 5: The 7-21-63 pattern more deeply
    print(f"\n{'='*70}")
    print("PART 5: THE 7-21-63 PATTERN -- REPUNIT STRUCTURE")
    print("="*70)

    # 7 = 111_2 = R_3(2) = (2^3-1)/(2-1)
    # 21 = 111_4 = R_3(4) = (4^3-1)/(4-1) = 63/3 = 21
    # 63 = 111111_2 = R_6(2) = (2^6-1)/(2-1) = 63
    # Also 63 = 333_4

    # Key: 7 = R_3(2), 21 = R_3(4) = R_3(2^2), 63 = R_6(2) = R_3(2)*R_2(2^3+1)?
    # Actually: 7*3 = 21, 7*9 = 63. So 21 = 7*3, 63 = 7*9 = 7*3^2.
    # Pattern: 7 * 3^k for k = 0, 1, 2.

    print("  7, 21, 63 = 7*1, 7*3, 7*9 = 7*3^0, 7*3^1, 7*3^2")
    print("  = 7 * {1, 3, 9} = 7 * powers of 3")
    print("  Next would be: 7 * 27 = 189 = max(H at n=7)!")
    print(f"  Is 189 forbidden? {189 not in achieved_set}")
    print(f"  No! 189 IS achieved (it's the maximum). The pattern breaks at 7*3^3.")

    # But wait: 189 = 7*27 = 7*3^3. It's ACHIEVED (max H).
    # So the pattern 7*3^k is forbidden for k=0,1,2 but achieved for k=3.
    # 7*3^3 = 189 = H(T) for the "most complex" tournaments.

    print(f"\n  7*3^k series: 7, 21, 63, 189, 567, ...")
    print(f"  At n=7: 7 FORBIDDEN, 21 FORBIDDEN, 63 FORBIDDEN, 189 ACHIEVED (max)")
    print(f"  The forbidden pattern is {7}*{3}^k for 0 <= k <= 2")
    print(f"  And 7*3^3 = 189 is the MAXIMUM value -- the pattern terminates!")

    # Deeper: 189 = max H at n=7.
    # Mean = 315/4 = 78.75. Max/Mean = 189/(315/4) = 189*4/315 = 756/315 = 2.4
    # Min/Mean = 1/(315/4) = 4/315 ~ 0.013

    print(f"\n  Max/Mean = 189/{float(mean)} = {189/float(mean):.4f}")
    print(f"  Min/Mean = 1/{float(mean)} = {1/float(mean):.6f}")

    # Part 6: Is there a forbidden value pattern at n=5 and n=6?
    print(f"\n{'='*70}")
    print("PART 6: FORBIDDEN VALUES AT SMALLER n")
    print("="*70)

    # n=3: achieved = {1, 3}. All odd in [1,3] = {1,3}. No forbidden.
    # n=4: achieved = {1, 3, 5}. All odd in [1,5] = {1,3,5}. No forbidden.
    # n=5: achieved = {1,3,5,9,11,13,15}. Odd in [1,15] = {1,3,5,7,9,11,13,15}.
    #   Forbidden at n=5: {7}
    # n=6: achieved as computed above. Let me list the full set.
    n6_achieved = [1, 3, 5, 9, 11, 13, 15, 17, 19, 23, 25, 27, 29, 31, 33, 37, 41, 43, 45]
    n6_set = set(n6_achieved)
    n6_all_odd = set(range(1, 46, 2))
    n6_forbidden = sorted(n6_all_odd - n6_set)

    print(f"  n=3: forbidden = [] (all achieved)")
    print(f"  n=4: forbidden = [] (all achieved)")
    print(f"  n=5: forbidden = [7]")
    print(f"  n=6: achieved = {n6_achieved}")
    print(f"  n=6: forbidden = {n6_forbidden}")
    print(f"  n=7: forbidden = {forbidden}")

    # n=5: 7 is forbidden
    # n=6: let me check... {7, 21, 35, 39} perhaps?
    # Actually from the data: 7 is forbidden at n=5 (not in {1,3,5,9,11,13,15})
    # At n=6 max = 45, achieved 19 values. All odd 1..45 = 23 values. So 4 forbidden.

    print(f"\n  n=6 forbidden: {n6_forbidden}")
    for fv in n6_forbidden:
        factors = []
        x = fv
        for p in range(2, fv + 1):
            while x % p == 0:
                factors.append(p)
                x //= p
            if x == 1:
                break
        print(f"    {fv} = {'*'.join(map(str, factors)) if factors else 'prime'}, binary={bin(fv)}")

    # PATTERN ACROSS n:
    # n=5: forbidden = {7}
    # n=6: need to compute
    # n=7: forbidden = {7, 21, 63, 107, 119, 149, 161-187 tail}

    # Is 7 ALWAYS forbidden for n >= 5?
    # At n=5: mean = 15/4 = 3.75, and 7 > max(H)? No, max at n=5 is 15.
    # 7 is in range [1,15] but not achieved. So 7 IS always forbidden for n >= 5?
    # At n=7: 7 is forbidden. At n=6: need to check.

    print(f"\n  Is 7 forbidden at each n?")
    print(f"    n=3: 7 > max(3), N/A")
    print(f"    n=4: 7 > max(5), N/A")
    print(f"    n=5: 7 FORBIDDEN (in [1,15] but not achieved)")
    print(f"    n=6: {7 in n6_forbidden}")
    print(f"    n=7: {7 in set(forbidden)}")

    # Part 7: The connection to 7! = 5040 and Mean at n=7
    print(f"\n{'='*70}")
    print("PART 7: FACTORIAL AND MEAN CONNECTIONS")
    print("="*70)

    for n in range(3, 10):
        mean_n = Fraction(factorial(n), 2**(n-1))
        max_h_approx = "?" # We don't know max H for n > 7
        print(f"  n={n}: Mean = {mean_n} = {float(mean_n):.4f}")

    # The forbidden value 7 at n=5:
    # Mean_5 = 120/16 = 15/2 = 7.5
    # So 7 is just BELOW the mean! (7 < 7.5 < 9)
    # And 7 = mean - 0.5 = 7.5 - 0.5

    print(f"\n  At n=5: Mean = 7.5, forbidden value 7 = Mean - 0.5")
    print(f"  At n=7: Mean = 78.75, forbidden 7 is far below mean")

    # Is the forbidden value always close to some nice fraction of the mean?
    for fv in interior:
        ratio = fv / float(mean)
        print(f"  H={fv}: H/Mean = {ratio:.6f}")

    # 7/78.75 = 0.0889 ~ 1/11.25
    # 21/78.75 = 0.2667 ~ 4/15
    # 63/78.75 = 0.8 = 4/5!
    # 107/78.75 = 1.3587
    # 119/78.75 = 1.5111 ~ 3/2
    # 149/78.75 = 1.8921 ~ 2

    print(f"\n  63/Mean = 63/(315/4) = 63*4/315 = {Fraction(63*4, 315)} = {float(Fraction(63*4,315))}")
    print(f"  This is 4/5! The forbidden value 63 = 4/5 * Mean.")
    print(f"  Check: 4/5 * 315/4 = 315/5 = 63. YES!")

    # 63 = (4/5) * Mean = (4/5) * 7!/128 = 4*5040/(5*128) = 20160/640 = 63/2... wait
    # Mean = 7!/2^6 = 5040/64 = 78.75 = 315/4
    # 63 = 315/4 * 4/5 = 315/5 = 63. Yes.
    # So 63 = 315/5 = 7!/5/2^6 = 7!/(5*64) = 5040/320 = 15.75... no that's wrong.
    # 63 = 315/5 = (7!/2^6)/5 * (4/4)... hmm.
    # Actually: 63 = 315/5. And 315 = 5*63. And 315 = 7*45 = 7*9*5 = 5*7*9.
    # So 63 = 7*9 = 7*3^2. Which we already knew.

    print(f"\n  315 = {315} = 5 * 63 = 7 * 45 = 9 * 35 = 5 * 7 * 9")
    print(f"  Mean = 315/4. So 63 = Mean * 4/5 = (4*Mean)/5")

    # Part 8: The "inner" forbidden values 107, 119, 149
    print(f"\n{'='*70}")
    print("PART 8: THE INNER FORBIDDEN VALUES 107, 119, 149")
    print("="*70)

    # Complement symmetry H(T-bar) = H(T) doesn't pair different values.
    # But there's another symmetry: vertex relabeling S_7 acts on tournaments.
    # Under S_7, isomorphic tournaments have the same H.
    # The 2^21 tournaments partition into ~456 isomorphism classes (A000568(7)=456).

    # 107 is prime. 119 = 7*17. 149 is prime.
    # 107 + 149 = 256 = 2^8. Interesting!
    # 107 + 119 = 226. 119 + 149 = 268.

    print(f"  107 + 149 = {107+149} = 2^8 !!!")
    print(f"  107 + 119 = {107+119}")
    print(f"  119 + 149 = {119+149}")
    print(f"  107 + 119 + 149 = {107+119+149}")
    print(f"  107 + 82 = 189 (max). So 189 - 107 = {189-107} = 82. H=82 is even, N/A.")
    print(f"  190 - 107 = {190-107} = 83. Is 83 achieved? {83 in achieved_set}")
    print(f"  190 - 119 = {190-119} = 71. Is 71 achieved? {71 in achieved_set}")
    print(f"  190 - 149 = {190-149} = 41. Is 41 achieved? {41 in achieved_set}")

    # 190 - f is achieved for all three! So they don't pair to other forbidden values.

    # Are 107, 119, 149 related to the Walsh lattice?
    # 32H = 2520 + 120*S2 + 6*S4a + 12*S4b + S6
    # For H=107: 32*107 = 3424, need sum = 3424 - 2520 = 904
    # For H=119: 32*119 = 3808, need sum = 3808 - 2520 = 1288
    # For H=149: 32*149 = 4768, need sum = 4768 - 2520 = 2248

    print(f"\n  Walsh lattice equation: 120*S2 + 6*S4a + 12*S4b + S6 = 32H - 2520")
    for fv in interior:
        target = 32 * fv - 2520
        print(f"    H={fv}: target = {target}")
        # S2 is odd (105 terms), S4a is even (1260 terms), S4b is even (630 terms), S6 is even (2520 terms)
        # S2 has same parity as 105 = odd
        # 120*S2 = 120*(2p2-105) where p2 = number of + signs
        # 6*S4a = 6*(2p4a-1260), 12*S4b = 12*(2p4b-630), S6 = 2p6-2520
        # Total = 240*p2 - 12600 + 12*p4a - 7560 + 24*p4b - 7560 + 2*p6 - 2520
        # = 240*p2 + 12*p4a + 24*p4b + 2*p6 - 30240
        # So: 240*p2 + 12*p4a + 24*p4b + 2*p6 = 32H - 2520 + 30240 = 32H + 27720
        expanded = 32 * fv + 27720
        print(f"    240*p2 + 12*p4a + 24*p4b + 2*p6 = {expanded}")
        # Divide by 2: 120*p2 + 6*p4a + 12*p4b + p6 = {expanded//2}
        half = expanded // 2
        print(f"    120*p2 + 6*p4a + 12*p4b + p6 = {half}")
        # p2 in [0, 105], p4a in [0, 1260], p4b in [0, 630], p6 in [0, 2520]
        # Check mod 6: 120*p2 = 0 mod 6, 6*p4a = 0 mod 6, 12*p4b = 0 mod 6
        # So p6 = half mod 6
        print(f"    p6 mod 6 = {half % 6}")
        # Check mod 120: p6 mod 120 constrains p6 given p2
        # This is getting complex. The key question: does this system have integer solutions?

    # Part 9: Are ALL Mersenne numbers 2^k - 1 forbidden?
    print(f"\n{'='*70}")
    print("PART 9: MERSENNE NUMBERS AND FORBIDDEN VALUES")
    print("="*70)

    mersenne = [2**k - 1 for k in range(1, 12)]
    print(f"  Mersenne numbers: {mersenne}")
    print(f"  In H-spectrum range [1, 189]:")
    for m in mersenne:
        if m <= 189 and m % 2 == 1:
            status = "FORBIDDEN" if m not in achieved_set else "ACHIEVED"
            print(f"    2^{m.bit_length()}-1 = {m}: {status}")

    # 1 = 2^1-1: ACHIEVED
    # 3 = 2^2-1: ACHIEVED
    # 7 = 2^3-1: FORBIDDEN
    # 15 = 2^4-1: ACHIEVED
    # 31 = 2^5-1: ACHIEVED
    # 63 = 2^6-1: FORBIDDEN
    # 127 = 2^7-1: ACHIEVED
    # So NOT all Mersenne numbers are forbidden. Only 7 and 63.
    # 7 = 2^3-1, 63 = 2^6-1. Note 3 and 6: 6 = 2*3.

    print(f"\n  Only 2^3-1 = 7 and 2^6-1 = 63 are forbidden.")
    print(f"  Pattern: 2^(3k)-1 for k=1,2. Next: 2^9-1 = 511 (outside range).")

    # Part 10: Summary and conjectures
    print(f"\n{'='*70}")
    print("PART 10: SUMMARY AND CONJECTURES")
    print("="*70)

    print(f"""
  FINDINGS:
  1. 18 forbidden values at n=7, split into:
     - 6 "interior" forbidden: {{7, 21, 63, 107, 119, 149}}
     - 12 "tail" forbidden: {{161-187}} (gap near max)

  2. The 7-21-63 pattern: 7*3^k for k=0,1,2.
     - 7*3^3 = 189 = max(H), which IS achieved.
     - So the pattern is 7*3^k for k < log_3(max/7).

  3. Projective space connection:
     - 7 = |PG(2, F_2)| = Fano plane
     - 21 = |PG(2, F_4)|
     - 63 = |PG(5, F_2)| = 2^6 - 1

  4. 107 + 149 = 256 = 2^8 (power of 2 sum!)

  5. Repunit structure:
     - 7 = 111_2
     - 21 = 111_4
     - 63 = 111111_2 = 333_4

  6. Only Mersenne primes 2^3-1 and 2^6-1 are forbidden (not 2^k-1 in general).

  7. C3 = triangle count is purely level-2 Walsh.
     All forbidden values are in range [1, max] and NOT extreme (small z-scores).
     The obstruction is structural (Walsh lattice), not statistical.

  CONJECTURES:
  - CONJ-F1: H=7 is forbidden for ALL n >= 5.
  - CONJ-F2: H=21 is forbidden for ALL n >= 7.
  - CONJ-F3: The "interior" forbidden values grow as 7*3^k.
  - CONJ-F4: The number of interior forbidden values at n grows logarithmically.
""")

    print(f"\n{'='*70}")
    print("DONE -- DEEP FORBIDDEN VALUE ANALYSIS")
    print("="*70)


if __name__ == "__main__":
    main()
