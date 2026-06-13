#!/usr/bin/env python3
"""
eisenstein_spectrum_connection_89.py -- opus-2026-03-14-S89

POTENTIAL CROWN JEWEL: The Eisenstein norm sequence N(F_n, F_{n+1})
matches the H-spectrum sizes!

H-spectrum sizes: 1, 1, 2, 3, 7, 19, 77 (n=1..7)
Eisenstein norms: 1, 1, 3, 7, 19, 49, 129 (n=0..6)

Overlapping subsequence: 1, 1, ?, 3/7, 7/19, 19/??, 77/??
Need to check if these are related or coincidental.
"""

from math import gcd


def main():
    print("="*70)
    print("EISENSTEIN NORM vs H-SPECTRUM SIZE CONNECTION")
    print("opus-2026-03-14-S89")
    print("="*70)

    # Fibonacci sequence
    fib = [0, 1]
    for _ in range(30):
        fib.append(fib[-1] + fib[-2])

    # Eisenstein norm: N(a,b) = a^2 - a*b + b^2
    print("\n  Eisenstein norms of consecutive Fibonacci pairs:")
    eis_norms = []
    for n in range(20):
        a, b = fib[n], fib[n+1]
        norm = a*a - a*b + b*b
        eis_norms.append(norm)
        print(f"    n={n:2d}: N(F_{n}, F_{n+1}) = N({a}, {b}) = {norm}")

    # H-spectrum sizes
    h_spec = {1: 1, 2: 1, 3: 2, 4: 3, 5: 7, 6: 19, 7: 77}
    print(f"\n  H-spectrum sizes:")
    for n in sorted(h_spec.keys()):
        print(f"    n={n}: |Spec| = {h_spec[n]}")

    # Compare
    print(f"\n  COMPARISON:")
    print(f"  {'n':>3} {'|Spec(n)|':>10} {'N(F_k,F_{k+1})':>15} {'match':>6}")

    # The sequence 1,1,3,7,19 appears in the Eisenstein norms at indices 0,1,2,3,4
    # The sequence 1,1,2,3,7,19,77 is the H-spectrum at n=1..7
    # Overlap: 1,1 match trivially. Then: 3 vs 2, 7 vs 3, 19 vs 7, 49 vs 19, 129 vs 77

    # Actually the match is SHIFTED: Eisenstein norms at n=2,3,4 = 3,7,19
    # H-spectrum at n=4,5,6 = 3,7,19
    # So the claim: |Spec(n+2)| = N(F_{n}, F_{n+1}) for n=2,3,4?

    print(f"\n  Test: |Spec(n+2)| = N(F_n, F_{{n+1}}):")
    for n in range(10):
        spec = h_spec.get(n+2, "?")
        eis = eis_norms[n]
        match = spec == eis if spec != "?" else "?"
        print(f"    n={n}: N(F_{n},F_{n+1})={eis}, |Spec({n+2})|={spec}, match={match}")

    # n=0: N=1, |Spec(2)|=1. Match!
    # n=1: N=1, |Spec(3)|=2. NO MATCH (1 vs 2)
    # n=2: N=3, |Spec(4)|=3. Match!
    # n=3: N=7, |Spec(5)|=7. Match!
    # n=4: N=19, |Spec(6)|=19. Match!
    # n=5: N=49, |Spec(7)|=77. NO MATCH (49 vs 77)

    # So: matches at n=0,2,3,4 but fails at n=1 and n=5.
    # Not an exact match. But 3 consecutive matches (3,7,19) is striking.

    # Maybe the relationship is: |Spec(n)| ~ Eisenstein norm for small n,
    # but diverges. Let's check other sequences that contain 2,3,7,19.

    print(f"\n  Subsequence 2, 3, 7, 19, 77:")
    # Ratios: 3/2=1.5, 7/3=2.33, 19/7=2.71, 77/19=4.05
    # Products: 2*3=6, 3*7=21, 7*19=133, 19*77=1463

    # Could it be a(n) = a(n-1)*a(n-2) - 1?
    # 3*2-1=5 (not 7). No.
    # a(n) = a(n-1)*a(n-2) + a(n-3)?
    # 7*3+2=23 (not 19). No.

    # What about: a(n) = a(n-1)^2 - a(n-2)?
    # 3^2-2=7. YES!
    # 7^2-3=46. Not 19.
    # Nope.

    # a(n) = a(n-1) + a(n-2)*(a(n-2)-1)?
    # 3 + 2*1 = 5 (not 7). No.

    # a(n) = a(n-1)*2 + a(n-2)?
    # 3*2+2=8 (not 7). No.

    # a(n) = a(n-1) + a(n-1)*a(n-2)?
    # 3+3*2=9 (not 7). No.

    # Let's look at DIFFERENCES: 2,3,7,19,77
    # Differences: 1, 4, 12, 58
    # Second differences: 3, 8, 46
    # Not obvious.

    # Ratios more carefully: 77/19 = 4.0526... This doesn't fit a simple pattern.
    # Maybe it's 2,3,7,19,77,... is its own new sequence.

    # OEIS: Let's check 2,3,7,19,77
    print(f"\n  Sequence 2,3,7,19,77: checking patterns")
    s = [2, 3, 7, 19, 77]
    for i in range(2, len(s)):
        # Quadratic: s[i] = A*s[i-1] + B*s[i-2]
        # 7 = 3A + 2B
        # 19 = 7A + 3B
        # From first: B = (7-3A)/2
        # Sub: 19 = 7A + 3(7-3A)/2 = 7A + 21/2 - 9A/2 = 14A/2 - 9A/2 + 21/2 = 5A/2 + 21/2
        # 38 = 5A + 21 -> 5A = 17 -> A = 17/5. Not integer!
        pass

    # Try: s[i] = s[i-1]*s[i-2] - s[i-3]?
    # Need s[0]=a, s[1]=b, s[2]=c.
    # Let's set s[-1] = some value and try s[i] = s[i-1]*s[i-2] - s[i-3]
    # s[3] = s[2]*s[1] - s[0] = 7*3 - 2 = 19. YES!
    # s[4] = s[3]*s[2] - s[1] = 19*7 - 3 = 130. Not 77.

    print("  Test s[i] = s[i-1]*s[i-2] - s[i-3]:")
    print(f"    s[3] = 7*3 - 2 = {7*3-2} (expected 19): {'MATCH' if 7*3-2==19 else 'FAIL'}")
    print(f"    s[4] = 19*7 - 3 = {19*7-3} (expected 77): {'MATCH' if 19*7-3==77 else 'FAIL'}")

    # s[4] = 130, not 77. Close but no cigar.

    # Try: s[i] = s[i-1]*s[i-2]/s[i-3] * k?
    # 77 = 19*7/3 * k -> k = 77*3/(19*7) = 231/133 = 33/19. Not clean.

    # Maybe: a(n) satisfies a non-standard recurrence.
    # 2, 3, 7, 19, 77
    # Look at: a(n)*a(n-2) - a(n-1)^2:
    # 7*2 - 9 = 5
    # 19*3 - 49 = 8
    # 77*7 - 361 = 178
    # Not obvious.

    # a(n+1)*a(n-1) - a(n)^2:
    # 3*2 - 4 = 2
    # 7*2 - 9 = 5
    # 19*3 - 49 = 8
    # 77*7 - 361 = 178
    # Hmm: 2, 5, 8... arithmetic for a bit?

    print(f"\n  Somos/determinant-like relations:")
    for i in range(1, len(s)-1):
        det = s[i+1]*s[i-1] - s[i]**2
        print(f"    s[{i+1}]*s[{i-1}] - s[{i}]^2 = {s[i+1]}*{s[i-1]} - {s[i]}^2 = {det}")

    # 2, 5, 8: differences 3, 3. Then 178: huge jump. Not arithmetic.

    # The 77 = 7 * 11. Maybe the sequence is n=7 dependent.
    # H-spectrum at n=8 would tell us more.

    # Let's also look at |Spec| / sum of first few Fibonacci products
    print(f"\n  |Spec(n)| vs Fibonacci-based formulas:")
    for n in range(3, 8):
        spec = h_spec[n]
        # F(2n-3)?
        f_idx = 2*n - 3
        if f_idx >= 0 and f_idx < len(fib):
            print(f"    n={n}: |Spec|={spec}, F({f_idx})={fib[f_idx]}")

    # F(3)=2, F(5)=5, F(7)=13, F(9)=34, F(11)=89
    # vs: 2, 3, 7, 19, 77. Only n=3 matches.

    # Lucas numbers: L(n) = F(n-1) + F(n+1)
    lucas = [2, 1]
    for _ in range(20):
        lucas.append(lucas[-1] + lucas[-2])
    print(f"\n  Lucas numbers: {lucas[:15]}")
    print(f"  |Spec|: 1, 1, 2, 3, 7, 19, 77")
    # L: 2, 1, 3, 4, 7, 11, 18, 29, 47, 76, 123
    # L(4)=7 and |Spec(5)|=7. L(9)=76 and |Spec(7)|=77.
    # Almost! 77 = L(9)+1 = 76+1.

    print(f"  L(4) = {lucas[4]} vs |Spec(5)| = 7: match = {lucas[4] == 7}")
    print(f"  L(9) = {lucas[9]} vs |Spec(7)| = 77: diff = {77 - lucas[9]}")

    # Hmm, L(9)=76 and |Spec|=77. Off by 1.
    # Could it be |Spec(n)| = L(2n-6) + epsilon?
    for n in range(3, 8):
        idx = 2*n - 6
        if idx >= 0:
            l_val = lucas[idx]
            spec = h_spec[n]
            print(f"    n={n}: L({idx})={l_val}, |Spec|={spec}, diff={spec-l_val}")

    # n=3: L(0)=2, Spec=2, diff=0 MATCH!
    # n=4: L(2)=3, Spec=3, diff=0 MATCH!
    # n=5: L(4)=7, Spec=7, diff=0 MATCH!
    # n=6: L(6)=18, Spec=19, diff=1
    # n=7: L(8)=47, Spec=77, diff=30

    # So for n=3,4,5: |Spec(n)| = L(2n-6) exactly!
    # n=6: off by 1. n=7: way off.
    # Not the right formula.

    # What if |Spec(n)| has a recursive structure based on the forbidden values?
    # |Spec(n)| = (max-1)/2 + 1 - |forbidden interior|
    # max at n=5: 15. (15-1)/2+1 = 8. 8-1(forbidden) = 7. YES!
    # max at n=6: 45. (45-1)/2+1 = 23. 23-4(forbidden at n=6) = 19. YES!
    # max at n=7: 189. (189-1)/2+1 = 95. 95-18(forbidden) = 77. YES!

    print(f"\n  |Spec(n)| = (max_H - 1)/2 + 1 - |forbidden|:")
    data = {
        3: (3, 0),    # max=3, 0 forbidden
        4: (5, 0),    # max=5, 0 forbidden
        5: (15, 1),   # max=15, 1 forbidden (7)
        6: (45, 4),   # max=45, 4 forbidden (7,21,35,39)
        7: (189, 18),  # max=189, 18 forbidden
    }
    for n in sorted(data.keys()):
        max_h, n_forb = data[n]
        total_odd = (max_h - 1) // 2 + 1  # odd numbers in [1, max_h]
        predicted = total_odd - n_forb
        actual = h_spec[n]
        print(f"    n={n}: max={max_h}, odds={total_odd}, forb={n_forb}, pred={predicted}, actual={actual}, match={predicted==actual}")

    # This is TAUTOLOGICAL (|Spec| = total odd - forbidden). But it tells us:
    # The real question is: what determines max_H and |forbidden|?

    # max_H at each n:
    print(f"\n  max(H) sequence: 1, 1, 3, 5, 15, 45, 189")
    maxes = [1, 1, 3, 5, 15, 45, 189]
    for i in range(1, len(maxes)):
        r = maxes[i] / maxes[i-1] if maxes[i-1] > 0 else 0
        print(f"    max({i+1})/max({i}) = {maxes[i]}/{maxes[i-1]} = {r:.4f}")

    # 3/1=3, 5/3=1.67, 15/5=3, 45/15=3, 189/45=4.2
    # The ratios: 3, 5/3, 3, 3, 4.2
    # For n >= 4: max(H) = max(n-1) * (n-1)?
    # 5*3 = 15. YES! (n=5: max = 5*3 = 15)
    # 15*3 = 45. YES! (n=6: max = 15*3 = 45)
    # 45*4.2 = 189. Hmm, 189/45 = 4.2 = 21/5. Not clean.

    # Actually: max H at n is known. For transitive tournament, H=1.
    # For the tournament with max H (all rotational/quadratic residue):
    # max H = (n-1)!! for odd n?
    # (n-1)!! for n=3: 2!!=2, but max=3. Nope.

    # From the data: 3, 5, 15, 45, 189
    # 3 = 3, 5 = 5, 15 = 3*5, 45 = 9*5 = 3*15, 189 = ?
    # 189 = 27*7 = 3^3 * 7. And 189/45 = 4.2 = 21/5.

    # Let me check: max H for n=3..7 vs n!/2^{n-1} (the mean)
    from math import factorial
    print(f"\n  max(H)/Mean(H):")
    means = {3: 3/2, 4: 3, 5: 15/2, 6: 45/2, 7: 315/4}
    for n in range(3, 8):
        r = maxes[n-1] / means[n]
        print(f"    n={n}: max/mean = {maxes[n-1]}/{means[n]} = {r:.4f}")
    # n=3: 3/1.5=2, n=4: 5/3=1.67, n=5: 15/7.5=2, n=6: 45/22.5=2, n=7: 189/78.75=2.4
    # Interesting: max/mean = 2 for n=3,5,6. And 2.4 for n=7.
    # max(H) = 2*Mean for even max (n=3: 3 = 2*1.5, n=5: 15 = 2*7.5, n=6: 45 = 2*22.5)

    print(f"\n  OBSERVATION: max(H) = 2*Mean for n=3,5,6 (exact)")
    print(f"  At n=7: max = 189 = 2.4*Mean, which is > 2*Mean = 157.5")

    print(f"\n{'='*70}")
    print("DONE -- EISENSTEIN SPECTRUM CONNECTION")
    print("="*70)


if __name__ == "__main__":
    main()
