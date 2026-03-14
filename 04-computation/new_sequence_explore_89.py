#!/usr/bin/env python3
"""
new_sequence_explore_89.py -- opus-2026-03-14-S89

The sequence a(n) = V(n) * n! / 2 where V(n) = Var(H)/Mean(H)^2
is NOT in OEIS. Explore its structure deeply.

a(n) = sum_{k=1}^{floor((n-1)/2)} (n-2k)^k * n! / P(n,2k)
     = sum_{k=1}^{floor((n-1)/2)} (n-2k)^k * (n-2k)!

Wait: n!/P(n,2k) = n!/[n!/(n-2k)!] = (n-2k)!
So a(n) = sum_{k=1}^{K} (n-2k)^k * (n-2k)!

This is BEAUTIFUL! The n!/2 scaling reveals:
  a(n) = sum_{k=1}^{K} (n-2k)^k * (n-2k)!

Each term is j^k * j! where j = n-2k, k ranges from 1 to K.
"""

from fractions import Fraction
from math import factorial, comb


def var_ratio(n):
    """Compute exact Var/Mean^2 using the Grand Theorem."""
    K = (n - 1) // 2
    total = Fraction(0)
    for k in range(1, K + 1):
        pn2k = 1
        for i in range(2 * k):
            pn2k *= (n - i)
        term = Fraction(2 * (n - 2*k)**k, pn2k)
        total += term
    return total


def a_n(n):
    """Compute a(n) = V(n) * n! / 2 = sum_{k=1}^K (n-2k)^k * (n-2k)!"""
    K = (n - 1) // 2
    total = 0
    for k in range(1, K + 1):
        j = n - 2 * k
        total += j**k * factorial(j)
    return total


def main():
    print("="*70)
    print("NEW INTEGER SEQUENCE: a(n) = V(n) * n! / 2")
    print("opus-2026-03-14-S89")
    print("="*70)

    # Verify the simplification
    print("\nPart 1: Verify a(n) = sum (n-2k)^k * (n-2k)!")
    print("="*70)
    for n in range(3, 15):
        vr = var_ratio(n)
        from_vr = vr * factorial(n) // 2
        from_formula = a_n(n)
        print(f"  n={n}: V*n!/2 = {from_vr}, formula = {from_formula}, match = {from_vr == from_formula}")

    # Display the sequence
    print("\nPart 2: The sequence a(n) for n = 3..30")
    print("="*70)
    for n in range(3, 31):
        val = a_n(n)
        print(f"  a({n:2d}) = {val}")

    # Term-by-term decomposition
    print("\nPart 3: Term-by-term: (n-2k)^k * (n-2k)!")
    print("="*70)
    for n in [5, 7, 9, 11, 13]:
        K = (n - 1) // 2
        print(f"\n  n={n} (K={K}):")
        for k in range(1, K + 1):
            j = n - 2 * k
            term = j**k * factorial(j)
            print(f"    k={k}: j={j}, j^k={j**k}, j!={factorial(j)}, term = {term}")
        print(f"    Total a({n}) = {a_n(n)}")

    # The structure j^k * j! = Gamma(j+1) * j^k
    # This looks like a moment: integral of x^k over Poisson or similar?

    # Related: sum_{j} j^k * j! is like evaluating a generating function
    # at specific points. Let's see the EGF.

    print("\nPart 4: EGF / OGF analysis")
    print("="*70)

    # a(n) = sum over k: (n-2k)^k * (n-2k)!
    # For odd n=2m+1: k goes from 1 to m, j = 2m+1-2k goes from 2m-1 down to 1
    #   a(2m+1) = sum_{k=1}^m (2m+1-2k)^k * (2m+1-2k)!
    #           = sum_{j=1,3,5,...,2m-1} j^{(2m+1-j)/2} * j!
    # For even n=2m: k goes from 1 to m-1, j = 2m-2k goes from 2m-2 down to 2
    #   a(2m) = sum_{k=1}^{m-1} (2m-2k)^k * (2m-2k)!
    #         = sum_{j=2,4,...,2m-2} j^{(2m-j)/2} * j!

    print("\n  For odd n=2m+1: a(n) = sum over ODD j in [1, n-2] of j^{(n-j)/2} * j!")
    print("  For even n=2m: a(n) = sum over EVEN j in [2, n-2] of j^{(n-j)/2} * j!")

    # This is remarkable: odd n sums over odd j, even n over even j!
    # The index j has the same parity as n.

    # Let's verify
    for n in range(3, 12):
        K = (n - 1) // 2
        terms = []
        for k in range(1, K + 1):
            j = n - 2 * k
            terms.append((j, k, j**k * factorial(j)))
        js = [t[0] for t in terms]
        parities = [j % 2 for j in js]
        n_parity = n % 2
        print(f"  n={n} (parity {n_parity}): j values = {js}, parities = {parities}")

    # Ratio a(n+1)/a(n)
    print("\nPart 5: Growth rate a(n+1)/a(n)")
    print("="*70)
    prev = a_n(3)
    for n in range(4, 25):
        val = a_n(n)
        ratio = val / prev
        # Compare with n
        print(f"  a({n})/a({n-1}) = {ratio:.6f}, n-1={n-1}, (n-1)*ratio/n = {(n-1)*ratio/n:.6f}")
        prev = val

    # The dominant term for large n is k=1: (n-2)^1 * (n-2)! = (n-2) * (n-2)! = (n-1)!
    # Wait: (n-2) * (n-2)! = (n-2)!*(n-2). Hmm, (n-1)! = (n-1)*(n-2)!, so (n-2)*(n-2)! < (n-1)!
    # Actually the dominant term IS (n-2)*(n-2)! for the k=1 term.
    # Next: k=2 gives (n-4)^2 * (n-4)!, which is much smaller.

    print("\nPart 6: Dominant term analysis")
    print("="*70)
    for n in range(3, 20):
        val = a_n(n)
        k1_term = (n-2) * factorial(n-2)  # k=1 term
        ratio = val / k1_term if k1_term > 0 else float('inf')
        print(f"  n={n}: a(n)={val}, (n-2)*(n-2)!={k1_term}, a(n)/dom = {ratio:.8f}")

    # So a(n) ~ (n-2)*(n-2)! * (1 + correction)
    # The correction = sum_{k>=2} (n-2k)^k * (n-2k)! / [(n-2)*(n-2)!]
    # For k=2: (n-4)^2 * (n-4)! / [(n-2)*(n-2)!] = (n-4)^2 / [(n-2)*(n-3)*(n-2)]
    #        = (n-4)^2 / [(n-2)^2 * (n-3)] ~ 1/n for large n.
    # So a(n) ~ (n-2)*(n-2)! * (1 + O(1/n))
    # And (n-2)*(n-2)! = (n-1)! - (n-2)! ... no, it's just (n-2)*(n-2)!.
    # Note: (n-2)! * (n-2) = (n-2)! * (n-2), not a standard factorial.

    # OEIS search: b(n) = (n-2)*(n-2)! = 1, 4, 18, 96, 600, 4320, 35280, ...
    # vs a(n) =                           1, 4, 19, 104, 655, 4720, 38443, ...
    # Difference: 0, 0, 1, 8, 55, 400, 3163, ...

    print("\nPart 7: Correction terms a(n) - (n-2)*(n-2)!")
    print("="*70)
    corrections = []
    for n in range(3, 20):
        val = a_n(n)
        dom = (n-2) * factorial(n-2)
        diff = val - dom
        corrections.append(diff)
        print(f"  n={n}: a(n)-dom = {diff}")

    # corrections: 0, 0, 1, 8, 55, 400, 3163, 27648, ...
    # Is THIS in OEIS? Let's check the ratios.
    print("\n  Correction ratios:")
    for i in range(2, len(corrections)):
        if corrections[i-1] > 0:
            r = corrections[i] / corrections[i-1]
            print(f"    c({i+3})/c({i+2}) = {r:.4f}")

    # The corrections are a(n) - (n-2)*(n-2)! = sum_{k>=2} (n-2k)^k * (n-2k)!
    # The k=2 term is (n-4)^2 * (n-4)!
    # For n=5: 1^2 * 1! = 1. Check: correction at n=5 is 1. YES!
    # For n=6: 2^2 * 2! = 8. Check: correction at n=6 is 8. YES!
    # For n=7: 3^2 * 3! + 1^3 * 1! = 54 + 1 = 55. Check: 55. YES!
    # For n=8: 4^2 * 4! + 2^3 * 2! = 384 + 16 = 400. Check: 400. YES!

    print("\n  Verification: correction = sum_{k>=2} (n-2k)^k * (n-2k)!")
    for n in range(5, 15):
        K = (n - 1) // 2
        correction = sum((n-2*k)**k * factorial(n-2*k) for k in range(2, K+1))
        expected = a_n(n) - (n-2) * factorial(n-2)
        print(f"    n={n}: correction = {correction}, expected = {expected}, match = {correction == expected}")

    # Part 8: Recurrence?
    print("\nPart 8: Looking for a recurrence")
    print("="*70)

    # Try a(n) = alpha*n * a(n-1) + beta * a(n-2) + ...
    # a(3)=1, a(4)=4, a(5)=19, a(6)=104, a(7)=655
    # If a(n) = c1*n*a(n-1) + c2*a(n-2):
    #   19 = 5*c1*4 + c2*1 = 20c1 + c2
    #   104 = 6*c1*19 + c2*4 = 114c1 + 4c2
    # From first: c2 = 19 - 20c1
    # Substitute: 104 = 114c1 + 4(19-20c1) = 114c1 + 76 - 80c1 = 34c1 + 76
    # 34c1 = 28, c1 = 14/17... not clean.

    # Try a(n) = A*(n-1)*a(n-1) + B*a(n-2):
    #   19 = 3A*4 + B*1 = 12A + B
    #   104 = 4A*19 + B*4 = 76A + 4B
    # From first: B = 19-12A
    # Sub: 104 = 76A + 4(19-12A) = 76A + 76 - 48A = 28A + 76
    # 28A = 28, A = 1
    # B = 19 - 12 = 7

    # Test: a(n) = (n-1)*a(n-1) + 7*a(n-2)?
    print("  Testing a(n) = (n-1)*a(n-1) + 7*a(n-2):")
    seq = [0, 0, 0, 1, 4]
    for n in range(5, 20):
        pred = (n-1) * seq[n-1] + 7 * seq[n-2]
        actual = a_n(n)
        seq.append(actual)
        print(f"    n={n}: predicted={pred}, actual={actual}, match={pred==actual}")

    # That won't work for all n. Let me try more general:
    # a(n) = f(n)*a(n-1) + g(n)*a(n-2)
    print("\n  Computing f(n) = [a(n) - g*a(n-2)] / a(n-1) for various g:")
    vals = {n: a_n(n) for n in range(3, 20)}

    # Try to find f(n), g(n) such that a(n) = f(n)*a(n-1) + g(n)*a(n-2)
    # This gives: f(n) = [a(n) - g(n)*a(n-2)] / a(n-1)
    # Two equations for n, n+1 with two unknowns f(n), f(n+1), g(n), g(n+1)
    # Too many unknowns unless f and g have simple forms.

    # Let's just compute the ratios directly:
    print("\n  Direct ratio a(n)/a(n-1):")
    for n in range(4, 20):
        r = Fraction(vals[n], vals[n-1])
        # Try to express as (An+B)/C
        print(f"    a({n})/a({n-1}) = {r} = {float(r):.6f}")

    # Part 9: Connection to subfactorials / derangements?
    print("\nPart 9: Connection to derangements and subfactorials")
    print("="*70)

    # D(n) = n! * sum_{k=0}^n (-1)^k/k! (derangements)
    def derangement(n):
        return sum((-1)**k * factorial(n) // factorial(k) for k in range(n+1))

    for n in range(3, 15):
        dn = derangement(n-2)
        ratio = Fraction(vals[n], dn) if dn > 0 else None
        print(f"  n={n}: a(n)={vals[n]}, D(n-2)={dn}, ratio={ratio}")

    # Part 10: EGF coefficients
    print("\nPart 10: EGF coefficients a(n)/n!")
    print("="*70)
    for n in range(3, 20):
        coeff = Fraction(vals[n], factorial(n))
        print(f"  a({n})/({n}!) = {coeff} = {float(coeff):.10f}")

    # Note: a(n)/n! = V(n)/2 = Var(H)/(2*Mean(H)^2)
    # This is half the normalized variance.

    print(f"\n{'='*70}")
    print("DONE -- NEW SEQUENCE EXPLORATION")
    print("="*70)


if __name__ == "__main__":
    main()
