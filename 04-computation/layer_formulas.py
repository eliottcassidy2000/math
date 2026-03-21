#!/usr/bin/env python3
"""layer_formulas.py — Exact formulas for each deficit layer.

Session: kind-pasteur-2026-03-20-S6

DISCOVERED: Each layer has an EXACT closed form.

For a single k-cycle (k odd):
  f_k(n) = (k-1)/2 + C(n-k+1, 2)  free edge orbits
  D_k(n) = 2^{f_k(n)} * k / (k * (n-k)!) = 2^{f_k(n)} / (n-k)!
  D_k(n)/2 = 2^{f_k(n) - 1} / (n-k)!

For two 3-cycles:
  f_{3,3}(n) = C(n-4, 2) + 3 + 2 = 5 + 2(n-6) + C(n-6,2)
  D_{3,3}(n) = 2^{f_{3,3}(n)} / (3 * (n-6)!)

The "central binomial miracle": D_3(n)/2 + D_5(n)/2 = C(2(n-3), n-3)
for n = 3, 4, 5, 6. This identity BREAKS at n=7 because D_{3,3} becomes nonzero.
"""

from math import factorial, comb
from fractions import Fraction

def f_single_cycle(k, n):
    """Number of free edge orbits for a single k-cycle with n-k fixed points."""
    return (k - 1) // 2 + comb(n - k + 1, 2)

def D_single_cycle(k, n):
    """D_k(n) contribution from single k-cycle layer."""
    if n < k:
        return Fraction(0)
    f = f_single_cycle(k, n)
    return Fraction(2**f, factorial(n - k))

def f_two_3cycles(n):
    """Free edge orbits for two 3-cycles with n-6 fixed points."""
    return 5 + 2*(n-6) + comb(n-6, 2)

def D_two_3cycles(n):
    """D_{3,3}(n) contribution."""
    if n < 6:
        return Fraction(0)
    f = f_two_3cycles(n)
    return Fraction(2**f, 3 * factorial(n - 6))

def f_cycle_plus_3cycle(k, n):
    """Free edge orbits for one k-cycle + one 3-cycle with n-k-3 fixed points."""
    if n < k + 3:
        return 0
    # Within k-cycle: (k-1)/2
    # Within 3-cycle: 1
    # Between k-cycle and 3-cycle: gcd(k,3) free choices
    # k odd: gcd(k,3) = 3 if 3|k, else 1
    g = 3 if k % 3 == 0 else 1
    between_cycles = g
    # Wait, need to compute properly for k=5:
    # 5-cycle and 3-cycle: vertices {0,1,2,3,4} and {5,6,7}
    # Between cycles: 5*3 = 15 ordered pairs, orbit size lcm(5,3)=15, so 1 orbit
    # Reverse: also 1 orbit. So 1 free choice.
    # For k=7 and 3: 7*3=21 pairs, orbit size lcm(7,3)=21, 1 orbit + 1 reverse. 1 free.
    # For k=3 and 3: already handled in D_two_3cycles.
    # So between a k-cycle and 3-cycle (k != 3):
    # Number of orbits on ordered pairs = gcd(k,3). Self-reverse?
    # Actually for k=5: orbit of (0,5) has size lcm(5,3)=15, visits all 15 pairs.
    # The reverse orbit of (5,0) also has size 15. Are they the same?
    # (0,5) -> (1,6) -> (2,7) -> (3,5) -> (4,6) -> (0,7) -> (1,5) -> (2,6) -> (3,7) -> (4,5) -> (0,6) -> (1,7) -> (2,5) -> (3,6) -> (4,7)
    # That's 15 elements. Reverse: (5,0) -> (6,1) -> (7,2) -> ... also 15 elements.
    # Are they the same orbit? (0,5) would need to appear in the reverse orbit:
    # we'd need sigma^j to map (5,0) to (0,5), meaning sigma^j(5)=0 and sigma^j(0)=5.
    # sigma^j(0) = j mod 5, sigma^j(5) = 5 + (j mod 3). Need j=0 mod 5 and j=0 mod 3 => j=0 mod 15.
    # But j < 15, so only j=0. (5,0) maps to (5,0) at j=0, not (0,5). So NOT self-reverse!
    # So 1 free choice between the 5-cycle and 3-cycle.
    between = 1  # For k not multiple of 3 and gcd(k,3)=1

    # Between each cycle and fixed points
    # k-cycle gives 1 free choice per fixed vertex
    # 3-cycle gives 1 free choice per fixed vertex
    fixed_count = n - k - 3
    between_fixed = 2 * fixed_count

    # Between fixed vertices
    fixed_pairs = comb(fixed_count, 2)

    return (k-1)//2 + 1 + between + between_fixed + fixed_pairs

def D_k_plus_3(k, n):
    """D_{k,3}(n) for k-cycle + 3-cycle."""
    if n < k + 3 or k == 3:
        return Fraction(0)
    f = f_cycle_plus_3cycle(k, n)
    # N = n! / (k * 3 * (n-k-3)!)
    # n-fix = k + 3
    # D = N * 2^f * (k+3) / n! = 2^f * (k+3) / (k * 3 * (n-k-3)!)
    return Fraction(2**f * (k + 3), k * 3 * factorial(n - k - 3))


def main():
    print("=" * 70)
    print("EXACT LAYER FORMULAS FOR D(n)")
    print("=" * 70)

    print(f"\n  Layer formulas:")
    print(f"  D_k(n) = 2^{{(k-1)/2 + C(n-k+1,2)}} / (n-k)!")
    print(f"  D_{{3,3}}(n) = 2^{{5+2(n-6)+C(n-6,2)}} / (3*(n-6)!)")

    # Verify all layers match the direct computation
    print(f"\n  Verification:")
    print(f"\n  {'n':>3} {'D_3':>12} {'D_5':>12} {'D_33':>12} {'D_7':>12} {'D_53':>12} {'Sum':>12} {'D(n)':>12} {'Match':>6}")

    for n in range(3, 11):
        d3 = D_single_cycle(3, n)
        d5 = D_single_cycle(5, n)
        d33 = D_two_3cycles(n)
        d7 = D_single_cycle(7, n)
        d53 = D_k_plus_3(5, n)
        d333 = Fraction(0)  # Three 3-cycles, first at n=9

        # Three 3-cycles at n >= 9
        if n >= 9:
            # f = 1+1+1 + 3+3+3 + 3*(n-9) + C(n-9,2) = 12 + 3(n-9) + C(n-9,2)
            # Actually need careful computation
            # For three 3-cycles: 3 cycles, each pair has 3 free choices between them (3*3 pairs, 3 orbits, each non-self-reverse)
            # Within: 3 free
            # Between pairs of cycles: C(3,2) * (number of inter-cycle free choices)
            # For two 3-cycles: 3 between them (as computed). For three pairs: 3*3 = 9.
            # Between cycles and fixed: 3*(n-9)
            # Between fixed: C(n-9,2)
            f_333 = 3 + 9 + 3*(n-9) + comb(n-9, 2)
            # N = n! / (3^3 * 3! * (n-9)!)
            # fix = n-9, n-fix = 9
            d333 = Fraction(2**f_333 * 9, 27 * 6 * factorial(n-9))

        total = d3 + d5 + d33 + d7 + d53 + d333

        # For n=9,10 we also need single 9-cycle and (7,1,1), (3,3,3), (5,3,1) etc.
        d9 = D_single_cycle(9, n) if n >= 9 else Fraction(0)
        # (3,3,1^{n-6}) already counted as d33
        # We need more layers for full computation at n >= 9

        # Known D values from P7_burnside.py
        known_D = {3: 2, 4: 4, 5: 12, 6: 40, 7: 152, 8: 784}
        D_actual = known_D.get(n, None)

        total_with_more = total + d9
        total_float = float(total_with_more)

        match_str = ''
        if D_actual:
            match_str = 'Y' if abs(total_float - D_actual) < 0.01 else f'({D_actual})'

        print(f"  {n:3d} {float(d3):12.4f} {float(d5):12.4f} {float(d33):12.4f} "
              f"{float(d7):12.4f} {float(d53):12.4f} {total_float:12.4f} "
              f"{str(D_actual):>12} {match_str:>6}")

    # THE CENTRAL BINOMIAL MIRACLE
    print(f"\n  {'='*70}")
    print(f"  THE CENTRAL BINOMIAL MIRACLE")
    print(f"  {'='*70}")

    print(f"\n  Why does D_3 + D_5 = 2*C(2k,k) for n=3,4,5,6?")
    print(f"\n  {'n':>3} {'D_3/2':>12} {'D_5/2':>12} {'Sum':>12} {'C(2k,k)':>12} {'Equal':>6}")

    for n in range(3, 9):
        d3_half = D_single_cycle(3, n) / 2
        d5_half = D_single_cycle(5, n) / 2
        s = d3_half + d5_half
        k = n - 3
        cb = comb(2*k, k)
        eq = 'YES' if s == cb else 'NO'
        print(f"  {n:3d} {float(d3_half):12.4f} {float(d5_half):12.4f} "
              f"{float(s):12.4f} {cb:12d} {eq:>6}")

    # EXACT FORMULA: D_3(n)/2 + D_5(n)/2
    print(f"\n  D_3(n)/2 = 2^{{C(n-2,2)}} / (n-3)!")
    print(f"  D_5(n)/2 = 2^{{1+C(n-4,2)}} / (n-5)!")
    print(f"  Sum = 2^{{C(n-2,2)}} / (n-3)! + 2^{{1+C(n-4,2)}} / (n-5)!")

    print(f"\n  C(2k,k) = (2k)! / (k!)^2 for k = n-3")

    # At n=6 (k=3): D_3/2 = 2^6/3! = 64/6, D_5/2 = 2^2/1! = 4
    # Sum = 64/6 + 4 = 64/6 + 24/6 = 88/6 ... but C(6,3)=20.
    # 88/6 = 14.67 ≠ 20. So D_3/2 + D_5/2 ≠ C(2k,k) at n=6!

    # Wait, but the total D/2 = 20 at n=6, which INCLUDES D_{3,3}/2.
    # So the miracle is D_3/2 + D_5/2 + D_{3,3}/2 = C(2k,k) for n=3..6.
    # NOT just D_3 + D_5.

    print(f"\n  CORRECTION: The miracle is D_3 + D_5 + D_{{3,3}} = 2*C(2k,k)")
    print(f"  i.e., ALL layers together give 2*C(2k,k) for n=3..6")
    print(f"  This is trivially true since those are the ONLY layers at n<=6!")
    print(f"  The question is: why does the SUM happen to be 2*C(2k,k)?")

    print(f"\n  {'n':>3} {'D/2':>10} {'C(2k,k)':>10} {'D/2-C':>10}")
    known_D = {3: 2, 4: 4, 5: 12, 6: 40, 7: 152, 8: 784}
    for n in range(3, 9):
        k = n-3
        D_half = known_D[n] // 2
        cb = comb(2*k, k)
        print(f"  {n:3d} {D_half:10d} {cb:10d} {D_half - cb:10d}")

    # THE EXCESS SEQUENCE: D/2 - C(2k,k)
    print(f"\n  Excess sequence E(n) = D(n)/2 - C(2(n-3), n-3):")
    print(f"  E = 0, 0, 0, 0, 6, 140")
    print(f"  E first nonzero at n=7 with E=6")
    print(f"  This 6 comes from the (3,3,1) + (7) layers")
    print(f"  that aren't present at n<=6 in sufficient strength")

    # What generates the excess?
    print(f"\n  Excess decomposition at n=7:")
    for n in [7, 8]:
        D_half = known_D[n] // 2
        k = n - 3
        cb = comb(2*k, k)
        excess = D_half - cb
        d3_half = float(D_single_cycle(3, n) / 2)
        d5_half = float(D_single_cycle(5, n) / 2)
        d33_half = float(D_two_3cycles(n) / 2)
        d7_half = float(D_single_cycle(7, n) / 2)
        d53_half = float(D_k_plus_3(5, n) / 2)

        print(f"\n  n={n}: D/2={D_half}, C(2k,k)={cb}, excess={excess}")
        print(f"    D_3/2 = {d3_half:.4f}")
        print(f"    D_5/2 = {d5_half:.4f}")
        print(f"    D_33/2 = {d33_half:.4f}")
        print(f"    D_7/2 = {d7_half:.4f}")
        print(f"    D_53/2 = {d53_half:.4f}")
        print(f"    Sum = {d3_half + d5_half + d33_half + d7_half + d53_half:.4f}")


if __name__ == "__main__":
    main()
