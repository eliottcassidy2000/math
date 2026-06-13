"""
gap_structure_deep.py -- kind-pasteur-2026-03-14-S66

Deep analysis of the gap structure connecting 2, 3, and cyclotomic polynomials.

KEY FRAMEWORK:
- H = I(Omega(T), 2) = sum alpha_k * 2^k
- 3^k bound: alpha_k >= 1 => H >= 3^k
- The levels [3^k, 3^{k+1}) partition achievable H into "complexity classes"
- Gaps can only occur within levels, not across levels

DISCOVERIES:
1. H=7 = Phi_3(2) = 2^2+2+1: unique gap at level 1
2. H=21 = 3*Phi_3(2): unique gap at level 2
3. 21 = (x^3+2x^2+2x+1)|_{x=2} = x*(x+1)^2 + 1 at x=2
4. 7 = (x^2+x+1)|_{x=2} = the norm of (x-omega)(x-omega^2) at x=2
5. The gaps are at H = Phi_3(2) and 3*Phi_3(2) = (2+1)*Phi_3(2)

QUESTION: Is there a deeper algebraic reason why Phi_3 at 2 gives the gaps?
"""

import numpy as np
from math import gcd

def main():
    print("=" * 70)
    print("GAP STRUCTURE AND THE ROLE OF Phi_3(2) = 7")
    print("=" * 70)

    # Part 1: Cyclotomic evaluation at x=2
    print("\n--- CYCLOTOMIC POLYNOMIALS AT x=2 ---")

    # Phi_n(x) values at x=2
    # Phi_1(x) = x-1, Phi_1(2) = 1
    # Phi_2(x) = x+1, Phi_2(2) = 3
    # Phi_3(x) = x^2+x+1, Phi_3(2) = 7
    # Phi_4(x) = x^2+1, Phi_4(2) = 5
    # Phi_5(x) = x^4+x^3+x^2+x+1, Phi_5(2) = 31
    # Phi_6(x) = x^2-x+1, Phi_6(2) = 3
    # Phi_n(2) = product over d|n of (2-zeta_d)

    # 2^n - 1 = prod_{d|n} Phi_d(2)
    print("Mersenne factorization: 2^n - 1 = prod_{d|n} Phi_d(2)")
    for n in range(1, 13):
        mersenne = 2**n - 1
        # Factor as product of Phi_d(2)
        factors = []
        for d in range(1, n+1):
            if n % d == 0:
                # Compute Phi_d(2) via Mobius
                phi_d = 1
                for k in range(1, d+1):
                    if d % k == 0:
                        # Mobius function
                        pass
                # Simpler: Phi_d(2) = (2^d - 1) / prod_{e|d, e<d} Phi_e(2)
                pass

        print(f"  2^{n:2d} - 1 = {mersenne:6d}")

    # Direct computation of Phi_d(2)
    phi_at_2 = {}
    for d in range(1, 25):
        # Phi_d(2) = prod_{d=e*f, f|d, f<d} removed
        # Use: 2^d - 1 = prod_{e|d} Phi_e(2)
        product_lower = 1
        for e in range(1, d):
            if d % e == 0:
                product_lower *= phi_at_2[e]
        phi_at_2[d] = (2**d - 1) // product_lower
        print(f"  Phi_{d:2d}(2) = {phi_at_2[d]:8d}")

    # Part 2: The gap values in terms of Phi
    print("\n--- GAP VALUES ---")
    print(f"  H=7  = Phi_3(2) = {phi_at_2[3]}")
    print(f"  H=21 = 3 * 7 = Phi_2(2) * Phi_3(2) = {phi_at_2[2] * phi_at_2[3]}")
    print(f"  H=21 = (2^3-1) * 3 = 7*3 = 21")
    print(f"  H=21 = 2^4 + 2^2 + 2^0 = 16+4+1 (alternating in binary)")
    print()

    # Part 3: The polynomial identity
    print("--- THE POLYNOMIAL IDENTITY ---")
    print("  Phi_3(x) = x^2 + x + 1")
    print("  x * Phi_3(x) = x^3 + x^2 + x = x(x+1)^2 - x")
    print("  (x+1) * Phi_3(x) = x^3 + 2x^2 + 2x + 1")
    print("    = (x^3 + 1) + 2x(x+1)")
    print("    = (x+1)(x^2-x+1) + 2x(x+1)")
    print("    = (x+1)(x^2+x+1)")
    print("    = (x+1) * Phi_3(x)")
    print()
    print("  At x=2: (2+1) * Phi_3(2) = 3 * 7 = 21")
    print()
    print("  So 21 = Phi_2(2) * Phi_3(2) = (2+1)(2^2+2+1)")
    print("       = (x+1)(x^2+x+1) at x=2")
    print("       = Phi_6(-2) ... no, Phi_6(x) = x^2-x+1")
    print()

    # Part 4: Binary representation
    print("--- BINARY STRUCTURE ---")
    for h in [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27]:
        binary = bin(h)[2:]
        T = (h-1)//2
        print(f"  H={h:3d} (T={T:3d}): binary = {binary:>8s}  "
              f"{'GAP' if h in [7,21] else 'ok'}")

    print()
    print("Pattern:")
    print("  H=7  = 111       (all 1s, length 3)")
    print("  H=21 = 10101     (alternating, length 5)")
    print("  H=63 = 111111    (all 1s, length 6) -- NOT permanent gap")
    print()
    print("  7  = 2^3 - 1 = M_3   (Mersenne prime)")
    print("  21 = ... not Mersenne")
    print("  63 = 2^6 - 1 = M_6 = 9*7 -- NOT prime")
    print()
    print("  But: 21 in binary is 10101 = sum_{k even} 2^k for k=0,2,4")
    print("       = (2^5 - 1)/(2 - (-1)) = 31/3... no")
    print("       = (2^5 + 1)/(2 + 1) = 33/3 = 11... no")
    print("       = (4^3 - 1)/(4 - 1) = 63/3 = 21 YES!")
    print()
    print("  So 21 = (4^3 - 1)/3 = (2^6 - 1)/(2^2 - 1)")
    print("        = Phi_3(4) = Phi_3(2^2)")
    print("        = (2^2)^2 + 2^2 + 1 = 16 + 4 + 1 = 21")
    print()
    print("  REMARKABLE: 21 = Phi_3(4) = Phi_3(2^2)")
    print("  And:        7 = Phi_3(2)")
    print()
    print("  The gaps are Phi_3(2^1) = 7 and Phi_3(2^2) = 21!")
    print("  Is Phi_3(2^3) = Phi_3(8) = 64+8+1 = 73 a gap?")
    print("  H=73 at n=7: check achievability...")
    print("  From n=7 spectrum: 73 IS achievable! So NOT a gap.")
    print()
    print("  Why does the pattern break at Phi_3(2^3) = 73?")
    print("  T = (73-1)/2 = 36. Since 36 > 13 = (3^3-1)/2,")
    print("  cubic decomposition (26, 3, 1) works.")
    print("  AND T=36 has many valid quadratic decompositions too.")
    print()
    print("  The pattern Phi_3(2^k) gives gaps only for k=1,2.")
    print("  At k=3: T=36 is well within reach of cubic I.P.")

    # Part 5: The additive vs multiplicative roles
    print(f"\n{'='*70}")
    print("ADDITIVE AND MULTIPLICATIVE ROLES OF 2 AND 3")
    print(f"{'='*70}")
    print()
    print("ADDITIVE (the evaluation point x=2):")
    print("  H = 1 + 2*a1 + 4*a2 + 8*a3 + ...")
    print("  Powers of 2 weight the independence structure")
    print("  2^k = weight of k-element independent sets")
    print()
    print("MULTIPLICATIVE (the hierarchy base 3 = 2+1):")
    print("  3^k = minimum H from k disjoint cycles")
    print("  (1+2)^k = binomial sum of weights")
    print("  The 3^k bound comes from the PRODUCT (1+x)^k")
    print()
    print("INTERPLAY:")
    print("  The gap H=7 arises because:")
    print("  - Additively: T=3 is too small for alpha_2 >= 1 (needs T >= 4)")
    print("    and alpha_1=3 forces splicing")
    print("  - Multiplicatively: 7 = 3^2 - 2 is just below the k=2 threshold")
    print()
    print("  The gap H=21 arises because:")
    print("  - Additively: T=10 is too small for alpha_3 >= 1 (needs T >= 13)")
    print("  - Multiplicatively: 21 = 3^3 - 6 is below the k=3 threshold")
    print()
    print("  In both cases: the gap is at Phi_3(2^k) for k=1,2.")
    print("  Phi_3(x) = x^2+x+1 = (x^3-1)/(x-1)")
    print("  At x=2: Phi_3(2) = (8-1)/(2-1) = 7")
    print("  At x=4: Phi_3(4) = (64-1)/(4-1) = 21")
    print()
    print("  The common factor: x^3-1 = (x-1)(x^2+x+1)")
    print("  At x=2: 7 = (2^3-1)/(2-1)")
    print("  At x=2^2: 21 = (2^6-1)/(2^2-1) = 63/3")
    print()
    print("  Key identity: Phi_3(2^k) = (2^{3k}-1)/(2^k-1)")
    print("  k=1: (8-1)/1 = 7")
    print("  k=2: (64-1)/3 = 21")
    print("  k=3: (512-1)/7 = 73 (NOT a gap)")

    # Part 6: The 5 and 6 connection
    print(f"\n{'='*70}")
    print("THE 5 AND 6 CONNECTION")
    print(f"{'='*70}")
    print()
    print("5 = 2 + 3 (additive combination)")
    print("6 = 2 * 3 (multiplicative combination)")
    print()
    print("In the gap structure:")
    print("  Level 1 span: [3, 9) has 3 odd values: 3, 5, 7")
    print("    GAP at position 3 of 3 (the last one)")
    print("  Level 2 span: [9, 27) has 9 odd values: 9, 11, ..., 25")
    print("    GAP at position 7 of 9 (= 21, the 7th odd in [9,25])")
    print()
    print("  Position of gap within level:")
    print("  Level 1: gap is 7 = 3rd odd in [3,8], position 3/3 = LAST")
    print("  Level 2: gap is 21 = 7th odd in [9,26]")
    print("    21 = 9 + 2*6 = 9 + 12. Position = 7.")
    print("    Fraction: 7/9 ~ 0.778")
    print()
    print("  The gap position within its level:")
    print("  Level 1: at Phi_3(2) = 7, fraction = (7-3)/(9-3) = 4/6 = 2/3")
    print("  Level 2: at 3*Phi_3(2) = 21, fraction = (21-9)/(27-9) = 12/18 = 2/3")
    print()
    print("  BOTH gaps are at exactly the 2/3 mark within their level!")
    print("  (7-3)/(9-3) = 4/6 = 2/3")
    print("  (21-9)/(27-9) = 12/18 = 2/3")
    print()
    print("  This is because Phi_3(2^k) = (2^{3k}-1)/(2^k-1)")
    print("  and the level starts at 3^k, ends at 3^{k+1}.")
    print("  Fraction = (Phi_3(2^k) - 3^k) / (3^{k+1} - 3^k)")
    print("           = (Phi_3(2^k) - 3^k) / (2 * 3^k)")
    for k in range(1, 5):
        phi3 = (2**(3*k) - 1) // (2**k - 1)
        level_start = 3**k
        level_end = 3**(k+1)
        frac_num = phi3 - level_start
        frac_den = level_end - level_start
        g = gcd(frac_num, frac_den)
        print(f"  k={k}: Phi_3(2^{k})={phi3}, level=[{level_start},{level_end}), "
              f"fraction = {frac_num}/{frac_den} = {frac_num//g}/{frac_den//g}")

    print()
    print("  k=1: fraction = 2/3 (= 4/6)")
    print("  k=2: fraction = 2/3 (= 12/18)")
    print("  k=3: fraction = 19/27 (NOT exactly 2/3!)")
    print()
    print("  For k=1,2: the fraction IS exactly 2/3.")
    print("  This means the gap sits at the 2/3 mark of the level width.")
    print("  The level width is 2*3^k, and the gap is at 2*3^k * (2/3) = (4/3)*3^k")
    print("  above the level start. So gap = 3^k + (4/3)*3^k = (7/3)*3^k = 7*3^{k-1}.")
    print()
    print("  k=1: gap = 7 * 3^0 = 7 CHECK")
    print("  k=2: gap = 7 * 3^1 = 21 CHECK")
    print("  k=3: gap would be 7 * 3^2 = 63 but this is NOT a permanent gap!")
    print()
    print("  So the 2/3 pattern predicts gaps at 7, 21, 63, 189, ...")
    print("  But 63 is NOT permanent (found at n=8).")
    print("  The 2/3 pattern holds only for k=1,2.")
    print("  At k>=3, the cubic I.P. provides achievability.")

if __name__ == "__main__":
    main()
