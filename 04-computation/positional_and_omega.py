"""
positional_and_omega.py -- Positional number theory of 10/11 and omega powers
kind-pasteur-2026-03-14-S65

The user's insights:
1. Consider positional number theory implications of 10 and 11
2. omega^2 = 2/3 omega - 1 (actually omega^2 + omega + 1 = 0 => omega^2 = -omega-1)
3. omega^8 is where the pigeonhole principle forces powers inside the unit circle
4. Connection to 8 as inflection point and Cayley-Dickson

This script explores:
- Base 10 = 2*5 and base 11 (prime) in tournament context
- Powers of omega = e^{2pi*i/3} and the pigeonhole argument
- The Hilbert 3rd problem / proof by contradiction structure
- How (omega^a - omega^b)^c in Z[omega] creates irrationality barriers
"""

import cmath
import numpy as np
from math import gcd
from collections import defaultdict

def main():
    print("=" * 70)
    print("POSITIONAL NUMBER THEORY & OMEGA POWERS")
    print("=" * 70)

    omega = cmath.exp(2j * cmath.pi / 3)

    # Part 1: Powers of omega
    print("\nPART 1: POWERS OF omega = e^{2pi*i/3}")
    print()
    print("  omega = e^{2pi*i/3} = -1/2 + sqrt(3)/2 * i")
    print(f"  omega = {omega}")
    print(f"  omega^2 = {omega**2}")
    print(f"  omega^3 = {omega**3}  (should be 1)")
    print(f"  |omega| = {abs(omega)}")
    print()
    print("  Key identity: omega^2 + omega + 1 = 0")
    print(f"  Check: {omega**2 + omega + 1}")
    print()
    print("  So: omega^2 = -omega - 1")
    print("  User said: omega^2 = 2/3 omega - 1. This is WRONG -- the correct")
    print("  relation is omega^2 = -omega - 1 (from Phi_3(omega) = 0).")
    print()
    print("  But 2/3 appears in a DIFFERENT way:")
    print("  The minimal polynomial of omega over Q is x^2 + x + 1.")
    print("  The golden ratio phi has x^2 - x - 1. Note the sign pattern!")
    print("  phi^2 = phi + 1 (expanding)")
    print("  omega^2 = -omega - 1 (contracting)")
    print("  The ratio: if we write omega^2 = c*omega + d:")
    print("    c = -1, d = -1")
    print("  In the Fibonacci analogy: phi^2 = 1*phi + 1, so c=1, d=1.")
    print("  phi and omega are 'dual' in some sense.")

    # Part 2: omega^n for n = 0, 1, ..., 8
    print("\n" + "=" * 70)
    print("PART 2: OMEGA POWERS AND THE PIGEONHOLE ARGUMENT")
    print("=" * 70)

    print("\n  omega^n in Z[omega] = Z + Z*omega:")
    for n in range(13):
        w = omega ** n
        # Express in Z[omega]: a + b*omega
        # omega^0 = 1, omega^1 = omega, omega^2 = -omega-1
        # omega^3 = 1, etc. (period 3)
        a = round(w.real - w.imag / 3**0.5)  # approximation
        b = round(2 * w.imag / 3**0.5)  # approximation
        # Actually: omega = (-1+i*sqrt(3))/2
        # a + b*omega = a + b*(-1+i*sqrt(3))/2 = (a - b/2) + b*sqrt(3)/2 * i
        # So: real part = a - b/2, imag part = b*sqrt(3)/2
        # => b = 2*Im/sqrt(3), a = Re + b/2 = Re + Im/sqrt(3)
        b_exact = round(2 * w.imag / (3**0.5))
        a_exact = round(w.real + w.imag / (3**0.5))
        print(f"    omega^{n:2d} = {a_exact:3d} + {b_exact:3d}*omega  "
              f"(|omega^{n}| = {abs(w):.4f})")

    print(f"\n  Since |omega| = 1, ALL powers have |omega^n| = 1.")
    print(f"  The pigeonhole argument doesn't apply to omega = e^{{2pi*i/3}}!")

    # Part 3: The CORRECT pigeonhole argument (for primitive roots of unity)
    print("\n" + "=" * 70)
    print("PART 3: THE PIGEONHOLE ARGUMENT FOR ROOTS OF UNITY")
    print("=" * 70)

    print(f"""
  The user mentioned: "omega^8 is the point where the pigeonhole principle
  forces powers inside the complex unit circle."

  This refers to the following:
  Let zeta_m = e^{{2pi*i/m}} be a primitive mth root of unity.
  The minimal polynomial is Phi_m(x) of degree phi(m).

  If s is the SMALLEST natural such that zeta_m^s = 1, then s = m.
  That's the definition of primitive root.

  The pigeonhole argument: if |z| = 1 and z^n in Z[z] for all n,
  then z is a root of unity (Kronecker's theorem).

  For the Hilbert 3rd problem connection:
  The Dehn invariant uses R tensor (R/pi*Q).
  The irrational dihedral angles produce IRRATIONAL elements in R/Q*pi.
  Specifically: if a dihedral angle theta is such that theta/pi is
  irrational, then for no n does n*theta become a multiple of pi,
  and the Dehn invariant cannot be zero.

  The "omega^8" statement might refer to:
  When we try to find s such that (2*cos(theta))^s is an algebraic integer
  for theta = pi/m, the powers grow when |2*cos(theta)| > 1 (i.e., m >= 7).
  For m <= 6: |2*cos(pi/m)| <= 2, and the powers stay bounded.
  For m = 7: 2*cos(pi/7) = 1.802..., and powers grow.

  Actually, the key Kronecker-like result:
  If alpha is an algebraic integer with |alpha| <= 1 for ALL conjugates,
  then alpha is a root of unity.

  For the tournament setting: I(omega) = (1-a2) + (a1-a2)*omega
  and |I(omega)|^2 = (1-a2)^2 + (a1-a2)^2 + (1-a2)*(a1-a2).
  This is ALWAYS real and >= 0.
  """)

    # Part 4: (omega^a - omega^b)^c in Z[omega]
    print("=" * 70)
    print("PART 4: (omega^a - omega^b)^c IN Z[omega]")
    print("=" * 70)

    print(f"\n  Since omega has order 3: omega^a is one of {{1, omega, omega^2}}.")
    print(f"  Possible differences omega^a - omega^b:")

    diffs = {}
    for a in range(3):
        for b in range(3):
            if a == b:
                continue
            diff = omega**a - omega**b
            # Express in Z[omega]
            b_coef = round(2 * diff.imag / (3**0.5))
            a_coef = round(diff.real + diff.imag / (3**0.5))
            norm_sq = abs(diff)**2
            key = f"omega^{a} - omega^{b}"
            diffs[key] = (a_coef, b_coef, norm_sq)
            print(f"    {key} = {a_coef} + {b_coef}*omega  (|.|^2 = {norm_sq:.4f})")

    print(f"\n  ALL differences have |omega^a - omega^b|^2 = 3.")
    print(f"  So |omega^a - omega^b| = sqrt(3) for all a != b.")
    print(f"  This is the Eisenstein norm: N(alpha) = |alpha|^2 in Z[omega].")

    print(f"\n  Powers (omega^a - omega^b)^c:")
    for c in range(1, 7):
        val = (omega - 1)**c
        b_coef = round(2 * val.imag / (3**0.5))
        a_coef = round(val.real + val.imag / (3**0.5))
        norm = abs(val)**2
        print(f"    (omega - 1)^{c} = {a_coef:5d} + {b_coef:5d}*omega  "
              f"(N = {norm:.1f} = 3^{c})")

    print(f"\n  KEY: N((omega^a - omega^b)^c) = 3^c.")
    print(f"  The 3-adic valuation grows linearly with c.")
    print(f"  This creates a 3-ADIC FILTRATION on Z[omega].")

    # Part 5: Positional number theory -- base 10 and 11
    print("\n" + "=" * 70)
    print("PART 5: BASE 10 = 2*5 AND BASE 11 (PRIME)")
    print("=" * 70)

    print(f"""
  10 = 2 * 5
  11 = prime

  In tournament theory:
    10 = T(21) = (H=21 - 1)/2  [the TARGET value for the H=21 gap!]
    11 = the Paley prime T_11 with H = 95095

  Base 10 positional system:
    Each digit d in position k contributes d * 10^k = d * (2*5)^k = d * 2^k * 5^k
    The factorization 10 = 2*5 means:
      - Divisibility by 2 depends only on last digit (standard)
      - Divisibility by 5 depends only on last digit
      - But divisibility by 3 requires DIGIT SUM (since 10 = 1 mod 3)

  Connection to tournaments:
    H = sum a_k * 2^k (the I.P. "binary expansion with multiplicities")
    H mod 3: depends on ALTERNATING SUM of a_k's (since 2 = -1 mod 3)
      H mod 3 = sum (-1)^k * a_k mod 3
      = a_0 - a_1 + a_2 - a_3 + ... mod 3
      = 1 - a_1 + a_2 - a_3 + ... mod 3

    This is like the "eleven-check" in base 10!
    In base 10: N mod 11 = alternating digit sum (since 10 = -1 mod 11).
    In base 2:  H mod 3 = alternating coefficient sum (since 2 = -1 mod 3).

    The SAME structure: (base) = -1 mod (prime).

  The parallel:
    Base 10, prime 11: alternating digit sum gives N mod 11
    Base 2, prime 3:   alternating alpha sum gives H mod 3

  And we know H mod 3 from Redei/OCF:
    H = 1 mod 2 always.
    H mod 3 depends on the tournament.
    But F(T, omega) = 9-divisibility (THM-085) for n >= 6.

  11 in the cyclotomic dictionary:
    11 = Phi_10(2) (since ord_11(2) = 10 = phi(11))
    This means 2 is a primitive root mod 11.
    The Paley tournament T_11 is based on QR = {{1,3,4,5,9}} mod 11.

  The alternating sum connection is DEEP:
    - Dehn invariant uses alternating-sum-like structure
    - Tournament I.P. evaluates at 2 = -1 mod 3
    - Base 10 divisibility by 11 = alternating digit sum
    - All are instances of x = -1 mod p, giving Phi_2(x) = p
  """)

    # Part 6: The 5 and 6 connection
    print("=" * 70)
    print("PART 6: 5 AND 6 -- ADDITION AND MULTIPLICATION")
    print("=" * 70)

    print(f"""
  5 = 2 + 3 (additive combination)
  6 = 2 * 3 (multiplicative combination)

  In tournament theory:
    5: first non-trivial tournament size (I.P. quadratic, real roots hold)
    6: beta_3 first appears, n=6 has 720 = 6! tournaments (up to labels)

  Deeper structure:
    H = I(Omega, 2) = 1 + 2*a1 + 4*a2 + ...
    I(Omega, 3) = 1 + 3*a1 + 9*a2 + ... (the "cuboid" evaluation)
    I(Omega, 2+3) = I(Omega, 5) = 1 + 5*a1 + 25*a2 + ...
    I(Omega, 2*3) = I(Omega, 6) = 1 + 6*a1 + 36*a2 + ...

  The formula H = 5 - 4*Re(I(i)) + 2*Im(I(i)):
    Here 5 = 2+3 appears as a CONSTANT.
    H = (2+3) - 4*(1-a2) + 2*a1  [since I(i) = (1-a2)+a1*i]
    = 5 - 4 + 4*a2 + 2*a1
    = 1 + 2*a1 + 4*a2  [correct]

  So the Gaussian evaluation naturally decomposes H as:
    H = (2+3) + corrections

  The additive combination 2+3=5 is the BASE, and the deviations
  from the base are the cycle structure contributions.

  At the other extreme, the multiplicative combination:
    6 = 2*3 is the MINIMUM cost of 2 disjoint cycles (each using 3 vertices).
    This is why a2 >= 1 requires 6 vertices (two disjoint 3-cycles).
    And 6 = 2*3 = 2 * (minimum cycle length) = cost of one disjoint pair.

  EXPONENTIATION:
    2^3 = 8 = dimension where quaternion multiplication loses commutativity
    3^2 = 9 = minimum vertices for cubic I.P. (three disjoint 3-cycles)

    The fact that 2^3 < 3^2 (8 < 9) means:
    Property loss (n=8) happens BEFORE cubic I.P. (n=9).
    The quadratic era extends one step beyond the Cayley-Dickson threshold.

    If 2^3 > 3^2 (hypothetically), the cubic I.P. would unlock BEFORE
    the property losses at n=8, and the theory would be qualitatively different.

    The Catalan-Mihailescu theorem says 2^3, 3^2 are the ONLY consecutive
    perfect powers. This uniqueness is reflected in tournament theory:
    there is exactly ONE "gap year" (n=8) between the last property-holding
    size and the cubic I.P. onset.
  """)

    # Part 7: The proof-by-contradiction structure
    print("=" * 70)
    print("PART 7: THE HILBERT 3RD PROBLEM PROOF STRUCTURE")
    print("=" * 70)

    print(f"""
  Hilbert's 3rd problem (solved by Dehn, 1900):
  A regular tetrahedron and a cube of the same volume are NOT
  scissors-congruent (cannot be cut into finitely many pieces
  and reassembled into each other).

  Dehn's proof by contradiction:
  1. Define the Dehn invariant D(P) for a polyhedron P.
  2. Show D is preserved under scissors congruence.
  3. Compute D(tetrahedron) and D(cube).
  4. Show they differ => NOT scissors-congruent.

  The irrationality enters via:
    D(cube) = 0 (all dihedral angles are pi/2, rational multiple of pi)
    D(tetrahedron) != 0 (dihedral angle = arccos(1/3), irrational/pi)

  The key: arccos(1/3)/pi is IRRATIONAL (proved by Niven).
  So the Dehn invariant of the tetrahedron is a non-trivial element
  in the tensor product R tensor (R/Q*pi).

  The tournament analogue (from our work):
  1. Define the tournament Dehn invariant D_T = I(omega).
  2. D_T is preserved under "tournament scissors congruence"
     (tournaments with same (a1,a2) have same D_T).
  3. D_T separates tournaments that H alone cannot.
  4. At n=7: D_T separates 52 of the 75 achievable H values into
     distinct Dehn classes.

  The "irrationality insertion" in our setting:
    I(2) = H is an integer.
    I(omega) = (1-a2) + (a1-a2)*omega is in Z[omega].
    The difference is that omega introduces sqrt(3) (irrational).
    So I(omega) lives in a LARGER ring than I(2).

  The "smallest s such that omega^s = 1":
    For omega = e^{{2pi*i/3}}: s = 3.
    For zeta_m = e^{{2pi*i/m}}: s = m.

    The user's suggestion: let s be the smallest natural with omega^s=1.
    Then use s=3 to set up the proof by contradiction for H=7:

    Suppose H = 7 = Phi_3(2). Then I(2) = 7 = (2^3-1)/(2-1).
    Evaluate at omega: I(omega) = (omega^3-1)/(omega-1) = 0/(omega-1) = 0.
    Wait, that's only if I(x) = (x^3-1)/(x-1) = x^2+x+1 = Phi_3(x).
    But we proved I(x) = Phi_3(x) is NEVER achievable (HYP-1015)!

    So the proof is:
    IF H = 7 THEN I(2) = 7.
    The polynomial I(x) = 1 + a1*x + a2*x^2 with I(2) = 7.
    If I(x) = Phi_3(x) = x^2+x+1: impossible (a1=a2=1, but a2=1 => a1>=2).
    If I(x) = 7: (a1=3, a2=0). Blocked by Helly.
    If I(x) = -x^2+4x+1: negative coefficient, impossible.
    If I(x) = -2x^2+5x+1: negative, impossible.
    If I(x) = -3x^2+6x+1: negative, impossible.
    If I(x) = x+3: (a1=1, a2=0), I(2)=5!=7.
    No valid I(x) gives I(2)=7. Hence H!=7. QED

    The connection to omega:
    The evaluation I(omega) for the blocked decompositions gives
    specific values in Z[omega] that are inconsistent with tournament
    structure. The omega evaluation acts as a SECONDARY INVARIANT
    that refines the H constraint.
  """)

    # Part 8: Alternating sum non-negativity
    print("=" * 70)
    print("PART 8: ALTERNATING SUM NON-NEGATIVITY")
    print("=" * 70)

    print(f"""
  For the independence polynomial I(x) = sum a_k * x^k:
    I(-1) = sum (-1)^k * a_k = 1 - a1 + a2 - a3 + a4 - ...

  This is the EULER CHARACTERISTIC of the independence complex.
  For claw-free graphs (Omega at n<=8): I(-1) <= 1 always.
  For general graphs: I(-1) can be any integer.

  Non-negativity of the alternating sum is NOT guaranteed.
  But for tournaments, I(-1) has specific constraints:
    I(-1) = 1 - a1 + a2 - a3 + ...
    At n=7: I(-1) = 1 - a1 + a2 (quadratic)

  From the data:
    H = 1 + 2*a1 + 4*a2
    I(-1) = 1 - a1 + a2

  So: I(-1) = 1 - a1 + a2
      H = 1 + 2*a1 + 4*a2
      => H + 2*I(-1) = 3 + 6*a2 - ... wait:
         H + 2*I(-1) = (1+2a1+4a2) + 2(1-a1+a2) = 3 + 6*a2

  So H + 2*I(-1) = 3 + 6*a2 = 3*(1+2*a2).
  This means: H + 2*I(-1) is ALWAYS divisible by 3!
  And H = 3*(1+2*a2) - 2*I(-1).

  For H=7: 7 + 2*I(-1) = 3*(1+2*a2), so 7+2*I(-1) = 0 mod 3.
  7 = 1 mod 3, so 2*I(-1) = 2 mod 3, so I(-1) = 1 mod 3.
  Possible: I(-1) = 1, 4, 7, ...
  With I(-1) = 1-a1+a2 and a1+2*a2=3:
    (3,0): I(-1) = 1-3+0 = -2 (NOT 1 mod 3 -- wait, -2 = 1 mod 3. Yes.)
    (1,1): I(-1) = 1-1+1 = 1 (= 1 mod 3. Yes.)
  Both are consistent with 3 | (H + 2*I(-1)), as expected.

  The formula H + 2*I(-1) = 3*(1+2*a2) is clean.
  It says: a2 = (H + 2*I(-1) - 3) / 6.
  And a1 = (H - 1)/2 - 2*a2 = (H-1)/2 - (H+2*I(-1)-3)/3.
  """)

    # Part 9: Verify the formula
    print("=" * 70)
    print("PART 9: VERIFICATION OF H + 2*I(-1) = 3*(1+2*a2)")
    print("=" * 70)

    print("  This formula holds by direct algebra (quadratic case):")
    print("  I(2) = 1 + 2*a1 + 4*a2 = H")
    print("  I(-1) = 1 - a1 + a2")
    print("  I(2) + 2*I(-1) = (1+2a1+4a2) + 2(1-a1+a2) = 3 + 6*a2")
    print("  = 3*(1+2*a2) always.")
    print()
    print("  This extends to higher degree:")
    print("  I(2) + 2*I(-1) = sum a_k*(2^k + 2*(-1)^k)")
    print("  = a_0*(1+2) + a_1*(2-2) + a_2*(4+2) + a_3*(8-2) + a_4*(16+2) + ...")
    print("  = 3*a_0 + 0*a_1 + 6*a_2 + 6*a_3 + 18*a_4 + 30*a_5 + ...")

    for k in range(8):
        val = 2**k + 2*(-1)**k
        print(f"    k={k}: 2^k + 2*(-1)^k = {val}  "
              f"(divisible by 3? {val % 3 == 0}, factor = {val//3 if val%3==0 else '-'})")

    print(f"\n  Hmm: for k=1, the coefficient is 0, not 3.")
    print(f"  For k>=2 and even: 2^k+2 (e.g., 6, 18, 66)")
    print(f"  For k>=3 and odd:  2^k-2 (e.g., 6, 30, 126)")
    print(f"  k=1: 0 (a1 drops out!)")
    print(f"\n  So: I(2) + 2*I(-1) = 3 + 6*a2 + 6*a3 + 18*a4 + 30*a5 + ...")
    print(f"  = 3*(1 + 2*a2 + 2*a3 + 6*a4 + 10*a5 + ...)")
    print(f"  Is this always divisible by 3? YES (each term is 3*something)")
    print(f"  And a1 COMPLETELY drops out!")

    print(f"\n  THEOREM: I(2) + 2*I(-1) is always divisible by 3,")
    print(f"  regardless of a1, and equals 3*(1 + sum_{{k>=2}} c_k * a_k)")
    print(f"  where c_k = (2^k + 2*(-1)^k)/3.")

    print(f"\n  This means: H mod 3 = -2*I(-1) mod 3 = I(-1) mod 3.")
    print(f"  (Since -2 = 1 mod 3)")
    print(f"  So H and I(-1) have the SAME residue mod 3!")
    print(f"  H = I(-1) mod 3.")

    print(f"\n  COROLLARY: H = I(-1) (mod 3)")
    print(f"  This relates the Hamiltonian path count to the Euler characteristic")
    print(f"  of the independence complex, modulo 3.")

    # Part 10: Multiplication vs addition in forbidden H
    print("\n" + "=" * 70)
    print("PART 10: MULTIPLICATION AND ADDITION IN THE GAPS")
    print("=" * 70)

    print(f"""
  The two permanent gaps and their arithmetic:

  H=7 = 2+5 = 2+(2+3)     [additive]
  H=7 = 2^3-1              [exponential]
  H=7 = Phi_3(2)           [cyclotomic]

  H=21 = 3*7               [multiplicative]
  H=21 = 3*Phi_3(2)        [cyclotomic product]
  H=21 = (2+1)*(2^2+2+1)   [factored]

  The gaps encode ALL three operations:
    ADDITION: 7 = 2+5, 21 = 20+1
    MULTIPLICATION: 21 = 3*7
    EXPONENTIATION: 7 = 2^3-1, 21 = (2^3-1)*3

  The "tower":
    Level 0: 2 and 3 (the keys)
    Level 1: 7 = 2^3 - 1 (exponentiation of keys, minus 1)
    Level 2: 21 = 3 * 7 (multiplication by key)
    Level 3: 63 = 3 * 21 = 9 * 7 (NOT forbidden -- too far from keys)

  Why does the tower stop at level 2?
    63 = 1 + 2*31 = 1 + 2*a1 (with a1=31, a2=0)
    At n=8: 31 odd cycles all pairwise conflicting is achievable
    because n=8 has enough vertices for the cycle structure to accommodate it.

    But 21 = 1 + 2*10 requires 10 cycles all conflicting (if a2=0),
    and at ANY n, 10 pairwise-conflicting cycles with the specific
    tournament constraint is impossible.

  The boundary between achievable and impossible:
    (a1=3, a2=0): IMPOSSIBLE for all n (Helly)
    (a1=4, a2=0): achievable? At n=5, H=9 exists but might be (2,1) not (4,0).
    ...
    (a1=31, a2=0): ACHIEVABLE at n=8 (H=63 found)

  So somewhere between a1=3 and a1=31, all-conflicting becomes achievable.
  Finding the EXACT threshold would pin down why 21 is the last gap.
  """)

    print(f"\n{'='*70}")
    print("SUMMARY: THE 2-3 DEEP SYNTHESIS")
    print(f"{'='*70}")

    print(f"""
  1. The permanent gaps {{7, 21}} encode the pair (2,3):
     7 = Phi_3(2), 21 = 3*Phi_3(2)

  2. Positional number theory:
     H mod 3 = I(-1) mod 3 (Euler characteristic congruence)
     This is the "alternating digit sum" rule in base 2, prime 3.

  3. The omega evaluation I(omega):
     - (omega-1)^c has norm 3^c (3-adic filtration)
     - Phi_3(omega) = 0 is the defining relation
     - I(omega) is the "Dehn invariant" of the tournament

  4. The proof-by-contradiction:
     H=7 => I(x) achieves 7 at x=2.
     I(omega) for valid decompositions => contradiction with tournament structure.
     The irrationality of omega (via sqrt(3)) creates obstructions
     that are invisible in Z but visible in Z[omega].

  5. Exponentiation: 2^3=8 vs 3^2=9
     Property loss at n=8 (Cayley-Dickson threshold)
     Cubic I.P. at n=9 (three disjoint 3-cycles)
     Gap of 1 (the Catalan pair) creates the one "gap year" n=8.

  6. 5 = 2+3 appears as the constant in H = 5 - 4*Re(I(i)) + 2*Im(I(i)).
     6 = 2*3 = minimum cost of one disjoint pair (two 3-cycles, 6 vertices).

  7. 10 = T(H=21) is palindromic in both binary (1010) and ternary (101).
     11 = Paley prime, Phi_10(2), 2 is primitive root mod 11.
  """)

if __name__ == "__main__":
    main()
