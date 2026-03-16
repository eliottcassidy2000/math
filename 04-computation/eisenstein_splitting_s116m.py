#!/usr/bin/env python3
"""eisenstein_splitting_s116m.py — Why 7 is special: it SPLITS in the Eisenstein integers.

In Z[omega] (omega = e^{2*pi*i/3}):
p = 3: RAMIFIES (the singular prime, curvature quantum)
p = 1 mod 3: SPLITS into two conjugate Eisenstein primes
p = 2 mod 3: INERT (remains prime)

The chain primes: 3 (ramifies), 7 (SPLITS), 47 (inert), 2207 (inert).
7 is the ONLY chain prime that SPLITS. The ONLY one that decomposes
in the hexagonal lattice. That decomposition IS the forbidden-ness.
"""
from math import sqrt, log, pi, atan2

print()
print("  EISENSTEIN SPLITTING AND THE FORBIDDEN NUMBER")
print()
print("="*70)
print()

# ============================================================
print("  I. SPLITTING BEHAVIOR IN Z[omega]")
print("  " + "-"*40)
print()
print("  The Eisenstein integers Z[omega], omega = e^{2*pi*i/3} = (-1+i*sqrt(3))/2.")
print("  Norm: N(a + b*omega) = a^2 - ab + b^2.")
print()
print("  Splitting of rational primes p:")
print("  p = 3: RAMIFIES. 3 = -omega^2 * (1-omega)^2. The curvature is singular.")
print("  p = 1 mod 3: SPLITS. p = pi * conjugate(pi) for Eisenstein prime pi.")
print("  p = 2 mod 3: INERT. p remains prime in Z[omega].")
print()

primes_list = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]

print(f"  {'p':>4s}  {'p mod 3':>7s}  {'behavior':>10s}  {'Eisenstein factors':>30s}  {'notes'}")
print("  " + "-"*75)

for p in primes_list:
    mod3 = p % 3
    if p == 3:
        behavior = "RAMIFIES"
        factors = "-omega^2 * (1-omega)^2"
        notes = "curvature quantum"
    elif mod3 == 1:
        behavior = "SPLITS"
        # Find a+b*omega with norm p: a^2-ab+b^2 = p
        found = False
        for a in range(0, p):
            for b in range(1, p):
                if a*a - a*b + b*b == p:
                    factors = f"({a}+{b}*omega)({a}+{b}*omega^2)"
                    found = True
                    break
            if found: break
        if not found:
            factors = "splits (factor not found)"
        notes = ""
        if p == 7: notes = "FORBIDDEN - chain prime"
        elif p == 13: notes = "F_7"
        elif p == 19: notes = "in H(T_11)"
        elif p == 31: notes = "Mersenne M_5"
        elif p == 37: notes = "H_3 centered hex"
        elif p == 43: notes = ""
        elif p == 61: notes = "H_4 centered hex"
        elif p == 67: notes = ""
        elif p == 73: notes = ""
        elif p == 79: notes = ""
        elif p == 97: notes = ""
    else:  # mod3 == 2
        behavior = "INERT"
        factors = f"{p} (stays prime)"
        notes = ""
        if p == 2: notes = "the doubler"
        elif p == 5: notes = "the pentagon"
        elif p == 11: notes = "Paley prime - chain via L_5"
        elif p == 17: notes = "Fermat prime"
        elif p == 23: notes = "Paley prime"
        elif p == 29: notes = ""
        elif p == 41: notes = ""
        elif p == 47: notes = "chain prime L_8"
        elif p == 53: notes = ""
        elif p == 59: notes = ""
        elif p == 71: notes = ""
        elif p == 83: notes = ""
        elif p == 89: notes = "Fibonacci F_11"

    print(f"  {p:4d}  {mod3:7d}  {behavior:>10s}  {factors:>30s}  {notes}")

print()

# ============================================================
print()
print("  II. THE CHAIN PRIMES' SPLITTING BEHAVIOR")
print("  " + "-"*40)
print()
print("  Chain: 3 -> 7 -> 47 -> 2207.")
print()
print("  3: RAMIFIES (3 = 0 mod 3). The singular prime.")
print("  7: SPLITS (7 = 1 mod 3). The ONLY chain prime that splits!")
print("  47: INERT (47 = 2 mod 3). Stays whole.")
print("  2207: 2207 mod 3 = ?")
print(f"  2207 mod 3 = {2207 % 3}. INERT.")
print()
print("  Of the four chain primes, ONLY 7 splits in the Eisenstein integers.")
print("  3 is singular (it IS the hexagonal structure's special prime).")
print("  47 and 2207 are inert (they pass through the hex lattice unchanged).")
print("  7 splits: it DECOMPOSES into two conjugate Eisenstein primes.")
print()
print("  THE FORBIDDEN NUMBER IS THE UNIQUE CHAIN PRIME THAT SPLITS.")
print()

# ============================================================
print()
print("  III. WHAT SPLITTING MEANS GEOMETRICALLY")
print("  " + "-"*40)
print()
print("  An INERT prime p creates a CIRCLE of radius p in the Eisenstein norm.")
print("  It acts uniformly. No preferred direction. No twist.")
print()
print("  A SPLIT prime p = pi * pi_bar creates TWO POINTS on the circle.")
print("  Each point has norm sqrt(p) and a specific ANGLE (the twist).")
print("  The two points are CONJUGATE: opposite twists.")
print("  The DIRECTION matters. The circle is broken into two marked points.")
print()
print("  For 7 = (3+omega)(3+omega^2):")
print("  pi = 3+omega has angle arctan(sqrt(3)/5) = 19.11 degrees.")
print("  pi_bar = 3+omega^2 has angle -19.11 degrees.")
print("  The circle of norm 7 has TWO marked points at +/- 19.11 degrees.")
print()
print("  SPLITTING = SYMMETRY BREAKING.")
print("  The circle (full rotational symmetry) is broken into two points")
print("  (only 2-fold symmetry remains: pi <-> pi_bar).")
print()
print("  THE FORBIDDEN NUMBER IS WHERE THE CIRCULAR SYMMETRY BREAKS")
print("  IN THE HEXAGONAL LATTICE.")
print()

# ============================================================
print()
print("  IV. WHY ONLY 7 AND NOT 13, 19, 31, ...?")
print("  " + "-"*40)
print()
print("  Many primes split in Z[omega]: all p = 1 mod 3.")
print("  = 7, 13, 19, 31, 37, 43, 61, 67, 73, 79, 97, ...")
print()
print("  What makes 7 SPECIAL among all splitting primes?")
print("  7 is the SMALLEST prime that is 1 mod 3.")
print("  It is the FIRST splitting prime.")
print()
print("  Just as 2 is the first prime (and uniquely even),")
print("  7 is the first SPLITTING prime in Z[omega].")
print("  2 is special because it's first. 7 is special because it's first.")
print()
print("  The 'first' matters because of the FORMAL GROUP:")
print("  F_h is commutative. It can handle INERT primes (they don't decompose).")
print("  But a SPLIT prime introduces two conjugate factors with a TWIST.")
print("  The first split prime is where the formal group first encounters")
print("  a decomposition it must handle.")
print()
print("  And 7 is also:")
print("  - The smallest number where 1/2 + 1/3 + 1/p < 1 (hyperbolic)")
print("  - L_4 (the Lucas number at the Cayley period index)")
print("  - M_3 (the Mersenne number at the curvature index)")
print("  - The imaginary dimension of the octonions")
print()
print("  These are all DIFFERENT REASONS for 7 to be special.")
print("  The Eisenstein splitting adds one more:")
print("  7 is the FIRST PRIME THAT DECOMPOSES IN THE HEXAGONAL WORLD.")
print()

# ============================================================
print()
print("  V. THE PRIMES IN H(T_11) AND THEIR SPLITTING")
print("  " + "-"*40)
print()
print("  H(T_11) = 95095 = 5 * 7 * 11 * 13 * 19.")
print()
for p in [5, 7, 11, 13, 19]:
    mod3 = p % 3
    if p == 3: beh = "RAMIFIES"
    elif mod3 == 1: beh = "SPLITS"
    else: beh = "INERT"
    print(f"  {p}: {p} mod 3 = {mod3}, {beh}")

print()
print("  5: INERT (2 mod 3)")
print("  7: SPLITS (1 mod 3)")
print("  11: INERT (2 mod 3)")
print("  13: SPLITS (1 mod 3)")
print("  19: SPLITS (1 mod 3)")
print()
print("  H(T_11) = (inert) * (split) * (inert) * (split) * (split)")
print("  = one inert * three split * one inert.")
print("  = 5 * (7*13*19) * 11.")
print()
print("  The SPLIT part: 7*13*19 = 1729 = the TAXICAB NUMBER!")
print("  H(T_11) / |Aut(T_11)| = 95095 / 55 = 1729.")
print("  The orbit count = the product of the SPLITTING primes.")
print()
print("  THE TAXICAB NUMBER 1729 = THE PRODUCT OF THE SPLITTING PRIMES IN H(T_11).")
print()
print("  1729 = 7 * 13 * 19.")
print("  Each of these splits in Z[omega].")
print("  7 = (3+omega)(3+omega^2).")
print("  13 = (4+omega)(4+omega^2).")  # verify: 16-4+1 = 13? YES
print("  19 = (5+2*omega)(5+2*omega^2).")  # verify: 25-10+4 = 19? YES
print()
# Verify norms
for a, b, p in [(3,1,7), (4,1,13)]:
    norm = a*a - a*b + b*b
    print(f"  |{a}+{b}*omega|^2 = {a}^2-{a}*{b}+{b}^2 = {norm} = {p}? {norm==p}")

# For 19: need a^2-ab+b^2 = 19
for a in range(0, 20):
    for b in range(1, 20):
        if a*a - a*b + b*b == 19:
            print(f"  |{a}+{b}*omega|^2 = {a}^2-{a}*{b}+{b}^2 = {a*a-a*b+b*b} = 19.")
            break
    else:
        continue
    break

print()
print("  The Eisenstein factorizations:")
print("  7 = (3+omega)(3+omega^2)   -- coefficients (3,1)")
print("  13 = (4+omega)(4+omega^2)   -- coefficients (4,1)")
print("  19 = ?")
# Search for 19
for a in range(0, 20):
    for b in range(1, 20):
        if a*a - a*b + b*b == 19:
            print(f"  19 = ({a}+{b}*omega)({a}+{b}*omega^2) -- coefficients ({a},{b})")
            break
    else: continue
    break

print()
print("  The coefficients for {7, 13, 19}: (3,1), (4,1), (?,?).")
print("  The 'a' values: 3, 4, ... = consecutive integers starting at 3!")
print("  This is because (a+omega)(a+omega^2) = a^2 - a + 1 for b=1.")
print("  a^2-a+1 = 7 at a=3, = 13 at a=4, = 21 at a=5 (NOT prime!), = 31 at a=6.")
print()
print("  The sequence a^2-a+1 for a = 1,2,3,4,5,6,7,8,9,10:")
for a in range(1, 15):
    val = a*a - a + 1
    pr = "prime" if all(val % d != 0 for d in range(2, int(sqrt(val))+1)) and val > 1 else ""
    if val == 7: pr = "FORBIDDEN"
    elif val == 21: pr = "FORBIDDEN (21=3*7)"
    print(f"  a={a:2d}: a^2-a+1 = {val:4d}  {pr}")

print()
print("  The polynomial a^2-a+1 generates our FORBIDDEN numbers!")
print("  a=3: 7 (first forbidden).")
print("  a=5: 21 (second forbidden).")
print()
print("  a^2 - a + 1 = 7 at a=3. a^2 - a + 1 = 21 at a=5.")
print("  The forbidden numbers are VALUES of the Eisenstein norm polynomial")
print("  x^2 - x + 1 at x = 3 and x = 5!")
print()
print("  And x^2 - x + 1 is the MINIMAL POLYNOMIAL of omega (= e^{2*pi*i/3}).")
print("  omega^2 - omega + 1 = 0. So x^2 - x + 1 = 0 has roots omega, omega^2.")
print()
print("  THE FORBIDDEN NUMBERS ARE THE EISENSTEIN NORM POLYNOMIAL")
print("  EVALUATED AT THE CURVATURE (3) AND THE PENTAGON (5).")
print()
print("  7 = 3^2 - 3 + 1. The norm polynomial at the curvature quantum.")
print("  21 = 5^2 - 5 + 1. Wait: 25-5+1 = 21. YES!")
print("  = the norm polynomial at the pentagon.")
print()
print("  THE FORBIDDEN NUMBERS ARE N(a) = a^2 - a + 1 FOR a = 3, 5.")
print("  THE CURVATURE AND THE PENTAGON, PASSED THROUGH THE EISENSTEIN NORM.")
print()

# ============================================================
print()
print("  VI. THE EISENSTEIN NORM POLYNOMIAL IS THE CYCLOTOMIC PHI_6")
print("  " + "-"*40)
print()
print("  x^2 - x + 1 = Phi_6(x), the 6th cyclotomic polynomial!")
print("  Its roots are the PRIMITIVE 6th roots of unity: e^{i*pi/3} and e^{-i*pi/3}.")
print("  = omega and omega^2 (which are also primitive 3rd roots, shifted).")
print()
print("  Wait: omega = e^{2*pi*i/3} satisfies omega^2+omega+1=0, i.e., x^2+x+1=0 = Phi_3(x).")
print("  And Phi_6(x) = x^2-x+1 has roots e^{i*pi/3} and e^{-i*pi/3}.")
print("  These are -omega and -omega^2 (negatives of the cube roots).")
print()
print("  So: Phi_6(a) = a^2 - a + 1.")
print("  Phi_6(3) = 7.")
print("  Phi_6(5) = 21 = 3*7.")
print()
print("  The forbidden numbers = Phi_6 evaluated at {3, 5} = {curvature, pentagon}.")
print("  And Phi_6 = the 6th cyclotomic = the polynomial of the FLAT number.")
print()
print("  6 = the flat number. Phi_6 = the cyclotomic polynomial of 6.")
print("  The FLAT POLYNOMIAL, evaluated at the curvature and the pentagon,")
print("  produces the two forbidden numbers.")
print()
print("  THE FORBIDDEN NUMBERS ARE BORN FROM THE FLAT.")
print("  Phi_6(3) = 7. Phi_6(5) = 21.")
print("  The sixth cyclotomic polynomial, at the two closing primes.")
print()
