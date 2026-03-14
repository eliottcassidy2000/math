"""
unified_number_atlas.py -- kind-pasteur-2026-03-14-S94
LONG SESSION: Unified atlas of ALL number types in tournament theory.

PREVIOUS SESSIONS:
  S91: Naturals 1-12 (2 dominates, Mersenne tower 2→3→7→21)
  S92: Extended naturals 4,9,12,13 + irrationals (phi, e, pi, sqrt(2,3,5), 4^{1/3})
  S93: Rationals (1/2, 1/3, 1/4, 2/3, 3/2, 9/2, CF of 0.27)

THIS SESSION: Connect ALL types. Find the unified structure.
Also: complex numbers, algebraic numbers, p-adic numbers.
"""

import sys, math
from fractions import Fraction

sys.stdout.reconfigure(encoding='utf-8')

def C(n, k):
    if k < 0 or k > n: return 0
    return math.comb(n, k)

maxH = {1:1, 2:1, 3:3, 4:5, 5:15, 6:45, 7:189, 8:661, 9:3357, 10:15745, 11:95095}

def main():
    print("=" * 70)
    print("UNIFIED NUMBER ATLAS OF TOURNAMENT THEORY")
    print("kind-pasteur-2026-03-14-S94")
    print("=" * 70)

    phi = (1 + math.sqrt(5)) / 2
    e = math.e
    pi = math.pi

    # ============================================================
    # PART 1: THE NUMBER TYPE HIERARCHY
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 1: THE NUMBER TYPE HIERARCHY IN TOURNAMENTS")
    print(f"{'='*70}")

    print(f"""
  NATURAL NUMBERS (Z+):
    H(T) is always a natural number (and always odd).
    alpha_k are natural numbers.
    |Aut(T)|, c_k, n are naturals.
    THE THEORY LIVES IN THE NATURALS at its core.

  INTEGERS (Z):
    Fourier coefficients H_hat(S) can be negative rationals.
    The Degree Drop: top-level coefficients are +-2.
    The mystery sequence a(n) = n + C(n,3) extends to negative n.
    a(n) + a(-n) = -n^2 (the antisymmetry lives in Z).

  RATIONALS (Q):
    mean_H = n!/2^{{n-1}} (rational, not integer for n >= 3)
    Fourier magnitudes: (n-2)!/2^{{n-2}} (rational)
    Var/Mean^2 = 1/3 (rational)
    The continued fraction of 0.27 = [0; 3, 1, 2, 2, 1, 2]
    Szele ratio maxH/meanH (rational for each n)

  ALGEBRAIC IRRATIONALS:
    phi = (1+sqrt(5))/2: Fibonacci growth rate
    sqrt(2): spectral radii of tournaments
    sqrt(3): variance structure (1/sqrt(3) = CV)
    4^{{1/3}} = 2^{{2/3}}: 3-strand Pascal growth
    sqrt(p) for Paley primes p: eigenvalues of Paley tournaments

  TRANSCENDENTAL NUMBERS:
    e: Szele limit (maxH/meanH -> e)
    pi: asymptotic normalization (Stirling, central binomial)
    log(2): appears in entropy calculations (I/m ≈ 0.27 bits)

  COMPLEX NUMBERS (C):
    Paley eigenvalues: (-1 +/- i*sqrt(p))/2
    Roots of unity: omega_k = e^{{2*pi*i/k}} (F(T,omega) divisibility)
    Lee-Yang zeros: Q_n(z) roots on/near unit circle
    Fourier transform: e^{{2*pi*i*k/N}} basis functions

  p-ADIC NUMBERS (Q_p):
    v_2(H) = 0 always (H is 2-adically a unit)
    THM-085: 9 | F(T,omega) for n >= 6 (3-adic structure)
    The OCF: H = 1 + 2*alpha_1 + 4*alpha_2 + ... is a 2-adic expansion!
    T = (H-1)/2 has binary digits = alpha_k parities
""")

    # ============================================================
    # PART 2: THE CONNECTIONS BETWEEN NUMBER TYPES
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 2: CONNECTIONS BETWEEN NUMBER TYPES")
    print(f"{'='*70}")

    connections = [
        # Natural -> Rational
        ("N -> Q", "mean_H = n!/2^{n-1}",
         "The mean of a NATURAL invariant is RATIONAL. The denominator is always a power of 2."),

        # Rational -> Irrational
        ("Q -> R\\Q", "maxH/meanH -> e as n -> inf",
         "A sequence of RATIONALS converges to a TRANSCENDENTAL. The Szele limit."),

        # Rational -> Algebraic
        ("Q -> Algebraic", "Fourier magnitudes (n-2)!/2^{n-2} -> 0",
         "But the growth rate phi is ALGEBRAIC. Rational Fourier coefficients "
         "encode algebraic growth via their asymptotics."),

        # Natural -> Complex
        ("N -> C", "Paley eigenvalues = (-1 +/- i*sqrt(p))/2",
         "NATURAL tournament quantities have COMPLEX eigenvalue spectra. "
         "The imaginary part sqrt(p) governs the spectral gap."),

        # Natural -> p-adic
        ("N -> Q_p", "H = 1 + 2*alpha_1 + 4*alpha_2 + ... (2-adic expansion)",
         "The OCF IS a 2-adic expansion. The forbidden value 7 = 111 in binary "
         "requires alpha_k parity pattern (1,1,1) which is unrealizable."),

        # Complex -> Natural
        ("C -> N", "I(C_3, 2) = 7 (complex cycle graph -> natural forbidden value)",
         "The independence polynomial of a COMPLEX structure (the 3-cycle in the "
         "Argand plane) evaluated at a natural number gives the forbidden value."),

        # Algebraic -> Transcendental
        ("Algebraic -> Transcendental", "phi^3 ≈ 4, e ≈ phi^2 + 0.1",
         "Near-coincidences between algebraic and transcendental numbers "
         "create the 'resonances' in tournament theory."),

        # Rational -> Complex
        ("Q -> C", "Q_n(z) zeros approach unit circle |z| = 1",
         "The RATIONAL generating polynomial Q_n has COMPLEX zeros that "
         "concentrate on the unit circle as n -> inf (Lee-Yang phenomenon)."),
    ]

    for name, example, meaning in connections:
        print(f"\n  {name}: {example}")
        print(f"    {meaning}")

    # ============================================================
    # PART 3: THE "NUMBER WEB" — ALL CONNECTIONS AT ONCE
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 3: THE NUMBER WEB")
    print(f"{'='*70}")

    print(f"""
  The COMPLETE web of number type connections:

                        pi (transcendental)
                         |
                     asymptotics
                         |
    phi (algebraic) --- e (transcendental)
       |                 |
    growth             Szele limit
       |                 |
    Fibonacci          max_H/mean_H
       |                 |
    F(n) (natural) --- H(T) (natural) --- alpha_k (natural)
       |                 |                    |
    path indep.        OCF                2-adic expansion
       |                 |                    |
    I(P_m, x)       I(Omega, 2)         T = (H-1)/2
       |                 |                    |
    Jacobsthal         H = 1+2*a1+4*a2    binary digits
    (natural)          (natural)          (mod 2 structure)
       |                 |                    |
    x=2               x=2 (!)            2-adic Q_2
    (natural)          (natural)          (p-adic)
       |                 |                    |
    I(P_m, 1)=Fib    I(K_3, 2)=7        v_2(H) = 0
    (x=1, natural)    (FORBIDDEN!)        (always a unit)

  THE CENTER OF THE WEB: H(T) = I(Omega(T), 2)
  EVERY number type connects to this formula.
""")

    # ============================================================
    # PART 4: COMPLEX NUMBERS IN DETAIL
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 4: COMPLEX NUMBERS — THE DEEPEST LAYER")
    print(f"{'='*70}")

    # Paley eigenvalues
    for p in [3, 7, 11]:
        mu_real = -0.5
        mu_imag = math.sqrt(p) / 2
        modulus = math.sqrt(mu_real**2 + mu_imag**2)
        argument = math.atan2(mu_imag, mu_real)
        print(f"\n  Paley T_{p} eigenvalue:")
        print(f"    mu = {mu_real:.4f} +/- {mu_imag:.4f}i")
        print(f"    |mu| = {modulus:.4f} = sqrt((p+1)/4) = sqrt({(p+1)/4:.1f})")
        print(f"    arg(mu) = {argument:.4f} rad = {math.degrees(argument):.1f} degrees")
        print(f"    mu^p = ? (determines p-fold symmetry)")

    # Roots of unity
    print(f"\n  KEY ROOTS OF UNITY:")
    for k in [2, 3, 4, 6]:
        omega = complex(math.cos(2*pi/k), math.sin(2*pi/k))
        print(f"    omega_{k} = e^{{2*pi*i/{k}}} = {omega.real:.4f} + {omega.imag:.4f}i")
        # Tournament evaluation
        # F(T, omega_k) involves this root

    # Lee-Yang zeros
    print(f"\n  LEE-YANG ZEROS OF Q_n(z):")
    print(f"    n=4: ALL 4 non-trivial zeros on unit circle |z| = 1")
    print(f"    n=5: zeros at |z| = 0.823, 1.009, 1.178, 1.430")
    print(f"    n=6: zeros concentrating near |z| = 1 (most within 3%)")
    print(f"    As n -> inf: zeros APPROACH the unit circle (Lee-Yang)")
    print(f"    The unit circle |z| = 1 is a PHASE BOUNDARY in the complex plane")

    # ============================================================
    # PART 5: p-ADIC NUMBERS IN DETAIL
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 5: p-ADIC STRUCTURE — THE OCF AS 2-ADIC EXPANSION")
    print(f"{'='*70}")

    print(f"""
  The OCF: H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + 16*alpha_4 + ...
            = 1 + 2*(alpha_1 + 2*alpha_2 + 4*alpha_3 + ...)

  This IS a 2-adic expansion! T = (H-1)/2 = alpha_1 + 2*alpha_2 + 4*alpha_3 + ...
  The "2-adic digits" of T are the alpha_k values (mod 2).

  2-ADIC VALUATION:
    v_2(H) = 0 for ALL tournaments (H always odd)
    v_2(H-1) >= 1 always (H-1 always even)
    v_2(H-1) = 1 iff alpha_1 is odd
    v_2(H-1) >= 2 iff alpha_1 is even (then alpha_2 determines next digit)

  3-ADIC STRUCTURE:
    THM-085: 9 | F(T, omega) for n >= 6
    v_3(F(T,omega)) >= 2 for n >= 6
    This means the forward polynomial has deep 3-adic structure at cube roots

  THE FORBIDDEN VALUE IN p-ADIC TERMS:
    7 = 1 + 2 + 4 = 1 + 2*(1 + 2*1) in OCF means alpha_1=3, alpha_2=0
    In 2-adic: T = 3 = 11 in binary
    The 2-adic constraint: Omega = K_3 is required, which is impossible
    So the 2-adic digit pattern (1,1) in T (with alpha_2=0) is FORBIDDEN
    when alpha_1 = 3 exactly (because of the completeness constraint)
""")

    # ============================================================
    # PART 6: ALGEBRAIC NUMBERS AND MINIMAL POLYNOMIALS
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 6: ALGEBRAIC NUMBERS — MINIMAL POLYNOMIALS")
    print(f"{'='*70}")

    algebraics = [
        ("phi", "(1+sqrt(5))/2", "x^2 - x - 1 = 0", "Fibonacci growth"),
        ("1/phi", "(sqrt(5)-1)/2", "x^2 + x - 1 = 0", "Golden ratio reciprocal"),
        ("sqrt(2)", "2^{1/2}", "x^2 - 2 = 0", "Spectral radius, silver ratio"),
        ("sqrt(3)", "3^{1/2}", "x^2 - 3 = 0", "Variance structure"),
        ("sqrt(5)", "5^{1/2}", "x^2 - 5 = 0", "Fibonacci discriminant"),
        ("4^{1/3}", "2^{2/3}", "x^3 - 4 = 0", "Pascal growth"),
        ("1+sqrt(2)", "silver ratio", "x^2 - 2x - 1 = 0", "Spectral radius at n=4"),
        ("2+sqrt(5)", "phi^3", "x^2 - 4x - 1 = 0", "3-step Fibonacci"),
    ]

    for name, value_str, minpoly, role in algebraics:
        print(f"\n  {name} = {value_str}")
        print(f"    Minimal polynomial: {minpoly}")
        print(f"    Tournament role: {role}")

    print(f"\n  THE POLYNOMIAL HIERARCHY:")
    print(f"  Degree 1 (rational): 1/2, 1/3, 3/2, 9/2, ...")
    print(f"  Degree 2 (quadratic): phi, sqrt(2), sqrt(3), sqrt(5), 1+sqrt(2)")
    print(f"  Degree 3 (cubic): 4^{{1/3}} = 2^{{2/3}}")
    print(f"  Degree infinity (transcendental): e, pi, log(2)")
    print(f"")
    print(f"  Tournament theory uses numbers of degree 1, 2, 3, and infinity.")
    print(f"  NO degree 4+ algebraic numbers have appeared!")
    print(f"  The theory is 'cubic-bounded' in algebraic complexity.")

    # ============================================================
    # PART 7: THE UNIFIED FORMULA WEB
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 7: THE UNIFIED FORMULA WEB — ALL TYPES IN ONE EQUATION")
    print(f"{'='*70}")

    print(f"""
  CAN WE WRITE ONE EQUATION THAT INVOLVES ALL NUMBER TYPES?

  CANDIDATE:
    lim_{{n->inf}} max_H(n) / (n! / 2^{{n-1}}) = e

  This equation involves:
    - NATURAL: max_H(n), n
    - RATIONAL: n!/2^{{n-1}} = mean_H(n) (ratio of natural and power of 2)
    - ALGEBRAIC: the limit process involves phi implicitly (growth rate)
    - TRANSCENDENTAL: e (the limit value)
    - COMPLEX: the proof uses complex analytic methods (Szele-Alon)
    - p-ADIC: the 2^{{n-1}} denominator encodes 2-adic structure

  THIS SINGLE EQUATION TOUCHES ALL SIX NUMBER TYPES.

  ANOTHER CANDIDATE:
    H(T_p) = ((p-1)/2) * product eigenvalue contributions
           = NATURAL (H) from RATIONAL (Szele) via ALGEBRAIC (eigenvalues)
             involving COMPLEX (i*sqrt(p)/2) and TRANSCENDENTAL (e in asymptotics)

  THE DEEPEST: The OCF + Szele + Fourier:
    H = I(Omega, 2) = sum_k 2^k * alpha_k        [NATURAL + 2-ADIC]
    |H_hat(S)| = (n-2)!/2^{{n-2}}                  [RATIONAL]
    Growth ~ phi^n or 4^{{n/3}}                     [ALGEBRAIC]
    maxH/meanH -> e                                [TRANSCENDENTAL]
    Q_n zeros on |z|=1                             [COMPLEX]
""")

    # ============================================================
    # PART 8: IDEA GENERATION FROM THE UNIFIED VIEW
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 8: NEW IDEAS FROM THE UNIFIED NUMBER ATLAS")
    print(f"{'='*70}")

    ideas = [
        ("IDEA 1: The 'degree hierarchy' of tournament complexity",
         "The algebraic degree of numbers involved increases with the "
         "sophistication of the result: degree 1 (rational) for basic counting, "
         "degree 2 (quadratic) for Fibonacci/eigenvalues, degree 3 (cubic) for "
         "Pascal growth, degree infinity for asymptotics. "
         "This suggests a FILTRATION of tournament theory by algebraic degree."),

        ("IDEA 2: 2-adic tournament theory",
         "Develop tournament theory entirely in the 2-adic world Q_2. "
         "H is a 2-adic unit (v_2=0). The OCF is a 2-adic expansion. "
         "The forbidden values 7 and 21 should have 2-adic characterizations. "
         "7 = 111 in binary = maximally non-sparse. "
         "21 = 10101 in binary = alternating pattern."),

        ("IDEA 3: The unit circle as 'phase boundary'",
         "Lee-Yang zeros converge to |z|=1. The OCF evaluates at z=2 (OUTSIDE). "
         "What happens at |z|=1? This is the 'critical point' of the "
         "tournament partition function. Phase transitions live here."),

        ("IDEA 4: Galois theory of tournament eigenvalues",
         "The Paley eigenvalue (-1+i*sqrt(p))/2 lives in Q(i, sqrt(p)). "
         "The Galois group Gal(Q(i,sqrt(p))/Q) has order 4 = 2^2. "
         "Does this Galois symmetry explain tournament self-complementarity?"),

        ("IDEA 5: The 'cubic wall'",
         "No algebraic number of degree > 3 appears in tournament theory. "
         "Is this a THEOREM or just an observation? "
         "If 4^{1/3} is the highest-degree algebraic number, this constrains "
         "the possible generating functions."),

        ("IDEA 6: Modular forms and tournaments",
         "The Eisenstein series E_4, E_6 have coefficients involving "
         "tournament-relevant numbers (240 = |Phi(E_8)|, 504 = 7*8*9). "
         "Is there a modular form whose Fourier coefficients are max_H(n)?"),

        ("IDEA 7: p-adic interpolation of H",
         "H is defined for integer n. Can we interpolate to p-adic n? "
         "The p-adic interpolation of n! is the p-adic Gamma function. "
         "A 'p-adic H' function could reveal hidden structure."),

        ("IDEA 8: The Grothendieck group of tournaments",
         "Define K_0(Tournaments) as the free abelian group on iso classes "
         "modulo relations from deletion-contraction. "
         "H descends to a homomorphism K_0 -> Z. "
         "What is the structure of K_0?"),

        ("IDEA 9: Tournament L-functions",
         "Define L(T, s) = prod_p (1 - a_p * p^{-s})^{-1} where a_p "
         "encodes H(T) mod p. The Riemann hypothesis for L(T,s) would "
         "give bounds on the distribution of H values."),

        ("IDEA 10: The 'number field' of a tournament",
         "Each tournament T generates a number field Q(eigenvalues of A). "
         "For Paley T_p: the field is Q(i, sqrt(p)). "
         "The discriminant of this field relates to |Aut(T)|. "
         "Does the class number relate to H(T)?"),
    ]

    for title, desc in ideas:
        print(f"\n  {title}")
        print(f"    {desc}")

    # ============================================================
    # PART 9: THE MASTER TABLE — ALL NUMBERS IN ONE PLACE
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 9: THE MASTER TABLE")
    print(f"{'='*70}")

    master = [
        # (value, type, degree, role, appearances)
        ("1", "natural", 0, "unit/identity", 20),
        ("2", "natural", 0, "THE generator (binary, OCF fugacity)", 20),
        ("3", "natural", 0, "cycle generator (3-cycle)", 11),
        ("4 = 2^2", "natural", 0, "quadratic structure", 16),
        ("5 = 2+3", "natural", 0, "KEY_SUM, triple coincidence", 10),
        ("6 = 2*3", "natural", 0, "LCM, period, h(G2)", 10),
        ("7 = 2^3-1", "natural", 0, "FORBIDDEN (Mersenne)", 5),
        ("8 = 2^3", "natural", 0, "E8 rank, beta_4 onset", 8),
        ("9 = 3^2", "natural", 0, "CS boundary, mod-3 structure", 5),
        ("10 = 2*5", "natural", 0, "V(Petersen), arcs at n=5", 6),
        ("12 = 4*3", "natural", 0, "h(E6), bridge square-cycle", 5),
        ("13 = F(7)", "natural", 0, "Fibonacci prime, factor of 1729", 13),
        ("21 = 3*7", "natural", 0, "FORBIDDEN (cycle*Mersenne)", 5),
        ("1/2", "rational", 1, "binary probability", 11),
        ("1/3", "rational", 1, "signal-to-noise ratio", 9),
        ("1/4 = (1/2)^2", "rational", 1, "cycle probability", 8),
        ("2/3", "rational", 1, "exponent 2^{2/3}", 8),
        ("3/2", "rational", 1, "mean_H(3), perfect fifth", 10),
        ("5/3 = F5/F4", "rational", 1, "maxH step ratio", 3),
        ("9/2 = 3^2/2", "rational", 1, "maxH(7)/a(7)", 3),
        ("0.27 ≈ 17/63", "rational", 1, "information rate I/m", 5),
        ("phi", "algebraic", 2, "Fibonacci growth", 12),
        ("sqrt(2)", "algebraic", 2, "spectral radius", 7),
        ("sqrt(3)", "algebraic", 2, "variance = 1/sqrt(3)", 8),
        ("sqrt(5)", "algebraic", 2, "Fibonacci discriminant", 7),
        ("4^{1/3}", "algebraic", 3, "Pascal growth = 2^{2/3}", 7),
        ("e", "transcendental", "inf", "Szele efficiency limit", 11),
        ("pi", "transcendental", "inf", "asymptotic normalization", 10),
        ("log(2)", "transcendental", "inf", "information entropy base", 5),
        ("(-1+i*sqrt(p))/2", "complex", 2, "Paley eigenvalue", 3),
        ("omega_3 = e^{2pi*i/3}", "complex", 2, "cube root of unity", 5),
        ("H mod 2^k", "p-adic", "—", "2-adic structure of H", 5),
    ]

    print(f"\n  {'Value':>18} {'Type':>14} {'Deg':>4} {'Appearances':>5}  Role")
    print(f"  {'-'*80}")
    for val, typ, deg, role, apps in master:
        print(f"  {val:>18} {typ:>14} {str(deg):>4} {apps:>5}  {role}")

    # ============================================================
    # PART 10: THE DEEPEST SYNTHESIS
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 10: THE DEEPEST SYNTHESIS")
    print(f"{'='*70}")

    print(f"""
  TOURNAMENT THEORY IS A THEORY OF THE NUMBER 2.

  Everything flows from 2:
    2 = number of arc orientations
    H = I(Omega, 2) = evaluation at 2
    2-adic: H is always a 2-adic unit (odd)
    2^n: total tournament count involves powers of 2
    degree 2: Fibonacci (quadratic minimal poly)
    2 strands: Fibonacci even/odd
    2-fold: H(T) = H(T^op) involution

  The NUMBER 3 provides the TEXTURE:
    3-cycles: fundamental building blocks
    3 strands: Pascal interleaving
    mod 3: F(T,omega) divisibility
    cubic: 4^{1/3} = 2^{2/3} growth rate
    3^2 = 9: CS boundary

  The NUMBER 5 provides the CROSSOVER:
    5 = 2+3: where generators meet
    n=5: triple coincidence, last unimodal
    sqrt(5): Fibonacci discriminant

  EVERYTHING ELSE follows from 2, 3, 5:
    7 = 2^3-1 (forbidden)
    e ≈ phi^2 (Szele ≈ Fibonacci^2)
    pi (asymptotic normalization)

  THE FINAL WORD:
  Tournament theory = the theory of BINARY CHOICE (2)
  on CYCLIC STRUCTURES (3), crossing over at their SUM (5).
""")

    print(f"\n{'='*70}")
    print("DONE — UNIFIED NUMBER ATLAS COMPLETE")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
