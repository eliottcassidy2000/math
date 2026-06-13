"""
rational_archaeology.py -- kind-pasteur-2026-03-14-S93
LONG SESSION: Rational number archaeology in tournament theory.

Rationals are the CONNECTIVE TISSUE between integers.
Every ratio between two tournament quantities is a rational number
with potential structural meaning.

FOCUS: Where do specific rational numbers appear, and what do they mean?
"""

import sys, math
from collections import Counter, defaultdict
from fractions import Fraction

sys.stdout.reconfigure(encoding='utf-8')

def C(n, k):
    if k < 0 or k > n: return 0
    return math.comb(n, k)

maxH = {1:1, 2:1, 3:3, 4:5, 5:15, 6:45, 7:189, 8:661, 9:3357, 10:15745, 11:95095}
meanH = {n: Fraction(math.factorial(n), 2**(n-1)) for n in range(1, 12)}

def main():
    print("=" * 70)
    print("RATIONAL NUMBER ARCHAEOLOGY — THE CONNECTIVE TISSUE")
    print("kind-pasteur-2026-03-14-S93")
    print("=" * 70)

    # ============================================================
    # PART 1: THE FRACTION 1/2
    # ============================================================
    print(f"\n{'='*70}")
    print("THE FRACTION 1/2")
    print(f"{'='*70}")
    items = [
        "mean_H = n!/2^{n-1}: every factor of 2 in the denominator",
        "H is always odd => H/2 is never an integer (half-integer obstruction)",
        "The Fourier level-2 magnitude at n=3,4: (n-2)!/2^{n-2} = 1/2",
        "Var(H)/Mean^2 = E_nonconst/E_0 = (1/2)^{...}... no, = 1/3",
        "P(T flip = GS) = 1/2^{m - dof} per GS constraint",
        "Each arc independently contributes 1/2 to forward probability",
        "F(T,0) + F(T,n-1) <= 2, so (F(0)+F(n-1))/H <= 2/H ≈ 0 for large H",
        "The 'half-Redei' number T = (H-1)/2 is always an integer",
        "Paley eigenvalue: mu = (-1 +/- i*sqrt(p))/2 — the 1/2 in spectral theory",
        "The strand ratio in 3-strand Pascal: s1/s0 approaches 2, so s0/s1 → 1/2",
        "1/2 = probability of each arc direction in a random tournament",
    ]
    for item in items:
        print(f"  - {item}")
    print(f"  ESSENCE: 1/2 is the probability of each binary choice.")
    print(f"  It's the MOST FUNDAMENTAL rational in tournament theory.")

    # ============================================================
    # PART 2: THE FRACTION 1/3
    # ============================================================
    print(f"\n{'='*70}")
    print("THE FRACTION 1/3")
    print(f"{'='*70}")
    items = [
        "Var(H)/Mean(H)^2 ≈ 1/3 = the UNIVERSAL RATIO",
        "Fourier: E_nonconst/E_0 = 1/3 at n=3,4 (EXACT), ≈1/3 at n=5+",
        "1/3 of all vertex triples form 3-cycles in a random tournament",
        "The probability of a random 3-vertex subtournament being a cycle: 1/4... no, 2/8 = 1/4",
        "Actually: P(3-cycle) = 2/8 = 1/4 for random, but expected c3 = C(n,3)/4",
        "c3/C(n,3) = 1/4 for random tournament (not 1/3!)",
        "But the OCF coefficient is 2: H = 1 + 2*alpha_1, so alpha_1/H ≈ 1/2",
        "The GS fraction: #GS/#total = 2^dof/2^m approaches 0",
        "1/3 = the coefficient of variation squared when Fourier energy is 25/75",
    ]
    for item in items:
        print(f"  - {item}")
    print(f"  ESSENCE: 1/3 is the signal-to-noise ratio of H.")
    print(f"  The 'noise floor' (level 0) is 3x larger than the 'signal' (level 2+).")

    # ============================================================
    # PART 3: THE FRACTION 3/2 (= 1.5)
    # ============================================================
    print(f"\n{'='*70}")
    print("THE FRACTION 3/2")
    print(f"{'='*70}")
    items = [
        "mean_H(3) = 3!/2^2 = 6/4 = 3/2",
        "Z_M(2) = 3/2 = perfect fifth in music (opus discovery)",
        "3/2 = s2/s0 ratio... no, that's 2. But s1/s0 at large n approaches 2, not 3/2",
        "3/2 = the ratio max_H(3)/mean_H(3)... mean=3/2, max=3, ratio=2",
        "Fourier level-2 at n=5: (n-2)!/2^{n-2} = 3!/2^3 = 6/8 = 3/4. Not 3/2",
        "3/2 appears as the score of the middle vertex in transitive n=4 tournament",
        "In the 3-strand Pascal: consecutive strand ratios approach 2, 2, 1",
        "3/2 = the average of 1 and 2 (the two generators)",
        "The Szele ratio starts at 2.0 (n=3,5,6) and approaches e ≈ 2.718",
        "3/2 = Fraction(maxH[3], meanH[3].denominator)... = 3/(3/2) = 2",
    ]
    for item in items:
        print(f"  - {item}")
    print(f"  ESSENCE: 3/2 = mean_H(3) = the simplest non-trivial tournament mean.")
    print(f"  Also = perfect fifth = the most consonant interval after octave (2/1).")

    # ============================================================
    # PART 4: THE FRACTION 2/3
    # ============================================================
    print(f"\n{'='*70}")
    print("THE FRACTION 2/3")
    print(f"{'='*70}")
    items = [
        "2/3 = 1 - 1/3 = the complement of the universal ratio",
        "E_0/E_total ≈ 3/4 = 1 - 1/4... actually 75/100 = 3/4, not 2/3",
        "H mod 3: P(H ≡ 0 mod 3) ≈ 1/3, P(H ≢ 0) ≈ 2/3",
        "The GS code rate: dof/m ≈ 2/3 for large n",
        "4^{1/3} = 2^{2/3}: the cube root of 4 IS 2 raised to 2/3!",
        "2/3 is the EXPONENT that links the generator 2 to the 3-fold structure",
        "2/3 appears in the asymptotic: C(2n,n) ~ 4^n / sqrt(pi*n) and 4 = 2^2",
        "The info rate 0.27 ≈ 1/3.7 ≈ 2/(3*e^{0.3})",
    ]
    for item in items:
        print(f"  - {item}")
    print(f"  ESSENCE: 2/3 is the EXPONENT of the 2-3 interaction.")
    print(f"  4^(1/3) = 2^(2/3) links the binary generator to the ternary structure.")

    # ============================================================
    # PART 5: THE FRACTION 9/2 (= 4.5)
    # ============================================================
    print(f"\n{'='*70}")
    print("THE FRACTION 9/2")
    print(f"{'='*70}")
    items = [
        "max_H(7)/a(7) = 189/42 = 9/2 (EXACT!)",
        "189 = 3^3 * 7, 42 = 6*7 = 2*3*7, so 189/42 = 3^3/(2*3) = 3^2/2 = 9/2",
        "9/2 = 3^2/2 = (cycle generator)^2 / (binary generator)",
        "This ratio says: the maximizer at n=7 exceeds the 'vertex+triangle count'",
        "by exactly the square of the cycle structure divided by the binary choice",
        "9/2 = mean_H(5)/... no. mean_H(5) = 120/16 = 15/2 = 7.5",
        "Actually mean_H(5) = 15/2 = 7.5. And mean_H(3) = 3/2 = 1.5",
        "Ratio of means: mean_H(5)/mean_H(3) = (15/2)/(3/2) = 5. Not 9/2.",
        "But mean_H(6)/mean_H(4) = (720/32)/(24/8) = 22.5/3 = 15/2. Still not 9/2.",
    ]
    for item in items:
        print(f"  - {item}")
    print(f"  ESSENCE: 9/2 = 3^2/2^1 = the 'second power of cycle over first power of binary'.")
    print(f"  It's the ratio that connects the Paley maximizer to the vertex+triangle count.")

    # ============================================================
    # PART 6: SYSTEMATIC RATIO SCAN
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 6: SYSTEMATIC RATIO SCAN — ALL PAIRS OF TOURNAMENT QUANTITIES")
    print(f"{'='*70}")

    # Compute key ratios between tournament quantities
    quantities = {}
    for n in range(3, 12):
        quantities[f"maxH({n})"] = maxH.get(n, None)
        quantities[f"meanH({n})"] = float(meanH.get(n, 0))
        quantities[f"C({n},2)"] = C(n, 2)
        quantities[f"C({n},3)"] = C(n, 3)
        quantities[f"n!/{n}"] = math.factorial(n)

    # Find all ratios that are "nice" fractions (small numerator/denominator)
    nice_ratios = []

    for n in range(3, 9):
        mh = maxH[n]
        mn = float(meanH[n])
        arcs = C(n, 2)
        tris = C(n, 3)
        nfact = math.factorial(n)

        pairs = [
            (f"maxH({n})/meanH({n})", Fraction(mh, 1) / meanH[n]),
            (f"maxH({n})/C({n},2)", Fraction(mh, arcs) if arcs > 0 else None),
            (f"maxH({n})/C({n},3)", Fraction(mh, tris) if tris > 0 else None),
            (f"maxH({n})/n!", Fraction(mh, nfact)),
            (f"maxH({n})/{n}", Fraction(mh, n)),
            (f"meanH({n})/C({n},2)", meanH[n] / arcs if arcs > 0 else None),
            (f"C({n},3)/C({n},2)", Fraction(tris, arcs) if arcs > 0 else None),
        ]

        for name, ratio in pairs:
            if ratio is not None and isinstance(ratio, Fraction):
                if ratio.denominator <= 20 and ratio.numerator <= 200:
                    nice_ratios.append((name, ratio))

    print(f"\n  'Nice' ratios (small denominator <= 20):")
    for name, ratio in sorted(nice_ratios, key=lambda x: float(x[1])):
        print(f"    {name} = {ratio} = {float(ratio):.6f}")

    # ============================================================
    # PART 7: THE "FRACTION LATTICE" OF TOURNAMENT THEORY
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 7: THE FRACTION LATTICE")
    print(f"{'='*70}")

    # Collect all "structural" fractions
    structural_fractions = set()

    # Fourier magnitudes
    for n in range(3, 8):
        mag = Fraction(math.factorial(n-2), 2**(n-2))
        structural_fractions.add((f"|H_hat level 2| at n={n}", mag))

    # Mean H
    for n in range(3, 8):
        structural_fractions.add((f"mean_H({n})", meanH[n]))

    # maxH/meanH = Szele ratio
    for n in range(3, 8):
        ratio = Fraction(maxH[n], 1) / meanH[n]
        structural_fractions.add((f"maxH/meanH at n={n}", ratio))

    # GS fraction
    for n in range(3, 10):
        m = C(n-1, 2)
        dof = (m + (n-1)//2) // 2
        frac = Fraction(2**dof, 2**m)
        structural_fractions.add((f"GS fraction at n={n}", frac))

    # Info rate ≈ 0.27
    for n in [3, 4, 5]:
        m = C(n, 2)
        # Approximate: number of H values
        nH = {3: 2, 4: 3, 5: 7}
        if n in nH:
            # I/m ≈ log2(nH)/m
            info_approx = Fraction(int(1000 * math.log2(nH[n]) / m), 1000)
            structural_fractions.add((f"I/m approx at n={n}", info_approx))

    print(f"\n  Structural fractions sorted by value:")
    for name, frac in sorted(structural_fractions, key=lambda x: float(x[1])):
        print(f"    {frac} = {float(frac):.6f} [{name}]")

    # ============================================================
    # PART 8: THE FRACTION 1/4 — THE CYCLE PROBABILITY
    # ============================================================
    print(f"\n{'='*70}")
    print("THE FRACTION 1/4")
    print(f"{'='*70}")
    items = [
        "P(3 random vertices form a 3-cycle) = 2/8 = 1/4",
        "Expected c3 in random tournament: C(n,3)/4",
        "1/4 = (1/2)^2: the probability of two specific arc orientations",
        "The Fourier energy at level 2: 25% = 1/4 of total at n=3,4",
        "alpha_1 = c3 + c5 + ..., and E[c3] = C(n,3)/4",
        "The 75/25 energy split: the 25% = 1/4 signal fraction",
        "1/4 = det(Paley T_p) coefficient: ((p+1)/4)^{(p-1)/2}",
        "The '4' in 4^{1/3} (Pascal growth) = (1/4)^{-1}",
    ]
    for item in items:
        print(f"  - {item}")
    print(f"  ESSENCE: 1/4 = (1/2)^2 = the probability of a directed 3-cycle.")
    print(f"  It's the square of the fundamental binary probability.")

    # ============================================================
    # PART 9: FIBONACCI RATIOS
    # ============================================================
    print(f"\n{'='*70}")
    print("FIBONACCI RATIOS IN TOURNAMENT THEORY")
    print(f"{'='*70}")

    def fib(n):
        if n <= 0: return 0
        a, b = 0, 1
        for _ in range(n): a, b = b, a + b
        return a

    print(f"\n  F(n+1)/F(n) converges to phi:")
    for n in range(2, 15):
        ratio = Fraction(fib(n+1), fib(n))
        print(f"    F({n+1})/F({n}) = {ratio} = {float(ratio):.6f}")

    # Tournament quantities as Fibonacci ratios?
    print(f"\n  max_H ratios as Fibonacci ratios:")
    for n in range(3, 10):
        if n+1 in maxH:
            ratio = Fraction(maxH[n+1], maxH[n])
            # Is this close to any F(k+1)/F(k)?
            for k in range(2, 15):
                fib_ratio = Fraction(fib(k+1), fib(k))
                if abs(float(ratio) - float(fib_ratio)) < 0.1:
                    print(f"    maxH({n+1})/maxH({n}) = {ratio} ≈ F({k+1})/F({k}) = {fib_ratio}")

    # ============================================================
    # PART 10: THE RATIO 0.27 AS A FRACTION
    # ============================================================
    print(f"\n{'='*70}")
    print("THE 0.27 CONSTANT AS A FRACTION")
    print(f"{'='*70}")

    # 0.27 ≈ what simple fraction?
    target = 0.27
    best_frac = None
    best_diff = 1
    for d in range(1, 100):
        for num in range(1, d):
            diff = abs(num/d - target)
            if diff < best_diff:
                best_diff = diff
                best_frac = Fraction(num, d)

    print(f"  0.27 ≈ {best_frac} = {float(best_frac):.6f} (diff = {best_diff:.6f})")

    # More approximations
    print(f"\n  Rational approximations to 0.27:")
    for d in range(2, 20):
        num = round(target * d)
        frac = Fraction(num, d)
        diff = abs(float(frac) - target)
        if diff < 0.01:
            print(f"    {frac} = {float(frac):.6f} (diff = {diff:.6f})")

    # opus suggested 27/100 = 3^3/100 = 3^3/(4*25) = 3^3/(2^2*5^2)
    print(f"\n  27/100 = 3^3/(2^2 * 5^2) = {float(Fraction(27,100)):.6f}")
    print(f"  This encodes the 2-3-5 structure!")
    print(f"  The information rate I/m ≈ 27/100 = (cycle generator)^3 / (binary^2 * sum^2)")

    # ============================================================
    # PART 11: IDEA GENERATION FROM RATIONALS
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 11: NEW IDEAS FROM RATIONAL STRUCTURE")
    print(f"{'='*70}")

    ideas = [
        ("IDEA 1: The '1/3 theorem'",
         "PROVE that Var(H)/Mean(H)^2 → 1/3 as n → infinity. "
         "This should follow from the Fourier energy being 75/25 at level 0 vs level 2, "
         "which in turn follows from the level-2 coefficient formula (n-2)!/2^{n-2}."),

        ("IDEA 2: The '9/2 bridge'",
         "max_H(7)/a(7) = 9/2 = 3^2/2. Are there other n where max_H(n)/a(n) "
         "is a clean power of 3 divided by a power of 2? "
         "At n=3: 3/4 = 3/2^2. At n=5: 15/15 = 1 = 3^0/2^0. "
         "Pattern: max_H(2k+1)/a(2k+1) = 3^k/2^k = (3/2)^k?"),

        ("IDEA 3: The '1/4 fixed point'",
         "The cycle probability 1/4 is a FIXED POINT of the 'square' map "
         "x → x^2 at x = 1/2. Does this fixed point structure explain "
         "why the Fourier energy split is exactly 75/25 = 3:1?"),

        ("IDEA 4: Rational generating function",
         "Define G(x) = sum_{n>=1} max_H(n) * x^n. Is G(x) a rational function? "
         "If so, its poles would reveal the growth rate structure."),

        ("IDEA 5: The 'Stern-Brocot' position of tournament rationals",
         "Place the key tournament fractions on the Stern-Brocot tree. "
         "1/3 and 2/3 are neighbors. 1/2 is the root. "
         "The path from 1/2 to 1/3 goes L,L = two left turns."),

        ("IDEA 6: Continued fraction of 0.27",
         "0.27 = [0; 3, 1, 2, 1, ...]. The first convergent is 1/4 = 0.25. "
         "The second is 1/3 = 0.333. The '0.27' lives BETWEEN 1/4 and 1/3, "
         "exactly in the interval defined by the two fundamental rationals!"),

        ("IDEA 7: The 'rational spectrum' of Omega graphs",
         "For each graph G, I(G, 2) is an integer. But I(G, x) for rational x "
         "gives a rational number. The set {I(G, p/q) : G is a tournament's Omega} "
         "is a subset of rationals with specific structure."),
    ]

    for title, desc in ideas:
        print(f"\n  {title}")
        print(f"    {desc}")

    # Verify IDEA 2
    print(f"\n  VERIFYING IDEA 2: max_H(2k+1)/a(2k+1) = (3/2)^k?")
    for k in range(4):
        n = 2*k + 3  # n = 3, 5, 7, 9
        if n in maxH:
            an = n + C(n, 3)
            ratio = Fraction(maxH[n], an)
            predicted = Fraction(3, 2) ** k
            print(f"    k={k}, n={n}: maxH/a = {ratio} = {float(ratio):.4f}, "
                  f"(3/2)^{k} = {predicted} = {float(predicted):.4f}, "
                  f"match = {ratio == predicted}")

    # Verify IDEA 6
    print(f"\n  VERIFYING IDEA 6: Continued fraction of 0.27:")
    # 0.27 = 27/100
    frac = Fraction(27, 100)
    # Continued fraction expansion
    cf = []
    a, b = frac.numerator, frac.denominator
    while b > 0:
        q, r = divmod(a, b)
        cf.append(q)
        a, b = b, r
    print(f"    27/100 = {cf}")
    print(f"    Convergents:")
    h_prev, h_curr = 0, 1
    k_prev, k_curr = 1, 0
    for i, c in enumerate(cf):
        h_prev, h_curr = h_curr, c * h_curr + h_prev
        k_prev, k_curr = k_curr, c * k_curr + k_prev
        conv = Fraction(h_curr, k_curr)
        print(f"      [{', '.join(str(x) for x in cf[:i+1])}] = {conv} = {float(conv):.6f}")

    print(f"\n{'='*70}")
    print("DONE — DEEP RATIONAL ARCHAEOLOGY")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
